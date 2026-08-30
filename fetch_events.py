"""Fetch this repo's Alaska teleseismic sample dataset from USGS + EarthScope.

Downloads 5 teleseismic earthquakes recorded by Alaska's `AK`/`AV`/`TA`
broadband networks (BHZ only), removes the instrument response, and
annotates event metadata and an initial P-pick. See PROVENANCE.md for the
event list, station-selection methodology, and window formula.

Re-run with `uv run python fetch_events.py` from this directory to
regenerate the dataset from scratch.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pysmo import MiniEvent, MiniStation
from pysmo.classes import SAC, StationXML
from pysmo.functions import detrend, taper
from pysmo.lib.io import (
    DEFAULT_REQUEST_RETRIES,
    DEFAULT_RETRY_DELAY_SECONDS,
    DEFAULT_TIMEOUT_SECONDS,
    http_get,
)
from pysmo.tools.azdist import azimuth, haversine
from pysmo.tools.signal import remove_response
from pysmo.tools.traveltime import travel_times
from pysmo.tools.web import fetch_stationxml

OUTPUT_DIR = Path(__file__).parent

# pysmo has no station-discovery wrapper, so this speaks raw FDSN text.
STATION_URL = "https://service.earthscope.org/fdsnws/station/1/query"
USGS_EVENT_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"

STATIONS_PER_EVENT = 40

# Stays clear of the P shadow zone (~100 deg onwards) so every station
# anchors on the same phase.
MAX_P_DISTANCE_DEG = 95.0

# P-coda window: covers pick/QC headroom before P and the coda after,
# well short of S and surface waves.
MARGIN_BEFORE = pd.Timedelta(minutes=2)
DURATION_AFTER = pd.Timedelta(minutes=3)


@dataclass(frozen=True)
class EventSpec:
    label: str
    usgs_eventid: str


EVENTS = [
    EventSpec("komandorskiye_ostrova", "us20009x42"),
    EventSpec("fiji_region", "usc000stdc"),
    EventSpec("nepal", "us20002926"),
    EventSpec("iraq", "us2000bmcg"),
    EventSpec("solomon_islands", "usc000phx5"),
]


def _fetch_event(eventid: str) -> MiniEvent:
    data = http_get(
        USGS_EVENT_URL,
        {"eventid": eventid, "format": "geojson"},
        timeout_seconds=DEFAULT_TIMEOUT_SECONDS,
        request_retries=DEFAULT_REQUEST_RETRIES,
        retry_delay_seconds=DEFAULT_RETRY_DELAY_SECONDS,
    )
    feature = json.loads(data)
    lon, lat, depth_km = feature["geometry"]["coordinates"]
    origin_ms = feature["properties"]["time"]
    return MiniEvent(
        latitude=lat,
        longitude=lon,
        depth=depth_km * 1000.0,
        time=pd.Timestamp(origin_ms, unit="ms", tz="UTC"),
    )


_IDEP_BY_UNITS = {"m": "disp", "m/s": "vel", "m/s**2": "acc"}


def _deconvolve(sac: SAC, station: MiniStation) -> str:
    """Remove the instrument response from `sac.seismogram`, in place.

    Returns the SAC `idep` value matching `response.input_units`, for the
    caller to set.
    """
    seismogram = sac.seismogram
    xml = fetch_stationxml(station=station)
    response = StationXML.from_bytes(xml, time=seismogram.begin_time).response

    nyquist = 0.5 / seismogram.delta.total_seconds()
    stage_nyquist = min(
        (stage.input_sample_rate / 2 for stage in response.stages), default=nyquist
    )
    f4 = 0.8 * min(nyquist, stage_nyquist)
    f3 = f4 * 0.9
    f1 = min(abs(pole) for pole in response.poles if pole != 0) / 10
    f2 = f1 * 10

    detrend(seismogram)
    taper(seismogram, 0.05)
    remove_response(seismogram, response, pre_filt=(f1, f2, f3, f4))
    return _IDEP_BY_UNITS.get(response.input_units.lower(), "unkn")


def _discover_stations(origin_time: pd.Timestamp) -> list[MiniStation]:
    """Discover unique AK/AV/TA BHZ stations operating at `origin_time`."""
    window_start = (origin_time - pd.Timedelta(days=1)).floor("s")
    window_end = (origin_time + pd.Timedelta(days=1)).floor("s")
    text = http_get(
        STATION_URL,
        {
            "net": "AK,AV,TA",
            "cha": "BHZ",
            "level": "channel",
            "format": "text",
            # Excludes "TA"'s earlier Lower-48 deployment; AK/AV are all >=53N.
            "minlatitude": 50,
            "starttime": window_start.isoformat(),
            "endtime": window_end.isoformat(),
        },
        timeout_seconds=DEFAULT_TIMEOUT_SECONDS,
        request_retries=DEFAULT_REQUEST_RETRIES,
        retry_delay_seconds=DEFAULT_RETRY_DELAY_SECONDS,
    ).decode("ascii")

    # One row per operational epoch; dedupe to the first per station.
    seen: dict[str, MiniStation] = {}
    for line in text.splitlines()[1:]:
        net, sta, loc, cha, lat, lon = line.split("|")[:6]
        key = f"{net}.{sta}"
        if key in seen:
            continue
        # dataselect rejects a blank location (HTTP 422); FDSN's "--" works.
        loc = loc.strip() or "--"
        seen[key] = MiniStation(
            name=sta,
            network=net,
            location=loc,
            channel=cha,
            latitude=float(lat),
            longitude=float(lon),
        )
    return list(seen.values())


def _select_subset(stations: list[MiniStation], count: int) -> list[MiniStation]:
    """Deterministically pick `count` stations spread across the array by latitude."""
    ordered = sorted(stations, key=lambda s: (s.latitude, s.network, s.name))
    if count >= len(ordered):
        return ordered
    # linspace, not a floor-divided stride: a stride degrades to 1 (no
    # spread at all) once the pool isn't much bigger than count.
    raw_indices = (round(i) for i in np.linspace(0, len(ordered) - 1, count))
    unique_indices = sorted(dict.fromkeys(raw_indices))
    return [ordered[i] for i in unique_indices]


def _select_common_stations(
    events: list[MiniEvent], count: int = STATIONS_PER_EVENT
) -> list[MiniStation]:
    """Pick one station set shared by every given event, not one per event.

    Selecting independently per event (each against its own differently
    shaped pool) barely overlaps station-to-station even when the same
    physical station is available for both — the linspace lands on
    different indices. Multiple events per station is valuable for this
    test dataset (ICCS/MCCC exercise the same station across events), so
    intersect first and spread by latitude once.
    """
    pools = [
        {f"{s.network}.{s.name}": s for s in _discover_stations(event.time)}
        for event in events
    ]
    common_keys = set.intersection(*(set(pool) for pool in pools))
    candidates = [
        station
        for key in common_keys
        for station in [next(pool[key] for pool in pools if key in pool)]
        if all(haversine(event, station) < MAX_P_DISTANCE_DEG for event in events)
    ]
    return _select_subset(candidates, count)


def _fetch_one_event(
    spec: EventSpec, event: MiniEvent, stations: list[MiniStation]
) -> None:
    event_dir = OUTPUT_DIR / (
        f"Event_{event.time:%Y.%m.%d.%H.%M.%S}.{event.time.microsecond // 1000:03d}"
    )
    event_dir.mkdir(exist_ok=True)

    manifest_rows: list[dict[str, object]] = []
    for station in stations:
        dist_deg = haversine(event, station)
        try:
            arrivals = travel_times(depth=event.depth, distance=dist_deg, phases=["P"])
            predicted_p = event.time + arrivals["P"]
            starttime = predicted_p - MARGIN_BEFORE
            endtime = predicted_p + DURATION_AFTER
            sac = SAC.fetch(station=station, starttime=starttime, endtime=endtime)
            # dataselect can return a short prefix instead of erroring outright
            # when the station has a data gap partway through the window.
            if sac.seismogram.end_time < endtime - pd.Timedelta(seconds=1):
                raise RuntimeError(
                    f"incomplete data: got {sac.seismogram.begin_time} to "
                    f"{sac.seismogram.end_time}, requested up to {endtime}"
                )
            idep = _deconvolve(sac, station)
        except Exception as exc:  # noqa: BLE001 - skip and log, don't abort the run
            print(f"skip {station.network}.{station.name}: {exc}")
            continue

        sac.native.idep = idep
        # v7 for double-precision o/t0 (float32 visibly truncates predicted_p).
        sac.native.nvhdr = 7
        sac.event.latitude = event.latitude
        sac.event.longitude = event.longitude
        sac.event.depth = event.depth
        sac.event.time = event.time
        assert sac.native.o is not None  # populated as soon as event.time is set

        sac.timestamps.t0 = predicted_p  # initial pick, for ICCS-style workflows

        loc = station.location.strip()
        filename = f"{station.network}.{station.name}.{loc}.BHZ"
        sac.write(event_dir / filename)

        manifest_rows.append(
            {
                "network": station.network,
                "station": station.name,
                "location": loc,
                "latitude": station.latitude,
                "longitude": station.longitude,
                "distance_deg": round(dist_deg, 4),
                "azimuth_deg": round(azimuth(event, station), 4),
                "samples": len(sac.seismogram.data),
                "o": round(sac.native.o, 3),
                "idep": idep,
            }
        )
        print(f"fetched {filename}: {starttime} to {endtime}")

    if not manifest_rows:
        raise RuntimeError(f"{spec.label}: no stations fetched, check geometry")
    manifest_path = event_dir / "manifest.csv"
    with manifest_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"{spec.label}: {len(manifest_rows)}/{len(stations)} stations fetched")


def main() -> None:
    fetched = [(spec, _fetch_event(spec.usgs_eventid)) for spec in EVENTS]

    shared_stations = _select_common_stations([event for _, event in fetched])
    for spec, event in fetched:
        _fetch_one_event(spec, event, shared_stations)


if __name__ == "__main__":
    main()
