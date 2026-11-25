"""Fetch this repo's Alaska teleseismic sample dataset from USGS + EarthScope.

Downloads 6 teleseismic earthquakes recorded by Alaska's `AK`/`AV`/`TA`
broadband networks, removes the instrument response, resamples to 20 Hz,
and annotates event metadata, channel orientation and labelled initial
picks: P on BHZ for five events (pP too for the deep Fiji event), S on
BHN/BHE for four of them, and SKS on BHN/BHE for a sixth event in the core
shadow. See PROVENANCE.md for the event list, station-selection
methodology, and window formula.

Re-run with `uv run --group fetch python fetch_events.py` from this
directory to regenerate the dataset from scratch. The `fetch` group brings
ObsPy, which computes pP and SKS.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from attrs import evolve

from pysmo import Event, MiniEvent, MiniStation, Station
from pysmo.classes import SAC, StationXML
from pysmo.exceptions import MissingValueError
from pysmo.functions import detrend, resample, taper
from pysmo.lib.io import (
    DEFAULT_REQUEST_RETRIES,
    DEFAULT_RETRY_DELAY_SECONDS,
    DEFAULT_TIMEOUT_SECONDS,
    http_get,
)
from pysmo.tools.azdist import azimuth, haversine
from pysmo.tools.signal import lowpass, remove_response
from pysmo.tools.web import fetch_stationxml

from traveltime_backends import arrival_times

OUTPUT_DIR = Path(__file__).parent

# pysmo has no station-discovery wrapper, so this speaks raw FDSN text.
STATION_URL = "https://service.earthscope.org/fdsnws/station/1/query"
USGS_EVENT_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"

STATIONS_PER_EVENT = 40

# The sites `_select_common_stations` chose in November 2026. Pinned
# because its result depends on what the station service lists at fetch
# time, which changes; PROVENANCE.md and testkit's fixtures name these
# sites. Re-run the selection only to change the set on purpose.
PINNED_SITES = (
    "AK.ATKA", "AK.BAL", "AK.BWN", "AK.CAPN", "AK.CAST", "AK.CCB", "AK.CHI",
    "AK.CHUM", "AK.EYAK", "AK.FIRE", "AK.FYU", "AK.GRNC", "AK.HIN", "AK.JIS",
    "AK.KAI", "AK.KHIT", "AK.KIAG", "AK.KLU", "AK.KNK", "AK.MCK", "AK.NICH",
    "AK.PAX", "AK.PIN", "AK.PNL", "AK.SAMH", "AK.SKN", "AK.YAH", "AV.AKRB",
    "AV.AU22", "AV.AUJA", "AV.MAPS", "AV.MSW", "AV.OKFG", "AV.RDDF", "AV.RDWB",
    "AV.SPCN", "AV.SPCP", "AV.SSBA", "TA.POKR", "TA.TOLK",
)  # fmt: skip

# Stays clear of the P shadow zone (~100 deg onwards) so every station of
# a P event anchors on the same phase.
MAX_P_DISTANCE_DEG = 95.0

# Window from MARGIN_BEFORE before an event's first phase to DURATION_AFTER
# after its last: pick/QC headroom before, the coda after.
MARGIN_BEFORE = pd.Timedelta(minutes=2)
DURATION_AFTER = pd.Timedelta(minutes=3)

# Sampling interval of every trace. 20 Hz is ample for P aligned at up to
# 2 Hz and for S and SKS below 0.2 Hz.
DELTA = pd.Timedelta(seconds=0.05)
# Below the new Nyquist (10 Hz): resample's Fourier method cuts the
# spectrum hard, so the low-pass keeps it from ringing.
RESAMPLE_LOWPASS_HZ = 8.0

# SAC header holding each phase's initial pick. The label goes in the
# matching `k` header. t1 to t3 stay free for AIMBAT v1.
PICK_HEADERS = {"P": "t0", "pP": "t5", "S": "t6", "SKS": "t7"}

# Phases picked on the horizontals; the others are picked on BHZ.
HORIZONTAL_PHASES = frozenset({"S", "SKS"})

VERTICAL = ("BHZ",)
THREE_COMPONENTS = ("BHZ", "BHN", "BHE")


@dataclass(frozen=True)
class EventSpec:
    label: str
    usgs_eventid: str
    phases: tuple[str, ...] = ("P",)
    channels: tuple[str, ...] = VERTICAL


EVENTS = [
    EventSpec("komandorskiye_ostrova", "us20009x42"),
    # 434 km deep: pP arrives 94 to 99 s after P.
    EventSpec(
        "fiji_region",
        "usc000stdc",
        phases=("P", "pP", "S"),
        channels=THREE_COMPONENTS,
    ),
    EventSpec("nepal", "us20002926", phases=("P", "S"), channels=THREE_COMPONENTS),
    EventSpec("iraq", "us2000bmcg", phases=("P", "S"), channels=THREE_COMPONENTS),
    EventSpec(
        "solomon_islands",
        "usc000phx5",
        phases=("P", "S"),
        channels=THREE_COMPONENTS,
    ),
    # 88 to 108 deg from the sites: most are in the core shadow, so no P.
    EventSpec("central_peru", "usb000s5x1", phases=("SKS",), channels=THREE_COMPONENTS),
]


def phases_on_channel(spec: EventSpec, channel: str) -> tuple[str, ...]:
    """The event's phases picked on `channel`."""
    horizontal = channel[-1] != "Z"
    return tuple(p for p in spec.phases if (p in HORIZONTAL_PHASES) == horizontal)


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


def _deconvolve(sac: SAC, station: MiniStation) -> tuple[str, StationXML]:
    """Remove the instrument response from `sac.seismogram`, in place.

    Returns the SAC `idep` value matching `response.input_units`, for the
    caller to set, and the channel epoch, which carries the orientation.
    """
    seismogram = sac.seismogram
    xml = fetch_stationxml(station=station)
    epoch = StationXML.from_bytes(xml, time=seismogram.begin_time)
    response = epoch.response

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
    return _IDEP_BY_UNITS.get(response.input_units.lower(), "unkn"), epoch


def _set_orientation(sac: SAC, epoch: StationXML) -> None:
    """Copy the epoch's orientation, leaving the headers unset if it has none."""
    try:
        azimuth_deg, dip_deg = epoch.azimuth, epoch.dip
    except MissingValueError:
        return
    sac.orientation.azimuth = azimuth_deg
    sac.orientation.dip = dip_deg


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


def _pinned_stations(events: list[MiniEvent]) -> list[MiniStation]:
    """Look up the `PINNED_SITES` among the stations operating at the events.

    Raises:
        RuntimeError: If a pinned site is not listed at any of the events.
    """
    pool: dict[str, MiniStation] = {}
    for event in events:
        for station in _discover_stations(event.time):
            pool.setdefault(f"{station.network}.{station.name}", station)
    missing = [site for site in PINNED_SITES if site not in pool]
    if missing:
        raise RuntimeError(f"pinned sites not listed: {', '.join(missing)}")
    return [pool[site] for site in PINNED_SITES]


def predicted_arrivals(
    spec: EventSpec, event: Event, station: Station
) -> dict[str, pd.Timestamp]:
    """Predicted arrival time of each of the event's phases at a station."""
    travel = arrival_times(
        depth=event.depth, distance=haversine(event, station), phases=spec.phases
    )
    return {phase: event.time + time for phase, time in travel.items()}


def set_picks(
    sac: SAC, arrivals: dict[str, pd.Timestamp], phases: tuple[str, ...]
) -> None:
    """Write the arrival of each of `phases` as a labelled initial pick (`PICK_HEADERS`)."""
    for phase in phases:
        time = arrivals[phase]
        header = PICK_HEADERS[phase]
        setattr(sac.timestamps, header, time)
        setattr(sac.native, f"k{header}", phase)


def _fetch_channel(
    station: MiniStation, starttime: pd.Timestamp, endtime: pd.Timestamp
) -> tuple[SAC, str]:
    """Fetch one channel, remove its response and resample it to `DELTA`.

    Returns the SAC and its `idep`.

    Raises:
        RuntimeError: If the data are incomplete.
    """
    sac = SAC.fetch(station=station, starttime=starttime, endtime=endtime)
    # dataselect can return a short prefix instead of erroring outright
    # when the station has a data gap partway through the window.
    if sac.seismogram.end_time < endtime - pd.Timedelta(seconds=1):
        raise RuntimeError(
            f"incomplete data: got {sac.seismogram.begin_time} to "
            f"{sac.seismogram.end_time}, requested up to {endtime}"
        )
    idep, epoch = _deconvolve(sac, station)
    lowpass(sac.seismogram, RESAMPLE_LOWPASS_HZ, corners=4, zerophase=True)
    resample(sac.seismogram, DELTA)
    _set_orientation(sac, epoch)
    return sac, idep


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
        # A site is kept only with every channel of the event, all sharing
        # one window.
        try:
            arrivals = predicted_arrivals(spec, event, station)
            starttime = min(arrivals.values()) - MARGIN_BEFORE
            endtime = max(arrivals.values()) + DURATION_AFTER
            fetched = [
                (
                    channel,
                    *_fetch_channel(
                        evolve(station, channel=channel), starttime, endtime
                    ),
                )
                for channel in spec.channels
            ]
        except Exception as exc:  # noqa: BLE001 - skip and log, don't abort the run
            print(f"skip {station.network}.{station.name}: {exc}")
            continue

        loc = station.location.strip()
        for channel, sac, idep in fetched:
            sac.native.idep = idep
            sac.event.latitude = event.latitude
            sac.event.longitude = event.longitude
            sac.event.depth = event.depth
            sac.event.time = event.time
            assert sac.native.o is not None  # populated as soon as event.time is set

            set_picks(sac, arrivals, phases_on_channel(spec, channel))

            filename = f"{station.network}.{station.name}.{loc}.{channel}"
            # v7 for double-precision o/t0 (float32 visibly truncates the picks).
            sac.write(event_dir / filename, nvhdr=7)

            manifest_rows.append(
                {
                    "network": station.network,
                    "station": station.name,
                    "location": loc,
                    "channel": channel,
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
    print(f"{spec.label}: {len(manifest_rows)} files fetched")


def main() -> None:
    fetched = [(spec, _fetch_event(spec.usgs_eventid)) for spec in EVENTS]

    # The P events locate the pinned sites; the SKS event uses the same ones.
    shared_stations = _pinned_stations(
        [event for spec, event in fetched if "P" in spec.phases]
    )
    for spec, event in fetched:
        _fetch_one_event(spec, event, shared_stations)


if __name__ == "__main__":
    main()
