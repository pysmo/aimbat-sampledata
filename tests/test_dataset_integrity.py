"""Integrity checks for the committed `Event_*/` dataset itself.

These don't hit the network — they validate what's already on disk, so they
catch regressions like stale files left over from a previous station
selection, or a truncated seismogram (a real FDSN data gap) slipping into a
commit.
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pysmo.classes import SAC
from pysmo.tools.azdist import backazimuth

from fetch_events import (
    DELTA,
    DURATION_AFTER,
    EVENTS,
    MARGIN_BEFORE,
    MAX_P_DISTANCE_DEG,
    OUTPUT_DIR,
    PICK_HEADERS,
    PINNED_SITES,
    EventSpec,
    predicted_arrivals,
)

EVENT_DIRS = sorted(OUTPUT_DIR.glob("Event_*"))
SAC_GLOB = "*.BH?"
MIN_WINDOW = MARGIN_BEFORE + DURATION_AFTER
TRUNCATION_TOLERANCE = pd.Timedelta(seconds=2)
EVENT_DIR_PATTERN = re.compile(
    r"Event_(\d{4})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{3})$"
)
MANIFEST_FIELDS = [
    "network",
    "station",
    "location",
    "channel",
    "latitude",
    "longitude",
    "distance_deg",
    "azimuth_deg",
    "samples",
    "o",
    "idep",
    "incidence_deg",
]


def _manifest_rows(event_dir: Path) -> list[dict[str, str]]:
    with (event_dir / "manifest.csv").open(newline="") as fh:
        return list(csv.DictReader(fh))


def _sac_paths(event_dir: Path) -> list[Path]:
    return sorted(event_dir.glob(SAC_GLOB))


def _origin_time_from_dirname(event_dir: Path) -> pd.Timestamp:
    """The directory name is derived from `event.time` at fetch time (to
    millisecond precision, matching USGS's own resolution), so it's an
    independent ground truth to check the SAC files' `event.time` against.
    """
    match = EVENT_DIR_PATTERN.match(event_dir.name)
    assert match, f"unexpected directory name: {event_dir.name}"
    year, month, day, hour, minute, second, ms = (int(g) for g in match.groups())
    return pd.Timestamp(
        year=year,
        month=month,
        day=day,
        hour=hour,
        minute=minute,
        second=second,
        microsecond=ms * 1000,
        tz="UTC",
    )


def test_one_directory_per_event() -> None:
    assert len(EVENT_DIRS) == len(EVENTS)


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_manifest_matches_files_on_disk(event_dir: Path) -> None:
    """No orphaned SAC files, and no manifest row missing its file."""
    rows = _manifest_rows(event_dir)
    manifest_files = {
        f"{r['network']}.{r['station']}.{r['location']}.{r['channel']}" for r in rows
    }
    disk_files = {p.name for p in _sac_paths(event_dir)}
    assert manifest_files == disk_files


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_manifest_has_expected_columns(event_dir: Path) -> None:
    rows = _manifest_rows(event_dir)
    assert rows, f"{event_dir.name}: manifest.csv has no rows"
    assert list(rows[0].keys()) == MANIFEST_FIELDS


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_only_pinned_sites(event_dir: Path) -> None:
    sites = {f"{r['network']}.{r['station']}" for r in _manifest_rows(event_dir)}
    assert sites <= set(PINNED_SITES), sorted(sites - set(PINNED_SITES))


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_every_site_has_every_channel(event_dir: Path) -> None:
    """A site is kept only with all of the event's output channels."""
    channels: defaultdict[str, set[str]] = defaultdict(set)
    for row in _manifest_rows(event_dir):
        channels[f"{row['network']}.{row['station']}"].add(row["channel"])
    expected = set(_spec(event_dir).outputs)
    for site, found in channels.items():
        assert found == expected, f"{site}: {sorted(found)}"


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_p_event_stations_within_shadow_zone(event_dir: Path) -> None:
    if "P" not in _spec(event_dir).phases:
        pytest.skip("no P phase")
    for row in _manifest_rows(event_dir):
        assert float(row["distance_deg"]) < MAX_P_DISTANCE_DEG


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_seismograms_are_not_truncated(event_dir: Path) -> None:
    """A short prefix instead of the full window (a real FDSN data gap) may
    not even reach the first arrival — useless for ICCS/MCCC and confusing
    in AIMBAT.
    """
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        duration = sac.seismogram.end_time - sac.seismogram.begin_time
        assert duration >= MIN_WINDOW - TRUNCATION_TOLERANCE, (
            f"{sac_path.name}: only {duration} of data, expected at least ~{MIN_WINDOW}"
        )


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_channels_of_a_site_share_the_window(event_dir: Path) -> None:
    windows: defaultdict[str, set[tuple[pd.Timestamp, int]]] = defaultdict(set)
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        site = sac_path.name.rsplit(".", 1)[0]
        begin = sac.seismogram.begin_time.round(DELTA)
        windows[site].add((begin, len(sac.seismogram.data)))
    for site, found in windows.items():
        assert len(found) == 1, f"{site}: {found}"


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_sampling_interval(event_dir: Path) -> None:
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        assert sac.seismogram.delta == DELTA, f"{sac_path.name}: {sac.seismogram.delta}"


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_orientation_is_set(event_dir: Path) -> None:
    """Z points up, T lies flat at the back azimuth + 270, and L and Q lie in
    the vertical plane through the ray at the back azimuth + 180, with dips
    `incidence - 90` and `incidence`.
    """
    incidences = {
        f"{r['network']}.{r['station']}.{r['location']}": r["incidence_deg"]
        for r in _manifest_rows(event_dir)
    }
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        channel = _channel(sac_path)
        if channel == "BHZ":
            assert sac.orientation.dip == pytest.approx(-90.0), sac_path.name
            continue
        baz = backazimuth(sac.event, sac.station)
        if channel == "BHT":
            expected = ((baz + 270) % 360, 0.0)
        else:
            incidence = float(incidences[sac_path.name.rsplit(".", 1)[0]])
            dip = incidence - 90 if channel == "BHL" else incidence
            expected = ((baz + 180) % 360, dip)
        found = (sac.orientation.azimuth, sac.orientation.dip)
        assert found == pytest.approx(expected, abs=1e-3), sac_path.name


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_q_is_rotated_from_z_and_r(event_dir: Path) -> None:
    """Q = -Z sin(i) + R cos(i), with R recovered from Q and L."""
    spec = _spec(event_dir)
    if "BHQ" not in spec.outputs:
        pytest.skip("no BHQ")
    for row in _manifest_rows(event_dir):
        if row["channel"] != "BHQ":
            continue
        site = f"{row['network']}.{row['station']}.{row['location']}"
        z, l_, q = (
            np.asarray(SAC.from_file(event_dir / f"{site}.{c}").seismogram.data)
            for c in ("BHZ", "BHL", "BHQ")
        )
        i = np.radians(float(row["incidence_deg"]))
        # L = Z cos i + R sin i, so R follows from Z and L.
        r = (l_ - z * np.cos(i)) / np.sin(i)
        expected = -z * np.sin(i) + r * np.cos(i)
        np.testing.assert_allclose(q, expected, atol=1e-4 * np.abs(q).max())


# SAC v7's double-precision footer still leaves a few ns of float64 rounding
# noise (measured: ~1ns) when converting a "seconds offset from a reference
# time" header back into a timestamp. That's expected and harmless. v6's
# float32 `o` header, by contrast, visibly truncates at the millisecond
# level or worse for offsets in this range — many orders of magnitude
# bigger than float64 noise, and exactly what nvhdr=7 exists to avoid.
EVENT_TIME_TOLERANCE = pd.Timedelta(microseconds=1)


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_event_time_is_read_back_precisely(event_dir: Path) -> None:
    """Guards against SAC v6's float32 `o` header silently regressing back
    in — that's why `fetch_events.py` writes v7 (`sac.native.nvhdr = 7`),
    which adds a double-precision footer.
    """
    expected = _origin_time_from_dirname(event_dir)
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        assert sac.native.nvhdr == 7, f"{sac_path.name}: not SAC v7"
        drift = abs(sac.event.time - expected)
        assert drift <= EVENT_TIME_TOLERANCE, (
            f"{sac_path.name}: event time {sac.event.time} != {expected} "
            f"(off by {drift})"
        )


def test_shared_station_overlap() -> None:
    """Regression guard for the bug this dataset redesign fixed: selecting
    stations independently per event barely overlapped station-to-station
    across events, defeating the point of a multi-event ICCS/MCCC dataset.
    """
    events_per_site: defaultdict[str, set[str]] = defaultdict(set)
    for event_dir in EVENT_DIRS:
        for row in _manifest_rows(event_dir):
            events_per_site[f"{row['network']}.{row['station']}"].add(event_dir.name)
    counts = Counter({site: len(events) for site, events in events_per_site.items()})

    shared = sum(1 for n in counts.values() if n >= 2)
    assert shared / len(counts) > 0.5, (
        "fewer than half of all stations are shared across multiple events"
    )


EVENT_LABELS = {
    "Event_2014.04.12.20.14.39.300": "solomon_islands",
    "Event_2014.08.24.23.21.45.520": "central_peru",
    "Event_2014.11.01.18.57.22.380": "fiji_region",
    "Event_2015.04.25.06.11.25.950": "nepal",
    "Event_2017.07.17.23.34.13.740": "komandorskiye_ostrova",
    "Event_2017.11.12.18.18.17.180": "iraq",
}
PICK_TOLERANCE = pd.Timedelta(milliseconds=100)
TAUP_PHASES = {"pP", "SKS"}


def _spec(event_dir: Path) -> EventSpec:
    return next(spec for spec in EVENTS if spec.label == EVENT_LABELS[event_dir.name])


def _channel(sac_path: Path) -> str:
    return sac_path.name.rsplit(".", 1)[1]


def test_every_event_directory_has_a_label() -> None:
    assert sorted(EVENT_LABELS) == [p.name for p in EVENT_DIRS]


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_every_phase_is_picked_on_every_channel(event_dir: Path) -> None:
    """Each channel of a site carries every pick of the event, and no other."""
    spec = _spec(event_dir)
    for sac_path in _sac_paths(event_dir):
        phases = spec.phases
        sac = SAC.from_file(sac_path)
        for phase, header in PICK_HEADERS.items():
            time = getattr(sac.timestamps, header)
            if phase in phases:
                assert time is not None, f"{sac_path.name}: no {phase} pick"
                label = getattr(sac.native, f"k{header}").strip()
                assert label == phase, f"{sac_path.name}: k{header} is {label!r}"
            else:
                assert time is None, f"{sac_path.name}: unexpected {phase} pick"


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_picks_track_their_backends(event_dir: Path) -> None:
    spec = _spec(event_dir)
    if TAUP_PHASES & set(spec.phases):
        pytest.importorskip("obspy")
    for sac_path in _sac_paths(event_dir):
        sac = SAC.from_file(sac_path)
        arrivals = predicted_arrivals(spec, sac.event, sac.station)
        for phase in spec.phases:
            pick = getattr(sac.timestamps, PICK_HEADERS[phase])
            assert abs(pick - arrivals[phase]) < PICK_TOLERANCE, (
                f"{sac_path.name}: {phase}"
            )
