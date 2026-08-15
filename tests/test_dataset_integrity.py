"""Integrity checks for the committed `Event_*/` dataset itself.

These don't hit the network — they validate what's already on disk, so they
catch regressions like stale files left over from a previous station
selection, or a truncated seismogram (a real FDSN data gap) slipping into a
commit.
"""

import csv
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import pytest

from pysmo.classes import SAC

from fetch_events import (
    DURATION_AFTER,
    EVENTS,
    MARGIN_BEFORE,
    MAX_P_DISTANCE_DEG,
    OUTPUT_DIR,
)

EVENT_DIRS = sorted(OUTPUT_DIR.glob("Event_*"))
EXPECTED_WINDOW = MARGIN_BEFORE + DURATION_AFTER
TRUNCATION_TOLERANCE = pd.Timedelta(seconds=2)
EVENT_DIR_PATTERN = re.compile(
    r"Event_(\d{4})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{2})\.(\d{3})$"
)
MANIFEST_FIELDS = [
    "network",
    "station",
    "location",
    "latitude",
    "longitude",
    "distance_deg",
    "azimuth_deg",
    "samples",
    "o",
    "idep",
]


def _manifest_rows(event_dir: Path) -> list[dict[str, str]]:
    with (event_dir / "manifest.csv").open(newline="") as fh:
        return list(csv.DictReader(fh))


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
    """No orphaned `.BHZ` files, and no manifest row missing its file."""
    rows = _manifest_rows(event_dir)
    manifest_files = {
        f"{r['network']}.{r['station']}.{r['location']}.BHZ" for r in rows
    }
    disk_files = {p.name for p in event_dir.glob("*.BHZ")}
    assert manifest_files == disk_files


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_manifest_has_expected_columns(event_dir: Path) -> None:
    rows = _manifest_rows(event_dir)
    assert rows, f"{event_dir.name}: manifest.csv has no rows"
    assert list(rows[0].keys()) == MANIFEST_FIELDS


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_stations_within_shadow_zone(event_dir: Path) -> None:
    for row in _manifest_rows(event_dir):
        assert float(row["distance_deg"]) < MAX_P_DISTANCE_DEG


@pytest.mark.parametrize("event_dir", EVENT_DIRS, ids=lambda p: p.name)
def test_seismograms_are_not_truncated(event_dir: Path) -> None:
    """A short prefix instead of the full window (a real FDSN data gap) may
    not even reach the P arrival — useless for ICCS/MCCC and confusing in
    AIMBAT.
    """
    for sac_path in sorted(event_dir.glob("*.BHZ")):
        sac = SAC.from_file(sac_path)
        duration = sac.seismogram.end_time - sac.seismogram.begin_time
        assert duration >= EXPECTED_WINDOW - TRUNCATION_TOLERANCE, (
            f"{sac_path.name}: only {duration} of data, expected ~{EXPECTED_WINDOW}"
        )


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
    for sac_path in sorted(event_dir.glob("*.BHZ")):
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
    counts: Counter[str] = Counter()
    for event_dir in EVENT_DIRS:
        for row in _manifest_rows(event_dir):
            counts[f"{row['network']}.{row['station']}"] += 1

    shared = sum(1 for n in counts.values() if n >= 2)
    assert shared / len(counts) > 0.5, (
        "fewer than half of all stations are shared across multiple events"
    )
