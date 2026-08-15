"""Unit tests for the pure (non-network) station-selection helpers."""

from pysmo import MiniStation

from fetch_events import _select_subset


def _station(name: str, lat: float) -> MiniStation:
    return MiniStation(
        name=name,
        network="XX",
        location="--",
        channel="BHZ",
        latitude=lat,
        longitude=0.0,
    )


def test_select_subset_returns_all_when_pool_not_larger_than_count() -> None:
    stations = [_station(f"S{i}", float(i)) for i in range(5)]
    result = _select_subset(stations, 10)
    assert sorted(s.name for s in result) == sorted(s.name for s in stations)


def test_select_subset_spreads_across_full_latitude_range() -> None:
    """The original bug this guards against: a floor-divided stride
    degrades to no spread at all once the pool isn't much bigger than the
    target count.
    """
    stations = [_station(f"S{i}", float(i - 50)) for i in range(100)]
    subset = _select_subset(stations, 10)

    assert len(subset) == 10
    lats = [s.latitude for s in subset]
    assert lats == sorted(lats)
    assert lats[0] == -50.0
    assert lats[-1] == 49.0
    assert max(b - a for a, b in zip(lats, lats[1:])) > 1


def test_select_subset_is_deterministic() -> None:
    stations = [_station(f"S{i}", float(i)) for i in range(37)]
    first = [s.name for s in _select_subset(stations, 15)]
    second = [s.name for s in _select_subset(stations, 15)]
    assert first == second
