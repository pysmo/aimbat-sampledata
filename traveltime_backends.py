"""Travel-time backends for the sample dataset's initial picks.

Each backend has the shape of `pysmo.tools.traveltime.TravelTimeBackend`.
`PHASE_BACKENDS` names the backend used for each phase. pysmo's built-in
solver covers P but not depth phases, so pP comes from ObsPy's TauP for
now. Moving pP to pysmo once its solver covers it means changing its entry
in `PHASE_BACKENDS`.

ObsPy is imported inside `obspy_taup`, so importing this module needs no
ObsPy. It is in the `fetch` dependency group only.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from pysmo.tools.traveltime import TravelTimeBackend, travel_times

MODEL = "iasp91"


def pysmo_builtin(
    *, depth: float, distance: float, phases: Sequence[str]
) -> Mapping[str, pd.Timedelta]:
    """Travel times from pysmo's built-in solver."""
    return travel_times(depth=depth, distance=distance, phases=phases, model=MODEL)  # type: ignore[arg-type]


def obspy_taup(
    *, depth: float, distance: float, phases: Sequence[str]
) -> Mapping[str, pd.Timedelta]:
    """First-arrival travel times from ObsPy's TauP.

    Needs the `fetch` dependency group.
    """
    from obspy.taup import TauPyModel

    arrivals = TauPyModel(model=MODEL).get_travel_times(
        source_depth_in_km=depth / 1000.0,
        distance_in_degree=distance,
        phase_list=list(phases),
    )
    times: dict[str, pd.Timedelta] = {}
    for arrival in arrivals:
        if arrival.name in phases and arrival.name not in times:
            times[arrival.name] = pd.Timedelta(seconds=arrival.time)
    return times


PHASE_BACKENDS: dict[str, TravelTimeBackend] = {
    "P": pysmo_builtin,
    "pP": obspy_taup,
}


def arrival_times(
    *, depth: float, distance: float, phases: Sequence[str]
) -> dict[str, pd.Timedelta]:
    """Travel time of each phase from the backend `PHASE_BACKENDS` names for it.

    Raises:
        KeyError: If a phase has no backend.
        ValueError: If a backend gives no arrival for its phase.
    """
    times: dict[str, pd.Timedelta] = {}
    for phase in phases:
        result = PHASE_BACKENDS[phase](depth=depth, distance=distance, phases=[phase])
        if phase not in result:
            raise ValueError(f"no {phase} arrival at {distance:.2f} deg")
        times[phase] = result[phase]
    return times
