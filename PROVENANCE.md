# Alaska teleseismic sample dataset

6 real teleseismic earthquakes, each recorded by a subset of Alaska's `AK`
(Alaska Earthquake Center), `AV` (Alaska Volcano Observatory, USGS/UAF-GI),
and `TA` (USArray Transportable Array, NSF-funded EarthScope/IRIS)
broadband networks. One event has BHZ only, the other five BHZ, BHN and
BHE, all resampled to 20 Hz. The mix is deliberate: real projects are not
uniform. (The previous, undocumented event this replaces is kept verbatim
on the `aimbat_v1` branch.)

## Events

| Label | USGS eventid | Origin time (UTC) | Location | Depth | Mag | Channels | Stations |
|---|---|---|---|---|---|---|---|
| `komandorskiye_ostrova` | [us20009x42](https://earthquake.usgs.gov/earthquakes/eventpage/us20009x42) | 2017-07-17T23:34:13.740Z | 54.4434, 168.857 (Commander Islands, Russia) | 10 km | 7.7 | Z | 29 |
| `fiji_region` | [usc000stdc](https://earthquake.usgs.gov/earthquakes/eventpage/usc000stdc) | 2014-11-01T18:57:22.380Z | -19.6903, -177.7587 (Fiji region) | 434 km | 7.1 | Z, N, E | 38 |
| `nepal` | [us20002926](https://earthquake.usgs.gov/earthquakes/eventpage/us20002926) | 2015-04-25T06:11:25.950Z | 28.2305, 84.7314 (Nepal, Gorkha earthquake) | 8.2 km | 7.8 | Z, N, E | 29 |
| `iraq` | [us2000bmcg](https://earthquake.usgs.gov/earthquakes/eventpage/us2000bmcg) | 2017-11-12T18:18:17.180Z | 34.9110, 45.9588 (Halabja, Iraq) | 19 km | 7.3 | Z, N, E | 36 |
| `solomon_islands` | [usc000phx5](https://earthquake.usgs.gov/earthquakes/eventpage/usc000phx5) | 2014-04-12T20:14:39.300Z | -11.3487, 162.0025 (Kirakira, Solomon Islands) | 22.6 km | 7.6 | Z, N, E | 25 |
| `central_peru` | [usb000s5x1](https://earthquake.usgs.gov/earthquakes/eventpage/usb000s5x1) | 2014-08-24T23:21:45.520Z | -14.598, -73.5714 (central Peru) | 101 km | 6.8 | Z, N, E | 34 |

Each is fetched fresh from the USGS event API by `eventid` at run time (see
`fetch_events.py`). Directory names follow `Event_YYYY.MM.DD.HH.MM.SS.SSS/`,
derived from each event's resolved origin time.

## Known data characteristics

Real data occasionally comes with real complications. Documented here so
they don't get mistaken for bugs, and deliberately left in — a dataset of
only clean best-case seismograms would be a worse fixture for exercising
AIMBAT's manual QC step than one with a few genuinely awkward examples.

- **`nepal`: local contamination at two Augustine Volcano stations.**
  `AV.AU22` and `AV.AUJA` (co-located, ~3 km apart) both show a sharp,
  impulsive local-event onset ~80s before the predicted Nepal P arrival,
  decaying over the next 30-60s — not present at any other station in the
  shared set. Augustine is a densely monitored, seismically active
  volcano (AVO network since 1993); this is very likely an unrelated
  local VT earthquake or similar that happened to land in the fetch
  window, not a processing error. Kept as-is: a good real-world case for
  exercising ICCS/MCCC's outlier rejection.

- **`komandorskiye_ostrova`: emergent (slow-onset) P arrival, elevated
  cycle-skip risk.** Median rise time to half-peak amplitude is ~16s
  (~10 wave cycles at this event's ~1.6s dominant period) versus under 1
  cycle for `fiji_region`, and this is consistent across most of the
  event's stations, not a handful of outliers — likely a genuine
  property of this M7.7 event's rupture (complex/slow nucleation), not a
  processing artifact; overall SNR is actually among the highest of the 5
  events, so this isn't an amplitude problem. In ICCS, start with a wider
  `window_pre`/`window_post` and a lower `bandpass_fmin` (both
  runtime-adjustable on the `ICCS` instance) to get a robust coarse
  alignment before narrowing either for a final refinement pass — a
  narrow window from the start is the most likely way to get a cycle
  skip on this event. After `run_mccc()`, check for high `errors`
  combined with high `cc_means` and low `cc_stds` — `mccc()`'s own
  docstring names that combination as the cycle-skip signature. Kept
  as-is deliberately: a good real-world case for testing that workflow.

- **`fiji_region`: pP at `AV.AKRB` cycle-skips in MCCC.** `AV.AKRB` has
  the weakest P of the event, and its pP aligns in ICCS but MCCC moves it
  by about 6 s, to a stack CC near zero. The pick is kept: a real case for
  the manual QC step between ICCS and MCCC, here on the second phase.

- **`nepal`: S energy arrives about 40 s after the predicted S.** At every
  station, the largest transverse S pulse comes about 40 s after the iasp91
  prediction, consistent with the long rupture of the Gorkha earthquake,
  which released most of its moment well after the origin time. The S
  picks are the predictions, so alignment needs a window that reaches the
  later pulse.

- **Orientation that is not 0/90, and that changes between events.**
  `AK.CCB`'s BHN and BHE point 6 deg and 96 deg at every event. `AK.GRNC`'s
  point 16 deg and 106 deg at `solomon_islands` (April 2014) and 0 deg and
  90 deg from `central_peru` (August 2014) onwards: one stream whose
  orientation depends on the event.

## Station selection

All events share the same ~40-station set (rather than reselecting
independently per event), so most of those stations have multi-event
coverage — valuable for this test dataset, since ICCS/MCCC workflows
exercise a station across events. The set is pinned in `fetch_events.py`
(`PINNED_SITES`): the station service's listings change over time, so
re-running the selection below gives a different set (half the sites
changed between two runs in November 2026), and these notes, testkit's
fixtures and the checks behind the picks all refer to this one. A fetch
looks the pinned sites up at the five P events and fails if one is no
longer listed. `central_peru`, beyond the P shadow, uses the same sites.

The set was chosen once by `_select_common_stations`, which still
documents the method:

1. Queries `fdsnws-station` (`net=AK,AV,TA`, `cha=BHZ`, `level=channel`,
   window = origin time ± 1 day, `minlatitude=50` to exclude `TA`'s earlier
   Lower-48 deployment) via raw FDSN text — `pysmo.tools.web` has no
   station-discovery function. This is repeated once per event, since
   operational stations vary over the 2014–2017 span of the set.
2. Dedupes each event's result to unique `net.sta`.
3. Intersects the 5 P events' station sets to the `net.sta` operating
   during *every one* of their windows, then drops any candidate at
   ≥95 degrees epicentral distance (`MAX_P_DISTANCE_DEG`) from *any* of
   those events — beyond that lies the direct-P shadow zone, and every
   station in an event must be anchored on the same phase for
   cross-correlation to make sense.
4. Selects ~40 stations evenly spread by latitude (`numpy.linspace` over
   the sorted candidates, not a floor-divided stride, which degrades to no
   spread at all once the pool isn't much bigger than the target count).
A fetch then:

5. Fetches the pinned station list for each of the 6 events, every channel
   of the event, skipping (and logging) a station if any of its channels
   has no waveform data for the window, a gap in it (dataselect returns
   several segments), no usable location code, or no deconvolvable
   response. A station is therefore kept with all of the event's channels
   or not at all.

Station counts in the table above are what survived per-event fetching —
no event was padded to hit exactly 40, and a station can fail for one
event while succeeding for another.

## Window formula

For each surviving station, `predicted_arrivals` in `fetch_events.py`
computes the predicted arrival of each of the event's phases (iasp91),
locally rather than from a web service.

- **Start** = first predicted arrival − 2 min
- **End** = last predicted arrival + 3 min

All channels of a station share the window. For `komandorskiye_ostrova`
that is a P-coda cut, as before; for the four events with S it runs from
before P to the S coda (about 15 min); for `central_peru` it surrounds
SKS and the S that follows it. The dataset is about 32 MB: 515 files at
40 unique stations, 13 of which appear in all 6 events, 13 more in 5, and
only 1 in as few as 2.

`_fetch_one_event` checks that each channel's fetched data actually spans
the full requested window (a real FDSN data gap can otherwise return a
short partial trace that doesn't even reach the first arrival), skipping
and logging any station that doesn't; `tests/test_dataset_integrity.py`
also guards against a truncated seismogram slipping into a commit.

## Instrument response removal

Raw dataselect output is uncalibrated digital counts. Each station's
seismogram is deconvolved (`_deconvolve` in `fetch_events.py`) following
`pysmo.tools.signal.remove_response`'s recipe:

1. Fetch the station's StationXML response (`fetch_stationxml`).
2. Detrend and taper (5%).
3. Deconvolve with frequency corners derived from the response itself
   (Nyquist- and pole-bounded), not hardcoded.

Output units follow `response.input_units` (`m/s` for every station here);
`idep` is set accordingly (`vel`) and recorded per file in `manifest.csv`.

## Resampling

Every trace is then low-passed at 8 Hz (zero-phase, 4 corners) and
resampled to 20 Hz (`DELTA`), down from 40 Hz (`TA`) or 50 Hz (`AK`,
`AV`). 20 Hz is ample for P aligned at up to 2 Hz and for S and SKS below
0.2 Hz, and it keeps the dataset small. pysmo's `resample` uses the
Fourier method, which cuts the spectrum hard; the low-pass keeps it from
ringing.

## Orientation

Each channel's azimuth and dip come from its StationXML epoch at the time
of the event (the one fetched for the response) and are written through
pysmo's `SAC.orientation` (`cmpaz`, `cmpinc`). Every channel here has
them; a channel whose StationXML had none would be written without them.

## Metadata annotation

FDSN dataselect has no concept of an earthquake event, so fetched files
have no `evla`/`evlo`/`evdp`/`o` set. These are annotated afterwards via
pysmo's own `SacEvent` setters directly in Python (unlike
`reference_event`'s stricter SAC-binary annotation in the pysmo repo,
which exists there to avoid round-tripping through the code under test —
not a concern here). `iztype` is left at its default (`unkn`).

Each file also gets labelled initial picks, the predicted arrivals of its
event's phases, on the channels each phase is picked on:

| Phase | Header | Label | Channels | Events |
|---|---|---|---|---|
| P | `t0` | `kt0 = "P"` | BHZ | all but `central_peru` |
| pP | `t5` | `kt5 = "pP"` | BHZ | `fiji_region` |
| S | `t6` | `kt6 = "S"` | BHN, BHE | `fiji_region`, `nepal`, `iraq`, `solomon_islands` |
| SKS | `t7` | `kt7 = "SKS"` | BHN, BHE | `central_peru` |

S and SKS are measured on the transverse and radial components, which
need rotation. Their picks are on both horizontals as the station's
arrival times, for the rotated components to take over. In AIMBAT before
rotation, `aimbat data add --phase P --phase pP` stores them without
promoting them. ScS, PP, Pdif and PKIKP were checked on record sections
and left out: none was clear enough at these distances.

P and S come from pysmo's tau-p solver (`travel_times`). pP and SKS come
from ObsPy's TauP, since pysmo's solver does not cover depth or core
phases. `traveltime_backends.py` names the travel-time backend for each
phase. At 434 km depth, Fiji's pP arrives 94 to 99 s after P. `t5` to
`t7` are clear of the `t1` to `t3` that AIMBAT v1 writes.
`tests/test_dataset_integrity.py` checks the labels and channels, and that
each pick still agrees with its backend within 100 ms (the pP and SKS
checks only with ObsPy installed).

Each event directory also has a `manifest.csv` (network, station,
location, channel, latitude, longitude, distance_deg, azimuth_deg,
samples, o, idep), one row per file.

## SAC version 7

Every file is written as v7 (`sac.write(..., nvhdr=7)`), which adds a
double-precision footer for time headers. Needed because `o` and the
picks are derived from a real floating-point computation, and v6's
float32 precision would visibly truncate them.

## Regenerating

```sh
uv sync --group fetch
uv run --group fetch python fetch_events.py
```

The `fetch` dependency group brings ObsPy, needed for pP and SKS.

`pyproject.toml` pins `pysmo` to a specific commit via `[tool.uv.sources]`
(`git = "https://github.com/pysmo/pysmo", rev = "<sha>"`), since
`fetch_sac`/`SAC.fetch`/`travel_times` aren't yet in a tagged
release. When iterating on a local `pysmo` checkout, override transiently
with `{ path = "../pysmo", editable = true }` — swap back to the git+rev
pin before committing, since the path form has no meaning in CI (a fresh
checkout has no sibling `pysmo` directory).

## Licence / attribution

- **`AK`** — Alaska Earthquake Center, University of Alaska Fairbanks.
- **`AV`** — Alaska Volcano Observatory, a joint USGS/University of
  Alaska Fairbanks Geophysical Institute program.
- **`TA`** — USArray Transportable Array, part of the NSF-funded
  EarthScope/IRIS effort.

Waveform and station metadata served via `service.earthscope.org`
(formerly IRIS DMC). Event parameters are public USGS FDSN
earthquake-catalogue information.
