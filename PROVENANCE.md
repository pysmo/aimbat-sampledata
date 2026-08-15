# Alaska teleseismic sample dataset

5 real teleseismic earthquakes, each recorded by a subset of Alaska's `AK`
(Alaska Earthquake Center), `AV` (Alaska Volcano Observatory, USGS/UAF-GI),
and `TA` (USArray Transportable Array, NSF-funded EarthScope/IRIS)
broadband networks — BHZ channel only. (The previous, undocumented event
this replaces is kept verbatim on the `aimbat_v1` branch.)

## Events

| Label | USGS eventid | Origin time (UTC) | Location | Depth | Mag | Stations |
|---|---|---|---|---|---|---|
| `komandorskiye_ostrova` | [us20009x42](https://earthquake.usgs.gov/earthquakes/eventpage/us20009x42) | 2017-07-17T23:34:13.740Z | 54.4434, 168.857 (Commander Islands, Russia) | 10 km | 7.7 | 29 |
| `fiji_region` | [usc000stdc](https://earthquake.usgs.gov/earthquakes/eventpage/usc000stdc) | 2014-11-01T18:57:22.380Z | -19.6903, -177.7587 (Fiji region) | 434 km | 7.1 | 39 |
| `nepal` | [us20002926](https://earthquake.usgs.gov/earthquakes/eventpage/us20002926) | 2015-04-25T06:11:25.950Z | 28.2305, 84.7314 (Nepal, Gorkha earthquake) | 8.2 km | 7.8 | 35 |
| `iraq` | [us2000bmcg](https://earthquake.usgs.gov/earthquakes/eventpage/us2000bmcg) | 2017-11-12T18:18:17.180Z | 34.9110, 45.9588 (Halabja, Iraq) | 19 km | 7.3 | 39 |
| `solomon_islands` | [usc000phx5](https://earthquake.usgs.gov/earthquakes/eventpage/usc000phx5) | 2014-04-12T20:14:39.300Z | -11.3487, 162.0025 (Kirakira, Solomon Islands) | 22.6 km | 7.6 | 25 |

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
  event's 29 stations, not a handful of outliers — likely a genuine
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

## Station selection

All 5 events share the same ~40-station set (rather than reselecting
independently per event), so most of those stations have multi-event
coverage — valuable for this test dataset, since ICCS/MCCC workflows
exercise a station across events. `fetch_events.py`:

1. Queries `fdsnws-station` (`net=AK,AV,TA`, `cha=BHZ`, `level=channel`,
   window = origin time ± 1 day, `minlatitude=50` to exclude `TA`'s earlier
   Lower-48 deployment) via raw FDSN text — `pysmo.tools.web` has no
   station-discovery function. This is repeated once per event, since
   operational stations vary over the 2014–2017 span of the set.
2. Dedupes each event's result to unique `net.sta`.
3. Intersects all 5 per-event station sets to the `net.sta` operating
   during *every one* of the events' windows, then drops any candidate at
   ≥95 degrees epicentral distance (`MAX_P_DISTANCE_DEG`) from *any* of
   those events — beyond that lies the direct-P shadow zone, and every
   station in an event must be anchored on the same phase for
   cross-correlation to make sense.
4. Selects ~40 stations evenly spread by latitude (`numpy.linspace` over
   the sorted candidates, not a floor-divided stride, which degrades to no
   spread at all once the pool isn't much bigger than the target count).
5. Fetches this same station list for each of the 5 events, skipping (and
   logging) any with no waveform data for that event's window, no usable
   location code, or no deconvolvable response.

Station counts in the table above are what survived per-event fetching —
no event was padded to hit exactly 40, and a station can fail for one
event while succeeding for another.

## Window formula

For each surviving station: `dist_deg = haversine(event, station)`,
`travel_times = fetch_travel_times(depth_km, dist_deg, ["P"])`,
`predicted_p = origin_time + travel_times["P"]`.

- **Start** = `predicted_p − 2min`
- **End** = `predicted_p + 3min`

A P-coda-focused cut, not a full-event recording: 2 minutes ahead of P for
pick/QC headroom, 3 minutes after for the coda, well short of S and
surface waves. Dataset stays small (~9.7 MB total across 5 events /
167 station-events, 40 unique stations — 17 of which appear in all 5
events, 14 more in 4 of 5, only 1 in as few as 2).

`_fetch_one_event` checks that each station's fetched data actually spans
the full requested window (a real FDSN data gap can otherwise return a
short partial trace that doesn't even reach the P arrival), skipping and
logging any that don't; `tests/test_dataset_integrity.py` also guards
against a truncated seismogram slipping into a commit.

## Instrument response removal

Raw dataselect output is uncalibrated digital counts. Each station's
seismogram is deconvolved (`_deconvolve` in `fetch_events.py`) following
`pysmo.tools.signal.remove_response`'s recipe:

1. Fetch the station's StationXML response (`fetch_stationxml`).
2. Detrend and taper (5%).
3. Deconvolve with frequency corners derived from the response itself
   (Nyquist- and pole-bounded), not hardcoded.

Output units follow `response.input_units` (`m/s` for every station here);
`idep` is set accordingly (`vel`) and recorded per-station in
`manifest.csv`.

## Metadata annotation

FDSN dataselect has no concept of an earthquake event, so fetched files
have no `evla`/`evlo`/`evdp`/`o` set. These are annotated afterwards via
pysmo's own `SacEvent` setters directly in Python (unlike
`reference_event`'s stricter SAC-binary annotation in the pysmo repo,
which exists there to avoid round-tripping through the code under test —
not a concern here). `iztype` is left at its default (`unkn`).

Each file also gets an initial phase-arrival pick: `sac.timestamps.t0` is
set to the same `predicted_p` used to anchor its fetch window.

Each event directory also has a `manifest.csv` (network, station,
location, latitude, longitude, distance_deg, azimuth_deg, samples, o,
idep).

## SAC version 7

Every file is written as v7 (`sac.native.nvhdr = 7`), which adds a
double-precision footer for time headers. Needed because `o`/`t0` are
derived from a real floating-point computation (`predicted_p`), and
v6's float32 precision would visibly truncate them.

## Regenerating

```sh
uv sync
uv run python fetch_events.py
```

`pyproject.toml` pins `pysmo` to a specific commit via `[tool.uv.sources]`
(`git = "https://github.com/pysmo/pysmo", rev = "<sha>"`), since
`fetch_sac`/`fetch_travel_times`/`SAC.fetch` aren't yet in a tagged
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
