# Alaska teleseismic sample dataset

5 real teleseismic earthquakes, each recorded by a subset of Alaska's `AK`
(Alaska Earthquake Center), `AV` (Alaska Volcano Observatory, USGS/UAF-GI),
and `TA` (USArray Transportable Array, NSF-funded EarthScope/IRIS)
broadband networks — BHZ channel only. (The previous, undocumented event
this replaces is kept verbatim on the `aimbat_v1` branch.)

## Events

| Label | USGS eventid | Origin time (UTC) | Location | Depth | Mag | Stations |
|---|---|---|---|---|---|---|
| `komandorskiye_ostrova` | [us20009x42](https://earthquake.usgs.gov/earthquakes/eventpage/us20009x42) | 2017-07-17T23:34:13.740Z | 54.4434, 168.857 (Commander Islands, Russia) | 10 km | 7.7 | 36 |
| `kumamoto_japan` | [us20005iis](https://earthquake.usgs.gov/earthquakes/eventpage/us20005iis) | 2016-04-15T16:25:06.220Z | 32.7906, 130.7543 (Kumamoto, Japan) | 10 km | 7.0 | 36 |
| `pinotepa_mexico` | [us2000d3km](https://earthquake.usgs.gov/earthquakes/eventpage/us2000d3km) | 2018-02-16T23:39:39.280Z | 16.3855, -97.9787 (Pinotepa de Don Luis, Mexico) | 22 km | 7.2 | 35 |
| `fiji_region` | [usc000stdc](https://earthquake.usgs.gov/earthquakes/eventpage/usc000stdc) | 2014-11-01T18:57:22.380Z | -19.6903, -177.7587 (Fiji region) | 434 km | 7.1 | 34 |
| `papua_new_guinea` | [us10007uph](https://earthquake.usgs.gov/earthquakes/eventpage/us10007uph) | 2017-01-22T04:30:22.960Z | -6.2464, 155.1718 (Papua New Guinea region) | 135 km | 7.9 | 28 |

Each is fetched fresh from the USGS event API by `eventid` at run time (see
`fetch_events.py`). Directory names follow `Event_YYYY.MM.DD.HH.MM.SS.SSS/`,
derived from each event's resolved origin time.

## Station selection

For each event, `fetch_events.py`:

1. Queries `fdsnws-station` (`net=AK,AV,TA`, `cha=BHZ`, `level=channel`,
   window = origin time ± 1 day, `minlatitude=50` to exclude `TA`'s earlier
   Lower-48 deployment) via raw FDSN text — `pysmo.tools.web` has no
   station-discovery function.
2. Dedupes to unique `net.sta`.
3. Drops candidates at ≥95 degrees epicentral distance (`MAX_P_DISTANCE_DEG`)
   — beyond that lies the direct-P shadow zone, and every station in an
   event must be anchored on the same phase for cross-correlation to make
   sense.
4. Selects ~40 stations evenly spread by latitude (`numpy.linspace` over
   the sorted candidates, not a floor-divided stride, which degrades to no
   spread at all once the pool isn't much bigger than the target count).
5. Fetches each selected station, skipping (and logging) any with no
   waveform data for the window, no usable location code, or no
   deconvolvable response. Station counts above are what survived this —
   no event was padded to hit exactly 40.

## Window formula

For each surviving station: `dist_deg = haversine(event, station)`,
`travel_times = fetch_travel_times(depth_km, dist_deg, ["P"])`,
`predicted_p = origin_time + travel_times["P"]`.

- **Start** = `predicted_p − 2min`
- **End** = `predicted_p + 3min`

A P-coda-focused cut, not a full-event recording: 2 minutes ahead of P for
pick/QC headroom, 3 minutes after for the coda, well short of S and
surface waves. Dataset stays small (~9.4 MB total across 5 events /
169 stations).

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

`pyproject.toml` depends on `pysmo` via a local editable path (`../pysmo`),
since `fetch_sac`/`fetch_travel_times`/`SAC.fetch` aren't yet in a tagged
release. Switch to a pinned `git+https://github.com/pysmo/pysmo@<sha>` if
this needs to run outside that layout.

## Licence / attribution

- **`AK`** — Alaska Earthquake Center, University of Alaska Fairbanks.
- **`AV`** — Alaska Volcano Observatory, a joint USGS/University of
  Alaska Fairbanks Geophysical Institute program.
- **`TA`** — USArray Transportable Array, part of the NSF-funded
  EarthScope/IRIS effort.

Waveform and station metadata served via `service.earthscope.org`
(formerly IRIS DMC). Event parameters are public USGS FDSN
earthquake-catalogue information.
