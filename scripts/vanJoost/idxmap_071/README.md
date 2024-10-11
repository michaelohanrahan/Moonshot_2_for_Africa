# Run Wflow with coarser forcing (via index mapping)

Allow Wflow to run with forcing that is on a different resolution than the Wflow model. Forcing
can be on a irregular spaced grid (and potentially different projection, as long as lat-lon
values are provided for each forcing pixel). An index-mapping file is required, that maps each
Wflow model grid cell to an index in the forcing file. This file can be created using the
Python script. This assumes a nearest neighbor style matching.

A lapse rate can be supplied to correct temperature for the original forcing and Wflow model
resolutions. The elevation (orography) is expected to be in meters, and on the same grid
(extend and resolution) as the forcing. The wflow_dem is used as the Wflow model elevation.

## Scripts
- `run_wflow_idxmap_071.jl`: custom `Wflow.run()` function, that allows to run Wflow with a forcing
  file that is not matching the Wflow model grid. As additional input, an index file need to be
  presented (at the model resolution) that describes which forcing pixel index needs to be used
  for each Wflow model pixel. This results in a on-the-fly nearest neighbor interpolation.
- `wflow_settings_base.toml`: example TOML file that uses the new `[forcing]` section that is
  required by the custom run function.
- `create_idx_file.py`: Python script that creates the index mapping file at the Wflow model
  resolution
- `overlap_orog.py`: additional Python script to ensure that the forcing orography file matches
  the extent of the forcing file.


## Add to toml

```toml
[forcing]
lapse_correction = true
lapse_rate = -0.0065
path_orography = "path/to/elevation.nc"
layer_name = "orog"
path_idx = "path/to/forcing_idx.nc"
lon_idx_name = "lon_idx"
lat_idx_name = "lat_idx"
```

## How to run
Running the custom Wflow function is similar to running a normal Wflow simulation. A slightly
different syntax is required, that calls the julia script, with the path to the TOML file as
input. An example: `julia run_wflow_idxmap_071.jl "path/to/settings.toml"` (options as number of
threads can still be passed).