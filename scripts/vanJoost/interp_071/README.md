# Run Wflow with coarser forcing (regular grid)

Allow Wflow to run with forcing that is on a different resolution than the Wflow model. Forcing
is required to be in the same projection as the Wflow model, and on a regular grid. The lat-lon
values are used to interpolate the forcing values to the model resolution, using nearest
neighbor interpolation.

A lapse rate can be supplied to correct temperature for the original forcing and Wflow model
resolutions. The elevation (orography) is expected to be in meters, and on the same grid
(extend and resolution) as the forcing. The wflow_dem is used as the Wflow model elevation.

## New section in the TOML

```toml
[forcing]
lapse_correction = true
lapse_rate = -0.0065
path_orography = "path/to/elevation.nc"
layer_name = "orog"
```