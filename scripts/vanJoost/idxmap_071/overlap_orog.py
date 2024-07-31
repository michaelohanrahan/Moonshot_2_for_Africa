import xarray as xr


def create_overlap_orog(forcing_fn, orography_fn, outfile_fn):
    # Read forcing dataset
    forc = xr.open_dataset(forcing_fn, chunks={"time": 10})

    # Read orography dataset
    orog = xr.open_dataset(orography_fn)

    # Find overlapping region (assuming orography is bigger than the forcing file)
    selection = orog.sel(rlat=forc.rlat, rlon=forc.rlon)
    selection.to_netcdf(outfile_fn)
    return selection


if __name__ == "__main__":
    forcing_fn = "./data/RACMO23_RESAMPLES_t2m_resample_future_Pd_definitive_biascorrected_RACMO.nc"
    orography_fn = "./data/orog_WEU-11_EC-Earth-Consortium-EC-Earth3_historical_r0i0p5f0_KNMI-RACMO23E_v2ds_fx.nc"
    outfile_fn = "./data/RACMO23_orog.nc"

    selection = create_overlap_orog(
        forcing_fn=forcing_fn, orography_fn=orography_fn, outfile_fn=outfile_fn
    )
