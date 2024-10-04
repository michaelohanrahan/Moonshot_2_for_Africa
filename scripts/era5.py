import xarray as xr
import icecream as ic

def main():
    # Open the dataset
    ds = xr.open_zarr(
        'gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3',
        chunks=None,
        storage_options=dict(token='anon'),
    )
    ic(ds)
    # Select only the variables you need
    variables = [
        '2m_temperature',
        '2m_dewpoint_temperature',
        '10m_u_component_of_wind',
        '10m_v_component_of_wind',
        'surface_net_solar_radiation'
    ]

    # Select the variables and the full time range
    subset = ds[variables].isel(time=0)

    # If you want to save this subset locally
    # subset.to_netcdf('era5_subset.nc')

    # If you want to work with the data directly, you can do so here
    # For example, to print information about the subset:
    print(subset)

if __name__ == "__main__":
    main()