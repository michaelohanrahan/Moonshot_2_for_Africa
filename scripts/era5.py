import xarray as xr
from icecream import ic
import traceback

import xarray as xr
import scipy.spatial
import numpy as np

def mirror_point_at_360(ds):
    extra_point = (
        ds.where(ds.longitude == 0, drop=True)
        .assign_coords(longitude=lambda x: x.longitude + 360)
    )
    return xr.concat([ds, extra_point], dim='values')

def build_triangulation(x, y):
    grid = np.stack([x, y], axis=1)
    return scipy.spatial.Delaunay(grid)

def interpolate(data, tri, mesh):
    indices = tri.find_simplex(mesh)
    ndim = tri.transform.shape[-1]
    T_inv = tri.transform[indices, :ndim, :]
    r = tri.transform[indices, ndim, :]
    c = np.einsum('...ij,...j', T_inv, mesh - r)
    c = np.concatenate([c, 1 - c.sum(axis=-1, keepdims=True)], axis=-1)
    result = np.einsum('...i,...i', data[:, tri.simplices[indices]], c)
    return np.where(indices == -1, np.nan, result)

def regrid(grid):
    longitude = np.linspace(0, 360, num=360*4+1)
    latitude = np.linspace(-90, 90, num=180*4+1)
    mesh = np.stack(np.meshgrid(longitude, latitude, indexing='ij'), axis=-1)
    
    

def main():
    # Open the dataset
    # ds = xr.open_zarr(
    #     'gs://gcp-public-data-arco-era5/co/single-level-reanalysis.zarr',
    #     chunks={'time':48},
    #     consolidated=True,
    #     storage_options=dict(token='anon'),
    # )
    
    ds = xr.open_zarr(
    "gs://gcp-public-data-arco-era5/ar/1959-2022-full_37-1h-0p25deg-chunk-1.zarr-v2",
    chunks={'time': 48},
    consolidated=True,
    )
    
    ic(ds)
    print(f"ATTRS::\n{ds.attrs}")
    print(f"TEMPORAL BOUNDS::\n{ds.time.values.min()} - {ds.time.values.max()}")
    for var in ds.data_vars:
        print(f"VAR:: {var}")
if __name__ == "__main__":
    main()