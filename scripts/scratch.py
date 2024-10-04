from hydromt import DataCatalog

d = DataCatalog(r"p:\wflow_global\hydromt\deltares_data.yml")

d.get_rasterdataset("era5_daily_zarr")
