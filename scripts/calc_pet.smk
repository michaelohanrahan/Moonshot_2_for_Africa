from pathlib import Path
import os 


os.chdir(Path('/p/moonshot2-casestudy/Wflow/africa'))

years = range(1980, 2020)
out_dir = Path("data/4-output/PET_global")
method = 'penman-monteith_tdew'
freq = 'daily'
temp_precip_forcing = 'era5_daily_zarr'
# bbox = '[-22.148438,-36.315125,52.382813,38.685510]'
bbox = '[0,-66,359.75,66]'

rule all:
    input:
        expand(Path(out_dir, f'{temp_precip_forcing}_{method}_{freq}_'+'{year}'+'.nc'), year=years)

rule calc_pet:
    input: 
        dc = Path('data/1-external/deltares_data_linux.yml')
    params:
        tpf = temp_precip_forcing,
        method = method,
        tmin = '{year}',
        tmax = '{year}'
    output:
        Path(out_dir, f'{temp_precip_forcing}_{method}_{freq}_'+'{year}'+'.nc')
    threads: 4
    resources:
        mem_mb = 32000
    shell:
        """
        pixi run python scripts/00_PET.py --tpf {params.tpf} --method {params.method} --tmin {params.tmin} --tmax {params.tmax} --bbox {bbox}
        """

    