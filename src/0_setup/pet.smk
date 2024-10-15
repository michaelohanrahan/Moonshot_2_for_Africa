from pathlib import Path
import numpy as np

out_dir = Path("/p/moonshot2-casestudy/Wflow/africa/data/3-input/global_era5_pet/")

years = np.arange(2000, 2025)
months = np.arange(1, 13)
months = [f"{month:02d}" for month in months]

tpf = "era5_hourly"
method = "debruin"
freq_data = "H"
freq_save = "M"
freq_out = "hourly"

#example p:\moonshot2-casestudy\Wflow\africa\data\3-input\global_era5_pet\era5_hourly_debruin_PET_hourly_2024_01_2024-01.nc

def define_pet_file(year, month, filetype):
    if not isinstance(year, list):
        year = [int(year)]
    pet_files = []
    for year in years:
        for month in months:
            pet_file = Path(out_dir, f"{tpf}_{method}_PET_{freq_out}_{year}_{month:02d}_{year}-{month:02d}.{filetype}")
            pet_files.append(pet_file)
    return pet_files

rule all:
    input:
        expand(Path(out_dir, f"{tpf}_{method}_PET_{freq_out}_"+"{year}_{month}_{year}-{month}.nc"), year=years, month=months),


rule calc_pet:
    params:
        tpf = tpf,
        method = method,
        freq_data = freq_data,
        freq_save = freq_save,
        year = "{year}",
        month = "{month}"
    output:
        output_file = Path(out_dir, f"{tpf}_{method}_PET_{freq_out}_"+"{year}_{month}_{year}-{month}.nc")
    group: "calc_pet"
    resources:
        mem_mb = 8000,
        time = "0-02:00:00",
        threads = 1
    shell:
        "pixi run python src/0_setup/00_PET.py \
        --tpf {params.tpf}\
        --method {params.method}\
        --freq-data {params.freq_data}\
        --freq-save {params.freq_save}\
        --year {params.year}\
        --month {params.month}"


            
