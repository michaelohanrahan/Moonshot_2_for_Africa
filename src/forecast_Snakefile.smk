import sys
import os
from glob import glob
# Add the parent directory of 'src' to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Use importlib for the module with numbers in its name
import importlib
create_forecast = importlib.import_module('2_build.2_create_forecast')
Jobs = create_forecast.Jobs
Config = create_forecast.Config
from icecream import ic
# Assuming 'scripts' is a directory at the same level as 'src'
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from helper import syscheck
import json

workdir: f"{syscheck()}/moonshot2-casestudy/Wflow/africa"

print(f"{'*'*10}\n::CWD:: {os.getcwd()}\n{'*'*10}")
base_dir = Path(os.getcwd()).as_posix()
#DIRECTORIES
log_dir = Path(base_dir, "data","0-log")
inter_dir = Path(base_dir,"data","2-interim")
input_dir = Path(base_dir,"data","3-input")
output_dir = Path(base_dir,"data","4-output")

def get_clusters(forecast):
    """
    Get the clusters for a given forecast.
    input:
        forecast: The forecast config file.
    Returns:
        forecast: The forecast name.
        cluster_ids: The cluster IDs.
    """
    jobs = Jobs(Path(base_dir,forecast).as_posix())  # Create a Jobs object
    jobs.prepare()  # This calls locate_clusters internally
    return forecast, jobs.cluster_ids


def prepare(base_dir):
    """
    Preparing helper files for the forecast workflow using the create_forecast module.
    """
    
    #Listing the files with forecast parameters
    CONFIGS = glob(str(Path(base_dir,"forecasts", "2024-10-16_MS2_workshop").as_posix())+"/*.yml")
    print("\nCONFIGS:")
    for config in CONFIGS:
        print(config)

    FORECASTS = [Path(config).name.split(".")[0] for config in CONFIGS]
    print("\nFORECASTS:")
    for forecast in FORECASTS:
        print(forecast)

    CLUSTERS = {fc: get_clusters(cfg)[1] for fc,cfg in zip(FORECASTS, CONFIGS)}

    print("\nCLUSTERS:")
    for fc, clusters in CLUSTERS.items():
        print(f"{fc}: {clusters}")

    model_cfg = "data/3-input/wflow_build/{}/wflow_sbm.toml"
    
    FC_DICT = {fc: [Path(base_dir, model_cfg.format(c)).as_posix() for c in clusters] for fc, clusters in CLUSTERS.items()}
    print("\nFC_DICT:")
    for fc, locations in FC_DICT.items():
        print(f"{fc}:")
        for location in locations:
            print(f"  {location}")
    
    #output key value pairs to a json file
    with open(Path(output_dir, 'FC_DICT.json').as_posix(), 'w') as f:
        json.dump(FC_DICT, f, indent=4)

    return CONFIGS, FORECASTS, CLUSTERS, FC_DICT
#CONFIGS: list of config files
#FORECASTS: list of forecast names
#CLUSTERS: dictionary with forecast names as keys and cluster IDs as values
#FC_DICT: dictionary with forecast names as keys and model locations as values
CONFIGS, FORECASTS, CLUSTERS, FC_DICT = prepare(base_dir)

"""
::: ALL :::
States that we want an output for each forecast, which is the postprocessed forecast.
We ALSO want an output per cluster, which is the forecast for that specific cluster.
Define the output files and get the wildcards later.
"""

def get_all_output_files(filename):
    return [str(output_dir)+f"/{forecast}/{cluster}/"+filename 
            for forecast in FORECASTS for cluster in CLUSTERS[forecast]]

rule all:
    input: 
        get_all_output_files("output_scalar.nc")
         
""" 
::: PREPARE FORECAST :::
Running Wflow for each cluster, for each forecast. 
"""

rule prepare_forecast_instate_config:
    output:
        file = str(output_dir)+"/{forecast}/{cluster}/output_scalar.nc"
    run:
        import os
        import subprocess
        os.makedirs(os.path.dirname(output.file), exist_ok=True)
        with open(output.file, "w") as f:
            f.write("output_scalar.nc")

# rule run_forecast:
#     input:
#         lambda wildcards: expand("{forecast}/wflow_forecast_{cluster}.toml", forecast=wildcards.forecast, cluster=get_clusters(wildcards.forecast))
#     output:
#         "{forecast}/wflow_output_{cluster}.nc"
#     run:
#         if os.path.isfile(f"{wildcards.forecast}/wflow_state_{wildcards.cluster}.toml"):
#             shell("julia run_script.jl {wildcards.forecast}/wflow_state_{wildcards.cluster}.toml")
#         shell("julia run_script.jl {wildcards.forecast}/wflow_forecast_{wildcards.cluster}.toml")

# rule postprocess_forecast:
#     output: "{forecast}/output.nc"  # one output per forecast
#     input:
#         lambda wildcards: expand("{forecast}/wflow_output_{cluster}.nc", forecast=wildcards.forecast, cluster=get_clusters(wildcards.forecast))
#     script: "postprocess.py"  # make sure this script is Snakemake compatible

