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
import yaml

workdir: f"{syscheck()}/moonshot2-casestudy/Wflow/africa"

print(f"{'*'*10}\n::CWD:: {os.getcwd()}\n{'*'*10}")
base_dir = Path(os.getcwd()).as_posix()
#DIRECTORIES
log_dir = Path(base_dir, "data","0-log")
inter_dir = Path(base_dir,"data","2-interim")
input_dir = Path(base_dir,"data","3-input")
output_dir = Path(base_dir,"data","4-output")

def get_clusters_and_state_files(forecast):
    """
    Get the clusters and state file names for a given forecast.
    """
    jobs = Jobs(Path(base_dir,forecast).as_posix())
    jobs.prepare()
    clusters = jobs.cluster_ids
    state_files = {cluster: jobs.get_state_file_name(forecast, cluster) for cluster in clusters}
    return forecast, clusters, state_files

def prepare(base_dir):
    """
    Preparing helper files for the forecast workflow using the create_forecast module.
    """
    
    #Listing the files with forecast parameters
    CONFIGS = glob(str(Path(base_dir,"forecasts", "2024-10-16_MS2_workshop").as_posix())+"/*.yml")
    print("\nCONFIGS:")
    for config in CONFIGS:
        print(config)

    FORECASTS = [yaml.safe_load(open(config))["name"] for config in CONFIGS]
    print("\nFORECASTS:")
    for forecast in FORECASTS:
        print(forecast)

    CLUSTERS = {fc: get_clusters_and_state_files(cfg)[1:][0] for fc,cfg in zip(FORECASTS, CONFIGS)}
    print("\nCLUSTERS:")
    for fc, clusters in CLUSTERS.items():
        print(f"{fc}: {clusters}")

    
    
    STATE_FILES = {fc: get_clusters_and_state_files(cfg)[2] for fc,cfg in zip(FORECASTS, CONFIGS)}
    print("\nSTATE_FILES:")
    for fc, state_files in STATE_FILES.items():
        print(f"{fc}: {state_files}")

    model_cfg = "data/3-input/wflow_forecast/{}/wflow_sbm.toml"
    
    FC_DICT = {fc: [Path(base_dir, model_cfg.format(c)).as_posix() for c in clusters] for fc, clusters in CLUSTERS.items()}
    print("\nFC_DICT:")
    for fc, locations in FC_DICT.items():
        print(f"{fc}:")
        for location in locations:
            print(f"  {location}")
    
    #output key value pairs to a json file
    with open(Path(output_dir, 'FC_DICT.json').as_posix(), 'w') as f:
        json.dump(FC_DICT, f, indent=4)

    return CONFIGS, FORECASTS, CLUSTERS, STATE_FILES, FC_DICT
#CONFIGS: list of config files
#FORECASTS: list of forecast names
#CLUSTERS: dictionary with forecast names as keys and cluster IDs as values
#FC_DICT: dictionary with forecast names as keys and model locations as values
CONFIGS, FORECASTS, CLUSTERS, STATE_FILES, FC_DICT = prepare(base_dir)


"""
::: ALL :::
States that we want an output for each forecast, which is the postprocessed forecast.
We ALSO want an output per cluster, which is the forecast for that specific cluster.
Define the output files and get the wildcards later.
"""

state_dir = "data/3-input/wflow_state/{}/{}" #will be formatted as forcing, date

def get_warmup_toml(FORECASTS, CLUSTERS):
    return [str(output_dir)+f"/{forecast}/{cluster}/warmup.toml" 
            for forecast in FORECASTS 
            for cluster in CLUSTERS[forecast]]


# def get_instate_files(FORECASTS, CLUSTERS, STATE_FILES):
#     return [str(output_dir)+f"/wflow_state/{cluster}/{file}",
#             for forecast in FORECASTS 
#             for cluster in CLUSTERS[forecast] 
#             for file in STATE_FILES[forecast][cluster]]

def get_all_output_files(filename, FORECASTS, CLUSTERS):
    return [str(output_dir)+f"/{forecast}/{cluster}/"+filename 
            for forecast in FORECASTS 
            for cluster in CLUSTERS[forecast]]

rule all:
    input: 
        get_warmup_toml(FORECASTS, CLUSTERS),
        get_all_output_files("output_scalar.nc", FORECASTS, CLUSTERS),
         
""" 
::: PREPARE FORECAST :::
Running Wflow for each cluster, for each forecast. 
"""

rule prepare_forecast_instate_config:
    output:
        file = [Path(output_dir, "{forecast}", "{cluster}", "warmup.toml").as_posix() 
                for forecast in FORECASTS 
                for cluster in CLUSTERS[forecast]]
    params:
        config = Path(CONFIGS["{forecast}"]).as_posix()
    localrule: True
    script:
        "2_build/2_create_forecast.py {params.config}"

"""
::: RUN INSTATE :::
Running the instate run for each forecast, cluster.
#TODO finish the create forecast  logic
"""

rule run_forecast:
    input:
        toml = str(output_dir)+"/{forecast}/{cluster}/wflow_sbm_forecast.toml"
    params: 
        warmup = str(output_dir)+"/{forecast}/{cluster}/wflow_sbm_warmup.toml"
    output:
        file = str(output_dir)+"/{forecast}/{cluster}/output_scalar.nc"
    run:
    
        if os.path.isfile(params.warmup):
            shell(
                """
                julia --project='{params.project}' run_script.jl {params.warmup}
                julia --project='{params.project}' run_script.jl {input.toml}
                """
            )
        else:
            shell(
                """
                julia --project='{params.project}' run_script.jl {input.toml}
                """
            )
