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

workdir: f"{syscheck()}/moonshot2-casestudy/Wflow/africa"

print(f"{'*'*10}\n::CWD:: {os.getcwd()}\n{'*'*10}")

base_dir = Path(os.getcwd()).as_posix()

def get_clusters(forecast):
    jobs = Jobs(Path(base_dir,forecast).as_posix())  # Create a Jobs object
    jobs.prepare()  # This calls locate_clusters internally
    return forecast, jobs.cluster_ids

CONFIGS = glob(str(Path(base_dir,"forecasts", "2024-10-16_MS2_workshop").as_posix())+"/*.yml")
print("CONFIGS:")
for config in CONFIGS:
    print(config)

model_location = "data/3-input/wflow_build/{}/staticmaps.nc"
CLUSTERS_n = {fc: get_clusters(fc)[1] for fc in CONFIGS}
print("CLUSTER NUMBERS:")
for fc, clusters in CLUSTERS_n.items():
    print(f"{fc}: {clusters}")

CLUSTERS_m = {fc: [Path(base_dir, model_location.format(c)).as_posix() for c in clusters] for fc, clusters in CLUSTERS_n.items()}
print("MODEL LOCATIONS:")
for fc, locations in CLUSTERS_m.items():
    print(f"{fc}:")
    for location in locations:
        print(f"  {location}")


"""
::: ALL :::
States that we want an output for each forecast, which is the postprocessed forecast.
We ALSO want an output per cluster, which is the forecast for that specific cluster.
"""
rule all:
    input: 
        # expand("{forecast}/wflow_output_{cluster}.nc", forecast=CONFIGS, cluster=get_clusters(forecast)),
        expand(Path("{forecast}","output.nc"), forecast=CONFIGS)

rule prepare_forecast:
    input:
        Path("{forecast}").as_posix()  # one input config
    output:
        Path("{forecast}","output.nc")
    script:
        Path("src","2-build","prepare.py").as_posix()  # make sure this script is Snakemake compatible

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

