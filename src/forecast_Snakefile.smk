import sys
import os

# Add the parent directory of 'src' to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
print("Python path:")
for path in sys.path:
    print(path)
import importlib
create_forecast = importlib.import_module('2_build.2_create_forecast')
Jobs = create_forecast.Jobs
Config = create_forecast.Config
from icecream import ic
from scripts.helper import syscheck

workdir: f"{syscheck()}/moonshot2-casestudy/Wflow/africa"

print(f"{'*'*10}\n::CWD:: {os.getcwd()}\n{'*'*10}")


def get_clusters(forecast):
    jobs = Jobs(f"forecast_config/{forecast}.yml")  # Create a Jobs object
    jobs.prepare()  # This calls locate_clusters internally
    return forecast, jobs.cluster_ids

CONFIGS = glob.glob("forecasts/2024-10-16_MS2_workshop/*.yml")
model_location = "data/3-input/wflow_build/{}/staticmaps.nc"
CLUSTERS = {fc: get_clusters(fc) for fc in CONFIGS}


"""
::: ALL :::
States that we want an output for each forecast, which is the postprocessed forecast.
We ALSO want an output per cluster, which is the forecast for that specific cluster.
"""
rule all:
    input: 
        # expand("{forecast}/wflow_output_{cluster}.nc", forecast=CONFIGS, cluster=get_clusters(forecast)),
        expand("{forecast}/output.nc", forecast=CONFIGS)

rule prepare_forecast:
    input:
        "config/"+"{forecast}"+".yml"  # one input config
    output:
        expand("{forecast}/output.nc", forecast=CONFIGS)
    script:
        "src/2-build/prepare.py"  # make sure this script is Snakemake compatible

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

