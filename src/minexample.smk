import os
from snakemake.utils import min_version

# Ensure minimum Snakemake version
min_version("6.0")

# Define samples and their specific conditions
SAMPLE_CONDITIONS = {
    "A": ["control", "treatment1", "treatment2"],
    "B": ["control", "treatment1"],
    "C": ["control", "treatment2", "treatment3"]
}

# Create output directory
os.makedirs("results", exist_ok=True)

def get_all_output_files():
    return [f"results/{sample}_{condition}_processed.txt" 
            for sample in SAMPLE_CONDITIONS 
            for condition in SAMPLE_CONDITIONS[sample]]

rule all:
    input:
        get_all_output_files(),
        "results/summary.txt"

rule process_sample:
    output:
        "results/{sample}_{condition}_processed.txt"
    params:
        threshold = 0.05
    run:
        if wildcards.condition in SAMPLE_CONDITIONS[wildcards.sample]:
            shell("echo 'Processing {wildcards.sample} under {wildcards.condition} condition' > {output}")
        else:
            raise ValueError(f"Invalid condition {wildcards.condition} for sample {wildcards.sample}")

def get_processed_files(wildcards):
    return [f"results/{sample}_{condition}_processed.txt" 
            for sample in SAMPLE_CONDITIONS 
            for condition in SAMPLE_CONDITIONS[sample]]

rule summarize:
    input:
        get_processed_files
    output:
        "results/summary.txt"
    run:
        with open(output[0], "w") as out:
            for file in input:
                out.write(f"Summarized: {file}\n")

rule clean:
    shell:
        "rm -rf results"