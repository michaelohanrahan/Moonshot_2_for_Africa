#= 
This script is used to run Wflow clusters in parallel on the cluster.
Example used from h7 wiki: https://publicwiki.deltares.nl/display/Deltareken/Example%3A+Using+Julia+in+parallel.
The batch script runs this script and provides an argument $SLURM_ARRAY_TASK_ID which corresponds to the i-th cluster that we want to run.
To run clusters i = 1 to N, we have to provide the keyword array=1-N+1.

Example code to run in the batch script for 4 clusters on the test nodes (2 threads per Wflow run):

#!/bin/bash
#SBATCH --job-name=ms2_cluster_$SLURM_ARRAY_TASK_ID
#SBATCH --cpus-per-task=2
#SBATCH --partition test
#SBATCH --array=1-5
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
  
julia -t $SLURM_CPUS_PER_TASK 03_run_cluster.jl $SLURM_ARRAY_TASK_ID
=#

# Paths for Windows and Linux
if Sys.iswindows()
    path_clusters = "p:/moonshot2-casestudy/Wflow/africa/src/3-model/wflow_build"
    path_run_wflow = "p:/moonshot2-casestudy/Wflow/africa/scripts/03_run_wflow_interp_080.jl"
elseif Sys.islinux()
    path_clusters = "/p/moonshot2-casestudy/Wflow/africa/src/3-model/wflow_build"
    path_run_wflow = "/p/moonshot2-casestudy/Wflow/africa/scripts/03_run_wflow_interp_080.jl"
end

# Get all clusters in the directory, filter out any non-directories.
# Convert to Int to sort by value, not lexicographically.
items = readdir(path_clusters)
folders = filter(item -> isdir(joinpath(path_clusters, item)), items)
clusters = sort(parse.(Int, folders))

# Find cluster i of 1 ... N to run
i = parse.(Int, only(ARGS))
cluster = clusters[i]
config = joinpath(path_clusters, string(cluster), "wflow_sbm.toml")
println("processing cluster $i with id $cluster and config $config")

# Run Wflow with on-the-fly forcing using script.
include(path_run_wflow)
run(config)