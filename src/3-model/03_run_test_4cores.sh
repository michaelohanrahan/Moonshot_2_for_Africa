#!/bin/bash
#SBATCH --job-name=ms2_test_4cores
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=4
#SBATCH --partition 4vcpu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sebastian.hartgring@deltares.nl

# trying 1st cluster
julia -t 4 03_run_cluster.jl 1