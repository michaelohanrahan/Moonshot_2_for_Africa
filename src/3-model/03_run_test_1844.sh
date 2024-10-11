#!/bin/bash
#SBATCH --job-name=ms2_test_1844
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition 16vcpu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=sebastian.hartgring@deltares.nl

# cluster with ID 1844 corresponds to the 61th cluster in the list  
julia -t 4 03_run_cluster.jl 61