#!/bin/bash -l
#SBATCH --job-name=pantanal_plot_inmaps                                    # Job name
#SBATCH --output="/p/11206499-pantanal/Climate change/h7/run_plot_%j.log"      # Standard output and error log
#SBATCH --time=0:30:00                                                                 # Job duration (hh:mm:ss)
#SBATCH --partition test #16vcpu
#SBATCH --ntasks=4                                                                   # Number of tasks (analyses) to run
#SBATCH --mail-user=sebastian.hartgring@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

conda init
conda activate hydromt-wflow
python /c/Users/hartgrin/OneDrive - Stichting Deltares/Projecten/Moonshot/Moonshot_2_for_Africa/scripts/02_hydromt_wflow_build.py
