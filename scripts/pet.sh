#!/bin/bash
#SBATCH --job-name=PET_MS2                                             # Job name
#SBATCH --output=/u/ohanrah/documents/Moonshot_2_for_Africa/data/0-log/h7/pet_%j.log # Standard output and error log
#SBATCH --time=2-00:00:01  # Job duration (hh:mm:ss)
#SBATCH --account=hyd
#SBATCH --partition 1vcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Set variables
# Go one directory up to set PWD
cd "/u/ohanrah/documents/Moonshot_2_for_Africa"

pixi run snakemake -c 1 -s "src/0_setup/pet.smk" --profile "./.config/forecast/" --unlock
echo "Unlocked directory"
pixi run snakemake -c 1 -s "src/0_setup/pet.smk" --profile "./.config/forecast/" --forceall -n
echo "Running snakemake with 1 core"
pixi run snakemake -c 1 -s "src/0_setup/pet.smk" --profile "./.config/forecast/" --forceall

