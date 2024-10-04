#!/bin/bash
#SBATCH --job-name=PET_global
#SBATCH --cpus-per-task=4
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1
#SBATCH --mail-type=all
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
  
if [ ! -f "pixi.toml" ]; then
    cp /p/moonshot2-casestudy/Wflow/africa/pixi.toml .
    pixi install
fi

echo "current working directory: $PWD"
echo "calculating PET"

pixi run snakemake -s "/p/moonshot2-casestudy/Wflow/africa/scripts/calc_pet.smk" --cores 4 --forceall --unlock
pixi run snakemake -s "/p/moonshot2-casestudy/Wflow/africa/scripts/calc_pet.smk" --cores 4 --forceall -n --quiet rules
pixi run snakemake -s "/p/moonshot2-casestudy/Wflow/africa/scripts/calc_pet.smk" --cores 4 --forceall --keep-incomplete --keep-going