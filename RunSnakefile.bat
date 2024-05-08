@echo off

echo Activating the pixi environment
echo Unlocking the directory
pixi run snakemake --unlock -s "snakefile" --configfile "config/snakeConfig.yaml"

echo Executing Snakefile "snakefile" with config "config/snakeConfig.yaml"
pixi run snakemake all -s "snakefile" -c 4 --configfile "config/snakeConfig.yaml"  --rerun-incomplete

pause