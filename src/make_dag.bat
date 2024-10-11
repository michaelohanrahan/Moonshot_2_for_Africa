@echo off

echo Activating the pixi environment
pixi run snakemake -s "snakefile" --dag | dot -Tsvg > dag.svg

pause