@echo off

pixi shell
snakemake --dag | dot -Tsvg > dag.svg

pause