@echo

echo Activating the pixi environment
pixi run snakemake --forceall --dag | dot -Tpdf > dag.pdf 

pause