Moonshot 2: Wflow for Africa
==============================

Creating a workflow to cluster a multipolygon vector layer into manageably sized Wflow model domains. These will beshelved models stored for later use in investigating compound flooding. Forcing will be implemented 'on-the-fly' with downscaling upon initiation of a model run to save on storage cost.

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── bin                 <- Your compiled model code can be stored here (not tracked by git)
    ├── config              <- Configuration files, e.g., for doxygen or for your model if needed
    ├── data                
    │   ├── 1-external      <- Data external to the project.
    │   ├── 2-interim       <- Intermediate data that has been altered.
    │   ├── 3-input         <- The processed data sets, ready for modeling.
    │   ├── 4-output        <- Data dump from the model.
    │   └── 5-visualization <- Post-processed data, ready for visualisation.
    ├── docs                <- Documentation, e.g., doxygen or scientific papers (not tracked by git)
    ├── scripts             <- Main branch python files
    ├── reports             <- For a manuscript source, e.g., LaTeX, Markdown, etc., or any project reports
    │   └── figures         <- Figures for the manuscript or reports
    └── src                 <- Source code for this project
        ├── 0-setup         <- Install necessary software, dependencies, pull other git projects, etc.
        ├── 1-prepare       <- Scripts and programs to process data, from 1-external to 2-interim.
        ├── 2-build         <- Scripts to create model specific input from 2-interim to 3-input. 
        ├── 3-model         <- Scripts to run model and convert or compress model results, from 3-input to 4-output.
        ├── 4-analyze       <- Scripts to post-process model results, from 4-output to 5-visualization.
        └── 5-visualize     <- Scripts for visualisation of your results, from 5-visualization to ./report/figures.

## GIT: [https://github.com/michaelohanrahan/Moonshot_2_for_Africa.git]

## Scripts

### 01_cluster_basins.py

This script takes a multipolygon file of independent basins and performs clustering. The clustering is performed
in a manner that reduces overlap of the resulting model domains. 

#TODO: determine if the resulting domains, which are large, should be split. 