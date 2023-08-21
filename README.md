# MYRIAD – Hazard Event Sets Algorithm (MYRIAD-HESA): 
This repository provides a python script that can be used to compile historically-based multi-hazard event sets from single hazard events that have occured on varying time, space, and intensity scales. 
## Python requirements/dependencies

* [Python](https://www.python.org/) 3.11.4
* [NumPy](https://numpy.org/_) 1.25.2
* [Pandas](https://pandas.pydata.org/) 2.0.3
* [Shapely](https://shapely.readthedocs.io/en/stable/manual.html) 2.0.1
* [GeoPandas](https://geopandas.org/en/stable/index.html) 0.13.2

All packages with the correct dependencies can be extracted from the YML file as an environment using [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

```bash

# Add conda-forge channel for extra packages
conda config --add channels conda-forge

# Create a conda environment for the project and install packages
conda env create -f environment.yml
conda activate MYRIAD-HESA

```
## How to cite
If you have used MYRIAD-HESA, please cite our paper:



### License
Copyright (C) 2023 Judith Claassen. All versions released under the [MIT license](LICENSE).
