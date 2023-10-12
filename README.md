# GLOF-simulation
The repository contains the GLOF methodology, code, literature review, and steps to produce a result.

![r.avaflow simulation (landslide-lake interaction)](https://github.com/iamtekson/GLOF-simulation/blob/main/RESULTS/r.avaflow_results/13_results/13_plots/13_hflow_map.gif "r.avaflow result")

## What is included in this repo?

**1. SOFTWARE:** The software folder includes r.avaflow version 2.4. The software was originally available on [r.avaflow official website](https://www.landslidemodels.org/r.avaflow/software.php).

**2. INSTALLATION GUIDELINE:** The r.avaflow installation guideline is explained in [`installation.md` file](https://github.com/iamtekson/GLOF-simulation/blob/main/Installation.md)

**3. SIMULATION CODE:** The code to reproduce the results are stored inside the `SIMULATION CODE` folder.

**4. INPUT DATA:** All the required input data to run the simulation are stored in the `INPUT DATA` folder. Inside this folder, there are the following folders:

  * Bathymetric data (`.csv` file)
  * kml (landslide scenarios shapes drawn in Google Earth Pro)
  * rasters (all the required input raster maps, i.e., landslide rasters, lake rasters, and elevation rasters)
  * shp (vector representation of the landslide and lake boundary)

**5. RESULTS:** All the results were generated from r.avaflow software and notebooks to generate the hydrographs for a report.
  * notebooks (`python` code for generating the graphs for the report)
  * r.avaflow_results (all the results generated from r.avaflow software It contains `.asc` files, `.txt` files, and some `.gif` files for animation of the landslide-lake interaction.

## Installation guide

For the GLOF simulation, we have used [Mergili et. al., 2014](https://gmd.copernicus.org/articles/8/4027/2015/) r.avaflow software. To install the software properly, please follow the [installation instruction](./Installation.md).

## Data Preparation

To prepare the data, open your grass GIS by typing `grass` on the terminal and follow the following steps:

1. Create a new location on the project directory with a `generic cartesian coordinate system`
2. Click on `PERMANENT` and `switch mapset`
3. Follow the following steps in the terminal:

```sh
# Change the directory
cd path/to/input/data/directory

# Change data files to grass format
r.in.gdal -o input=elevation.tif output=elevation

# Check the region 
g.region -p

# Change the default region
g.region -s rast=elevation

# Change all other raster data to grass format using r.in.gdal command
r.in.gdal -o input=hrelease.tif output=hrelease
```

## Simulation using r.avaflow

To start the simulation with r.avaflow, run the following command in the terminal:

```sh
# Basic run
r.avaflow prefix=1 phases=m elevation=elevation hrelease1=hrelease1 hrelease3=hrelease3

# Run with some additional parameters
r.avaflow prefix=3 phases=m elevation=elevation hrelease1=hrelease1 hrelease3=hrelease3 friction=35,20,3,0,0,0,0.05 time=10,300
```

To know more about the r.avaflow input, check the [official documentation of r.avaflow](https://www.landslidemodels.org/r.avaflow/manual.php).
