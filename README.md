# GLOF-simulation
The repository contain the GLOF methodology, code, literature review and steps to produce a result.

## What is included in this repo?

**1. SOFTWARE:** The software folder includes the r.avaflow version 2.4. The software originally available in [r.avaflow official website](https://www.landslidemodels.org/r.avaflow/software.php).

**2. INSTALLATION GUIDLINE:** The r.avaflow installation guidline is explained in [`installation.md` file](https://github.com/iamtekson/GLOF-simulation/blob/main/Installation.md)

**3. SIMULATION CODE:** The code to reproduce the results are stored inside `SIMULATION CODE` folder.

**4. INPUT DATA:** All the required input data to run the simulation are stored inside `INPUT DATA` folder. Inside this folder there are following folders,

  * Bathymetric data (`.csv` file)
  * kml (landslide scenarios shapes drawn in google earth pro)
  * rasters (all the required input raster maps, i.e., landslide rasters, lake raster, elevation rasters)
  * shp (vector representation of the landslide and lake boundary)

**5. RESULTS:** All the results generated from r.avaflow software and notebooks to generate the hydrographs for report.
  * notebooks (`python` code for generating the graphs for report)
  * r.avaflow_results (all the result generated from r.avaflow software. It contain `.asc` files, `.txt` files and some `.gif` files for animation of the landslide-lake interaction)

## Installation guide

For the GLOF simulation, we have used [Mergili et. al., 2014](https://gmd.copernicus.org/articles/8/4027/2015/) r.avaflow software. To install the software properly, please follow the [installation instruction](./Installation.md).

## Data preparation

To prepare the data, open your grass GIS by typing `grass` on terminal and follow the following steps,

1. Create new location on the project directory with `generic cartesian coordinate system`
2. Click on `PERMENANT` and `switch mapset`
3. Follow the following steps in terminal,

```sh
# Change the directory
cd path/to/input/data/directory

# Change data files to grass format
r.in.gdal -o input=elevation.tif output=elevation

# Check the region 
g.region -p

# Change the default region
g.region -s rast=elevation

# Change all other rasater data to grass format using r.in.gdal command
r.in.gdal -o input=hrelease.tif output=hrelease
```

## Simulation using r.avaflow

To start the simulation with r.avaflow, run the following command in terminal,

```sh
# Basic run
r.avaflow prefix=1 phases=m elevation=elevation hrelease1=hrelease1 hrelease2=hrelease2 

# Run with some additional parameters
r.avaflow elevation=elevation hrelease1=hrelease1 hrelease3=hrelease2 phases=m prefix=3 friction=35,20,3,0,0,0,0.05 time=10,300
```

To know more about the r.avaflow input, check the [official documentation of r.avaflow](https://www.landslidemodels.org/r.avaflow/manual.php).
