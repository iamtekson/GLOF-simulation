# GLOF-simulation
The repository contain the GLOF methodology, code, literature review and steps to produce a result.

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
r.avaflow elevation=elevation hrelease1=hrelease1 hrelease3=hrelease2 phases=m prefix=3 friction=35,20,3,0,0,0,0 control=0,0,1,0,0,0 time=10,300
```

To know more about the r.avaflow input, check the [official documentation of r.avaflow](https://www.landslidemodels.org/r.avaflow/manual.php).