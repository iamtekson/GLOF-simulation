# Grass GIS

![grass data structure](https://github.com/iamtekson/GLOF-simulation/assets/39838116/b5d0bfa3-dccb-430e-9d25-796571ecc606)

* g. is the general command
* r. is raster command
* v. is the vector command

```shell

# list of raster and vector datasets
g.list rast
g.list vect

# raster, vector information
r.info raster_layer
v.info vector_layer

# create new location
g.proj -c location=miller_lake -t datumtrans=1

# change the new location
g.mapset mapset=PERMANENT location=miller_lake

# Check the region 
g.region -p

# Change the default region
g.region -s rast=elevation
```

