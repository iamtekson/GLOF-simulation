# ArcGIS pro raster calculator

Change the environments as below:

* Processing extends to same as AOI
* Snap raster to same `alos_palsar_aoi_remove_acnopy_hpp.tif`
* Save output as `aoi_lake_removed.tif`

```shell
# remove lake from main AOI dem
Con(IsNull("lake.tif"),  "alos_palsar_aoi_remove_canopy_hpp.tif",  "alos_palsar_aoi_remove_canopy_hpp.tif" -  "lake.tif")

# remove ls14 and save it as elevation14.tif
Con(IsNull("ls14.tif"),  "aoi_lake_removed.tif",  "aoi_lake_removed.tif" -  "ls14.tif")
```
## Code simulation-related learning

* Please don't use the `layers=1` in simulation; it will end up with the `Numeric Failure` error
* Use `Control > Entrainment control = 1` rather than 2, If you put 2, the result will be strange as below image

![Simplified Pudasaini and Krautblatter (2021) entrainment and deposition model](https://github.com/iamtekson/GLOF-simulation/assets/39838116/c2d232ac-b42e-495e-8fd4-84d3074a86c1)

