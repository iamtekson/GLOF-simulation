# ArcGIS pro raster calculator

Change the environments as below:

* Processing extends to same as AOI
* Snap raster to same `alos_palsar_aoi_remove_acnopy_hpp.tif`
* Save output as `aoi_lake_removed.tif`

```shell
# remove lake from main aoi dem
Con(IsNull("lake.tif"),  "alos_palsar_aoi_remove_canopy_hpp.tif",  "alos_palsar_aoi_remove_canopy_hpp.tif" -  "lake.tif")

# remove ls14 and save it as elevation14.tif
Con(IsNull("ls14.tif"),  "aoi_lake_removed.tif",  "aoi_lake_removed.tif" -  "ls14.tif")
```
