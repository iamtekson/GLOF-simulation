r.in.gdal -o --overwrite input=elevation1.tif output=elevation1
r.in.gdal -o --overwrite input=elevation2.tif output=elevation2
r.in.gdal -o --overwrite input=elevation3.tif output=elevation3
r.in.gdal -o --overwrite input=elevation4.tif output=elevation4
r.in.gdal -o --overwrite input=elevation5.tif output=elevation5
r.in.gdal -o --overwrite input=elevation6.tif output=elevation6
r.in.gdal -o --overwrite input=elevation7.tif output=elevation7


r.in.gdal -o --overwrite input=ls1_final.tif output=ls1
r.in.gdal -o --overwrite input=ls2_final.tif output=ls2
r.in.gdal -o --overwrite input=ls3_final.tif output=ls3
r.in.gdal -o --overwrite input=ls4_final.tif output=ls4
r.in.gdal -o --overwrite input=ls5_final.tif output=ls5
r.in.gdal -o --overwrite input=ls6_final.tif output=ls6
r.in.gdal -o --overwrite input=ls7_final.tif output=ls7

r.in.gdal -o --overwrite input=hrelease3.tif output=hrelease3

g.region -s rast=elevation1

g.region -p

r.avaflow prefix=1 phases=m elevation=elevation1 hrelease1=ls1 hrelease3=hrelease3 orthophoto=orthophoto.tif hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180
r.avaflow prefix=2 phases=m elevation=elevation2 hrelease1=ls2 hrelease3=hrelease3 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180
r.avaflow prefix=3 phases=m elevation=elevation3 hrelease1=ls3 hrelease3=hrelease3 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180
r.avaflow prefix=4 phases=m elevation=elevation4 hrelease1=ls4 hrelease3=hrelease3 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180
r.avaflow prefix=5 phases=m elevation=elevation5 hrelease1=ls5 hrelease3=hrelease3 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180

r.avaflow prefix=6 phases=m elevation=elevation6 hrelease1=ls6 hrelease3=hrelease3 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180
r.avaflow prefix=7 phases=m elevation=elevation7 hrelease1=ls7 hrelease3=hrelease3 time=10,1000 hydrocoords=500851,5581936,100,180,501253,5581693,100,180,501830,5581550,100,180

