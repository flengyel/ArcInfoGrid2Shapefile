aig2shp.py: ArcInfo Grid ASCII to ESRI Shapefile conversion
=====================
Vectorization of ArcInfo Grid ASCII files as ESRI Shapefiles. 
Python program to convert an ArcInfo Grid ASCII raster file to an ESRI shapefile with 
a single layer containing oriented grid square polygons centered at the pixels 
of the grid, each containing an attribute value equal to the value of the corresponding 
ArcInfo raster pixel. The ArcInfo Grid ASCII header <code>cellsize</code> parameter
determines the size of the grid squares. 

The program was written to upload raster data in a format useable by 
[CartoDB](http://www.cartodb.com). See the [following correspondence](https://groups.google.com/d/msg/cartodb/fbjRhgO-AMo/x8Mfy_Z_8DgJ) on the [CartoDB google group](https://groups.google.com/forum/?fromgroups=#!forum/cartodb).

## Notes on coordinates ##
The following was determined from ESRI documentation of ArcInfo Grid ASCII files online, and 
was verified by examination of the output of the <code>gdalinfo</code> GDAL utility applied to 
ArcInfo Grid ASCII files, and by examination of converted files within GIS systems.
* The ArcInfo Grid ASCII raster format stores values in Lon, Lat <code>(x, y)</code> order, but the ESRI shapefile format stores coordinates in Lat, Lon <code>(y, x)</code> order. 
* The origin of the raster is defined to be the <code>(x, y)</code> coordinates of the upper-left corner of the upper left grid square of size <code>cellsize</code>, namely <code>(x0, y0) = (xllcorner, yllcorner + nrows * cellsize)</code>. 
* The coordinates of the center of the grid square at the origin are then <code>(x0 + cellsize/2, y0 - cellsize/2)</code>. 
* In general, the geographic coordinates <code>(x, y)</code> of the pixel at 
<code>(row, col)</code> are given by <code>(x, y) = (xllcorner + (col + 1/2) * cellsize, yllcorner + (nrows - row - 1/2) * cellsize)</code>.

    
## Usage ##
```
usage: aig2shp.py [-h] [-a attribute] [-e minX minY maxX maxY] [-l LAYER]
                     [--version] [--wgs84] [-x EXCLUDE]
                     grid_ASCII_file ESRI_shapefile

Create ESRI Shapefile grid poly coverage from Arcinfo Grid ASCII raster.

   positional arguments:
     grid_ASCII_file       ArcInfo Grid ASCII input file.
     ESRI_shapefile        ESRI shapefile output file.

   optional arguments:
     -h, --help            show this help message and exit
     -a attribute, --attr attribute
                           Name of attribute for ArcInfo grid values. Defaults to
                           "value."
     -e minX minY maxX maxY, --extent minX minY maxX maxY
                           Bounding box in geographic coordinates.
     -l LAYER, --layer LAYER
                           Shapefile layer name string.
     --version             show program's version number and exit
     --wgs84               Set spatial reference to WGS84/EPSG:4326 in shapefile
                           layer. Projection file (.prj) is written out.
     -x EXCLUDE, --exclude EXCLUDE
                           Exclude zero.

   (CC BY-NC-SA 3.0 US) 2013 Florian Lengyel, CUNY Environmental CrossRoads
   Initiative, Advanced Science Research Center, The City College of New York.
   Contact: gmail/skype/twitter florianlengyel.
```

## Author ##
Florian Lengyel, [CUNY Environmental CrossRoads Initiative](http://asrc.cuny.edu/crossroads), 
[Advanced Science Research Center](http://asrc.cuny.edu/crossroads),
[The City College of New York](http://www.ccny.cuny.edu), [CUNY](http://www.cuny.edu).
Contact: gmail/skype/twitter florianlengyel 

## License ##
[Attribution-NonCommercial-ShareAlike 3.0 United States](http://creativecommons.org/licenses/by-nc-sa/3.0/us/)
Full license in License.txt
