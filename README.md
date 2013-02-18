aig2shp.py: ArcInfo Grid ASCII to ESRI Shapefile conversion
=====================
Vectorization of ArcInfo Grid ASCII files as ESRI Shapefiles. 
Python program to convert an ArcInfo Grid ASCII raster file to an ESRI shapefile with 
a single layer containing oriented grid square polygons centered at the pixels 
of the grid, each containing an attribute value equal to the value of the corresponding 
ArcInfo raster pixel. The ArcInfo Grid ASCII header <code>cellsize</code> parameter
determines the size of the grid squares. 

In contrast to the [ArcGIS Raster to Polygon](http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//001200000008000000) and the [gdal_polygonize.py](http://www.gdal.org/gdal_polygonize.html) 
raster to polygon feature layer conversion utilities, <code>aig2shp.py</code> does not require 
float raster values to be discretized into (nonnegative) integer values before conversion.  Your 
raster data is preserved, which may have some advantages for geospatial analysis.  One disadvantage 
is that data isn't compressed, since a discrete classification isn't used.  Discretization of raster 
values can facilitate symbology assignment for cartographic display, and may be offered in future 
revisions.  Taking unions of adjacent grid squares with the same value isn't attempted. However, 
this leaves open the possibility of simultaneous vectorization of multiple raster files so that a 
vector grid square (at the highest resolution) is associated with multiple attributes, each from 
different rasters. (Not supported so far.)

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
usage: aig2shp.py [-h] [-a attribute] [-e minX minY maxX maxY] [-l LAYER] [-n]
                  [-v] [--version] [--wgs84]
                  grid_ASCII_file ESRI_shapefile

Create ESRI Shapefile grid poly coverage from ArcInfo Grid ASCII raster.

positional arguments:
  grid_ASCII_file       ArcInfo Grid ASCII input file.
  ESRI_shapefile        ESRI shapefile output file.

optional arguments:
  -h, --help            show this help message and exit
  -a attribute, --attr attribute
                        Name of attribute for ArcInfo grid values. Defaults to
                        "value."
  -e minX minY maxX maxY, --extent minX minY maxX maxY
                        Bounding box of subset of raster in geographic
                        coordinates.
  -l LAYER, --layer LAYER
                        Shapefile layer name string.
  -n, --nonzero         Exclude zero values.
  -v, --verbose         Display verbose output.
  --version             Show program version number and exit.
  --wgs84               Set spatial reference to WGS84/EPSG:4326 in shapefile
                        layer. Projection file (.prj) is written out.

Software is released under The MIT License (c) 2013 Florian Lengyel. All other
original work is licenced (CC BY-NC-SA 3.0 US) 2013 Florian Lengyel, CUNY
Environmental CrossRoads Initiative, Advanced Science Research Center, The
City College of New York. Contact: gmail/skype/twitter florianlengyel.
```

## Example ##
In this example, cropland for a region including Africa was subsetted from the Ramankutty
cropland raster data set [1] in ArcInfo Grid ASCII format, and used to produce a corresponding 
shapefile. 
``` 
./aig2shp.py -e -34.892837 -17.338675  37.428152 57.845763 --wgs84 -n -v ramankutty_cropland2000_frac_5m.asc vector_squares.shp
```
The shapefile was uploaded to CartoDB; a screenshot is shown below.
[<img src="https://raw.github.com/flengyel/ArcInfoGrid2Shapefile/master/AfricaCropland.png">](https://raw.github.com/flengyel/ArcInfoGrid2Shapefile/master/AfricaCropland.png)

## Dependencies ##
* [ProgressBar](http://code.google.com/p/python-progressbar/)
* GDAL, OGR python libraries
* numpy
* argparse

## Author ##
Florian Lengyel, [CUNY Environmental CrossRoads Initiative](http://asrc.cuny.edu/crossroads), 
[Advanced Science Research Center](http://asrc.cuny.edu/crossroads),
[The City College of New York](http://www.ccny.cuny.edu), [CUNY](http://www.cuny.edu).
Contact: gmail/skype/twitter florianlengyel 

## License ##
Original software is licensed under the [MIT License](http://opensource.org/licenses/MIT): (c) 2013 Florian Lengyel. All other original work is licensed under an [Attribution-NonCommercial-ShareAlike 3.0 United States](http://creativecommons.org/licenses/by-nc-sa/3.0/us/) 
Creative Commons License: (CC BY-NC-SA 3.0 US) 2013 Florian Lengyel. Derivative work is licensed accordingly.  Full license in License.txt.

## References ##
[1] [Chad Monfreda, Navin Ramankutty and Jonathan A. Foley. Farming the planet: 2. Geographic distribution of crop areas, yields, physiological types, and net primary production in the year 2000. Global Biogeochemical Cycles. Volume 22, Issue 1, March 2008.](http://dx.doi.org/10.1029/2007GB002947)
