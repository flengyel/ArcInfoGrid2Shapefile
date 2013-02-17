ArcInfoGrid2Shapefile
=====================
Vectorization of Arc Info Grid ASCII files as ESRI Shapefiles. 

## Author ##
Florian Lengyel, Advanced Science Research Center, CUNY

### Usage ###
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
     --   version             show     ram's version number and exit
     --wgs84               Set spatial reference to WGS84/EPSG:4326 in shapefile
                           layer. Projection file (.prj) is written out.
     -x EXCLUDE, --exclude EXCLUDE
                           Exclude zero.

   (CC BY-NC-SA 3.0 US) 2013 Florian Lengyel, CUNY Environmental CrossRoads
   Initiative, Advanced Science Research Center, The City College of New York.
   Contact: gmail/twitter florianlengyel.

