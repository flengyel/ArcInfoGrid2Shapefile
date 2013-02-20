#!/usr/bin/python
# aig2shp.py
# Author: Florian Lengyel
# Date: February 17, 2013
# Contact: gmail/skype/twitter florianlengyel
# License: MIT License (c) 2013 Florian Lengyel
 
from __future__ import division
try:
  from osgeo import ogr
  from osgeo import osr
except:
  import ogr
  import osr
import os, sys
import numpy as np
import argparse as arg
from progressbar import ProgressBar, Percentage, Bar

# Define the arguments first
descript = "Create ESRI Shapefile grid poly coverage from ArcInfo Grid ASCII raster."
epistr   = "Software is released under The MIT License (c) 2013 Florian Lengyel, CUNY Environmental CrossRoads Initiative, Advanced Science Research Center, The City College of New York. Contact: gmail/skype/twitter florianlengyel."
parser  = arg.ArgumentParser( description = descript, epilog=epistr )
# value of dest derived from first long opt
parser.add_argument('-a', '--attr',
                    metavar='attribute',
		    default='value',
		    help='Name of attribute for ArcInfo grid values. Defaults to "value."')
parser.add_argument('-e', '--extent', 
		    nargs=4, 
		    type=float,
		    metavar=('minX', 'minY', 'maxX', 'maxY'),
		    help='Bounding box of subset of raster in geographic coordinates.')
parser.add_argument('-l', '--layer', 
                    default='grid_value',
                    help='Shapefile layer name string.')
parser.add_argument('-m', '--multiplier',
		    type=int,
		    help='Multiply attribute column by the multiplier and take integer part. Useful in conjunction with QGIS dissolve.')
parser.add_argument('-n', '--nonzero', 
		    action="store_true",
		    help='Exclude zero values.')
parser.add_argument('-q', '--quiet', 
		    action="store_true",
		    help='Suppress progress bar.')
parser.add_argument('-v', '--verbose', 
		    action="store_true",
		    help='Display verbose output.')
parser.add_argument('--version', 
		    action='version', 
		    version='%(prog)s 0.3',
		    help='Show program version number and exit.')
parser.add_argument('--wgs84', 
		    action="store_true", 
		    help='Set spatial reference to WGS84/EPSG:4326 in shapefile layer. Projection file (.prj) is written out.')
parser.add_argument('infile', 
		    metavar='grid_ASCII_file',
		    help='ArcInfo Grid ASCII input file.')
parser.add_argument('outfile',
		    metavar='ESRI_shapefile',
		    help='ESRI shapefile output file.')



class ArcInfoGridASCII(object):
  """Header corresponding to Arc Info Grid ASCII file"""

  def getField(self, fieldname, typeconv):
    line = self.file.readline()
    field = line.split()
    if len(field) != 2:
      raise ValueError, "Expected field/value pair."
    name = field[0]
    if name  != fieldname:
      raise ValueError, "fieldname {0} missing.".format(fieldname)
    try:
      value = typeconv(field[1])
    except:
      raise ValueError, "Integer conversion of {0} failed.".format(field[1])
    if self.args.verbose:
      print '{0}: {1}'.format(fieldname, field[1])
    return value

  def __init__(self, args):
    self.args = args
    try:
      self.file = open(args.infile, "r")
    except IOError as e:
      raise IOError, "I/O error({0}): {1}".format(e.errno, e.strerror)
    self.ncols  = self.getField("ncols", int)
    self.nrows  = self.getField("nrows", int)
    self.xll    = self.getField("xllcorner", float)
    self.yll    = self.getField("yllcorner", float)
    self.cell   = self.getField("cellsize", float)
    self.nodata = self.getField("NODATA_value", float)
    # ESRI documentation states that the origin of the grid is
    # the upper left hand corner. From gdalinfo of the ascii grid file
    # of our test, it is evident that the coordinates of the  
    # center of the grid square whose upper left vertex is at the origin
    # are (xOrigin, yOrigin) as computed below.
    self.mid    = self.cell / 2
    self.xOrigin = self.xll + self.mid 
    self.yOrigin = (self.yll + self.nrows * self.cell) - self.mid

    # Finally, for the bounding box consistency check, compute the coordinates 
    # of the upper right corner of the upper right grid square
    self.xur    = self.xll + self.ncols * self.cell
    self.yur    = self.yll + self.nrows * self.cell

  def cart2geo(self, row, col):  
    """Convert Cartesian row and column to geographic coordinates of upper left
       corner of the cell. Returns (x, y)."""
    # ESRI documentation states that the origin of the grid is at upper left
    # and the terminus at the lower right. The coordinates of the upper left
    # are then (self.xll, -self.yll)
    x = self.xOrigin + col * self.cell
    y = self.yOrigin - row * self.cell
    return (x, y)

  def createGridSquare(self, lon, lat): 
    """Return wkb polygon with coordinates of grid square. Associate attribute"""
    ring = ogr.Geometry(ogr.wkbLinearRing)  # geometry for grid square
    minX = lon - self.mid
    maxX = lon + self.mid
    minY = lat - self.mid
    maxY = lat + self.mid

    #you may need to project coordinates
    #if args.wgs84:
    # we might want to project to WGS84
    #wgs84 = proj.Proj("+init=EPSG:4326") #LatLon with WGS84 datum used by Google
    #This seems unnecessary -- although this is probably because the grid files
    #I have been using were already in EPSG:4326.
    
    # ESRI shapefiles store points in lat/lon order, whereas ArcInfo Grid ASCII
    # files store points in lon/lat order.
    ring.AddPoint(minY, minX)
    ring.AddPoint(maxY, minX)
    ring.AddPoint(maxY, maxX)
    ring.AddPoint(minY, maxX)
    ring.AddPoint(minY, minX) # close the ring
    poly = ogr.Geometry(ogr.wkbPolygon) # create a new polygon
    poly.AddGeometry(ring) # add ring to polygon
    # add the attribute to the polygon
    return poly
    
class ExtentHandler(object):
  """Handle bounding boxes within ArcInfo Grid ASCII extents. Sets comparison
     function depending on arguments supplied to command line"""
  def __True__(self, lat, lon):  # use this comparison when -e, --extent is absent
     return True

  def __cmpFun__(self, lon, lat):
    return (self.minX <= lon and lon <= self.maxX and 
	    self.minY <= lat and lat <= self.maxY) 	  

  def __init__(self, hdr, args):
    self.hdf = hdr

    if args.extent == None:
      self.cmpFun = self.__True__
    else:  
      # set the extent
      self.minX, self.minY, self.maxX, self.maxY = args.extent
      # verify that the extent is within the bounds    
      valid = (hdr.xll <= self.minX and hdr.yll <= self.minY and
              self.maxX <= hdr.xur and self.maxY <= hdr.yur)
      if not valid:
        extError = 'Extent out of bounds: {0} {1} < {2} {3}'
        if self.minX < hdr.xll:
          raise ValueError, extError.format('minX', self.minX, hdr.xll, 'xllcorner')
        if self.minY < hdr.yll:
          raise ValueError, extError.format('minY', self.minY, hdr.yll, 'yllcorner')
        if hdr.xur < self.maxX:
          raise ValueError, extError.format('xurcorner', hdr.xul, self.maxX, 'maxX')
        if hdr.yur < self.maxY:
          raise ValueError, extError.format('yurcorner', hdr.xur, self.maxY, 'maxY')
      # verify that the extent defines a box 
      consistent = self.minX <= self.maxX and self.minY <= self.maxY
      if not consistent:
        incError = 'Extent inconsistent: {0} {1} < {2} {3}'
        if self.maxX < self.minX:
          raise ValueError, incError.format('maxX', self.maxX, self.minX, 'minX')
        if self.maxY < self.minY:
          raise ValueError, incError.format('maxY', self.maxY, self.minY, 'minY')
      # our exacting standards have been met
      self.cmpFun = self.__cmpFun__

  def compare(self, lat, lon):
    return self.cmpFun(lat, lon)

args = parser.parse_args()  # parse command line arguments

if args.verbose:
  print "Reading header..."

hdr = ArcInfoGridASCII(args)
ext = ExtentHandler(hdr, args)

if args.verbose:
  print "Reading array..."

grid1D = np.fromfile(hdr.file, sep = " \n")

# verify that the array can be reshaped
items = grid1D.shape[0]
if hdr.ncols * hdr.nrows != items:
  errorStr = "Number of items read in is {0} instead of {1}={2}*{3}"
  raise IOError, errorStr.format(items, hdr.nrows*hdr.ncols,
		                  hdr.nrows, hdr.ncols)

# reshape the array
if args.verbose:
  print "Reshaping array to grid..."
grid = np.reshape( grid1D, (hdr.nrows, hdr.ncols) )

# create the shapefile
driverName = "ESRI Shapefile"
drv = ogr.GetDriverByName( driverName )
if drv is None:
  raise ValueError, "{0} driver not available.".format(driverName)

shpFile = args.outfile
if os.path.exists( shpFile ):
  drv.DeleteDataSource( shpFile )

ds = drv.CreateDataSource( shpFile )
if ds is None:
  raise IOError, "Creation of output file {0} failed.".format(shpFile)

spatialReference = None
if args.wgs84:
  spatialReference = osr.SpatialReference()
  spatialReference.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')
layer = ds.CreateLayer( args.layer , spatialReference, ogr.wkbPolygon )
if layer is None:
  raise ValueError, "Layer creation failed."

# define the attribute at the centroid of the grid square
# a future version should set the type of the field from the command line
fieldef = ogr.FieldDefn( args.attr , ogr.OFTReal )
if layer.CreateField ( fieldef ) != 0:
  raise ValueError, "OGR field definition failed."

if args.verbose:
  print 'Converting to shapefile...'

if not args.quiet: # show the progress bar unless instructed otherwise
  pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval = hdr.nrows).start()

for row  in range(0, hdr.nrows):
  for col in range(0, hdr.ncols):
    v = grid[row][col]
    if v != hdr.nodata and (not args.nonzero or v != 0):
      lat, lon = hdr.cart2geo(row, col)
      if ext.compare(lon, lat):
        poly = hdr.createGridSquare(lon, lat)
        feature = ogr.Feature( layer.GetLayerDefn() )
	if args.multiplier != None:
	  v = int(v * args.multiplier)
        feature.SetField(args.attr, v)
        feature.SetGeometry(poly) # set the attribute
        if layer.CreateFeature(feature):
          raise ValueError, "Could not create feature in shapefile."
        feature.Destroy()
  if not args.quiet:
    pbar.update(row+1)

if not args.quiet:
  pbar.finish()

ds.Destroy()

