#!/usr/bin/python
# -*- coding: utf-8 -*-
# DO NOT REMOVE THE SECOND LINE!
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
parser.add_argument('-c', '--coords',
		    nargs=2,
		    type=int,
		    metavar=('row', 'col'),
		    help='Debugging input. Raster coordinates of starting pixel.')
parser.add_argument('-d', '--dissolve',
		    action='store_true',
		    help='Dissolve Arc Info ASCII Grid in (row, col) space before converting to shapefile.')
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
parser.add_argument('-v', '--verbosity', 
		    action="count", default=0,
		    help="Display verbose message output. Each additional 'v' increases message verbosity: -vv is very verbose, and -vvv is very very verbose.")
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
    if self.args.verbosity >= 1:
      print '{0}: {1}'.format(fieldname, field[1])
    return value

  def __init__(self, args):
    if args.verbosity >= 1:
      print "Reading header..."
    try:
      self.file = open(args.infile, "r")
    except IOError as e:
      raise IOError, "I/O error({0}): {1}".format(e.errno, e.strerror)
    self.args   = args  # may not need to save this!
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

  def geo2cart(self, x, y):  # occasionally not correct -- don't convert back!
    """Inverse transform of cart2geo."""
    col = int( ( x - self.xOrigin ) / self.cell )
    row = int( ( self.yOrigin - y ) / self.cell )
    return (row, col)

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


# a StackExchange creation
def enum(**enums):
  return type('Enum', (), enums)

# vertex directions
vx = enum(z=0, r=1, d=2, l=3, u=4)


class Dissolver(object):
  """Dissolve polygons in raster space. Uses box coordinate representation.

     The box coordinates of the pixel at (i, j) meaning row i, column j,
     are denoted [r, c] and are related by r,c = 2*i+1, 2*j+1. The
     meaning of upper left UL, upper right UR, lower left LL, lower right
     LR, top edge, bot edge, left edge and right edge are define once
     and for all in the following diagram. The definitions are intended
     to comport with the orientation of raster space.

      UL→   Top   UR↓       [r-1,c-1]  [r-1,c] [r-1,c+1]
      Left  [r,c] Right     [r,c-1]    [r,c]   [r,c+1]
      LL↑   Bot   LR←       [r+1,c-1]  [r+1,c] [r+1,c+1]

      Note: Left, Top, Right and Bot aren't used. A little wasteful.
  """

  def __init__(self, args, hdr, ext, raster):
    """Create a Dissolver, using 
       args -- parsed arguments 
       hdr  -- Grid ASCII header 
       ext  -- extent of grid in coordinate space
       raster -- reshaped grid"""

    self.args = args
    self.hdr  = hdr
    self.ext  = ext
    self.raster = raster  # this is a raster object

    # This is the polygon region number. 
    self.region = 0  # The box at [r,c] is marked by the region number    

    # Create the box coordinate space (i, j) -> [2*i+1, 2*j+1]
    # Square brackets denote box coordinates, parentheses denote pixel coordinates.
    self.boxRows  = 2*hdr.nrows + 1
    self.boxCols  = 2*hdr.ncols + 1
    self.box   = np.zeros(shape = (self.boxRows, self.boxCols), dtype = np.int8) 
    if self.args.coords == None:
      self.i =  0  # CAUTION: the nextValid() function starts from (i,j+1)
      self.j = -1  # and must be called to set the current valid pixel
    else:
      self.i, self.j = self.args.coords
      self.j -= 1  # back up since nextValid() starts from (i, j+1)

    # define maps to add to [r,c] vertex coordinates
    self.goR    = { vx.r   :  0,  vx.d :  2, vx.l   :  0,  vx.u : -2}
    self.goC    = { vx.r   :  2,  vx.d :  0, vx.l   : -2,  vx.u :  0}

    # clockwise and counter-clockwise direction changes
    self.cw     = { vx.r : vx.d, vx.d : vx.l, vx.l : vx.u, vx.u : vx.r }
    self.ccw    = { vx.r : vx.u, vx.u : vx.l, vx.l : vx.d, vx.d : vx.r }

    # relative coordinates of inside and outside squares, 
    # given a vertex [r, c] and a direction vector vx.
    self.insideR = { vx.r : 1, vx.d : 1, vx.l :-1, vx.u :-1 }
    self.insideC = { vx.r : 1, vx.d :-1, vx.l :-1, vx.u : 1 }

    self.outsideR = { vx.r :-1, vx.d : 1, vx.l : 1, vx.u :-1 }
    self.outsideC = { vx.r : 1, vx.d : 1, vx.l :-1, vx.u :-1 }

    # used for diagnostics
    self.diag = { vx.r : 'r', vx.d : 'd', vx.l : 'l', vx.u : 'u' }

  def move(self, r, c, v):
    """recompute box coordinates using the mov enum"""
    if self.args.verbosity >= 7:
      r1 = r + self.goR[v]
      c1 = c + self.goC[v]
      edr = (r + r1) >> 1
      edc = (c + c1) >> 1
      self.box[edr][edc] -= 1
      if self.box[edr][edc] < -2:
        raise ValueError, 'Edge [{0},{1}] visited more than twice.'.format(edr,edc)
    return (r + self.goR[v], c + self.goC[v])

  def isBoxCoord(self, r, c):
    """True iff [r,c] is a valid box coordinate"""
    return ((0 <= r) and (r <  self.boxRows) and 
		  (0 <= c) and c <  (self.boxCols))


  def pixel2box(self, i, j):
    """Box coordinates [r, c] of pixel at raster coordinates (i, j)"""
    return ((i<<1)+1, (j<<1)+1)

  def box2pixel(self, r, c):
    """Return pixel coordinates (i,j) of box at box coordinates [r,c]"""
    return ((r-1)>>1, (c-1)>>1)
  
  def nextValid(self):
    """Find the next potential starting point [r,c] for an outer loop or an inner 
       annulus ring. Note that some squares may be "deep inside" another 
       region, but cannot be the starting point of any ring. The traverse() 
       algorithm should recognize these and mark them with the appropriate
       region number--though it could skip over them without marking."""

    def isValid(self, i, j): # this local function isn't used elsewhere
      """True iff the pixel at (i, j) is not nodata and the box at [r, c] has
      no region number assigned. Note the side effect--hence this is local."""
      self.r, self.c = self.pixel2box(i, j)
      #print i, j, self.raster[i][j]
      pixel = self.raster.grid[i][j]
      return ( pixel != self.hdr.nodata and 
               (not self.args.nonzero or pixel != 0) and 
	       self.box[self.r][self.c] == 0 )

    # box coordinates [r, c] of the next valid pixel
    i = self.i  # exhaust the current row
    # (i,j) was valid; start at (i, j+1) 
    for j in range(self.j+1, hdr.ncols): # (i,j+1) to (i, ncols-1)
      if isValid(self, i, j): # note that isValid is local and needs self
	 self.i, self.j = i, j
	 return True
    # continue with remaining rows
    for i in range(self.i+1, hdr.nrows):
      for j in range(0, hdr.ncols):
        if isValid(self, i, j):
	 self.i, self.j = i, j
	 return True
    return False

  def markBox(self, r, c, v, region):
    """Mark the inside  box at vertex [r, c] in direction v. If annulus
       is true, you should not mark any boxes.
       CAUTION: not marking the first box visited is an error!
    """
    r += self.insideR[v]
    c += self.insideC[v]
    if self.args >= 4:
      region0 = self.box[r][c]	    
      print '[{0},{1}] in region {2} marking with {3}'.format(r, c, region0, region)
    self.box[r][c] = region

  def isInOut(self, isIn, r, c, v, cls):
    """Returns True if the square clockwise (ccw if isIn is False) 
       from v at [r, c] exists and is in the class cls; false otherwise."""
    if isIn:   
      r += self.insideR[v]
      c += self.insideC[v]
    else:
      r += self.outsideR[v]
      c += self.outsideC[v]

    if not self.isBoxCoord(r, c): 
      return False
    # [r, c] coordinates valid
    i, j = self.box2pixel(r, c)  # get corresponding raster coordinates
    if self.raster.grid[i][j] == cls:
      return True
    return False


  def getRegion(self, r, c, v):
   """If v == vx.r, get region of box "above" (outside) vertex [r,c] looking right
      If v == vx.l, get region of box "left" (outside) vertex [r,c] looking left.
      The perspective is that [r,c] is the upper left corner of box at [r+1,c+1]
   """
   if self.args.verbosity >= 5:
     r0 = r+self.outsideR[v]
     c0 = c+self.outsideC[v]
     print 'getRegion: box outside [{0},{1}] in direction {2} is [{3},{4}]'.format(
		     r, c, self.diag[v], r0, c0)
   return self.box[r+self.outsideR[v]][c+self.outsideC[v]]

     
  def traverse(self): 
    """Traverse and build polygons with rings.
       The state is ([r, c], v, v0, annulus): the current vertex and the current direction
    """	  
    i, j = self.box2pixel(self.r, self.c)
    cls = self.raster.grid[i][j]    # get the classification of the raster

    # The first step is to find the region number of the box at [self.r, self.c]
    # move to upper left vertex
    r = self.r-1
    c = self.c-1
    v = vx.r   # Look outside right (east) of [r, c]

    if self.args.verbosity >= 5:
      print 'In traverse box [{0},{1}], vertex [{2},{3}] vector {4}.'.format(
	self.r, self.c, r, c, v)

    if self.isInOut(False, r, c, v, cls):
      # now it gets interesting. This box is in the same region as the box above.
      # find the region and set it
      region = self.getRegion(r, c, v)
      print 'Box above vertex [{0},{1}] is in the same region {2}.'.format(r, c, region)
      if region == 0:
        print 'ERROR ERROR: zero region zero region zero region'
	exit(-1)
      if self.isInOut(False, r, c, vx.l, cls): # look "outside" left (west) of [r, c]
	# [r, c] is an interior box. Mark with region and exit.
	print 'marking the box, and exiting'
	self.markBox(r, c, vx.r, region)
	return
      else:
	print 'Go counter clockwise -- this is an annulus.'
	print 'Do not mark boxes.'
	v = vx.l
	annulus = True
    else:
      # this is a new polygonal region. The orientation is clockwise
      self.region += 1 
      region = self.region
      annulus = False

    # the direction and the region are known
    r0, c0 = r, c  # remember the initial vertex [r0, c0]


    if self.args.verbosity >= 6:
      print 'Writing state ([{0}, {1}], {2}, {3}, {4})'.format(r,c, self.diag[v], 
		                                               cls, region)
    r, c = self.move(r, c, v)  # move in direction v
    v0 = v                     # save the previous direction

    if not self.isInOut(True, r, c, v, cls):
      v = self.cw[v]           # Inside is cw from [r, c] if not annulus
    else:
      if self.isInOut(False, r, c, v, cls):
        v = self.ccw[v]        # Outside is ccw from [r, c] if not annulus

    # mark the box after you move and have a direction
    # the boxes will be marked when you return
    #if not annulus:
    self.markBox(r, c, v, region) 

    if v != v0: # if you changed direction, output vertex
      if self.args.verbosity >= 6:
        print 'Writing state ([{0}, {1}], {2}, {3}, {4})'.format(r, c, self.diag[v], 
			                                         cls, region)

    while r != r0 or c != c0:
      r, c = self.move(r, c, v)  # move in direction v
      v0 = v                     # save the previous direction
      # define the next direction
      if not self.isInOut(True, r, c, v, cls):
        v = self.cw[v]                # cw from [r, c] if not annulus
      else:
        if self.isInOut(False, r, c, v, cls):
          v = self.ccw[v]             # ccw from [r, c] if not annulus

      # you moved and have a direction -- mark the inside box
      # the same box will be marked for each move...oh well
      #if not annulus:
      self.markBox(r, c, v, region)

      if v != v0:
        if self.args.verbosity >= 6:
          print 'Writing state ([{0}, {1}], {2}, {3}, {4})'.format(r,c,self.diag[v], 
			                                      cls, region)


class Raster(object):
  """Raster objects hold the grid of pixes, which can be passed to Dissolve objects
     without copying the entire numpy array"""
  def __init__(self, args, hdr):
    	  
    if args.verbosity >= 1:
      print "Reading array..."

    grid1D = np.fromfile(hdr.file, sep = " \n")

    # verify that the array can be reshaped
    items = grid1D.shape[0]
    if hdr.ncols * hdr.nrows != items:
      errorStr = "Number of items read in is {0} instead of {1}={2}*{3}"
      raise IOError, errorStr.format(items, hdr.nrows*hdr.ncols,
		                  hdr.nrows, hdr.ncols)

    # reshape the array
    if args.verbosity >= 1:
      print "Reshaping array to grid..."
    self.grid = np.reshape( grid1D, (hdr.nrows, hdr.ncols) )



args = parser.parse_args()  # parse command line arguments


hdr = ArcInfoGridASCII(args)
ext = ExtentHandler(hdr, args)
raster = Raster(args, hdr)


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

if args.verbosity >= 1:
  print 'Converting to shapefile...'

if not args.quiet: # show the progress bar unless instructed otherwise
  pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval = hdr.nrows).start()

if not args.dissolve:
  for row  in range(0, hdr.nrows):
    for col in range(0, hdr.ncols):
      v = raster.grid[row][col]
      if v != hdr.nodata and (not args.nonzero or v != 0):
        lon, lat = hdr.cart2geo(row, col) # Note reversal!
        if ext.compare(lat, lon):
          poly = hdr.createGridSquare(lat, lon)
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
else:
  dis = Dissolver(args, hdr, ext, raster) 
  if args.verbosity >= 5:
    print raster.grid
  while dis.nextValid(): # find the next valid box in raster/box coordinates
    dis.traverse() # traverse boundary
    if args.verbosity >= 5:
      print dis.box

if not args.quiet:
  pbar.finish()

ds.Destroy()

