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
descript = "Create ESRI Shapefile from ArcInfo Grid ASCII raster."
epistr   = "Software is released under The MIT License (c) 2013 Florian Lengyel, CUNY Environmental CrossRoads Initiative, Advanced Science Research Center, The City College of New York. Contact: gmail/skype/twitter florianlengyel."
parser  = arg.ArgumentParser( description = descript, epilog=epistr )
# value of dest derived from first long opt
parser.add_argument('-a', '--attr',
                    metavar='attribute',
		    default='value',
		    help='Name of attribute for ArcInfo grid values. Defaults to "value."')
parser.add_argument('-d', '--dissolve',
		    action='store_true',
		    help='Dissolve ArcInfo ASCII Grid in raster space before converting to shapefile.')
parser.add_argument('-e', '--extent', 
		    nargs=4, 
		    type=float,
		    metavar=('minX', 'minY', 'maxX', 'maxY'),
		    help='Bounding box of subset of raster in geographic coordinates.')
parser.add_argument('-l', '--layer', 
                    default='grid_value',
                    help='Shapefile layer name string. Default is grid_value.')
parser.add_argument('-n', '--nonzero', 
		    action="store_true",
		    help='Exclude zero values.')
parser.add_argument('-O', '--opt',
		    action='store_true',
		    help='Enable greedy cell marking optimization.')
parser.add_argument('-q', '--quiet', 
		    action="store_true",
		    help='Suppress progress bar.')
parser.add_argument('-v', '--verbosity', 
		    action="count", default=0,
		    help="Display verbose message output. Each additional 'v' increases message verbosity: -vv is very verbose, and -vvv is very very verbose.")
parser.add_argument('--version', 
		    action='version', 
		    version='%(prog)s 0.5',
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
    
    ring.AddPoint(minX, minY)
    ring.AddPoint(minX, maxY)
    ring.AddPoint(maxX, maxY)
    ring.AddPoint(maxX, minY)
    ring.AddPoint(minX, minY) # close the ring
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
          raise ValueError, extError.format('xurcorner', hdr.xur, self.maxX, 'maxX')
        if hdr.yur < self.maxY:
          raise ValueError, extError.format('yurcorner', hdr.yur, self.maxY, 'maxY')
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

  def compare(self, lon, lat):
    return self.cmpFun(lon, lat)

# a StackExchange creation
def enum(**enums):
  return type('Enum', (), enums)

# vertex directions
vx = enum(z=0, r=1, d=2, l=3, u=4)

class PolygonDB(object):
  """List of polygons and potential holes created by Dissolver object"""

  def __init__(self, args):
    """Create the initial table of polygons. Each polygon as a region number,
    starting coordinates (r, c) and a list of potential holes with their
    starting coordinates."""

    self.db = dict()  # create empty dictionary of polygon lists

  def addPolygon(self, region, r, c):
    """Add polygon with region number 'region' and top left coordinates [r, c] to list"""
    self.db[region] = [(r-1, c-1)]

  def dump(self):
    for key in self.db:
      polyList = self.db[key]
      print 'region {0} start {1} holes {2}'.format(key, polyList[0], polyList[1:])

  def addHole(self, region, r, c):
    """Add the coordinates of the potiential hole associated to region
    The coordinates are the top left  of the box at [r,c].  Traverse counter 
    clockwise such that the outside is  the 'not region' region"""
    self.db[region].append((r-1, c-1))


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

  def __init__(self, args, hdr, ext, raster, polyDB):
    """Create a Dissolver, using 
       args   -- parsed arguments 
       hdr    -- Grid ASCII header 
       ext    -- extent of grid in coordinate space
       raster -- reshaped grid
       polyDB -- polygon/hole database"""

    self.args = args
    self.hdr  = hdr
    self.ext  = ext
    self.raster = raster  # this is a raster object
    self.polyDB = polyDB


    # Create the box coordinate space (i, j) -> [2*i+1, 2*j+1]
    # Square brackets denote box coordinates, parentheses denote pixel coordinates.
    self.boxRows  = 2*hdr.nrows + 1
    self.boxCols  = 2*hdr.ncols + 1



    # use 32 bit integers 
    # First refactoring: box coordinates are identical to
    # pixel coordinates. Translation is done at the last step.
    # This will save 4 times the space at the expense of some
    # computation initially. Subsequent refactorings will remove
    # some computation
    self.nrows = self.hdr.nrows
    self.ncols = self.hdr.ncols
    self.box   = np.zeros(shape = (self.nrows, self.ncols), dtype = np.int32) 

    # current position within the box 
    self.i =  0  # CAUTION: the nextValid() function starts from (i,j+1)
    self.j = -1  # and must be called to set the current valid pixel

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
    self.diag = { vx.z : 'z', vx.r : 'r', vx.d : 'd', vx.l : 'l', vx.u : 'u' }

    # Simply connected region to class map
    # The region outside the raster is 0 by definition. Its classification is
    # np.nan (the numpy library not a number), which is never equal to anything.
    # This is precisely the behavior one wants. The upper left hand corner is
    # in region 1, by definition. All other region numbers are defined as
    # the program runs.
    self.region2class = { 0L : np.nan }  # the outside region

    # Global locally simply connected region number
    self.region = 0L  # The box at [r,c] is marked by the region number    

  def move(self, r, c, v):
    """recompute box coordinates using the mov enum"""
    return (r + self.goR[v], c + self.goC[v])

  # Should be is VertexCoord
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
    """Find the next starting point [r,c] of bounding polygon 
       or potential annulus ring."""

    def isValid(self, i, j): # this local function isn't used elsewhere
      """True iff the pixel at (i, j) (box at [r, c]) has not been
      assigned a nonzero region number. Handling special data values goes
      elsewhere."""
      self.r, self.c = self.pixel2box(i, j)
      return self.box[i][j] == 0

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

  def isSucc(self, r, c):
    """Check if [r, c] is the successor of [self.r, self.c]"""
    #Case 1: [r, c] == [self.r, c+2]
    if self.r == r and self.c + 2 == c:
      return True
    #Case 2: [r, c] == [self.r+2, 1] and self.c == self.boxCols-2
    if r == self.r+2 and c == 1 and self.c == self.boxCols-2:
      return True
    return False

  def advCrnt(self, r, c):
    """Advance the last marked box [self.r, self.c] to [r, c].
       Use with isSuccessor()."""
    if self.isSuccessor(r, c):
      self.r, self.c = r,c
      self.i, self.j = self.box2pixel(r, c)

  def getInsideRegion(self, r, c, v):
    """Get the region of the box at vertex [r, c] clockwise 
       from direction v"""
    i, j = self.box2pixel(r+self.insideR[v], c+self.insideC[v])
    return (self.box[i][j])

  def stampBox(self, r, c, v, region):
    """Like markBox below, only without checking.
       Used by pathFix()"""
    i, j = self.box2pixel(r+self.insideR[v], c+self.insideC[v])
    self.box[i][j] = region

  def markBox(self, r, c, v, region):
    """Mark the inside  box [r0, c0] at vertex [r, c] in direction v. 
       Return True if the self.box[r0][c0] is 0 or is already marked
       region; return False if the box is nonzero and marked with
       some other region
    """
    i, j = self.box2pixel(r+self.insideR[v], c+self.insideC[v])
    reg = self.box[i][j]

    if (reg != 0 and reg != region):
      return False	# don't set the region! This is the stopping point
    self.box[i][j] = region
    # This seems to add nothing in practice
    #if self.args.opt and self.isSucc(r, c):
    #  print 'optimize!'
    #  self.advCrnt(r, c)
    return True

  def isInOut(self, isIn, r, c, v, cls):
    """Returns True if the square clockwise (ccw if isIn is False) 
       from [r, c] in  direction v exists and is in the class cls; 
       false otherwise."""
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
    return (self.raster.grid[i][j] == cls)


  # more optimized version for second pass, after all regions identified
  def isInOutReg(self, isIn, r, c, v, reg):
    """Returns True if the square clockwise (ccw if isIn is False) 
       from [r, c] in  direction v exists and is in the region reg; 
       false otherwise."""
    if isIn:   
      r += self.insideR[v]
      c += self.insideC[v]
    else:
      r += self.outsideR[v]
      c += self.outsideC[v]
    if not self.isBoxCoord(r, c): 
      return False
    # [r, c] coordinates valid
    i, j = self.box2pixel(r, c)
    return (self.box[i][j] == reg)


  # Simply Connected Region handling. The upper left-hand corner has 
  # region 1, valid or not, by  definition. If it is a nodata region, 
  # this is either treated like every other region, or no polygon is 
  # generated for it. Provisionally nodata  regions are treated the 
  # same way as all others. For now, nodata is treated by
  # the algorithm as if it were a valid region.

  def getCls(self, r, c):
    """Get class of box at [r, c]"""
    i, j = self.box2pixel(r, c)
    return self.raster.grid[i][j]

  def boxRegion(self):
    """Determine the region of the new box. Inductive assumption:
       The boxes left and above have assigned regions. The upper 
       left has region 1.  The region2class[] map is defined for the 
       boxes left and above.

       This routine needs to check whether the traversal can proceed, and in
       which direction.  Check the box above and to your left.    
       
       Left      Above    Box
       (0,nan)   (0, nan) (1, vx.r)
       (0,nan)   (R, cls) if Box.cls == cls then (R, vx,0) don't go anywhere!
                          you'll cross an edge!
			  if Box.cls != cls, then (new R, vx.r) 
       (R, cls)  (0, nan) if Box.cls == cls, then error -- you should have 
                             been here
       (R, x)    (R, x)   If Box.cls == x, then (R, vx.0) -- don't traverse
       (R, x)    (S, y)   if Box.cls == x, then (R, vx.r) you can go right
                          if Box.cls == y  then (S, vx.0) you cannot move
			  if Box.cls is neither, (new R, vx.r)
    """			  
    # left square
    rLeft = self.r
    cLeft = self.c-2
    
    rTop  = self.r-2
    cTop  = self.c

    # get i,j coordinates of these for later refactoring
    i, j = self.box2pixel(self.r, self.c)
    iLeft, jLeft = self.box2pixel(rLeft, cLeft)
    iTop, jTop = self.box2pixel(rTop, cTop)

    myCls = self.getCls(self.r, self.c)

    if not self.isBoxCoord(rLeft, cLeft):
      # left invalid, look above
      if not self.isBoxCoord(rTop, cTop):
	self.region += 1
        self.polyDB.addPolygon(self.region, self.r, self.c)
	self.region2class[self.region] = myCls
        self.box[i][j] = self.region  
        return (self.region, vx.r)    
      else:
	# Left undefined, Top is defined -- get cls
	if self.getCls(rTop, cTop) == myCls:
	  region = self.box[iTop][jTop]
          self.box[i][j] = region
	  # you cannot move right -- this would mean crossing an edge
	  return (region, vx.z)
        else:
	  # my class is new at the left edge.
	  self.region += 1
          self.polyDB.addPolygon(self.region, self.r, self.c)
	  self.region2class[self.region] = myCls  # add the new region
          self.box[i][j] = self.region
	  return (self.region, vx.r)

    # Left is defined. Suppose the top is not defined. This is a sanity check

    leftCls = self.getCls(rLeft, cLeft)
    if not self.isBoxCoord(rTop, cTop):
      if leftCls == myCls:
	# This is an error. The traversal should have found this
	raise ValueError, 'Box [{0},{1] has same class {2} as [{3},{4}]'.format(
			                         self.r, self.c, myCls, rLeft, cLeft)
      else:
      # New class
        self.region += 1
        self.polyDB.addPolygon(self.region, self.r, self.c)
	self.region2class[self.region] = myCls  # add the new region
        self.box[i][j] = self.region
	return (self.region, vx.r)

    # finally, the top and left are defined 
    topRegion  = self.box[iTop][jTop]
    leftRegion = self.box[iLeft][jLeft]

    try:
      if myCls == self.region2class[topRegion]:
        self.box[i][j] = topRegion
        # you cannot move right -- crossing deleted edge
        return (topRegion, vx.z)
    except KeyError:
      print 'EXCEPTION: box[{0}][{1}] topRegion {2} top box [{3},{4}] not in map, myCls {5}'.format( self.r, self.c, topRegion, rTop, cTop, myCls)
      print 'self.box[{0}][{1}] = {2} next available region {3}'.format(
		      rTop,cTop,self.box[iTop][jTop], self.region)
      print self.region2class
      exit(-2)

    # my class differs from the top region  

    try:
      if myCls == self.region2class[leftRegion]:
        self.box[i][j] = leftRegion
        return (leftRegion, vx.r)

    except KeyError:
      print 'EXCEPTION: box[{0}][{1}] topRegion {2} top box [{3},{4}] not in map, myCls {5}'.format( self.r, self.c, topRegion, rTop, cTop, myCls)
      print 'self.box[{0}][{1}] = {2} next available region {3}'.format(
		      rTop,cTop,self.box[iTop][jTop], self.region)
      print self.region2class
      exit(-3)


    # my class starts a new region -- this could join with others.
    self.region += 1
    self.polyDB.addPolygon(self.region, self.r, self.c)
    #if self.region2class[leftRegion] == self.region2class[topRegion]:
    if leftRegion == topRegion: # left and top region are defined
      # to be sure, find the diagonally oppostite region at [i-1, j-1]
      if leftRegion == self.box[i-1, j-1]:
        # this is potentially a hole -- you may have to undo
        self.polyDB.addHole(leftRegion, self.r, self.c)  
    self.region2class[self.region] = myCls  # add the new region
    self.box[i][j] = self.region
    return (self.region, vx.r)



  def pathFix(self, r0, c0, r1, c1, region, cls):
    """Your optimism was unfounded. Update the boxes
       along the path from [r0, c0] to [r1, c1] with 
       the correct region 'region'. Remove the entry
       from the polygon database and revert the region
       number.
    """
    v = vx.r  # you went right in this case at [r0, c0]

    r, c = r0, c0  # start of path
    self.stampBox(r, c, v, region) # fix boxRegion()

    r, c = self.move(r, c, v)  # move in direction v
    v0 = v                     # save the previous direction

    if not self.isInOut(True, r, c, v, cls):
      v = self.cw[v]           # Inside is cw from [r, c] 
    else:
      if self.isInOut(False, r, c, v, cls):
        v = self.ccw[v]        # Outside is ccw from [r, c] 
	# stamp the box -- forget about collisions
	self.stampBox(r, c, v0, region)  
          
    if v != v0: # if you changed direction, output vertex
      if self.args.verbosity >= 5:
	print 'Turn (pre while): ([{0}, {1}], {2}, {3}, {4})'.format(r, c, 
			                                         self.diag[v], 
			                                         cls, region)
    # mark the box inside in the direction v0 you came from [r,c]
    # check this!!! v0 or v?
    self.stampBox(r, c, v, region)

    while r != r1 or c != c1:
      r, c = self.move(r, c, v)  # move in direction v
      v0 = v                     # save the previous direction
      # define the next direction
      if not self.isInOut(True, r, c, v, cls):
        v = self.cw[v]           # cw from [r, c] if not annulus
      else:
        if self.isInOut(False, r, c, v, cls):
          v = self.ccw[v]        # ccw from [r, c] if not annulus
	  self.stampBox(r, c, v0, region)  

      if v != v0:
        if self.args.verbosity >= 5:
	  print 'Turn (while): ([{0}, {1}], {2}, {3}, {4})'.format(
			  r, c, self.diag[v], cls, region)
      # mark the box inside in the direction v from [r,c]
      self.stampBox(r, c, v, region)
    
    # Reset the region number -- this should be
    # handled as a "region number service."
    # Set the polygon to an invalid entry
    self.polyDB.db[self.region] = [(-1,-1)] # make empty
    # del self.polyDB.db[self.region]
    self.region -= 1




  def traverse(self): 
    """Traverse and build polygons with rings. The state is 
       ([r, c], v, v0): the current vertex and the current and previous 
       directions. The previous direction is needed to mark regions 
       during a ccw turn.

       The boxRegion() algorithm is optimistic and will sometimes identify
       new regions that should be combined with previously defined regions.
       The pathFix() algorithm will undo a misadventure of this kind. Also,
       the boxRegion() algorithm will sometimes misidentify holes that 
       are outside the associated region. The isHole() algorithm identifies
       these by checking if the path coincides with the boundary (it will 
       pass through the origin of the path). The alternative to optimistic 
       traversal and marking is to completely identify a region.
    """	  

    # get the classification of the new box
    (region, v) = self.boxRegion()  # All boxes visited will be marked with region 
    if region < 0:
      print 'FATAL ERROR: region # {0} < 0. Too many regions.'.format(region)
      exit(-1)

    # the class is needed to determine direction
    cls =  self.getCls(self.r, self.c)  

    # move to upper left corner
    r = self.r-1
    c = self.c-1

    verbosity = self.args.verbosity
    if verbosity >= 5:
      print 'Traverse box [{0},{1}], vertex [{2},{3}] vector {4} region {5} class {6}.'.format(
	self.r, self.c, r, c, self.diag[v], region, cls)

    if v == vx.z:  # can't proceed
      if verbosity >=7:
	print 'Traverse box: exiting, cannot go right'
      return


    # the direction and the region are known
    r0, c0 = r, c  # remember the initial vertex [r0, c0]


    r, c = self.move(r, c, v)  # move in direction v
    v0 = v                     # save the previous direction

    if not self.isInOut(True, r, c, v, cls):
      v = self.cw[v]           # Inside is cw from [r, c] 
    else:
      if self.isInOut(False, r, c, v, cls):
        v = self.ccw[v]        # Outside is ccw from [r, c] 
	# mark box -- keep track of collisions
	if not self.markBox(r, c, v0, region):  
          if verbosity >= 7:
	    print 'Mark box:collision at [{0},{1}] should be {2}'.format(
			  r, c, self.getInsideRegion(r, c, v0))
	    print 'Fixup path from [{0},{1}] to [{2},{3}]'.format(
			  r0, c0, r, c)
	  self.pathFix(r0, c0, r, c, self.getInsideRegion(r, c, v0), cls)
	  return
          
    if v != v0: # if you changed direction, output vertex
      if verbosity >= 5:
	print 'Turn (pre while): ([{0}, {1}], {2}, {3}, {4})'.format(r, c, 
			                                         self.diag[v], 
			                                         cls, region)
    # mark the box inside in the direction v0 you came from [r,c]
    if not self.markBox(r, c, v, region):
      if verbosity >= 7:
        print 'Mark box:collision at [{0},{1}] should be {2}'.format(
			  r, c, self.getInsideRegion(r, c, v))
        print 'Fixup path from [{0},{1}] to [{2},{3}]'.format(
			  r0, c0, r, c)
      self.pathFix(r0, c0, r, c, self.getInsideRegion(r, c, v), cls)
      return

    while r != r0 or c != c0:
      r, c = self.move(r, c, v)  # move in direction v
      v0 = v                     # save the previous direction
      # define the next direction
      if not self.isInOut(True, r, c, v, cls):
        v = self.cw[v]           # cw from [r, c] if not annulus
      else:
        if self.isInOut(False, r, c, v, cls):
          v = self.ccw[v]        # ccw from [r, c] if not annulus
	  verbosity = self.args.verbosity
	  if verbosity >= 6:
	    print 'While: going ccw from {0} to {1}'.format(
			    self.diag[v0], self.diag[v])
	  if not self.markBox(r, c, v0, region):  
	    if verbosity >= 7:
              print 'Mark box: ccw collision at [{0},{1}] should be {2}'.format(
			  r, c, self.getInsideRegion(r, c, v0))
              print 'Fixup path from [{0},{1}] to [{2},{3}]'.format(
			  r0, c0, r, c)
            self.pathFix(r0, c0, r, c, self.getInsideRegion(r, c, v0), cls)
	    return

      if v != v0:
        if self.args.verbosity >= 5:
	  print 'Turn (while): ([{0}, {1}], {2}, {3}, {4})'.format(
			  r, c, self.diag[v], cls, region)
      # mark the box inside in the direction v from [r,c]
      if not self.markBox(r, c, v, region):
	if verbosity >= 7:
          print 'Mark box: cw/0 collision at [{0},{1}] should be {2}'.format(
			  r, c, self.getInsideRegion(r, c, v))
          print 'Fixup path from [{0},{1}] to [{2},{3}]'.format(
			  r0, c0, r, c)
        self.pathFix(r0, c0, r, c, self.getInsideRegion(r, c, v),cls)
	return

# Let us pray

  def isHole(self, r, c, region, r1, c1): 
    """Returns True iff [r, c] are the coordinates of a vertex (specifically, 
       the upper right vertex of the upper left most pixel) of an inner ring 
       of the region r. Traverse hole counterclockwise, using the potential 
       hole list from the first pass.  The state is ([r, c], v, v0): the 
       current vertex and the current and previous directions. The previous 
       direction is used to verify regions during a ccw turn.

       By moving counterclockwise, the same traversal algorithm as traverse() 
       above will traverse a "hole". This time a check needs to be performed 
       at each counterclockwise turn. Checks for the region are unnecessary 
       when going "straight" or clockwise, since the isInOutReg() predicates 
       have checked the region number.  The starting configuration is the 
       upper right vertex [r, c] of the upper rightmost pixel of a potential 
       hole, in which the top, diagonally upper left and left pixels all have 
       the same region R.

       A gotcha: traversal of a potential hole could lead you outside the 
       region.  This happens if you move to the coordinates of any vertex 
       of the bounding polygon.  Since the polygon database maintains the 
       start vertex of the region, this can be checked.

       R   R  [r, c]
       R   *
    """	  

    v = vx.d # go counterclockwise from upper left corner
    # the direction and the region are known
    r0, c0 = r, c  # remember the initial vertex [r0, c0]

    r, c = self.move(r, c, v)  # move in direction v
    if r == r1 and c == c1:    # check if on bounding polygon of region
      return False
    v0 = v                     # save the previous direction

    if not self.isInOutReg(True, r, c, v, region):
      v = self.cw[v]           # Inside is cw from [r, c] 
    else:
      if self.isInOutReg(False, r, c, v, region):
        v = self.ccw[v]        # Outside is ccw from [r, c] 
	#  After a ccw move, verify the inside and outside boxes
	#  from [r, c] in the previous direction
	if ((not self.isInOutReg(True, r, c, v0, region)) or
		(not self.isInOutReg(False, r, c, v0, region))):
	  return False

    if v != v0: # if you changed direction, output vertex
      if self.args.verbosity >= 4:
	print 'isHole: Changed direction: ([{0}, {1}], {2}, {3})'.format(
                                                 r, c, self.diag[v], region )
    while r != r0 or c != c0:
      r, c = self.move(r, c, v)  # move in direction v
      if r == r1 and c == c1:    # outside since you passed through [r1, c1]
	return False
      v0 = v                     # save the previous direction
      # define the next direction
      if not self.isInOutReg(True, r, c, v, region):
        v = self.cw[v]           # cw from [r, c] if not annulus
      else:
        if self.isInOutReg(False, r, c, v, region):
          v = self.ccw[v]        # ccw from [r, c] if not annulus
	  #  After a ccw move, verify the inside and outside boxes
	  #  from [r, c] in the previous direction
	  if ((not self.isInOutReg(True, r, c, v0, region)) or
		(not self.isInOutReg(False, r, c, v0, region))):
	    return False

      if v != v0:
        if self.args.verbosity >= 4:
	  print 'isHole: Changed direction: ([{0}, {1}], {2}, {3})'.format(
                                r, c, self.diag[v], region)
    # you made it!	  
    return True	  

  def box2geo(self, r, c):
    """Convert vertex at [r, c] to geographic coordinates (x, y) = (lon, lat). 
       Assumes that r and c are both even--which is true for vertices of boxes."""
    return ( self.hdr.xll + self.hdr.cell * (c >> 1),
             self.hdr.yll + self.hdr.cell * (self.hdr.nrows - (r >> 1)) )


  def addRing(self, poly, region, r0, c0, v):
    """Add ring with starting coordinates [r0, c0] to polygon poly.
       Traverse in the direction v. (This cries out for a general traversal 
       class.)"""
    
    if self.args.verbosity >= 5:
      print 'AddRing'

    ring = ogr.Geometry(ogr.wkbLinearRing)

    # add the origin
    lon, lat = self.box2geo(r0, c0)
    ring.AddPoint(lon, lat)    # note reversal

    r, c = r0, c0  # start of path
    v0 = v                     # save the previous direction

    r, c = self.move(r, c, v)  # move in direction v

    if not self.isInOutReg(True, r, c, v, region):
      v = self.cw[v]           # Inside is cw from [r, c] 
    else:
      if self.isInOutReg(False, r, c, v, region):
        v = self.ccw[v]        # Outside is ccw from [r, c] 
          
    if v != v0: # if you changed direction, output vertex
      lon, lat = self.box2geo(r, c)
      ring.AddPoint(lon, lat)    # note reversal

    while r != r0 or c != c0:
      r, c = self.move(r, c, v)  # move in direction v
      v0 = v                     # save the previous direction
      # define the next direction
      if not self.isInOutReg(True, r, c, v, region):
        v = self.cw[v]           # cw from [r, c] if not annulus
      else:
        if self.isInOutReg(False, r, c, v, region):
          v = self.ccw[v]        # ccw from [r, c] if not annulus

      if v != v0:
        lon, lat = self.box2geo(r, c)
        ring.AddPoint(lon, lat)    # note reversal

    ring.CloseRings()
    # Add the ring to the polygon
    poly.AddGeometry( ring )



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


if __name__ == '__main__':  
  args = parser.parse_args()  # parse command line arguments
  hdr = ArcInfoGridASCII(args)
  ext = ExtentHandler(hdr, args)
  raster = Raster(args, hdr)
  verbosity = args.verbosity


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

  if verbosity >= 1:
    print 'Converting to shapefile...'


  if not args.dissolve:
    if not args.quiet: # show the progress bar unless instructed otherwise
      pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval = hdr.nrows).start()
    for row  in range(0, hdr.nrows):
      for col in range(0, hdr.ncols):
        v = raster.grid[row][col]
        if v != hdr.nodata and (not args.nonzero or v != 0):
          lon, lat = hdr.cart2geo(row, col) # Note reversal!
          if ext.compare(lon, lat):
            poly = hdr.createGridSquare(lon, lat)
            feature = ogr.Feature( layer.GetLayerDefn() )
            feature.SetField(args.attr, v)
            feature.SetGeometry(poly) # set the attribute
            if layer.CreateFeature(feature):
              raise ValueError, "Could not create feature in shapefile."
            feature.Destroy()
      if not args.quiet:
        pbar.update(row+1)
    if not args.quiet:
      pbar.finish()
  else:
    polyDB = PolygonDB(args) # create empty polygon database

    # create a polygon dissolver
    dis = Dissolver(args, hdr, ext, raster, polyDB) 

    if verbosity >= 6:
      print raster.grid

    # this is the dissolve algorithm
    while dis.nextValid(): # find the next valid box in raster/box coordinates
      dis.traverse() # traverse boundary

      if verbosity >= 7:
        print dis.box

    # dump the polygons
    for region in polyDB.db:
      polyList = polyDB.db[region]
      r0, c0 = polyList[0]  # obtain start coordinates of polygon

      if r0 == -1 and c0 == -1:
        continue  # ignore bad region
      
      if verbosity >= 6:
        print 'region {0} start {1} holes {2}'.format(
                                    region, polyList[0], polyList[1:])

      # otherwise obtain the class (bin number) of the region
      cls = dis.region2class[region]

      # ignore nodata polygons for now
      if cls == hdr.nodata:
        if verbosity >= 6:
	  print 'NODATA region'
      else:
        # create the polygon
        # first poly
        poly = ogr.Geometry(ogr.wkbPolygon)

        # add the outer clockwise ring starting at [r0, c0]
        dis.addRing(poly, region, r0, c0, vx.r)

        for vertex in polyList[1:]:
	  r, c = vertex

	  if dis.isHole(r, c, region, r0, c0):
            if verbosity >= 6:
	      print '[{0},{1}] is a hole of [{2},{3}]'.format( r, c, r0, c0 )

            # add inner counterclockwise ring at [r, c]
            dis.addRing(poly, region, r, c, vx.d)
      
        # add the class as an attribute
        feature = ogr.Feature( layer.GetLayerDefn() )
        feature.SetField(args.attr, cls)
        feature.SetGeometry( poly ) # set the attribute
        if layer.CreateFeature( feature ):
          raise ValueError, "Could not create feature in shapefile."
        feature.Destroy()

    if args.verbosity >= 7:
      print dis.box

  # write out the shapefile
  ds.Destroy()
