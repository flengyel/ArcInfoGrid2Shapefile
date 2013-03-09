#!/usr/bin/python
# -*- coding: utf-8 -*-
# DO NOT REMOVE THE SECOND LINE!
# aig2shp.py
# Author: Florian Lengyel
# Date: February 17, 2013
# Contact: gmail/skype/twitter florianlengyel
# License: MIT License (c) 2013 Florian Lengyel
 
from __future__ import division
import os, sys
import numpy as np
from  scipy import stats, special
import argparse as arg



class ArcInfoGridASCII(object):
  """Header corresponding to Arc Info Grid ASCII file"""

  def getField(self, fieldname, typeconv):
    """Expect fieldname + ':' and data of type typeconv. 
       Convert data to typeconv and return to caller."""
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

  def putField(self, fd, fieldname, val):
    """Write field to file as text"""
    fd.write('{0:<14}{1!s}\n'.format(fieldname, val))

  def writeReclass(self, args, ext, raster, endpoints):
    """Write the reclassified file"""	    
    with open(args.outfile, 'w') as fd:
      # Write header
      self.putField(fd, 'ncols', self.ncols)
      self.putField(fd, 'nrows', self.nrows)
      self.putField(fd, 'xllcorner', self.xll)
      self.putField(fd, 'yllcorner', self.yll)
      self.putField(fd, 'cellsize', self.cell)
      self.putField(fd, 'NODATA_value', int(self.nodata)) # must be int
      # Write reclassified data based on endpoints
      for row  in range(0, self.nrows):
        for col in range(0, self.ncols):
          v = raster.grid[row][col]
	  if args.extent != None:
	    lon, lat = self.cart2geo(row, col)
	    if not ext.compare(lon, lat):
	      v = self.nodata
          if v != self.nodata:
	    if args.scaling > 1:
	      v = int(v * args.scaling)
            cat = np.searchsorted(endpoints, v, side='right')
	  else:
	    cat = int(self.nodata)
          fd.write('%i ' % cat)
        fd.write('\n')  	


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
    

#handle extents by converting to rows and columns, and limiting ranges.
#But for now, don't handle extents.

class ExtentHandler(object):
  """Handle bounding boxes within ArcInfo Grid ASCII extents. Sets comparison
     function depending on arguments supplied to command line"""
  def __True__(self, lon, lat):  # use this comparison when -e, --extent is absent
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

  def compare(self, lon, lat):
    return self.cmpFun(lon, lat)

if __name__ == '__main__':
  from progressbar import ProgressBar, Percentage, Bar

  class Raster(object):
    """Raster objects hold the grid of pixels, which can be passed to Dissolve objects
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

  # Define the arguments first
  descript = "Create reclassified ArcInfo Grid ASCII raster. Produces an integer valued ArcInfo Grid ASCII file."
  epistr   = "Software is released under The MIT License (c) 2013 Florian Lengyel, CUNY Environmental CrossRoads Initiative, Advanced Science Research Center, The City College of New York. Contact: gmail/skype/twitter florianlengyel."
  parser  = arg.ArgumentParser( description = descript, epilog=epistr )
  # value of dest derived from first long opt
  parser.add_argument('-b', '--bins',
                    type=int,
		    default=20,
		    help='Number of bins. The max is 127.')
  parser.add_argument('-e', '--extent', 
		    nargs=4, 
		    type=float,
		    metavar=('minX', 'minY', 'maxX', 'maxY'),
		    help='Bounding box of subset of raster in geographic coordinates.')
  parser.add_argument('-p','--prt',
		     action='store_true',
		     help='Print bin intervals and exit.')
  parser.add_argument('-q', '--quiet', 
		    action='store_true',
		    help='Suppress progress bar.')
  parser.add_argument('-r', '--reclass',
		      nargs=1,
		      choices=['eq', 'hist'],
		      default='eq',
		      help='Reclassification. Equal area histogram (quantiles), equal division binning (histogram).')
  parser.add_argument('-s', '--scaling',
		      type=int,
		      default=1,
		      help='Scale raster values by multiplier (usually a factor of 10) and take integer part.')
  parser.add_argument('-v', '--verbosity', 
		    action="count", default=0,
		    help="Display verbose message output. Each additional 'v' increases message verbosity: -vv is very verbose, and -vvv is very very verbose.")
  parser.add_argument('--version', 
		    action='version', 
		    version='%(prog)s 0.1',
		    help='Show program version number and exit.')
  parser.add_argument('infile', 
		    metavar='input_ArcInfo_Grid_ASCII_file',
		    help='ArcInfo Grid ASCII input file.')
  parser.add_argument('outfile',
		    metavar='output_ArcInfo_Grid_ASCII_file',
		    help='output ArcInfo Grid ASCII file')

  args = parser.parse_args()  # parse command line arguments
  hdr = ArcInfoGridASCII(args)
  ext = ExtentHandler(hdr, args)
  raster = Raster(args, hdr)

  if not args.quiet: # show the progress bar unless instructed otherwise
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval = hdr.nrows).start()

  # pass one: compute the bins. May use interesting algorithms.
  a = [] # list of valid values.
  for row  in range(0, hdr.nrows):
    for col in range(0, hdr.ncols):
      v = raster.grid[row][col]
      if args.extent != None:
        lon, lat = hdr.cart2geo(row, col)
	if not ext.compare(lon, lat):
          v = hdr.nodata # act as if nodata
      if v != hdr.nodata:
	if args.scaling > 1:
	  v = int(v * args.scaling)
        a.append(v)
    if not args.quiet:
      pbar.update(row+1)

  x = np.array(a) # make numpy array
  x.sort()

  endpoints = []
  if args.reclass[0] == 'eq':
    # get quantiles.
    for i in range(1, args.bins):
      endpoints.append(stats.scoreatpercentile(x, i * int(100/args.bins))) 
    endpoints.append(stats.scoreatpercentile(x, 100))

  if args.reclass[0] == 'hist':
    counts,endpoints = np.histogram(x, bins=args.bins)

  if args.prt:
    print endpoints
    exit(0)

  hdr.writeReclass(args, ext, raster, endpoints)

  if not args.quiet:
    pbar.finish()


