#!/usr/bin/python
# Convert BC-Consulting QGIS Color Table plugin syntax to CartoDB.CSS
# Author: Florian Lengyel
# Date:   Feb 18, 2013
# License: BC-Consult's code is GPL'd, hence so is this.

import sys
import argparse as arg

descript = "Convert ArcGIS or BC-Consulting QGIS Color Table plugin syntax to CartoDB.CSS"
epistr   = "Software is released under the GPL (c) 2013 Florian Lengyel, CUNY Environmental CrossRoads Initiative, Advanced Science Research Center, The City College of New York. Contact: gmail/skype/twitter florianlengyel."
parser  = arg.ArgumentParser( description = descript, epilog=epistr )
parser.add_argument('--clr',
		    help='Use ArcInfo .clr format',
		    action='store_true')
parser.add_argument('--ge', 
		    help='Use [column_name >= value] instead of >',
		    action='store_true')
parser.add_argument( '--polygon_opacity',
		     help='Opacity of layer. Default is 0.8',
		     default=0.8,
		     type=float )
parser.add_argument( '--line_opacity',
		     help='Opacity of line. Default is 0',
		     default=0,
		     type=float )
parser.add_argument( '-c', '--column_name', 
		     help='Name of CartoDB column',
		     default='value' )
parser.add_argument( '-s', '--scaling',
		     type=int,
		     default=1,
		     help='Scale values by int multiple; e.g., a power of 10.')
parser.add_argument( 'table_name', 
		     help='Name of CartoDB table')
parser.add_argument( 'color_table', 
		     type=arg.FileType('r'),
		     help='Name of BC-Consult Color Table file' )


class CartoCSS(object):
  def __init__(self, args):
    self.openBraces = 1  # start with one openbrace
    print '#{0} {{'.format(args.table_name)
    print '  polygon-opacity:{0};'.format(args.polygon_opacity)
    print '  line-opacity:{0};'.format(args.line_opacity)

  def entry(self, args, map):
    gt = '>'
    if args.ge:
      gt = '>='
    print '    [{0} {1} {2}] {{'.format(args.column_name, gt, map['value'])
    print '      polygon-fill:rgb({0},{1},{2});'.format( map['red'], 
	                                                 map['green'], 
							 map['blue'])
    print '       polygon-clip:false;'
    self.openBraces += 1  # close the braces at the end.

  def closeBraces(self):
    for i in range(0, self.openBraces):
      sys.stdout.write('}')

args = parser.parse_args()
css = CartoCSS(args)

# discard the first line
# Format of BC consult includes a header line
if not args.clr:
  args.color_table.readline()

# if ArcInfo .clr format, the separator is a space
if args.clr:
  sep = ' '
else:
  sep = ','

for line in args.color_table:
  field = line.split(sep)
  map = {  'value': str(int(field[0].strip())*args.scaling), 
	   'red'  : field[1].strip(),
	   'green': field[2].strip(),
	   'blue' : field[3].strip()}
  css.entry(args, map)

css.closeBraces()

