#!/bin/bash

./arcinfo.py -s 1000 -b 20 -r eq \
	-e -17.338675 -34.892837  57.845763  37.428152 \
	ramankutty_cropland2000_frac_5m.asc croplandeq20.asc

./aig2shp.py -d -vv --wgs84 croplandeq20.asc x.shp

