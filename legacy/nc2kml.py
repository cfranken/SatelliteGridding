#!/usr/bin/env python

# 2012 Christian Frankenberg, christian.frankenberg@jpl.nasa.gov
#
# Version 1 - 10-July-2012 Initial first try
#

import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from matplotlib.ticker import NullFormatter
import matplotlib.mlab as mlab
import scipy.stats as st
import os
#import kmlcircle as k
from simplekml import *
from optparse import OptionParser
#import tai64n
import datetime
from datetime import timedelta

# Dictionary defining what L2 fields to read:
dict_l2 = {'var': "/ch4", 'lat_': "/lat",'lon_': "/lon", 'n':"n"}

#fields to use in the KML file
dict_fields = {'var':'xCH4'}
  
# Paul Tol's colorbar:
# Linear colormap (rainbow)
cols = [(0,0,0)]
for x in np.linspace(0,1, 255):
  rcol = ((0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3))
  gcol = (0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6)
  bcol = (1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5))
# print rcol, gcol, bcol
  cols.append((rcol, gcol, bcol))
  
# rcol = 255*((0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3))
# gcol = 255*(0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6)
# bcol = 255*(1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5))
# cols.append((rcol, gcol, bcol))

def bytsclGen(x, minvalue=None, maxvalue=None, top=255.0):
  if x<minvalue: return 0
  elif x>maxvalue: return top
  else: return int((x-minvalue)*top/(maxvalue-minvalue))

def plotOrbit(options, args):
        
    kml = Kml(name='Level 3 viewer')
    # doc = kml.newdocument(name="ACOS-GOSAT Level 3 XCO2")
    fol = kml.newfolder()
    fol2 = kml.newfolder()
    fol.name = "L3 gridded averages"
    fol.description = "lat/lon gridded averages"
   
    schema = kml.newschema(name='Stats')
    for i in dict_fields:
        schema.newsimplefield(name=i, type='float', displayname=dict_fields[i])

  
    for file in args:
        print("Analyzing file " + file)
        l2_obj = h5py.File(file, "r")
  
        try:
            year = float(file[-10:-6])
            month = float(file[-5:-3])
            time_str = file[-10:-6] + '-'+ file[-5:-3]
            useTime=True
        except:
            year = 0
            useTime=False
            time_str =''
        print(year, time_str)
              # Read all necessart datasets
        process = True
        for i in dict_l2:
            cmd = i + '=l2_obj["' +dict_l2[i]+'"][:]'
            try:
                exec(cmd)
                print(cmd)
            except:
                   print('Error executing ' + cmd + 'in ' + file)
                   process = False
        # Shortcut:
        var=l2_obj["/ch4"][:]
        lat_=l2_obj["/lat"][:]
        lon_=l2_obj["/lon"][:]
        n=l2_obj["n"][:]

        if process:
            print(var.shape)
            nz,ny,nx = var.shape
            dx = lon_[2]-lon_[1]
            dy = lat_[2]-lat_[1]
            print('Dimensions (lat*lon): ' +  str(nx) + '  ' +  str(ny))
        else:
            nx = 0
            ny = 0
            dx = 0
            dy = 0
        wo = np.where((n>0.5))
        print((wo[0]), " ", (wo[1]), " ",(wo[2]) )
        #print(size(n))
        n_val = len(wo[0])
        for i in range(n_val):
            #print(wo[2][i])
            #print(wo[1][i])
            lon = lon_[wo[1][i]]-dx/2.
            lat = lat_[wo[0][i]]-dx/2.
            #print(var.shape)
            # print ix, iy, lon, lat
            if lon>options.lonMin and lon<options.lonMax and lat>options.latMin and lat<options.latMax and np.isfinite(var[wo[0][i],wo[1][i],wo[2][i]]) and var[wo[0][i],wo[1][i],wo[2][i]]!=-999.9:
                pol = fol.newpolygon(name='grid box')
                pol.outerboundaryis=[(lon,lat),(lon+dx, lat), (lon+dx, lat+dy), (lon, lat+dy),(lon,lat) ]
          #print np.isfinite(var[wo[0][i],wo[1][i],wo[2][i]])
                cIndex = bytsclGen(var[wo[0][i],wo[1][i],wo[2][i]], minvalue=options.min, maxvalue=options.max, top=254)
                cColor = 'dd%02x%02x%02x' % (int(cols[cIndex][2]*254), int(cols[cIndex][1]*254), int(cols[cIndex][0]*254))
                pol.polystyle.color = cColor  # Transparent red
                pol.polystyle.outline = 0
                pol.extendeddata.schemadata.schemaurl = schema.id
              #for j in dict_fields:
          # cmd = 'pol.extendeddata.schemadata.newsimpledata("'+j+'",'+j+'[ix,iy])'
          # exec cmd
             
              #if useTime:
                 
              #    pol.timestamp = time_str
    print('Saving to ' + options.outFile + ', can take some time')
    kml.save(options.outFile)
def standalone_main():
    # Parsing options
    parser = OptionParser(usage="usage: %prog -o outfile.kml level2_file(s) ")
    parser.add_option(  "-o", "--output", dest="outFile",
                       metavar="FILE",
                       default='./orbit_result.kmz',
                       help="Name of output figure file to be used in this run")
    parser.add_option(  "--min", dest="min",
                       type='float',
                       default=0,
                       help="Lower range of colorbar")
    parser.add_option(  "--max", dest="max",
                       type='float',
                       default=1.5,
                       help="Upper range of colorbar")
    parser.add_option(  "--latMin", dest="latMin",
                       type='float',
                       default=-90,
                       help="Lower range of latitude")
    parser.add_option(  "--latMax", dest="latMax",
                       type='float',
                       default=90.0,
                       help="Upper range of latitude")
    parser.add_option(  "--lonMin", dest="lonMin",
                       type='float',
                       default=-179,
                       help="Lower range of longitude")
    parser.add_option(  "--lonMax", dest="lonMax",
                       type='float',
                       default=179.0,
                       help="Upper range of longitude")
    parser.add_option( "--all", dest="plotall", action="store_true",default=False,help="set this if you want ALL sounding ids to be plotted (can be lots of files!)")

    # Parse command line arguments
    (options, args) = parser.parse_args()
    # call the function defined on top
    plotOrbit(options, args)
        
if __name__ == "__main__":
    standalone_main()



#no guarantee, Christian Frankenberg
