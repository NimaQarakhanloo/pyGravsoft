#!/usr/bin/python
# $Id: geoip.py 286 2009-07-26 15:36:32Z tjansson $
"""Geoip,
Thomas R. N. Jansson,
tjansson@fys.ku.dk"""

#Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *    
    # This is the new name in python 3.

import os
from SimpleDialog import SimpleDialog
import geomodule

## Constants
programname = "GEOIP - Grid interpolation"
version = "0.3"
jobfile = "geoip.inp"
logfile = "geoip.log"
execute = "geoip"
description = """GEOIP

Gravsoft manual can be found in doc/gravsoft-manual.pdf
"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## Functions -- these needs to be defined before they are used
######################################################

        def write_to_file():
            file = open(jobfile, 'w')
            if file:
                wordlist = [\
                gfile.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                mode.textentry.get()," ",nsp.textentry.get()," ",
                rkm.textentry.get()," "]

                if latlon.textentry.get() == "0.0 0.0 0.0 0.0":
                    wordlist.extend([" f ","\n",])
                else:
                    wordlist.extend([" t ","\n", latlon.textentry.get(), "\n"])

                if mode.textentry.get() != "4":
                    wordlist.extend([pointfile.textentry.get(), "\n"])
                if int(mode.textentry.get()) > 10:
                    wordlist.extend([idno.textentry.get(), "\n"])
                if mode.textentry.get() == "4":
                    wordlist.extend([LINT.state.get(), " ",
                    griddef.textentry.get(), "\n"])

                if iellzone.textentry.get() != "0 0":
                    wordlist.extend([iellzone.textentry.get(), "\n"])

                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+
                " succesfully","normal")
            else:
                statusbar.config("ERROR writing to "+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The first section
######################################################

        #LSEL  = geomodule.NewRadioButton(frame, root,"Select only data within a subregion:", "select", "nocommand", "", "Select only data wihin a subregion (lat1-lat2, lon1-lon2) this option only works with point lists.")
        LINT = geomodule.NewRadioButton(frame, root,
        "Should grid values be integers:", "select",
        "nocommand", "", """
Specifying if the output should be in integer grid format (only needed for mode 3 and 4)""")

        #mainLabel = geomodule.mainline(frame, programname, version)
        gfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Input grid filename:", "grid.gri", """
The input file is the grid file containing data grid in standard format, i.e.
scanned in e-w lines from n to s, initiated by label
(lat1,lat2,lon1,lon2,dlat,dlon) describing grid limits.utm grids may be used,
in this case lat/lon should be northing/easting in meter, followed by
ellipsoid number (1:wgs84, 2:ed50, 3:nad27) and utm zone.""")
        mode = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Operation mode:", "1", """
1: prediction point list in geographic coordinates (degrees)
2: do, with lat and lon given in degrees, minutes, seconds
3: prediction point list in utm coordinates
4: predictions wanted in grid (geographic or utm)
5: individual prediction points in lat, lon (degrees)
6: do, with lat, lon in degrees, minutes, seconds
7: individual prediction points in utm
10: like 1, with a data value in file written after predictions
11: like 1, predictions  s u b t r a c t e d  from values
    given in file (data value must follow after the height)
12: do, but predictions  a d d e d  to values in file
13: like 11, with a second data value for each point not changed
14: like 12, do
15: 'pointfile' contains a grid, from which the
    interpolated values from 'gridfile' are  s u b t r a c t e d,
    i.e. 'outfile' = 'pointfile' - 'gridfile'
16: 'pointfile' contains a grid, to which the
    interpolated values from gridfile are  a d d e d.
17: 'pointfile' contains a grid which defines the interpolation
    points. the given grid values are only used in two-height
    interpolation mode (117), see below.
18: 'pointfile' contains a grid with unknown values (9999).
    the unknown values are interpolated from the gridfile,
    other grid values left untouched.
19: list of terrain corrections converted to rtm anomalies
    through bouguer reduction to reference level 'gfile'
20: list of free-air anomaly data converted to rtm-reduced
    data using a bouguer plate approximation only
    additional input: density
21: conversion of free-air data to Bouguer using grid
    additional input: density
22: migration of ERS-1 data over ice caps. grid is slope
    grid in degrees
31, 32, ..: like 11, 12, .. for utm coordinates in 'pointfile'

- if mode is negative, mode = abs mode and  t w o  grid
interpolation is performed (the two grids must be in
the same file, e.g. deflections of the vertical from
program 'geofour')
- if mode > 100 then point interpolation is done between
two grids in different heights using mode = mode-100.
the levels to which the two grids refer must be input.
if mode=117 then 'pointfile' is assumed
to be a height grid defining the interpolation level.""")
        nsp = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal",
        "Spline windows size:", "0", """
spline window size. 0 means bilinear interpolation,
1 is equivalent to a 8x8 spline or 8""")
        rkm = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal",
        "Minimum distance:", "0", """
Minimum required distance to closest edge (km)
for interpolation to take place.""")

        geomodule.seperator_line(frame, 10,
        "Pointfile definition (mode = 1-3, 5-7 and mode > 10)")
        pointfile = geomodule.FileSelector(frame, root, 14,
        geomodule.textwidth, "normal",
        "Point filename:", "pointfile.dat", "For mode 1,2,3 and >10")
        idno = geomodule.NewEntry(frame, root, 15,
        geomodule.textwidth, "normal",
        "Data number in record:", "1", """
Data number in line (statno,lat,lon,height,d(1),d(2)..).
For mode 11-14 and 19-22.""")


        geomodule.seperator_line(frame, 20, "Gridfile definition (mode = 4)")
        griddef = geomodule.NewEntry(frame, root, 21,
        geomodule.textwidth, "normal",
        "Output grid definition:", "1.0 5.0 2.0 6.0 0.1 0.1", "nohelp")
        #LSEL.draw(24)
        iellzone = geomodule.NewEntry(frame, root, 22,
        geomodule.textwidth, "normal",
        "    Elipsoide and zone numbers:", "0 0",
        "(0 0) indicates geographical coordinates.")
        latlon = geomodule.NewEntry(frame, root, 25,
        geomodule.textwidth, "normal", "Subgrid coordinates:",
        "0.0 0.0 0.0 0.0", """
Subregion definition lat1,lat2, lon1,lon2. This option only
works with point lists. (0.0 0.0 0.0 0.0) for no subgrid.""")
        LINT.draw(26)

        geomodule.seperator_line(frame, 30, "Height definition (mode > 100)")
        h1h2 = geomodule.NewEntry(frame, root, 31,
        geomodule.textwidth, "normal",
        "Altitude of grids:", "0 1000", "nohelp")

        geomodule.seperator_line(frame, 90, "Running options. Working in "+geomodule.getcwd)
        resultfile = geomodule.FileSelectorSave(frame, root, 93,
        geomodule.textwidth, "normal",
        "Name of file to hold interpolated values:" , "result.gri", """
Here are output on grid format is written. For point prediction both
predicted values and error estimates are written here.""")

######################################################
##The Bottom line
######################################################
        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        geomodule.quitwriterun(frame, geomodule.maxrow,
        lambda:[frame.quit()],
        lambda:[write_to_file()],
        lambda:[run_program()],
        description, root)

######################################################
## Initiate the program and start program loop
######################################################

root = Tk()
app = App(root)
root.title(programname)
root.mainloop()
