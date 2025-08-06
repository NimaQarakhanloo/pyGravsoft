#!/usr/bin/python
# $Id: geofour.py 286 2009-07-26 15:36:32Z tjansson $
"""Geofour,
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
programname = "GEOFOUR - Planar FFT for gravity field modelling"
version = "0.3"
jobfile = "geofour.inp"
logfile = "geofour.log"
execute = "geofour"
description = """GEOFOUR

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
                hfile.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                mode.textentry.get()," ",attkm.textentry.get(),"\n",
                swcorner.textentry.get()," ",inne.textentry.get()," ",
                iwndow.textentry.get(),"\n"]

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

        #mainLabel = geomodule.mainline(frame, programname, version)

        gfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Input data filename:",
        "grid.gri", """
The input file is the grid file containing data grid in standard format,
i.e. scanned in e-w lines from n to s, initiated by label
(lat1,lat2,lon1,lon2,dlat,dlon) describing grid limits. utm grids may be used,
in this case lat/lon should be northing/easting in meter, followed by
ellipsoid number (1:wgs84, 2:ed50, 3:nad27) and utm zone.""")
        hfile = geomodule.FileSelector(frame, root, 3,
        geomodule.textwidth, "normal", "Input height grid file:",
        "hgrid.gri", """
The height grid file is a grid file on text format corresponding to the input
file. It is not needed in the simple mode (dummy name must be specified). The
height grid must have the same spacings and relative position as the input
file grid, and must have at least the wanted area in common.""")
        mode = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal", "Operation mode:", "1", """
1  conversion gravity to geoid (stokes formula).
2  gravity to deflections (vening meinesz).
3  conversion geoid to gravity.
4  gravity to tzz (eotvos)
5  gravity to plumbline curvatures tyz,txz (eotvos)
6  upward continuation to h (km) (NB! point innerzone pt)
7  Defections (ksi,eta) to gravity
8  ksi (arcsec) to gravity
9  eta (arcsec) to gravity
10 downward continuation of gravity data to sea level.
11 gravity to geoid, using harmonic continuation to a mean
   height reference level, followed by stokes and continuation
   of computed geoid to surface of topography.
12 gravity to deflections using harmonic continuation.
   continuation of deflections from reference level done by
   using plumbline curvature parameters txz and tyz.
0  nothing, gravity wiener filtering only if attkm<0 a
   high-pass filtering opreation is done""")
        attkm = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal",
        "Wiener filter resolution:", "0.0", """
Wiener filtering resolution for noisy data (kaula rule).
This option should be used only for high-pass filtering
operations (modes.ge.3).\nIt specifies the resolution (km)
where the wiener attenuation filter is 0.5. If equal to
0.0 it means no attenuation is used. Values < 0: special
filtering, should only be used in when mode equal 0.""")
        swcorner = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal", "South-West corner position:",
        "0.0 0.0", """
(degrees (or m for utm)) sw corner of wanted
subgrid if both coordinates are equal to 0 then the sw-corner
of the input grid is used.""")
        inne = geomodule.NewEntry(frame, root, 7,
        geomodule.textwidth, "normal", "Number of points in subgrid:",
        "200 256", """
Number of points of subgrid in northern and easteren direction
(must be even numbers) if the wanted subgrid is bigger than the
actual grid the transform grid is padded with zeros. on output
only the non-padded part of the grid is written.""")
        iwndow = geomodule.NewEntry(frame, root, 8,
        geomodule.textwidth, "normal", "Tapering window width:", "5",
        "Width of cosine-tapered window zone in grid points")

        resultfile = geomodule.FileSelectorSave(frame, root, 95,
        geomodule.textwidth, "normal", "Name of file to hold predictions:" ,
        "result.gri", """
Here are output on grid format is written. For point prediction
both predicted values and error estimates are written here. """)

        geomodule.seperator_line(frame, 20, "Output parameters")

######################################################
##The Bottom line
######################################################

        geomodule.seperator_line(frame,geomodule.maxrow-2, "Running options. Working in "+geomodule.getcwd)
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
