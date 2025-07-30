#!/usr/bin/python
# $Id: geogrid.py 286 2009-07-26 15:36:32Z tjansson $
"""Geogrid,
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
programname = "GEOGRID - Gridding or interpolation of irregualar distributed data"
version = "0.3"
jobfile = "geogrid.inp"
logfile = "geogrid.log"
execute = "geogrid"
description = """GEOGRID

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
                inputfile.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                errorfile.textentry.get(), "\n",
                nd.textentry.get()," ",idno.textentry.get()," ",
                nqmax.textentry.get()," ",itrend.textentry.get()," ",
                ipred.textentry.get(),"\n",
                predpar.textentry.get(), "\n",
                mode.textentry.get()," ", rkm.textentry.get(), "\n",
                predpoints.textentry.get(), "\n"]

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

        inputfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal",
        "Input data filename:", "observations.dat", """
The file contains data points in standard gi-format, i.e.
given as lines 'stat-id, lat, lon, height, data(1), .. ,
data(nd)', where 'idno' specifies which data to be used.
data values <= -9999 or >= 9999 signals missing data and
not used.""")
        nd = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal",
        "Number of data values:", "2", "Number of columns after altitude.")
        idno = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal",
        "Input position of data element: :", "1",
        "Column number after altitude, if altitude is used input 0.")
        resultfile = geomodule.FileSelectorSave(frame, root, 5,
        geomodule.textwidth, "normal",
        "File to hold predictions:" , "result.gri", """
Here are output on grid format is written. For point
prediction both predicted values and error estimates are written here. """)
        errorfile = geomodule.FileSelectorSave(frame, root, 6,
        geomodule.textwidth, "normal",
        "File to hold errorestimates:" , "errors.gri", """
For grid the error estimates are written here. For weighted
means interpolation it contains distances in kilometers to
the closets data points.""")

        geomodule.seperator_line(frame, 10)

        nqmax = geomodule.NewEntry(frame, root, 11,
        geomodule.textwidth, "normal",
        "Number of closest points:", "1", """
Number of closest points per quadrant used in predictions.
if less or equal to zero all points are used. """)
        itrend = geomodule.NewEntry(frame, root, 12,
        geomodule.textwidth, "normal",
        "Trend surface removal method", "0", """
0 - none (just grid data)
1 - remove mean
2 - remove linear function in x and y
3 - remove 2nd polynomium
4 - remove 3nd polynomium
5 - remove 4 parameter datumshift""")
        ipred = geomodule.NewEntry(frame, root, 13,
        geomodule.textwidth, "normal",
        "Select prediction method:", "1", """
Prediction method
0 - none (detrend only)
1 - collocation (kriging)
2 - weigthed means""")
        predpar = geomodule.NewEntry(frame, root, 14,
        geomodule.textwidth, "normal", "Prediction variables",
        "10.0 1.0", """
Followed by if collocation half corellation distance (km) and
RMS of data noise. If weighted means prediction power.

Example (10.0 1.0)
1 means collocation, with 10.0 half corellation
distance and 1.0 RMS of data noise.""")

        mode = geomodule.NewEntry(frame, root, 15,
        geomodule.textwidth, "normal",
        "Mode  number:", "1", """
1: grid in geographical coordinates
2, 3: grid in utm projection
4: profile interpolation
5: point interpolation, data points given in <efile> in
   standard gi format, i.e. as gi-id, lat, lon, height,
   data1, data2, ...
6: do, with input data assumed to be bouguer anomalies
    (grs67/2.67). output data will be gravity values,
    derived from interpolated anomalies and prediction
    point heights.
7: fill-in of unknown (9999) grid values. In this case
   <ifile> must contain a grid rather than a point list.
   Only the 9999-values are actually predicted. negative:
   Same as positive value, except only positive prediction
   values allowed (all negative values set to zero). Use
   for DTM interpolation only.""")
        predpoints = geomodule.NewEntry(frame, root, 16,
        geomodule.textwidth, "normal",
        "Specify prediction points:", "54.5 57.5 7.0 13.0 0.1 0.1", """
1: fi1, fi2, la1, la2, dfi, dla (geographic grid)
2: n1, n2, e1, e2, dn, de / ie, zone (utm grid)
3: fi1, la1, dn, de, nn, ne / ie, zone (utm grid
   with sw-corner lat/lon, spacing, and number of points)
4: fi1, la1, fi2, la2, n (profile with 'n' points)
5, 6, 7: (nothing)

The utm projections are defined by their zone no and
ellipsoid no ('ie'): ie = 1: wgs84/nad83, ie = 2:
hayf/ed50, ie = 3: clarke/nad27.""")
        rkm = geomodule.NewEntry(frame, root, 17,
        geomodule.textwidth, "normal",
        "Margin for data selection area:", "10.0", """
Margin for data selection area surrounding the wanted
grid (or retangualr envelope of the profile or wanted
points for mode > 3)""")

######################################################
##The Bottom line
######################################################

        geomodule.seperator_line(frame,geomodule.maxrow-2,
                                 "Running options. Working in "+
                                 geomodule.getcwd)
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
