#!/usr/bin/python
# $Id: fitgeoid.py 286 2009-07-26 15:36:32Z tjansson $
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
programname = "FITGEOID - Fit surface to GPS levelling"
version = "0.1"
jobfile = "fitgeoid.inp"
logfile = "fitgeoid.log"
execute = "fitgeoid"
description = """FITGEOID

This program calls geogrid and geoip to grid differences between
a gravimetric geoid and a set of gps geoid observations.

For further consult the source code.
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
                wordlist = [geoidfile.textentry.get(), "\n",
                gpsfile.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                nqmax.textentry.get()," ",
                itrend.textentry.get()," ",
                ipred.textentry.get(),"\n",
                predpar.textentry.get(), "\n",
                griddef.textentry.get(),"\n"]

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

       # mainLabel = geomodule.mainline(frame, programname, version)

        geoidfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal",
        "Input geoid grid filename:", "geoid.gri", """
The file contains data points in standard grid-format, i.e.
given as lines 'stat-id, lat, lon, height, data(1), .. ,
data(nd)', where 'idno' specifies which data to be used.
data values <= -9999 or >= 9999 signals missing data and
not used.""")
        gpsfile = geomodule.FileSelector(frame, root, 3,
        geomodule.textwidth, "normal",
        "Input GPS filename:", "gps.dat", "Standard gravsoft format.")

        resultfile = geomodule.FileSelectorSave(frame, root, 4,
        geomodule.textwidth, "normal",
        "File to hold final geoid values:" , "result.gri", """
Here are output on grid format written.""")

        griddef = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal",
        "Input grid definition:", "0.0 8.0 90.0 98.0 0.1 0.1", """
Min and max lattiude, min and max logntitude, lattiude and longtitude spacing""")


        geomodule.seperator_line(frame, 10)

        nqmax = geomodule.NewEntry(frame, root, 11,
        geomodule.textwidth, "normal",
        "Number of closets points:", "10", """
Number of closets points per quadrant used in predictions.
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
