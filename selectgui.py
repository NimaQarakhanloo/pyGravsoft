#!/usr/bin/python
# $Id: selectgui.py 286 2009-07-26 15:36:32Z tjansson $
"""Select GUI,
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
programname = "SELECT - Thin and/or average data"
version = "0.2"
jobfile = "select.inp"
logfile = "select.log"
execute = "select"
description = """SELECT

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
                mode.textentry.get()," ",iang.textentry.get()," ",
                datanumber.textentry.get(),"\n",
                pixeldef.textentry.get(), "\n"]

                if rejection.textentry.get() != "0":
                    wordlist.extend([rejection.textentry.get(),"\n",])

                if winspec.textentry.get() != "0.0 0.0 0.0 0.0":
                    wordlist.extend([winspec.textentry.get(),"\n",])

                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+
                " succesfully","normal")
            else:
                statusbar.config("ERROR writing to "+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile,logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The first section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)

        inputfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Input data file:",
        "data.dat",
        "Data may either be a point file or a grid file in free format.")
        mode = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Operation mode:", "1","""
0: reformat all data
1: select data closest to knots (if dfi=0 all data in area selected)
2: average all data in grid cell and output in list format
3: do, output in grid format with 9999.0 in cells without data
4: primitive screen plot
5: as 1, but data with estimated error larger than 'rejlev' are rejected
6: as 1, but a window inside the area is excluded
7: as 1, but select minimal value in cell
8: as 1, but with random, normal distributed noise added. noise standard
   deviation must be input after grid specification""")
        iang = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal", "Code for coordinates and format:",
        "1", """
1: no, lat/lon in degrees, height, data1, data2 ...
2: do, with lat/lon in degrees and minutes
3: do, with lat/lon in deg, min, sec
4: altimetry format (n=1 gives ssh only, 2 ssh+error, 3 ssh+error+t),
5: binary format,
6: gravity data in kms 80-char format
   (n=1 fa, n=2 fa+err, n=3 g,fa,ba, n=-3: fa,err,tc)
7: grid format (if n=0 a height grid of integers assumed)""")
        datanumber = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal", "Data cloumn number:", "1", """
This is the number of data following no, lat, lon,
elev (max 3) in standard format, or number of
wanted data in data line format. For averaging
only data number 'n' is used and output.""")

        geomodule.seperator_line(frame, 10, "Pixel definition (mode > 0)")
        pixeldef = geomodule.NewEntry(frame, root, 11,
        geomodule.textwidth, "normal", "Pixel definition:",
        "1.0 5.0 2.0 6.0 0.1 0.1", """
(south boundary, north boundary, east boundary, west
boundary, pixel distance northern, pixel distance
easteren, degrees or m) (iell, izone - only when
limits are northing and easting in meter)""")

        geomodule.seperator_line(frame, 20,
        "Rejection level (mode 5 and 7 only)")
        rejection = geomodule.NewEntry(frame, root, 21,
        geomodule.textwidth, "normal", "Rejection level:",
        "0", "0 means no rejection.")

        geomodule.seperator_line(frame, 30,
        "Window specification (mode 6 and 7)")
        winspec = geomodule.NewEntry(frame, root, 31,
        geomodule.textwidth, "normal",
        "Windows specification:", "0.0 0.0 0.0 0.0", """
Window definition lat1,lat2, lon1,lon2. (0.0 0.0 0.0 0.0) for no window.""")

        geomodule.seperator_line(frame, 90,"Running options. Working in "+geomodule.getcwd)
        resultfile = geomodule.FileSelectorSave(frame, root, 93,
        geomodule.textwidth, "normal", "Name of file to hold ouput:",
        "result.gri", """
Here are output on grid format is written. For point
prediction both predicted values and error estimates
are written here.""")

######################################################
##The Bottom line
######################################################

        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        geomodule.quitwriterun(frame, geomodule.maxrow,
        lambda:[frame.quit()],
        lambda:[write_to_file()],
        lambda:[run_program()],
        description,root)

#        geomodule.menubar(frame,
#        lambda:[frame.quit()],
#        lambda:[write_to_file()],
#        lambda:[run_program()],
#        description,
#        root)

######################################################
## Initiate the program and start program loop
######################################################

root = Tk()
app = App(root)
root.title(programname)
root.mainloop()
