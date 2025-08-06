#!/usr/bin/python
# $Id: tcgrid.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for tcgrid
Thomas R. N. Jansson
tjansson@fys.ku.dk
minor addition to help 2009-01-25 by cct.
"""

#Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *    
    # This is the new name in python 3.

import tkFont
import tkFileDialog

import os
import geomodule

## Constants
programname = "TCGRID - DTM grids and mean terrain surfaces for RTM method"
version = "0.2"
jobfile = "tcgrid.inp"
logfile = "tcgrid.log"
execute = "tcgrid"
description = """TCGRID

Program to produce a mean  elevation grid of a file containing
a digital terrain model (standard format or ngs modified).
The mean grid is obtained by simple averaging. areas with no
elevations are given value 9999. if an individual grid
element is covered partly by elevations, only these will be
averaged. if no data is available the unknown flag 9999 will be written.

If wanted, the mean elevation grid may be low-pass filtered
using a moving-average window of 'iffi' x 'ifla' cells.
if 'ifla' = -1 the longitude resolution is changed to match
the latitude variation.

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
                inputfile.textentry.get(), " \n",
                resultfile.textentry.get(), " \n",
                itype.textentry.get()," ", area.textentry.get()," ",
                idfila.textentry.get()," ", iffila.textentry.get(), " \n"]
                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+
                " succesfully","normal")
            else:
                statusbar.config("ERROR writing to"+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphicalfirst section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)
        inputfile = geomodule.FileSelector(frame, root, 1,
        geomodule.textwidth, "normal",
        "Name of grid file:", "dtm.gri", "nohelp")
        itype = geomodule.NewEntry(frame, root, 2,
        geomodule.textwidth, "normal",
        "Specification of grid file format:", "0",
        "0,1 standard GRAVSOFT format. 2 NGS")
        area = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal",
        "Wanted area boundaries:", "54.0 58.0 10.0 20.0",
        "Area boundaries in degrees. (0 0 0 0 indicates whole area)")
        idfila = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal",
        "Grid cell averages:", "5 5", "Integers only.")
        iffila = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal",
        "Moving average window size:", "3 3",
        "Odd numbers only. (= 0 or 1 meand no low pass filtering).")

        resultfile = geomodule.FileSelectorSave(frame, root, 20,
        geomodule.textwidth, "normal",
        "Name of file to hold result:", "result.dat", "nohelp")

        geomodule.seperator_line(frame,70,
                                 "Running options. Working in "
                                 +geomodule.getcwd)
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
