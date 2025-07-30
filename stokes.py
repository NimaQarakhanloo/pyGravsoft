#!/usr/bin/python
# $Id: stokes.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for Stokes
Thomas R. N. Jansson
tjansson@fys.ku.dk
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
programname = "STOKES - Space domain integration for geoid or deflection of the vertical"
version = "0.4"
jobfile = "stokes.inp"
logfile = "stokes.log"
execute = "stokes"
description = """STOKES
program for integration of gravity data to geoid or deflectios
of the vertical. gravity data must be given in gridded format,
e.g. produced by 'geogrid'. the integrals are evaluated at given
station locations. in a 3 x 3 innerzone around each computation point
a local spline densification of the given data is made, analogous
to the terrain effect integration programme 'tc'.
gravity grid values 9999 or larger signals unknown values. such
cells are not taking part in the integration, but 9999's must not be
found in the innerzone (7 x 7 subgrid around computation points).

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
                gridfile.textentry.get(), " \n",
                stationfile.textentry.get(), " \n",
                resultfile.textentry.get(), " \n",
                mode.textentry.get()," ", LMEAN.state.get()," ",
                psi12.textentry.get(),"\n"]

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
        LMEAN = geomodule.NewRadioButton(frame, root,"Remove mean?","select",
        "nocommand", "", "nohelp")

        #mainLabel = geomodule.mainline(frame, programname, version)
        gridfile = geomodule.FileSelector(frame, root, 1,
        geomodule.textwidth, "normal",
        "Grid file:", "gravi.gri", "Grid with gravity anomalies")
        stationfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Station file:",
        "stations.dat", "Location of computation points.")
        mode = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Operation mode:",
        "1", "1: Stokes integration.\n2: Vening-Meinesz")
        LMEAN.draw(4)
        psi12 = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal", "Cap size range [degrees]:",
        "0 9.0", "nohelp")

        geomodule.seperator_line(frame,30,
                                 "Running options. Working in "
                                 +geomodule.getcwd)
        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        resultfile = geomodule.FileSelectorSave(frame, root, 31,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "result.dat", "nohelp")
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
