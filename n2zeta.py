#!/usr/bin/python
# $Id: n2zeta.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for N2ZETA
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
programname = "N2ZETA - Transformation of geoide heights to height anomalies"
version = "0.1"
jobfile = "n2zeta.inp"
logfile = "n2zeta.log"
execute = "n2zeta"
description = """N2ZETA

The program computes the difference between the orthometric
height and the normal height as being equal to the Bouger-anomaly
multiplied by the altitude and divided by normal gravity.
a file with geoid heights may be converted to heightanomalies
(or the reverse) using the program.
For the conversion a fixed density of 2.67 g/cm**3 is used,
corresponding to a Bouger factor of 0.1119 mgal/m.

Input files:
file with free-air gravity anomalies of the points to be corrected.
if corrections are to be applied, file with geoid heights or height
anomalies (m).

Output file:
file to hold corrections.
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
                "t\n",
                gravityfile.textentry.get(), " \n",
                outputfile.textentry.get(), " \n",
                ndat1.textentry.get(), " \n",
                "t\n",
                inputfile.textentry.get(), " \n",
                ndat2.textentry.get(), " \n"]

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
        geomodule.textwidth, "normal", "Data file:",
        "data/nmgpslev.dat", "File with geoid heights.")
        gravityfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Gravity anomaly file:",
        "data/nmgpslev_grav.dat", "File with gravity anomalies in same points as data files.")
        ndat1 = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal",
        "Data column number of geoid heights:", "1", "nohelp")
        ndat2 = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal",
        "Data column number of gravity anomalies:", "3", "nohelp")

        geomodule.seperator_line(frame,30,
                                 "Running options. Working in "
                                 +geomodule.getcwd)
        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        outputfile = geomodule.FileSelectorSave(frame, root, 31,
        geomodule.textwidth, "normal", "Height anomaly file:",
        "zeta.dat", "File with height anomalies ")
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
