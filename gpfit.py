#!/usr/bin/python
# $Id: gpfit.py 222 2008-09-15 11:59:34Z tjansson $
""" Fitting flat-earth covariance function to gravity data,
Thomas R. N. Jansson,
tjansson@fys.ku.dk"""

##Load graphical tool kit
from Tkinter import *
import os
from SimpleDialog import SimpleDialog
import geomodule

## Constants
programname = "GPFIT - Fitting flat-earth covariance function to gravity data"
version = "0.1"
jobfile = "gpfit.inp"
logfile	= "gpfit.log"
execute	= "gpfit"
description	= """ GPFIT

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
                wordlist = [inputfile.textentry.get(), "\n", \
                sampleintervalsize.textentry.get(), " ", \
                maximumsampledistance.textentry.get(), " ",
				d12.textentry.get(), " ",
				t12.textentry.get(), "\n"]

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
## The graphical section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)
        inputfile = geomodule.FileSelector(frame, root, 2, geomodule.textwidth,
        "normal", "Gravity file:", "gravity.dat",
        "Data must be in first column after height.")
        sampleintervalsize = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Sample intervalsize [km]:",
        "10.0", "nohelp")
        maximumsampledistance = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal", "Maximal sampling distance [km]:",
        "50.0", "nohelp")

        geomodule.seperator_line(frame,10, "Configure parameters")
        d12 = geomodule.NewEntry(frame, root, 11,
        geomodule.textwidth, "normal",
        "Search range of d parameter [km] : ", "1.0 10.0",
        "d is a variable in logarithmic covariance function.")
        t12 = geomodule.NewEntry(frame, root, 12,
        geomodule.textwidth, "normal",
        "Search range of t parameter [km]: ", "10.0 100.0",
        "t is a variable in logarithmic covariance function.")

        geomodule.seperator_line(frame,80,
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
