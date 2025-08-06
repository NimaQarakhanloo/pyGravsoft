#!/usr/bin/python
# $Id: covfft.py 286 2009-07-26 15:36:32Z tjansson $
# vim: tabstop=4 expandtab shiftwidth=4
"""Estimation of 2D covariance functions using FFT,
Thomas R. N. Jansson,
tjansson@fys.ku.dk
"""

##Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *
    # This is the new name in python 3.


import os
from SimpleDialog import SimpleDialog
import geomodule

## Constants
programname = "COVFFT - Estimation of 2D covariance functions using FFT"
version = "0.1"
jobfile = "covfft.inp"
logfile = "covfft.log"
execute = "covfft"
description = """COVFFT

PROGRAM FOR COMPUTING COVARIANCE FUNCTIONS, POWER SPECTRA,
AND POTENTIAL DEGREE VARIANCES FOR GRIDDED GRAVITY DATA
BY FFT. POWER SPECTRUM WILL BE GIVEN IN UNITS OF
MGAL**2*DEGREE**2, AND EXPRESSED IN DB (10LOG10).
THE PROGRAM ALSO OUTPUTS ANISOTROPY INDEX AND SLOPE OF
OF POWER SPECTRUM

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
                wordlist = [inputfile.textentry.get(),"\n",
                resultfile.textentry.get(), "\n",
                gridfile.textentry.get(), "\n",
                rfila.textentry.get(), "\n",
                inne.textentry.get(), "\n",
                iwndow.textentry.get()]
                file.writelines(wordlist)
                file.close

                statusbar.config("Data writtten to "+jobfile+" succesfully",
                                 "normal")
            else:
                statusbar.config("ERROR writing to "+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphical section
######################################################

#        mainLabel = geomodule.mainline(frame, programname, version)

        inputfile = geomodule.FileSelector(frame,root,10,
        geomodule.textwidth, "normal",
        "Name of grid file: ",
        "grav.gri", "nohelp")

        rfila = geomodule.NewEntry(frame,root,11,
        geomodule.textwidth, "normal",
        "SW corner of wanted subgrid (deg): ",
        "56.0 10.0", "nohelp")

        inne = geomodule.NewEntry(frame,root,12,
        geomodule.textwidth, "normal",
        "Number of points in subgrid (even)",
        "10 10", "nohelp")

        iwndow = geomodule.NewEntry(frame,root,13,
        geomodule.textwidth, "normal",
        "Width of cosine-tapered windows zone in grid points",
        "2", "nohelp")

        resultfile = geomodule.FileSelectorSave(frame,root,14,
        geomodule.textwidth, "normal",
        "Name of result file:",
        "covariance.txt", "nohelp")

        gridfile = geomodule.FileSelector(frame,root,15,
        geomodule.textwidth, "normal",
        "Name of 2D function gridfile:",
        "cov.gri", "nohelp")

        geomodule.seperator_line(frame, 60, "Running options. Working in "+geomodule.getcwd)
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
