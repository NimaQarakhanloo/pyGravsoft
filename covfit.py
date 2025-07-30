#!/usr/bin/python
# $Id: covfit.py 286 2009-07-26 15:36:32Z tjansson $
# vim: tabstop=4 expandtab shiftwidth=4
"""Covariance fitting GUI,
Thomas R. N. Jansson,
tjansson@fys.ku.dk
"""

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
programname = "COVFIT - Covariance fitting"
version = "0.5"
jobfile = "covfit.inp"
logfile = "covfit.log"
execute = "covfit16"
description = """COVFIT

http://cct.gfy.ku.dk/covfit16v1.pdf
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
                wordlist = ["t\n",
                "4 f f t f \n",
                "2\n",
                "4\n",
                covariancemodel.textentry.get(), " f\n",
                "-1 2 ", scale.textentry.get(), " f\n",
                errordegreefile.textentry.get(), "\n",
                numberofiterations.textentry.get(), " ",
                aprioriweights.textentry.get(), "\n",
                "1\n",
                numberofvalues.textentry.get(), " ",
                datatype.textentry.get()," ",
                datatype.textentry.get(), " ",
                meanaltitude.textentry.get(), " ",
                meanaltitude.textentry.get(), " 1 1.0 t\n",
                datavariance.textentry.get(), " ",
                datavariance.textentry.get(), " ",
                areapecification.textentry.get(), "\n",
                inputfile.textentry.get(), "\n"]
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

       # mainLabel = geomodule.mainline(frame, programname, version)
        inputfile = geomodule.FileSelector(frame,root,10,
        geomodule.textwidth, "normal",
        "Name of file with empirical covariances: ",
        "table.txt", "nohelp")
        numberofvalues = geomodule.NewEntry(frame,root,11,
        geomodule.textwidth, "normal",
        "Input number of values in table: ", "10", "nohelp")
        datatype = geomodule.NewEntry(frame,root,12,
        geomodule.textwidth, "normal",
        "Input code for observations: ", "3",
        """
        1 - HEIGHT-ANOMALY OR GEOID UNDULATION
        2 - GRAVITY DISTURBANCE
        3 - GRAVITY ANOMALY
        5 - VERTICAL GRAVITY DISTURBANCE GRADIENT""" )
        meanaltitude = geomodule.NewEntry(frame,root,13,
        geomodule.textwidth, "normal",
        "Input the mean altitude (m): ", "0.0", """
This value may be obtained using TC on a digital terrain model.
        """)
        datavariance = geomodule.NewEntry(frame,root,14,
        geomodule.textwidth, "normal",
        "Input data variance at mean altitude: ", "400", "nohelp")
        areapecification = geomodule.NewEntry(frame,root,15,
        geomodule.textwidth, "normal",
        "Input data area specification :", "54.5 57.5 0.1 7.0 13.0 0.2","""
        Min, Max latitude, mean lat.spacing,
        Min, Max longtitude, mean long.spacing
        The values must reflect the area-size and mean 
        data-spacing of the data used to estimate the 
        empirical covariance function.
        """)

        geomodule.seperator_line(frame, 20, "Model parameters")
        covariancemodel = geomodule.NewEntry(frame,root,21,
        geomodule.textwidth, "normal",
        "Input convariance model parameters: ", "-1.0 400.0 360",
        """
        (-1.0) THE DEPTH TO THE BJERHAMMAR SPHERE IN KM.
        (400) VARIANCE OF GRAVITY ANOMALIES AT ZERO ALTITUDE.
        (360) MAXIMAL DEGREE FOR EMPIRICAL DEGREE-VARIANCES. """ )
        scale = geomodule.NewEntry(frame,root,22,
        geomodule.textwidth,"normal",
        "Input error degree variance scale factor :", "1.0", """
The resulting value in the output AA must always be
positive. If not there is an error in other input parameters.
        """)
        errordegreefile	= geomodule.FileSelector(frame,root,23,
        geomodule.textwidth,"normal",
        "Input name of error degree variance file: ",
        "data/egm96.edg", "nohelp")

        geomodule.seperator_line(frame, 40, "Itteration paramters")
        numberofiterations = geomodule.NewEntry(frame,root,41,
        geomodule.textwidth,"normal",
        "Input number of iterations: ", "10",
        "Iterations for nonlinar adjustment.")
        aprioriweights = geomodule.NewEntry(frame,root,42,
        geomodule.textwidth,"normal",
        "Input three weights: ", "1.0 1.0 1.0",
        "Weights for nonlinar adjustment.")

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
