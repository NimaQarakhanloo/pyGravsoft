#!/usr/bin/python
# $Id: empcov.py 286 2009-07-26 15:36:32Z tjansson $
"""Emperical Covariance Estimation,
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
programname = "EMPCOV - Emperical Covariance Estimation"
version = "0.5"
jobfile = "empcov.inp"
logfile	= "empcov.log"
execute	= "empcov"
description	= """ EMPCOV

http://cct.gfy.ku.dk/empcov.pdf
"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## Initiate all the button objects - they will first be drawn when
## foo.draw(linenumber) is called.
######################################################
        LMEAN	= geomodule.NewRadioButton(frame, root,
        "Should mean value be subtracted: " , "unselect",
        "nocommand", "", "nohelp")
        LAREA	=	geomodule.NewRadioButton(frame, root,
        "Should data in subarea be used: " ,"unselect",
        lambda:[areaspecification.enable()],
        lambda:[areaspecification.disable()], "nohelp" )

######################################################
## Functions -- these needs to be defined before they are used
######################################################

        def write_to_file():
            file = open(jobfile, 'w')
            if file:
                wordlist = ["Emperical Covarince Estimation\n", \
                sampleintervalsize.textentry.get(), " ",
                numberofsamplinginterval.textentry.get(), " 3 F T ",
                LMEAN.state.get(), " \n", resultfile.textentry.get(), "\n", \
                "100000 9 T 3 1 0 ", histogrambinsize.textentry.get(), " ",
                LAREA.state.get(), " \n", \
                str(int(position.textentry.get())+1), " ", str(int(position.textentry.get())+1), " 0 ", "\n",
                inputfile.textentry.get(), "\n"]
                if LAREA.state.get() == "t":
                    wordlist.extend([areaspecification.textentry.get(), "\n"])
                wordlist.extend(["T\n"])

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
        inputfile = geomodule.FileSelector(frame, root, 2, geomodule.textwidth,
        "normal", "Input data filename:", "observations.dat", "nohelp")
        position = geomodule.NewEntry(frame, root, 3, geomodule.textwidth,
        "normal", "Input position of data element: :", "1",
        "Column number after  altitude, if altitude is used input 0.")
        sampleintervalsize = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal", "Input sample intervalsize (arcmin):",
        "10.0", "nohelp")
        numberofsamplinginterval = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal", "Input number of sampling intervals:",
        "20", "nohelp")

        geomodule.seperator_line(frame,10, "Configure parameters")
        LMEAN.draw(11)
        LAREA.draw(12)
        areaspecification = geomodule.NewEntry(frame, root, 67,
        geomodule.textwidth, "normal",
        "     Input area boundaries: ", "54.5 57.5 7.0 13.0",
        "Min, Max latitude, Min, Max longtitude")
        areaspecification.disable()
        histogrambinsize = geomodule.NewEntry(frame, root, 74,
        geomodule.textwidth, "normal", 	"Input histogram bin size ",
        "5.0", "nohelp")
        resultfile = geomodule.FileSelectorSave(frame, root, 79,
        geomodule.textwidth, "normal", "Name of file to hold result:",
        "table.txt", "nohelp")

        geomodule.seperator_line(frame,80,
                                 "Running options. Working in "+geomodule.getcwd)
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
