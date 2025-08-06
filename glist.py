#!/usr/bin/python
# $Id: glist.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for GCOMB
Thomas R. N. Jansson,
tjansson@fys.ku.dk
"""

#Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *    
    # This is the new name in python 3.

import tkFileDialog
import tkFont

import os
import geomodule

## Constants
programname = "GLIST - Converts grid file to list file"
version = "0.1"
jobfile = "glist.inp"
logfile = "glist.log"
execute = "glist"
description = """GLIST - Converts grid file to list file

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
                resultfile.textentry.get(), " \n"]

                #################################
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
        inputfile = geomodule.FileSelector(frame, root, 71,
        geomodule.textwidth, "normal",
        "Name of grid file:" , "data/geoid.gri", "nohelp")
        resultfile  = geomodule.FileSelectorSave(frame, root, 72,
        geomodule.textwidth, "normal",
        "Name of list file:" , "geoid.dat", "nohelp")

        geomodule.seperator_line(frame,90,
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
