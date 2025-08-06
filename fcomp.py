#!/usr/bin/python
# $Id: fcomp.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for fcomp
Thomas R. N. Jansson
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
programname = "FCOMP - File Comparison"
version = "0.5"
jobfile = "fcomp.inp"
logfile = "fcomp.log"
execute = "fcomp"
description = """ FCOMP

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
                inputfile1.textentry.get(), " \n",
                inputfile2.textentry.get(), " \n",
                resultfile.textentry.get(), " \n"]
                if operation.list.get(ACTIVE) == "Addition":
                    wordlist.extend([ "2 "])
                if operation.list.get(ACTIVE) == "Subtraction":
                    wordlist.extend([ "1 "])
                wordlist.extend([numberofcolumns.textentry.get(), " ",
                binsize.textentry.get(), "\n"])
                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+
                jobfile+" succesfully","normal")
            else:
                statusbar.config("ERROR writing to "+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphicalfirst section
######################################################

       # mainLabel = geomodule.mainline(frame, programname, version)
        inputfile1 = geomodule.FileSelector(frame, root, 11,
        geomodule.textwidth, "normal",
        "Name of input file 1:",
        "file1.dat", "nohelp")
        inputfile2 = geomodule.FileSelector(frame, root, 12,
        geomodule.textwidth, "normal",
        "Name of input file 2:",
        "file2.dat", "nohelp")
        numberofcolumns = geomodule.NewEntry(frame, root, 13,
        geomodule.textwidth, "normal",
        "Number of columns after alt. column:",
        "2", "nohelp")
        binsize = geomodule.NewEntry(frame, root, 14,
        geomodule.textwidth, "normal",
        "Histogram binsize:",
        "1.0", "nohelp")
        # When editing the list in operation be carefull to
        # change in the write_to_file
        operation = geomodule.liste(frame, 15,
        "Column operation", ["Addition", "Subtraction"])
        resultfile = geomodule.FileSelectorSave(frame, root, 20,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "result.dat", "nohelp")

        geomodule.seperator_line(frame,70, "Running options. Working in "+geomodule.getcwd)
        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        geomodule.quitwriterun(frame, geomodule.maxrow,
        lambda:[frame.quit()],
        lambda:[write_to_file()],
        lambda:[run_program()],description, root)

######################################################
## Initiate the program and start program loop
######################################################

root = Tk()
app = App(root)
root.title(programname)
root.mainloop()
