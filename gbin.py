#!/usr/bin/python
# $Id: gbin.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for gbin
Thomas R. N. Jansson
tjansson@fys.ku.dk
"""

#Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *    
    # This is the new name in python 3.

import os
import geomodule

import tkFont
import tkFileDialog

## Constants
programname = "GBIN - Convert grid to binary or reverse"
version = "0.5"
jobfile = "gbin.inp"
logfile = "gbin.log"
execute = "gbin"
description = """GBIN

Gravsoft manual can be found in doc/gravsoft-manual.pdf
"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

###################
## Functions
##################

        def write_to_file():
            file = open(jobfile, 'w')
            if file:
                wordlist = [\
                inputfile.textentry.get(), " \n",
                resultfile.textentry.get(), " \n"]

                if choices.list.get(ACTIVE) == "ASCII to binary":
                    wordlist.extend([ "1 "])
                if choices.list.get(ACTIVE) == "Binary to ASCII":
                    wordlist.extend([ "2 "])
                if choices.list.get(ACTIVE) == "Binary to ASCII integers":
                    wordlist.extend([ "3 "])

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
        # When editing the list in choices be carefull to
        # change in the write_to_file
        choices = geomodule.liste(frame, 70, "Select conversion:",
        ["ASCII to binary","Binary to ASCII", "Binary to ASCII integers"])
        inputfile = geomodule.FileSelector(frame, root, 71,
        geomodule.textwidth, "normal",
        "Name of file to be converted:", "inputfile", "nohelp")
        resultfile = geomodule.FileSelectorSave(frame, root, 72,
        geomodule.textwidth, "normal",
        "Name of converted file:", "outputfile","nohelp")

        geomodule.seperator_line(frame,90, "Running options. Working in "+geomodule.getcwd)
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
