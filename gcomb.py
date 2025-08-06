#!/usr/bin/python
# $Id: gcomb.py 286 2009-07-26 15:36:32Z tjansson $
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

import os
import geomodule

import tkFont
import tkFileDialog

## Constants
programname = "GCOMB - Combining two grids"
version = "0.5"
jobfile = "gcomb.inp"
logfile = "gcomb.log"
execute = "gcomb"
description = """GCOMB

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
                resultfile.textentry.get(), " \n",
                mode.textentry.get()," ", iout.textentry.get()," \n"]

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
        inputfile1 = geomodule.FileSelector(frame, root, 71,
        geomodule.textwidth, "normal",
        "Name of grid file 1:" , "file1.gri", "nohelp")
        inputfile2 = geomodule.FileSelector(frame, root, 72,
        geomodule.textwidth, "normal",
        "Name of grid file 2:" , "file2.gri", "nohelp")
        mode = geomodule.NewEntry(frame, root, 73,
        geomodule.textwidth, "normal",
        "Write grid combinationmode: " , "2", """
mode = 1: subtract 'grid1' minus 'grid2'
mode = 2: add 'grid1' plus 'grid2'
mode = 5: grid overwrite values in 'grid2' is written on top of 'grid1',
except when 9999 is encountered, then the grid1-values are kept.
mode = 6: grid select values in 'grid1' are written only when there
is no data in 'grid2'.
mode = 9: Grid 1 is a Bouguer anomaly grid, grid 2 a height grid.
Output N - zeta (difference geoid minus quasigeoid) in linear approximation.
If grid 1 is a free-air grid, the difference zeta* - zeta is obtained.
        """)
        iout = geomodule.NewEntry(frame, root, 74,
        geomodule.textwidth, "normal",
        "Output file type:" , "2", """
1 = Binary form.
2 = Txt-file with reals.
3 = Txt-file with integers. """)

        resultfile  = geomodule.FileSelectorSave(frame, root, 95,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "result.gri", "nohelp")

        geomodule.seperator_line(frame,90,"Running options. Working n "
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
