#!/usr/bin/python
# $Id: geomain.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for geomain
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
programname = "GEOMAIN - 2 points: Distance and azimuth or reverse"
version = "0.2"
jobfile = "geomain.inp"
logfile = "geomain.log"
execute = "geomain2"
description = """GEOMAIN

THE PROGRAM SOLVES THE 2 MAIN GEODETIC PROBLEMS FOR THE ELLIPSOID.

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
                "t", "\n",
                ispher.textentry.get(), " \n",
                "3", " \n",
                "f", " \n",
                inputfile.textentry.get(), "\n",
                "t", "\n",
                resultfile.textentry.get(), "\n",
                LDIREC.state.get()]

                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+" succesfully",
                                     "normal")
            else:
                statusbar.config("ERROR writing to"+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphicalfirst section
######################################################
        LDIREC = geomodule.NewRadioButton(frame, root,"Direct solution?",
        "select", "nocommand", "", """
From latitude, longitude, azimuth and distance to another
point. If NO from latitude and longtitude of two points
to azimuth and distance. All angles in degrees and
distances in meters.""")

        #mainLabel = geomodule.mainline(frame, programname, version)
        inputfile = geomodule.FileSelector(frame, root, 1,
        geomodule.textwidth, "normal", "Input file:" ,
        "data/geomain_dir.dat", """
Input either station number, latitude and longtitude of
 first and second point (see data/geomain_dir.dat)
 or stationnumber, latitude,
 longtitude, azimuth and distance (see data/geomain-rev.dat) .""")
        LDIREC.draw(2)
        ispher = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Select spheriod:" ,
        "4", """
1: SPHERE,
2: CLARKE 1866,
3: HAYFORD 1909 (INTERNATIONAL),
4: GRS 1980,
5: CLARK1880,
6: BESSEL 1841,
7: KRASOVSKY 1940,
8: WGS 1972,
9: AUSTRALIAN 1965,
10: AIRY1849,
11: EVEREST 1830,
12: HOUGH 1956,:
13: FISHER 1960,
14&15: SPHERE""")

        geomodule.seperator_line(frame,185, "Running options. Working in "+geomodule.getcwd)
        resultfile = geomodule.FileSelectorSave(frame, root, 186,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "result.dat", "nohelp")

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
