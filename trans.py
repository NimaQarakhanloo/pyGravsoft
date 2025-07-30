#!/usr/bin/python
# $Id: trans.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for trans13
Thomas R. N. Jansson
tjansson@fys.ku.dk
last change 2009-03-11 by cct.
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
programname = "TRANS - Transformation of coordinates to or from a 2D or 3D system"
version = "0.2"
jobfile = "trans.inp"
logfile = "trans.log"
execute = "trans13"
description ="""TRANS13

PROGRAM FOR THE TRANSFORMATION OF FILES OF COORDINATES TO
OR FROM GEOGRAPHICAL COORDINATES FROM/TO A PLANE SYSTEM.

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
                "Transformation input:", " \n",
                LFORWA.state.get(), " \n",
                ispher.textentry.get(), " \n",
                iproj.textentry.get(), "\n"]

                if (int(iproj.textentry.get()) == 2) or \
                (int(iproj.textentry.get()) == 4) or \
                (int(iproj.textentry.get()) == 5):
                    wordlist.extend([cesfa.textentry.get(), "\n"])
                if (int(iproj.textentry.get()) == 2) or \
                (int(iproj.textentry.get()) == 3) or \
                (int(iproj.textentry.get()) == 5):
                    wordlist.extend([cmer.textentry.get(), "\n"])
                if int(iproj.textentry.get()) == 2:
                    wordlist.extend([xconst.textentry.get(), "\n"])
                if (int(iproj.textentry.get()) == 2) or \
                (int(iproj.textentry.get()) == 3) or \
                (int(iproj.textentry.get()) == 4):
                    wordlist.extend([bpar.textentry.get(), "\n"])
                if int(iproj.textentry.get()) == 4 :
                    wordlist.extend([spar12.textentry.get(), "\n"])
                if int(iproj.textentry.get()) == 1 :
                    wordlist.extend([izone.textentry.get(), "\n"])
                if int(iproj.textentry.get()) == 9 :
                    if LSHIFT.state.get() == "t":
                        wordlist.extend(["t", "\n",
                        dxyz.textentry.get()," ","\n",
                        s1.textentry.get()," ","\n",
                        eps321.textentry.get(), "\n"])
                    else:
                        wordlist.extend(["f", "\n"])
                wordlist.extend([nd.textentry.get()," ",
                nd.textentry.get(), "\n",
                 "3", "\n",
                 "f", "\n",
                 "f", "\n",
                  inputfile.textentry.get(), "\n",
                    "t", "\n",
                    resultfile.textentry.get()])

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
## The graphicalfirst section
######################################################
        LFORWA = geomodule.NewRadioButton(frame, root,
        "Convert from geographical to 2D or 3D?", "select",
        "nocommand", "", """
If No is select the transformation is from 2D or 3D to
geographical coordinates.""")
        LSHIFT = geomodule.NewRadioButton(frame, root,
        "Is datumshift necessary?", "select",
        lambda: [dxyz.enable(), s1.enable(),eps321.enable()],
        lambda: [dxyz.disable(), s1.disable(),eps321.disable()],
        "This may be needed between different reference systems.")

        
        #mainLabel = geomodule.mainline(frame, programname, version)
        inputfile = geomodule.FileSelector(frame, root, 1,
        geomodule.textwidth, "normal", "Point datafile:",
        "data/position.dat", """
Point data on GRAVSOFT standard format (int stationnumber,
latitude,longitude [degrees], altitude [m])""")
        nd = geomodule.NewEntry(frame, root, 2,
        geomodule.textwidth, "normal",
        "Number of datacolumns in point file:" , "1",
        "The last data column is reproduced in output.")
        ispher = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal", "Select spheriod:", "4","""
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
12: HOUGH 1956,
13: FISHER 1960,
14&15: SPHERE""")
        iproj = geomodule.NewEntry(frame, root, 4,
        geomodule.textwidth, "normal",
        "Transformation type:" , "1","""
1: U T M,
2: TRANSVERSE MERCATOR,
3: MERCATOR',
4: LAMBERT CONFORMAL CONIC (2 STD. PARAL.),
5: POLAR STEREOGRAPHIC (AZIMUTHAL-),
6: SYSTEM 34 JYLLAND,
7: SYSTEM 34 SJAELLAND,
8: SYSTEM 45 BORNHOLM,
9: CARTESIAN (X,Y,Z).""")
        LFORWA.draw(5)

        geomodule.seperator_line(frame,20, "Transformation specifications")
        izone = geomodule.NewEntry(frame, root, 21,
        geomodule.textwidth, "normal",
        "[Type 1]\t\t UTM zone:", "34", "nohelp")
        xconst = geomodule.NewEntry(frame, root, 22,
        geomodule.textwidth, "normal",
        "[Type 2]\t\t Abscissa constant [m]:" , "100000", "nohelp")
        cesfa = geomodule.NewEntry(frame, root, 23,
        geomodule.textwidth, "normal",
        "[Type 2,4]\t Central scale factor:" , "1.0", "nohelp")
        cmer = geomodule.NewEntry(frame, root, 24,
        geomodule.textwidth, "normal",
        "[Type 2,3,4]\t Longitude of central meridian:",
        "15.0", "Longtitude in degrees.")
        bpar = geomodule.NewEntry(frame, root, 25,
        geomodule.textwidth, "normal",
        "[Type 2,3,4]\t Latitude of base parallel [deg]:","56.0", "nohelp")
        spar12 = geomodule.NewEntry(frame, root, 26,
        geomodule.textwidth, "normal",
        "[Type 4]\t\t Latitude of first and second parallel [deg]:",
        "54.0 58.0", "nohelp")

        geomodule.seperator_line(frame,80,
        "For transformation between geographical and 3D coordinates (Type 9)")
        LSHIFT.draw(81)
        dxyz = geomodule.NewEntry(frame, root, 82,
        geomodule.textwidth, "normal",
        "     Translation vector [m]:", "102.0 102.0 129.0", "nohelp")
        s1 = geomodule.NewEntry(frame, root, 83,
        geomodule.textwidth, "normal",
        "     Scale factor - 1.0:", "0.0000025", "nohelp")
        eps321 = geomodule.NewEntry(frame, root, 84,
        geomodule.textwidth, "normal",
        "     Rotations around x,y,z axis [arcsec]:", "0.4 -0.2 0.4", "nohelp")

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
