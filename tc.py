#!/usr/bin/python
# $Id: tc.py 286 2009-07-26 15:36:32Z tjansson $
"""TC GUI,
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
programname = "TC - Compute terain effect on gravimetric quantities"
version = "0.3"
jobfile = "tc.inp"
logfile = "tc.log"
execute = "tc"
description = """TC

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
                statfile.textentry.get(), "\n",
                dtmfile1.textentry.get(), "\n",
                dtmfile2.textentry.get(), "\n",
                dtmfile3.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                itype.textentry.get()," ",ikind.textentry.get()," ",
                izcode.textentry.get()," ", istype.textentry.get(), " ",
                rho.textentry.get(),"\n",
                winspec.textentry.get(), "\n",
                r1.textentry.get()," ",r2.textentry.get(), "\n",
                nd.textentry.get(), "\n"]

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
        statfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal",
        "Station list file:","stations.dat", "nohelp")
        dtmfile1 = geomodule.FileSelector(frame, root, 3,
        geomodule.textwidth, "normal",
        "Detailed elevation grid file:", "dtmfile1.gri", "nohelp")
        dtmfile2 = geomodule.FileSelector(frame, root, 4,
        geomodule.textwidth, "normal",
        "Coarse elevation grid file:", "dtmfile2.gri", "nohelp")
        dtmfile3 = geomodule.FileSelector(frame, root, 5,
        geomodule.textwidth, "normal",
        "Reference elevation grid file:", "dtmfile3.gri", "nohelp")

        itype = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal",
        "Data type:", "5", """
1  gravity disturbance (dg) - mgal
2  deflections (ksi, eta) - arcsec
3  height anomaly (ha) - meter
4  dg, ksi, eta
5  gravity anomaly (dg - ind.eff.)
6  anomaly, ksi, eta, ha
7  tzz - vertical gravity gradient (z positive up)
8  txx, tyy, tzz
9  all gradients.""")
        ikind = geomodule.NewEntry(frame, root, 7,
        geomodule.textwidth, "normal",
        "Type of effect:", "4", """
1  topographic effect
2  isostatic effec
3  terrain corrections
4  residual terrain effetcs
5  do, using precomputed terrain corrections to dist 'rtc' km
   (if tc-values are missing, ordinary rtm effects are computed)""")
        izcode = geomodule.NewEntry(frame, root, 8,
        geomodule.textwidth, "normal",
        "Placement of station:", "1", """
0  station on terrain, change station elevation
1          do        , change terrain model
2          do        , change terrain model in land points only
3  station free
4  station free, no spline densification
   (dma special - if topography is above computation
   level, the terrain effects will be computed at the
   topography level in both mode 3 and 4)""")
        istype = geomodule.NewEntry(frame, root, 9,
        geomodule.textwidth, "normal",
        "Type of operation:", "1", """
0  no statfile, compute in grid
1  compute effects in statfile
2  add effects to value in statfile
3  subtract effects
4  statfile with 80-char KMS gravity recs (80-char output
   for ikind=3, otherwise normal output format) neg  as 1,2 or
   3 but UTM data in statfile.""")
        nd = geomodule.NewEntry(frame, root, 10,
        geomodule.textwidth, "normal",
        "   Data column (operation 2 or 3):", "2",
        "If type of operation equal 2 or 3 the data columnnumber must be input.")
        rho = geomodule.NewEntry(frame, root, 11, geomodule.textwidth,
        "normal", "Density:", "2.67", "Density in grams per cubic centimeter.")
        winspec = geomodule.NewEntry(frame, root, 12, geomodule.textwidth,
        "normal", "Maximum window:",
        "40.0 60.0 0.0 10.0", """
(minimum, maximal northern coordinate, minimum,
maximum eastern coordinate defining fixed maximum
area for which terrain effect is computed,
irrespectively of the specified computational
radii. unit: deg (utm: m).""")
        r1 = geomodule.NewEntry(frame, root, 13,
        geomodule.textwidth, "normal",
        "Minimum computation distance of inner grid:", "40.0", """
Minimum computation distance (r1) of inner grid (km).
The inner grid is actually used in the smallest
'subsquare' of the coarse grid, covering a circle
of radius r1. """)
        r2 = geomodule.NewEntry(frame, root, 14,
        geomodule.textwidth, "normal",
        "Maximal radius of computation", "100.0", """
Maximal radius, r2, of computation (km). if r2 = 0 no outer
grid is used. if a fixed-area computation is wanted 'r2'
must be sufficiently large.""")

        geomodule.seperator_line(frame, 90, "Running options. Working in "+geomodule.getcwd)
        resultfile = geomodule.FileSelectorSave(frame, root, 93,
        geomodule.textwidth, "normal",
        "Name of file to hold ouput:" , "result.dat", "nohelp")

######################################################
##The Bottom line
######################################################

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
