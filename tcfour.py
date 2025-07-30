#!/usr/bin/python
# $Id: tcfour.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Graphical interface for TCFOUR
Thomas R. N. Jansson
tjansson@fys.ku.dk
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
programname = "TCFOUR - Terrain effect computation by FFT"
version = "0.2"
jobfile = "tcfour.inp"
logfile = "tcfour.log"
execute = "tcfour"
description = """TCFOUR

Program for fft analysis of digital terrain models.
space domain expressions of the integral kernels are
transformed rather than using the fourier domain analytical
kernels, in order to allow control over integration
radii.

A subgrid of the grid in 'dtmfile1' is analyzed, specified
by its sw-corner (fic,lac) and number of
points. both 'in' and 'ie' must be  e v e n .
a reference height grid may be subtracted (for rtm effects).

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
                dtm1file.textentry.get(), " \n",
                dtm2file.textentry.get(), " \n",
                resultfile.textentry.get(), " \n",
                mode.textentry.get()]
                if LREF.state.get()  == "t":
                    wordlist.extend([" t " ])
                else:
                    wordlist.extend([" f " ])
                wordlist.extend([dist12.textentry.get(),"\n",
                filac1.textentry.get()," ", ine1.textentry.get(),"\n"])
                if int(mode.textentry.get()) == 2:
                    wordlist.extend([filac2.textentry.get()," ",
                    ine2.textentry.get(),"\n"])

                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+
                " succesfully","normal")
            else:
                statusbar.config("ERROR writing to"+jobfile,"warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphicalfirst section
######################################################
        LREF = geomodule.NewRadioButton(frame, root,
        "Subtract reference grid?","select","nocommand", "", "nohelp")

        #mainLabel = geomodule.mainline(frame, programname, version)
        dtm1file = geomodule.FileSelector(frame, root, 1,
        geomodule.textwidth, "normal", "DTM file:", "dtm.gri", "nohelp")
        dtm2file = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Reference DTM file:",
        "dtm.ref", "nohelp")
        mode = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal",
        "Operation mode:", "4", """
0: simple filtering. 'dist1' gives the wavelength
   for cut-off in units of km.
   'dist2' > 0: low pass, < 0: high pass.
1: covariance function and power spectra,
   (2-d cov.fct. output as dtm with c(0,0)=1000)
   'dist1, dist2' is not used.
2: terrain corrections, using grid1 to distance 'dist1',
   grid2 to distance 'dist2' (both in km).
   (the detailed grid is used out to a square of half side-
   length 'dist1', border points in both inner and outer
   grid being attenuated to compensate edge effects.)
3: terrain corrections from 'dist1' to 'dist2' (km).
4: rtm gravity effect, computed as terrain correction
   to distance 'dist1' and bouguer reduction to reference
   level.
5: geoid effect for residual topography, condensation
   approximation. only one grid, second grid specification
   is the reference grid (lref = false means no ref used).
   effects will be computed in elevation 'dist2' (km) above
   the topography (i.e., above the condensation level).
   computed to distance 'dist1'.
6: isostatic geoid effects, second order expansion, in
   elevation 'dist2' above the geoid. (mixed continental and
   oceanic area using equivalent sea/rock root conversion).
   computed to distance 'dist1' (km).
7: gravity effect of the isostatic compensation,
   computed to distance 'dist1' (km).
8: rtm deflections of the vertical (condensation approximation)
   (written on outfile as two consecutive grids in unit arcsec)
   effect computed between distances 'dist1' to 'dist2' (km).""")
        LREF.draw(4)
        dist12 = geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal",
        "Distances of computation [km]:", "0 100", "See operation mode help.")
        filac1 = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal",
        "South-West corner of grid:", "55 10", "Degrees or meters for UTM.")
        ine1 = geomodule.NewEntry(frame, root, 7,
        geomodule.textwidth, "normal",
        "Number of points in North and East:", "120 120",
        "Numbers must be even.")

        geomodule.seperator_line(frame,20, "Operation mode 2 only")
        filac2 = geomodule.NewEntry(frame, root, 21,
        geomodule.textwidth, "normal",
        "South-West corner of subgrid:", "55.2 10.2",
        "Degrees or meters for UTM.")
        ine2 = geomodule.NewEntry(frame, root, 22,
        geomodule.textwidth, "normal",
        "Number of points in North and East:", "6 6", "Numbers must be even.")
        resultfile = geomodule.FileSelectorSave(frame, root, 31,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "result.dat", "nohelp")

        geomodule.seperator_line(frame,30,
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
