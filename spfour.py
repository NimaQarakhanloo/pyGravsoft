#!/usr/bin/python
# $Id: spfour.py 286 2009-07-26 15:36:32Z tjansson $
"""SPFOUR GUI,
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
programname = "SPFOUR - Spherical multiband FFT for gravimetric computations"
version = "0.4"
jobfile = "spfour.inp"
logfile = "spfour.log"
execute = "spfour"
description = """SPFOUR

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
                gridfile.textentry.get(), "\n",
                resultfile.textentry.get(), "\n",
                mode.textentry.get()," ",LMEAN.state.get()," ",
                nref.textentry.get()," ", dist1.textentry.get()," ",
                dist2.textentry.get(), "\n"]
                if mode.textentry.get() == "1":
                    wordlist.extend([nmod1.textentry.get()," ",
                    nmod2.textentry.get(),"\n",])
                if (mode.textentry.get() == "3") or (mode.textentry.get() == "11"):
                    wordlist.extend([height.textentry.get(),"\n",])
                if mode.textentry.get() == "13":
                    wordlist.extend([ref.textentry.get(),"\n",])
                wordlist.extend([ficlac.textentry.get()," ",
                inne.textentry.get()," ",iwn.textentry.get(), "\n",])

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

        LMEAN = geomodule.NewRadioButton(frame, root,
        "Should mean values be computed " , "select",
        "nocommand", "", "nohelp")

        #mainLabel = geomodule.mainline(frame, programname, version)
        gridfile = geomodule.FileSelector(frame, root, 2,
        geomodule.textwidth, "normal", "Input grid file:",
        "grid.gri", "Grid on standard GRAVSOFT format.")
        mode = geomodule.NewEntry(frame, root, 3,
        geomodule.textwidth, "normal",
        "Operation mode:", "1", """
1   geoid prediction from gravity data
2   (deflections - not implemented yet)
3   harmonic upward continuation of gravity data
10  gravity effect of isostasy - 2nd order expansion
    (sign: positive on land. add to Bouguer anomalies)
11  bouguer effect in airborne gravimetry at 'height'
12  geoid terrain effect (mass layer approximation)
13  geoid RTM terrain effect (3rd order approximation)
    reffile must be exactly same grid parameters as the basic grid""")
        LMEAN.draw(4)
        nref= geomodule.NewEntry(frame, root, 5,
        geomodule.textwidth, "normal", "Number of reference parallels:",
        "2", """
This is the number of reference parallels used.
if nref = 1 just one fft is done, nref = 2
yields two ffts with reference on north and
south boundary of subgrid, nref = 3 yields
three FFT's with reference latitudes on
north, central and south parallel, nref = 4 ...""")
        dist1 = geomodule.NewEntry(frame, root, 6,
        geomodule.textwidth, "normal",
        "Inner range of kernel:", "0", "nohelp")
        dist2 = geomodule.NewEntry(frame, root, 7,
        geomodule.textwidth, "normal",
        "Outer range of kernel:", "999", "nohelp")
        nmod1 = geomodule.NewEntry(frame, root, 8,
        geomodule.textwidth, "normal",
        "Lower range of kernel modification:", "0", "nohelp")
        nmod2 = geomodule.NewEntry(frame, root, 9,
        geomodule.textwidth, "normal",
        "Upper range of kernel modification:", "0", "nohelp")
        height = geomodule.NewEntry(frame, root, 10,
        geomodule.textwidth, "normal",
        "Height [m]:", "0", "nohelp")

        geomodule.seperator_line(frame, 20, "Output grid definition")
        ficlac = geomodule.NewEntry(frame, root, 21,
        geomodule.textwidth, "normal", "Geographical coordinates (SW):",
        "55 10",
        "This is geographical coordinates (fic, lac) of the wanted SW-corner.")
        inne = geomodule.NewEntry(frame, root, 22,
        geomodule.textwidth, "normal",
        "Number of points in the tranformed grid:", "40 40", """
This is the number of points (inn, ine) actually transformed
(may be a subgrid). If the implied subgrid by fic,lac,inn,ine
is greater than the actual grid zero padding will be done.
The reference latitude will not be affected by the zero
padding, and the padded values will not be output (the output
grid might thus be smaller than the wanted number of points)""")
        iwn = geomodule.NewEntry(frame, root, 23,
        geomodule.textwidth, "normal", "Tapering windowing width:",
        "0", """
if the width (iwn) is greater than 0 the 'iwn' points closed
to the edge of the internal selected grid will be windowed by
a cosine-taper window. note: the windowing is done after taking
the data mean!""")
        geomodule.seperator_line(frame, 30, "Geoid RTM effects (mode = 13)")
        ref = geomodule.FileSelector(frame, root, 31, geomodule.textwidth,
        "normal", "RTM file:", "geoid.rtm",
        "This file must include the exactly the same grid as the grid file.")

        geomodule.seperator_line(frame, 90, "Running options. Working in "+geomodule.getcwd)
        resultfile = geomodule.FileSelectorSave(frame, root, 93,
        geomodule.textwidth, "normal",
        "File to hold ouput:" , "dtm_ref.gri", """
Here are output on grid format is written. For point
prediction both predicted values and error estimates
are written here. """)

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
