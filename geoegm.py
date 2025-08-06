#!/usr/bin/python
# $Id: geoegm.py 286 2009-07-26 15:36:32Z tjansson $
"""Gravity Model Evaluation GUI,
Thomas R. N. Jansson,
tjansson@fys.ku.dk
Change 2008-03-10 by cct. """

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
programname = "GEOEGM - Gravity Model Evaluation"
version = "1.1.6"
jobfile = "geoegm.inp"
logfile = "geoegm.log"
execute = "geocol17"
description = """GEOEGM - Gravity Model Evaluation

http://cct.gfy.ku.dk/geocol17_rev1.pdf
"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## Initiate all the button objects - they will first be drawn when foo.draw(linenumber) is called.
######################################################

        inputformat = geomodule.NewEntry(frame,root, 16,
        geomodule.textwidth, "normal", "    Input format",
        "(2I4,2D20.12)", "nohelp")
        LFORM = geomodule.NewRadioButton(frame,root,
        "Are the coefficients formatted?", "select",
        lambda: [inputformat.enable()],
        lambda: [inputformat.disable()],
         "nohelp")
        LNCOL = geomodule.NewRadioButton(frame,root,
        "Should collocation be used?" , "select","nocommand", "", "nohelp")
        LERR = geomodule.NewRadioButton(frame,root,
        "Should error estimates be computed ", "unselect",
        "nocommand", "", "nohelp")
        LCOMP = geomodule.NewRadioButton(frame,root,
        "Should computed values be subtracted from observed " ,
        "unselect", 
        lambda: [nd.enable()],
        lambda: [nd.disable()],
        "nohelp")
      #  LMAP = geomodule.NewRadioButton(frame,root,
      #  "Should primitive map be output ", "unselect" ,
      #  "nocommand", "", "nohelp")
        LPUNCH = geomodule.NewRadioButton(frame,root,
        "Output to file ", "select",
        lambda: [resultfile.enable()],
        lambda: [resultfile.disable()], "nohelp")
        #LMEAN  = geomodule.NewRadioButton(frame,root,
        #"Should mean values be computed " , "unselect",
        #"nocommand", "", "nohelp")
        LGRID = geomodule.NewRadioButton(frame,root,
        "Should a grid be used in computations ", "select",
        lambda:    [datainput.disable(),
        gridspecification.enable(),
        height.enable()],
        lambda: [datainput.enable(),
        gridspecification.disable(),
        height.disable()],
        "nohelp")
        LPOT = geomodule.NewRadioButton(frame,root,
        "Should spherical harmonic expansion be used?", "select",
        lambda: [LFORM.enable(), LGRID.enable()],
        lambda: [LFORM.disable(),LGRID.disable()], "")
        LSTAT  = geomodule.NewRadioButton(frame,root,
        "Should statistics be output ", "unselect",
        lambda: [histogrambinsize.enable()],
        lambda: [histogrambinsize.disable()], "nohelp" )

######################################################
## Functions -- these needs to be defined before they are used
######################################################

        def write_to_file():
            file = open(jobfile, 'w')
            if file:
                wordlist = ["f\n", \
                "f f ", LPOT.state.get(), " f f f ",
                LNCOL.state.get(), " f\n", \
                refsysListbox.get(ACTIVE), "\n" ,\
                "EGM", "\n", \
                semimajor.textentry.get(), " 0.0 ",
                maximaldegree.textentry.get(), " f f ",
                LFORM.state.get(), " f f \n"]
                if LFORM.state.get() == "t":
                    wordlist.extend([inputformat.textentry.get(), "\n"])
                wordlist.extend([inputfile.textentry.get(), "\n", \
                LGRID.state.get(), " ", LERR.state.get(), " ",
                LCOMP.state.get(), " f \n"])
                if LGRID.state.get() == "t":
                    wordlist.extend([gridspecification.textentry.get(), " ",
                    datatype.textentry.get(), " -1 ",
                    height.textentry.get(), " ", LPUNCH.state.get(), " ",
                    LPUNCH.state.get(), " f\n"])
                    if LPUNCH.state.get() == "t":
                        wordlist.extend([resultfile.textentry.get(), "\n"])
                else:
                    if LCOMP.state.get() == "t":
                        wordlist.extend([" -1 2 3 3 4 ", \
                        str(int(nd.textentry.get())+4), " 0 ", \
                        datatype.textentry.get(), " -1 0.0 ", \
                        LPUNCH.state.get(), " f f f f f ",
                        LSTAT.state.get(), " f f t \n",\
                        datainput.textentry.get(), "\n", \
                        "25 \n"])
                    else:
                        wordlist.extend([" -1 2 3 3 4 0 0 ",
                        datatype.textentry.get(), " -1 0.0 ",
                        LPUNCH.state.get(), " f f f f f f f f t \n",\
                        datainput.textentry.get(), "\n", \
                        "25 \n"])
                    if LPUNCH.state.get() == "t":
                        wordlist.extend([resultfile.textentry.get(), "\n"])
                    if LSTAT.state.get() == "t":
                        wordlist.extend([histogrambinsize.textentry.get(),
                        "\n"])
                wordlist.extend([ "t\n"])
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
        refsysListbox = Listbox(frame, height="2",
        selectmode="SINGLE", takefocus=0, width=geomodule.textwidth)
        choices2 = ["5 - GRS80", "7 - Best current",]
        for item in choices2:
            refsysListbox.insert(END, item)
        refsysListbox.selection_set(0)
        refsysLabel = Label(frame,text="Select reference system ")
        refsysLabel.grid(row=11, column=0, columnspan=2, sticky=W)
        refsysListbox.grid(row=11, column=1, columnspan=2, sticky=W)
        inputfile = geomodule.FileSelector(frame,root,14,
        geomodule.textwidth, "normal",
        "Input gravity model filepath:", "data/EGM96",
        "Data must be in column 5 of the input file.")
        LFORM.draw(15)
        semimajor = geomodule.NewEntry(frame,root,19,
        geomodule.textwidth, "normal",
        "Input GM, semi-major axis (M): ",
        "3.986004415D14 6378136.3", "nohelp")
        maximaldegree     = geomodule.NewEntry(frame,root,20,
        geomodule.textwidth, "normal",
        "Input maximal degree: ", "360", "nohelp")

######################################################
## Second section
######################################################

        geomodule.seperator_line(frame,64,"Configure parameters")
        datatype = geomodule.NewEntry(frame,root,65,
        geomodule.textwidth, "normal",
        "Input datatype code: ", "11", """
11 - HEIGHT-ANOMALY OR GEOID UNDULATION
03 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
04 - DEFLECTION OF THE VERTICAL, PRIME VERTI.
12 - GRAVITY DISTURBANCE
13 - GRAVITY ANOMALY
15 - VERTICAL GRAVITY DISTURBANCE GRADIENT""" )
        LGRID.draw(66)
        gridspecification = geomodule.NewEntry(frame,root,67,
        geomodule.textwidth, "normal", "     Input grid specification :",
        "54.5 57.5 7.0 13.0 0.1 0.2","""
The numbers show examples of:
Min and Max latitude, Min and Max longtitude, latitude spacing, longtitude spacing

Notice that numbers must be seperated by space.
""")
        height = geomodule.NewEntry(frame,root,68,
        geomodule.textwidth, "normal", "     Input grid altitude (m) :",
        "0.0", "nohelp")
        datainput = geomodule.FileSelector(frame,root,71,
        geomodule.textwidth, "normal",
        "Input name of datafile (Gravsoft format): ",
        "coordinates.dat", "nohelp")
        datainput.disable()
        LCOMP.draw(72)
        nd = geomodule.NewEntry(frame,root, 73,
        geomodule.textwidth, "normal",
        "     Data column number:", "1", "Column number after altitude.")
        nd.disable()
        LSTAT.draw(74)
        histogrambinsize = geomodule.NewEntry(frame,root,75,
        geomodule.textwidth, "normal", "     Input histogram bin size ",
        "5.0", "nohelp")
        histogrambinsize.disable()
        LPUNCH.draw(78)
        resultfile = geomodule.FileSelectorSave(frame,root,79,
        geomodule.textwidth, "normal",
        "     Name of file to hold result:" , "geoid.gri", "nohelp")

######################################################
##The Bottom line
######################################################

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
