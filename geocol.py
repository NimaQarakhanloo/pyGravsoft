#!/usr/bin/python
# $Id: geocol.py 286 2009-07-26 15:36:32Z tjansson $
# vim: tabstop=4 expandtab shiftwidth=4
"""Geodetic Collocation GUI,
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
programname = "GEOCOL - Geodectic Collocation"
version = "0.7.5"
jobfile = "geocol.inp"
logfile = "geocol.log"
execute = "geocol17"
description = """ GEOCOL

Documentation can found in doc/geocol17_rev1.pdf
"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## Initiate all the button objects - they will first
# be drawn when foo.draw(linenumber) is called.
######################################################

        #inputformat = geomodule.NewEntry(frame,root, 16,
        #geomodule.textwidth, "normal",
        #"Input format", "(2I4,2D20.12)", "nohelp")
        #LFORM  = geomodule.NewRadioButton(frame, root,
        #"Are the coefficients formatted?", "select", "nocommand", "", "nohelp")
        LNCOL  = geomodule.NewRadioButton(frame, root,
        "Should collocation be used?" , "select", "nocommand", "", "nohelp")
        #LERR   = geomodule.NewRadioButton(frame, root,
        #"Should error estimates be computed ",
        #"unselect", "nocommand", "", "nohelp")
        LCOMP  = geomodule.NewRadioButton(frame, root,
        "Should computed values be subtracted from observed ",
        "unselect",
        lambda: [nd3.enable()],
        lambda: [nd3.disable()],
        "Data must be in column 5 of the input file.")
        LMAP   = geomodule.NewRadioButton(frame, root,
        "Should primitive map be output ", "select" ,
        "nocommand", "", "nohelp")
        LPUNCH = geomodule.NewRadioButton(frame, root,
        "Output to file ", "select", "nocommand", "", "nohelp")
        LMEAN  = geomodule.NewRadioButton(frame, root,
        "Should mean values be computed ", "unselect",
        "nocommand", "", "nohelp")
        LGRID  = geomodule.NewRadioButton(frame, root,
        "Should a grid be used in computations ", "select",
        lambda: [gridspecification.enable(), height.enable(), predictioninput.disable()],
        lambda: [gridspecification.disable(), height.disable(), predictioninput.enable()] , "nohelp")
        LPOT   = geomodule.NewRadioButton(frame, root,
        "Should spherical harmonic expansion be used?", "unselect",
        lambda: [LGRID.enable()],
        lambda: [LGRID.disable()], "")
        LSTAT  = geomodule.NewRadioButton(frame, root,
        "Should statistics be output ", "select",
        lambda: [histogrambinsize.enable()],
        lambda: [histogrambinsize.disable()], "nohelp" )

######################################################
## Functions -- these needs to be defined before they are used
######################################################

        def write_to_file():
            file = open(jobfile, 'w')
            if file:

                LPARAM = "f"
                if abs(int(datatype.textentry.get())) == 11:
                    LPARAM = "t"
                if datatype2.textentry.get() != "" and abs(int(datatype2.textentry.get())) == 11:
                    LPARAM = "t"

                wordlist = ["f\n", \
                "t f ", LPOT.state.get(), " f f ", LPARAM ," f f\n", \
                "1\n", \
                "neq\n", \
                "21 2000\n", \
                refsysListbox.get(ACTIVE), "\n" ]
                #"EGM", "\n"]
                #semimajor.textentry.get(), " 0.0 ", \
                #maximaldegree.textentry.get(), " f f ", \
                #" f f f \n"]
                #if LFORM.state.get() == "t":
                #    wordlist.extend([inputformat.textentry.get(), "\n"])
                #wordlist.extend([inputfile.textentry.get(), "\n", \
                wordlist.extend(["2\n", \
                "4\n", \
                covariancemodel.textentry.get(), " f f t f\n", \
                "-1 2 ", scale.textentry.get(), "\n", \
                errordegreefile.textentry.get(), "\n"])

                if LPARAM == "t":
                     wordlist.extend(["f\n"])

                ################### Section 3
                wordlist.extend([" -1 2 3 3 4 ", \
                str(int(nd.textentry.get())+4)," 0 ", \
                datatype.textentry.get(), " -1 0.0 ", \
                " f f f t f f f f f t \n",\
                datainput.textentry.get(), "\n", \
                "25 \n"])

                if LPARAM == "t":
                    if abs(int(datatype.textentry.get())) == 13: # also -13
                        wordlist.extend(["t 0\n"])
                    if abs(int(datatype.textentry.get())) == 11: # also -11
                        wordlist.extend(["t 1\n11111\n"])

                wordlist.extend([obserror.textentry.get(), "\n",])

                ### Second optional dataset
                if datainput2.textentry.get() == "":
                    if saneq.textentry.get() == "0":
                        wordlist.extend([ "t f \nf f f \n"])
                    else:
                        wordlist.extend([ "t t\nt ",
                        saneq.textentry.get(),"\nf f f \n"])
                else:
                    wordlist.extend([ "f f \n"])
                    wordlist.extend([" -1 2 3 3 4 ",
                    str(int(nd2.textentry.get())+4)," 0 ",
                    datatype2.textentry.get(), " -1 0.0 ",
                    " f f f t f f f f f t \n",\
                    datainput2.textentry.get(), "\n", \
                    "26 \n"])

                    if LPARAM == "t":
                        if abs(int(datatype2.textentry.get())) == 13:
                            wordlist.extend(["t 0\n"])
                        if abs(int(datatype2.textentry.get())) == 11:
                            wordlist.extend(["t 1\n11111\n"])

                    wordlist.extend([obserror2.textentry.get(), "\n",])
                    if saneq.textentry.get() == "0":
                        wordlist.extend([ "t f \nf f f \n"])
                    else:
                        wordlist.extend([ "t t\nt ",
                        saneq.textentry.get(),"\nf f f\n"])

                ################### Section 4
                wordlist.extend([LGRID.state.get(), " ",
                " t ", LCOMP.state.get(), " f\nf\n"])
                #if LERR.state.get() == "t":
                #    wordlist.extend(["f\n"])
                if LGRID.state.get() == "t":
                    wordlist.extend([gridspecification.textentry.get(), " ",
                    predictiontype.textentry.get(), " -1 ", height.textentry.get(),
                    " ", LMAP.state.get(), " ", LPUNCH.state.get(), " ",
                    LMEAN.state.get(), "\n"])
                    if LPUNCH.state.get() == "t":
                        wordlist.extend([resultfile.textentry.get(), "\n",
                        resultfile.textentry.get(), ".err", "\n"
                        "dcova\n"])

                    if LPARAM == "t":
                        if abs(int(predictiontype.textentry.get())) == 13: # also -13
                            wordlist.extend(["0\n"])
                        if abs(int(predictiontype.textentry.get())) == 11: # also -11
                            wordlist.extend(["1 1111\n"])

                else:
                    if LCOMP.state.get() == "t":
                        wordlist.extend([" 1 2 3 3 4 ",
                        str(int(nd3.textentry.get())+4), " 0 ",
                        predictiontype.textentry.get(), " -1 0.0 ",
                        LPUNCH.state.get(), " f f t f f ",
                        LSTAT.state.get(), " f f t \n",\
                        predictioninput.textentry.get(), "\n", \
                        "27 \n"])
                    else:
                        wordlist.extend([" 1 2 3 3 4 0 0 ",
                        predictiontype.textentry.get(),
                        " -1 0.0 ", LPUNCH.state.get(),
                        " f f t f f f f f t \n",\
                        predictioninput.textentry.get(), "\n", \
                        "27 \n"])

                    if LPUNCH.state.get() == "t":
                        wordlist.extend([resultfile.textentry.get(), "\n"])
                    if LSTAT.state.get() == "t":
                        wordlist.extend([histogrambinsize.textentry.get(),
                        "\n"])
                    if LPARAM == "t":
                        if abs(int(predictiontype.textentry.get())) == 13: # also -13
                            wordlist.extend(["t 0\n"])
                        if abs(int(predictiontype.textentry.get())) == 11: # also -11
                            wordlist.extend(["t 1\n11111\n"])

                    if LCOMP.state.get() == "t":
                        wordlist.extend(["0.0 \n"])                        
                        if LSTAT.state.get() == "t" and grrfile.textentry.get() != "":
                            wordlist.extend(["4 \n"])
                            wordlist.extend([grrfile.textentry.get(),"\n"])
                        else:
                            wordlist.extend(["-1 \n"])
                        
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

        #inputfile = geomodule.FileSelector(frame,root,14,
        #geomodule.textwidth, "normal",
        #"Input gravity model filepath:", "data/EGM96",
        #"Coefficients must be in free format.")
        #LFORM.draw(15)
        #semimajor = geomodule.NewEntry(frame,root,19,
        #geomodule.textwidth, "normal",
        #"Input GM, semi-major axis (M): ",
        #"3.986004415D14 6378136.3", "nohelp")
        #maximaldegree = geomodule.NewEntry(frame,root,20,
        #geomodule.textwidth, "normal",
        #"Input maximal degree: ", "360", "nohelp")

######################################################
## The second section
######################################################

        geomodule.seperator_line(frame, 24,
        'Analytic covariance function definition')
        covariancemodel = geomodule.NewEntry(frame,root,25,
        geomodule.textwidth, "normal",
        "Input convariance model parameters: ", "-1.0 400.0 360", """
(-1.0) THE DEPTH TO THE BJERHAMMAR SPHERE IN KM. Value (RE-RB)/1000 from COVFIT OUTPUT 

(400)  VARIANCE OF GRAVITY ANOMALIES AT ZERO ALTITUDE. Value VARG from COVFIT output.

(360)  MAXIMAL DEGREE FOR EMPIRICAL DEGREE-VARIANCES. Must correspond to maximal degree of EGM used in remove-step (GEOEGM).
 """)
        scale = geomodule.NewEntry(frame,root,26,
        geomodule.textwidth,"normal",
        "Input error degree variance scale factor :", "1.0", """
Vaule is AA in COVFIT output.        
        """)
        errordegreefile    = geomodule.FileSelector(frame,root,27,
        geomodule.textwidth,"normal",
        "Input name of error degree variance file: ",
        "data/egm96.edg", "nohelp")

######################################################
## The third section
######################################################

        geomodule.seperator_line(frame, 35,'Observation dataset parameters')
        datatype = geomodule.NewEntry(frame,root,36,
        geomodule.textwidth, "normal",
        "Input code for observations: ", "13", """

03 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
04 - DEFLECTION OF THE VERTICAL, PRIME VERTI.
11 - HEIGHT-ANOMALY OR GEOID UNDULATION
12 - GRAVITY DISTURBANCE
13 - GRAVITY ANOMALY
15 - VERTICAL GRAVITY DISTURBANCE GRADIENT""" )
        datainput = geomodule.FileSelector(frame,root, 37,
        geomodule.textwidth, "normal",
        "Input name of datafile (Gravsoft format): ",
        "data/nmfa-egm96-tc.dat", "nohelp")
        #datainput.disable()
        obserror = geomodule.NewEntry(frame,root, 38,
        geomodule.textwidth, "normal",
        "Observation error:", "0.1", "nohelp")
        nd = geomodule.NewEntry(frame,root, 39,
        geomodule.textwidth, "normal",
        "Data column number:", "1", "Column number after altitude.")

######################################################
## The third section
######################################################

        geomodule.seperator_line(frame, 45,
        'Second observation dataset parameters (optional)')
        datatype2 = geomodule.NewEntry(frame,root,46,
        geomodule.textwidth, "normal",
        "Input code for observations: ","","""

03 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
04 - DEFLECTION OF THE VERTICAL, PRIME VERTI.
11 - HEIGHT-ANOMALY OR GEOID UNDULATION
12 - GRAVITY DISTURBANCE
13 - GRAVITY ANOMALY
15 - VERTICAL GRAVITY DISTURBANCE GRADIENT""" )
        datainput2 = geomodule.FileSelector(frame,root, 47,
        geomodule.textwidth, "normal",
        "Input name of datafile (Gravsoft format): ", "", "nohelp")
        #datainput2.disable()
        obserror2 = geomodule.NewEntry(frame,root, 48,
        geomodule.textwidth, "normal",
        "Observation error:", "", "nohelp")
        nd2    = geomodule.NewEntry(frame,root, 49,
        geomodule.textwidth, "normal",
        "Data column number:", "", "Column number after altitude.")

######################################################
## The fourth section
######################################################

        geomodule.seperator_line(frame, 60,'Prediction type definition')
        saneq = geomodule.NewEntry(frame,root,61,
        geomodule.textwidth, "normal",
        "Number of already reduced equations: ","0", """
Already reduced normal equations maybe used again if number
of equations is different from 0.""" )
        predictiontype = geomodule.NewEntry(frame,root,63,
        geomodule.textwidth, "normal",
        "Input code for predictions: ", "11","""

03 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
04 - DEFLECTION OF THE VERTICAL, PRIME VERTI.
11 - HEIGHT-ANOMALY OR GEOID UNDULATION
12 - GRAVITY DISTURBANCE
13 - GRAVITY ANOMALY
15 - VERTICAL GRAVITY DISTURBANCE GRADIENT""" )
        #LERR.draw(65)
        LGRID.draw(66)
        gridspecification = geomodule.NewEntry(frame,root,67,
        geomodule.textwidth, "normal",
        "     Input grid specification :",
        "54.5 57.5 7.0 13.0 0.1  0.2",
        "Min, Max latitude, Min, Max longtitude, lat.spacing, long.spacing")
        height = geomodule.NewEntry(frame,root,68,
        geomodule.textwidth, "normal",
        "     Input grid altitude (m) :", "0.0", "nohelp")
        predictioninput = geomodule.FileSelector(frame,root,69,
        geomodule.textwidth, "normal",
        "Input name of predictionfile: ",
        "data/nmzeta.dat", """
The file must be in Gravsoft format if we
predict in points. """)
        predictioninput.disable()
        LCOMP.draw(72)
        nd3    = geomodule.NewEntry(frame,root,73,
        geomodule.textwidth, "normal",
        "     Data column number", "2", "nohelp")
        nd3.disable() 
        LSTAT.draw(74)
        histogrambinsize = geomodule.NewEntry(frame,root,75,
        geomodule.textwidth, "normal",
        "     Input histogram bin size ", "5.0", "nohelp")
        #LPUNCH.draw(78)
        grrfile = geomodule.FileSelectorSave(frame,root,78,
        geomodule.textwidth, "normal",
        "     File to hold suspected gross errors:" , "",
        "File with suspected gross errors is needed if both files are compared and statistics are outputed. Rejection level is fixed to 4. If no name is written no suspected gross errors be looked for.")
        resultfile = geomodule.FileSelectorSave(frame,root,79,
        geomodule.textwidth, "normal",
        "File to hold result:" , "grid11.ex4", """
If prediction in grid a grid 
filename is needed for predicted values. Error 
estimates will be written to a file with the 
same name but with a suffix .err""")

######################################################
##The Bottom line
######################################################

        geomodule.seperator_line(frame, geomodule.maxrow-2,
                                 "Running options. Working in "+geomodule.getcwd)
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
