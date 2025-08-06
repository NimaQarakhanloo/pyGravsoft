#!/usr/bin/python
# $Id: gpcol1.py 286 2009-07-26 15:36:32Z tjansson $
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
programname = "GPCOL1 - Flat-earth collocation"
version = "0.7.5"
jobfile = "gpcol1.inp"
logfile = "gpcol1.log"
execute = "gpcol1"
description = """ GPCOL1

Gravsoft manual can be found in doc/gravsoft-manual.pdf
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

        LFORM  = geomodule.NewRadioButton(frame, root,
                "Are the coefficients formatted?", "select", "nocommand", "", "nohelp")
        LBLOCK = geomodule.NewRadioButton(frame, root,
                "Use blocked computations?", "select",
                lambda: [fila.enable(), dlatlon.enable(),nmin.enable(), blatlon.enable()],
                lambda: [fila.disable(), dlatlon.disable(),nmin.disable(), blatlon.disable()], "nohelp")
        LNCOL  = geomodule.NewRadioButton(frame, root,
                "Should collocation be used?" , "select", "nocommand", "", "nohelp")
        LGRID  = geomodule.NewRadioButton(frame, root,
                "Should a grid be used in computations ", "select",
                lambda: [gridspecification.enable(), height.enable()],
                lambda: [gridspecification.disable(), height.disable()], "nohelp")
        LERR   = geomodule.NewRadioButton(frame, root,
                "Should error estimates be computed ", "unselect", 
                lambda: [errorspecification.enable()],
                lambda: [errorspecification.disable()], "Determines whether \
                        or not error estimates be computed")
        LCOMP  = geomodule.NewRadioButton(frame, root,
                "Should computed values be subtracted from observed ",
                "unselect",    "nocommand", "",
                "Data must be in column 5 of the input file.")
        LMAP   = geomodule.NewRadioButton(frame, root,
                "Should primitive map be output ", "unselect" ,
                "nocommand", "", "nohelp")
        LPUNCH = geomodule.NewRadioButton(frame, root,
                "Output to file ", "select", "nocommand", "", "nohelp")
        LMEAN  = geomodule.NewRadioButton(frame, root,
                "Should mean values be computed ", "unselect",
                "nocommand", "", "nohelp")
        LPOT   = geomodule.NewRadioButton(frame, root,
                "Should spherical harmonic expansion be used?", "select",
                lambda: [LFORM.enable(), LGRID.enable()],
                lambda: [LFORM.disable(),LGRID.disable()], "")
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
                wordlist = [""]
            
                if datainput2.textentry.get() == "":
                    wordlist.extend([ "1\n"])
                else: 
                    wordlist.extend([ "2\n"])
            
                wordlist.extend([ covariancemodel.textentry.get(), "\n"])
                wordlist.extend([ datatype.textentry.get(), " ", 
                obserror.textentry.get(), "\n",  
                datainput.textentry.get(), "\n"])
                
                if datainput2.textentry.get() != "":
                    wordlist.extend([datatype2.textentry.get(), 
                    " ", 
                    obserror2.textentry.get(),"\n", 
                    datainput2.textentry.get(), 
                    "\n" ])
                    
                mode = "1"
                if dfile.textentry.get() != "":
                    mode = "2"
                if LGRID.state.get() == "f":
                    mode = "3"
                
                wordlist.extend([ predictiontype.textentry.get(), " ", 
                    mode, " ", 
                    LBLOCK.state.get(), "\n"])
                
                if int(mode) >= 2:
                    if dfile.textentry.get() == "":
                        wordlist.extend([dfile1.textentry.get(), "\n"])
                    else:
                        wordlist.extend([dfile.textentry.get(), "\n"])
                else: 
                    wordlist.extend([gridspecification.textentry.get(),
                    " ", height.textentry.get(), "\n"])
                 
                if LBLOCK.state.get() == "t":
                    if int(mode) == 3:
                        wordlist.extend([fila.textentry.get(), "\n"])
                    wordlist.extend([dlatlon.textentry.get(), " ",
                    blatlon.textentry.get(), " ",
                    nmin.textentry.get(), "\n" ])
                    
                wordlist.extend([resultfile.textentry.get(), "\n"])
                
                if LBLOCK.state.get() == "f":
                    wordlist.extend([errorspecification.textentry.get(), "\n"])
                    
                file.writelines(wordlist)
                file.close
                statusbar.config("Data writtten to "+jobfile+
                " succesfully","normal")
            else:
                statusbar.config("ERROR writing to "+jobfile,"warning")

        def run_program():
            write_to_file()
            statusbar.config("Data send to "+execute,"normal")
            geomodule.runprogram(execute, jobfile, logfile)

######################################################
## The first section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)

        geomodule.seperator_line(frame, 24,
        'Logarithmic covariance function definition')
        covariancemodel = geomodule.NewEntry(frame,root,25,
        geomodule.textwidth, "normal",
        "Convariance model parameters: ", "7.0 8.0 30.0", """
(10.0)  Standard deviation of gravity (mgal).
(20.0)  Covariance parameter D [km].
(20.0)  Covariance paramater T [km]
Determined by running gpfit.""")

######################################################
## The second section
######################################################

        geomodule.seperator_line(frame, 35,'Observation dataset parameters')
        datatype = geomodule.NewEntry(frame,root,36,
        geomodule.textwidth, "normal",
        "Input code for observations: ", "03", """
01 - HEIGHT-ANOMALY OR GEOID UNDULATION
02 - PAIR OF DEFLECTIONS OF THE VERTICAL
03 - GRAVITY ANOMALY
05 - GEOID HEIGHT USED AS SLOPES WITH DLIM=10.0 KM 
06 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
07 - DEFLECTION OF THE VERTICAL, PRIME VERTI. """ )
        datainput = geomodule.FileSelector(frame,root, 37,
        geomodule.textwidth, "normal",
        "Input name of datafile (Gravsoft format): ",
        "coordinates.dat", "Data must be in first column after height.")
        #datainput.disable()
        obserror = geomodule.NewEntry(frame,root, 38,
        geomodule.textwidth, "normal",
        "Observation error:", "2.0", "nohelp")
    
######################################################
## The third section
######################################################

        geomodule.seperator_line(frame, 45,
        'Second observation dataset parameters (optional)')
        datatype2 = geomodule.NewEntry(frame,root,46,
        geomodule.textwidth, "normal",
        "Input code for observations: ","03","""
01 - HEIGHT-ANOMALY OR GEOID UNDULATION
02 - PAIR OF DEFLECTIONS OF THE VERTICAL
03 - GRAVITY ANOMALY
05 - GEOID HEIGHT USED AS SLOPES WITH DLIM=10.0 KM 
06 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
07 - DEFLECTION OF THE VERTICAL, PRIME VERTI. """ )
        datainput2 = geomodule.FileSelector(frame,root, 47,
        geomodule.textwidth, "normal",
        "Input name of datafile (Gravsoft format): ", "", 
        "Data must be in first column after height.")
        #datainput2.disable()
        obserror2 = geomodule.NewEntry(frame,root, 48,
        geomodule.textwidth, "normal",
        "Observation error:", "1.5", "nohelp")

######################################################
## The fourth section
######################################################

        geomodule.seperator_line(frame, 60,'Prediction type definition')
        predictiontype = geomodule.NewEntry(frame,root,63,
        geomodule.textwidth, "normal",
        "Input code for predictions: ", "03","""
01 - HEIGHT-ANOMALY OR GEOID UNDULATION
02 - PAIR OF DEFLECTIONS OF THE VERTICAL
03 - GRAVITY ANOMALY
05 - GEOID HEIGHT USED AS SLOPES WITH DLIM=10.0 KM 
06 - DEFLECTION OF THE VERTICAL, MERIDIAN COMP.
07 - DEFLECTION OF THE VERTICAL, PRIME VERTI. """ )
        LERR.draw(65)
        errorspecification = geomodule.FileSelector(frame,root,66,
        geomodule.textwidth, "normal",
        "     Input file for error estimates:",
        "efile","nohelp")
        LGRID.draw(67)
        gridspecification = geomodule.NewEntry(frame,root,68,
        geomodule.textwidth, "normal",
        "     Input grid specification :",
        "2 4 100 103 0.0166667 0.0166667",
        "Min, Max latitude, Min, Max longtitude, lat.spacing, long.spacing")
        height = geomodule.NewEntry(frame,root,69,
        geomodule.textwidth, "normal",
        "     Input grid altitude (m) :", "0", 
        "If no value is given a grid with variable heights will be used. ")
        dfile = geomodule.FileSelector(frame,root,70,
        geomodule.textwidth, "normal",
        "     Grid file with varying heights:", "", 
        "This field should be filled in if grid altitude is blank. Gravsoft grid file.")
        dfile1 = geomodule.FileSelector(frame,root,71,
        geomodule.textwidth, "normal",
        "Input name of prediction point file (Gravsoft format): ",
        "predictions.dat", "nohelp")
      
        
        geomodule.seperator_line(frame, 80,'Blocked computations')
        LBLOCK.draw(81)
        fila = geomodule.NewEntry(frame,root,82,
        geomodule.textwidth, "normal",
        "     Block boundaries:", "2 4 100 103", 
        "nohelp")
        dlatlon = geomodule.NewEntry(frame,root,83,
        geomodule.textwidth, "normal",
        "     Size of blocks in latitude and longtitude:", "1 1", 
        "nohelp")
        nmin = geomodule.NewEntry(frame,root,84,
        geomodule.textwidth, "normal",
        "     Minimum number of obs. in central block:", "3", """
Minimum number of obs in central block 
for a solution (the central block 
expanded with blat/3,blon/3)""")
        blatlon = geomodule.NewEntry(frame,root,85,
        geomodule.textwidth, "normal",
        "     Borders in latitude and longtitude:", "0.6 0.6", 
        "nohelp")

######################################################
##Invoke some buttons
######################################################
        LERR.invokeno()
        
######################################################
##The Bottom line
######################################################

        geomodule.seperator_line(frame, geomodule.maxrow-3,
                                 "Running options. Working in "
                                 +geomodule.getcwd)
        statusbar = geomodule.statusbar(frame, geomodule.maxrow-1)
        resultfile = geomodule.FileSelectorSave(frame,root,geomodule.maxrow-2,
        geomodule.textwidth, "normal",
        "Name of file to hold result:" , "ofile", "nohelp")

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
