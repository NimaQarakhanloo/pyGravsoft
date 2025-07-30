#!/usr/bin/python
# $Id: geoid.py 286 2009-07-26 15:36:32Z tjansson $
# vim: tabstop=4 expandtab shiftwidth=4
"""GEOID - ,
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
from SimpleDialog import SimpleDialog
import geomodule

## Constants
programname = "GEOID - Binary geoid grid interpolation"
version = "0.1"
jobfile = "geoid.inp"
logfile = "geoid.log"
execute = "geoid"
description = """GEOID

Linear interpolation and transformation program for binary
geoid gridfiles. This program is an alternative to 'geoip', and
differs a.o. in the way that subgrid are are not stored in memory,
but accessed everytime from disc.

The program may also be used for transformation between UTM, geographic
and cartesian coordinate systems, may convert ellipsoidal to orthometric
heights and vice versa. An option to fit geoid residuals to known control
may be activated (not in this version),
Input to the program is interactive.

"""

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## Functions -- these needs to be defined before they are used
######################################################

        LTREND  = geomodule.NewRadioButton(frame, root,
        "Fit the geoid to local heights (mode 2)" , "select", "nocommand", "", """Heights of the local points must be
given in a separate file of form
1 = statno, lat, lon, h, H (degrees)
2 = statno, lat, lon, h, H (deg,min,sec)
3 = statno, lat, lon, H, geoid height (degrees)
 (where h = ellipsoidal height, H = orthometric height)""")
 
        LINP  = geomodule.NewRadioButton(frame, root,
        "Data points from file?" , "select", "nocommand", "", "nohelp")

        def write_to_file():
            file = open(jobfile, 'w')
            if file:
                wordlist = [itask.textentry.get(),"\n"]
                
                if int(itask.textentry.get()) == 4:
                    wordlist.extend([ids.textentry.get(),"\n"])
                 
                wordlist.extend([gfile.textentry.get(), "\n"])

                if int(itask.textentry.get()) == 2:
                    wordlist.extend([LTREND.state.get(),"\n"])
                    if LTREND.state.get() == "t":
                        wordlist.extend([ffile.textentry.get(),"\n"])
                        wordlist.extend([iform.textentry.get(),"\n"])
                        wordlist.extend([itrend.textentry.get(),"\n"])
                
                wordlist.extend([LINP.state.get(),"\n"])
                wordlist.extend([dfile.textentry.get(),"\n"])
                wordlist.extend([ofile.textentry.get(),"\n"])
                wordlist.extend([itypi.textentry.get(),"\n"])
                
                if int(itypi.textentry.get()) == 3:
                    wordlist.extend([idati.textentry.get(),"\n"])
                if int(itypi.textentry.get()) == 4:
                    wordlist.extend([iell.textentry.get(),"\n"])
                if int(itask.textentry.get()) == 0:
                    wordlist.extend([itypo.textentry.get(),"\n"])
                if int(itypo.textentry.get()) == 3:
                    wordlist.extend([idato.textentry.get(),"\n"])
                if int(itypo.textentry.get()) == 4:
                    wordlist.extend([iello.textentry.get(),"\n"])
                                        
                file.writelines(wordlist)
                file.close
                
                statusbar.config("Data writtten to "+jobfile+" succesfully",
                                 "normal")
            else:
                statusbar.config("ERROR writing to "+jobfile, "warning")

        def run_program():
            write_to_file()
            geomodule.runprogram(execute, jobfile, logfile)
            statusbar.config("Data send to "+execute,"normal")

######################################################
## The graphical section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)
        
        itask = geomodule.NewEntry(frame,root,10,
        geomodule.textwidth, "normal",
        "Mode number: ",
        "1", """
0: just coordinate conversion (i.e., geoid height assumed to be zero)
1: interpolate geoid heights in same system as in file
2: convert ellipsoidal to "normal" heights
3: convert "normal" heights to ellipsoidal heights
4: interpolate geoid height and convert to other datum (e.g. ED50 geoid)
5: deflections of the vertical (simple linear differentiation of geoid)
6: deflections of the vertical in local system (e.g. ED50 deflections)""")

        ids = geomodule.NewEntry(frame,root,11,
        geomodule.textwidth, "normal",
        "Datum for geoid predictions (mode 4):",
        "1", """
        Input wanted datum for geoid prediction: 
2 = ED50, 
3 = NAD27, 
4 = Qornoq:""")        
        
        gfile = geomodule.FileSelector(frame,root,12,
        geomodule.textwidth, "normal",
        "Binary gridfile name:",
        "geoid2_fit.bin", "nohelp")
        
        LTREND.draw(13)   

        geomodule.seperator_line(frame, 20,"Trend options")
        
        ffile = geomodule.FileSelector(frame,root,21,
        geomodule.textwidth, "normal",
        "File with data for trend estimation",
        "geoid.dat", "nohelp")

        iform = geomodule.NewEntry(frame,root,22,
        geomodule.textwidth, "normal",
        "   File format",
        "1", """
Heights of the local points must be
given in a separate file of form
1 = statno, lat, lon, h, H (degrees)
2 = statno, lat, lon, h, H (deg,min,sec)
3 = statno, lat, lon, H, geoid height (degrees)
(where h = ellipsoidal height, H = orthometric height)""")
        
        itrend = geomodule.NewEntry(frame,root,23,
        geomodule.textwidth, "normal",
        "   Select type of trendfunction",
        "2", """
1 = bias fit (requires at least 1 known point)
2 = fit linear trend (requires min. 3 points or more)
3 = fit second order polynomial in x and y (min. 6 points)
4 = fit third order polynomial (min. 10 points)
5 = fit datumshift parameters dX,dY,dZ (min. 3 points)
6 = fit datumshift dX,dY,dZ,scale (min. 4 points)""")

#        LINP.draw(16)


        
        geomodule.seperator_line(frame, 30,"Input options")
        
        dfile = geomodule.FileSelector(frame,root,31,
        geomodule.textwidth, "normal",
        "Input file name:",
        "points_inp.dat", "nohelp")
        
        itypi = geomodule.NewEntry(frame,root,32,
        geomodule.textwidth, "normal",
        "Type of input:",
        "1", """If mode<=1 :
1 = statno, lat, lon (degrees)
2 = statno, lat, lon (deg,min,sec)
3 = statno, X, Y, Z (meter)
4 = statno, N, E (UTM, meter)

If mode>1:
1 = statno, lat, lon, height (degrees)
2 = statno, lat, lon, height (deg,min,sec)
3 = statno, X, Y, Z (meter)
4 = statno, N, E, height (meter)""")

        idati = geomodule.NewEntry(frame,root,33,
        geomodule.textwidth, "normal",
        "   Select datumshift for XYZ input (type 3):",
        "1", """1 = None (WGS84 to WGS84)
2 = WGS84 to ED50
3 = WGS84 to NAD27
4 = WGS84 to Qornoq
5 = NWL9D to WGS84""")

        iell = geomodule.NewEntry(frame,root,34,
        geomodule.textwidth, "normal",
        "   UTM ellipsoid and zone number (type 4):",
        "1 33", """1 = WGS84
2 = Hayford (ED50)
3 = Clarke (NAD27)
4 = Bessel""")

        
        
        geomodule.seperator_line(frame, 40,"Ouput options")

        itypo = geomodule.NewEntry(frame,root,41,
        geomodule.textwidth, "normal",
        "Type of output:",
        "1", """1 = statno, lat, lon (degrees)
2 = statno, lat, lon (deg,min,sec)
3 = statno, X, Y, Z (meter)
4 = statno, N, E (UTM, meter)""")

        idato = geomodule.NewEntry(frame,root,42,
        geomodule.textwidth, "normal",
        "   Datumshift for XYZ output (type 3):",
        "1", """1 = None (WGS84 to WGS84)
2 = ED50 to WGS84
3 = NAD27 to WGS84
4 = Qornoq to WGS84
5 = WGS84 to NWL9D""")

        iello = geomodule.NewEntry(frame,root,43,
        geomodule.textwidth, "normal",
        "   UTM ellipsoid and zone number (type 4):",
        "1 33", """1 = WGS84
2 = Hayford (ED50)
3 = Clarke (NAD27)
4 = Bessel""")

       
        
        geomodule.seperator_line(frame, 60,
                                 "Running options. Working in "
                                 +geomodule.getcwd)
                                 
        ofile = geomodule.FileSelectorSave(frame,root,61,
        geomodule.textwidth, "normal",
        "Output file name:",
        "out.dat", "nohelp")
                                         
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
