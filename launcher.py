#!/usr/bin/python
# $Id: launcher.py 286 2009-07-26 15:36:32Z tjansson $
"""\
Launcher
Thomas R. N. Jansson
tjansson@fys.ku.dk"""

#Load graphical tool kit
try:
    from Tkinter import *
except ImportError:
    from tkinter import *    
    # This is the new name in python 3.

import os
import geomodule

geomodule.maxwidth  = 140
# Since the interface now has to columns of buttons this must be wider

## Constants
programname = "PyGravsoft Launcher - Gravity field programs"
version = "0.6"

## Starting the program and creating classes
class App:
    def __init__(self, master):
        frame = Frame(master)
        frame.grid()

######################################################
## The first section
######################################################

        #mainLabel = geomodule.mainline(frame, programname, version)


        geomodule.seperator_line(frame,10, "3D Applications")
        covfit = geomodule.LauncherButton(frame,11,8,0,"COVFIT",
        lambda:[covfit.run("covfit.py")], "Covariance Fitting")
        empcov = geomodule.LauncherButton(frame,11,8,1,"EMPCOV", lambda:
        [empcov.run("empcov.py")],	"Empirical Covariance Estimation")
        geocol	= geomodule.LauncherButton(frame,12,8,0,"GEOCOL", lambda:
        [geocol.run("geocol.py")], "Geodectic Collocation")
        geogrid = geomodule.LauncherButton(frame,12,8,1,"GEOGRID",
        lambda: [geogrid.run("geogrid.py")],
        "Gridding or Interpolation of Irregular Distributed Data")
        geoip = geomodule.LauncherButton(frame,13,8,0,"GEOIP", lambda:
        [geoip.run("geoip.py")], "Grid Interpolation")
        geoegm = geomodule.LauncherButton(frame,13,8,1,"GEOEGM",
        lambda:[geoegm.run("geoegm.py")], "Gravity Model Evaluation")
        stokes = geomodule.LauncherButton(frame,14,8,0,"STOKES", lambda:
        [stokes.run("stokes.py")],
        "Space Domain Integration for Geoid or Deflections of the Vertical")
        tc  = geomodule.LauncherButton(frame,14,8,1,"TC", lambda:
        [tc.run("tc.py")], "Gravimetric Terrain Effects")


        geomodule.seperator_line(frame,40, "2D Applications")
        spfour = geomodule.LauncherButton(frame,41,8,0,"SPFOUR", lambda:
        [spfour.run("spfour.py")],
        "Spherical Multiband FFT for Gravimetric Computations")
        tcfour = geomodule.LauncherButton(frame,41,8,1,"TCFOUR",
        lambda: [tcfour.run("tcfour.py")],
        "Terrain effect computation by FFT")
        covfft = geomodule.LauncherButton(frame,42,8,0,"COVFFT", lambda:
        [covfft.run("covfft.py")],
        "Estimation of 2D Covariance Functions Using FFT")
        geofour = geomodule.LauncherButton(frame,42,8,1,"GEOFOUR",
        lambda: [geofour.run("geofour.py")],
        "Planar FFT for Gravity Field Modelling")
        gpfit = geomodule.LauncherButton(frame,43,8,0,"GPFIT", lambda:
        [gpfit.run("gpfit.py")],
        "Fitting Flat-earth Covariance Function to Gravity Data")
        gpcol1 = geomodule.LauncherButton(frame,43,8,1,"GPCOL1", lambda:
        [gpcol1.run("gpcol1.py")],
        "Flat-earth Collocation")
        fitgeoid = geomodule.LauncherButton(frame,44,8,0,"FITGEOID", lambda:
        [gpcol1.run("fitgeoid.py")],
        "Fit surface to GPS levelling")
        tcgrid = geomodule.LauncherButton(frame,44,8,1,"TCGRID", lambda:
        [tcgrid.run("tcgrid.py")],
        "DTM Grids and Mean Terrain Surfaces for RTM Method")


        geomodule.seperator_line(frame,70, "Service Programs")
        selectgui = geomodule.LauncherButton(frame,71,8,0,"SELECT", lambda:
        [selectgui.run("selectgui.py")], "Select, Thin and/or Average Data")
        geomain	= geomodule.LauncherButton(frame,71,8,1, "GEOMAIN",
        lambda: [geomain.run("geomain.py")],
        "2 points: Distance and Azimuth or Reverse")
        trans = geomodule.LauncherButton(frame,72,8,0,"TRANS",
        lambda: [trans.run("trans.py")],
        "Transformation of Coordinates to or from a 2D or 3D System")
        n2zeta = geomodule.LauncherButton(frame,72,8,1,"N2ZETA",
        lambda: [n2zeta.run("n2zeta.py")],
        "Transformation of Geoid Heights to Height Anomalies")
        fcomp = geomodule.LauncherButton(frame,73,8,0,"FCOMP", lambda:
        [fcomp.run("fcomp.py")], "File Comparison")
        g2sur = geomodule.LauncherButton(frame,73,8,1,"G2SUR", lambda:
        [g2sur.run("g2sur.py")],
        "Conversion of GRAVSOFT Grids to SURFER Format")
        gbin = geomodule.LauncherButton(frame,74,8,0,"GBIN",
        lambda: [gbin.run("gbin.py")], "Convert Grid to Binary or Reverse")
        g2gmt = geomodule.LauncherButton(frame,74,8,1,"G2GMT", lambda:
        [g2gmt.run("g2gmt.py")],
        "Conversion of GRAVSOFT Grids to GMT Format")       
        gcomb = geomodule.LauncherButton(frame,75,8,0,"GCOMB", lambda:
        [gcomb.run("gcomb.py")], "Combining Two Grids")        
        glist = geomodule.LauncherButton(frame,75,8,1,"GLIST", lambda:
        [glist.run("glist.py")],
        "Converts grid file to list file")

        geomodule.seperator_line(frame,geomodule.maxrow-2)
        button = Button(frame, text="QUIT", width=6,
        command=frame.quit)
        button.grid(row=geomodule.maxrow, column=0, sticky=W)

######################################################
## Initiate the program and start program loop
######################################################

root = Tk()
app = App(root)
root.title(programname)
root.mainloop()
