#!/usr/bin/env python3

#############################################################################
#
# AUTHORS:      Martin Mergili and Shiva P. Pudasaini
# CONTRIBUTORS: Massimiliano Alvioli, Matthias Benedikt, Emmanuel Delage,
#               Wolfgang Fellin, Jan-Thomas Fischer, Sigridur S. Gylfadottir,
#               Andreas Huber, Ivan Marchesini, Markus Metz, Markus Neteler,
#               Alexander Ostermann, Matthias Rauter
#
# PURPOSE:      The simulation model for avalanche and debris flows
#
# COPYRIGHT:    (c) 2013 - 2021 by the authors
#               (c) 2020 - 2021 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 2014 - 2021 by the University of Bonn
#               (c) 2000 - 2021 by the GRASS Development Team
#               (c) 1993 - 2021 by the R Development Core Team
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################

#%module
#% description: The mass flow simulation tool
#% keywords: Raster
#% keywords: Landslide
#% keywords: Numerical simulation
#%end

#%flag
#% key: a
#% description: Map plots of pressure and kinetic energy
#% guisection: flags
#%end

#%flag
#% key: e
#% description: Model execution
#% guisection: flags
#%end

#%flag
#% key: k
#% description: Keep result GRASS raster maps
#% guisection: flags
#%end

#%flag
#% key: m
#% description: Multiple model runs
#% guisection: flags
#%end

#%flag
#% key: t
#% description: Map plots of impact wave or tsunami height
#% guisection: flags
#%end

#%flag
#% key: v
#% description: Evaluation and visualization
#% guisection: flags
#%end

#%option
#% key: prefix
#% type: string
#% description: Prefix for output files and folders
#% required: no
#% multiple: no
#%end

#%option
#% key: cores
#% type: string
#% description: Number of cores to be used for multiple model runs
#% required: no
#% multiple: no
#%end

#%option
#% key: cellsize
#% type: string
#% description: Cell size in metres
#% required: no
#% multiple: no
#%end

#%option
#% key: phases
#% type: string
#% description: Phase(s) to be considered (maximum 3, comma-separated, x = mixture, s = solid, fs = fine solid, f = fluid, m = multi-phase)
#% required: no
#% multiple: yes
#%end

#%option
#% key: phasetext
#% type: string
#% description: Names of phases to be considered (three comma-separated values for multi-phase, one value otherwise)
#% required: no
#% multiple: yes
#%end

#%option
#% key: gravity
#% type: string
#% description: Gravity (m/s2)
#% required: no
#% multiple: no
#%end

#%option
#% key: limiter
#% type: string
#% description: Numerical limiter (1 = Minmod, 2 = Superbee, 3 = Woodward, 4 = van Leer)
#% required: no
#% multiple: no
#%end

#%option
#% key: controls
#% type: string
#% description: Controls (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: aoicoords
#% type: string
#% description: Set of coordinates delineating area of interest (N,S,E,W)
#% required: no
#% multiple: yes
#%end

#%option
#% key: elevation
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of elevation raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: rhrelease1
#% type: string
#% description: Ratio(s) of PHASE 1 release height (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: vhrelease
#% type: string
#% description: Variation of release height (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hrelease1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hrelease3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: trelease
#% type: string
#% gisprompt: old,raster,dcell
#% description: Release time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: trelstop
#% type: string
#% gisprompt: old,raster,dcell
#% description: Release hydrograph stop time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: tslide
#% type: string
#% gisprompt: old,raster,dcell
#% description: Time of initial sliding raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: tstop
#% type: string
#% gisprompt: old,raster,dcell
#% description: Stopping time raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: vinx3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release velocity in x direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: viny3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 3 release velocity in y direction raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of entrainment raster map (comma-separated)
#% required: no
#% multiple: no
#%end

#%option
#% key: rhentrmax1
#% type: string
#% description: Ratio(s) of maximum height of PHASE 1 entrainment (comma-separated)
#% required: no
#% multiple: no
#%end

#%option
#% key: vhentrmax
#% type: string
#% description: Variation of maximum height of entrainment (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hentrmax1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 1 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 2 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hentrmax3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of maximum height of PHASE 3 entrainment raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: density
#% type: string
#% description: Density parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: friction
#% type: string
#% description: Friction parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: viscosity
#% type: string
#% description: Viscosity parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: basal
#% type: string
#% description: Basal surface parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: transformation
#% type: string
#% description: Phase transformation parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: dynfric
#% type: string
#% description: Parameters for the dynamic adaptation of friction (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: slidepar
#% type: string
#% description: Parameters for initial sliding (search radius, exponent for weighting, coefficient for deformation)
#% required: no
#% multiple: yes
#%end

#%option
#% key: special
#% type: string
#% description: Special parameters (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: phi1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: phi2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: phi3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of internal friction angle of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: delta3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction angle of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: tufri
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of turbulent friction coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny1
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 1 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny2
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 2 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ny3
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of viscosity of PHASE 3 raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ambdrag
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of ambient drag coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: flufri
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of fluid friction coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: centr
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of entrainment coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: cvshear
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of shear velocity coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: deltab
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of basal friction difference raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans12
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1-PHASE 2 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans13
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 1-PHASE 3 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: ctrans23
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of PHASE 2-PHASE 3 transformation coefficient raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: zones
#% type: string
#% gisprompt: old,raster,cell
#% description: Zones raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: impactarea
#% type: string
#% gisprompt: old,raster,cell
#% description: Name of observed impact area raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hdeposit
#% type: string
#% gisprompt: old,raster,dcell
#% description: Name of observed deposit height raster map
#% required: no
#% multiple: no
#%end

#%option
#% key: hydrograph
#% type: string
#% description: Path(es) to hydrograph file(s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: hydrocoords
#% type: string
#% description: Pairs of coordinates of hydrograph(s) (x1, y1, x2, y2, ...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: frictiograph
#% type: string
#% description: Path(es) to frictiograph file(s)
#% required: no
#% multiple: no
#%end

#%option
#% key: transformograph
#% type: string
#% description: Path(es) to transformograph file(s)
#% required: no
#% multiple: no
#%end

#%option
#% key: sampling
#% type: string
#% description: Type of parameter sampling (positive number = random, 0 = controlled, negative number = OAT)
#% required: no
#% multiple: no
#%end

#%option
#% key: cfl
#% type: string
#% description: CFL criterion, initial time step length (s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: time
#% type: string
#% description: Time interval (s), stop for writing output (s)
#% required: no
#% multiple: yes
#%end

#%option
#% key: slomo
#% type: string
#% description: Factor for slow motion
#% required: no
#% multiple: no
#%end

#%option
#% key: thresholds
#% type: string
#% description: Lower thresholds for display of flow height, kinetic energy and pressure (comma-separated)
#% required: no
#% multiple: yes
#%end

#%option
#% key: profile
#% type: string
#% description: Pairs of coordinates for profile plots (x1,y1,x2,y2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: ctrlpoints
#% type: string
#% description: Pairs of coordinates for control points (x1,y1,x2,y2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: reftime
#% type: string
#% description: Reference time of reach for each control point (t1,t2,...)
#% required: no
#% multiple: yes
#%end

#%option
#% key: phexagg
#% type: string
#% description: Exaggeration of flow height in profile plots
#% required: no
#% multiple: no
#%end

#%option
#% key: orthophoto
#% type: string
#% description: Path to orthophoto tiff image
#% required: no
#% multiple: no
#%end


# Importing libraries

import grass.script as grass
from grass.script import core as grasscore
import math
from PIL import Image
import os
import queue
import random
import shutil
import subprocess
import sys
import threading
import time


# Defining fundamental functions, classes, and variables

queueLock = threading.Lock()  # queue for multi-core processing
workQueue = queue.Queue()

exitFlag = 0  # exit flag
ambvars = grass.gisenv()  # path to GRASS data
locpath = ambvars.GISDBASE + "/" + ambvars.LOCATION_NAME  # path to GRASS location
temppath = (
    ambvars.GISDBASE
    + "/"
    + ambvars.LOCATION_NAME
    + "/"
    + ambvars.MAPSET
    + "/.tmp/rtemp"
)  # path to temporary directory
scriptpath = "$HOME/.grass7/addons/etc/r.avaflow.rcode"  # path to R scripts
scriptpath2 = "$HOME/.grass7/addons/scripts"  # path to r.avaflow.mult
mainmapset = ambvars.MAPSET  # name of main mapset


def which(file):  # function for reading GRASS directory:
    for path in os.environ["PATH"].split(":"):
        if os.path.exists(path + "/" + file):
            return path + "/" + file
    return None


for i in range(70, 80):
    if not which("grass" + str(i)) is None:
        bingrass = which("grass" + str(i))  # reading grass binary


def ErrorMessage(specify):  # function for error message:
    grass.message(" ")
    grass.error("Please revise the " + specify + ".")
    grass.message(" ")
    sys.exit()


class myThread(threading.Thread):  # class for threading:
    def __init__(self, name, q):
        threading.Thread.__init__(self)
        self.name = name
        self.q = q

    def run(self):
        process_data(self.name, self.q)


def StartBatch(jid):  # function for start of multi-core processing:
    print("Executing model run %s" % jid)
    execute = "bash " + temppath + "/tmp" + str(jid) + "/batch" + str(jid)
    os.system(execute + " < /dev/null > " + temppath + "/out" + str(jid))
    print("Model run %s completed." % jid)
    return


def process_data(threadName, q):  # function supporting start of multi-core processing:
    while not exitFlag:
        queueLock.acquire()
        if not workQueue.empty():
            time.sleep(0.25)
            data = q.get()
            queueLock.release()
            StartBatch(data)
        else:
            queueLock.release()


mstring = [
    "_hflow",
    "_tflow",
    "_pflow",
    "_basechange",
    "_vflow",
    "_tsun",
    "_treach",
]  # list of names of raster maps


def corrasc(corrname):  # function for correcting output ascii rasters

    os.rename(corrname + ".asc", corrname + "o.asc")

    fpressureo = open(corrname + "o.asc", "r")
    fpressure = open(corrname + ".asc", "w")

    fdata = ""
    for i in range(0, 5):
        fdatah = fpressureo.readline()
        fdata = fdata + "\t".join(fdatah.split()) + "\n"

    fdatah = fpressureo.readline()
    fdatai = fdatah.split()
    if fdatai[0] == "NODATA_value":
        fdata = fdata + fdatah
        fdataj = fdatai[1]
    else:
        fdata = fdata + "NODATA_value\t-9999\n" + fdatah
        fdataj = "-9999"

    fdata = fdata + "".join(fpressureo.readlines())
    fdata = fdata.replace("corner", "center").replace(fdataj, "-9999")

    fpressure.write(fdata)

    fpressureo.close()
    fpressure.close()

    os.system("rm -f " + corrname + "o.asc")
    os.system("rm -f " + corrname + ".asc.aux.xml")


def main():  # starting main function

    # Setting flags and parameters

    aflag = flags["a"]
    eflag = flags["e"]
    kflag = flags["k"]
    mflag = flags["m"]
    tflag = flags["t"]
    vflag = flags["v"]
    pf = options["prefix"]
    cores = options["cores"]
    cellsize = options["cellsize"]
    phases = options["phases"]
    phasetext0 = options["phasetext"]
    gravity = options["gravity"]
    limiter = options["limiter"]
    controls = options["controls"]
    aoicoords = options["aoicoords"]
    elevation = options["elevation"]
    hrelease = options["hrelease"]
    rhrelease1 = options["rhrelease1"]
    vhrelease = options["vhrelease"]
    hrelease1 = options["hrelease1"]
    hrelease2 = options["hrelease2"]
    hrelease3 = options["hrelease3"]
    trelease = options["trelease"]
    trelstop = options["trelstop"]
    tslide = options["tslide"]
    stoptime = options["tstop"]
    vinx = options["vinx"]
    viny = options["viny"]
    vinx1 = options["vinx1"]
    viny1 = options["viny1"]
    vinx2 = options["vinx2"]
    viny2 = options["viny2"]
    vinx3 = options["vinx3"]
    viny3 = options["viny3"]
    hentrmax = options["hentrmax"]
    rhentrmax1 = options["rhentrmax1"]
    vhentrmax = options["vhentrmax"]
    hentrmax1 = options["hentrmax1"]
    hentrmax2 = options["hentrmax2"]
    hentrmax3 = options["hentrmax3"]
    phi1 = options["phi1"]
    phi2 = options["phi2"]
    phi3 = options["phi3"]
    delta1 = options["delta1"]
    delta2 = options["delta2"]
    delta3 = options["delta3"]
    tufri = options["tufri"]
    ny1 = options["ny1"]
    ny2 = options["ny2"]
    ny3 = options["ny3"]
    ambdrag = options["ambdrag"]
    flufri = options["flufri"]
    centr = options["centr"]
    cvshear = options["cvshear"]
    deltab = options["deltab"]
    ctrans12 = options["ctrans12"]
    ctrans13 = options["ctrans13"]
    ctrans23 = options["ctrans23"]
    zones = options["zones"]
    impactarea = options["impactarea"]
    hdeposit = options["hdeposit"]
    hydrograph = options["hydrograph"]
    hydrocoords = options["hydrocoords"]
    density = options["density"]
    friction = options["friction"]
    viscosity = options["viscosity"]
    basal = options["basal"]
    transformation = options["transformation"]
    special = options["special"]
    dynfric = options["dynfric"]
    frictiograph = options["frictiograph"]
    transformograph = options["transformograph"]
    sampling = options["sampling"]
    slidepar = options["slidepar"]
    cfl = options["cfl"]
    times = options["time"]
    slomo = options["slomo"]
    thresholds = options["thresholds"]
    profile = options["profile"]
    ctrlpoints = options["ctrlpoints"]
    reftime = options["reftime"]
    phexagg = options["phexagg"]
    orthophoto = options["orthophoto"]

    # Prefix
    if not pf:
        pf = "avf"

    print()
    print()
    print(
        "#############################################################################"
    )
    print("#")
    print("# MODULE:       r.avaflow 2.4 Pre-release")
    print("#")
    print("# AUTHORS:      Martin Mergili and Shiva P. Pudasaini")
    print("# CONTRIBUTORS: Massimiliano Alvioli, Matthias Benedikt, Emmanuel Delage,")
    print("#               Wolfgang Fellin, Jan-Thomas Fischer, Sigridur S. Gylfadottir,")
    print("#               Andreas Huber, Ivan Marchesini, Markus Metz, Markus Neteler,")
    print("#               Alexander Ostermann, Matthias Rauter")
    print("#")
    print("# PURPOSE:      The simulation model for avalanche and debris flows")
    print("# VERSION:      20210728 (28 July 2021)")
    print("#")
    print("# COPYRIGHT:    (c) 2013 - 2021 by the authors")
    print("#               (c) 2020 - 2021 by the University of Graz")
    print("#               (c) 2013 - 2021 by the BOKU University, Vienna")
    print("#               (c) 2015 - 2020 by the University of Vienna")
    print("#               (c) 2014 - 2021 by the University of Bonn")
    print("#               (c) 2000 - 2021 by the GRASS Development Team")
    print("#               (c) 1993 - 2021 by the R Development Core Team")
    print("#")
    print("#               This program is free software under the GNU General Public")
    print("#               License (>=v2). Read the file COPYING that comes with GRASS")
    print("#               for details.")
    print("#")
    print("# PREFIX:       %s" % pf)
    print("#")
    print(
        "#############################################################################"
    )
    print()
    print()

    # Making sure that all parameters are correctly defined

    # Flags
    if not eflag and not vflag:

        eflag = "1"
        vflag = "1"

    if mflag:
        aflag = "1"

    # Profile
    if not aflag:
        profile = None
    if profile:
        if not phexagg:
            phexagg = "1"
        try:
            phexagg = float(phexagg)
        except ValueError:
            ErrorMessage("profile exaggeration value")
        phexagg = str(phexagg)

    if eflag:  # for model execution mode:

        # Raster maps
        for inrast in [
            elevation,
            hrelease,
            hrelease1,
            hrelease2,
            hrelease3,
            hentrmax,
            hentrmax1,
            hentrmax2,
            hentrmax3,
            vinx,
            viny,
            vinx1,
            viny1,
            vinx2,
            viny2,
            vinx3,
            viny3,
            phi1,
            phi2,
            phi3,
            deltab,
            tufri,
            delta1,
            delta2,
            delta3,
            ny1,
            ny2,
            ny3,
            ambdrag,
            flufri,
            centr,
            cvshear,
            ctrans12,
            ctrans13,
            ctrans23,
            zones,
            impactarea,
            hdeposit,
        ]:  # !!!CHECK not yet working

            if inrast and (
                not grasscore.find_file(inrast, element="cell", mapset=".")
                and not grasscore.find_file(inrast, element="fcell", mapset=".")
                and not grasscore.find_file(inrast, element="dcell", mapset=".")
            ):
                ErrorMessage("GRASS raster maps")

        # Mesh spacing
        nmesh = "1"

        # Number and types of phases
        if not phases:
            phases = "m"
        if phases == "m":
            phases = "s,fs,f"
        phases = list(map(str, phases.split(",")))

        if phasetext0:
            phasetext = list(map(str, phasetext0.split(",")))
        if phasetext0 and not len(phasetext) == len(phases):
            ErrorMessage("phase descriptions")

        if len(phases) > 3:
            ErrorMessage("number of phases")
        else:
            if len(phases) == 1:
                model = 1
            elif len(phases) == 2:
                ErrorMessage(
                    "number of phases. Please use the multi-phase model (three phases) and provide no release data for the phase which is not needed."
                )
            elif len(phases) == 3:
                model = 7
            phasesnum = []
            if not phasetext0:
                phasetext = []

        for i in range(0, len(phases)):
            if (
                not (len(phases) == 1 and phases[i] == "x")
                and not phases[i] == "s"
                and not phases[i] == "fs"
                and not phases[i] == "f"
            ):
                ErrorMessage("types of phases")
            elif phases[i] == "x":
                phasesnum.append(1)
                model = 0
                if not phasetext0:
                    phasetext.append("Mixture")
            elif phases[i] == "s":
                phasesnum.append(1)
                if not phasetext0:
                    phasetext.append("Solid")
            elif phases[i] == "fs":
                phasesnum.append(2)
                if not phasetext0:
                    phasetext.append("Fine")
            elif phases[i] == "f":
                phasesnum.append(3)
                if not phasetext0:
                    phasetext.append("Fluid")

        if model <= 3:
            phasetext.append("NA")
            phasetext.append("NA")

        if not gravity:
            gravity = "9.81"

        if not limiter:
            limiter = "1"
        elif (
            not limiter == "1"
            and not limiter == "2"
            and not limiter == "3"
            and not limiter == "4"
        ):
            ErrorMessage("limiter")

        # Elevation
        if not elevation:
            ErrorMessage("elevation raster map")

        # Release
        if model <= 3 and not hydrograph and not hrelease:
            ErrorMessage("definition of flow release")
        elif (
            model == 7
            and not hydrograph
            and not hrelease1
            and not hrelease2
            and not hrelease3
            and not (hrelease and rhrelease1)
        ):
            ErrorMessage("definition of flow release")

        # Entrainable height
        if model == 7 and (hentrmax and not rhentrmax1):
            ErrorMessage("definition of entrainable height")

        # Hydrographs
        if hydrograph and not hydrocoords:
            ErrorMessage("input hydrograph coordinates")

        # Controls
        if not controls:
            controls = "0,1,0,0,0,0"
        controls = list(map(str, controls.split(",")))
        if not len(controls) == 6:
            ErrorMessage("number of controls")

        for econtrols in controls:
            try:
                econtrols = int(econtrols)
            except ValueError:
                ErrorMessage("control parameter values")

        # Thresholds
        if not thresholds:
            thresholds = "0.1,10000,10000,0.001"
        thresholds = list(map(str, thresholds.split(",")))
        if not len(thresholds) == 4:
            ErrorMessage("number of thresholds")
        thresholdsc = thresholds[3]

        for ethresholds in thresholds:
            try:
                ethresholds = float(ethresholds)
            except ValueError:
                ErrorMessage("threshold parameter values")

        # Parameter sampling
        if mflag:
            if not sampling:
                sampling = "100"
            try:
                sampling = int(sampling)
            except ValueError:
                ErrorMessage("value provided for sampling")
            sampling=str(sampling)

        # Parameters for initial sliding
        if not slidepar:
            slidepar = "0,0,0"
        slidepar = list(map(str, slidepar.split(",")))
        if not len(slidepar) == 3:
            ErrorMessage("number of parameters for initial sliding")

        for eslidepar in slidepar:
            try:
                eslidepar = float(eslidepar)
            except ValueError:
                ErrorMessage("initial sliding parameter values")

        # Flow parameters

        if not mflag:

            if model == 0:

                if not density:
                    density = "2000"
                if not friction:
                    friction = "35,20,3.0"
                if not basal:
                    basal = "-7.0,0.05,0.0,0.0"
                if not special:
                    special = "2,0"
                if not dynfric:
                    dynfric = "0.0,-6.0"
                lmax = 12

            elif model == 1:

                if not density:
                    density = "2700"
                if not friction:
                    friction = "35,20,0.05"
                if not viscosity:
                    viscosity = "-9999,-9999"
                if not basal:
                    basal = "-7.0,0.05,0.0,0.0"
                if not special:
                    special = "0.0,1,10,1,2,0"
                if not dynfric:
                    dynfric = "0.0,-6.0"
                lmax = 18

            elif model == 7:

                if not density:
                    density = "2700,1800,1000"
                if not friction:
                    friction = "35,20,0,0,0,0,0.05"
                if not viscosity:
                    viscosity = "-9999,-9999,-3.0,-9999,-3.0,0.0"
                if not basal:
                    basal = "-7.0,0.05,0.0,0.333,0.0"
                if not transformation:
                    transformation = "0.0,0.0,0.0"
                if not special:
                    special = "0.0,10,0.12,1,1,1,3,1,0.1,1,1,1,1,1,0,0,0,1,1,1,10,0,1,1,1,2,0"
                if not dynfric:
                    dynfric = "0.0,-6.0,0.0"
                lmax = 54

        elif int(sampling) > 0:

            if model == 0:

                if not density:
                    density = "2000,2000"
                if not friction:
                    friction = "35,35,20,20,3.0,3.0"
                if not basal:
                    basal = "-7.0,-7.0,0.05,0.05,0.0,0.0,0.0,0.0"
                if not special:
                    special = "2,2,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,-6.0,-6.0"
                lmax = 12

            elif model == 1:

                if not density:
                    density = "2700,2700"
                if not friction:
                    friction = "35,35,20,20,0.05,0.05"
                if not viscosity:
                    viscosity = "-9999,-9999,-9999,-9999"
                if not basal:
                    basal = "-7.0,-7.0,0.05,0.05,0.0,0.0,0.0,0.0"
                if not special:
                    special = "0.0,0.0,1,1,10,10,1,1,2,2,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,-6.0,-6.0"
                lmax = 18

            elif model == 7:

                if not density:
                    density = "2700,2700,1800,1800,1000,1000"
                if not friction:
                    friction = "35,35,20,20,0,0,0,0,0,0,0,0,0.05,0.05"
                if not viscosity:
                    viscosity = "-9999,-9999,-9999,-9999,-3.0,-3.0,-9999,-9999,-3.0,-3.0,0.0,0.0"
                if not basal:
                    basal = "-7.0,-7.0,0.05,0.05,0.0,0.0,0.333,0.333,0.0,0.0"
                if not transformation:
                    transformation = "0.0,0.0,0.0,0.0,0.0,0.0"
                if not special:
                    special = "0.0,0.0,10,10,0.12,0.12,1,1,1,1,1,1,3,3,1,1,0.1,0.1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,10,10,0,0,1,1,1,1,1,1,2,2,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,-6.0,-6.0,0.0,0.0"
                lmax = 54
        else:

            if model == 0:

                if not density:
                    density = "2000,2000,0"
                if not friction:
                    friction = "35,35,0,20,20,0,3.0,3.0,0"
                if not basal:
                    basal = "-7.0,-7.0,0,0.05,0.05,0,0.0,0.0,0,0.0,0.0,0"
                if not special:
                    special = "2,2,0,0,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,0,-6.0,-6.0,0"
                lmax = 12

            elif model == 1:

                if not density:
                    density = "2700,2700,0"
                if not friction:
                    friction = "35,35,0,20,20,0,0.05,0.05,0"
                if not viscosity:
                    viscosity = "-9999,-9999,0,-9999,-9999,0"
                if not basal:
                    basal = "-7.0,-7.0,0,0.05,0.05,0,0.0,0.0,0,0.0,0.0,0"
                if not special:
                    special = "0.0,0.0,0,1,1,0,10,10,0,1,1,0,2,2,0,0,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,0,-6.0,-6.0,0"
                lmax = 18

            elif model == 7:

                if not density:
                    density = "2700,2700,0,1800,1800,0,1000,1000,0"
                if not friction:
                    friction = "35,35,0,20,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0.05,0"
                if not viscosity:
                    viscosity = "-9999,-9999,0,-9999,-9999,0,-3.0,-3.0,0,-9999,-9999,0,-3.0,-3.0,0,0.0,0.0,0"
                if not basal:
                    basal = "-7.0,-7.0,0,0.05,0.05,0,0.0,0.0,0,0.333,0.333,0,0.0,0.0,0"
                if not transformation:
                    transformation = "0.0,0.0,0,0.0,0.0,0,0.0,0.0,0"
                if not special:
                    special = "0.0,0.0,0,10,10,0,0.12,0.12,0,1,1,0,1,1,0,1,1,0,3,3,0,1,1,0,0.1,0.1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,10,10,0,0,0,0,1,1,0,1,1,0,1,1,0,2,2,0,0,0,0"
                if not dynfric:
                    dynfric = "0.0,0.0,0,-6.0,-6.0,0,0.0,0.0,0"
                lmax = 54

        if model == 0:
            flowparam = (
                density + "," + friction + "," + basal + "," + special + "," + dynfric
            )
        elif model == 1:
            flowparam = (
                density
                + ","
                + friction
                + ","
                + viscosity
                + ","
                + basal
                + ","
                + special
                + ","
                + dynfric
            )
        else:
            flowparam = (
                density
                + ","
                + friction
                + ","
                + viscosity
                + ","
                + basal
                + ","
                + transformation
                + ","
                + special
                + ","
                + dynfric
            )
        flowparam = list(map(str, flowparam.split(",")))

        if not mflag and not len(flowparam) == lmax:
            ErrorMessage("number of flow parameters")
        elif mflag and int(sampling) > 0 and not len(flowparam) == 2 * lmax:
            ErrorMessage("number of flow parameters")
        elif mflag and int(sampling) <= 0 and not len(flowparam) == 3 * lmax:
            ErrorMessage("number of flow parameters")

        for eflowparam in flowparam:
            try:
                eflowparam = float(eflowparam)
            except ValueError:
                ErrorMessage("parameter values")

        # CFL criterion and initial time step length
        if not cfl:
            cfl = "0.40,0.001"
        cfl = list(map(str, cfl.split(",")))
        if not len(cfl) == 2:
            ErrorMessage("number of cfl parameters")

        for ecfl in cfl:
            try:
                ecfl = float(ecfl)
            except ValueError:
                ErrorMessage("cfl values")

        # Time interval and stop time
        if not times:
            times = "10,300"
        times = list(map(str, times.split(",")))
        tint = times[0]
        tstop = times[1]

        if mflag:
            tint = tstop

        for etimes in times:
            try:
                etimes = float(etimes)
            except ValueError:
                ErrorMessage("times")

        # Control point information
        if ctrlpoints:  # if control point coordinates are defined:

            ctrlpoints = list(
                map(str, ctrlpoints.split(","))
            )  # splitting list of control point coordinates
            lctrlpts = int(len(ctrlpoints) / 2)  # number of control points

            xctrl = []
            yctrl = []
            for ictrlpts in range(0, lctrlpts):  # extracting control point coordinates:
                xctrl.append(ctrlpoints[2 * ictrlpts])
                yctrl.append(ctrlpoints[2 * ictrlpts + 1])
        else:
            lctrlpts = 0

        if reftime:  # if reference times of reach are defined:

            reftime = list(
                map(str, reftime.split(","))
            )  # splitting list of reference times of reach

            if not len(reftime) == lctrlpts:
                ErrorMessage("number of reference times of reach")

        # Factor for slow motion
        if not slomo or not (phases == "fs" or phases == "f"):
            slomo = "1"

        # Number of cores
        if mflag:
            if not cores:
                cores = "8"
            try:
                cores = int(cores)
            except ValueError:
                ErrorMessage("number of cores")
            cores=str(cores)

    # Preparing environment

    print("1. PREPARING ENVIRONMENT.")

    os.system("rm -rf " + temppath)  # if necessary, deleting temporary directory
    os.system("mkdir " + temppath)  # creating temporary directory

    os.environ["GRASS_VERBOSE"] = "-1"
    grass.run_command(
        "g.gisenv", set="GRASS_VERBOSE=-1"
    )  # suppressing errors and warnings (overruled for some functions)

    ascpath = pf + "_results/" + pf + "_ascii/"  # path to directory with ascii rasters
    filepath = (
        pf + "_results/" + pf + "_files/"
    )  # path to directory with result text files

    if eflag:  # for model execution mode:

        if os.path.exists(pf + "_results"):
            os.system(
                "rm -rf " + pf + "_results"
            )  # if result directory already exists, deleting it
        os.mkdir(pf + "_results")  # creating directory for results
        os.mkdir(
            pf + "_results/" + pf + "_files"
        )  # creating directory for result files
        os.mkdir(
            pf + "_results/" + pf + "_plots"
        )  # creating directory for result plots
        if not mflag:
            os.mkdir(
                pf + "_results/" + pf + "_plots/" + pf + "_maps_timesteps"
            )  # creating directory for result maps of all time steps
        if not mflag:
            os.mkdir(pf + "_results/" + pf + "_plots/" + pf + "_profiles_timesteps")
        # creating directory for result profiles of all time steps
        os.mkdir(ascpath)  # creating directory for result ascii rasters

        # Defining GRASS region

        grass.run_command("g.region", flags="d")  # setting default region
        if cellsize:
            grass.run_command("g.region", flags="a", res=cellsize)  # updating cell size
        if aoicoords:  # updating extent:
            aoicoords = list(map(str, aoicoords.split(",")))
            if not len(aoicoords) == 4:
                ErrorMessage("number of coordinates of the area of interest")
            grass.run_command(
                "g.region",
                flags="a",
                n=aoicoords[0],
                s=aoicoords[1],
                w=aoicoords[2],
                e=aoicoords[3],
            )

        grass.mapcalc('"_aoi"=pow("%s",0)' % elevation, overwrite=True, quiet=True)
        cellsize = grass.raster_info("_aoi")["nsres"]  # reading cell size
        rnorth = grass.raster_info("_aoi")["north"]  # reading boundaries of aoi map
        rsouth = grass.raster_info("_aoi")["south"]
        rwest = grass.raster_info("_aoi")["west"]
        reast = grass.raster_info("_aoi")["east"]

        grass.run_command("g.region", flags="d")  # setting default region
        grass.run_command(
            "g.region", flags="a", n=rnorth, s=rsouth, w=rwest, e=reast
        )  # updating bounds
        grass.run_command("g.region", flags="a", res=cellsize)  # updating cell size

        # Pre-processing input raster maps

        grass.run_command(
            "g.copy", rast=(elevation, pf + "_elev"), overwrite=True, quiet=True
        )  # elevation map
        grass.run_command(
            "r.out.gdal",
            input=pf + "_elev",
            output=ascpath + pf + "_elev.asc",
            format="AAIGrid",
            overwrite=True,
        )  # exporting to ascii
        corrasc(ascpath + pf + "_elev")

        if model <= 3:

            if hrelease:
                grass.mapcalc(
                    '"%s_hflow0000"=if(isnull("%s")==1,0,"%s")'
                    % (pf, hrelease, hrelease),
                    overwrite=True,
                    quiet=True,
                )
            else:
                grass.mapcalc(
                    '"%s_hflow0000"=0' % pf, overwrite=True, quiet=True
                )  # release height map

            if not controls[3] == "0":

                if hentrmax:
                    grass.run_command(
                        "g.copy",
                        rast=(hentrmax, pf + "_hentrmax"),
                        overwrite=True,
                        quiet=True,
                    )  # maximum height of entrainment map
                else:
                    grass.mapcalc(
                        '"%s_hentrmax"=10000' % pf, overwrite=True, quiet=True
                    )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_hentrmax",
                    output=ascpath + pf + "_hentrmax.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_hentrmax")

            if vinx:
                grass.run_command(
                    "g.copy",
                    rast=(vinx, pf + "_vflowx0000"),
                    overwrite=True,
                    quiet=True,
                )  # x velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowx0000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowx0000",
                    output=ascpath + pf + "_vflowx0000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowx0000")

            if viny:
                grass.run_command(
                    "g.copy",
                    rast=(viny, pf + "_vflowy0000"),
                    overwrite=True,
                    quiet=True,
                )  # y velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowy0000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowy0000",
                    output=ascpath + pf + "_vflowy0000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowy0000")

        elif model == 7:

            if hrelease1:
                grass.run_command(
                    "g.copy",
                    rast=(hrelease1, pf + "_hflow10000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 1 release height map
            else:
                grass.mapcalc('"%s_hflow10000"=0' % pf, overwrite=True, quiet=True)
            grass.run_command(
                "r.out.gdal",
                input=pf + "_hflow10000",
                output=ascpath + pf + "_hflow10000.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting to ascii
            corrasc(ascpath + pf + "_hflow10000")

            if hrelease2:
                grass.run_command(
                    "g.copy",
                    rast=(hrelease2, pf + "_hflow20000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 2 release height map
            else:
                grass.mapcalc('"%s_hflow20000"=0' % pf, overwrite=True, quiet=True)
            grass.run_command(
                "r.out.gdal",
                input=pf + "_hflow20000",
                output=ascpath + pf + "_hflow20000.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting to ascii
            corrasc(ascpath + pf + "_hflow20000")

            if vinx1:
                grass.run_command(
                    "g.copy",
                    rast=(vinx1, pf + "_vflowx10000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 1 x velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowx10000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowx10000",
                    output=ascpath + pf + "_vflowx10000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowx10000")

            if viny1:
                grass.run_command(
                    "g.copy",
                    rast=(viny1, pf + "_vflowy10000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 1 y velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowy10000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowy10000",
                    output=ascpath + pf + "_vflowy10000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowy10000")

            if vinx2:
                grass.run_command(
                    "g.copy",
                    rast=(vinx2, pf + "_vflowx20000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 2 x velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowx20000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowx20000",
                    output=ascpath + pf + "_vflowx20000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowx20000")

            if viny2:
                grass.run_command(
                    "g.copy",
                    rast=(viny2, pf + "_vflowy20000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 2 y velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowy20000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowy20000",
                    output=ascpath + pf + "_vflowy20000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowy20000")

            grass.mapcalc(
                '"%s_hflow0000"="%s_hflow10000"+"%s_hflow20000"' % (pf, pf, pf),
                overwrite=True,
                quiet=True,
            )  # total release height map

            if not controls[3] == "0":

                if hentrmax1 or hentrmax2 or hentrmax3:

                    if hentrmax1:
                        grass.run_command(
                            "g.copy",
                            rast=(hentrmax1, pf + "_hentrmax1"),
                            overwrite=True,
                            quiet=True,
                        )  # maximum height of PHASE 1 entrainment map
                    else:
                        grass.mapcalc(
                            '"%s_hentrmax1"=0' % pf, overwrite=True, quiet=True
                        )
                else:
                    grass.mapcalc(
                        '"%s_hentrmax1"=10000' % pf, overwrite=True, quiet=True
                    )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_hentrmax1",
                    output=ascpath + pf + "_hentrmax1.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_hentrmax1")

                if hentrmax2 or hentrmax1 or hentrmax3:

                    if hentrmax2:
                        grass.run_command(
                            "g.copy",
                            rast=(hentrmax2, pf + "_hentrmax2"),
                            overwrite=True,
                            quiet=True,
                        )  # maximum height of PHASE 2 entrainment map
                    else:
                        grass.mapcalc(
                            '"%s_hentrmax2"=0' % pf, overwrite=True, quiet=True
                        )
                else:
                    grass.mapcalc(
                        '"%s_hentrmax2"=10000' % pf, overwrite=True, quiet=True
                    )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_hentrmax2",
                    output=ascpath + pf + "_hentrmax2.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_hentrmax2")

            if hrelease3:
                grass.run_command(
                    "g.copy",
                    rast=(hrelease3, pf + "_hflow30000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 3 release height map
            else:
                grass.mapcalc('"%s_hflow30000"=0' % pf, overwrite=True, quiet=True)
            grass.run_command(
                "r.out.gdal",
                input=pf + "_hflow30000",
                output=ascpath + pf + "_hflow30000.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting to ascii
            corrasc(ascpath + pf + "_hflow30000")

            grass.mapcalc(
                '"%s_hflow0000"="%s_hflow0000"+"%s_hflow30000"' % (pf, pf, pf),
                overwrite=True,
                quiet=True,
            )  # total release height map

            if vinx3:
                grass.run_command(
                    "g.copy",
                    rast=(vinx3, pf + "_vflowx30000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 3 x velocity map
            elif aflag:
                grass.mapcalc('"%s_vflowx30000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowx30000",
                    output=ascpath + pf + "_vflowx30000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowx30000")

            if viny3:
                grass.run_command(
                    "g.copy",
                    rast=(viny3, pf + "_vflowy30000"),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 3 y velocity map
            else:
                grass.mapcalc('"%s_vflowy30000"=0' % pf, overwrite=True, quiet=True)
            if aflag:
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_vflowy30000",
                    output=ascpath + pf + "_vflowy30000.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_vflowy30000")

            if not controls[3] == "0":

                if hentrmax3 or hentrmax1 or hentrmax2:

                    if hentrmax3:
                        grass.run_command(
                            "g.copy",
                            rast=(hentrmax3, pf + "_hentrmax3"),
                            overwrite=True,
                            quiet=True,
                        )  # maximum height of PHASE 2 entrainment map
                    else:
                        grass.mapcalc(
                            '"%s_hentrmax3"=0' % pf, overwrite=True, quiet=True
                        )
                else:
                    grass.mapcalc(
                        '"%s_hentrmax3"=10000' % pf, overwrite=True, quiet=True
                    )
                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_hentrmax3",
                    output=ascpath + pf + "_hentrmax3.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(ascpath + pf + "_hentrmax3")

        if impactarea:
            grass.run_command(
                "g.copy",
                rast=(impactarea, pf + "_impactarea"),
                overwrite=True,
                quiet=True,
            )  # observed impact area map

            maximpactarea = grass.raster_info(impactarea)[
                "max"
            ]  # maximum value of impact area
            if maximpactarea > 0:
                impdef = "1"  # setting control for impact area map
            else:
                impdef = "0"

            grass.run_command(
                "r.out.gdal",
                input=pf + "_impactarea",
                output=ascpath + pf + "_impactarea.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting to ascii
            corrasc(ascpath + pf + "_impactarea")
        else:
            impdef = "0"  # setting control for impact area map to negative

        if hdeposit:
            grass.run_command(
                "g.copy", rast=(hdeposit, pf + "_hdeposit"), overwrite=True, quiet=True
            )  # observed height of deposit map

            maxhdeposit = grass.raster_info(hdeposit)[
                "max"
            ]  # maximum value of height of deposit
            if maxhdeposit > 0:
                depdef = "1"  # setting control for deposit map
            else:
                depdef = "0"

            grass.run_command(
                "r.out.gdal",
                input=pf + "_hdeposit",
                output=ascpath + pf + "_hdeposit.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting to ascii
            corrasc(ascpath + pf + "_hdeposit")
        else:
            depdef = "0"  # setting control for height of deposit map to negative

        if (
            model == 7 and hrelease and rhrelease1
        ):  # if release heights are defined by total height and ratio of PHASE 1:

            hrelease1 = pf + "_hflow10000"
            hrelease2 = pf + "_hflow20000"
            hrelease3 = pf + "_hflow30000"

            if mflag == 0:  # for single model run:

                grass.mapcalc(
                    '"%s"="%s"*%s'
                    % (hrelease1, hrelease.replace(".asc", ""), rhrelease1),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 3 release height map
                grass.mapcalc(
                    '"%s"=0' % hrelease2, overwrite=True, quiet=True
                )  # PHASE 2 release height map
                grass.mapcalc(
                    '"%s"="%s"*(1-%s)'
                    % (hrelease3, hrelease.replace(".asc", ""), rhrelease1),
                    overwrite=True,
                    quiet=True,
                )  # PHASE 3 release height map

                grass.run_command(
                    "r.out.gdal",
                    input=hrelease1,
                    output=ascpath + hrelease1 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting PHASE 1 release height map to ascii
                grass.run_command(
                    "r.out.gdal",
                    input=hrelease2,
                    output=ascpath + hrelease2 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting PHASE 2 release height map to ascii
                grass.run_command(
                    "r.out.gdal",
                    input=hrelease3,
                    output=ascpath + hrelease3 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting PHASE 3 release height map to ascii

                corrasc(ascpath + hrelease1)
                corrasc(ascpath + hrelease2)
                corrasc(ascpath + hrelease3)

            grass.mapcalc(
                '"%s_hflow0000"="%s"' % (pf, hrelease.replace(".asc", "")),
                overwrite=True,
                quiet=True,
            )  # release height map

        grass.mapcalc(
            '"%s_hflow_max0000"="%s_hflow0000"' % (pf, pf), overwrite=True, quiet=True
        )  # maximum release height map
        grass.run_command(
            "r.out.gdal",
            input=pf + "_hflow0000",
            output=ascpath + pf + "_hflow0000.asc",
            format="AAIGrid",
            overwrite=True,
        )  # exporting to ascii
        corrasc(ascpath + pf + "_hflow0000")
        grass.run_command(
            "r.out.gdal",
            input=pf + "_hflow_max0000",
            output=ascpath + pf + "_hflow_max0000.asc",
            format="AAIGrid",
            overwrite=True,
        )  # exporting to ascii
        corrasc(ascpath + pf + "_hflow_max0000")

        if (
            model and hentrmax and rhentrmax1 and not controls[3] == "0"
        ):  # if maximum heights of entrainment are defined by total height and ratio of PHASE 1:

            hentrmax1 = pf + "_hentrmax1"
            hentrmax2 = pf + "_hentrmax2"
            hentrmax3 = pf + "_hentrmax3"

            if mflag == 0:  # for single model run:
                grass.mapcalc(
                    '"%s"="%s"*%s'
                    % (hentrmax1, hentrmax.replace(".asc", ""), rhentrmax1),
                    overwrite=True,
                    quiet=True,
                )  # maximum PHASE 1 height of entrainment map
                grass.mapcalc(
                    '"%s"=0' % hentrmax2, overwrite=True, quiet=True
                )  # maximum PHASE 2 height of entrainment map
                grass.mapcalc(
                    '"%s"="%s"*(1-%s)'
                    % (hentrmax3, hentrmax.replace(".asc", ""), rhentrmax1),
                    overwrite=True,
                    quiet=True,
                )  # maximum PHASE 3 height of entrainment map

                grass.run_command(
                    "r.out.gdal",
                    input=hentrmax1,
                    output=ascpath + hentrmax1 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting maximum height of PHASE 1 entrainment map to ascii
                grass.run_command(
                    "r.out.gdal",
                    input=hentrmax2,
                    output=ascpath + hentrmax2 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting maximum height of PHASE 2 entrainment map to ascii
                grass.run_command(
                    "r.out.gdal",
                    input=hentrmax3,
                    output=ascpath + hentrmax3 + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting maximum height of PHASE 3 entrainment map to ascii

                corrasc(ascpath + hentrmax1)
                corrasc(ascpath + hentrmax2)
                corrasc(ascpath + hentrmax3)

            grass.mapcalc(
                '"%s_hentrmax"="%s"' % (pf, hentrmax.replace(".asc", "")),
                overwrite=True,
                quiet=True,
            )  # maximum height of entrainment map

        time.sleep(1.0)  # waiting to make sure that all rasters have been built

        if (
            model <= 3
        ):  # controls for release masses, release velocities, and entrainment:

            if hrelease:
                releasemass = 1
            else:
                releasemass = 0

            if vinx or viny:
                releasevelocity = 1
            else:
                releasevelocity = 0

            if not controls[3] == "0" and hentrmax:
                entrainment = 1
            else:
                entrainment = 0

        elif model == 7:

            if hrelease1 or hrelease2 or hrelease3:
                releasemass = 1
            else:
                releasemass = 0

            if vinx1 or viny1:
                releasevelocity = 1
            else:
                releasevelocity = 0

            if vinx2 or viny3:
                releasevelocity2 = 1
            else:
                releasevelocity2 = 0

            if vinx3 or viny3:
                releasevelocity3 = 1
            else:
                releasevelocity3 = 0

            if not controls[3] == "0" and (hentrmax1 or hentrmax2 or hentrmax3):
                entrainment = 1
            else:
                entrainment = 0

        # Pre-processing hydrograph information

        if hydrograph:  # if hydrograph file is defined:

            hydrograph = list(
                map(str, hydrograph.split(","))
            )  # splitting list of hydrographs
            hyd_nin = len(hydrograph)

            hydlines = []
            for hydi in range(0, hyd_nin):
                fhyd = open(hydrograph[hydi], "r")  # opening hydrograph file
                hydlines.append(sum(1 for line in fhyd) - 2)  # counting number of lines
                fhyd.close()  # closing file
        else:
            hyd_nin = 0

        if hydrocoords:  # if hydrograph coordinates are defined:

            hydrocoords = list(
                map(str, hydrocoords.split(","))
            )  # splitting hydrograph coordinate input
            hyd_nout = int(len(hydrocoords) / 4 - hyd_nin)

        # Pre-processing frictiograph information

        if frictiograph:  # if frictiograph file is defined:

            ffri = open(frictiograph, "r")  # opening frictiograph file
            frilines = sum(1 for line in ffri) - 2  # counting number of lines
            ffri.close()  # closing file

        # Pre-processing transformograph information

        if transformograph:  # if transformograph file is defined:

            ftra = open(transformograph, "r")  # opening transformograph file
            tralines = sum(1 for line in ftra) - 2  # counting number of lines
            ftra.close()  # closing file

        if not mflag:  # for single model run:

            # Writing model parameters to file

            p1file = open(temppath + "/param1.txt", "w")  # opening parameter files

            gt = [float(flowparam[0])]
            for l in range(1, lmax):
                gt.append(float(flowparam[l]))

            print(mainmapset, file=p1file)  # name of main mapset
            print("0", file=p1file)  # identifier for single model run
            print(str(nmesh), file=p1file)  # mesh spacing
            print(str(model), file=p1file)  # model
            print(len(phasesnum), file=p1file)  # number of phases
            for l in range(0, len(phasesnum)):
                print(str(phasesnum[l]), file=p1file)  # phases
            if aflag or profile:
                print("1", file=p1file)  # control for additional output
            else:
                print("0", file=p1file)
            print(gravity, file=p1file)  # gravity
            print(limiter, file=p1file)  # numerical limiter

            print(controls[0], file=p1file)  # control for correction of flow height
            print(controls[1], file=p1file)  # control for diffusion control
            print(controls[2], file=p1file)  # control for surface control
            print(controls[3], file=p1file)  # control for entrainment
            print(controls[4], file=p1file)  # control for deposition and stopping
            print(controls[5], file=p1file)  # control for dynamic variation of friction

            print(str(releasemass), file=p1file)  # control for (PHASE 1) release height
            if model == 7:
                print(
                    str(releasemass), file=p1file
                )  # control for PHASE 2 release height
            if model == 7:
                print(
                    str(releasemass), file=p1file
                )  # control for PHASE 3 release height
            print(
                str(releasevelocity), file=p1file
            )  # control for (PHASE 1) release velocity
            if model == 7:
                print(
                    str(releasevelocity2), file=p1file
                )  # control for PHASE 2 release velocity
            if model == 7:
                print(
                    str(releasevelocity3), file=p1file
                )  # control for PHASE 3 release velocity
            print(
                str(entrainment), file=p1file
            )  # control for maximum height of (PHASE 1) entrainment
            if model == 7:
                print(
                    str(entrainment), file=p1file
                )  # control for maximum height of PHASE 2 entrainment
            if model == 7:
                print(
                    str(entrainment), file=p1file
                )  # control for maximum height of PHASE 3 entrainment
            if zones:
                print("1", file=p1file)  # control for zones map
            else:
                print("0", file=p1file)
            if centr:
                print("1", file=p1file)  # control for entrainment coefficient map
            else:
                print("0", file=p1file)
            if cvshear:
                print("1", file=p1file)  # control for shear velocity coefficient map
            else:
                print("0", file=p1file)
            if phi1:
                print(
                    "1", file=p1file
                )  # control for internal friction angle of PHASE 1 map
            else:
                print("0", file=p1file)
            if phi2:
                print(
                    "1", file=p1file
                )  # control for internal friction angle of PHASE 2 map
            else:
                print("0", file=p1file)
            if phi3:
                print(
                    "1", file=p1file
                )  # control for internal friction angle of PHASE 3 map
            else:
                print("0", file=p1file)
            if deltab:
                print(
                    "1", file=p1file
                )  # control for basal friction difference map
            else:
                print("0", file=p1file)
            if tufri:
                print(
                    "1", file=p1file
                )  # control for turbulent friction coefficient map
            else:
                print("0", file=p1file)
            if delta1:
                print(
                    "1", file=p1file
                )  # control for basal friction angle of PHASE 1 map
            else:
                print("0", file=p1file)
            if delta2:
                print(
                    "1", file=p1file
                )  # control for basal friction angle of PHASE 2 map
            else:
                print("0", file=p1file)
            if delta3:
                print(
                    "1", file=p1file
                )  # control for basal friction angle of PHASE 3 map
            else:
                print("0", file=p1file)
            if ny1:
                print("1", file=p1file)  # control for viscosity of (PHASE 1) map
            else:
                print("0", file=p1file)
            if ny2:
                print("1", file=p1file)  # control for viscosity of PHASE 2 map
            else:
                print("0", file=p1file)
            if ny3:
                print("1", file=p1file)  # control for viscosity of PHASE 3 map
            else:
                print("0", file=p1file)
            if ambdrag:
                print("1", file=p1file)  # control for ambient drag coefficient map
            else:
                print("0", file=p1file)
            if flufri:
                print("1", file=p1file)  # control for PHASE 2 friction coefficient map
            else:
                print("0", file=p1file)
            if ctrans12:
                print(
                    "1", file=p1file
                )  # control for PHASE 1-PHASE 2 transformation coefficient map
            else:
                print("0", file=p1file)
            if ctrans13:
                print(
                    "1", file=p1file
                )  # control for PHASE 1-PHASE 3 transformation coefficient map
            else:
                print("0", file=p1file)
            if ctrans23:
                print(
                    "1", file=p1file
                )  # control for PHASE 2-PHASE 3 transformation coefficient map
            else:
                print("0", file=p1file)
            if trelease:
                print("1", file=p1file)  # control for release time map
            else:
                print("0", file=p1file)
            if trelstop:
                print("1", file=p1file)  # control for release hydrograph stop time map
            else:
                print("0", file=p1file)
            if stoptime:
                print("1", file=p1file)  # control for stopping time map
            else:
                print("0", file=p1file)
            if tslide:
                print("1", file=p1file)  # control for time of initial sliding map
            else:
                print("0", file=p1file)
            if hydrograph or hydrocoords:
                print("1", file=p1file)  # control for hydrograph
            else:
                print("0", file=p1file)
            if frictiograph:
                print("1", file=p1file)  # control for frictiograph
            else:
                print("0", file=p1file)
            if transformograph:
                print("1", file=p1file)  # control for transformograph
            else:
                print("0", file=p1file)
            print(elevation, file=p1file)  # name of elevation map
            if model <= 3 and releasemass == 1:
                print(pf + "_hflow0000", file=p1file)  # name of release height map
            elif model == 7 and releasemass == 1:
                print(
                    pf + "_hflow10000", file=p1file
                )  # name of PHASE 1 release height map
            if model == 7 and releasemass == 1:
                print(
                    pf + "_hflow20000", file=p1file
                )  # name of PHASE 2 release height map
            if model == 7 and releasemass == 1:
                print(
                    pf + "_hflow30000", file=p1file
                )  # name of PHASE 3 release height map
            if model <= 3 and releasevelocity == 1:
                print(
                    pf + "_vflowx0000", file=p1file
                )  # name of release velocity in x direction map
                print(
                    pf + "_vflowy0000", file=p1file
                )  # name of release velocity in y direction map
            elif model == 7 and releasevelocity == 1:
                print(
                    pf + "_vflowx10000", file=p1file
                )  # name of PHASE 1 release velocity in x direction map
                print(
                    pf + "_vflowy10000", file=p1file
                )  # name of PHASE 1 release velocity in y direction map
            if model == 7 and releasevelocity2 == 1:
                print(
                    pf + "_vflowx20000", file=p1file
                )  # name of PHASE 2 release velocity in x direction map
                print(
                    pf + "_vflowy20000", file=p1file
                )  # name of PHASE 2 release velocity in y direction map
            if model == 7 and releasevelocity3 == 1:
                print(
                    pf + "_vflowx30000", file=p1file
                )  # name of PHASE 3 release velocity in x direction map
                print(
                    pf + "_vflowy30000", file=p1file
                )  # name of PHASE 3 release velocity in y direction map
            if model <= 3 and entrainment == 1:
                print(
                    pf + "_hentrmax", file=p1file
                )  # name of maximum height of entrainment map
            elif model == 7 and entrainment == 1:
                print(
                    pf + "_hentrmax1", file=p1file
                )  # name of maximum height of PHASE 1 entrainment map
            if model == 7 and entrainment == 1:
                print(
                    pf + "_hentrmax2", file=p1file
                )  # name of maximum height of PHASE 2 entrainment map
            if model == 7 and entrainment == 1:
                print(
                    pf + "_hentrmax3", file=p1file
                )  # name of maximum height of PHASE 3 entrainment map
            if zones:
                print(zones, file=p1file)  # name of zones map
            if centr:
                print(centr, file=p1file)  # name of entrainment coefficient map
            if cvshear:
                print(cvshear, file=p1file)  # name of shear velocity coefficient map
            if phi1:
                print(
                    phi1, file=p1file
                )  # name of internal friction angle of PHASE 1 map
            if phi2:
                print(
                    phi2, file=p1file
                )  # name of internal friction angle of PHASE 2 map
            if phi3:
                print(
                    phi3, file=p1file
                )  # name of internal friction angle of PHASE 3 map
            if deltab:
                print(
                    deltab, file=p1file
                )  # name of basal friction difference map
            if tufri:
                print(tufri, file=p1file)  # name of turbulent friction coefficient map
            if delta1:
                print(
                    delta1, file=p1file
                )  # name of basal friction angle of PHASE 1 map
            if delta2:
                print(
                    delta2, file=p1file
                )  # name of basal friction angle of PHASE 2 map
            if delta3:
                print(
                    delta3, file=p1file
                )  # name of basal friction angle of PHASE 3 map
            if ny1:
                print(ny1, file=p1file)  # name of viscosity of (PHASE 1) map
            if ny2:
                print(ny2, file=p1file)  # name of viscosity of PHASE 2 map
            if ny3:
                print(ny3, file=p1file)  # name of viscosity of PHASE 3 map
            if ambdrag:
                print(ambdrag, file=p1file)  # name of ambient drag coefficient map
            if flufri:
                print(flufri, file=p1file)  # name of PHASE 2 friction coefficient map
            if ctrans12:
                print(
                    ctrans12, file=p1file
                )  # name of PHASE 1-PHASE 2 transformation coefficient map
            if ctrans13:
                print(
                    ctrans13, file=p1file
                )  # name of PHASE 1-PHASE 3 transformation coefficient map
            if ctrans23:
                print(
                    ctrans23, file=p1file
                )  # name of PHASE 2-PHASE 3 transformation coefficient map
            if trelease:
                print(trelease, file=p1file)  # name of release time map
            if trelstop:
                print(trelstop, file=p1file)  # name of release hydrograph stop time map
            if stoptime:
                print(stoptime, file=p1file)  # name of stopping time map
            if tslide:
                print(tslide, file=p1file)  # name of time of initial sliding map

            if hydrograph or hydrocoords:
                print(hyd_nin, file=p1file)  # number of input hydrographs
                print(hyd_nout, file=p1file)  # number of output hydrographs
                for hydi in range(0, hyd_nin):
                    print(hydrograph[hydi], file=p1file)  # path to hydrograph file
                    print(
                        str(hydlines[hydi]), file=p1file
                    )  # number of hydrograph time steps
                for hydi in range(0, hyd_nin + hyd_nout):
                    print(
                        hydrocoords[4 * hydi], file=p1file
                    )  # x coordinate of hydrograph
                    print(
                        hydrocoords[4 * hydi + 1], file=p1file
                    )  # y coordinate of hydrograph
                    print(
                        hydrocoords[4 * hydi + 2], file=p1file
                    )  # length of hydrograph profile
                    print(
                        hydrocoords[4 * hydi + 3], file=p1file
                    )  # aspect of hydrograph profile
            else:
                hyd_nin = 0
                hyd_nout = 0

            if frictiograph:
                print(frictiograph, file=p1file)  # path to frictiograph file
                print(str(frilines), file=p1file)  # number of frictiograph time steps

            if transformograph:
                print(transformograph, file=p1file)  # path to transformograph file
                print(
                    str(tralines), file=p1file
                )  # number of transformograph time steps

            print(lmax, file=p1file)  # number of flow parameters
            for l in range(0, lmax):
                print(round(gt[l], 10), file=p1file)  # flow parameters

            print(
                thresholdsc, file=p1file
            )  # threshold of flow height (for computation)
            print(thresholds[0], file=p1file)  # threshold of flow height (for display)
            print(thresholds[1], file=p1file)  # threshold of flow kinetic energy
            print(thresholds[2], file=p1file)  # threshold of flow pressure
            print(tint, file=p1file)  # time for writing output
            print(tstop, file=p1file)  # process duration at which to stop
            print(slomo, file=p1file)  # factor for slow motion

            print(slidepar[0], file=p1file)  # search radius for initial sliding
            print(
                slidepar[1], file=p1file
            )  # exponent for weighting for initial sliding
            print(slidepar[2], file=p1file)  # coefficient for deformation

            print(cfl[0], file=p1file)  # cfl criterion
            print(cfl[1], file=p1file)  # maximum length of time step

            print("%s_" % pf, file=p1file)  # prefix
            print(
                "%s_results/%s_ascii/" % (pf, pf), file=p1file
            )  # path and prefix for storing output maps
            print(
                "%s_results/%s_files/" % (pf, pf), file=p1file
            )  # path and prefix for storing output files

            p1file.close()  # closing parameter file

            shutil.copyfile(
                temppath + "/param1.txt",
                pf + "_results/" + pf + "_files/" + pf + "_param.txt",
            )  # copying parameter file to results folder

            # Routing flow

            print("2. ROUTING FLOW.")

            start = time.time()  # storing time (start of main computation)

            os.environ["XINT"] = "1"  # exporting id of model run
            grass.run_command("r.avaflow.main")  # executing r.avaflow

            # Finalizing computation

            ftimesteps = open(
                pf + "_results/" + pf + "_files/" + pf + "_nout1.txt", "r"
            )  # opening file with number of time steps and success
            ntimesteps = ftimesteps.readline()  # reading number of time steps
            ntimesteps = int(ntimesteps.replace("\n", ""))  # removing newline
            csuccess = (
                ftimesteps.readline()
            )  # reading control for success of simulation
            csuccess = int(csuccess.replace("\n", ""))  # removing newline
            maxvflow = ftimesteps.readline()  # reading maximum flow velocity
            maxvflow = float(maxvflow.replace("\n", ""))  # removing newline
            vol_entr = ftimesteps.readline()  # reading entrained PHASE 1 volume
            vol_entr = float(vol_entr.replace("\n", ""))  # removing newline
            vol_entr2 = ftimesteps.readline()  # reading entrained PHASE 2 volume
            vol_entr2 = float(vol_entr2.replace("\n", ""))  # removing newline
            vol_entr3 = ftimesteps.readline()  # reading entrained PHASE 3 volume
            vol_entr3 = float(vol_entr3.replace("\n", ""))  # removing newline
            ftimesteps.close()  # closing file with number of time steps

            if aflag:
                mstringa = ["_basechange", "_tflow", "_pflow", "_vflow"]
                acount = 4
            else:
                mstringa = ["_basechange"]
                acount = 1

            for j in range(0, acount):  # loop over impact parameters:

                if model <= 3:

                    grass.mapcalc(
                        '"%s%s0000"=0' % (pf, mstringa[j]), overwrite=True, quiet=True
                    )  # total map

                    grass.run_command(
                        "r.out.gdal",
                        input=pf + mstringa[j] + "0000",
                        output=ascpath + pf + mstringa[j] + "0000.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting map to ascii

                    corrasc(ascpath + pf + mstringa[j] + "0000")

                elif model == 7:

                    grass.mapcalc(
                        '"%s%s10000"=0' % (pf, mstringa[j]), overwrite=True, quiet=True
                    )
                    grass.mapcalc(
                        '"%s%s20000"=0' % (pf, mstringa[j]), overwrite=True, quiet=True
                    )

                    grass.run_command(
                        "r.out.gdal",
                        input=pf + mstringa[j] + "10000",
                        output=ascpath + pf + mstringa[j] + "10000.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting PHASE 1 map to ascii
                    grass.run_command(
                        "r.out.gdal",
                        input=pf + mstringa[j] + "20000",
                        output=ascpath + pf + mstringa[j] + "20000.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting PHASE 2 map to ascii

                    corrasc(ascpath + pf + mstringa[j] + "10000")
                    corrasc(ascpath + pf + mstringa[j] + "20000")

                    grass.mapcalc(
                        '"%s%s30000"=0' % (pf, mstringa[j]), overwrite=True, quiet=True
                    )
                    grass.mapcalc(
                        '"%s%s0000"="%s%s10000"+"%s%s20000"+"%s%s30000"'
                        % (
                            pf,
                            mstringa[j],
                            pf,
                            mstringa[j],
                            pf,
                            mstringa[j],
                            pf,
                            mstringa[j],
                        ),
                        overwrite=True,
                        quiet=True,
                    )

                    grass.run_command(
                        "r.out.gdal",
                        input=pf + mstringa[j] + "30000",
                        output=ascpath + pf + mstringa[j] + "30000.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting PHASE 3 map to ascii
                    grass.run_command(
                        "r.out.gdal",
                        input=pf + mstringa[j] + "0000",
                        output=ascpath + pf + mstringa[j] + "0000.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting total map to ascii

                    corrasc(ascpath + pf + mstringa[j] + "30000")
                    corrasc(ascpath + pf + mstringa[j] + "0000")

            minval = grass.raster_info(pf + "_basechange_fin")[
                "min"
            ]  # minimum value of basal change
            maxval = grass.raster_info(pf + "_basechange_fin")[
                "max"
            ]  # maximum value of basal change

            if minval == 0 and maxval == 0:
                basechange = 0
            else:
                basechange = 1

            if basechange == 1:
                grass.mapcalc(
                    '"%s_hflow_dep"=if("%s_basechange_fin">0,"%s_basechange_fin",0)'
                    % (pf, pf, pf),
                    overwrite=True,
                    quiet=True,
                )
            else:
                grass.mapcalc(
                    '"%s_hflow_dep"="%s_hflow_fin"' % (pf, pf),
                    overwrite=True,
                    quiet=True,
                )  # simulated height of deposition map

            grass.mapcalc(
                '"%s_elev_mod"="%s_elev"+"%s_basechange_fin"' % (pf, pf, pf),
                overwrite=True,
                quiet=True,
            )  # updated elevation map

            grass.run_command(
                "r.out.gdal",
                input=pf + "_elev_mod",
                output=ascpath + pf + "_elev_mod.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting updated elevation map to ascii
            corrasc(ascpath + pf + "_elev_mod")

            # Evaluating result

            if profile or impdef == "1" or depdef == "1" or not controls[3] == "0":

                fval = open(
                    pf + "_results/" + pf + "_files/" + pf + "_evaluation.txt", "w"
                )  # opening file for evaluation parameters for writing
                print("Criterion\tShortcut\tValue", file=fval)  # writing file header

            if model <= 3 and not controls[3] == "0":

                print(
                    "Entrained volume\tVE\t" + str(vol_entr) + " cubic metres",
                    file=fval,
                )

            elif model == 7 and not controls[3] == "0":

                print(
                    "Entrained PHASE 1 volume\tVE1\t" + str(vol_entr) + " cubic metres",
                    file=fval,
                )
                print(
                    "Entrained PHASE 2 volume\tVE2\t"
                    + str(vol_entr2)
                    + " cubic metres",
                    file=fval,
                )
                print(
                    "Entrained PHASE 3 volume\tVE3\t"
                    + str(vol_entr3)
                    + " cubic metres",
                    file=fval,
                )

            if profile or impdef == "1" or depdef == "1" or not controls[3] == "0":
                print(" ", file=fval)

            if profile:  # if profile is defined:

                grass.run_command(
                    "r.profile",
                    input=elevation,
                    output=temppath + "/felev.txt",
                    coordinates=profile,
                    null="-9999",
                    quiet=True,
                    overwrite=True,
                )
                # profile of elevation

                fprof = open(
                    temppath + "/felev.txt", "r"
                )  # opening profile of elevation for reading
                flines = sum(1 for line in fprof)  # counting number of lines
                fprof.close()  # closing file

                flowline_start = (
                    0.0
                )  # setting initial value for start to start point of profile
                if releasemass == 1:  # if release height is defined:

                    grass.run_command(
                        "r.profile",
                        input=pf + "_hflow0000",
                        output=temppath + "/fhrelease.txt",
                        coordinates=profile,
                        null="-9999",
                        quiet=True,
                        overwrite=True,
                    )
                    # profile of release height

                    fprof = open(
                        temppath + "/fhrelease.txt"
                    )  # again, opening profile file for reading
                    for j in range(1, flines):  # loop over all lines:
                        fline = fprof.readline()  # reading and splitting line
                        fline = list(map(str, fline.split(" ")))
                        if float(fline[2].replace("\n", "")) > 0:
                            flowline_start = float(
                                fline[1]
                            )  # identifying profile position at start of flow
                            break
                    fprof.close()  # closing file

                if (
                    impdef == "1"
                ):  # identifying profile position at terminus of observed flow:

                    mstringj = 0  # initializing counter
                    grass.run_command(
                        "r.profile",
                        input=pf + "_impactarea",
                        output=temppath + "/fimp.txt",
                        coordinates=profile,
                        null="-9999",
                        quiet=True,
                        overwrite=True,
                    )  # profile of map

                    flowline_stopimp = None
                    fprof = open(temppath + "/fimp.txt")  # opening profile for reading
                    for j in range(1, flines):  # loop over all lines:
                        fline = fprof.readline()  # reading and splitting line
                        fline = list(map(str, fline.split(" ")))
                        if float(fline[2].replace("\n", "")) == 1 or (
                            not flowline_stopimp and j == flines - 1
                        ):
                            flowline_stopimp = float(
                                fline[1]
                            )  # identifying profile position at terminus of observed flow
                    fprof.close()  # closing file
                    maxlength_imp = flowline_stopimp - flowline_start  # flow length
                    print(
                        "Flow length to terminus of observed impact area\tLio\t%s m"
                        % maxlength_imp,
                        file=fval,
                    )

                if (
                    depdef == "1"
                ):  # identifying profile position at terminus of observed flow:

                    mstringj = 0  # initializing counter
                    grass.run_command(
                        "r.profile",
                        input=pf + "_hdeposit",
                        output=temppath + "/fdep.txt",
                        coordinates=profile,
                        null="-9999",
                        quiet=True,
                        overwrite=True,
                    )  # profile of map

                    flowline_stopdep = None
                    fprof = open(temppath + "/fdep.txt")  # opening profile for reading
                    for j in range(1, flines):  # loop over all lines:
                        fline = fprof.readline()  # reading and splitting line
                        fline = list(map(str, fline.split(" ")))
                        if float(fline[2].replace("\n", "")) >= float(
                            thresholds[0]
                        ) or (not flowline_stopdep and j == flines - 1):
                            flowline_stopdep = float(
                                fline[1]
                            )  # identifying profile position at terminus of observed flow
                    fprof.close()  # closing file
                    maxlength_dep = flowline_stopdep - flowline_start  # flow length
                    print(
                        "Flow length to terminus of observed deposit\tLdo\t%s m"
                        % maxlength_dep,
                        file=fval,
                    )

                mstringj = 0  # initializing counter
                for mstringi in ["_hflow_max", "_hflow_dep"]:

                    grass.run_command(
                        "r.profile",
                        input=pf + mstringi,
                        output=temppath + "/f" + mstringi + ".txt",
                        coordinates=profile,
                        null="-9999",
                        quiet=True,
                        overwrite=True,
                    )  # profile of map

                    flowline_stop = None
                    fprof = open(
                        temppath + "/f" + mstringi + ".txt"
                    )  # opening profile for reading
                    for j in range(1, flines):  # loop over all lines:
                        fline = fprof.readline()  # reading and splitting line
                        fline = list(map(str, fline.split(" ")))
                        if float(fline[2].replace("\n", "")) >= float(
                            thresholds[mstringj]
                        ) or (not flowline_stop and j == flines - 1):
                            flowline_stop = float(
                                fline[1]
                            )  # identifying profile position at termination of flow
                    fprof.close()  # closing file
                    maxlength = flowline_stop - flowline_start  # flow length
                    if mstringi == "_hflow_dep":
                        print(
                            "Flow length to terminus of simulated deposit\tLdm\t%s m"
                            % maxlength,
                            file=fval,
                        )
                    elif mstringi == "_hflow_max":
                        print(
                            "Flow length to terminus of simulated impact area\tLim\t%s m"
                            % maxlength,
                            file=fval,
                        )
                    # writing flow length to text file
                    if impdef == "1" and mstringi == "_hflow_max":
                        flowline_diff = flowline_stop - flowline_stopimp
                        if maxlength_imp > 0:
                            lratio_imp = flowline_diff / maxlength_imp
                        else:
                            lratio = -9999
                        print(
                            "Difference between simulated and observed flow length (impact area)\tdLi\t%s m"
                            % str(flowline_diff),
                            file=fval,
                        )
                        print(
                            "Ratio between simulated and observed flow length (impact area)\trLi\t%s"
                            % str(lratio_imp),
                            file=fval,
                        )
                    elif (
                        depdef == "1" and mstringi == "_hflow_dep"
                    ):  # differences between modelled and observed flow lengths:
                        flowline_diff = flowline_stop - flowline_stopdep
                        if maxlength_dep > 0:
                            lratio_dep = flowline_diff / maxlength_dep
                        else:
                            lratio = -9999
                        print(
                            "Difference between simulated and observed flow length (deposit)\tdLd\t%s m"
                            % str(flowline_diff),
                            file=fval,
                        )
                        print(
                            "Ratio between simulated and observed flow length (deposit)\trLd\t%s"
                            % str(lratio_dep),
                            file=fval,
                        )

                    if not mstringi == "_hflow_max":
                        mstringj = mstringj + 1  # updating counter

                print("", file=fval)

            if impdef == "1" and depdef == "1":
                mstringlist = [
                    "_hflow_max",
                    "_hflow_dep",
                ]  # raster maps to be processed
            elif impdef == "1":
                mstringlist = ["_hflow_max"]
            elif depdef == "1":
                mstringlist = ["_hflow_dep"]
            else:
                mstringlist = []

            mstringj = 0
            for mstringi in mstringlist:  # loop over all maps to be processed:

                if mstringi == "_hflow_max":
                    mstringref = "_impactarea"
                    grass.mapcalc(
                        '"_obsbin"=if("%s%s">0,1,if("%s%s"==0,0,null()))'
                        % (pf, mstringref, pf, mstringref),
                        overwrite=True,
                        quiet=True,
                    )  # processed map of observed impact area

                elif mstringi == "_hflow_dep":
                    mstringref = "_hdeposit"
                    grass.mapcalc(
                        '"_obsbin"=if("%s%s">=%s,1,if("%s%s"==0,0,null()))'
                        % (pf, mstringref, thresholds[mstringj], pf, mstringref),
                        overwrite=True,
                        quiet=True,
                    )  # processed map height of observed deposit

                grass.mapcalc(
                    '"_simbin"=if("%s%s">=%s,1,0)'
                    % (pf, mstringi, thresholds[mstringj]),
                    overwrite=True,
                    quiet=True,
                )
                # processed map of simulated impact area or deposit
                grass.mapcalc(
                    '"_eval"="_simbin"+10*"_obsbin"', overwrite=True, quiet=True
                )  # combination of observed and simulated deposits
                grass.run_command(
                    "r.stats",
                    flags="c" + "n",
                    separator="x",
                    input="_eval",
                    output="%s/evaluation.txt" % temppath,
                    quiet=True,
                    overwrite=True,
                )  # scores tp=11, tn=0, fp=1, fn=10

                valfile = open(
                    temppath + "/evaluation.txt", "r"
                )  # opening raw file with evaluation scores
                valfile_lines = sum(1 for line in valfile)  # reading number of lines
                valfile.close()  # closing file

                if (
                    valfile_lines < 4
                ):  # if number of lines is smaller than 4 (at least one of the scores is 0):
                    valfile = open(
                        temppath + "/evaluation.txt", "r"
                    )  # opening file again
                    lr = valfile.readline()  # reading first line
                    l1 = lr.find("0x")  # checking for tn score
                    if l1 > -1:
                        lr = valfile.readline()  # if found, reading next line
                    l2 = lr.find("1x")  # checking for fp score
                    if l2 > -1:
                        lr = valfile.readline()  # if found, reading next line
                    l3 = lr.find(
                        "10x"
                    )  # checking for fn score (if not found, -1 is returned)
                    if l3 > -1:
                        lr = valfile.readline()  # if found, reading next line
                    l4 = lr.find("11x")  # checking for tp score
                    valfile.close()  # closing file again
                else:
                    l1 = (
                        1
                    )  # if number of lines is 4, setting all controls to positive (any value larger than -1 would be fine)
                    l2 = 1
                    l3 = 1
                    l4 = 1

                valfile = open(
                    temppath + "/evaluation.txt", "r"
                )  # once more opening temporary file with evaluation scores
                if l1 == -1:
                    tn = float(0)  # if no tn score was identified before, setting to 0
                else:
                    tn = valfile.readline()  # else, reading line
                    tn = float(
                        tn.replace("0x", "").replace("\n", "")
                    )  # extracting tn score
                if l2 == -1:
                    fp = float(0)  # if no fp score was identified before, setting to 0
                else:
                    fp = valfile.readline()  # else, reading line
                    fp = float(
                        fp.replace("1x", "").replace("\n", "")
                    )  # extracting fp score
                if l3 == -1:
                    fn = float(0)  # if no fn score was identified before, setting to 0
                else:
                    fn = valfile.readline()  # else, reading line
                    fn = float(
                        fn.replace("10x", "").replace("\n", "")
                    )  # extracting fn score
                if l4 == -1:
                    tp = float(0)  # if no tp score was identified before, setting to 0
                else:
                    tp = valfile.readline()  # else, reading line
                    tp = float(
                        tp.replace("11x", "").replace("\n", "")
                    )  # extracting tp score
                valfile.close()  # closing file for the last time

                tn = 5 * (tp + fn) - fp  # normalizing tn score
                if tp + fn > 0:
                    foc = (tp + fp) / (tp + fn)  # conservativeness factor
                else:
                    foc = -9999

                if tp + fn > 0 and tn + fp > 0:

                    total = tp + tn + fp + fn  # sum of scores

                    tpr = 100 * tp / total  # tp rate in percent
                    tnr = 100 * tn / total  # tn rate in percent
                    fpr = 100 * fp / total  # fp ratte in percent
                    fnr = 100 * fn / total  # fn rate in percent
                    if (tp + fn) > 0:
                        tlr = (
                            100 * tp / (tp + fn)
                        )  # rate of true positive predictions out of all positive observations
                    else:
                        tlr = -9999
                    if (tn + fp) > 0:
                        tnlr = (
                            100 * tn / (tn + fp)
                        )  # rate of true negative predictions out of all negative observations
                    else:
                        tnlr = -9999
                    mpr = tpr + fpr  # total rate of positive predictions
                    mnr = tnr + fnr  # total rate of negative predictions
                    tr = tnr + tpr  # total rate of true predictions
                    fr = fnr + fpr  # total rate of false predictions
                    opr = tpr + fnr  # rate of positive observations
                    onr = tnr + fpr  # rate of negative observations

                    if not tpr + fpr + fnr == 0:
                        csi = tpr / (tpr + fpr + fnr)  # critical success index
                    else:
                        csi = -9999
                    if not (tpr + fnr) * (fnr + tnr) + (tpr + fpr) * (fpr + tnr) == 0:
                        hss = (2 * (tpr * tnr) - (fpr * fnr)) / (
                            (tpr + fnr) * (fnr + tnr) + (tpr + fpr) * (fpr + tnr)
                        )  # Heidke skill score
                    else:
                        hss = -9999

                    tprr = tpr / opr  # rate of true positive rate
                    fprr = fpr / onr  # rate of false positive rate
                    aucroc = (
                        (1 - fprr) * tprr
                        + fprr * tprr / 2
                        + (1 - fprr) * (1 - tprr) / 2
                    )  # area under the roc curve
                    d2pc = pow(
                        pow(1 - tprr, 2) + pow(fprr, 2), 0.5
                    )  # distance to perfect classification

                else:
                    tpr = -9999
                    tnr = -9999
                    fpr = -9999
                    fnr = -9999
                    tlr = -9999
                    tnlr = -9999
                    mpr = -9999
                    mnr = -9999
                    tr = -9999
                    fr = -9999
                    opr = -9999
                    onr = -9999
                    csi = -9999
                    hss = -9999
                    aucroc = -9999
                    d2pc = -9999

                if profile:
                    if mstringi == "_hflow_max":
                        lratio = lratio_imp
                    else:
                        lratio = lratio_dep
                else:
                    lratio = -1
                if lratio == -1:
                    lratio = -1.1
                if foc == 0:
                    foc = -1

                spi = (
                    max(0, min(0.2, 0.2 * csi))
                    + max(0, min(0.2, 0.2 * hss))
                    + max(0, min(0.2, 0.2 * (1 - d2pc)))
                    + max(
                        0,
                        min(
                            0.2,
                            0.2
                            * pow(
                                math.sin(
                                    math.pi / 2 * min(lratio + 1, 1 / (lratio + 1))
                                ),
                                5,
                            ),
                        ),
                    )
                    + max(
                        0,
                        min(
                            0.2, 0.2 * pow(math.sin(math.pi / 2 * min(foc, 1 / foc)), 5)
                        ),
                    )
                )  # synthetic performance index

                if lratio == -1.1:
                    lratio = -9999
                if foc == -1:
                    foc = -9999

                if mstringi == "_hflow_dep":
                    print("Evaluation scores for observed deposit", file=fval)
                elif mstringi == "_hflow_max":
                    print("Evaluation scores for observed impact area", file=fval)

                print("Rate of positive observations\trOP\t%.1f%%" % opr, file=fval)
                print("Rate of negative observations\trON\t%.1f%%" % onr, file=fval)
                print(
                    "True positive rate\trTP\t%.1f%%" % tpr, file=fval
                )  # writing rates to file
                print("True negative rate\trTN\t%.1f%%" % tnr, file=fval)
                print("False positive rate\trFP\t%.1f%%" % fpr, file=fval)
                print("False negative rate\trFN\t%.1f%%" % fnr, file=fval)
                print(
                    "True positive rate out of all positive observations\trTP/rOP\t%.1f%%"
                    % tlr,
                    file=fval,
                )
                print(
                    "True negative rate out of all negative observations\trTN/rON\t%.1f%%"
                    % tnlr,
                    file=fval,
                )
                print(
                    "Rate of total positive predictions\trMP\t%.1f%%" % mpr, file=fval
                )
                print(
                    "Rate of total negative predictions\trMN\t%.1f%%" % mnr, file=fval
                )
                print("Rate of true predictions\trT\t%.1f%%" % tr, file=fval)
                print("Rate of false predictions\trF\t%.1f%%" % fr, file=fval)
                print("Critical success index\tCSI\t%.3f" % csi, file=fval)
                print("Heidke Skill Score\tHSS\t%.3f" % hss, file=fval)
                print("Area under the ROC curve\tAUROC\t%.3f" % aucroc, file=fval)
                print(
                    "Distance to perfect classification\tD2PC\t%.3f" % d2pc, file=fval
                )
                print("Factor of conservativeness\tFoC\t%.3f" % foc, file=fval)
                print("Synthetic performance index\tSPI\t%.3f" % spi, file=fval)
                print("", file=fval)

                if not mstringi == "_hflow_max":
                    mstringj = mstringj + 1  # updating counter

            if profile or impdef == "1" or depdef == "1":
                fval.close()  # closing file

            # Creating profile file

            if profile:  # if profile is defined:

                jprof = 0
                jheader = ["elev"]
                grass.run_command(
                    "r.profile",
                    input=pf + "_elev",
                    output=temppath + "/fprof" + str(jprof) + ".txt",
                    coordinates=profile,
                    null="-9999",
                    quiet=True,
                    overwrite=True,
                )  # elevation

                if hdeposit:

                    jprof = jprof + 1
                    jheader.append("hdeposit")
                    grass.run_command(
                        "r.profile",
                        input=pf + "_hdeposit",
                        output=temppath + "/fprof" + str(jprof) + ".txt",
                        coordinates=profile,
                        null="-9999",
                        quiet=True,
                        overwrite=True,
                    )  # height of observed deposit

                for step in range(0, ntimesteps + 1):  # loop over all time steps:

                    if step < 10:
                        fill = "000" + str(step)  # formatting time step string
                    elif step < 100:
                        fill = "00" + str(step)
                    elif step < 1000:
                        fill = "0" + str(step)
                    else:
                        fill = str(step)

                    for nprof in (
                        "hflow",
                        "vflow",
                        "tflow",
                        "pflow",
                        "basechange",
                    ):  # loop over all parameters:

                        if model <= 3:

                            jprof = jprof + 1
                            jheader.append(nprof + fill)
                            grass.run_command(
                                "r.profile",
                                input=pf + "_" + nprof + fill,
                                output=temppath + "/fprof" + str(jprof) + ".txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of parameter

                        elif model == 7:

                            jprof = jprof + 1
                            jheader.append(nprof + "1" + fill)
                            grass.run_command(
                                "r.profile",
                                input=pf + "_" + nprof + "1" + fill,
                                output=temppath + "/fprof" + str(jprof) + ".txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of PHASE 1 parameter

                            jprof = jprof + 1
                            jheader.append(nprof + "2" + fill)
                            grass.run_command(
                                "r.profile",
                                input=pf + "_" + nprof + "2" + fill,
                                output=temppath + "/fprof" + str(jprof) + ".txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of PHASE 2 parameter

                            jprof = jprof + 1
                            jheader.append(nprof + "3" + fill)
                            grass.run_command(
                                "r.profile",
                                input=pf + "_" + nprof + "3" + fill,
                                output=temppath + "/fprof" + str(jprof) + ".txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of PHASE 3 parameter

                fprof2 = open(
                    pf + "_results/" + pf + "_files/" + pf + "_profile.txt", "w"
                )  # opening combined profile file for writing
                strprof = []  # string of file content
                for j in range(0, jprof + 1):  # loop over all parameters:
                    fprof = open(
                        temppath + "/fprof" + str(j) + ".txt", "r"
                    )  # opening profile file for reading
                    strprof.append(
                        fprof.readlines()
                    )  # writing entire file content to string
                    fprof.close()  # closing file

                    if j == 0:
                        fprof2.write(jheader[j])
                    else:
                        fprof2.write("\t" + jheader[j])
                fprof2.write("\n")

                for row in range(len(strprof[0])):  # loop over all time steps:
                    strprof2 = (
                        strprof[0][row].strip().split(" ")[0]
                    )  # coordinate values
                    for j in range(0, jprof + 1):  # loop over all parameters:
                        strprof2 = (
                            strprof2 + "\t" + strprof[j][row].strip().split(" ")[1]
                        )  # parameter values
                    fprof2.write(strprof2 + "\n")

                fprof2.close()  # closing file

            # Creating control point information file

            if ctrlpoints:  # if control points are defined:

                fctrlpts = open(
                    pf + "_results/" + pf + "_files/" + pf + "_ctrlpoints.txt", "w"
                )  # opening control points information file for writing
                fctrlpts.write(
                    "id\tt\tx\ty\thflow1\thflow2\thflow3\thentr1\thentr2\thentr3\n"
                )  # writing header to file

                for ictrl in range(0, lctrlpts):  # loop over all control points:

                    for step in range(0, ntimesteps + 1):  # loop over all time steps:

                        if step < 10:
                            fill = "000" + str(step)  # formatting time step string
                        elif step < 100:
                            fill = "00" + str(step)
                        elif step < 1000:
                            fill = "0" + str(step)
                        else:
                            fill = str(step)

                        if model <= 3:

                            ctrlstring1 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                                "\n", ""
                            )

                            ctrlstring2 = "0"
                            ctrlstring3 = "0"

                            ctrlstring4 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring4 = (
                                ctrlstring4.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring5 = "0"
                            ctrlstring6 = "0"

                        elif model <= 3:

                            ctrlstring1 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow1" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                                "\n", ""
                            )

                            ctrlstring2 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow2" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring2 = (
                                ctrlstring2.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring3 = "0"

                            ctrlstring4 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange1" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring4 = (
                                ctrlstring4.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring5 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange2" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring5 = (
                                ctrlstring5.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring6 = "0"

                        else:

                            ctrlstring1 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow1" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                                "\n", ""
                            )

                            ctrlstring2 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow2" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring2 = (
                                ctrlstring2.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring3 = grass.read_command(
                                "r.what",
                                map=pf + "_hflow3" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring3 = (
                                ctrlstring3.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring4 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange1" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring4 = (
                                ctrlstring4.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring5 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange2" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring5 = (
                                ctrlstring5.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                            ctrlstring6 = grass.read_command(
                                "r.what",
                                map=pf + "_basechange3" + fill,
                                separator="\t",
                                coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                                overwrite=True,
                                quiet=True,
                            )
                            ctrlstring6 = (
                                ctrlstring6.replace("\t\t", "\t")
                                .replace("\n", "")
                                .replace(str(xctrl[ictrl]) + "\t", "")
                                .replace(str(yctrl[ictrl]) + "\t", "")
                            )

                        fctrlpts.write(
                            "C"
                            + str(ictrl + 1)
                            + "\t"
                            + str(step * float(tint))
                            + "\t"
                            + ctrlstring1
                            + "\t"
                            + ctrlstring2
                            + "\t"
                            + ctrlstring3
                            + "\t"
                            + ctrlstring4
                            + "\t"
                            + ctrlstring5
                            + "\t"
                            + ctrlstring6
                            + "\n"
                        )

                    if model <= 3:

                        ctrlstring1 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                            "\n", ""
                        )

                        ctrlstring2 = "0"
                        ctrlstring3 = "0"

                    elif model <= 3:

                        ctrlstring1 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow1_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                            "\n", ""
                        )

                        ctrlstring2 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow2_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring2 = (
                            ctrlstring2.replace("\t\t", "\t")
                            .replace("\n", "")
                            .replace(str(xctrl[ictrl]) + "\t", "")
                            .replace(str(yctrl[ictrl]) + "\t", "")
                        )

                        ctrlstring3 = "0"

                    else:

                        ctrlstring1 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow1_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring1 = ctrlstring1.replace("\t\t", "\t").replace(
                            "\n", ""
                        )

                        ctrlstring2 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow2_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring2 = (
                            ctrlstring2.replace("\t\t", "\t")
                            .replace("\n", "")
                            .replace(str(xctrl[ictrl]) + "\t", "")
                            .replace(str(yctrl[ictrl]) + "\t", "")
                        )

                        ctrlstring3 = grass.read_command(
                            "r.what",
                            map=pf + "_hflow3_max",
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )
                        ctrlstring3 = (
                            ctrlstring3.replace("\t\t", "\t")
                            .replace("\n", "")
                            .replace(str(xctrl[ictrl]) + "\t", "")
                            .replace(str(yctrl[ictrl]) + "\t", "")
                        )

                    fctrlpts.write(
                        "C"
                        + str(ictrl + 1)
                        + "\tmax\t"
                        + ctrlstring1
                        + "\t"
                        + ctrlstring2
                        + "\t"
                        + ctrlstring3
                        + "\t"
                        + ctrlstring4
                        + "\t"
                        + ctrlstring5
                        + "\t"
                        + ctrlstring6
                        + "\n"
                    )
                    # writing point number and coordinates to file

                fctrlpts.close()  # closing file

            stop = time.time()  # storing time (end of main computation)
            comptime = stop - start  # storing computational time in seconds

            timefile = open(filepath + pf + "_time.txt", "w")
            timefile.write(str(comptime))  # writing computational time to file
            timefile.close()

        # Computing impact indicator indices and writing input files for AIMEC

        else:  # for multiple model runs:

            if int(sampling) > 0:
                nruns = int(sampling)  # number of model runs
            elif int(sampling) < 0:
                nruns = int(sampling) * -1
            ncores = int(cores)  # number of processors to be used

            os.mkdir(
                pf + "_results/" + pf + "_aimec"
            )  # creating directory for AIMEC input
            os.mkdir(
                pf + "_results/" + pf + "_aimec/depth"
            )  # creating directory for maximum flow height maps for AIMEC input
            os.mkdir(
                pf + "_results/" + pf + "_aimec/pressure"
            )  # creating directory for maximum kinetic energy maps for AIMEC input

            fparamt = open(
                temppath + "/fparamt.txt", "w"
            )  # opening temporary input parameter file for writing
            faimecp = open(
                pf + "_results/" + pf + "_aimec/" + pf + "_aimec_params.txt", "w"
            )  # opening parameter input file for AIMEC for writing

            if model == 0:
                print(
                    "nrunp\tvhr\tvhem\trho1\tphi1\tdelta1\ttufri1\tcentr\tcvshear\tdeltab\tcstop\thydadd\tdecel\tfrimin\tekincoef\n",
                    end=" ",
                    file=fparamt,
                )  # writing header to temporary input parameter file
                print(
                    "nrun, vhr, vhem, rho1, phi1, delta1, tufri1, centr, cvshear, deltab, cstop, hydadd, decel, frimin, ekincoef",
                    end=" ",
                    file=faimecp,
                )  # writing header to parameter input file for AIMEC

            elif model == 1:
                print(
                    "nrunp\tvhr\tvhem\trho1\tphi1\tdelta1\tflufri\tny1\ttauy1\tcentr\tcvshear\tdeltab\tcstop\tambdrag\tchi1\try\tkpmax1\thydadd\tdecel\tfrimin\tekincoef\n",
                    end=" ",
                    file=fparamt,
                )
                print(
                    "nrun, vhr, vhem, rho1, phi1, delta1, flufri, ny1, tauy1, centr, cvshear, deltab, cstop, ambdrag, chi1, ry, kpmax1, hydadd, decel, frimin, ekincoef",
                    end=" ",
                    file=faimecp,
                )

            elif model == 7:
                print(
                    "nrunp\tvhr\trhrs\tvhem\trhes\trho1\trho2\trho3\tphi1\tdelta1\tphi2\tdelta2\tphi3\tdelta3\tflufri\tny1\ttauy1\tny2\ttauy2\tny3\ttauy3\tcentr\tcvshear\tdeltab\tthetas\tcstop\tctp1p2\tctp1p3\tctp2p3\tambdrag\tnnvm\tlvm\tnvm\tfvm\tkdrag\tmdrag\tndrag\tut\trep\tj\tchi1\tchi2\tchi3\txi1\txi2\txi3\tap2p1\tap3p1\tap3p2\try\tjny\tkpmaxp1\tkpmaxp2\tkpmaxp3\thydadd\tdecel\tfrimin\tekincoef\tjpd\n",
                    end=" ",
                    file=fparamt,
                )
                print(
                    "nrun, vhr, rhrs, vhem, rhes, rho1, rho2, rho3, phi1, delta1, phi2, delta2, phi3, delta3, flufri, ny1, tauy1, ny2, tauy2, ny3, tauy3, centr, cvshear, deltab, thetas, cstop, ctp1p2, ctp1p3, ctp2p3, ambdrag, nnvm, lvm, nvm, fvm, kdrag, mdrag, ndrag, ut, rep, j, chi1, chi2, chi3, xi1, xi2, xi3, ap2p1, ap3p1, ap3p2, ry, jny, kpmaxp1, kpmaxp2, kpmaxp3, hydadd, decel, frimin, ekincoef, jpd",
                    end=" ",
                    file=faimecp,
                )

            if impdef == "1":
                grass.mapcalc(
                    '"_cimpactarea"=if("%s_impactarea">0,1,0)' % pf,
                    overwrite=True,
                    quiet=True,
                )  # binary map of observed impact area
                grass.run_command(
                    "r.out.gdal",
                    input="_cimpactarea",
                    output=pf + "_results/" + pf + "_aimec/" + pf + "_impactarea.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(pf + "_results/" + pf + "_aimec/" + pf + "_impactarea")

            if depdef == "1":
                grass.mapcalc(
                    '"_cdepositarea"=if("%s_hdeposit">0,1,0)' % pf,
                    overwrite=True,
                    quiet=True,
                )  # binary map of observed height of deposit
                grass.run_command(
                    "r.out.gdal",
                    input="_cdepositarea",
                    output=pf + "_results/" + pf + "_aimec/" + pf + "_depositarea.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting to ascii
                corrasc(pf + "_results/" + pf + "_aimec/" + pf + "_depositarea")

            if profile:  # creating profile file in directory for AIMEC input:

                fpaimec = open(
                    pf + "_results/" + pf + "_aimec/" + pf + "_profile.xyz", "w"
                )  # opening text file for writing
                profcoord = list(
                    map(str, profile.split(","))
                )  # splitting profile coordinate string
                proflength = len(profcoord)  # reading number of profile coordinates
                for profnum in range(0, proflength, 2):  # loop over all coordinates:
                    profilex = profcoord[profnum]  # extracting x coordinate
                    profiley = profcoord[profnum + 1]  # extracting y coordinate
                    profstring = grass.read_command(
                        "r.what",
                        map=pf + "_elev",
                        separator="\t",
                        coordinates=profilex + "," + profiley,
                        overwrite=True,
                        quiet=True,
                    )  # reading elevation value associated to profile coordinates
                    profstring = profstring.replace("\t\t", "\t").replace("\n", "")
                    print(
                        "%s" % profstring, file=fpaimec
                    )  # writing profile string to text file
                fpaimec.close()  # closing file

            # Preparing variation of input raster maps

            if model == 7 and hrelease:
                rhrls = list(
                    map(str, rhrelease1.split(","))
                )  # splitting string with ratio(s) of PHASE 1 release height
                rhrlsmin = rhrls[0]  # minimum
                if len(rhrls) == 3:  # maximum and third parameter
                    rhrlsmax = rhrls[1]
                    rhrls3 = rhrls[2]
                elif len(rhrls) == 2:
                    rhrlsmax = rhrls[1]
                elif len(rhrls) == 1:
                    rhrlsmax = rhrls[0]
                else:
                    ErrorMessage()  # else, exiting with error message

            if (hrelease or hrelease1 or hrelease2 or hrelease3) and vhrelease:
                vhrl = list(
                    map(str, vhrelease.split(","))
                )  # splitting string with variation of release height
                vhrlmin = vhrl[0]  # minimum
                if len(vhrl) == 3:  # maximum and third parameter
                    vhrlmax = vhrl[1]
                    vhrl3 = vhrl[2]
                elif len(vhrl) == 2:
                    vhrlmax = vhrl[1]
                elif len(vhrl) == 1:
                    vhrlmax = vhrl[0]
                else:
                    ErrorMessage()  # else, exiting with error message
            else:
                vhrlmin = "1"
                vhrlmax = "1"
                vhrl3 = "1"

            if model == 7 and hentrmax:
                rhems = list(
                    map(str, rhentrmax1.split(","))
                )  # splitting string with ratio(s) of maximum height of entrainment
                rhemsmin = rhems[0]  # minimum
                if len(rhems) == 3:  # maximum and third parameter
                    rhemsmax = rhems[1]
                    rhems3 = rhems[2]
                elif len(rhems) == 2:
                    rhemsmax = rhems[1]
                elif len(rhems) == 1:
                    rhemsmax = rhems[0]
                else:
                    ErrorMessage()  # else, exiting with error message

            if (hentrmax or hentrmax1 or hentrmax2 or hentrmax3) and vhentrmax:
                vhem = list(
                    map(str, vhentrmax.split(","))
                )  # splitting string with variation of maximum height of entrainment
                vhemmin = vhem[0]  # minimum
                if len(vhem) == 2:  # maximum and third parameter
                    vhemmax = vhem[1]
                    vhem3 = vhem[2]
                elif len(vhem) == 2:
                    vhemmax = vhem[1]
                elif len(vhem) == 1:
                    vhemmax = vhem[0]
                else:
                    ErrorMessage()  # else, exiting with error message
            else:
                vhemmin = "1"
                vhemmax = "1"
                vhem3 = "1"

            if sampling == "0":  # if controlled sampling is applied:

                njids = 2
                nruns0 = []
                nruns = 1

                if float(vhrlmin) == float(vhrlmax):
                    nruns0.append(1)
                else:
                    nruns0.append(int(vhrl3))
                nruns *= nruns0[0]

                if float(vhemmin) == float(vhemmax):
                    nruns0.append(1)
                else:
                    nruns0.append(int(vhem3))
                nruns *= nruns0[1]

                if model == 7 and hrelease and rhrelease1:

                    njids += 1
                    if float(rhrlsmin) == float(rhrlsmax):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(rhrls3))
                    nruns *= nruns0[2]

                if model == 7 and hentrmax and rhentrmax1:

                    njids += 1
                    if float(rhemsmin) == float(rhemsmax):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(rhems3))
                    nruns *= nruns0[njids - 1]

                for l in range(0, lmax):

                    if float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                        nruns0.append(1)
                    else:
                        nruns0.append(int(flowparam[3 * l + 2]))
                    nruns *= nruns0[l + njids]

                jids = [0]
                for l in range(1, lmax + njids):
                    jids.append(0)

            elif int(sampling) < 0:  # if OAT sampling is applied:

                nruns00 = nruns
                nruns = 0
                nvar = []
                njids = 2

                if not float(vhrlmin) == float(vhrlmax):
                    nvar.append(1)
                    nruns += nruns00
                else:
                    nvar.append(0)

                if not float(vhemmin) == float(vhemmax):
                    nvar.append(1)
                    nruns += nruns00
                else:
                    nvar.append(0)

                if model == 7 and hrelease and rhrelease1:

                    njids += 1
                    if not float(rhrlsmin) == float(rhrlsmax):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                if model == 7 and hentrmax and rhentrmax1:

                    njids += 1
                    if not float(rhemsmin) == float(rhemsmax):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                for l in range(0, lmax):
                    if not float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                        nvar.append(1)
                        nruns += nruns00
                    else:
                        nvar.append(0)

                lnrun = 0
                ltest = 0

            for jid in range(
                1, nruns + 1
            ):  # loop over predefined number of randomized parameter combinations:

                if jid < 10:
                    jfill = "00000" + str(jid)  # formatting model run string
                elif jid < 100:
                    jfill = "0000" + str(jid)
                elif jid < 1000:
                    jfill = "000" + str(jid)
                elif jid < 10000:
                    jfill = "00" + str(jid)
                elif jid < 100000:
                    jfill = "0" + str(jid)
                else:
                    jfill = str(jid)

                if sampling == "0":  # if controlled sampling is applied:

                    jidctrl = 0
                    for l in range(0, lmax + njids):

                        ll = l
                        if ll == 0 or jids[l] == 0 or jidctrl == 1:
                            jids[l] += 1
                        if jids[l] > nruns0[l]:
                            jids[l] = 1
                            jidctrl = 1
                        else:
                            jidctrl = 0

                # Varying input raster maps

                ipar = []
                if sampling == "0":  # for controlled sampling:

                    if float(vhrlmin) == float(vhrlmax):
                        vhrl = float(vhrlmin)
                    else:
                        vhrl = float(vhrlmin) + float(jids[0] - 1) / float(
                            nruns0[0] - 1
                        ) * (float(vhrlmax) - float(vhrlmin))
                        ipar.append(1)

                elif int(sampling) > 0:  # for random sampling:

                    vhrl = round(
                        random.uniform(float(vhrlmin), float(vhrlmax)), 2
                    )  # variation of release height

                else:  # if OAT sampling is applied:

                    if ltest < lmax + njids:
                        while nvar[ltest] == 0 and ltest < lmax + njids:
                            ltest += 1
                            lnrun = 0

                    if float(vhrlmin) == float(vhrlmax) or not ltest == 0:
                        vhrl = float(vhrl3)
                    else:
                        vhrl = float(vhrlmin) + float(lnrun) / float(nruns00 - 1) * (
                            float(vhrlmax) - float(vhrlmin)
                        )

                vhrl = str(vhrl)

                if sampling == "0":  # for controlled sampling:

                    if float(vhemmin) == float(vhemmax):
                        vhem = float(vhemmin)
                    else:
                        vhem = float(vhemmin) + float(jids[1] - 1) / float(
                            nruns0[1] - 1
                        ) * (float(vhemmax) - float(vhemmin))
                        if model <= 3:
                            ipar.append(2)
                        else:
                            ipar.append(3)

                elif int(sampling) > 0:  # for random sampling:

                    vhem = round(
                        random.uniform(float(vhemmin), float(vhemmax)), 2
                    )  # variation of maximum height of entrainment

                else:  # if OAT sampling is applied:

                    if float(vhemmin) == float(vhemmax) or not ltest == 1:
                        vhem = float(vhem3)
                    else:
                        vhem = float(vhemmin) + float(lnrun) / float(nruns00 - 1) * (
                            float(vhemmax) - float(vhemmin)
                        )

                vhem = str(vhem)

                nhrelease = pf + "_hflow0000" + str(jid)
                nhentrmax = pf + "_hentrmax" + str(jid)
                if (
                    model <= 3
                ):  # variables for names of input raster maps for model run:
                    nhentrmax = pf + "_hentrmax" + str(jid)
                elif model == 7:
                    nhrelease1 = pf + "_hflow10000" + str(jid)
                    nhrelease2 = pf + "_hflow20000" + str(jid)
                    nhentrmax1 = pf + "_hentrmax1" + str(jid)
                    nhentrmax2 = pf + "_hentrmax2" + str(jid)
                    nhrelease3 = pf + "_hflow30000" + str(jid)
                    nhentrmax3 = pf + "_hentrmax3" + str(jid)

                if model <= 3 and hrelease:  # if release height is given:
                    grass.mapcalc(
                        '"%s"="%s_hflow0000"*%s' % (nhrelease, pf, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # release height map
                elif model <= 3:
                    grass.mapcalc('"%s"=0' % nhrelease, overwrite=True, quiet=True)

                if model == 7 and hrelease1:  # if release height of PHASE 1 is given:
                    grass.mapcalc(
                        '"%s"="%s_hflow10000"*%s' % (nhrelease1, pf, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # release height map of PHASE 1
                elif model == 7:
                    grass.mapcalc('"%s"=0' % nhrelease1, overwrite=True, quiet=True)

                if model == 7 and hrelease2:  # if release height of PHASE 2 is given:
                    grass.mapcalc(
                        '"%s"="%s_hflow20000"*%s' % (nhrelease2, pf, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # release height map of PHASE 2
                elif model >= 3:
                    grass.mapcalc('"%s"=0' % nhrelease2, overwrite=True, quiet=True)

                if model == 7 and hrelease3:  # if release height of PHASE 3 is given:
                    grass.mapcalc(
                        '"%s"="%s_hflow30000"*%s' % (nhrelease3, pf, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # release height map of PHASE 3
                elif model == 7:
                    grass.mapcalc('"%s"=0' % nhrelease3, overwrite=True, quiet=True)

                if model == 7:
                    grass.mapcalc(
                        '"%s"="%s"+"%s"+"%s"'
                        % (nhrelease, nhrelease1, nhrelease2, nhrelease3),
                        overwrite=True,
                        quiet=True,
                    )

                if (
                    model == 7 and hentrmax and not controls[3] == "0"
                ):  # if maximum height of entrainment is given:
                    grass.mapcalc(
                        '"%s"="%s_hentrmax"*%s' % (nhentrmax, pf, vhem),
                        overwrite=True,
                        quiet=True,
                    )
                    # maximum height of entrainment map
                elif model <= 3:
                    grass.mapcalc('"%s"=0' % nhentrmax, overwrite=True, quiet=True)

                if (
                    model == 7 and hentrmax1 and not controls[3] == "0"
                ):  # if maximum height of entrainment of PHASE 1 is given:
                    grass.mapcalc(
                        '"%s"="%s_hentrmax1"*%s' % (nhentrmax1, pf, vhem),
                        overwrite=True,
                        quiet=True,
                    )
                    # maximum height of entrainment map of PHASE 1
                elif model == 7:
                    grass.mapcalc('"%s"=0' % nhentrmax1, overwrite=True, quiet=True)

                if (
                    model == 7 and hentrmax2 and not controls[3] == "0"
                ):  # if maximum height of entrainment of PHASE 2 is given:
                    grass.mapcalc(
                        '"%s"="%s_hentrmax2"*%s' % (nhentrmax2, pf, vhem),
                        overwrite=True,
                        quiet=True,
                    )
                    # maximum height of entrainment map of PHASE 2
                elif model == 7:
                    grass.mapcalc('"%s"=0' % nhentrmax2, overwrite=True, quiet=True)

                if (
                    model == 7 and hentrmax1 and not controls[3] == "0"
                ):  # if maximum height of entrainment of PHASE 3 is given:
                    grass.mapcalc(
                        '"%s"="%s_hentrmax3"*%s' % (nhentrmax3, pf, vhem),
                        overwrite=True,
                        quiet=True,
                    )
                    # maximum height of entrainment map of PHASE 3
                elif model == 7:
                    grass.mapcalc('"%s"=0' % nhentrmax3, overwrite=True, quiet=True)

                if model == 7:
                    grass.mapcalc(
                        '"%s"="%s"+"%s"+"%s"'
                        % (nhentrmax, nhentrmax1, nhentrmax2, nhentrmax3),
                        overwrite=True,
                        quiet=True,
                    )

                if (
                    model == 7 and hrelease and rhrelease1
                ):  # if release heights are defined by total height and ratio of PHASE 1:

                    if sampling == "0":  # for controlled sampling:

                        if float(rhrlsmin) == float(rhrlsmax):
                            rhrls = float(rhrlsmin)
                        else:
                            rhrls = float(rhrlsmin) + float(jids[2] - 1) / float(
                                nruns0[2] - 1
                            ) * (float(rhrlsmax) - float(rhrlsmin))
                            ipar.append(2)

                    elif int(sampling) > 0:  # for random sampling:

                        rhrls = random.uniform(
                            float(rhrlsmin), float(rhrlsmax)
                        )  # randomizing ratio of PHASE 1 release height

                    else:  # if OAT sampling is applied:

                        if float(rhrlsmin) == float(rhrlsmax) or not ltest == 2:
                            rhrls = float(rhrls3)
                        else:
                            rhrls = float(rhrlsmin) + float(lnrun) / float(
                                nruns00 - 1
                            ) * (float(rhrlsmax) - float(rhrlsmin))

                    rhrls = str(rhrls)  # constraining ratio and converting to string

                    grass.mapcalc(
                        '"%s"="%s_hflow0000"*%s*%s' % (nhrelease1, pf, rhrls, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # PHASE 1 release height map
                    grass.mapcalc(
                        '"%s"=0' % nhrelease2, overwrite=True, quiet=True
                    )  # PHASE 2 release height map
                    grass.mapcalc(
                        '"%s"="%s_hflow0000"*(1-%s)*%s' % (nhrelease3, pf, rhrls, vhrl),
                        overwrite=True,
                        quiet=True,
                    )  # PHASE 3 release height map

                elif model == 7:
                    grass.mapcalc(
                        '"%s"="%s_hflow10000"' % (nhrelease1, pf),
                        overwrite=True,
                        quiet=True,
                    )
                    grass.mapcalc(
                        '"%s"="%s_hflow20000"' % (nhrelease2, pf),
                        overwrite=True,
                        quiet=True,
                    )
                    grass.mapcalc(
                        '"%s"="%s_hflow30000"' % (nhrelease3, pf),
                        overwrite=True,
                        quiet=True,
                    )
                    rhrls = "-9999"

                if model == 7 and hentrmax and rhentrmax1:
                    # if maximum heights of entrainment are defined by total height and ratio of PHASE 1:

                    if sampling == "0":  # for controlled sampling:

                        if float(rhemsmin) == float(rhemsmax):
                            rhems = float(rhemsmin)
                        else:
                            rhems = float(rhemsmin) + float(
                                jids[njids - 1] - 1
                            ) / float(nruns0[njids - 1] - 1) * (
                                float(rhemsmax) - float(rhemsmin)
                            )
                            ipar.append(4)

                    elif int(sampling) > 0:  # for random sampling:

                        rhems = random.uniform(
                            float(rhemsmin), float(rhemsmax)
                        )  # randomizing ratio of maximum PHASE 1 height of entrainment

                    else:  # for OAT sampling:

                        if float(rhemsmin) == float(rhemsmax) or not ltest == njids - 1:
                            rhems = float(rhems3)
                        else:
                            rhems = float(rhemsmin) + float(lnrun) / float(
                                nruns00 - 1
                            ) * (float(rhemsmax) - float(rhemsmin))

                    rhems = str(rhems)  # constraining ratio and converting to string

                    grass.mapcalc(
                        '"%s"="%s_hentrmax"*%s*%s' % (nhentrmax1, pf, rhems, vhem),
                        overwrite=True,
                        quiet=True,
                    )  # maximum height of PHASE 1 entrainment map
                    grass.mapcalc(
                        '"%s"=0' % nhentrmax2, overwrite=True, quiet=True
                    )  # maximum height of PHASE 2 entrainment map
                    grass.mapcalc(
                        '"%s"="%s_hentrmax"*(1-%s)*%s' % (nhentrmax3, pf, rhems, vhem),
                        overwrite=True,
                        quiet=True,
                    )  # maximum height of PHASE 3 entrainment map

                elif model == 7:
                    if entrainment == 1:
                        grass.mapcalc(
                            '"%s"="%s_hentrmax1"' % (nhentrmax1, pf),
                            overwrite=True,
                            quiet=True,
                        )
                    if entrainment == 1:
                        grass.mapcalc(
                            '"%s"="%s_hentrmax2"' % (nhentrmax2, pf),
                            overwrite=True,
                            quiet=True,
                        )
                    rhems = "-9999"

                if model == 7:
                    if entrainment == 1:
                        grass.mapcalc(
                            '"%s"="%s_hentrmax3"' % (nhentrmax3, pf),
                            overwrite=True,
                            quiet=True,
                        )
                    rhems = "-9999"

                if model <= 3:

                    if releasemass == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhrelease,
                            output=ascpath + nhrelease + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting release height map to ascii
                    if entrainment == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhentrmax,
                            output=ascpath + nhentrmax + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting maximum height of entrainment map to ascii

                    if releasemass == 1:
                        corrasc(ascpath + nhrelease)
                    if entrainment == 1:
                        corrasc(ascpath + nhentrmax)

                elif model == 7:

                    if releasemass == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhrelease1,
                            output=ascpath + nhrelease1 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting PHASE 1 release height map to ascii
                    if releasemass == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhrelease2,
                            output=ascpath + nhrelease2 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting PHASE 2 release height map to ascii
                    if entrainment == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhentrmax1,
                            output=ascpath + nhentrmax1 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting maximum height of PHASE 1 entrainment map to ascii
                    if entrainment == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhentrmax2,
                            output=ascpath + nhentrmax2 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting maximum height of PHASE 2 entrainment map to ascii

                    if releasemass == 1:
                        corrasc(ascpath + nhrelease1)
                    if releasemass == 1:
                        corrasc(ascpath + nhrelease2)
                    if entrainment == 1:
                        corrasc(ascpath + nhentrmax1)
                    if entrainment == 1:
                        corrasc(ascpath + nhentrmax2)

                if model == 7:

                    if releasemass == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhrelease3,
                            output=ascpath + nhrelease3 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting PHASE 2 release height map to ascii
                    if entrainment == 1:
                        grass.run_command(
                            "r.out.gdal",
                            input=nhentrmax3,
                            output=ascpath + nhentrmax3 + ".asc",
                            format="AAIGrid",
                            overwrite=True,
                        )  # exporting maximum height of PHASE 2 entrainment map to ascii

                    if releasemass == 1:
                        corrasc(ascpath + nhrelease3)
                    if entrainment == 1:
                        corrasc(ascpath + nhentrmax3)

                # Writing flow parameters to file

                p1file = open(
                    temppath + "/param" + str(jid) + ".txt", "w"
                )  # opening parameter files
                gt = []

                if sampling == "0":  # if controlled sampling is applied:

                    for l in range(0, lmax):
                        if float(flowparam[3 * l]) == float(flowparam[3 * l + 1]):
                            gt.append(float(flowparam[3 * l]))
                        else:
                            gt.append(
                                float(flowparam[3 * l])
                                + float(jids[l + njids] - 1)
                                / float(nruns0[l + njids] - 1)
                                * (
                                    float(flowparam[3 * l + 1])
                                    - float(flowparam[3 * l])
                                )
                            )
                            if model <= 3:
                                ipar.append(3 + l)
                            else:
                                ipar.append(5 + l)

                elif int(sampling) > 0:  # if random sampling is applied:

                    for l in range(0, lmax):
                        gt.append(
                            random.uniform(
                                float(flowparam[2 * l]), float(flowparam[2 * l + 1])
                            )
                        )
                    # randomizing flow and basal surface parameters

                else:  # if OAT sampling is applied:

                    for l in range(0, lmax):
                        if (
                            float(flowparam[3 * l]) == float(flowparam[3 * l + 1])
                            or not ltest == l + njids
                        ):
                            gt.append(float(flowparam[3 * l + 2]))
                        else:
                            gt.append(
                                float(flowparam[3 * l])
                                + float(lnrun)
                                / float(nruns00 - 1)
                                * (
                                    float(flowparam[3 * l + 1])
                                    - float(flowparam[3 * l])
                                )
                            )

                    lnrun += 1
                    if lnrun == nruns00:
                        lnrun = 0
                        ltest += 1

                print(mainmapset, file=p1file)  # name of main mapset
                print("1", file=p1file)  # identifier for multiple model runs
                print(str(nmesh), file=p1file)  # mesh spacing
                print(model, file=p1file)  # model
                print(len(phasesnum), file=p1file)  # number of phases
                for l in range(0, len(phasesnum)):
                    print(str(phasesnum[l]), file=p1file)  # phases
                if aflag or profile:
                    print("1", file=p1file)  # control for additional output
                else:
                    print("0", file=p1file)
                print(gravity, file=p1file)  # gravity
                print(limiter, file=p1file)  # numerical limiter

                print(controls[0], file=p1file)  # control for correction of flow height
                print(controls[1], file=p1file)  # control for diffusion control
                print(controls[2], file=p1file)  # control for surface control
                print(controls[3], file=p1file)  # control for entrainment
                print(controls[4], file=p1file)  # control for deposition and stopping
                print(
                    controls[5], file=p1file
                )  # control for dynamic variation of friction

                print(
                    str(releasemass), file=p1file
                )  # control for (PHASE 1) release height
                if model == 7:
                    print(
                        str(releasemass), file=p1file
                    )  # control for PHASE 2 release height
                if model == 7:
                    print(
                        str(releasemass), file=p1file
                    )  # control for PHASE 3 release height
                print(
                    str(releasevelocity), file=p1file
                )  # control for (PHASE 1) release velocity
                if model == 7:
                    print(
                        str(releasevelocity2), file=p1file
                    )  # control for PHASE 2 release velocity
                if model == 7:
                    print(
                        str(releasevelocity3), file=p1file
                    )  # control for PHASE 3 release velocity
                print(
                    str(entrainment), file=p1file
                )  # control for maximum height of (PHASE 1) entrainment
                if model == 7:
                    print(
                        str(entrainment), file=p1file
                    )  # control for maximum height of PHASE 2 entrainment
                if model == 7:
                    print(
                        str(entrainment), file=p1file
                    )  # control for maximum height of PHASE 3 entrainment
                if zones:
                    print("1", file=p1file)  # control for zones map
                else:
                    print("0", file=p1file)
                if centr:
                    print("1", file=p1file)  # control for entrainment coefficient map
                else:
                    print("0", file=p1file)
                if cvshear:
                    print(
                        "1", file=p1file
                    )  # control for shear velocity coefficient map
                else:
                    print("0", file=p1file)
                if phi1:
                    print(
                        "1", file=p1file
                    )  # control for internal friction angle of PHASE 1 map
                else:
                    print("0", file=p1file)
                if phi2:
                    print(
                        "1", file=p1file
                    )  # control for internal friction angle of PHASE 2 map
                else:
                    print("0", file=p1file)
                if phi3:
                    print(
                        "1", file=p1file
                    )  # control for internal friction angle of PHASE 3 map
                else:
                    print("0", file=p1file)
                if deltab:
                    print(
                        "1", file=p1file
                    )  # control for basal friction difference map
                else:
                    print("0", file=p1file)
                if tufri:
                    print(
                        "1", file=p1file
                    )  # control for turbulent friction coefficient map
                else:
                    print("0", file=p1file)
                if delta1:
                    print(
                        "1", file=p1file
                    )  # control for basal friction angle of PHASE 1 map
                else:
                    print("0", file=p1file)
                if delta2:
                    print(
                        "1", file=p1file
                    )  # control for basal friction angle of PHASE 2 map
                else:
                    print("0", file=p1file)
                if delta3:
                    print(
                        "1", file=p1file
                    )  # control for basal friction angle of PHASE 3 map
                else:
                    print("0", file=p1file)
                if ny1:
                    print("1", file=p1file)  # control for viscosity of (PHASE 1) map
                else:
                    print("0", file=p1file)
                if ny2:
                    print("1", file=p1file)  # control for viscosity of PHASE 2 map
                else:
                    print("0", file=p1file)
                if ny3:
                    print("1", file=p1file)  # control for viscosity of PHASE 3 map
                else:
                    print("0", file=p1file)
                if ambdrag:
                    print("1", file=p1file)  # control for ambient drag coefficient map
                else:
                    print("0", file=p1file)
                if flufri:
                    print(
                        "1", file=p1file
                    )  # control for PHASE 2 friction coefficient map
                else:
                    print("0", file=p1file)
                if ctrans12:
                    print(
                        "1", file=p1file
                    )  # control for PHASE 1-PHASE 2 transformation coefficient map
                else:
                    print("0", file=p1file)
                if ctrans13:
                    print(
                        "1", file=p1file
                    )  # control for PHASE 1-PHASE 3 transformation coefficient map
                else:
                    print("0", file=p1file)
                if ctrans23:
                    print(
                        "1", file=p1file
                    )  # control for PHASE 2-PHASE 3 transformation coefficient map
                else:
                    print("0", file=p1file)
                if trelease:
                    print("1", file=p1file)  # control for release time map
                else:
                    print("0", file=p1file)
                if trelstop:
                    print(
                        "1", file=p1file
                    )  # control for release hydrograph stop time map
                else:
                    print("0", file=p1file)
                if stoptime:
                    print("1", file=p1file)  # control for stopping time map
                else:
                    print("0", file=p1file)
                if tslide:
                    print("1", file=p1file)  # control for time of initial sliding map
                else:
                    print("0", file=p1file)
                if hydrograph or hydrocoords:
                    print("1", file=p1file)  # control for hydrograph
                else:
                    print("0", file=p1file)
                if frictiograph:
                    print("1", file=p1file)  # control for frictiograph
                else:
                    print("0", file=p1file)
                if transformograph:
                    print("1", file=p1file)  # control for transformograph
                else:
                    print("0", file=p1file)
                print(elevation, file=p1file)  # name of elevation map
                if model <= 3 and releasemass == 1:
                    print(nhrelease, file=p1file)  # name of release height map
                elif model == 7 and releasemass == 1:
                    print(nhrelease1, file=p1file)  # name of PHASE 1 release height map
                if model == 7 and releasemass == 1:
                    print(nhrelease2, file=p1file)  # name of PHASE 2 release height map
                if model == 7 and releasemass == 1:
                    print(nhrelease3, file=p1file)  # name of PHASE 3 release height map
                if model <= 3 and releasevelocity == 1:
                    print(
                        pf + "_vflowx0000", file=p1file
                    )  # name of release velocity in x direction map
                    print(
                        pf + "_vflowy0000", file=p1file
                    )  # name of release velocity in y direction map
                elif model == 7 and releasevelocity == 1:
                    print(
                        pf + "_vflowx10000", file=p1file
                    )  # name of PHASE 1 release velocity in x direction map
                    print(
                        pf + "_vflowy10000", file=p1file
                    )  # name of PHASE 1 release velocity in y direction map
                if model == 7 and releasevelocity2 == 1:
                    print(
                        pf + "_vflowx20000", file=p1file
                    )  # name of PHASE 2 release velocity in x direction map
                    print(
                        pf + "_vflowy20000", file=p1file
                    )  # name of PHASE 2 release velocity in y direction map
                if model == 7 and releasevelocity3 == 1:
                    print(
                        pf + "_vflowx30000", file=p1file
                    )  # name of PHASE 3 release velocity in x direction map
                    print(
                        pf + "_vflowy30000", file=p1file
                    )  # name of PHASE 3 release velocity in y direction map
                if model <= 3 and entrainment == 1:
                    print(
                        nhentrmax, file=p1file
                    )  # name of maximum height of entrainment map
                elif model == 7 and entrainment == 1:
                    print(
                        nhentrmax1, file=p1file
                    )  # name of maximum height of PHASE 1 entrainment map
                if model == 7 and entrainment == 1:
                    print(
                        nhentrmax2, file=p1file
                    )  # name of maximum height of PHASE 2 entrainment map
                if model == 7 and entrainment == 1:
                    print(
                        nhentrmax3, file=p1file
                    )  # name of maximum height of PHASE 3 entrainment map
                if zones:
                    print(zones, file=p1file)  # name of zones map
                if centr:
                    print(centr, file=p1file)  # name of entrainment coefficient map
                if cvshear:
                    print(
                        cvshear, file=p1file
                    )  # name of shear velocity coefficient map
                if phi1:
                    print(
                        phi1, file=p1file
                    )  # name of internal friction angle of PHASE 1 map
                if phi2:
                    print(
                        phi2, file=p1file
                    )  # name of internal friction angle of PHASE 2 map
                if phi3:
                    print(
                        phi3, file=p1file
                    )  # name of internal friction angle of PHASE 3 map
                if deltab:
                    print(
                        deltab, file=p1file
                    )  # name of basal friction difference map
                if tufri:
                    print(
                        tufri, file=p1file
                    )  # name of turbulent friction coefficient map
                if delta1:
                    print(
                        delta1, file=p1file
                    )  # name of basal friction angle of PHASE 1 map
                if delta2:
                    print(
                        delta2, file=p1file
                    )  # name of basal friction angle of PHASE 2 map
                if delta3:
                    print(
                        delta3, file=p1file
                    )  # name of basal friction angle of PHASE 3 map
                if ny1:
                    print(ny1, file=p1file)  # name of viscosity of (PHASE 1) map
                if ny2:
                    print(ny2, file=p1file)  # name of viscosity of PHASE 2 map
                if ny3:
                    print(ny3, file=p1file)  # name of viscosity of PHASE 3 map
                if ambdrag:
                    print(ambdrag, file=p1file)  # name of ambient drag coefficient map
                if flufri:
                    print(
                        flufri, file=p1file
                    )  # name of PHASE 2 friction coefficient map
                if ctrans12:
                    print(
                        ctrans12, file=p1file
                    )  # name of PHASE 1-PHASE 2 transformation coefficient map
                if ctrans13:
                    print(
                        ctrans13, file=p1file
                    )  # name of PHASE 1-PHASE 3 transformation coefficient map
                if ctrans23:
                    print(
                        ctrans23, file=p1file
                    )  # name of PHASE 2-PHASE 3 transformation coefficient map
                if trelease:
                    print(trelease, file=p1file)  # name of release time map
                if trelstop:
                    print(
                        trelstop, file=p1file
                    )  # name of release hydrograph stop time map
                if stoptime:
                    print(stoptime, file=p1file)  # name of stopping time map
                if tslide:
                    print(tslide, file=p1file)  # name of time of initial sliding map

                if hydrograph or hydrocoords:
                    print(hyd_nin, file=p1file)  # number of input hydrographs
                    print(hyd_nout, file=p1file)  # number of output hydrographs
                    for hydi in range(0, hyd_nin):
                        print(hydrograph[hydi], file=p1file)  # path to hydrograph file
                        print(
                            str(hydlines[hydi]), file=p1file
                        )  # number of hydrograph time steps
                    for hydi in range(0, hyd_nin + hyd_nout):
                        print(
                            hydrocoords[4 * hydi], file=p1file
                        )  # x coordinate of hydrograph
                        print(
                            hydrocoords[4 * hydi + 1], file=p1file
                        )  # y coordinate of hydrograph
                        print(
                            hydrocoords[4 * hydi + 2], file=p1file
                        )  # length of hydrograph profile
                        print(
                            hydrocoords[4 * hydi + 3], file=p1file
                        )  # aspect of hydrograph profile
                else:
                    hyd_nin = 0
                    hyd_nout = 0

                if frictiograph:
                    print(frictiograph, file=p1file)  # path to frictiograph file
                    print(
                        str(frilines), file=p1file
                    )  # number of frictiograph time steps

                if transformograph:
                    print(transformograph, file=p1file)  # path to transformograph file
                    print(
                        str(tralines), file=p1file
                    )  # number of transformograph time steps

                print(lmax, file=p1file)  # number of flow parameters
                for l in range(0, lmax):
                    print(round(gt[l], 10), file=p1file)  # flow parameters

                print(
                    thresholdsc, file=p1file
                )  # threshold of flow height (for computation)
                print(
                    thresholds[0], file=p1file
                )  # threshold of flow height (for display)
                print(thresholds[1], file=p1file)  # threshold of flow kinetic energy
                print(thresholds[2], file=p1file)  # threshold of flow pressure

                print(tint, file=p1file)  # time for writing output
                print(tstop, file=p1file)  # process duration at which to stop
                print(slomo, file=p1file)  # factor for slow motion

                print(slidepar[0], file=p1file)  # search radius for initial sliding
                print(
                    slidepar[1], file=p1file
                )  # exponent for weighting for initial sliding
                print(slidepar[2], file=p1file)  # coefficient for deformation

                print(cfl[0], file=p1file)  # cfl criterion
                print(cfl[1], file=p1file)  # maximum length of time step

                print("%s_" % pf, file=p1file)  # prefix
                print(
                    "%s_results/%s_ascii/" % (pf, pf), file=p1file
                )  # path and prefix for storing output maps
                print(
                    "%s_results/%s_files/" % (pf, pf), file=p1file
                )  # path and prefix for storing output files

                p1file.close()  # closing parameter file

                fparamtstring = jfill
                if model <= 3:
                    fparamtstring = fparamtstring + "\t" + vhrl + "\t" + vhem
                else:
                    fparamtstring = (
                        fparamtstring
                        + "\t"
                        + vhrl
                        + "\t"
                        + rhrls
                        + "\t"
                        + vhem
                        + "\t"
                        + rhems
                    )
                for l in range(0, lmax):
                    fparamtstring = fparamtstring + "\t" + str(gt[l])
                print(
                    fparamtstring, file=fparamt
                )  # writing parameters to temporary input parameter file

                # Creating batch file

                os.mkdir(temppath + "/tmp%s" % jid)  # creating directory for batch file
                strtmp = temppath + "/tmp%s/batch" % jid  # file name for batch file
                out = open(strtmp + str(jid), "w")  # creating batch file
                os.system(
                    "rm -rf " + locpath + "/map%s" % jid
                )  # removing old mapset for model run
                grass.run_command(
                    "g.mapset", flags="c", mapset="map%s" % jid
                )  # creating new mapset for model run
                grass.run_command(
                    "g.mapsets", mapset=ambvars.MAPSET, operation="add"
                )  # making original mapset active

                grass.run_command(
                    "g.mapset", mapset=ambvars.MAPSET
                )  # switching back to original mapset
                os.system(
                    "mkdir " + locpath + "/map%s/.tmp/rtemp" % jid
                )  # creating directory in mapset for model run
                os.environ["PATH"] += (
                    os.pathsep + os.path.join(temppath + "/tmp%s") % jid
                )  # adding path to batch file

                print(
                    """#!/bin/bash
export jid=%s
export SHELL=\"/bin/bash\"
export cellsize=%s
export rnorth=%s
export rsouth=%s
export rwest=%s
export reast=%s
export rtemp=%s
export GRASS_BATCH_JOB=%s/r.avaflow.mult
cp %s/*.txt %s/map%s/.tmp/rtemp/
%s -text %s/map%s
unset GRASS_BATCH_JOB"""
                    % (
                        jid,
                        cellsize,
                        rnorth,
                        rsouth,
                        rwest,
                        reast,
                        temppath,
                        scriptpath2,
                        temppath,
                        locpath,
                        jid,
                        bingrass,
                        locpath,
                        jid,
                    ),
                    file=out,
                )
                # creating batch file

                out.close()  # closing batch file
                fd = os.open(strtmp + str(jid), os.O_RDONLY)  # opening batch file
                os.fchmod(fd, 0o755)  # making batch file executable
                os.close(fd)  # closing batch file

            fparamt.close()  # closing input parameter file

            # Executing batch processing

            start_batch = time.time()  # storing time (start of multi-core processing)

            neff = min(ncores, nruns)
            threadList = list(range(1, neff + 1))
            nameList = list(range(1, nruns + 1))
            threads = []
            for tName in threadList:
                thread = myThread(tName, workQueue)
                thread.start()
                threads.append(thread)
            queueLock.acquire()
            for word in nameList:
                workQueue.put(word)
            queueLock.release()
            while not workQueue.empty():
                pass
            global exitFlag
            exitFlag = 1
            for t in threads:
                t.join()
            print()
            print("Batch processing completed.")
            print()

            for jid in range(1, nruns + 1):
                os.system(
                    "rm -rf " + locpath + "/map" + str(jid)
                )  # removing mapsets for all model runs

            fval = open(
                pf + "_results/" + pf + "_files/" + pf + "_evaluation.txt", "w"
            )  # opening evaluation file for writing

            fvalstring = "nrune\tHmax\tHmaxd\tV1max\tV2max\tV3max\tTmax\tPmax\tHefin\tVE1\tVE2\tVE3"  # adding header to string for evaluation file:
            fvalstring = fvalstring + "\tLi\tLd\tLio\trLi\tLdo\trLd"
            fvalstring = (
                fvalstring
                + "\trTPi\trTNi\trFPi\trFNi\tCSIi\tHSSi\tAUROCi\tD2PCi\tFoCi\tSPIi"
            )
            fvalstring = (
                fvalstring
                + "\trTPd\trTNd\trFPd\trFNd\tCSId\tHSSd\tAUROCd\tD2PCd\tFoCd\tSPId"
            )
            for ictrl in range(0, lctrlpts):
                fvalstring = fvalstring + "\ttreach" + str(ictrl + 1)
                if reftime:
                    fvalstring = (
                        fvalstring
                        + "\ttreachref"
                        + str(ictrl + 1)
                        + "\ttreachrat"
                        + str(ictrl + 1)
                    )

            if impdef == "1":
                grass.mapcalc(
                    '"_obsbini"=if("%s_impactarea">0,1,if("%s_impactarea"==0,0,null()))'
                    % (pf, pf),
                    overwrite=True,
                    quiet=True,
                )  # processed map of observed impact area

            if depdef == "1":
                grass.mapcalc(
                    '"_obsbind"=if("%s_hdeposit">=%s,1,if("%s_hdeposit"==0,0,null()))'
                    % (pf, thresholds[0], pf),
                    overwrite=True,
                    quiet=True,
                )  # processed map of height of observed deposit

            print(fvalstring, file=fval)  # writing header to evaluation file

            # Impact and deposition indicator indices and evaluation file

            fparamt = open(
                temppath + "/fparamt.txt", "r"
            )  # opening temporary input parameter file for reading
            fparami = open(
                pf + "_results/" + pf + "_files/" + pf + "_params.txt", "w"
            )  # opening input parameter file for reading
            print(
                fparamt.readline().replace("\n", ""), file=fparami
            )  # writing header of temporary input parameter file to input parameter file

            for j in range(
                0, 3
            ):  # loop over all output parameters for wich impact indicator index should be derived:
                grass.mapcalc(
                    '"%s_iii%s"=0' % (pf, mstring[j]), overwrite=True, quiet=True
                )  # initializing impact indicator index map
                grass.mapcalc(
                    '"%s_dii"=0' % pf, overwrite=True, quiet=True
                )  # initializing deposit indicator index map

            nsuccess = 0  # intitalizing counter for number of successful simulations
            for jid in nameList:  # loop over all model runs:

                ftimesteps = open(
                    pf
                    + "_results/"
                    + pf
                    + "_files/"
                    + pf
                    + "_nout"
                    + str(jid)
                    + ".txt",
                    "r",
                )
                # opening file with number of time steps and success
                ftimesteps.readline()  # number of time steps is not needed
                csuccess = (
                    ftimesteps.readline()
                )  # reading control for success of simulation
                csuccess = int(csuccess.replace("\n", ""))  # removing newline
                maxvflow = ftimesteps.readline()  # reading maximum flow velocity
                maxvflow = float(maxvflow.replace("\n", ""))  # removing newline
                vol_entr = ftimesteps.readline()  # reading entrained (PHASE 1) volume
                vol_entr = float(vol_entr.replace("\n", ""))  # removing newline
                vol_entr2 = ftimesteps.readline()  # reading entrained PHASE 2 volume
                vol_entr2 = float(vol_entr2.replace("\n", ""))  # removing newline
                vol_entr3 = ftimesteps.readline()  # reading entrained PHASE 3 volume
                vol_entr3 = float(vol_entr3.replace("\n", ""))  # removing newline
                ftimesteps.close()  # closing file with number of time steps

                if jid < 10:
                    jfill = "00000" + str(jid)  # formatting model run string
                elif jid < 100:
                    jfill = "0000" + str(jid)
                elif jid < 1000:
                    jfill = "000" + str(jid)
                elif jid < 10000:
                    jfill = "00" + str(jid)
                elif jid < 100000:
                    jfill = "0" + str(jid)
                else:
                    jfill = str(jid)

                if csuccess == 1:
                    fvalstring = jfill  # writing number of model run to string

                if model <= 3:
                    mstringlist = [
                        "_hflow_max",
                        "_hflow_fin",
                        "_vflow_max",
                        "none",
                        "none",
                        "_tflow_max",
                        "_pflow_max",
                        "_basechange_fin",
                        "_treach",
                    ]
                elif model == 7:
                    mstringlist = [
                        "_hflow_max",
                        "_hflow_fin",
                        "_vflow1_max",
                        "_vflow2_max",
                        "_vflow3_max",
                        "_tflow_max",
                        "_pflow_max",
                        "_basechange_fin",
                        "_treach",
                    ]

                if csuccess == 1:

                    for mstringi in mstringlist:

                        if not mstringi == "none":

                            grass.run_command(
                                "r.in.gdal",
                                input=ascpath + pf + mstringi + str(jid) + ".asc",
                                output=pf + mstringi + str(jid),
                                overwrite=True,
                            )  # importing map

                    minval = grass.raster_info(pf + "_basechange_fin" + str(jid))[
                        "min"
                    ]  # minimum value of basal change
                    maxval = grass.raster_info(pf + "_basechange_fin" + str(jid))[
                        "max"
                    ]  # maximum value of basal change

                    if minval == 0 and maxval == 0:
                        basechange = 0
                    else:
                        basechange = 1

                    if basechange == 1:
                        grass.mapcalc(
                            '"%s_hflow_dep%s"=if("%s_basechange_fin%s">0,"%s_basechange_fin%s",0)'
                            % (pf, str(jid), pf, str(jid), pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                    else:
                        grass.mapcalc(
                            '"%s_hflow_dep%s"="%s_hflow_fin%s"'
                            % (pf, str(jid), pf, str(jid)),
                            overwrite=True,
                            quiet=True,
                        )
                    # simulated height of deposition map

                    if model <= 3:
                        mstringlist = [
                            "_hflow_max",
                            "_hflow_dep",
                            "_vflow_max",
                            "none",
                            "none",
                            "_tflow_max",
                            "_pflow_max",
                            "_basechange_fin",
                        ]
                    elif model == 7:
                        mstringlist = [
                            "_hflow_max",
                            "_hflow_dep",
                            "_vflow1_max",
                            "_vflow2_max",
                            "_vflow3_max",
                            "_tflow_max",
                            "_pflow_max",
                            "_basechange_fin",
                        ]

                    for mstringi in mstringlist:

                        if not mstringi == "none":

                            if mstringi == "_basechange_fin":
                                maxval = (
                                    grass.raster_info(pf + mstringi + str(jid))["min"]
                                    * -1
                                )
                            else:
                                maxval = grass.raster_info(pf + mstringi + str(jid))[
                                    "max"
                                ]  # maximum value

                            maxval = grass.raster_info(pf + mstringi + str(jid))[
                                "max"
                            ]  # maximum value
                            fvalstring = (
                                fvalstring + "\t" + str(maxval)
                            )  # adding maximum value to string
                            if mstringi == "_tflow_max" and maxval == 0:
                                csuccess = (
                                    0
                                )  # setting control for success to negative, if necessary

                        else:
                            fvalstring = (
                                fvalstring + "\t-9999"
                            )  # adding maximum value to string

                if csuccess == 1:

                    fvalstring = (
                        fvalstring + "\t" + str(vol_entr)
                    )  # adding entrained/deposited (PHASE 1) volume to string
                    if model == 7:
                        fvalstring = (
                            fvalstring + "\t" + str(vol_entr2)
                        )  # adding entrained/deposited PHASE 2 volume to string
                    else:
                        fvalstring = fvalstring + "\t-9999"
                    if model == 7:
                        fvalstring = (
                            fvalstring + "\t" + str(vol_entr3)
                        )  # adding entrained/deposited PHASE 3 volume to string
                    else:
                        fvalstring = fvalstring + "\t-9999"

                flinet = (
                    fparamt.readline()
                )  # reading line of temporary input parameter file

                if csuccess == 1:  # if simulation was successful:

                    nsuccess = nsuccess + 1  # updating number of successful simulations

                    grass.mapcalc(
                        '"%s_iis%s"=if("%s_hflow_max%s">%s,1,0)'
                        % (pf, jfill, pf, str(jid), thresholds[0]),
                        overwrite=True,
                        quiet=True,
                    )  # impact indicator score map

                    grass.mapcalc(
                        '"%s_dis%s"=if("%s_hflow_dep%s">%s,1,0)'
                        % (pf, jfill, pf, str(jid), thresholds[0]),
                        overwrite=True,
                        quiet=True,
                    )  # deposit indicator score map

                    for j in range(
                        0, 3
                    ):  # loop over all output parameters for wich impact indicator index should be derived:

                        grass.mapcalc(
                            '"%s_iii%s"=if("%s%s_max%s">%s,"%s_iii%s"+1,"%s_iii%s")'
                            % (
                                pf,
                                mstring[j],
                                pf,
                                mstring[j],
                                str(jid),
                                thresholds[j],
                                pf,
                                mstring[j],
                                pf,
                                mstring[j],
                            ),
                            overwrite=True,
                            quiet=True,
                        )  # updating impact indicator index map

                    grass.mapcalc(
                        '"%s_dii"=if("%s_hflow_dep%s">%s,"%s_dii"+1,"%s_dii")'
                        % (pf, pf, str(jid), thresholds[0], pf, pf),
                        overwrite=True,
                        quiet=True,
                    )  # updating deposit indicator index map

                    if profile:  # if profile is defined:

                        rstring = "_hflow0000"
                        grass.run_command(
                            "r.profile",
                            input=pf + rstring + str(jid),
                            output=temppath + "/fhrelease.txt",
                            coordinates=profile,
                            null="-9999",
                            quiet=True,
                            overwrite=True,
                        )
                        # profile of release height

                        fprof = open(
                            temppath + "/fhrelease.txt", "r"
                        )  # opening profile of release height for reading
                        flines = sum(1 for line in fprof)  # counting number of lines
                        fprof.close()  # closing file

                        flowline_start = (
                            0.0
                        )  # setting initial value for start to start point of profile
                        if releasemass == 1:  # if release height is defined:

                            fprof = open(
                                temppath + "/fhrelease.txt"
                            )  # again, opening profile file for reading
                            for j in range(1, flines):  # loop over all lines:
                                fline = fprof.readline()  # reading and splitting line
                                fline = list(map(str, fline.split(" ")))
                                if float(fline[2].replace("\n", "")) > 0:
                                    flowline_start = float(
                                        fline[1]
                                    )  # identifying profile position at start of flow
                                    break
                            fprof.close()  # closing file

                        if (
                            impdef == "1"
                        ):  # identifying profile position at terminus of observed flow:

                            mstringj = 0  # initializing counter
                            grass.run_command(
                                "r.profile",
                                input=pf + "_impactarea",
                                output=temppath + "/fimp.txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of map

                            flowline_stopimp = None
                            fprof = open(
                                temppath + "/fimp.txt"
                            )  # opening profile for reading
                            for j in range(1, flines):  # loop over all lines:
                                fline = fprof.readline()  # reading and splitting line
                                fline = list(map(str, fline.split(" ")))
                                if float(fline[2].replace("\n", "")) == 1 or (
                                    not flowline_stopimp and j == flines - 1
                                ):
                                    flowline_stopimp = float(
                                        fline[1]
                                    )  # identifying profile position at terminus of observed flow
                            fprof.close()  # closing file

                        if (
                            depdef == "1"
                        ):  # identifying profile position at terminus of observed flow:

                            mstringj = 0  # initializing counter
                            grass.run_command(
                                "r.profile",
                                input=pf + "_hdeposit",
                                output=temppath + "/fobs.txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of map

                            flowline_stopdep = None
                            fprof = open(
                                temppath + "/fobs.txt"
                            )  # opening profile for reading
                            for j in range(1, flines):  # loop over all lines:
                                fline = fprof.readline()  # reading and splitting line
                                fline = list(map(str, fline.split(" ")))
                                if float(fline[2].replace("\n", "")) >= float(
                                    thresholds[0]
                                ) or (not flowline_stopdep and j == flines - 1):
                                    flowline_stopdep = float(
                                        fline[1]
                                    )  # identifying profile position at terminus of observed flow
                            fprof.close()  # closing file

                        mstringj = 0  # initializing counter
                        for mstringi in ["_hflow_max", "_hflow_dep"]:

                            grass.run_command(
                                "r.profile",
                                input=pf + mstringi + str(jid),
                                output=temppath + "/f" + mstringi + ".txt",
                                coordinates=profile,
                                null="-9999",
                                quiet=True,
                                overwrite=True,
                            )  # profile of map

                            flowline_stop = None
                            fprof = open(
                                temppath + "/f" + mstringi + ".txt"
                            )  # opening profile for reading
                            for j in range(1, flines):  # loop over all lines:
                                fline = fprof.readline()  # reading and splitting line
                                fline = list(map(str, fline.split(" ")))
                                if float(fline[2].replace("\n", "")) >= float(
                                    thresholds[mstringj]
                                ) or (not flowline_stop and j == flines - 1):
                                    flowline_stop = float(
                                        fline[1]
                                    )  # identifying profile position at termination of flow
                            fprof.close()  # closing file
                            maxlength = flowline_stop - flowline_start  # flow length
                            fvalstring = (
                                fvalstring + "\t" + str(maxlength)
                            )  # adding flow length to string
                            if impdef == "1" and mstringi == "_hflow_max":
                                maxlengthimp = (
                                    flowline_stopimp - flowline_start
                                )  # observed movement length
                                flowline_diffimp = flowline_stop - flowline_stopimp
                                if maxlengthimp > 0:
                                    lratio_imp = flowline_diffimp / maxlengthimp
                                else:
                                    lratio_imp = -9999
                            elif depdef == "1" and mstringi == "_hflow_dep":
                                maxlengthdep = (
                                    flowline_stopdep - flowline_start
                                )  # observed movement length
                                flowline_diffdep = flowline_stop - flowline_stopdep
                                if maxlengthdep > 0:
                                    lratio_dep = flowline_diffdep / maxlengthdep
                                else:
                                    lratio_dep = -9999
                            if not mstringi == "_hflow_max":
                                mstringj = mstringj + 1  # updating counter

                        if impdef == "1":
                            fvalstring = (
                                fvalstring
                                + "\t"
                                + str(maxlengthimp)
                                + "\t"
                                + str(lratio_imp)
                            )
                        else:
                            fvalstring = fvalstring + "\t-9999\t-9999"
                        if depdef == "1":
                            fvalstring = (
                                fvalstring
                                + "\t"
                                + str(maxlengthdep)
                                + "\t"
                                + str(lratio_dep)
                            )
                        else:
                            fvalstring = fvalstring + "\t-9999\t-9999"

                    else:
                        fvalstring = (
                            fvalstring + "\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999"
                        )
                        lratio_imp = -9999
                        lratio_dep = -9999

                    if impdef == "1" and depdef == "1":
                        mstringlist2 = ["_hflow_max", "_hflow_dep"]
                    # raster maps to be processed for evaluation scores
                    elif impdef == "1":
                        mstringlist2 = ["_hflow_max"]
                    elif depdef == "1":
                        mstringlist2 = ["_hflow_dep"]
                    else:
                        mstringlist2 = []

                    if not impdef == "1":
                        fvalstring = (
                            fvalstring
                            + "\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999"
                        )

                    mstringj = 0
                    for mstringi in mstringlist2:  # loop over all maps to be processed:

                        grass.mapcalc(
                            '"_simbin"=if("%s%s%s">=%s,1,0)'
                            % (pf, mstringi, str(jid), thresholds[mstringj]),
                            overwrite=True,
                            quiet=True,
                        )  # processed map of simulated impact area or deposit

                        if mstringi == "_hflow_max":
                            grass.mapcalc(
                                '"_eval"="_simbin"+10*"_obsbini"',
                                overwrite=True,
                                quiet=True,
                            )
                        # combination of observed and simulated impact areas
                        elif mstringi == "_hflow_dep":
                            grass.mapcalc(
                                '"_eval"="_simbin"+10*"_obsbind"',
                                overwrite=True,
                                quiet=True,
                            )
                        # combination of observed and simulated deposits

                        grass.run_command(
                            "r.stats",
                            flags="c" + "n",
                            separator="x",
                            input="_eval",
                            output="%s/evaluation.txt" % temppath,
                            quiet=True,
                            overwrite=True,
                        )  # scores tp=11, tn=0, fp=1, fn=10

                        valfile = open(
                            temppath + "/evaluation.txt", "r"
                        )  # opening raw file with evaluation scores
                        valfile_lines = sum(
                            1 for line in valfile
                        )  # reading number of lines
                        valfile.close()  # closing file

                        if (
                            valfile_lines < 4
                        ):  # if number of lines is smaller than 4 (at least one of the scores is 0):
                            valfile = open(
                                temppath + "/evaluation.txt", "r"
                            )  # opening file again
                            lr = valfile.readline()  # reading first line
                            l1 = lr.find("0x")  # checking for tn score
                            if l1 > -1:
                                lr = valfile.readline()  # if found, reading next line
                            l2 = lr.find("1x")  # checking for fp score
                            if l2 > -1:
                                lr = valfile.readline()  # if found, reading next line
                            l3 = lr.find(
                                "10x"
                            )  # checking for fn score (if not found, -1 is returned)
                            if l3 > -1:
                                lr = valfile.readline()  # if found, reading next line
                            l4 = lr.find("11x")  # checking for tp score
                            valfile.close()  # closing file again
                        else:
                            l1 = (
                                1
                            )  # if number of lines is 4, setting all controls to positive (any value larger than -1 would be fine)
                            l2 = 1
                            l3 = 1
                            l4 = 1

                        valfile = open(
                            temppath + "/evaluation.txt", "r"
                        )  # once more opening temporary file with evaluation scores
                        if l1 == -1:
                            tn = float(
                                0
                            )  # if no tn score was identified before, setting to 0
                        else:
                            tn = valfile.readline()  # else, reading line
                            tn = float(
                                tn.replace("0x", "").replace("\n", "")
                            )  # extracting tn score
                        if l2 == -1:
                            fp = float(
                                0
                            )  # if no fp score was identified before, setting to 0
                        else:
                            fp = valfile.readline()  # else, reading line
                            fp = float(
                                fp.replace("1x", "").replace("\n", "")
                            )  # extracting fp score
                        if l3 == -1:
                            fn = float(
                                0
                            )  # if no fn score was identified before, setting to 0
                        else:
                            fn = valfile.readline()  # else, reading line
                            fn = float(
                                fn.replace("10x", "").replace("\n", "")
                            )  # extracting fn score
                        if l4 == -1:
                            tp = float(
                                0
                            )  # if no tp score was identified before, setting to 0
                        else:
                            tp = valfile.readline()  # else, reading line
                            tp = float(
                                tp.replace("11x", "").replace("\n", "")
                            )  # extracting tp score
                        valfile.close()  # closing file for the last time

                        tn = 5 * (tp + fn) - fp  # normalizing tn score
                        if tp + fn > 0:
                            foc = (tp + fp) / (tp + fn)  # conservativeness factor
                        else:
                            foc = -9999

                        if tn + fp > 0 and tp + fn > 0:

                            total = tp + tn + fp + fn  # sum of scores

                            tpr = 100 * tp / total  # tp rate in percent
                            tnr = 100 * tn / total  # tn rate in percent
                            fpr = 100 * fp / total  # fp rate in percent
                            fnr = 100 * fn / total  # fn rate in percent
                            opr = tpr + fnr  # rate of positive observations
                            onr = tnr + fpr  # rate of negative observations

                            if not tpr + fpr + fnr == 0:
                                csi = tpr / (tpr + fpr + fnr)  # critical success index
                            else:
                                csi = -9999
                            if (
                                not (tpr + fnr) * (fnr + tnr)
                                + (tpr + fpr) * (fpr + tnr)
                                == 0
                            ):
                                hss = (2 * (tpr * tnr) - (fpr * fnr)) / (
                                    (tpr + fnr) * (fnr + tnr)
                                    + (tpr + fpr) * (fpr + tnr)
                                )  # Heidke skill score
                            else:
                                hss = -9999

                            tprr = tpr / opr  # rate of true positive rate
                            fprr = fpr / onr  # rate of false positive rate
                            aucroc = (
                                (1 - fprr) * tprr
                                + fprr * tprr / 2
                                + (1 - fprr) * (1 - tprr) / 2
                            )  # area under the roc curve
                            d2pc = pow(
                                pow(1 - tprr, 2) + pow(fprr, 2), 0.5
                            )  # distance to perfect classification

                        else:
                            tpr = -9999
                            tnr = -9999
                            fpr = -9999
                            fnr = -9999
                            csi = -9999
                            hss = -9999
                            aucroc = -9999
                            d2pc = -9999

                        if profile:
                            if mstringi == "_hflow_max":
                                lratio = lratio_imp
                            else:
                                lratio = lratio_dep
                        else:
                            lratio = -1
                        if lratio == -1:
                            lratio = -1.1
                        if foc == 0:
                            foc = -1

                        spi = (
                            max(0, min(0.2, 0.2 * csi))
                            + max(0, min(0.2, 0.2 * hss))
                            + max(0, min(0.2, 0.2 * (1 - d2pc)))
                            + max(
                                0,
                                min(
                                    0.2,
                                    0.2
                                    * pow(
                                        math.sin(
                                            math.pi
                                            / 2
                                            * min(lratio + 1, 1 / (lratio + 1))
                                        ),
                                        5,
                                    ),
                                ),
                            )
                            + max(
                                0,
                                min(
                                    0.2,
                                    0.2
                                    * pow(math.sin(math.pi / 2 * min(foc, 1 / foc)), 5),
                                ),
                            )
                        )  # synthetic performance index

                        if lratio == -1.1:
                            lratio = -9999
                        if foc == -1:
                            foc = -9999

                        fvalstring = (
                            fvalstring
                            + "\t"
                            + str(tpr)
                            + "\t"
                            + str(tnr)
                            + "\t"
                            + str(fpr)
                            + "\t"
                            + str(fnr)
                        )
                        fvalstring = (
                            fvalstring
                            + "\t"
                            + str(csi)
                            + "\t"
                            + str(hss)
                            + "\t"
                            + str(aucroc)
                            + "\t"
                            + str(d2pc)
                            + "\t"
                            + str(foc)
                            + "\t"
                            + str(spi)
                        )
                        # adding evaluation results to string

                        if not mstringi == "_hflow_max":
                            mstringj = mstringj + 1  # updating counter

                    if not depdef == "1":
                        fvalstring = (
                            fvalstring
                            + "\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999\t-9999"
                        )

                    for ictrl in range(0, lctrlpts):  # loop over all control points:

                        ctrlstring = grass.read_command(
                            "r.what",
                            map=pf + "_treach" + str(jid),
                            separator="\t",
                            coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                            overwrite=True,
                            quiet=True,
                        )  # reading impact indicator index at control point
                        ctrlstring = (
                            ctrlstring.replace("\t\t", "\t")
                            .replace("\n", "")
                            .replace(str(xctrl[ictrl]) + "\t", "")
                            .replace(str(yctrl[ictrl]) + "\t", "")
                        )
                        if ctrlstring == "*":
                            ctrlstring = "9999"
                        fvalstring = fvalstring + "\t" + ctrlstring
                        if reftime and not ctrlstring == "9999":
                            fvalstring = fvalstring + "\t" + reftime[ictrl]
                            fvalstring = (
                                fvalstring
                                + "\t"
                                + str(float(ctrlstring) / float(reftime[ictrl]))
                            )
                        elif reftime:
                            fvalstring = fvalstring + "\t" + reftime[ictrl]
                            fvalstring = fvalstring + "\t" + "-9999"

                    print(
                        fvalstring, file=fval
                    )  # writing string for model run to text file

                    # Input for AIMEC

                    if jid < 10:
                        jfill = "00000" + str(jid)  # formatting model run string
                    elif jid < 100:
                        jfill = "0000" + str(jid)
                    elif jid < 1000:
                        jfill = "000" + str(jid)
                    elif jid < 10000:
                        jfill = "00" + str(jid)
                    elif jid < 100000:
                        jfill = "0" + str(jid)
                    else:
                        jfill = str(jid)

                    flinep = flinet.replace("\t", ";").replace(
                        "\n", ""
                    )  # preparing line for writing to parameter input file  for AIMEC
                    print(
                        flinep, file=faimecp
                    )  # writing line of input parameter file to parameter input file for AIMEC

                    grass.mapcalc(
                        '"_hflowii"=if("%s_hflow_max%s">=%s,"%s_hflow_max%s",0)'
                        % (pf, str(jid), thresholds[0], pf, str(jid)),
                        overwrite=True,
                        quiet=True,
                    )  # constraining maximum flow height map
                    grass.mapcalc(
                        '"_tflowii"=if("%s_tflow_max%s">=%s,"%s_tflow_max%s",0)'
                        % (pf, str(jid), thresholds[1], pf, str(jid)),
                        overwrite=True,
                        quiet=True,
                    )  # constraining maximum flow kinetic energy map

                    grass.run_command(
                        "r.out.gdal",
                        input="_hflowii",
                        output=pf + "_results/" + pf + "_aimec/depth/" + jfill + ".asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting maximum flow height raster to ascii
                    grass.run_command(
                        "r.out.gdal",
                        input="_tflowii",
                        output=pf
                        + "_results/"
                        + pf
                        + "_aimec/pressure/"
                        + jfill
                        + ".asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting maximum flow kinetic energy raster to ascii

                    corrasc(pf + "_results/" + pf + "_aimec/depth/" + jfill)
                    corrasc(pf + "_results/" + pf + "_aimec/pressure/" + jfill)

                    flinei = flinet.replace(
                        "\n", ""
                    )  # preparing line for writing to input parameter file

                else:

                    grass.mapcalc(
                        '"%s_iis%s"=0' % (pf, jfill), overwrite=True, quiet=True
                    )  # impact indicator score map
                    grass.mapcalc(
                        '"%s_dis%s"=0' % (pf, jfill), overwrite=True, quiet=True
                    )  # deposit indicator score map
                    flinei = "*" + flinet.replace(
                        "\n", ""
                    )  # preparing line for writing to input parameter file

                print(
                    flinei, file=fparami
                )  # writing line of temporary input parameter file to input parameter file

            fval.close()  # closing evaluation file
            fparamt.close()  # closing temporary input parameter file
            fparami.close()  # closing input parameter file
            faimecp.close()  # closing parameter input file for AIMEC

            # Finalizing impact indicator index maps

            for j in range(0, 3):

                grass.mapcalc(
                    '"%s_iii%s"=float("%s_iii%s")/float(%s)'
                    % (pf, mstring[j], pf, mstring[j], str(nsuccess)),
                    overwrite=True,
                    quiet=True,
                )

                grass.run_command(
                    "r.out.gdal",
                    input=pf + "_iii" + mstring[j],
                    output=ascpath + pf + "_iii" + mstring[j] + ".asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # exporting impact indicator index map to ascii
                corrasc(ascpath + pf + "_iii" + mstring[j])

            # Finalizing deposit indicator index maps

            grass.mapcalc(
                '"%s_dii"=float("%s_dii")/float(%s)' % (pf, pf, str(nsuccess)),
                overwrite=True,
                quiet=True,
            )

            grass.run_command(
                "r.out.gdal",
                input=pf + "_dii",
                output=ascpath + pf + "_dii.asc",
                format="AAIGrid",
                overwrite=True,
            )  # exporting deposit indicator index map to ascii
            corrasc(ascpath + pf + "_dii")

            # Finalizing multi-core processing

            stop_batch = time.time()  # storing time (stop of multi-core processing)
            comptime_batch = (
                stop_batch - start_batch
            )  # storing computational time for batch processing in seconds

            timefile = open(filepath + pf + "_time.txt", "w")
            timefile.write(
                str(comptime_batch)
            )  # writing computational time for batch processing to file
            timefile.close()

            # Creating control point information file

            if ctrlpoints:  # if control points are defined:

                fctrlpts = open(
                    pf + "_results/" + pf + "_files/" + pf + "_ctrlpoints.txt", "w"
                )  # opening control points information file for writing
                fctrlpts.write("Id\tx\ty\tiii\n")  # writing header to file

                for ictrl in range(0, lctrlpts):  # loop over all control points:

                    ctrlstring = grass.read_command(
                        "r.what",
                        map=pf + "_iii_hflow",
                        separator="\t",
                        coordinates=xctrl[ictrl] + "," + yctrl[ictrl],
                        overwrite=True,
                        quiet=True,
                    )  # reading impact indicator index at control point
                    ctrlstring = ctrlstring.replace("\t\t", "\t").replace("\n", "")

                    fctrlpts.write("C" + str(ictrl + 1) + "\t" + ctrlstring + "\n")
                    # writing point number and coordinates to file

                fctrlpts.close()  # closing file

        # Writing documentation file for flags and options

        documentation = open(
            pf + "_results/" + pf + "_files/" + pf + "_documentation.txt", "w"
        )  # opening documentation file for writing

        if aflag:
            print("Map plots of pressure and kinetic energy\tTRUE", file=documentation)
        else:
            print("Map plots of pressure and kinetic energy\tFALSE", file=documentation)
        if eflag:
            print("Model execution\tTRUE", file=documentation)
        else:
            print("Model execution\tFALSE", file=documentation)
        if kflag:
            print("Keeping GRASS rasters of results\tTRUE", file=documentation)
        else:
            print("Keeping GRASS rasters of results\tFALSE", file=documentation)
        if mflag:
            print("Multiple model runs\tTRUE", file=documentation)
        else:
            print("Multiple model runs\tFALSE", file=documentation)
        if tflag:
            print(
                "Map plots of impact wave or tsunami height\tTRUE", file=documentation
            )
        else:
            print(
                "Map plots of impact wave or tsunami height\tFALSE", file=documentation
            )
        if vflag:
            print("Evaluation and visualization\tTRUE", file=documentation)
        else:
            print("Evaluation and visualization\tFALSE", file=documentation)

        print("Prefix\t%s" % pf, file=documentation)

        print("Cell size\t%s" % cellsize, file=documentation)
        print("Northern boundary\t%s" % rnorth, file=documentation)
        print("Southern boundary\t%s" % rsouth, file=documentation)
        print("Western boundary\t%s" % rwest, file=documentation)
        print("Eastern boundary\t%s" % reast, file=documentation)
        print("Model\t%s" % str(model), file=documentation)

        print("Number of phases\t%s" % str(len(phases)), file=documentation)
        for i in range(0, len(phases)):
            print("Phase %i\t%s" % (i + 1, phasetext[i]), file=documentation)

        print("Gravity\t%s" % gravity, file=documentation)
        print("Numerical limiter\t%s" % limiter, file=documentation)
        print("Flow height conversion\t%s" % controls[0], file=documentation)
        print("Diffusion control\t%s" % controls[1], file=documentation)
        print("Surface control\t%s" % controls[2], file=documentation)
        print("Entrainment\t%s" % controls[3], file=documentation)
        print("Deposition and stopping\t%s" % controls[4], file=documentation)
        print("Dynamic adaptation of friction\t%s" % controls[5], file=documentation)

        if cores:
            print("Number of cores\t%s" % cores, file=documentation)
        else:
            print("Number of cores\tFALSE", file=documentation)
        if mflag:
            print("Sampling strategy\t%s" % sampling, file=documentation)
        else:
            print("Sampling strategy\tFALSE", file=documentation)

        if orthophoto:
            print("Orthophoto\t%s" % orthophoto, file=documentation)
        else:
            print("Orthophoto\tFALSE", file=documentation)
        if impactarea:
            print("Observed impact area\t%s" % impactarea, file=documentation)
        else:
            print("Observed impact area\tFALSE", file=documentation)
        if hdeposit:
            print("Height of observed deposit\t%s" % hdeposit, file=documentation)
        else:
            print("Height of observed deposit\tFALSE", file=documentation)
        print("Elevation\t%s" % elevation, file=documentation)

        if hrelease:
            print("Release height\t%s" % hrelease, file=documentation)
        else:
            print("Release height\tFALSE", file=documentation)
        if model == 7:
            if hrelease1:
                print("Solid release height\t%s" % hrelease1, file=documentation)
            else:
                print("Solid release height\tFALSE", file=documentation)
            if hrelease2:
                print("PHASE 2 release height\t%s" % hrelease2, file=documentation)
            else:
                print("PHASE 2 release height\tFALSE", file=documentation)
            if rhrelease1:
                print(
                    "Fraction of PHASE 1s within release height\t%s" % rhrelease1,
                    file=documentation,
                )
            else:
                print(
                    "Fraction of PHASE 1s within release height\tFALSE",
                    file=documentation,
                )
            if hrelease3:
                print("PHASE 3 release height\t%s" % hrelease3, file=documentation)
            else:
                print("PHASE 3 release height\tFALSE", file=documentation)
        if mflag and vhrelease:
            print("Variation of release height\t%s" % vhrelease, file=documentation)
        else:
            print("Variation of release height\tFALSE", file=documentation)
        print("Use of release mass\t%s" % str(releasemass), file=documentation)

        if hentrmax:
            print("Maximum entrainable height\t%s" % hentrmax, file=documentation)
        else:
            print("Maximum entrainable height\tFALSE", file=documentation)
        if model == 7:
            if hentrmax1:
                print(
                    "Maximum entrainable PHASE 1 height\t%s" % hentrmax1,
                    file=documentation,
                )
            else:
                print("Maximum entrainable PHASE 1 height\tFALSE", file=documentation)
            if hentrmax2:
                print(
                    "Maximum entrainable PHASE 2 height\t%s" % hentrmax1,
                    file=documentation,
                )
            else:
                print("Maximum entrainable PHASE 2 height\tFALSE", file=documentation)
            if rhentrmax1:
                print(
                    "Fraction of PHASE 1s within entrainable height\t%s" % rhentrmax1,
                    file=documentation,
                )
            else:
                print(
                    "Fraction of PHASE 1s within entrainable height\tFALSE",
                    file=documentation,
                )
            if hentrmax3:
                print(
                    "Maximum entrainable PHASE 3 height\t%s" % hentrmax3,
                    file=documentation,
                )
            else:
                print("Maximum entrainable PHASE 3 height\tFALSE", file=documentation)
        if mflag and vhentrmax:
            print("Variation of entrainable height\t%s" % vhentrmax, file=documentation)
        else:
            print("Variation of entrainable height\tFALSE", file=documentation)

        if model <= 3:
            if releasevelocity == 1:
                print("Release velocity in x direction\t%s" % vinx, file=documentation)
                print("Release velocity in y direction\t%s" % viny, file=documentation)
            else:
                print("Release velocity in x direction\tFALSE", file=documentation)
                print("Release velocity in y direction\tFALSE", file=documentation)
        elif model == 7:
            if releasevelocity == 1:
                print(
                    "PHASE 1 release velocity in x direction\t%s" % vinx1,
                    file=documentation,
                )
                print(
                    "PHASE 1 release velocity in y direction\t%s" % viny1,
                    file=documentation,
                )
            else:
                print(
                    "PHASE 1 release velocity in x direction\tFALSE", file=documentation
                )
                print(
                    "PHASE 1 release velocity in y direction\tFALSE", file=documentation
                )
            if releasevelocity2 == 1:
                print(
                    "PHASE 2 release velocity in x direction\t%s" % vinx2,
                    file=documentation,
                )
                print(
                    "PHASE 2 release velocity in y direction\t%s" % viny3,
                    file=documentation,
                )
            else:
                print(
                    "PHASE 2 release velocity in x direction\tFALSE", file=documentation
                )
                print(
                    "PHASE 2 release velocity in y direction\tFALSE", file=documentation
                )
            if releasevelocity3 == 1:
                print(
                    "PHASE 3 release velocity in x direction\t%s" % vinx3,
                    file=documentation,
                )
                print(
                    "PHASE 3 release velocity in y direction\t%s" % viny3,
                    file=documentation,
                )
            else:
                print(
                    "PHASE 3 release velocity in x direction\tFALSE", file=documentation
                )
                print(
                    "PHASE 3 release velocity in y direction\tFALSE", file=documentation
                )

        if phi1:
            print("PHASE 1 internal friction angle map\t%s" % phi1, file=documentation)
        else:
            print("PHASE 1 internal friction angle map\tFALSE", file=documentation)
        if phi2:
            print("PHASE 2 internal friction angle map\t%s" % phi2, file=documentation)
        else:
            print("PHASE 2 internal friction angle map\tFALSE", file=documentation)
        if phi3:
            print("PHASE 3 internal friction angle map\t%s" % phi3, file=documentation)
        else:
            print("PHASE 3 internal friction angle map\tFALSE", file=documentation)
        if deltab:
            print(
                "Basal friction difference map\t%s" % deltab,
                file=documentation,
            )
        else:
            print(
                "Basal friction difference map\tFALSE",
                file=documentation,
            )
        if tufri:
            print("Turbulent friction coefficient map\t%s" % tufri, file=documentation)
        else:
            print("Turbulent friction coefficient map\tFALSE", file=documentation)
        if delta1:
            print("PHASE 1 basal friction angle map\t%s" % delta1, file=documentation)
        else:
            print("PHASE 1 basal friction angle map\tFALSE", file=documentation)
        if delta2:
            print("PHASE 2 basal friction angle map\t%s" % delta2, file=documentation)
        else:
            print("PHASE 2 basal friction angle map\tFALSE", file=documentation)
        if delta3:
            print("PHASE 3 basal friction angle map\t%s" % delta3, file=documentation)
        else:
            print("PHASE 3 basal friction angle map\tFALSE", file=documentation)
        if ny1:
            print("PHASE 1 viscosity map\t%s" % ny1, file=documentation)
        else:
            print("PHASE 1 viscosity map\tFALSE", file=documentation)
        if ny2:
            print("PHASE 2 viscosity map\t%s" % ny2, file=documentation)
        else:
            print("PHASE 2 viscosity map\tFALSE", file=documentation)
        if ny3:
            print("PHASE 3 viscosity map\t%s" % ny3, file=documentation)
        else:
            print("PHASE 3 viscosity map\tFALSE", file=documentation)
        if ambdrag:
            print("Ambient drag coefficient map\t%s" % ambdrag, file=documentation)
        else:
            print("Ambient drag coefficient map\tFALSE", file=documentation)
        if flufri:
            print("Fluid friction coefficient map\t%s" % flufri, file=documentation)
        else:
            print("Fluid friction coefficient map\tFALSE", file=documentation)
        if ctrans12:
            print(
                "PHASE 1-PHASE 2 transformation coefficient map\t%s" % ctrans12,
                file=documentation,
            )
        else:
            print(
                "PHASE 1-PHASE 2 transformation coefficient map\tFALSE",
                file=documentation,
            )
        if ctrans13:
            print(
                "PHASE 1-PHASE 3 transformation coefficient map\t%s" % ctrans13,
                file=documentation,
            )
        else:
            print(
                "PHASE 1-PHASE 3 transformation coefficient map\tFALSE",
                file=documentation,
            )
        if ctrans23:
            print(
                "PHASE 2-PHASE 3 transformation coefficient map\t%s" % ctrans23,
                file=documentation,
            )
        else:
            print(
                "PHASE 2-PHASE 3 transformation coefficient map\tFALSE",
                file=documentation,
            )
        if zones:
            print("Zones map\t%s" % zones, file=documentation)
        else:
            print("Zones map\tFALSE", file=documentation)
        if centr:
            print("Entrainment coefficient map\t%s" % centr, file=documentation)
        else:
            print("Entrainment coefficient map\tFALSE", file=documentation)
        if cvshear:
            print("Shear velocity coefficient map\t%s" % cvshear, file=documentation)
        else:
            print("Shear velocity coefficient map\tFALSE", file=documentation)
        if trelease:
            print("Release time map\t%s" % trelease, file=documentation)
        else:
            print("Release time map\tFALSE", file=documentation)
        if trelstop:
            print("Release hydrograph stop time map\t%s" % trelstop, file=documentation)
        else:
            print("Release hydrograph stop time map\tFALSE", file=documentation)
        if stoptime:
            print("Stopping time map\t%s" % stoptime, file=documentation)
        else:
            print("Stopping time map\tFALSE", file=documentation)
        if tslide:
            print("Time of initial sliding map\t%s" % tslide, file=documentation)
        else:
            print("Time of initial sliding map\tFALSE", file=documentation)

        if hydrograph:
            print("Hydrograph\t%s" % hydrograph, file=documentation)
        else:
            print("Hydrograph\tFALSE", file=documentation)
        if hydrocoords:
            print(
                "Hydrograph profile coordinates\t%s" % hydrocoords, file=documentation
            )
            print("Number of input hydrographs\t%i" % hyd_nin, file=documentation)
            print("Number of output hydrographs\t%i" % hyd_nout, file=documentation)
        else:
            print("Hydrograph profile coordinates\tFALSE", file=documentation)

        if frictiograph:
            print("Frictiograph\t%s" % frictiograph, file=documentation)
        else:
            print("Frictiograph\tFALSE", file=documentation)
        if transformograph:
            print("Transformograph\t%s" % transformograph, file=documentation)
        else:
            print("Transformograph\tFALSE", file=documentation)

        print("Flow parameters\t%s" % flowparam, file=documentation)

        print("Time interval for output\t%s" % tint, file=documentation)
        print("Time at which to stop simulation\t%s" % tstop, file=documentation)
        print("Factor for slow motion\t%s" % slomo, file=documentation)

        print("Search radius for initial sliding\t%s" % slidepar[0], file=documentation)
        print(
            "Exponent for weighting for initial sliding\t%s" % slidepar[1],
            file=documentation,
        )
        print(
            "Coefficient for deformation during initial sliding\t%s" % slidepar[2],
            file=documentation,
        )

        print("CFL criterion\t%s" % cfl[0], file=documentation)
        print("Minimum time step length\t%s" % cfl[1], file=documentation)

        print(
            "Minimum flow height for simulation\t%s" % thresholdsc, file=documentation
        )
        print("Minimum flow height for display\t%s" % thresholds[0], file=documentation)
        print(
            "Minimum flow kinetic energy for display\t%s" % thresholds[1],
            file=documentation,
        )
        print(
            "Minimum flow pressure for display\t%s" % thresholds[2], file=documentation
        )

        if ctrlpoints:
            print("Coordinates of control points\t%s" % ctrlpoints, file=documentation)
        else:
            print("Coordinates of control points\tFALSE", file=documentation)
        if profile:
            print("Profile coordinates\t%s" % profile, file=documentation)
        else:
            print("Profile coordinates\tFALSE", file=documentation)
        if phexagg:
            print(
                "Vertical exaggeration factor for profile display\t%s" % phexagg,
                file=documentation,
            )
        else:
            print(
                "Vertical exaggeration factor for profile display\tFALSE",
                file=documentation,
            )

        if not mflag:
            print("Number of output time steps\t%i" % ntimesteps, file=documentation)
        else:
            print("Number of output time steps\tFALSE", file=documentation)
        if mflag:
            print("Success of simulation\t%i" % nsuccess, file=documentation)
        else:
            print("Success of simulation\t%i" % csuccess, file=documentation)
        if mflag:
            print("Maximum flow velocity\t-9999", file=documentation)
        else:
            print("Maximum flow velocity\t%s" % maxvflow, file=documentation)

        documentation.close()  # closing documentation file

    if vflag:  # evaluation and visualization mode:

        print()
        print("3. EVALUATION AND VISUALIZATION")
        print()

        # Reading documentation file for flags and options

        if not eflag:  # if mode only includes visualization:

            documentation = open(
                pf + "_results/" + pf + "_files/" + pf + "_documentation.txt", "r"
            )  # opening documentation file

            documentation.readline()
            documentation.readline()
            documentation.readline()
            kflag = "TRUE"
            mflag = (
                documentation.readline()
                .replace("Multiple model runs\t", "")
                .replace("\n", "")
            )
            if mflag == "FALSE":
                mflag = None
            documentation.readline()
            documentation.readline()

            pf = documentation.readline().replace("Prefix\t", "").replace("\n", "")

            cellsize = (
                documentation.readline().replace("Cell size\t", "").replace("\n", "")
            )
            rnorth = (
                documentation.readline()
                .replace("Northern boundary\t", "")
                .replace("\n", "")
            )
            rsouth = (
                documentation.readline()
                .replace("Southern boundary\t", "")
                .replace("\n", "")
            )
            rwest = (
                documentation.readline()
                .replace("Western boundary\t", "")
                .replace("\n", "")
            )
            reast = (
                documentation.readline()
                .replace("Eastern boundary\t", "")
                .replace("\n", "")
            )

            model = int(
                documentation.readline().replace("Model\t", "").replace("\n", "")
            )

            nphases = int(
                documentation.readline()
                .replace("Number of phases\t", "")
                .replace("\n", "")
            )

            phasetext = []
            for i in range(0, nphases):
                if i == 0:
                    phasetext.append(
                        documentation.readline()
                        .replace("Phase 1\t", "")
                        .replace("\n", "")
                    )
                if i == 1:
                    phasetext.append(
                        documentation.readline()
                        .replace("Phase 2\t", "")
                        .replace("\n", "")
                    )
                if i == 2:
                    phasetext.append(
                        documentation.readline()
                        .replace("Phase 3\t", "")
                        .replace("\n", "")
                    )

            documentation.readline()
            limiter = (
                documentation.readline()
                .replace("Numerical limiter\t", "")
                .replace("\n", "")
            )
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()

            documentation.readline()
            documentation.readline()

            orthophoto = (
                documentation.readline().replace("Orthophoto\t", "").replace("\n", "")
            )
            if orthophoto == "FALSE":
                orthophoto = None
            impactarea = (
                documentation.readline()
                .replace("Observed impact area\t", "")
                .replace("\n", "")
            )
            if impactarea == "FALSE":
                impactarea = None
            hdeposit = (
                documentation.readline()
                .replace("Height of observed deposit\t", "")
                .replace("\n", "")
            )
            if hdeposit == "FALSE":
                hdeposit = None
            elevation = (
                documentation.readline().replace("Elevation\t", "").replace("\n", "")
            )

            hrelease = (
                documentation.readline()
                .replace("Release height\t", "")
                .replace("\n", "")
            )
            if model == 7:
                hrelease1 = (
                    documentation.readline()
                    .replace("PHASE 1 release height\t", "")
                    .replace("\n", "")
                )
                hrelease2 = (
                    documentation.readline()
                    .replace("PHASE 2 release height\t", "")
                    .replace("\n", "")
                )
                documentation.readline()
                hrelease3 = (
                    documentation.readline()
                    .replace("PHASE 3 release height\t", "")
                    .replace("\n", "")
                )
            documentation.readline()
            releasemass = int(
                documentation.readline()
                .replace("Use of release mass\t", "")
                .replace("\n", "")
            )

            hentrmax = (
                documentation.readline()
                .replace("Maximum entrainable height\t", "")
                .replace("\n", "")
            )
            if model == 7:
                hentrmax1 = (
                    documentation.readline()
                    .replace("Maximum entrainable PHASE 1 height\t", "")
                    .replace("\n", "")
                )
                hentrmax2 = (
                    documentation.readline()
                    .replace("Maximum entrainable PHASE 2 height\t", "")
                    .replace("\n", "")
                )
                documentation.readline()
                hentrmax3 = (
                    documentation.readline()
                    .replace("Maximum entrainable PHASE 3 height\t", "")
                    .replace("\n", "")
                )
            documentation.readline()

            documentation.readline()
            documentation.readline()

            if model == 7:
                documentation.readline()
                documentation.readline()
                documentation.readline()
                documentation.readline()

            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()

            hydrograph = (
                documentation.readline().replace("Hydrograph\t", "").replace("\n", "")
            )
            if hydrograph == "FALSE":
                hydrograph = None
            hydrocoords = (
                documentation.readline()
                .replace("Hydrograph profile coordinates\t", "")
                .replace("\n", "")
            )
            if not hydrocoords == "FALSE":
                hyd_nin = int(
                    documentation.readline()
                    .replace("Number of input hydrographs\t", "")
                    .replace("\n", "")
                )
                hyd_nout = int(
                    documentation.readline()
                    .replace("Number of output hydrographs\t", "")
                    .replace("\n", "")
                )
            else:
                hydrocoords = None
                hyd_nin = 0
                hyd_nout = 0

            documentation.readline()
            documentation.readline()
            documentation.readline()

            tint = (
                documentation.readline()
                .replace("Time interval for output\t", "")
                .replace("\n", "")
            )
            tstop = (
                documentation.readline()
                .replace("Time at which to stop simulation\t", "")
                .replace("\n", "")
            )
            slomo = (
                documentation.readline()
                .replace("Factor for slow motion\t", "")
                .replace("\n", "")
            )
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            documentation.readline()
            thresholds = [
                documentation.readline()
                .replace("Minimum flow height for display\t", "")
                .replace("\n", "")
            ]
            thresholds.append(
                documentation.readline()
                .replace("Minimum flow kinetic energy for display\t", "")
                .replace("\n", "")
            )
            thresholds.append(
                documentation.readline()
                .replace("Minimum flow pressure for display\t", "")
                .replace("\n", "")
            )
            thresholds.append("0.001")

            ctrlpoints = (
                documentation.readline()
                .replace("Coordinates of control points\t", "")
                .replace("\n", "")
            )
            if ctrlpoints == "FALSE":
                ctrlpoints = None
            profile = (
                documentation.readline()
                .replace("Profile coordinates\t", "")
                .replace("\n", "")
            )
            if profile == "FALSE":
                profile = None
            phexagg = (
                documentation.readline()
                .replace("Vertical exaggeration factor for profile display\t", "")
                .replace("\n", "")
            )

            ntimesteps = (
                documentation.readline()
                .replace("Number of output time steps\t", "")
                .replace("\n", "")
            )
            success = int(
                documentation.readline()
                .replace("Success of simulation\t", "")
                .replace("\n", "")
            )
            if not ntimesteps == "FALSE":
                ntimesteps = int(ntimesteps)
            maxvflow = float(
                documentation.readline()
                .replace("Maximum flow velocity\t", "")
                .replace("\n", "")
            )

            documentation.close()  # closing documentation file

            if model <= 3:
                phasetext.append("NA")
            if model <= 3:
                phasetext.append("NA")

            grass.run_command("g.region", flags="d")  # setting default region
            grass.run_command(
                "g.region", n=rnorth, s=rsouth, w=rwest, e=reast
            )  # updating region
            grass.run_command("g.region", res=cellsize)  # updating cell size

            if impactarea:  # if observed impact area is defined:
                grass.run_command(
                    "g.copy",
                    rast=(impactarea, "_impactarea"),
                    overwrite=True,
                    quiet=True,
                )  # temporary impact area map

                maximpactarea = grass.raster_info(impactarea)[
                    "max"
                ]  # maximum value of impact area
                if maximpactarea > 0:
                    impdef = "1"  # setting control for impact area map
                else:
                    impdef = "0"
            else:
                impdef = (
                    "0"
                )  # if observed impact area map is not defined, setting control to negative

            if hdeposit:  # if height of deposit map is defined:
                grass.run_command(
                    "g.copy", rast=(hdeposit, "_hdeposit"), overwrite=True, quiet=True
                )  # temporary height of deposit map

                maxhdeposit = grass.raster_info(hdeposit)[
                    "max"
                ]  # maximum value of height of deposit
                if maxhdeposit > 0:
                    depdef = "1"  # setting control for deposit map
                else:
                    depdef = "0"
            else:
                depdef = (
                    "0"
                )  # if height of deposit map is not defined, setting control to negative

        else:  # success of simulation:
            if mflag:
                success = nsuccess
            else:
                success = csuccess

        if not orthophoto:
            orthophoto = "0"  # orthophoto for visualization

        # Profile plots

        if profile and not mflag and success > 0:

            maxhflow = grass.raster_info(pf + "_hflow_max")[
                "max"
            ]  # maximum flow height

            for j in [
                1,
                2,
                4,
            ]:  # loop over all parameters to be displayed as bar plots:

                if model <= 3:
                    maxsflow1 = grass.raster_info(pf + mstring[j] + "_max")[
                        "max"
                    ]  # maximum flow parameter
                    maxsflow2 = 0.0
                    maxsflow3 = 0.0
                elif model == 7:
                    maxsflow1 = grass.raster_info(pf + mstring[j] + "1_max")[
                        "max"
                    ]  # maximum PHASE 1 flow parameter
                    maxsflow2 = grass.raster_info(pf + mstring[j] + "2_max")[
                        "max"
                    ]  # maximum PHASE 2 flow parameter
                    if model == 7:
                        maxsflow3 = grass.raster_info(pf + mstring[j] + "3_max")[
                            "max"
                        ]  # maximum PHASE 3 flow parameter
                    else:
                        maxsflow3 = 0.0

                d_hdisp = max(
                    maxsflow1, maxsflow2, maxsflow3
                )  # unit conversion factor for scaling of plot:
                if not d_hdisp == 0:
                    mconv = math.log10(d_hdisp)
                else:
                    mconv = 3

                if mconv < 3:
                    mconv2 = "1"
                elif mconv < 6:
                    mconv2 = "0.001"
                elif mconv < 9:
                    mconv2 = "0.000001"
                elif mconv < 12:
                    mconv2 = "0.000000001"
                elif mconv < 15:
                    mconv2 = "0.000000000001"
                else:
                    mconv2 = "0.000000000000001"

                maxsflow1 = maxsflow1 * float(
                    mconv2
                )  # converting units of maxima for second y axis
                maxsflow2 = maxsflow2 * float(mconv2)
                maxsflow3 = maxsflow3 * float(mconv2)

                pnames = []  # initializing list of map images
                pnamesc = []  # initializing list of compressed profile images
                pwidthc = 400  # defining width of compressed profile images

                for step in range(0, ntimesteps + 1):  # loop over all time steps:

                    if step < 10:
                        fill = "000" + str(step)  # formatting time step string
                    elif step < 100:
                        fill = "00" + str(step)
                    elif step < 1000:
                        fill = "0" + str(step)
                    else:
                        fill = str(step)

                    subprocess.call(
                        "Rscript %s/r.avaflow.profile.R %s %s %i %s %s %s %s %s %s %s %s %s %s %s %s --slave --quiet"
                        % (
                            scriptpath,
                            pf,
                            str(model),
                            step,
                            cellsize,
                            maxhflow,
                            maxsflow1,
                            maxsflow2,
                            maxsflow3,
                            phexagg,
                            fill,
                            tint,
                            tstop,
                            depdef,
                            str(j),
                            mconv2,
                        ),
                        shell=True,
                    )  # creating profiles with R

                    pimg = Image.open(
                        pf
                        + "_results/"
                        + pf
                        + "_plots/"
                        + pf
                        + "_profiles_timesteps/"
                        + pf
                        + mstring[j]
                        + fill
                        + ".png"
                    )
                    # opening profile image
                    pheightc = int(
                        float(pimg.size[1]) * float(pwidthc / float(pimg.size[0]))
                    )  # height of compressed profile image
                    pimgc = pimg.resize(
                        (pwidthc, pheightc), Image.BILINEAR
                    )  # compressed profile image
                    pimgc.save(
                        temppath + "/" + mstring[j] + "c" + fill + ".png"
                    )  # saving compressed profile image

                    pnames.append(
                        pf
                        + "_results/"
                        + pf
                        + "_plots/"
                        + pf
                        + "_profiles_timesteps/"
                        + pf
                        + mstring[j]
                        + fill
                        + ".png"
                    )
                    # updating list of profile images
                    pnamesc.append(
                        temppath + "/" + mstring[j] + "c" + fill + ".png"
                    )  # updating list of profile images

                images = []

                for nnn in pnames:
                    frame = Image.open(nnn)
                    frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                    images.append(frame)

                # Save the frames as an animated GIF
                images[0].save(
                    pf
                    + "_results/"
                    + pf
                    + "_plots/"
                    + pf
                    + mstring[j]
                    + "_profile.gif",
                    save_all=True,
                    append_images=images[1:],
                    duration=200,
                    loop=0,
                )

                images = []

                for nnn in pnamesc:
                    frame = Image.open(nnn)
                    frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                    images.append(frame)

                # Save the frames as an animated GIF
                images[0].save(
                    pf
                    + "_results/"
                    + pf
                    + "_plots/"
                    + pf
                    + mstring[j]
                    + "_profilec.gif",
                    save_all=True,
                    append_images=images[1:],
                    duration=200,
                    loop=0,
                )

        # Hydrograph plots

        if (hydrograph or hydrocoords) and not mflag and success > 0:

            for ihyd in range(0, hyd_nin + hyd_nout):  # loop over all hydrographs:

                subprocess.call(
                    "Rscript %s/r.avaflow.hydrograph.R %s %s %s %s %s --slave --quiet"
                    % (
                        scriptpath,
                        pf,
                        str(model),
                        str(ihyd + 1),
                        str(hyd_nin),
                        str(tint),
                    ),
                    shell=True,
                )  # creating hydrograph plot with R

        # Map plots

        if (
            profile and success > 0
        ):  # if profile is defined, writing coordinates to file:

            profcoord = list(
                map(str, profile.split(","))
            )  # splitting profile coordinate string
            proflength = len(profcoord)  # reading number of profile coordinates
            for profnum in range(0, proflength, 2):  # loop over all coordinates:
                if profnum == 0:
                    profilex = profcoord[profnum]  # extracting x coordinate
                    profiley = profcoord[profnum + 1]  # extracting y coordinate
                else:
                    profilex = (
                        profilex + "\t" + profcoord[profnum]
                    )  # extracting x coordinate
                    profiley = (
                        profiley + "\t" + profcoord[profnum + 1]
                    )  # extracting y coordinate

            fmapprof = open(
                temppath + "/fmapprof.txt", "w"
            )  # opening text file for writing
            print(
                "%s\n%s" % (profilex, profiley), file=fmapprof
            )  # writing profile string to text file
            fmapprof.close()  # closing file
            profdef = "1"  # control for profile

        elif success > 0:
            profdef = "0"

        if ctrlpoints and success > 0:
            ctrlpts = str(lctrlpts)  # control for control points
        elif success > 0:
            ctrlpts = "0"

        if (hydrograph or hydrocoords) and success > 0:
            hydr = "1"  # control for hydrographs
        elif success > 0:
            hydr = "0"

        if impdef == "1" and success > 0:  # if observed impact area map is defined:

            grass.mapcalc(
                '"_limpactarea"=if("%s_impactarea"==1,1,null())' % pf,
                overwrite=True,
                quiet=True,
            )  # binary map of impact area
            grass.run_command(
                "r.to.vect",
                input="_limpactarea",
                output="d_limpactarea",
                type="area",
                flags="v" + "s",
                quiet=True,
                overwrite=True,
            )  # line vector of observed impact area
            grass.run_command(
                "v.out.ogr",
                input="d_limpactarea",
                output=temppath,
                format="ESRI_Shapefile",
                type="boundary",
                flags="c",
                quiet=True,
                overwrite=True,
            )  # exporting line vector of observed impact area to shape file

        if (
            depdef == "1" and success > 0
        ):  # if height of observed deposit map is defined:

            grass.mapcalc(
                '"_ldeposit"=if("%s_hdeposit">0,1,null())' % pf,
                overwrite=True,
                quiet=True,
            )  # binary map of observed deposit
            grass.run_command(
                "r.to.vect",
                input="_ldeposit",
                output="d_ldeposit",
                type="area",
                flags="v" + "s",
                quiet=True,
                overwrite=True,
            )  # line vector of observed deposit
            grass.run_command(
                "v.out.ogr",
                input="d_ldeposit",
                output=temppath,
                format="ESRI_Shapefile",
                type="boundary",
                flags="c",
                quiet=True,
                overwrite=True,
            )  # exporting line vector of observed deposit to shape file

        releasemass1 = 0
        if (
            releasemass == 1 and success > 0
        ):  # if (PHASE 1) release height map is defined:

            if model <= 3:  # name of (PHASE 1) release height map:
                if not mflag:
                    hrlss = pf + "_hflow0000"
                else:
                    hrlss = pf + "_hflow00001"
            elif model == 7:
                if not mflag:
                    hrlss = pf + "_hflow10000"
                else:
                    hrlss = pf + "_hflow100001"

            if grass.raster_info(hrlss)["max"] > 0:
                releasemass1 = 1
                grass.mapcalc(
                    '"_lrelease"=if("%s">0,1,null())' % hrlss,
                    overwrite=True,
                    quiet=True,
                )  # binary map of (PHASE 1) release height
                grass.run_command(
                    "r.to.vect",
                    input="_lrelease",
                    output="d_lrelease",
                    type="area",
                    flags="v" + "s",
                    quiet=True,
                    overwrite=True,
                )  # line vector of (PHASE 1) release height
                grass.run_command(
                    "v.out.ogr",
                    input="d_lrelease",
                    output=temppath,
                    format="ESRI_Shapefile",
                    type="boundary",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )  # exporting line vector of (PHASE 1) release height

        releasemass2 = 0
        if (
            model == 7 and releasemass == 1 and success > 0
        ):  # if PHASE 2 release height map is defined:

            if not mflag:
                hrlsf = pf + "_hflow20000"  # name of PHASE 2 release height map:
            else:
                hrlsf = pf + "_hflow200001"
            if grass.raster_info(hrlsf)["max"] > 0:
                releasemass2 = 1
                grass.mapcalc(
                    '"_lrelease2"=if("%s">0,1,null())' % hrlsf,
                    overwrite=True,
                    quiet=True,
                )  # binary map of PHASE 2 release height
                grass.run_command(
                    "r.to.vect",
                    input="_lrelease2",
                    output="d_lrelease2",
                    type="area",
                    flags="v" + "s",
                    quiet=True,
                    overwrite=True,
                )  # line vector of PHASE 2 release height
                grass.run_command(
                    "v.out.ogr",
                    input="d_lrelease2",
                    output=temppath,
                    format="ESRI_Shapefile",
                    type="boundary",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )  # exporting line vector of PHASE 2 release height

        releasemass3 = 0
        if (
            model == 7 and releasemass == 1 and success > 0
        ):  # if PHASE 3 release height map is defined:

            if not mflag:
                hrlsw = pf + "_hflow30000"  # name of PHASE 3 release height map:
            else:
                hrlsw = pf + "_hflow300001"
            if grass.raster_info(hrlsw)["max"] > 0:
                releasemass3 = 1
                grass.mapcalc(
                    '"_lrelease3"=if("%s">0,1,null())' % hrlsw,
                    overwrite=True,
                    quiet=True,
                )  # binary map of PHASE 3 release height
                grass.run_command(
                    "r.to.vect",
                    input="_lrelease3",
                    output="d_lrelease3",
                    type="area",
                    flags="v" + "s",
                    quiet=True,
                    overwrite=True,
                )  # line vector of PHASE 3 release height
                grass.run_command(
                    "v.out.ogr",
                    input="d_lrelease3",
                    output=temppath,
                    format="ESRI_Shapefile",
                    type="boundary",
                    flags="c",
                    quiet=True,
                    overwrite=True,
                )  # exporting line vector of PHASE 3 release height

        if not mflag and success > 0:  # for single model run:

            maxhentr = grass.raster_info(pf + "_basechange_fin")[
                "max"
            ]  # maximum height of entrainment
            maxhdep = (
                grass.raster_info(pf + "_basechange_fin")["min"] * -1
            )  # maximum height of deposition

            jrange = [0]
            if aflag:
                jrange.append(1)
                jrange.append(2)
            if maxhentr > 0 or maxhdep > 0:
                jrange.append(3)
            if model == 7 and tflag:
                jrange.append(5)
            jrange.append(6)

            for j in jrange:  # loop over all sets of maps to be created

                if j == 6:
                    d_hdisp = grass.raster_info(pf + "_treach")[
                        "max"
                    ]  # absolute maximum raster value
                    d_hdispn = 0  # absolute minimum raster value
                elif j < 3:
                    d_hdisp = grass.raster_info(pf + mstring[j] + "_max")[
                        "max"
                    ]  # absolute maximum raster value
                    d_hdispn = 0  # absolute minimum raster value
                elif j == 3:
                    d_hdisp = grass.raster_info(pf + mstring[j] + "_fin")["max"]
                    d_hdispn = grass.raster_info(pf + mstring[j] + "_fin")["min"]
                else:
                    grass.mapcalc(
                        '"%s_tsun_max"="%s"-"%s"'
                        % (pf, pf + "_hflow3_max", pf + "_hflow30000"),
                        overwrite=True,
                        quiet=True,
                    )
                    d_hdisp = grass.raster_info(pf + "_tsun_max")["max"]
                    d_hdispn = d_hdisp * -1

                if j == 6:
                    mconv = (
                        1
                    )  # unit conversion and number of digits for legend according to maximum value:
                elif not max(d_hdisp, -d_hdispn) == 0:
                    mconv = math.log10(max(d_hdisp, -d_hdispn))
                else:
                    mconv = 3

                if mconv < 3:
                    mconv2 = "1"
                    if mconv < 0:
                        mdig = 3
                    elif mconv < 1:
                        mdig = 2
                    elif mconv < 2:
                        mdig = 1
                    else:
                        mdig = 0
                elif mconv < 6:
                    mconv2 = "0.001"
                    if mconv < 4:
                        mdig = 2
                    elif mconv < 5:
                        mdig = 1
                    else:
                        mdig = 0
                elif mconv < 9:
                    mconv2 = "0.000001"
                    if mconv < 7:
                        mdig = 2
                    elif mconv < 8:
                        mdig = 1
                    else:
                        mdig = 0
                elif mconv < 12:
                    mconv2 = "0.000000001"
                    if mconv < 10:
                        mdig = 2
                    elif mconv < 11:
                        mdig = 1
                    else:
                        mdig = 0
                elif mconv < 15:
                    mconv2 = "0.000000000001"
                    if mconv < 13:
                        mdig = 2
                    elif mconv < 14:
                        mdig = 1
                    else:
                        mdig = 0
                else:
                    mconv2 = "0.000000000000001"
                    if mconv < 16:
                        mdig = 2
                    elif mconv < 17:
                        mdig = 1
                    else:
                        mdig = 0
                d_hdisp = d_hdisp * float(mconv2)

                if model <= 3 and j == 3:  # for entrainment and deposition map:
                    d_hdispn = d_hdispn * float(mconv2)
                    jimp = 0  # index of impact parameter
                    d_htrsh5 = "0"
                    d_htrsh4 = str(
                        "{:.20f}".format(round(d_hdisp, mdig))
                    )  # breaks for raster maps
                    d_htrsh3 = str("{:.20f}".format(round(d_hdisp * 10 / 20, mdig)))
                    d_htrsh2 = str("{:.20f}".format(round(d_hdisp * 6 / 20, mdig)))
                    d_htrsh1 = str("{:.20f}".format(round(d_hdisp * 3 / 20, mdig)))
                    if d_hdisp == 0:
                        d_htrshm = str(round(0, mdig))
                    else:
                        d_htrshm = str(
                            min(
                                float(thresholds[jimp]) * float(mconv2),
                                float(d_htrsh1) * 0.75,
                            )
                        )
                    # minimum value displayed

                    d_htrsh5n = "0"
                    d_htrsh4n = str("{:.20f}".format(round(d_hdispn, mdig)))
                    d_htrsh3n = str("{:.20f}".format(round(d_hdispn * 10 / 20, mdig)))
                    d_htrsh2n = str("{:.20f}".format(round(d_hdispn * 6 / 20, mdig)))
                    d_htrsh1n = str("{:.20f}".format(round(d_hdispn * 3 / 20, mdig)))
                    if d_hdispn == 0:
                        d_htrshmn = str(round(0, mdig))
                    else:
                        d_htrshmn = str(
                            max(
                                -1 * float(thresholds[jimp]) * float(mconv2),
                                float(d_htrsh1n) * 0.75,
                            )
                        )
                    # minimum value displayed

                elif model == 7 and j == 5:  # for tsunami map:
                    d_hdispn = d_hdispn * float(mconv2)
                    d_htrsh5 = "0"
                    d_htrsh4 = str(
                        "{:.20f}".format(round(d_hdisp, mdig))
                    )  # breaks for raster maps
                    d_htrsh3 = str("{:.20f}".format(round(d_hdisp * 15 / 20, mdig)))
                    d_htrsh2 = str("{:.20f}".format(round(d_hdisp * 10 / 20, mdig)))
                    d_htrsh1 = str("{:.20f}".format(round(d_hdisp * 5 / 20, mdig)))
                    d_htrshm = str("{:.20f}".format(round(d_hdisp * 1 / 25, mdig)))

                    d_htrsh5n = "0"
                    d_htrsh4n = str("{:.20f}".format(round(d_hdispn, mdig)))
                    d_htrsh3n = str("{:.20f}".format(round(d_hdispn * 15 / 20, mdig)))
                    d_htrsh2n = str("{:.20f}".format(round(d_hdispn * 10 / 20, mdig)))
                    d_htrsh1n = str("{:.20f}".format(round(d_hdispn * 5 / 20, mdig)))
                    d_htrshmn = str("{:.20f}".format(round(d_hdispn * 1 / 25, mdig)))

                elif j == 6:  # for time of reach map:
                    jimp = j
                    d_htrsh5 = str(
                        "{:.20f}".format(round(d_hdisp, mdig))
                    )  # breaks for raster maps
                    d_htrsh4 = str("{:.20f}".format(round(d_hdisp * 10 / 20, mdig)))
                    d_htrsh3 = str("{:.20f}".format(round(d_hdisp * 7 / 20, mdig)))
                    d_htrsh2 = str("{:.20f}".format(round(d_hdisp * 4 / 20, mdig)))
                    d_htrsh1 = str("{:.20f}".format(round(d_hdisp * 2 / 20, mdig)))
                    d_htrshm = str(round(0, mdig))

                    d_htrsh5n = "0"
                    d_htrsh4n = "0"
                    d_htrsh3n = "0"
                    d_htrsh2n = "0"
                    d_htrsh1n = "0"
                    d_htrshmn = "0"

                else:  # for all other maps:
                    jimp = j
                    d_htrsh5 = str("{:.20f}".format(round(d_hdisp, mdig)))  # breaks
                    d_htrsh4 = str("{:.20f}".format(round(d_hdisp * 12 / 20, mdig)))
                    d_htrsh3 = str("{:.20f}".format(round(d_hdisp * 8 / 20, mdig)))
                    d_htrsh2 = str("{:.20f}".format(round(d_hdisp * 5 / 20, mdig)))
                    d_htrsh1 = str("{:.20f}".format(round(d_hdisp * 2 / 20, mdig)))
                    if d_hdisp == 0:
                        d_htrshm = str(round(0, mdig))
                    else:
                        d_htrshm = str(
                            min(
                                float(thresholds[jimp]) * float(mconv2),
                                float(d_htrsh1) * 0.75,
                            )
                        )
                    # minimum value displayed

                    if model == 7 and j == 3:

                        d_hdispn = d_hdispn * float(mconv2)
                        d_htrsh5n = str(
                            "{:.20f}".format(round(d_hdispn, mdig))
                        )  # breaks (negative values)
                        d_htrsh4n = str(
                            "{:.20f}".format(round(d_hdispn * 12 / 20, mdig))
                        )
                        d_htrsh3n = str(
                            "{:.20f}".format(round(d_hdispn * 8 / 20, mdig))
                        )
                        d_htrsh2n = str(
                            "{:.20f}".format(round(d_hdispn * 5 / 20, mdig))
                        )
                        d_htrsh1n = str(
                            "{:.20f}".format(round(d_hdispn * 2 / 20, mdig))
                        )
                        if d_hdisp == 0:
                            d_htrshmn = str(round(0, mdig))
                        else:
                            d_htrshmn = str(
                                max(
                                    -float(thresholds[jimp]) * float(mconv2),
                                    float(d_htrsh1n) * 0.75,
                                )
                            )
                        # minimum value displayed

                    else:

                        d_htrsh5n = "0"
                        d_htrsh4n = "0"
                        d_htrsh3n = "0"
                        d_htrsh2n = "0"
                        d_htrsh1n = "0"
                        d_htrshmn = "0"

                mnames = []  # initializing list of map images
                mnamesc = []  # initializing list of compressed map images
                mwidthc = 640  # defining width of compressed map images

                for i in range(
                    0, ntimesteps + 2
                ):  # loop over all time steps plus one for maps of maximum values:

                    if i < 10:
                        fill = "000" + str(i)  # formatting time step string
                    elif i < 100:
                        fill = "00" + str(i)
                    elif i < 1000:
                        fill = "0" + str(i)
                    else:
                        fill = str(i)

                    if i == 0 and model <= 3:
                        grass.mapcalc(
                            '"_elmod"="%s_elev"+"%s_hflow0000"' % (pf, pf),
                            overwrite=True,
                            quiet=True,
                        )  # elevation map
                    elif i == 0:
                        grass.mapcalc(
                            '"_elmod"="%s_elev"+"%s_hflow10000"' % (pf, pf),
                            overwrite=True,
                            quiet=True,
                        )
                    elif i == ntimesteps + 1:
                        grass.mapcalc(
                            '"_elmod"="%s_elev"+"%s_basechange_fin"' % (pf, pf),
                            overwrite=True,
                            quiet=True,
                        )
                    elif model <= 3:
                        grass.mapcalc(
                            '"_elmod"="%s_elev"+"%s_hflow%s"+"%s_basechange%s"'
                            % (pf, pf, fill, pf, fill),
                            overwrite=True,
                            quiet=True,
                        )
                    else:
                        grass.mapcalc(
                            '"_elmod"="%s_elev"+"%s_hflow1%s"+"%s_basechange%s"'
                            % (pf, pf, fill, pf, fill),
                            overwrite=True,
                            quiet=True,
                        )
                    grass.run_command(
                        "r.relief",
                        input="_elmod",
                        output="_hillshade",
                        azimuth="315",
                        quiet=True,
                        overwrite=True,
                    )  # hillshade
                    grass.run_command(
                        "r.out.gdal",
                        input="_hillshade",
                        output=temppath + "/hillshade.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # export hillshade map to ascii

                    if (
                        i <= ntimesteps and model == 7 and j == 5
                    ):  # strings for names of maps:
                        mstringt = pf + "_tsun_" + fill
                        grass.mapcalc(
                            '"%s"="%s"-"%s"'
                            % (mstringt, pf + "_hflow3" + fill, pf + "_hflow30000"),
                            overwrite=True,
                            quiet=True,
                        )
                    elif i <= ntimesteps and not j == 6:
                        mstringt = pf + mstring[j] + fill
                        mstrings = pf + mstring[j] + "1" + fill  # !!!CHECK
                        if model == 7:
                            mstringf = pf + mstring[j] + "2" + fill
                        if model == 7:
                            mstringw = pf + mstring[j] + "3" + fill
                    elif j < 3:
                        mstringt = pf + mstring[j] + "_max"
                        mstrings = pf + mstring[j] + "1_max"
                        if model == 7:
                            mstringf = pf + mstring[j] + "2_max"
                        if model == 7:
                            mstringw = pf + mstring[j] + "3_max"
                    elif j == 3:
                        mstringt = pf + mstring[j] + "_fin"
                        mstrings = pf + mstring[j] + "1_fin"
                        if model == 7:
                            mstringf = pf + mstring[j] + "2_fin"
                        if model == 7:
                            mstringw = pf + mstring[j] + "3_fin"
                    elif model == 7 and j == 5:
                        mstringt = pf + "_tsun_max"
                    elif i == ntimesteps and j == 6:
                        mstringt = pf + "_treach"

                    if model == 7 and j == 5:  # for tsunami map:

                        grass.mapcalc(
                            '"_a2"=if("%s"*%s<%s,1,if("%s"*%s<%s,2,if("%s"*%s<%s,3,if("%s"*%s<%s,4,if("%s"*%s<=%s,5,\
                            if("%s"*%s<=%s,6,if("%s"*%s<=%s,7,if("%s"*%s<=%s,8,9))))))))'
                            % (
                                mstringt,
                                mconv2,
                                d_htrsh3n,
                                mstringt,
                                mconv2,
                                d_htrsh2n,
                                mstringt,
                                mconv2,
                                d_htrsh1n,
                                mstringt,
                                mconv2,
                                d_htrshmn,
                                mstringt,
                                mconv2,
                                d_htrshm,
                                mstringt,
                                mconv2,
                                d_htrsh1,
                                mstringt,
                                mconv2,
                                d_htrsh2,
                                mstringt,
                                mconv2,
                                d_htrsh3,
                            ),
                            overwrite=True,
                            quiet=True,
                        )  # reclassifying map

                        grass.mapcalc(
                            '"_a2c"="_a2"', overwrite=True, quiet=True
                        )  # raster map for contour creation

                        contstep = 1  # contour interval
                        contmin = 0  # minimum contour
                        contmax = 9  # maximum contour

                    elif model <= 3 or (i == ntimesteps and j == 6):

                        if j == 3:  # for entrainment and deposition map:

                            grass.mapcalc(
                                '"_a2"=if("%s"*%s<%s,1,if("%s"*%s<%s,2,if("%s"*%s<%s,3,if("%s"*%s<%s,4,if("%s"*%s<=%s,5,\
                                if("%s"*%s<=%s,6,if("%s"*%s<=%s,7,if("%s"*%s<=%s,8,9))))))))'
                                % (
                                    mstringt,
                                    mconv2,
                                    d_htrsh3n,
                                    mstringt,
                                    mconv2,
                                    d_htrsh2n,
                                    mstringt,
                                    mconv2,
                                    d_htrsh1n,
                                    mstringt,
                                    mconv2,
                                    d_htrshmn,
                                    mstringt,
                                    mconv2,
                                    d_htrshm,
                                    mstringt,
                                    mconv2,
                                    d_htrsh1,
                                    mstringt,
                                    mconv2,
                                    d_htrsh2,
                                    mstringt,
                                    mconv2,
                                    d_htrsh3,
                                ),
                                overwrite=True,
                                quiet=True,
                            )  # reclassifying map

                            grass.mapcalc(
                                '"_a2c"="_a2"', overwrite=True, quiet=True
                            )  # raster map for contour creation

                            contstep = 1  # contour interval
                            contmin = 0  # minimum contour
                            contmax = 9  # maximum contour

                        else:  # for all other maps:

                            grass.mapcalc(
                                '"_a2"=if("%s"*%s<%s,0,if("%s"*%s<%s,1,if("%s"*%s<%s,2,if("%s"*%s<%s,3,if("%s"*%s<%s,4,5)))))'
                                % (
                                    mstringt,
                                    mconv2,
                                    d_htrshm,
                                    mstringt,
                                    mconv2,
                                    d_htrsh1,
                                    mstringt,
                                    mconv2,
                                    d_htrsh2,
                                    mstringt,
                                    mconv2,
                                    d_htrsh3,
                                    mstringt,
                                    mconv2,
                                    d_htrsh4,
                                ),
                                overwrite=True,
                                quiet=True,
                            )  # reclassifying map
                            grass.mapcalc(
                                '"_a2c"="_a2"', overwrite=True, quiet=True
                            )  # raster map for contour creation

                            contstep = 1  # contour interval
                            contmin = 1  # minimum contour
                            contmax = 5  # maximum contour

                    elif (
                        model == 7
                    ):  # create composite map for Pudasaini two-phase model:

                        if (
                            j == 3 and not model == 7
                        ):  # for entrainment and deposition map:

                            grass.mapcalc(
                                '"_a21"="%s"*%s' % (mstringt, mconv2),
                                overwrite=True,
                                quiet=True,
                            )
                            # sum of all phases, unit conversion

                            if model == 7:  # for three-phase model:

                                grass.mapcalc(
                                    '"_a22"=if("%s"==0,0,if("%s"+"%s"==0,1,"%s"/("%s")))'
                                    % (
                                        mstrings,
                                        mstringf,
                                        mstringw,
                                        mstrings,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 1 component
                                grass.mapcalc(
                                    '"_a23"=if("%s"==0,0,if("%s"+"%s"==0,1,"%s"/("%s")))'
                                    % (
                                        mstringf,
                                        mstrings,
                                        mstringw,
                                        mstringf,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 2 component
                                grass.mapcalc(
                                    '"_a24"=if("%s"==0,0,if("%s"+"%s"==0,1,"%s"/("%s")))'
                                    % (
                                        mstringw,
                                        mstrings,
                                        mstringf,
                                        mstringw,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 3 component

                                grass.mapcalc(
                                    '"_a2"=if(abs("_a21")<=%s,0,if(abs("_a21")<=%s,1,if(abs("_a21")<=%s,2,if(abs("_a21")<=%s,3,if(abs("_a21")<=%s,4,5)))))'
                                    % (
                                        d_htrshm,
                                        d_htrsh1,
                                        d_htrsh2,
                                        d_htrsh3,
                                        d_htrsh4,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for magnitude
                                grass.mapcalc(
                                    '"_a2s"=if("_a2"==0,0,if("_a22"<0.2,1,if("_a22"<0.4,2,if("_a22"<0.6,3,\
                                    if("_a22"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 1 ratio
                                grass.mapcalc(
                                    '"_a2f"=if("_a2"==0,0,if("_a23"<0.2,1,if("_a23"<0.4,2,if("_a23"<0.6,3,\
                                    if("_a23"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 2 ratio
                                grass.mapcalc(
                                    '"_a2w"=if("_a2"==0,0,if("_a24"<0.2,1,if("_a24"<0.4,2,if("_a24"<0.6,3,\
                                    if("_a24"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 3 ratio
                                grass.mapcalc(
                                    '"_a2"=("_a2"-1)*125+("_a2s"-1)*25+("_a2f"-1)*5+"_a2w"',
                                    overwrite=True,
                                    quiet=True,
                                )  # combined ratio
                                grass.mapcalc(
                                    '"_a2c"="_a2"+124', overwrite=True, quiet=True
                                )  # raster map for contour creation

                                contstep = 125  # contour interval
                                contmin = 125  # minimum contour
                                contmax = 625  # maximum contour

                        elif not j == 6:  # for all other maps:

                            grass.mapcalc(
                                '"_a21"="%s"*%s' % (mstringt, mconv2),
                                overwrite=True,
                                quiet=True,
                            )
                            # sum of all phases, unit conversion

                            if model == 7:  # for three-phase model:

                                grass.mapcalc(
                                    '"_a22"=if("%s"==0,0,if("%s"+"%s"==0,1,abs("%s"/("%s"))))'
                                    % (
                                        mstrings,
                                        mstringf,
                                        mstringw,
                                        mstrings,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 1 component
                                grass.mapcalc(
                                    '"_a23"=if("%s"==0,0,if("%s"+"%s"==0,1,abs("%s"/("%s"))))'
                                    % (
                                        mstringf,
                                        mstrings,
                                        mstringw,
                                        mstringf,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 2 component
                                grass.mapcalc(
                                    '"_a24"=if("%s"==0,0,if("%s"+"%s"==0,1,abs("%s"/("%s"))))'
                                    % (
                                        mstringw,
                                        mstrings,
                                        mstringf,
                                        mstringw,
                                        mstringt,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # ratio of PHASE 3 component

                                grass.mapcalc(
                                    '"_a2"=if("_a21"<=%s,0,if("_a21"<=%s,1,if("_a21"<=%s,2,if("_a21"<=%s,3,if("_a21"<=%s,4,5)))))'
                                    % (
                                        d_htrshm,
                                        d_htrsh1,
                                        d_htrsh2,
                                        d_htrsh3,
                                        d_htrsh4,
                                    ),
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for magnitude

                                if j == 3:
                                    grass.mapcalc(
                                        '"_a2n"=if("_a21">=%s,0,if("_a21">=%s,1,if("_a21">=%s,2,if("_a21">=%s,3,if("_a21">=%s,4,5)))))'
                                        % (
                                            d_htrshmn,
                                            d_htrsh1n,
                                            d_htrsh2n,
                                            d_htrsh3n,
                                            d_htrsh4n,
                                        ),
                                        overwrite=True,
                                        quiet=True,
                                    )
                                else:
                                    grass.mapcalc(
                                        '"_a2n"=0', overwrite=True, quiet=True
                                    )
                                grass.mapcalc(
                                    '"_a2x"=max("_a2","_a2n")',
                                    overwrite=True,
                                    quiet=True,
                                )

                                grass.mapcalc(
                                    '"_a2s"=if("_a2"==0,0,if("_a22"<0.2,1,if("_a22"<0.4,2,if("_a22"<0.6,3,\
                                    if("_a22"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 1 ratio
                                grass.mapcalc(
                                    '"_a2f"=if("_a2"==0,0,if("_a23"<0.2,1,if("_a23"<0.4,2,if("_a23"<0.6,3,\
                                    if("_a23"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 2 ratio
                                grass.mapcalc(
                                    '"_a2w"=if("_a2"==0,0,if("_a24"<0.2,1,if("_a24"<0.4,2,if("_a24"<0.6,3,\
                                    if("_a24"<0.8,4,5)))))',
                                    overwrite=True,
                                    quiet=True,
                                )  # reclass for PHASE 3 ratio

                                grass.mapcalc(
                                    '"_a2"=("_a2"-1)*125+("_a2s"-1)*25+("_a2f"-1)*5+"_a2w"',
                                    overwrite=True,
                                    quiet=True,
                                )  # combined ratio
                                grass.mapcalc(
                                    '"_a2c"="_a2"+124', overwrite=True, quiet=True
                                )  # raster map for contour creation

                                if j == 3:

                                    grass.mapcalc(
                                        '"_a2sn"=if("_a2n"==0,0,if("_a22"<0.2,1,if("_a22"<0.4,2,if("_a22"<0.6,3,\
                                        if("_a22"<0.8,4,5)))))',
                                        overwrite=True,
                                        quiet=True,
                                    )  # reclass for PHASE 1 ratio
                                    grass.mapcalc(
                                        '"_a2fn"=if("_a2n"==0,0,if("_a23"<0.2,1,if("_a23"<0.4,2,if("_a23"<0.6,3,\
                                        if("_a23"<0.8,4,5)))))',
                                        overwrite=True,
                                        quiet=True,
                                    )  # reclass for PHASE 2 ratio
                                    grass.mapcalc(
                                        '"_a2wn"=if("_a2n"==0,0,if("_a24"<0.2,1,if("_a24"<0.4,2,if("_a24"<0.6,3,\
                                        if("_a24"<0.8,4,5)))))',
                                        overwrite=True,
                                        quiet=True,
                                    )  # reclass for PHASE 3 ratio

                                    grass.mapcalc(
                                        '"_a2n"=("_a2n"-1)*125+("_a2sn"-1)*25+("_a2fn"-1)*5+"_a2wn"',
                                        overwrite=True,
                                        quiet=True,
                                    )  # combined ratio
                                    grass.mapcalc(
                                        '"_a2cn"="_a2n"+124', overwrite=True, quiet=True
                                    )  # raster map for contour creation
                                    grass.mapcalc(
                                        '"_a2"=("_a2x"-1)*125+("_a2s"+"_a2sn"-1)*25+("_a2f"+"_a2fn"-1)*5+"_a2w"+"_a2wn"',
                                        overwrite=True,
                                        quiet=True,
                                    )  # combined ratio

                                contstep = 125  # contour interval
                                contmin = 125  # minimum contour
                                contmax = 625  # maximum contour

                    grass.run_command(
                        "r.out.gdal",
                        input="_a2",
                        output=temppath + "/a2.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # export composite map to ascii raster

                    maxhstep = grass.raster_info("_a2")[
                        "max"
                    ]  # maximum value of composite map
                    minhstep = grass.raster_info("_a2")[
                        "min"
                    ]  # maximum value of composite map
                    if not maxhstep == 0 or not minhstep == 0:

                        grass.run_command(
                            "r.contour",
                            input="_a2c",
                            output="d_contours",
                            step=contstep,
                            minlevel=contmin,
                            maxlevel=contmax,
                            quiet=True,
                            overwrite=True,
                        )
                        # creating contour lines of magnitude
                        grass.run_command(
                            "v.out.ogr",
                            flags="c",
                            input="d_contours",
                            output=temppath,
                            type="line",
                            format="ESRI_Shapefile",
                            quiet=True,
                            overwrite=True,
                        )  # exporting contour lines to shape file

                        if model == 7 and j == 3:

                            grass.run_command(
                                "r.contour",
                                input="_a2cn",
                                output="d_contours_neg",
                                step=contstep,
                                minlevel=contmin,
                                maxlevel=contmax,
                                quiet=True,
                                overwrite=True,
                            )
                            # creating contour lines of magnitude (entrainment)
                            grass.run_command(
                                "v.out.ogr",
                                flags="c",
                                input="d_contours_neg",
                                output=temppath,
                                type="line",
                                format="ESRI_Shapefile",
                                quiet=True,
                                overwrite=True,
                            )  # exporting contour lines to shape file (entrainment)

                    if not (i <= ntimesteps and j == 6):
                        subprocess.call(
                            "Rscript %s/r.avaflow.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                            " %s %s %s %s %s %s %s %s %s %s"
                            % (
                                scriptpath,
                                temppath,
                                pf,
                                str(model),
                                str(i),
                                fill,
                                contstep,
                                d_htrsh5,
                                d_htrsh4,
                                d_htrsh3,
                                d_htrsh2,
                                d_htrsh1,
                                d_htrshm,
                                d_htrsh5n,
                                d_htrsh4n,
                                d_htrsh3n,
                                d_htrsh2n,
                                d_htrsh1n,
                                d_htrshmn,
                                rnorth,
                                rsouth,
                                rwest,
                                reast,
                                cellsize,
                                tint,
                                tstop,
                                impdef,
                                depdef,
                                str(maxvflow),
                                profdef,
                                hydr,
                                str(releasemass1 + 2 * releasemass2 + 4 * releasemass3),
                                str(hyd_nin),
                                str(hyd_nout),
                                mstring[j],
                                mconv2,
                                str(mdig),
                                str(ntimesteps),
                                ctrlpts,
                                orthophoto,
                                phasetext[0],
                                phasetext[1],
                                phasetext[2],
                            ),
                            shell=True,
                        )  # creating maps with R

                    if (
                        i <= ntimesteps and not j == 6
                    ):  # for all time steps (not maps of maximum value), prepare maps for animated gifs

                        mimg = Image.open(
                            pf
                            + "_results/"
                            + pf
                            + "_plots/"
                            + pf
                            + "_maps_timesteps/"
                            + pf
                            + mstring[j]
                            + fill
                            + ".png"
                        )  # opening map image
                        mheightc = int(
                            float(mimg.size[1]) * float(mwidthc / float(mimg.size[0]))
                        )  # height of compressed map image
                        mimgc = mimg.resize(
                            (mwidthc, mheightc), Image.BILINEAR
                        )  # compressed map image
                        mimgc.save(
                            temppath + "/" + mstring[j] + "c" + fill + ".png"
                        )  # saving compressed map image

                        mnames.append(
                            pf
                            + "_results/"
                            + pf
                            + "_plots/"
                            + pf
                            + "_maps_timesteps/"
                            + pf
                            + mstring[j]
                            + fill
                            + ".png"
                        )
                        # updating list of map images
                        mnamesc.append(
                            temppath + "/" + mstring[j] + "c" + fill + ".png"
                        )  # updating list of compressed map images

                if not j == 6:

                    images = []

                    for nnn in mnames:
                        frame = Image.open(nnn)
                        frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                        images.append(frame)

                    # Save the frames as an animated GIF
                    images[0].save(
                        pf
                        + "_results/"
                        + pf
                        + "_plots/"
                        + pf
                        + mstring[j]
                        + "_map.gif",
                        save_all=True,
                        append_images=images[1:],
                        duration=200,
                        loop=0,
                    )

                    images = []

                    for nnn in mnamesc:
                        frame = Image.open(nnn)
                        frame = frame.convert("P", palette=Image.ADAPTIVE, colors=256)
                        images.append(frame)

                    # Save the frames as an animated GIF
                    images[0].save(
                        pf
                        + "_results/"
                        + pf
                        + "_plots/"
                        + pf
                        + mstring[j]
                        + "_mapc.gif",
                        save_all=True,
                        append_images=images[1:],
                        duration=200,
                        loop=0,
                    )

        elif success > 0:  # for multiple model runs:

            d_htrsh5 = "1.00"  # breaks of impact indicator index
            d_htrsh4 = "0.80"
            d_htrsh3 = "0.60"
            d_htrsh2 = "0.40"
            d_htrsh1 = "0.20"
            d_htrshm = "0.01"

            grass.run_command(
                "r.relief",
                input=pf + "_elev",
                output="_hillshade",
                azimuth="315",
                quiet=True,
                overwrite=True,
            )  # hillshade
            grass.run_command(
                "r.out.gdal",
                input="_hillshade",
                output=temppath + "/hillshade.asc",
                format="AAIGrid",
                overwrite=True,
            )  # export hillshade map to ascii

            for j in range(0, 4):  # loop over all sets of maps to be created

                if j < 3:
                    mstringx = "_iii" + mstring[j]
                else:
                    mstringx = "_dii"

                grass.mapcalc(
                    '"_a2"=if("%s%s"<%s,0,if("%s%s"<%s,1,if("%s%s"<%s,2,if("%s%s"<%s,3,if("%s%s"<%s,4,5)))))'
                    % (
                        pf,
                        mstringx,
                        d_htrshm,
                        pf,
                        mstringx,
                        d_htrsh1,
                        pf,
                        mstringx,
                        d_htrsh2,
                        pf,
                        mstringx,
                        d_htrsh3,
                        pf,
                        mstringx,
                        d_htrsh4,
                    ),
                    overwrite=True,
                    quiet=True,
                )  # reclassifying impact indicator index

                grass.run_command(
                    "r.out.gdal",
                    input="_a2",
                    output=temppath + "/a2.asc",
                    format="AAIGrid",
                    overwrite=True,
                )  # export reclassified map to ascii raster

                grass.mapcalc(
                    '"_a2c"="_a2"', overwrite=True, quiet=True
                )  # raster map for contour creation

                contstep = 1  # contour interval
                contmin = 1  # minimum contour
                contmax = 5  # maximum contour
                maxhstep = grass.raster_info("_a2")[
                    "max"
                ]  # maximum value of composite map
                if maxhstep >= contstep:

                    grass.run_command(
                        "r.contour",
                        input="_a2c",
                        output="d_contours",
                        step=contstep,
                        minlevel=contmin,
                        maxlevel=contmax,
                        quiet=True,
                        overwrite=True,
                    )  # creating contour lines

                    grass.run_command(
                        "v.out.ogr",
                        flags="c",
                        input="d_contours",
                        output=temppath,
                        type="line",
                        format="ESRI_Shapefile",
                        quiet=True,
                        overwrite=True,
                    )  # exporting contour lines to shape file

                subprocess.call(
                    "Rscript %s/r.avaflow.map.R %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"
                    " %s %s %s %s %s %s %s %s %s"
                    % (
                        scriptpath,
                        temppath,
                        pf,
                        str(model),
                        "1",
                        "iscore",
                        contstep,
                        d_htrsh5,
                        d_htrsh4,
                        d_htrsh3,
                        d_htrsh2,
                        d_htrsh1,
                        d_htrshm,
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        "0",
                        rnorth,
                        rsouth,
                        rwest,
                        reast,
                        cellsize,
                        tint,
                        tstop,
                        impdef,
                        depdef,
                        "0",
                        profdef,
                        hydr,
                        str(releasemass1 + 2 * releasemass2),
                        str(hyd_nin),
                        str(hyd_nout),
                        mstringx,
                        "1",
                        "2",
                        "0",
                        ctrlpts,
                        orthophoto,
                        phasetext[0],
                        phasetext[1],
                        phasetext[2],
                    ),
                    shell=True,
                )  # creating map with R

                # Evaluation through ROC plots

                if j == 0 and impdef == "1":

                    grass.run_command(
                        "r.out.gdal",
                        flags="c",
                        input=impactarea,
                        output=temppath + "/observed.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting map of observed impact area map to ascii raster
                    grass.run_command(
                        "r.out.gdal",
                        flags="c",
                        input=pf + mstringx,
                        output=temppath + "/xindex.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting impact indicator index map to ascii raster

                    with open(temppath + "/xindex.asc", "r") as xfin:
                        with open(temppath + "/index.asc", "w") as xfout:
                            for line in xfin:
                                xfout.write(line.replace("-nan", "0.0"))

                    xfin.close()
                    xfout.close()

                    subprocess.call(
                        "Rscript %s/r.avaflow.roc.R %s %s %s %s %s --slave --quiet"
                        % (scriptpath, temppath, pf, "1", "1", "1"),
                        shell=True,
                    )
                    # ROC plot relating iii to observation (without normalization)

                    subprocess.call(
                        "Rscript %s/r.avaflow.roc.R %s %s %s %s %s --slave --quiet"
                        % (scriptpath, temppath, pf, "2", "2", "1"),
                        shell=True,
                    )
                    # ROC plot relating iii to observation (with normalization)

                if j == 3 and depdef == "1":

                    grass.run_command(
                        "r.out.gdal",
                        flags="c",
                        input=hdeposit,
                        output=temppath + "/observed.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting map of observed deposition height to ascii raster
                    grass.run_command(
                        "r.out.gdal",
                        flags="c",
                        input=pf + mstringx,
                        output=temppath + "/xindex.asc",
                        format="AAIGrid",
                        overwrite=True,
                    )  # exporting deposition indicator map to ascii raster

                    with open(temppath + "/xindex.asc", "r") as xfin:
                        with open(temppath + "/index.asc", "w") as xfout:
                            for line in xfin:
                                xfout.write(line.replace("-nan", "0.0"))

                    xfin.close()
                    xfout.close()

                    subprocess.call(
                        "Rscript %s/r.avaflow.roc.R %s %s %s %s %s --slave --quiet"
                        % (scriptpath, temppath, pf, "3", "1", "2"),
                        shell=True,
                    )
                    # ROC plot relating dii to observation (without normalization)

                    subprocess.call(
                        "Rscript %s/r.avaflow.roc.R %s %s %s %s %s --slave --quiet"
                        % (scriptpath, temppath, pf, "4", "2", "2"),
                        shell=True,
                    )
                    # ROC plot relating dii to observation (with normalization)

                # Producing graphics summarizing the evaluation results

                if eflag and mflag and impdef == "1":

                    if sampling == "0" and len(ipar) == 2:

                        subprocess.call(
                            "Rscript %s/r.avaflow.multval.R %s %s %s %s %s --slave --quiet"
                            % (scriptpath, pf, str(ipar[0]), str(ipar[1]), "i", "0"),
                            shell=True,
                        )

                if eflag and mflag and depdef == "1":

                    if sampling == "0" and len(ipar) == 2:

                        subprocess.call(
                            "Rscript %s/r.avaflow.multval.R %s %s %s %s %s --slave --quiet"
                            % (scriptpath, pf, str(ipar[0]), str(ipar[1]), "d", "0"),
                            shell=True,
                        )

                if ctrlpoints and reftime:

                    if sampling == "0" and len(ipar) == 2:

                        for ictrlpoint in range(1, len(ctrlpoints) / 2 + 1):

                            subprocess.call(
                                "Rscript %s/r.avaflow.multval.R %s %s %s %s %s --slave --quiet"
                                % (
                                    scriptpath,
                                    pf,
                                    str(ipar[0]),
                                    str(ipar[1]),
                                    "t",
                                    str(ictrlpoint),
                                ),
                                shell=True,
                            )

                if mflag:

                    subprocess.call(
                        "Rscript %s/r.avaflow.merge.R %s --slave --quiet"
                        % (scriptpath, pf),
                        shell=True,
                    )  # merging parameter and evaluation tables

                    fin = open(
                        pf + "_results/" + pf + "_files/" + pf + "_allparams0.txt", "r"
                    )
                    fout = open(
                        pf + "_results/" + pf + "_files/" + pf + "_allparams.txt", "w"
                    )
                    for line in fin:
                        fout.write(line.replace("nrunp", "nrunp\tnrune"))
                    os.remove(
                        pf + "_results/" + pf + "_files/" + pf + "_allparams0.txt"
                    )

    # Cleaning file system and exiting

    os.system("rm -rf " + temppath)  # removing temporary directory
    if eflag and basechange == 0:
        os.system(
            "rm -rf " + ascpath + pf + "_basechange*"
        )  # removing obsolete ascii rasters

    grass.run_command(
        "g.remove", flags="f", type="rast", pattern="_*", quiet=True
    )  # removing temporary input and result raster maps
    grass.run_command(
        "g.remove", flags="f", type="vect", pattern="d_*", quiet=True
    )  # removing temporary vector maps

    if not kflag:
        grass.run_command(
            "g.remove", flags="f", type="rast", pattern=pf + "_*", quiet=True
        )  # removing result raster maps

    grass.run_command("g.region", flags="d")  # resetting default region

    if eflag and not mflag:

        if csuccess == 1:
            print()
            print(
                "Completed successfully in %.1f seconds (net computing time excluding visualization)."
                % comptime
            )
            print("Please find the collected results in the directory %s_results." % pf)
            print()

        elif csuccess == 0:
            print()
            print("The simulation was interrupted due to numerical failure.")
            print("Please find the collected results in the directory %s_results." % pf)
            print()

    elif eflag and mflag:

        print()
        print(
            "Completed in %.1f seconds (net computing time excluding visualization)."
            % comptime_batch
        )
        print("%s out of %s simulations were successful." % (str(nsuccess), str(nruns)))
        print("Please find the collected results in the directory %s_results." % pf)
        print()

    elif vflag and not eflag:

        print()
        print("Completed successfully.")
        print("Please find the collected results in the directory %s_results." % pf)
        print()

    sys.exit()  # exit


if __name__ == "__main__":
    options, flags = grass.parser()
    main()
