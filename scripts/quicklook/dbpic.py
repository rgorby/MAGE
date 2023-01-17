#!/usr/bin/env python


"""Plot the ground magnetic field perturbations from a magnetosphere run.

Plot the ground magnetic field perturbations from a magnetosphere run.

Author
------
Kareem Sorathia (kareem.sorathia@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import datetime
import os

# Import 3rd-party modules.
import cartopy.crs as ccrs
from cdasws import CdasWs
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

# Import project-specific modules.
import kaipy.cmaps.kaimaps as kmaps
import kaipy.gamera.deltabViz as dbViz
import kaipy.gamera.gampp as gampp
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as ktools
import kaipy.kaiViz as kv
import kaipy.satcomp.scutils as scutils


# Program constants and defaults

# Program description.
description = "Plot the ground magnetic field perturbations for a MAGE magnetosphere run."

# Default identifier for results to read.
default_runid = "msphere"

# Plot the last step by default.
default_step = -1

# Default vertical layer to plot.
default_k0 = 0

# Do not show anomalous currents by default.
default_Jr = False

# Default output filename.
default_output_filename = "qkdbpic.png"

# Size of figure in inches (width x height).
figSz = (12, 6)

# Valid spacecraft for position plotting.
valid_spacecraft = ('ACE')


def create_command_line_parser():
    """Create the command-line argument parser.
    
    Create the parser for command-line arguments.

    Parameters
    ----------
    None

    Returns
    -------
    parser : argparse.ArgumentParser
        Command-line argument parser for this script.
    """
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-d", type=str, metavar="directory", default=os.getcwd(),
        help="Directory containing data to read (default: %(default)s)"
    )
    parser.add_argument(
        "-n", type=int, metavar="step", default=default_step,
        help="Time slice to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-Jr", action="store_true", default=default_Jr,
        help="Show radial component of anomalous current (default: %(default)s)."
    )
    parser.add_argument(
        '-k0', type=int, metavar="layer", default=default_k0,
        help="Vertical layer to plot (default: %(default)s)")
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Name of spacecraft to plot (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    return parser


def datetime_to_float(d):
    epoch = datetime.datetime.utcfromtimestamp(0)
    total_seconds =  (d - epoch).total_seconds()
    # total_seconds will be in decimals (millisecond precision)
    return total_seconds


def float_to_datetime(fl):
    return datetime.datetime.fromtimestamp(fl)


def fetch_satellite_position(spacecraft, when):
    """Fetch the position of a satellite at a specified time.

    Fetch the position of a satellite at a specified time.
    Data is fetched from CDAWeb

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    when : datetime.datetime
        datetime for position fetch.

    Returns
    -------
    sc_lon, sc_lat : float
        Geographic longitude and latitude of spacecraft (degrees)
    """
    # <TEST>
    sc_lon, sc_lat = 90.0, 45.0
    # </TEST>

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Fetch the satellite position from CDAWeb.
    # Use the specified time as the start time, and nudge it by adding 1 hour
    # to get the end time.
    t0 = when.strftime("%Y-%m-%dT%H:%M:%SZ")
    t_end = when + datetime.timedelta(0, 3600)
    t1 = t_end.strftime("%Y-%m-%dT%H:%M:%SZ")
    status, data = cdas.get_data("GE_K0_MGF", "POS", t0, t1)

    # The ephemeris is in GSE X, Y, Z (km), so convert to geographic
    # longitude and latitude.
    xyz = data["POS"][0]
    scpos = Coords(xyz, "GSE", "car", use_irbem=False)
    scpos.ticks = Ticktock(t0)
    smpos = scpos.convert('GEO', 'sph')
    sc_lon = smpos.data[0][2]
    sc_lat = smpos.data[0][1]

    # Return the spacecraft position.
    return sc_lon, sc_lat


if __name__ == "__main__":
    """Plot the ground magnetic field perturbations."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    fdir = args.d
    runid = args.id
    nStp = args.n
    k0 = args.k0
    doJr = args.Jr
    spacecraft = args.spacecraft
    verbose = args.verbose
    if debug:
        print("args = %s" % args)

    # Fetch constants.
    bLin = dbViz.dbLin
    bMag = dbViz.dbMag
    jMag = dbViz.jMag

    # Compute the name of the file containing the ground magnetic field perturbations.
    ftag = runid + ".deltab"
    if debug:
        print("ftag = %s" % ftag)

    # Read the ground magnetic field perturbations.
    fname = os.path.join(fdir, ftag + ".h5")
    if debug:
        print("fname = %s" % fname)
    dbdata = gampp.GameraPipe(fdir, ftag)
    if debug:
        print("dbdata = %s" % dbdata)
    print("---")

    # Get the ID of the coordinate system, and the Earth radius.
    CoordID, Re = dbViz.GetCoords(fname)
    print("Found %s coordinate data ..." % CoordID)
    if debug:
        print("CoordID = %s" % CoordID)
        print("Re = %s" % Re)

    # If the last simulation step was requested, get the step number.
    if nStp < 0:
        nStp = dbdata.sFin
        print("Using Step %d" % nStp)

    # Check the vertical level.
    Z0 = dbViz.CheckLevel(dbdata, k0, Re)
    if debug:
        print("Z0 = %s" % Z0)

    # If currents were requested, read them. Otherwise, read the ground
    # magnetic field perturbations.
    if (doJr):
        print("Reading Jr ...")
        Jr = dbdata.GetVar("dbJ", nStp, doVerb=False)[:, :, k0]
        Q = Jr
    else:
        dBn = dbdata.GetVar("dBn", nStp, doVerb=True)[:, :, k0]
        Q = dBn

    # Convert MJD to UT.
    MJD = kh5.tStep(fname, nStp, aID="MJD")
    if debug:
        print("MJD = %s" % MJD)
    utS = ktools.MJD2UT([MJD])
    if debug:
        print("utS = %s" % utS)
    utDT= utS[0]
    if debug:
        print("utDT = %s" % utDT)

    # Create the mapping grid.
    crs = ccrs.PlateCarree()
    if debug:
        print("ccrs = %s" % ccrs)
    LatI, LonI, LatC, LonC = dbViz.GenUniformLL(dbdata, k0)
    if debug:
        print("LatI = %s" % LatI)
        print("LonI = %s" % LonI)
        print("LatC = %s" % LatC)
        print("LonC = %s" % LonC)

    # Fetch the color map.
    cmap = kmaps.cmDiv
    if debug:
        print("cmap = %s" % cmap)

    # Determine color bar settings.
    if (doJr):
        vQ = kv.genNorm(jMag)
        cbStr = "Anomalous current"
    else:
        vQ = kv.genNorm(bMag, doSymLog=True, linP=bLin)
        cbStr = r"$\Delta B_N$ [nT]"
    if debug:
        print("vQ = %s" % vQ)
        print("cbStr = %s" % cbStr)

    # Create plot in memory.
    mpl.use("Agg")

    # Create the figure to hold the plot.
    fig = plt.figure(figsize=figSz)

    # Specify the grid for the subplots.
    gs = gridspec.GridSpec(3, 1, height_ratios=[20, 1.0, 1.0], hspace=0.025)

    # Create the subplots.
    AxM = fig.add_subplot(gs[0, 0], projection=crs)
    AxCB = fig.add_subplot(gs[-1, 0])

    # Make the plot.
    AxM.pcolormesh(LonI, LatI, Q, norm=vQ, cmap=cmap)

    # If requested, overlay the spacecraft position.
    if spacecraft:
        print("Overplotting position of %s." % spacecraft)
        # Fetch the spacecraft position from CDAWeb.
        sc_lon, sc_lat = fetch_satellite_position(spacecraft, utS[0])
        if debug:
            print("sc_lon, sc_lat = %s, %s" % (sc_lon, sc_lat))

        # Plot a labelled dot at the location of the spacecraft position.
        AxM.plot(sc_lon, sc_lat, 'o')
        lon_nudge = 2.0
        lat_nudge = 2.0
        AxM.text(sc_lon + lon_nudge, sc_lat + lat_nudge, spacecraft)

    # Make the colorbar.
    kv.genCB(AxCB, vQ, cbStr, cM=cmap)

    # Add labels and other decorations.
    tStr = dbViz.GenTStr(AxM, fname, nStp)
    if debug:
        print("tStr = %s" % tStr)
    dbViz.DecorateDBAxis(AxM, crs, utDT)

    # Save the figure.
    fOut = default_output_filename
    if debug:
        print("fOut = %s" % fOut)
    kv.savePic(fOut)
