#!/usr/bin/env python


"""Create north and south REMIX plots.

Create north and south REMIX plots.

Author
------
Kareem Sorathia (kareem.sorathia@jhuapl.edu)
Eric Winter (eric.winter@jhuapl.edu)
"""


# Import standard modules.
import argparse
import sys
import os

# Import supplemental modules.
from astropy.time import Time
# import numpy as np
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import h5py as h5
from alive_progress import alive_it

# Import project-specific modules.
import kaipy.cdaweb_utils as cdaweb_utils
import kaipy.kaiH5 as kaiH5
import kaipy.kaiTools as ktools
import kaipy.remix.remix as remix
import kaipy.kaiViz as kv
import kaipy.kdefs as kd


# Program constants and defaults

# Program description.
description = """Creates simple multi-panel REMIX figure for a GAMERA magnetosphere run.
Top Row - FAC (with potential contours overplotted), Pedersen and Hall Conductances
Bottom Row - Joule heating rate, particle energy and energy flux
"""

# Default identifier for results to read.
default_runid = "msphere"

# Plot the last step by default.
default_step = -1

# Color to use for magnetic footprint positions.
FOOTPRINT_COLOR = 'red'


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
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--debug", action="store_true", default=False,
        help="Print debugging output (default: %(default)s)."
    )
    parser.add_argument(
        "-id", type=str, metavar="runid", default=default_runid,
        help="Run ID of data (default: %(default)s)"
    )
    parser.add_argument(
        "-n", type=int, metavar="step", default=default_step,
        help="Time slice to plot (default: %(default)s)"
    )
    parser.add_argument(
        '-nflux', action='store_true', default=False,
        help="Show number flux instead of energy flux (default: %(default)s)"
    )
    parser.add_argument(
        '-print', action='store_true', default=False,
        help="Print list of all steps and time labels, then exit (default: %(default)s)"
    )
    parser.add_argument(
        "--spacecraft", type=str, metavar="spacecraft", default=None,
        help="Names of spacecraft to plot magnetic footprints, separated by commas (default: %(default)s)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Print verbose output (default: %(default)s)."
    )
    parser.add_argument(
        '-GTYPE', action='store_true', default=False,
        help="Show RCM grid type in the eflx plot (default: %(default)s)"
    )
    parser.add_argument(
        '-PP', action='store_true', default=False,
        help="Show plasmapause (10/cc) in the eflx/nflx plot (default: %(default)s)"
    )
    parser.add_argument(
        '--vid', action='store_true', default=False,
        help="Make a video and store in mixVid directory (default: %(default)s)"
    )
    return parser

def makePlot_single(remixFile, nStp, h):

    with h5.File(remixFile, 'r') as f5:
        foundT = f5['Step#'+str(nStp)].attrs['MJD']
    #foundT = T[idxStp]
    if debug:
        print("foundT = %s" % foundT)
    print('Found time:', Time(foundT, format='mjd').iso)
    utS = ktools.MJD2UT(foundT)
    if debug:
        print("utS = %s" % utS)


    # Read the data into the remix object.
    ion = remix.remix(remixFile, nStp)
    if debug:
        print("ion = %s" % ion)

    # Add a label.
    plt.figtext(
        0.5, 0.94, 'MIX (' + h + ')\n' + Time(foundT, format='mjd').iso,
        fontsize=12, multialignment='center', horizontalalignment='center'
    )

    # Initialize the remix object based on the current hemisphere.
    ion.init_vars(h)

    # Create the plot layout for the current hemisphere.
    gs = gridspec.GridSpec(
        2, 3, figure=fig, left=0.03, right=0.97, top=0.9, bottom=0.03
    )

    # Create the individual plots for the current hemisphere.
    axs = [None]*6
    axs[0] = ion.plot('current', gs=gs[0, 0])
    axs[1] = ion.plot('sigmap', gs=gs[0, 1])
    axs[2] = ion.plot('sigmah', gs=gs[0, 2])
    axs[3] = ion.plot('joule', gs=gs[1, 0])
    axs[4] = ion.plot('energy', gs=gs[1, 1])
    if do_nflux:
        axs[5] = ion.plot('flux', gs=gs[1, 2],doGTYPE=do_GTYPE,doPP=do_PP)
    else:
        axs[5] = ion.plot('eflux', gs=gs[1, 2],doGTYPE=do_GTYPE,doPP=do_PP)

    # If requested, plot the magnetic footprints for the specified
    # spacecraft.
    if spacecraft:
        for sc in spacecraft:
            if verbose:
                print("Overplotting %s magnetic footprint for %s." % (h, sc))

            # Fetch the footprint position for this hemisphere.
            if h.lower() == 'north':
                fp_lat, fp_lon = cdaweb_utils.fetch_satellite_magnetic_northern_footprint_position(
                    sc, utS
                )
            else:
                fp_lat, fp_lon = cdaweb_utils.fetch_satellite_magnetic_southern_footprint_position(
                    sc, utS
                )
            if debug:
                print("fp_lat, fp_lon = %s, %s" % (fp_lat, fp_lon))

            # Skip if no footprint found.
            if fp_lat is None:
                print("No %s footprint found for spacecraft %s." % (h, sc))
                continue

            # The footprint locations are in geographic (GEO) coordinates.
            # They must be converted to Solar Magnetic (SM) coordinates
            # for plotting.

            # Convert the footprint position to the coordinate system used
            # by these plots, which show contours at the surface of the
            # ionosphere, about 122 km above the surface ofn the Earth.
            # Note that this adjustment assumes the field lines impinging
            # on the magnetic footprint descend vertically at the
            # footprint point, which is not technically accurate.
            fp_lat_rad = np.radians(fp_lat)
            fp_lon_rad = np.radians(fp_lon)
            fp_x = np.cos(fp_lat_rad)*np.cos(fp_lon_rad)
            fp_y = np.cos(fp_lat_rad)*np.sin(fp_lon_rad)
            fp_theta = np.arctan2(fp_y, fp_x)  # [-pi, pi]
            fp_r = np.sqrt(fp_x**2 + fp_y**2)

            # Plot a labelled dot at the location of each footprint.
            # Skip if no footprint position found.
            for ax in axs:
                ax.plot(fp_theta, fp_r, 'o', c=FOOTPRINT_COLOR)
                theta_nudge = 0.0
                r_nudge = 0.0
                ax.text(fp_theta + theta_nudge, fp_r + r_nudge, sc)


def makeAndSaveFig(remixFile, nStp, h, dir, outFname):
    # Make sure figure is clean
    plt.clf()

    # Draw to figure
    makePlot_single(remixFile, nStp, h)

    # Save to file
    outPath = os.path.join(dir, outFname)
    kv.savePic(outPath, dpiQ=300)


if __name__ == "__main__":
    """Plot the ground magnetic field perturbations."""

    # Set up the command-line parser.
    parser = create_command_line_parser()

    # Parse the command-line arguments.
    args = parser.parse_args()
    debug = args.debug
    runid = args.id
    nStp = args.n
    do_nflux = args.nflux
    do_print = args.print
    spacecraft = args.spacecraft
    verbose = args.verbose
    do_GTYPE = args.GTYPE
    do_PP = args.PP
    do_vid = args.vid
    if debug:
        print("args = %s" % args)

    # Construct the name of the REMIX results file.
    remixFile = runid + '.mix.h5'
    if debug:
        print("remixFile = %s" % remixFile)

    # Split the original string into a list of spacecraft IDs.
    if spacecraft:
        spacecraft = spacecraft.split(',')
        if debug:
            print("spacecraft = %s" % spacecraft)

    # Enumerate the steps in the results file.
    nsteps, sIds = kaiH5.cntSteps(remixFile)
    sIds = sorted(sIds)
    if debug:
        print("nsteps = %s" % nsteps)
        print("sIds = %s" % sIds)

    # Check that the requested step exists.
    if nStp >= 0 and not nStp in sIds:  # ANY nStp<0 gets last step.
        raise TypeError(f"Step #{nStp} not found in {remixFile}!")

    # Get the times from the result file.
    if do_print:
        T = kaiH5.getTs(remixFile, sIds, aID='MJD')
        if debug:
            print("T = %s" % T)
        for i, tt in enumerate(T):
            print('Step#%06d: ' % sorted(sIds)[i], Time(tt, format='mjd').iso)
        sys.exit(0)

    
    # Create the plots in a memory buffer.
    mpl.use('Agg')

    # Set global plot font options.
    mpl.rc('mathtext', fontset='stixsans', default='regular')
    mpl.rc('font', size=10)

    # Init figure
    fig = plt.figure(figsize=(12, 7.5))

    if not do_vid:  # Then we are making a single image, keep original functionality
        # Find the time for the specified step.
        if nStp == -1:
            # Take the last step.
            nStp = sorted(sIds)[-1]
        if debug:
            print("nStp = %s" % nStp)
        makeAndSaveFig(remixFile, nStp, 'NORTH', '.', 'remix_n.png')
        makeAndSaveFig(remixFile, nStp, 'SOUTH', '.', 'remix_s.png')
    
    else:  # Then we make a video, i.e. series of images saved to mixVid
        outDir = 'mixVid'
        kaiH5.CheckDirOrMake(outDir)

        # How many 0's do we need for filenames?
        n_pad = int(np.log10((len(sIds)))) + 1

        # Do NORTH first
        for i, nStp in enumerate(alive_it(sIds,title="North".ljust(kd.barLab),length=kd.barLen,bar=kd.barDef)):
            
            filename = "{}.{:0>{n}d}.png".format("remix_n", i, n=n_pad)
            makeAndSaveFig(remixFile, nStp, 'NORTH', outDir, filename)

        # Do SOUTH second
        for i, nStp in enumerate(alive_it(sIds,title="South".ljust(kd.barLab),length=kd.barLen,bar=kd.barDef)):
            filename = "{}.{:0>{n}d}.png".format("remix_s", i, n=n_pad)
            makeAndSaveFig(remixFile, nStp, 'SOUTH', outDir, filename)