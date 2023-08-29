"""Standard plots and tools for gamhelio results.

This module provides standard plots and tools for gamhelio results.

Author
------
Eric Winter (eric.winter62@gmail.com)
"""


# Import standard modules.

# Import supplemental modules.
import numpy as np

# Import project modules.
import kaipy.kaiViz as kv
from kaipy.kdefs import Tsolar


# Solar wind speed limits and color map.
VMax = 800.0   # km/s
VMin = 300.0   # km/s
MagVCM = "inferno"

# Normalized number density plot limits and color map.
DMax = 150.0   # #/cm**3
DMin = 2000.0  # #/cm**3
DCM = "copper_r"

# Normalized temperature plot limits and color map.
TMin = 0.2  # megakelvin (MK)
TMax = 2.0  # megakelvin (MK)
TCM = "copper"

# Radial magnetic field plot limits and color map.
BMax = 150.0   # nT
BMin = -150.0  # nT
BCM = "coolwarm"

# 1 AU number density plot limits and color map.
D0Max = 1.0   # #/cm**3
D0Min = 15.0  # #/cm**3
D0CM = "copper_r"

# 1 AU temperature plot limits.
T0Min = 0.01  # MK
T0Max = 0.25  # MK

# 1 AU radial magnetic field plot limits.
B0Min = -4.0
B0Max = 4.0

# Color to use for profile plots.
colorProf = "tab:orange"


def GetSizeBds(pic):
    """Return plot domain based on plot type.

    Return the limits of the plotting domain based on the plot type.

    Parameters
    ----------
    pic : str
        String identifying plot type.

    Returns
    -------
    xyBds : list of 4 float
        Plot boundaries in the form (xmin, xmax, ymin, ymax).

    Raises
    ------
    TypeError
        If an invalid plot type is specified.
    """
    if pic == "pic1" or pic == "pic2":
        # Inner heliosphere - (xmin, xmax, ymin, ymax), in units of Rsun.
        xyBds = [-220.0, 220.0, -220.0, 220.0]
    elif pic == "pic3":
        # Solar lon/lat plot - (lonmin, lonmax, latmin, latmax), in degrees.
        xyBds = [0.0, 360.0, -75.0, 75.0]
    elif pic == "pic4":
        # Global solar lon/lat plot - (lonmin, lonmax, latmin, latmax),
        # in degrees.
        xyBds = [0.0, 360.0, -90.0, 90.0]
    elif pic == "pic5":
        xyBds = []
    else:
        raise TypeError("No pic type specified.")
    return xyBds


def PlotEqMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot solar wind speed in the equatorial plane.

    Plot solar wind speed in the equatorial plane.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    MagV : np.ndarray of float, shape (Nr, Np)
        Equatorial solar wind speeds (km/s)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, "Speed [km/s]", cM=MagVCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the equatorial solar wind speed for the specified step.
    MagV = gsph.eqMagV(nStp)

    # Plot as a pcolormesh.
    Ax.pcolormesh(gsph.xxi, gsph.yyi, MagV, cmap=MagVCM, norm=vMagV)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel("X [R_S]")
        Ax.set_ylabel("Y [R_S]")

    # Return the equatorial solar wind speed.
    return MagV


def PlotEqD(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot scaled number density in the equatorial plane.

    Plot scaled number density in the equatorial plane.
    The values are scaled as n*(r/r0)^2.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    NormD : np.ndarray of float, shape (Nr, Np)
        Scaled number density (#/cm**3)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, r"Density n(r/r$_0)^2$ [cm$^{-3}$]", cM=DCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled equatorial solar wind number density
    # for the specified step.
    NormD = gsph.eqNormD(nStp)

    # Plot as a pcolormesh.
    Ax.pcolormesh(gsph.xxi, gsph.yyi, NormD, cmap=DCM, norm=vD)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Y [R_S]')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the scaled equatorial solar wind number density.
    return NormD


def PlotEqTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot scaled temperature in equatorial plane.

    Plot scaled temperature in equatorial plane.
    The values are scaled as T*(r/r0), in MK (megakelvin)

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Temp : np.ndarray of float, shape (Nr, Np)
        Scaled temperature (MK)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, r"Temperature T(r/r$_0$) [MK]", cM=TCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled equatorial temperature for the specified step.
    Temp = gsph.eqTemp(nStp)

    # Plot as a pcolormesh.
    Ax.pcolormesh(gsph.xxi, gsph.yyi, Temp, cmap=TCM, norm=vT)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Y [R_S]')

    # Return the scaled equatorial solar wind temperature.
    return Temp


def PlotEqBr(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot the scaled radial magnetic field in equatorial plane.

    Plot the scaled radial magnetic field in equatorial plane.
    The values are scaled as Br*(r/r0)**2 (nT).

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Br : np.ndarray of float, shape (Nr, Np)
        Radial magnetic field (nT)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, r'Radial MF B$_r$(r/r$_0)^2$ [nT]', cM=BCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled equatorial radial magnetic field for the
    # specified step.
    Br = gsph.eqNormBr(nStp)

    # Plot as a pcolormesh.
    Ax.pcolormesh(gsph.xxi, gsph.yyi, Br, cmap=BCM, norm=vB)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Y [R_S]')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the scaled equatorial radial magnetic field.
    return Br


def PlotMerMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot solar wind speed in the meridional plane.

    Plot solar wind speed in the meridional plane (y = 0) containing Earth.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Vr, Vl : np.ndarray of float, shape (???)
        Meridional solar wind speeds (km/s)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, "Speed [km/s]", cM=MagVCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the solar wind speed for the specified step.
    Vr, Vl = gsph.MerMagV(nStp)

    # Get the coordinates of the cell corners.
    xr, zr, xl, zl, r = gsph.MeridGridHalfs()

    # Plot the speeds on the right (r) and left (l) sides.
    Ax.pcolormesh(xr, zr, Vr, cmap=MagVCM, norm=vMagV)
    Ax.pcolormesh(xl, zl, Vl, cmap=MagVCM, norm=vMagV)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel("X [R_S]")
        Ax.set_ylabel("Z [R_S]")

    # Return the meridional solar wind speed.
    return Vr, Vl


def PlotMerDNorm(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot scaled number density in the meridional plane.

    Plot scaled number density in the meridional (y = 0) plane containing
    Earth.
    The values are scaled as n*(r/r0)^2.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Dr, Dl : np.ndarray of float, shape (???)
        Scaled number density (#/cm**3)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vD = kv.genNorm(DMin, DMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, r"Density n$(r/r_0)^2$ [cm$^{-3}$]", cM=DCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled meridional solar wind number density
    # for the specified step.
    Dr, Dl = gsph.MerDNrm(nStp)

    # Get the coordinates of the cell corners.
    xr, zr, xl, zl, r = gsph.MeridGridHalfs()

    # Plot the density on the right (r) and left (l) sides.
    Ax.pcolormesh(xr, zr, Dr, cmap=DCM, norm=vD, shading='auto')
    Ax.pcolormesh(xl, zl, Dl, cmap=DCM, norm=vD, shading='auto')

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Z [R_S]')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the scaled meridional solar wind number density.
    return Dr, Dl


def PlotMerTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot scaled temperature in meridional plane.

    Plot scaled temperature in meridional (y = 0) plane, containing Earth.
    The values are scaled as T*(r/r0), in MK (megakelvin)

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Tempr, Templ : np.ndarray of float, shape (???)
        Scaled temperature (MK)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vT = kv.genNorm(TMin, TMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, r'Temperature T(r/r$_0$) [MK]', cM=TCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled meridional temperature for the specified step.
    Tempr, Templ = gsph.MerTemp(nStp)

    # Get the coordinates of the cell corners.
    xr, zr, xl, zl, r = gsph.MeridGridHalfs()

    # Plot the temperature on the right (r) and left (l) sides.
    Ax.pcolormesh(xr, zr, Tempr, cmap=TCM, norm=vT)
    Ax.pcolormesh(xl, zl, Templ, cmap=TCM, norm=vT)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Z [R_S]')

    # Return the scaled meridional solar wind temperature.
    return Tempr, Templ


# #Plot normalized Br Br(r/r0)^2 in meridional plane Y=0
def PlotMerBrNorm(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot the scaled radial magnetic field in meridional plane.

    Plot the scaled radial magnetic field in meridional plane.
    The values are scaled as Br*(r/r0)**2 (nT).

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Br_r, Br_l : np.ndarray of float, shape (???)
        Radial magnetic field (nT)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, r'Radial MF B$_r$(r/r$_0)^2$ [nT]', cM=BCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Compute the scaled meridional radial magnetic field for the
    # specified step.
    Br_r, Br_l = gsph.MerBrNrm(nStp)

    # Get the coordinates of the cell corners.
    xr, zr, xl, zl, r = gsph.MeridGridHalfs()

    # Plot the radial magnetic field on the right (r) and left (l) sides.
    Ax.pcolormesh(xr, zr, Br_r, cmap=BCM, norm=vB, shading='auto')
    Ax.pcolormesh(xl, zl, Br_l, cmap=BCM, norm=vB, shading='auto')

    # Cell-centered coordinates first.
    xr_c = 0.25*(xr[:-1, :-1] + xr[:-1, 1:] + xr[1:, :-1] + xr[1:, 1:])
    zr_c = 0.25*(zr[:-1, :-1] + zr[:-1, 1:] + zr[1:, :-1] + zr[1:, 1:])
    xl_c = 0.25*(xl[:-1, :-1] + xl[:-1, 1:] + xl[1:, :-1] + xl[1:, 1:])
    zl_c = 0.25*(zl[:-1, :-1] + zl[:-1, 1:] + zl[1:, :-1] + zl[1:, 1:])

    # Plot as a contour map, with a black line at the 0 level (this is the
    # location of the heliospheric current sheet).
    Ax.contour(xr_c, zr_c, Br_r, [0.], colors='black')
    Ax.contour(xl_c, zl_c, Br_l, [0.], colors='black')

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('X [R_S]')
        Ax.set_ylabel('Z [R_S]')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the scaled meridional radial magnetic field.
    return Br_r, Br_l


# #Plot Speed at 1 AU
def PlotiSlMagV(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot solar wind speed at 1 AU.

    Plot solar wind speed at 1 AU.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    V : np.ndarray of float, shape (??)
        Solar wind speeds (km/s)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vMagV = kv.genNorm(VMin, VMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vMagV, "Speed [km/s]", cM=MagVCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Extract the solar wind speeds at 1 AU.
    V = gsph.iSliceMagV(nStp)

    # Get a solar lat/lon grid.
    lat, lon = gsph.iSliceGrid()

    # Plot as a pcolormesh.
    Ax.pcolormesh(lon, lat, V, cmap=MagVCM, norm=vMagV)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')

    # Return the solar wind speed.
    return V


def PlotiSlD(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot number density at 1 AU.

    Plot number density at 1 AU.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    D : np.ndarray of float, shape (???)
        Number density (#/cm**3)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vD = kv.genNorm(D0Min, D0Max, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vD, "Number density [cm-3]", cM=D0CM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the number density at 1 AU.
    D = gsph.iSliceD(nStp)

    # Get a solar lat/lon grid.
    lat, lon = gsph.iSliceGrid()

    # Plot as a pcolormesh.
    Ax.pcolormesh(lon, lat, D, cmap=D0CM, norm=vD)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')

    # Return the number density.
    return D


def PlotiSlTemp(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot temperature at 1 AU.

    Plot temperature at 1 AU.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Temp : np.ndarray of float, shape (???)
        Temperature (MK)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vT = kv.genNorm(T0Min, T0Max, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vT, "Temperature [MK]", cM=TCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the temperature at 1 AU.
    Temp = gsph.iSliceT(nStp)

    # Get a solar lat/lon grid.
    lat, lon = gsph.iSliceGrid()

    # Plot as a pcolormesh.
    Ax.pcolormesh(lon, lat, Temp, cmap=TCM, norm=vT)

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')

    # Return the temperature.
    return Temp


# #Plot Br and current sheet (Br=0) at 1 AU
def PlotiSlBr(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot radial magnetic field at 1 AU.

    Plot radial magnetic field at 1 AU.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Br : np.ndarray of float, shape (???)
        Radial magnetic field (nT)

    Raises
    ------
    None
    """
    # Create a Normalize object for matplotlib.
    vB = kv.genNorm(B0Min, B0Max, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, "Radial magnetic field [nT]", cM=BCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the radial magnetic field at 1 AU.
    Br = gsph.iSliceBr(nStp)

    # Get a solar lat/lon grid.
    lat, lon = gsph.iSliceGrid()

    # Compute cell-centered lon lat coordinates (for contour overlay).
    lon_c = 0.25*(lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
    lat_c = 0.25*(lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])

    # Plot as a pcolormesh, with a contour line at Br = 0 (current sheet).
    Ax.pcolormesh(lon, lat, Br, cmap=BCM, norm=vB)
    Ax.contour(lon_c, lat_c, Br, [0.], colors='black')

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')
        # for pic4
        Ax.set_aspect('equal')

    # Return the radial magnetic field
    return Br


def PlotiSlBrRotatingFrame(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot the radial magnetic field at the inner grid boundary.

    Plot the radial magnetic field at the inner grid boundary.
    The values are scaled as Br*(r/r0)**2 (nT).

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    Br : np.ndarray of float, shape (Nr, Np)
        Radial magnetic field (nT)

    Raises
    ------
    None
    """
    # Set magnetic field limits
    BMin = -5.0  # nT
    BMax = 5.0   # nT

    # Create a Normalize object for matplotlib.
    vB = kv.genNorm(BMin, BMax, doLog=False, midP=None)

    # Add the colorbar (optional).
    if AxCB:
        AxCB.clear()
        kv.genCB(AxCB, vB, "Radial magnetic field [nT]", cM=BCM, Ntk=7)

    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the radial magnetic field at the inner grid boundary.
    Br = gsph.iSliceBrBound(nStp)

    # Get a solar lat/lon grid.
    lat, lon = gsph.iSliceGrid()

    # Transform the data into a rotating frame.
    jd0 = gsph.MJDs.min()   # Julian date of the initial map
    jd_c = gsph.MJDs[nStp]  # Julian date of the current step

    # Compute the elapsed days.
    time_days = jd_c - jd0

    # Compute the angular speed.
    omega = 2*180.0/Tsolar

    # Compute cell-centered lon lat coordinates for contour.
    lon_c = 0.25*(lon[:-1, :-1] + lon[:-1, 1:] + lon[1:, :-1] + lon[1:, 1:])
    lat_c = 0.25*(lat[:-1, :-1] + lat[:-1, 1:] + lat[1:, :-1] + lat[1:, 1:])

    # Compute the starting and current angles.
    phi = lon_c[0, :]
    phi_prime = (phi - omega*time_days) % (2*180.)

    if np.where(np.ediff1d(phi_prime) < 0)[0].size != 0:
        # For the first map size =0, for other maps size=1
        ind0 = np.where(np.ediff1d(phi_prime) < 0)[0][0] + 1
    else:
        # First map
        ind0 = 0

    # Rotate the data by the computed angle.
    Br = np.roll(Br, -ind0, axis=1)

    # Plot as a pcolormesh, with a contour line for the current sheet.
    Ax.pcolormesh(lon, lat, Br, cmap=BCM, norm=vB)
    Ax.contour(lon_c, lat_c, Br, [0.], colors='black')

    # Set the plot bounds.
    kv.SetAx(xyBds, Ax)

    # Add axis labels and tick marks (optional).
    if doDeco:
        Ax.set_xlabel('Longitude')
        Ax.set_ylabel('Latitude')
        Ax.yaxis.tick_right()
        Ax.yaxis.set_label_position('right')
        # For pic4
        Ax.set_aspect('equal')

    # Return the radial magnetic field.
    return Br


def PlotDensityProf(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot a radial number density profile.

    Plot a radial number density profile.

    Nr = # radius grid cells in gamera result file
    Np = # phi grid cells in gamera result file

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    D : np.ndarray of float, shape (???)
        Number density (#/cm**3)

    Raises
    ------
    None
    """
    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the radial density profile.
    D = gsph.RadProfDen(nStp)

    # Fetch the radial grid.
    rad = gsph.RadialProfileGrid()

    # Plot as a line plot.
    Ax.plot(rad, D, colorProf)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('Density [cm-3]')
        Ax.set_ylim(250.0, 450.0)
        Ax.set_xlim(20.0, 220.0)

    # Return the number density.
    return D


def PlotSpeedProf(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot a radial speed profile.

    Plot a radial speed profile.

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    V : np.ndarray of float, shape (???)
        Radial speed (km/s)

    Raises
    ------
    None
    """
    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the radial speed profile.
    V = gsph.RadProfSpeed(nStp)

    # Fetch the radius grid.
    rad = gsph.RadialProfileGrid()

    # Plot as a line plot.
    Ax.plot(rad, V, colorProf)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('Speed [km/s]')
        # Ax.set_ylim(600.0, 750.0)
        Ax.set_xlim(20.0, 220.0)

    # Return the speed profile.
    return V


def PlotFluxProf(
        gsph, nStp, xyBds, Ax, AxCB=None, doClear=True, doDeco=True
):
    """Plot a radial flux profile.

    Plot a radial flux profile.

    Parameters
    ----------
    gsph : GameraPipe
        Pipe interface to simulation results.
    nStp : int
        Step number of simulation step to plot
    xyBds : list of 4 float
        List of plot limits [xmin, xmax, ymin, ymax]
    Ax : Axes
        Axes object to use for the plot
    AxCB : Axes, default None
        Axes object to use for Colorbar.
    doClear : bool, default True
        True to clear the existing plot before plotting,
    doDeco : bool, default True
        True to add axis labels.

    Returns
    -------
    F : np.ndarray of float, shape (???)
        Flux (???)

    Raises
    ------
    None
    """
    # Clear the plot (optional).
    if doClear:
        Ax.clear()

    # Fetch the radial flux profile.
    F = gsph.RadProfFlux(nStp)

    # Fetch the radial grid.
    rad = gsph.RadialProfileGrid()

    # Plot as a line plot.
    Ax.plot(rad, F, colorProf)

    # Add axis labels (optional).
    if doDeco:
        Ax.set_xlabel('Radial distance [R_sun]')
        Ax.set_ylabel('RhoVr^2')
        Ax.set_ylim(180000.0, 280000.0)
        Ax.set_xlim(20.0, 220.0)

    # Return the flux.
    return F


if __name__ == "__main__":
    """Module self-test code"""
