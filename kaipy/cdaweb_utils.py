"""Simple utilitiles for fetching data from CDAWeb.

This module provides utility functions which perform simple queries to CDAWeb.

This module currently support all of the spacecraft that are described in
the YAML file used by the satellite comparison script msphSatComp.py.

Author
------
Eric Winter (eric.winter62@gmail.com)
"""


# Import standard modules.
import datetime

# Import supplemental modules.
from cdasws import CdasWs
import numpy as np
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

# Import project modules.
import kaipy.satcomp.scutils as scutils


# Program constants

# Format string for CDAWeb datetime strings.
CDAWEB_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%SZ"


def fetch_satellite_geographic_position(spacecraft, when):
    """Fetch the geographic position of a satellite at a specified time.

    Fetch the geographic position of a satellite at a specified time.
    Data is fetched from CDAWeb. The first returned position is assumed to
    correspond to the requested value of "when".

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    when : datetime.datetime
        datetime for position fetch.

    Returns
    -------
    sc_rad, sc_lat, sc_lon : float
        Geographic longitude and latitude of spacecraft (degrees) and
        radius (km).
    """
    # Read the CDAWeb spacecraft database.
    sc_info = scutils.getScIds()

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Format the start and end time strings.
    t0 = when.strftime(CDAWEB_DATETIME_FORMAT)
    # <HACK>
    # Use the specified time as the start time, and nudge it by adding 1 minute
    # (in seconds) to get the end time.
    one_minute = 60
    t_end = when + datetime.timedelta(0, one_minute)
    # </HACK>
    t1 = t_end.strftime(CDAWEB_DATETIME_FORMAT)

    # Fetch the satellite position from CDAWeb.
    dataset = sc_info[spacecraft]['Ephem']['Id']
    variable = sc_info[spacecraft]['Ephem']['Data']
    status, data = cdas.get_data(dataset, variable, t0, t1)

    # Return if no data found.
    if data is None:
        return None, None, None

    # The ephemeris is in Cartesian coordinates.
    # Assume the first/only position returned is closest position at t0.
    # GSM coordinates have r = SQRT(x**2 + y**2 + z**2) as a 4th element,
    # and we don't need that, so only copy the first 3 elements of any
    # position.
    if data[sc_info[spacecraft]['Ephem']['Data']].shape == (3,):
        # Only 1 position was returned, so copy its first 3 elements.
        xyz = data[sc_info[spacecraft]['Ephem']['Data']][:3]
    else:
        # More than 1 position was returned, so copy the first 3 elements of
        # the first position.
        xyz = data[sc_info[spacecraft]['Ephem']['Data']][0, :3]  # >1 position returned

    # Create a Coords object for the Cartesian position in the coordinate system used by
    # the spacecraft, at the start time.
    scpos = Coords(
        xyz,
        sc_info[spacecraft]['Ephem']['CoordSys'],
        'car', use_irbem=False
    )
    scpos.ticks = Ticktock(t0)

    # Convert the spacecraft coordinates to geographic spherical coordinates.
    smpos = scpos.convert('GEO', 'sph')
    sc_rad = smpos.data[0][0]
    sc_lat = smpos.data[0][1]
    sc_lon = smpos.data[0][2]

    # Return the spacecraft radius, longitude and latitude.
    return sc_rad, sc_lat, sc_lon


def fetch_satellite_SM_position(spacecraft, when):
    """Fetch the Solar Magnetic position of a satellite at a specified time.

    Fetch the Solar Magnetic position of a satellite at a specified time, in
    SM cartesian coordinates. Data is fetched from CDAWeb. The first returned
    position is assumed to correspond to the requested value of "when".

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    when : datetime.datetime
        datetime for position fetch.

    Returns
    -------
    sc_x, sc_y, sc_z : float
        Cartesian SM coordinates of spacecraft (units of kilometers).
    """
    # Read the CDAWeb spacecraft database.
    sc_info = scutils.getScIds()

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Format the start and end time strings.
    t0 = when.strftime(CDAWEB_DATETIME_FORMAT)
    # <HACK>
    # Use the specified time as the start time, and nudge it by adding 1 minute
    # (in seconds) to get the end time.
    one_minute = 60
    t_end = when + datetime.timedelta(0, one_minute)
    # </HACK>
    t1 = t_end.strftime(CDAWEB_DATETIME_FORMAT)

    # Fetch the satellite position from CDAWeb.
    dataset = sc_info[spacecraft]['Ephem']['Id']
    variable = sc_info[spacecraft]['Ephem']['Data']
    status, data = cdas.get_data(dataset, variable, t0, t1)

    # Return if no data found.
    if data is None:
        return None, None, None

    # The ephemeris is in Cartesian coordinates.
    # Assume the first/only position returned is closest position at t0.
    # GSM coordinates have r = SQRT(x**2 + y**2 + z**2) as a 4th element,
    # and we don't need that, so only copy the first 3 elements of any
    # position.
    if data[sc_info[spacecraft]['Ephem']['Data']].shape == (3,):
        # Only 1 position was returned, so copy its first 3 elements.
        xyz = data[sc_info[spacecraft]['Ephem']['Data']][:3]
    else:
        # More than 1 position was returned, so copy the first 3 elements of
        # the first position.
        xyz = data[sc_info[spacecraft]['Ephem']['Data']][0, :3]  # >1 position returned

    # Create a Coords object for the Cartesian position in the coordinate system used by
    # the spacecraft, at the start time.
    scpos = Coords(
        xyz,
        sc_info[spacecraft]['Ephem']['CoordSys'],
        'car', use_irbem=False
    )
    scpos.ticks = Ticktock(t0)

    # Convert the spacecraft coordinates to SM cartesian coordinates.
    smpos = scpos.convert('SM', 'car')
    sc_x = smpos.data[0][0]
    sc_y = smpos.data[0][1]
    sc_z = smpos.data[0][2]

    # Return the spacecraft cartesian SM coordinates.
    return sc_x, sc_y, sc_z


def fetch_spacecraft_SM_trajectory(spacecraft, t_start, t_end):
    """Fetch the Solar Magnetic trajectory of a spacecraft in a time range.

    Fetch the trajectory of a spacecraft in a time range, in Solar Magnetic
    Cartesian coordinates. Data is fetched from CDAWeb.

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    t_start, t_end : datetime.datetime
        Start and end datetime for trajectory fetch.

    Returns
    -------
    XXX : XXX
        XXX
    """
    # Read the CDAWeb spacecraft database.
    sc_info = scutils.getScIds()

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Format the start and end time strings.
    t0 = t_start.strftime(CDAWEB_DATETIME_FORMAT)
    t1 = t_end.strftime(CDAWEB_DATETIME_FORMAT)

    # Fetch the satellite position from CDAWeb.
    dataset = sc_info[spacecraft]['Ephem']['Id']
    variable = sc_info[spacecraft]['Ephem']['Data']
    coord_sys = sc_info[spacecraft]['Ephem']['CoordSys']
    status, data = cdas.get_data(dataset, variable, t0, t1)

    # Return if no data found.
    if data is None:
        return None, None, None

    # Extract the trajectory times as datetime objects.
    # <HACK>
    # Add code to handle this in the YAML file.
    if 'Epoch_state' in data.keys():
        t_trajectory = data['Epoch_state']
    elif 'Epoch' in data.keys():
        t_trajectory = data['Epoch']
    # </HACK>

    # The ephemeris is in Cartesian coordinates in an arbitrary frame.
    # GSM coordinates have r = SQRT(x**2 + y**2 + z**2) as a 4th element,
    # and we don't need that, so only copy the first 3 elements of any
    # position.
    if len(data[variable].shape) == 1:
        # Only 1 position was returned, so copy its first 3 elements
        # into a 2-D array.
        xyz_trajectory = np.array([data[variable][:3]])
    else:
        # More than 1 position was returned, so copy the first 3 elements of
        # each position.
        xyz_trajectory = data[variable][:, :3]  # >1 position

    # Create a Coords object for the Cartesian trajectory.
    trajectory = Coords(xyz_trajectory, coord_sys, 'car', use_irbem=False)
    trajectory.ticks = Ticktock(t_trajectory)

    # Convert the spacecraft coordinates to SM cartesian coordinates.
    trajectory_sm = trajectory.convert('SM', 'car')

    # Return the spacecraft cartesian SM coordinates.
    return trajectory_sm.x, trajectory_sm.y, trajectory_sm.z


def fetch_satellite_magnetic_northern_footprint_position(spacecraft, when):
    """Fetch the position of the northern magnetic footprint.

    Fetch the positions of the northern magnetic footprint for a spacecraft
    at a specified time. Data is fetched from CDAWeb. The first returned
    positions are assumed to correspond to the requested value of "when".

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    when : datetime.datetime
        datetime for position fetch.

    Returns
    -------
    fp_lon, fp_lat : float
        Geographic longitude and latitude (degrees) of northern magnetic
        footprint.
    """
    # Initialize the footprint position.
    fp_lat = None
    fp_lon = None

    # Read the CDAWeb spacecraft database.
    sc_info = scutils.getScIds()

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Format the start and end time strings.
    t0 = when.strftime(CDAWEB_DATETIME_FORMAT)
    # <HACK>
    # Use the specified time as the start time, and nudge it by adding 1 minute
    # (in seconds) to get the end time.
    one_minute = 60
    t_end = when + datetime.timedelta(0, 3600)
    # </HACK>
    t1 = t_end.strftime(CDAWEB_DATETIME_FORMAT)

    # Fetch the footprint position from CDAWeb.
    dataset = sc_info[spacecraft]['MagneticFootprintNorth']['Id']
    variable = sc_info[spacecraft]['MagneticFootprintNorth']['Data']
    coordsys = sc_info[spacecraft]['MagneticFootprintNorth']['CoordSys']
    status, data = cdas.get_data(dataset, variable, t0, t1)

    # Return if no data found.
    if data is None:
        return fp_lat, fp_lon

    # Convert the position if needed.
    if coordsys == 'GEO':
        if type(variable) is list:
            # <HACK>
            # Assume index 0 is latitude, index 1 is longitude.
            fp_lat = data[variable[0]][0]
            fp_lon = data[variable[1]][0]
            # </HACK>
        else:
            raise TypeError("Unexpected variable type!")
    elif coordsys == 'GSM':
        # The position is in cartesian GSM coordinates (kilometers).
        if data[variable].shape == (3,):
            # Only 1 position was returned, so copy its first 3 elements.
            xyz = data[variable][:3]
        else:
            # More than 1 position was returned, so copy the first 3 elements of
            # the first position.
            xyz = data[variable][0, :3]  # >1 position returned

        # Create a Coords object for the Cartesian position in the coordinate system used by
        # the spacecraft, at the start time.
        fp_pos_orig = Coords(xyz, coordsys, 'car', use_irbem=False)
        fp_pos_orig.ticks = Ticktock(t0)

        # Convert the footprint coordinates to geographic spherical coordinates.
        fp_pos_sph = fp_pos_orig.convert('GEO', 'sph')
        fp_rad = fp_pos_sph.data[0][0]
        fp_lat = fp_pos_sph.data[0][1]
        fp_lon = fp_pos_sph.data[0][2]
    else:
        raise TypeError("Unexpected coordinate system: %s" % coordsys)

    # Canonicalize longitude to [-180, +180] range.
    if fp_lon > 180:
        fp_lon = fp_lon - 360

    # Return the spacecraft longitude and latitude.
    return fp_lat, fp_lon


def fetch_satellite_magnetic_southern_footprint_position(spacecraft, when):
    """Fetch the position of the southern magnetic footprint.

    Fetch the positions of the southern magnetic footprint for a spacecraft
    at a specified time. Data is fetched from CDAWeb. The first returned
    positions are assumed to correspond to the requested value of "when".

    Parameters
    ----------
    spacecraft : str
        CDAWeb-compliant spacecraft ID.
    when : datetime.datetime
        datetime for position fetch.

    Returns
    -------
    fp_lon, fp_lat : float
        Geographic longitude and latitude (degrees) of southern magnetic
        footprint.
    """
    # Initialize the footprint position.
    fp_lat = None
    fp_lon = None

    # Read the CDAWeb spacecraft database.
    sc_info = scutils.getScIds()

    # Create the CDAWeb connection.
    cdas = CdasWs()

    # Format the start and end time strings.
    t0 = when.strftime(CDAWEB_DATETIME_FORMAT)
    # <HACK>
    # Use the specified time as the start time, and nudge it by adding 1 minute
    # (in seconds) to get the end time.
    one_minute = 60
    t_end = when + datetime.timedelta(0, 3600)
    # </HACK>
    t1 = t_end.strftime(CDAWEB_DATETIME_FORMAT)

    # Fetch the footprint position from CDAWeb.
    dataset = sc_info[spacecraft]['MagneticFootprintSouth']['Id']
    variable = sc_info[spacecraft]['MagneticFootprintSouth']['Data']
    coordsys = sc_info[spacecraft]['MagneticFootprintSouth']['CoordSys']
    status, data = cdas.get_data(dataset, variable, t0, t1)

    # Return if no data found.
    if data is None:
        return fp_lat, fp_lon

    # Convert the position if needed.
    if coordsys == 'GEO':
        if type(variable) is list:
            # <HACK>
            # Assume index 0 is latitude, index 1 is longitude.
            fp_lat = data[variable[0]][0]
            fp_lon = data[variable[1]][0]
            # </HACK>
        else:
            raise TypeError("Unexpected variable type!")
    elif coordsys == 'GSM':
        # The position is in cartesian GSM coordinates (kilometers).
        if data[variable].shape == (3,):
            # Only 1 position was returned, so copy its first 3 elements.
            xyz = data[variable][:3]
        else:
            # More than 1 position was returned, so copy the first 3 elements of
            # the first position.
            xyz = data[variable][0, :3]  # >1 position returned

        # Create a Coords object for the Cartesian position in the coordinate system used by
        # the spacecraft, at the start time.
        fp_pos_orig = Coords(xyz, coordsys, 'car', use_irbem=False)
        fp_pos_orig.ticks = Ticktock(t0)

        # Convert the footprint coordinates to geographic spherical coordinates.
        fp_pos_sph = fp_pos_orig.convert('GEO', 'sph')
        fp_rad = fp_pos_sph.data[0][0]
        fp_lat = fp_pos_sph.data[0][1]
        fp_lon = fp_pos_sph.data[0][2]
    else:
        raise TypeError("Unexpected coordinate system: %s" % coordsys)

    # Canonicalize longitude to [-180, +180] range.
    if fp_lon > 180:
        fp_lon = fp_lon - 360

    # Return the spacecraft longitude and latitude.
    return fp_lat, fp_lon
