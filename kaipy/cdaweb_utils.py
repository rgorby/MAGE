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
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

# Import project modules.
import kaipy.satcomp.scutils as scutils


# Program constants

# Format string for CDAWeb datetime strings.
CDAWEB_DATETIME_FORMAT = "%Y-%m-%dT%H:%M:%SZ"


def fetch_satellite_geographic_position(spacecraft, when):
    """Fetch the position of a satellite at a specified time.

    Fetch the position of a satellite at a specified time.
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
    sc_lon, sc_lat, sc_rad : float
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
    status, data = cdas.get_data(
        sc_info[spacecraft]['Ephem']['Id'],
        sc_info[spacecraft]['Ephem']['Data'],
        t0, t1
    )

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

    # Return the spacecraft longitude and latitude.
    return sc_rad, sc_lat, sc_lon
