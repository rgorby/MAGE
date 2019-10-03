#!/usr/bin/env python

import argparse
from argparse import RawTextHelpFormatter
import datetime
import julian

fmt='%m/%d/%Y, %H:%M:%S'

if __name__ == "__main__":
    t0="1/1/2010, 00:00:00"

    MainS = """ Returns MJD (modified Julian date) from a given UT
                UT: UT string, MM/DD/YYYY, HH:MM:SS format
                ut2mjd.py "5/2/2010, 5:00:00"
    """

    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('UT',type=str,metavar="UT",default=t0,help='UT string to convert (use quotes) (default: %(default)s)')

    #Finalize parsing
    args = parser.parse_args()

    utStr = args.UT
    ut = datetime.datetime.strptime(utStr,fmt)
    mjd = julian.to_jd(ut,fmt='mjd')

    print("%15s (UT) => %f (MJD)"%(utStr,mjd))
