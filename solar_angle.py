import numpy as np
import argparse
#from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from astropy.coordinates import get_sun
from astropy.time import Time
#from astropy.coordinates import Angle

parser = argparse.ArgumentParser(description="Compute Solar angle")
parser.add_argument('-c', '--ra_dec', type=str, nargs=2,
                    help="RA Dec of the pulsar (format, e.g.: 10h22m57.9992s +10d01m52.78s)")
parser.add_argument('-m', '--mjd', type=float, nargs=1,
                    help="MJD at which computing the Solar angle")

args = parser.parse_args()

if args.ra_dec:
    ra = args.ra_dec[0]
    dec = args.ra_dec[1]
if args.mjd:
    mjd = args.mjd



#coo_sun = SkyCoord('09h30m23s', '+14d45m09s', frame=FK5)
coo_psr = SkyCoord(ra,dec, frame=FK5)

t = Time(mjd, format='mjd')
sunpos = get_sun(t)
sep = sunpos.separation(coo_psr)#.to(u.degree)
#a = Angle(sep, u.degree)

print "Solar angle is",sep.degree[0],"degrees"

