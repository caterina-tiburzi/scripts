import os, sys, argparse
import datetime, jdcal
import numpy as np
from astropy.time import Time

parser = argparse.ArgumentParser(description="")

parser.add_argument('-u', '--utc', type=str, nargs=1,
                    help="UTC in YYYY-MM-DDThh:mm:ss")

parser.add_argument('-p', '--psr', type=str, nargs=1,
                    help="Pulsar Jname")

parser.add_argument('-o', '--observatory', type=str, nargs=1,
                    help="Name of the observatory")

args = parser.parse_args()

UTC=args.utc[0]#"2016-06-03T16:09:00.0"#utc is a variant of ut, we can use it as ut (I think)
psr = args.psr[0]
obs = args.observatory[0]

pulsars = np.genfromtxt("radec.dat", names='ind, jpsr, rajd, decjd', dtype=None)
index_psr = list(pulsars['jpsr']).index(psr)
raj = pulsars['rajd'][index_psr]
decj = pulsars['decjd'][index_psr]


latsites = {"DE601": 50.5230, "DE602": 48.501455, "DE603": 50.980111, "DE604": 52.438, "DE605": 50.897, "DE609": 53.698454}
lonsites = {"DE601": 6.8837, "DE602": 11.287696, "DE603": 11.711167, "DE604": 13.016, "DE605": 6.424, "DE609": 9.969239}

lat = latsites[obs]
lon = lonsites[obs]

print "UTC:",UTC


mjd_j2000 = 51544.


yr=np.float((UTC.split("T")[0]).split("-")[0])
mo=np.float((UTC.split("T")[0]).split("-")[1])
da=np.float((UTC.split("T")[0]).split("-")[2])
hh=np.float((UTC.split("T")[1]).split(":")[0])
mm=np.float((UTC.split("T")[1]).split(":")[1])
ss=np.float((UTC.split("T")[1]).split(":")[2])

ut=hh+mm/60.+ss/3600.#ut in decimal hours

fmt = '%Y.%m.%d'
s = '%s.%s.%s'%((UTC.split("T")[0]).split("-")[0],(UTC.split("T")[0]).split("-")[1],(UTC.split("T")[0]).split("-")[2] )
dt = datetime.datetime.strptime(s, fmt)

tt = dt.timetuple()
nthdayyear = tt.tm_yday
mjd = sum(jdcal.gcal2jd(dt.year, dt.month, dt.day)) - 2400000.5

print "MJD:",mjd




#de609 = EarthLocation(lat='31d57.5m', lon='-111d35.8m', height=2096*u.m)
#c = SkyCoord(ra=raj*u.degree, dec=decj*u.degree, frame='icrs')

#from http://www.stargazing.net/kepler/altaz.html

A, B, C = 100.46, 0.985647, 15.

mjd = mjd + ut/24.
d = mjd - mjd_j2000

lst = A + B * d + lon + C*ut #lst in degrees
turns = np.floor(lst/360.)

if (turns != 0):
    lst = lst - turns * 360.#to bring back lst between 0 and 360 (why??)

print "LST:",lst

ha = lst - raj #angular distance btw the local meridian and the meridian passing the pulsar

if ha<0:
    ha = ha + 360.

alt = np.degrees(np.arcsin(np.sin(np.radians(decj))*np.sin(np.radians(lat))+np.cos(np.radians(decj))*np.cos(np.radians(lat))*np.cos(np.radians(ha))))
print "Elevation:", alt





#def date_to_nth_day(date, format='%Y%m%d'):
#    date = pd.to_datetime(date, format=format)
#    new_year_day = pd.Timestamp(year=date.year, month=1, day=1)
#    return (date - new_year_day).days + 1
