import healpy as hp
import astropy.coordinates
import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
maxhetdec = 74
minhetdec = -12
minhetdec_rad = (90-maxhetdec)*np.pi/180
maxhetdec_rad = (90-minhetdec)*np.pi/180
def prob_observable(m, header, time, plot = False):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable with HET at a particular time. Needs hetpix.dat,
    pixels and lonlat of HET pupil, in the directory.

    """

    # Determine resolution of sky map
    mplot = np.copy(m)
    npix = len(m)
    nside = hp.npix2nside(npix)
    # Get time now and Local Sidereal Time
    # time = astropy.time.Time.now()
    # Or at the time of the gravitational-wave event...
    # time = astropy.time.Time(header['MJD-OBS'], format='mjd')
    # Or at a particular time...
    # time = astropy.time.Time(time)

    # Geodetic coordinates of MacDonald Obs

    HET_loc = (-104.01472,30.6814,2025)
    hetpupil = np.loadtxt('hetpix.dat')
    hetfullpix = hp.query_strip(nside, minhetdec_rad, \
                            maxhetdec_rad)

    observatory = astropy.coordinates.EarthLocation(
        lat=HET_loc[1]*u.deg, lon=HET_loc[0]*u.deg, height=HET_loc[2]*u.m)

    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg
    HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
    HETtheta = (90-hetpupil[:,2])*np.pi/180
    newpix = hp.ang2pix(nside, HETtheta, HETphi)
    newpixp = newpix


    # Alt/az reference frame at the observatory, in this time
    frame = astropy.coordinates.AltAz(obstime=t, location=observatory)

    # Look up (celestial) spherical polar coordinates of HEALPix grid.
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    # Convert to RA, Dec.
    radecs = astropy.coordinates.SkyCoord(
        ra=phi*u.rad, dec=(0.5*np.pi - theta)*u.rad)
    # Transform grid to alt/az coordinates at observatory, in this time
    altaz = radecs.transform_to(frame)

    #Get RA,DEC of the sun in this time
    sun = astropy.coordinates.get_sun(time)
    # Where is the sun in the Texas sky, in this time?
    sun_altaz = sun.transform_to(frame)

    delta_time = np.linspace(0, 24, 1000)*u.hour
    times24 = t + delta_time
    frames24 = astropy.coordinates.AltAz(obstime = times24, location=observatory)
    sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frames24)
    timetilldark = 0*u.hour
    timetillbright = 0*u.hour

    nightstart = times24[sunaltazs24.alt<-18*u.deg][0]
    nighttimemask = np.array((sunaltazs24.alt<-18*u.deg))*1
    if (sun_altaz.alt > -18*u.deg):
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][1]
    else:
        nightend = times24[(np.roll(nighttimemask, 1) - nighttimemask) != 0][0]
    nightime = nightend - nightstart
    nightime.format = 'sec'
    #Moving to start of the night if in daytime
    if (sun_altaz.alt > -18*u.deg):
        timetilldark = (nightstart-t)
        timetilldark.format = 'sec'
        LST = nightstart.sidereal_time('mean').deg
        HETphi = ((hetpupil[:,1]+LST)%360)*np.pi/180
        newpix = hp.ang2pix(nside, HETtheta, HETphi)
    else:
        timetillbright = (nightend-t)
        timetillbright.format = 'sec'

    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, in this time and within 24 hours? 
    # Demand that it falls in the HETDEX pupil, that the sun is at least 18 
    # degrees below the horizon and that the airmass (secant of zenith angle 
    # approximation) is at most 2.5.

    msortedpix = np.flipud(np.argsort(m))
    cumsum = np.cumsum(m[msortedpix])
    cls = np.empty_like(m)
    cls[msortedpix] = cumsum*100
    p90i = np.where(cls <= 90)
    if plot:

        #SUN CIRCLE OF 18 DEGREES
        radius = 18
        phis = Angle(sun.ra).radian
        thetas = 0.5*np.pi-Angle(sun.dec).radian
        radius = np.deg2rad(radius)
        xyz = hp.ang2vec(thetas, phis)
        ipix_sun = hp.query_disc(nside, xyz, radius)

        #Coloring the plot, order important here!
        mplot[altaz.secz > 2.5] = 0.5
        mplot[altaz.alt < 0] = 0.4
        mplot[newpixp] = 0.8
        mplot[p90i] = 0.9
        mplot[ipix_sun] = 1
        hp.mollview(mplot, coord='C', cbar=False, max=1, title='HET NOW',)
        hp.graticule(local=True)

        plt.savefig('MOLL_GWHET_%s.pdf'%header['GraceID'])
        #plt.show()
    theta90, phi90 = hp.pix2ang(nside, p90i)
    #mask skymap pixels by hetdex accesible region
    theta90HETi = (theta90 > minhetdec_rad)*(theta90 < maxhetdec_rad)
    print(theta90HETi.sum(),theta90.min(),theta90.max())
    theta90HET = theta90[theta90HETi]
    phi90HET = phi90[theta90HETi]
    timetill90 = 0
    #if the region doesn't intersect HET now
    if len(np.intersect1d(p90i,newpix)) == 0:
        #if the region doesn't intersect HET at all
        if len(np.intersect1d(p90i,hetfullpix)) == 0:
            return 0 , 0 , -99, 0
        hetedge = np.loadtxt('hetedge.dat')
        hetedgef = lambda x: np.interp(x,hetedge[:,0],hetedge[:,1])
        y = theta90HET*180/np.pi #DEC SKYMAP
        x = (phi90HET*180/np.pi-LST)%360 #RA SKYMAP ZEROING HET PUPIL
        wsecs = np.min(x - hetedgef(y))*3600*12/180
        if wsecs < 0:
            #it's inside the pupil...
            hetedgef2 = lambda x: np.interp(x,hetedge[:,0],hetedge[:,2])
            y = (hetedgef2(y)+180)%360 -180
            x = (x+180)%360-180
            wsecs = np.min(x - y)*3600*12/180

        if timetilldark == 0:
            if wsecs > timetillbright.value:
                return 0 , 0 , -99, 0
        else:
            if wsecs > nightime.value:
                return 0 , 0 , -99, 0
        timetill90 = (wsecs+timetilldark.value)/3600
    elif timetilldark.value > 0:
        timetill90 = timetilldark.value/3600

    mask_arraynow = np.zeros(len(m), dtype=int)
    mask_arraynow[newpixp] = 1
    mask_arraynow *= (altaz.secz <= 2.5)&(sun_altaz.alt <= -18*u.deg)

    prob = m[mask_arraynow > 0].sum()
    probfull = m[np.intersect1d(p90i,hetfullpix)].sum()
    m[np.setdiff1d(np.arange(len(m)),np.intersect1d(p90i,hetfullpix),assume_unique=True)]=m.min()
    #hp.orthview(m)
    #plt.show()
    #plt.savefig('MOLL_GWHET_%s.pdf'%header['GraceID'])
    # Done!
    return prob, probfull, timetill90, m

def main():
        skymap, header = hp.read_map(sys.argv[1],
                                     h=True, verbose=False)
        header = {'GraceID': 'TEST'}
        time = astropy.time.Time.now()
        prob, probfull, timetill90, m = prob_observable(skymap, header, time, plot=True)
        print(timetill90)
        return timetill90
if __name__=="__main__":
    main()
