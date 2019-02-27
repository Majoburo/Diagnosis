import healpy as hp
import astropy.coordinates
import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
import matplotlib.pyplot as plt


def prob_observable(m, header, time, plot = False):
    """
    Determine the integrated probability contained in a gravitational-wave
    sky map that is observable with HET at a particular time. Needs hetpix.dat,
    pixels and lonlat of HET pupil, in the directory.
    
    """
    

    # Determine resolution of sky map
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
    observatory = astropy.coordinates.EarthLocation(
        lat=HET_loc[1]*u.deg, lon=HET_loc[0]*u.deg, height=HET_loc[2]*u.m)
    
    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg
    print(max(hetpupil[:,2]),min(hetpupil[:,2]))
    hetfullpix = hp.query_strip(nside, (90-max(hetpupil[:,2]))*np.pi/180, \
                            (90-min(hetpupil[:,2]))*np.pi/180)
    
    newpix = hp.ang2pix(nside,(90-hetpupil[:,2])*np.pi/180,
                        (hetpupil[:,1]+LST)*np.pi/180)

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
    sun_altaz = sun.transform_to(altaz)
    
    #delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    #midnight = astropy.time.Time('2012-7-13 00:00:00') - utcoffset
    #times24 = midnight + delta_midnight
    #frames24 = AltAz(obstime = times24, location=observatory)
    #sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frame24)
    #sunaltazs24.alt < -18*u.deg
    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, in this time and within 24 hours? 
    # Demand that it falls in the HETDEX pupil, that the sun is at least 18 
    # degrees below the horizon and that the airmass (secant of zenith angle 
    # approximation) is at most 2.5.
    mask_arraynow = np.zeros(len(m), dtype=int)
    mask_arraynow[newpix] = 1
    mask_arraynow *= (altaz.secz <= 2.5)&(sun_altaz.alt <= -18*u.deg)

    mask_arrayfull = np.zeros(len(m), dtype=int)
    mask_arrayfull[hetfullpix] = 1
    #mask_arrayfull *= (altaz.secz <= 2.5)#&(sun_altaz.alt <= -18*u.deg)

    
    prob = m[mask_arraynow > 0].sum()
    probfull = m[mask_arrayfull > 0].sum()
    
    if plot:

        msortedpix = np.flipud(np.argsort(m)) 
        cumsum = np.cumsum(m[msortedpix]) 
        cls = np.empty_like(m) 
        cls[msortedpix] = cumsum*100
        p90i =[]
        for i in range(len(m)):
            if cls[i] <= 90:
                p90i.append(i)
        p90i = np.array(p90i)
        
        #SUN CIRCLE OF 18 DEGREES
        radius = 18
        phi = Angle(sun.ra).radian
        theta = 0.5*np.pi-Angle(sun.dec).radian
        radius = np.deg2rad(radius)
        xyz = hp.ang2vec(theta, phi)
        ipix_sun = hp.query_disc(nside, xyz, radius)

        #Coloring the plot, order important here!
        m[altaz.secz > 2.5] = 0.5
        m[altaz.alt < 0] = 0.4
        m[newpix] = 0.8
        m[p90i] = 0.9
        m[ipix_sun] = 1
        hp.mollview(m, coord='C', cbar=False, max=1, title='HET',)
        hp.graticule(local=True)

        plt.savefig('MOLL_GWHET.pdf')
        #plt.show()

    # Done!
    return prob, probfull


