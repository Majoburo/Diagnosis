import healpy as hp
import astropy.coordinates
import astropy.time
import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def no_gaps(area1,area2):
    #IN RADIANS
    #ensure no coordinate gaps in source RAs
    shift = np.pi/6*(1+np.sqrt(5))/2 # use golden ratio so we're guaranteed to terminate eventually (unless area covers whole sky)
    i=0
    while (np.max(area1) - np.min(area1)) > 359*np.pi/180 or (np.max(area2) - np.min(area2)) > 359*np.pi/180:
        area1 = ((area1 + shift) % 2*np.pi - shift)
        area2 = ((area2 + shift) % 2*np.pi - shift)
        i+=1
        if i > 100:
            print("Couldn't calculate time until observation. Time will probably be negative. Check things manually...")
            return area1, area2
    return area1, area2

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
    hetfullpix = hp.query_strip(nside, (90-max(hetpupil[:,2]))*np.pi/180, \
                            (90-min(hetpupil[:,2]))*np.pi/180)

    observatory = astropy.coordinates.EarthLocation(
        lat=HET_loc[1]*u.deg, lon=HET_loc[0]*u.deg, height=HET_loc[2]*u.m)
    
    # Find pixels of HET pupil in this time
    t = astropy.time.Time(time,scale='utc',location=HET_loc)
    LST = t.sidereal_time('mean').deg

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
    #import pdb
    #pdb.set_trace()
    
    delta_time = np.linspace(0, 24, 1000)*u.hour
    #midnight = astropy.time.Time('2012-7-13 00:00:00') - utcoffset
    times24 = t + delta_time
    frames24 = astropy.coordinates.AltAz(obstime = times24, location=observatory)
    sunaltazs24 = astropy.coordinates.get_sun(times24).transform_to(frames24)
    timetilldark = 0*u.hour

    nightstart = times24[sunaltazs24.alt<-18*u.deg][0]
    nightend = times24[sunaltazs24.alt<-18*u.deg][-1]
    nightime = nightend - nightstart
    nightime.format = 'sec'
    #Moving to start of the night if in daytime
    if (sun_altaz.alt > -18*u.deg):
        timetilldark = (nightstart-t)
        timetilldark.format = 'sec'
        LST = nightstart.sidereal_time('mean').deg


        
    # How likely is it that the (true, unknown) location of the source
    # is within the area that is visible, in this time and within 24 hours? 
    # Demand that it falls in the HETDEX pupil, that the sun is at least 18 
    # degrees below the horizon and that the airmass (secant of zenith angle 
    # approximation) is at most 2.5.
    
    msortedpix = np.flipud(np.argsort(m)) 
    cumsum = np.cumsum(m[msortedpix]) 
    cls = np.empty_like(m) 
    cls[msortedpix] = cumsum*100
    p90i =[]
    for i in range(len(m)):
        if cls[i] <= 90:
            p90i.append(i)
    p90i = np.array(p90i)
    theta90, phi90 = hp.pix2ang(nside, p90i)

    #if not crossing
    #first find if distribution crosses 0 ra
    
    #ensure no coordinate gaps in hetpupil RAs
    #hetpupil[:,1] = ((hetpupil[:,1] + 180) % 360 - 180)
    HETtheta = (90-hetpupil[:,2])*np.pi/180
    HETphi = (hetpupil[:,1]+LST)*np.pi/180
    newpix = hp.ang2pix(nside, HETtheta, HETphi)
    
    timetill90 = 0
    #import pdb
    #pdb.set_trace()
    if not len(np.intersect1d(p90i,newpix))>0:
        if not len(np.intersect1d(p90i,hetfullpix))>0:
            return 0 , 0 , -99
        #ensure no coordinate gaps in source RAs
        #if (np.max(phi90) - np.min(phi90)) > 359*(np.pi/180):
        #    phi90 = ((phi90 + 180) % 360 - 180)
        phi90, HETphi = no_gaps(phi90, HETphi)
        iminth = np.argmin(HETtheta - theta90[np.argmin(phi90)])
        #at the right declination
        wsecs = (np.min(phi90)+np.pi-np.max(HETphi[iminth]))*12*3600/np.pi 
        if wsecs > nightime.value:
            return 0 , 0 , -99
        else:
            timetill90 = (wsecs+timetilldark.value)/3600
            #print('{:.1f} hours till you can observe the 90 %% prob region.'.format(
            #    (wsecs+timetilldark)/3600))
        #else if np.min(phi90) > np.max(HETphi):a
    elif timetilldark.value > 0:
       # print('{:.1f} hours till you can observe the 90 %% prob region.'.format(
        #        timetilldark.value/3600))
        timetill90 = timetilldark.value/3600
        
    
    mask_arraynow = np.zeros(len(m), dtype=int)
    mask_arraynow[newpix] = 1
    mask_arraynow *= (altaz.secz <= 2.5)&(sun_altaz.alt <= -18*u.deg)

    
    prob = m[mask_arraynow > 0].sum()
    probfull = m[np.intersect1d(p90i,hetfullpix)].sum()
    
    if plot:
        #Return to LST now for plotting
        LST = t.sidereal_time('mean').deg
        HETphi = (hetpupil[:,1]+LST)*np.pi/180
        newpix = hp.ang2pix(nside, HETtheta, HETphi)
        
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
        hp.mollview(m, coord='C', cbar=False, max=1, title='HET NOW',)
        hp.graticule(local=True)

        plt.savefig('MOLL_GWHET_%s.pdf'%header['GraceID'])
        #plt.show()

    # Done!
    return prob, probfull, timetill90


