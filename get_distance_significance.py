import healpy as hp # for working with HEALPix files
import numpy as np # needed for vector operations
from scipy.stats import norm # probability functions
from astropy.utils.data import download_file
from astropy.io import ascii
import argparse
from astropy.table import Table
from astropy.table import Column
from astroquery.vizier import Vizier
from scipy.special import gammaincinv
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import astropy.constants as c
import pandas as pd
from ligo.skymap.distance import conditional_pdf
import matplotlib.pyplot as plt
def parseargs():

    class GetLoc(argparse.Action):
        def __init__(self, option_strings, dest, nargs=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs not allowed")
            super(GetLoc, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            url = (values)
            filename = download_file(url, cache=True)
            setattr(namespace, self.dest, filename)

    parser = argparse.ArgumentParser(description='GET VOLUME PROBABILITY OF RA DEC DIST')
    parser.add_argument('--http', dest='fits', default='https://dcc.ligo.org/public/0146/G1701985/001/LALInference_v2.fits.gz', action=GetLoc, help='HTTPS link to LIGO event localization. It will download the file if not cached.')
    parser.add_argument('coords',nargs = 3,type = float, help= 'RA DEC DIST in degrees and Mpc')
    args = parser.parse_args()

    return args

def cdf(pdf):
    #Calculate contour in probability
    sortedpix = np.flipud(np.argsort(pdf))
    cumsum = np.cumsum(pdf[sortedpix])
    cls = np.empty_like(pdf)
    cls[sortedpix] = cumsum*100
    return cls

def vol_prob(ra,dec,dist,fits):
# Reading in the skymap prob and header
    locinfo, header = hp.read_map(fits, field=range(4), h=True)
    probb, distmu, distsigma, distnorm = locinfo

# Getting healpix resolution and pixel area in deg^2
    npix = len(probb)
    nside = hp.npix2nside(npix)
# Area per pixel in steradians
    pixarea = hp.nside2pixarea(nside)
# Get the catalog
    cat1 = pd.read_csv("./GLADE2.3HETd.csv", sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','d','B_Abs','K_Abs'],header=0,dtype=pd.np.float64)
    theta = 0.5*np.pi - dec*np.pi/180
    phi = ra*np.pi/180
    cls = cdf(probb)
    ipix = hp.ang2pix(nside, theta, phi)
    cls = cls[ipix]
    logdp_dV = np.log(probb[ipix]) + np.log(conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()) - np.log(pixarea)
    dp_dD = conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()
    print("[SIGMA DIST]: %.1f"%((dist-distmu[ipix])/distsigma[ipix]))
    r = np.linspace(0, 500)
    p = norm(distmu[ipix], distsigma[ipix]).pdf(r) * distnorm[ipix] * r**2
    plt.plot(r, p)
    plt.axvline(dist)
    plt.xlabel('Distance (Mpc)')
    plt.ylabel('Probability / Mpc')
    plt.savefig('pg.png')

    return dp_dD


def main():

    args = parseargs()
    dp_dD = vol_prob(args.coords[0],args.coords[1],args.coords[2],args.fits)
if __name__== "__main__":
    main()
