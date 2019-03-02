#Mainly a simplified copy of HET_obs.py
#https://github.com/sjanowiecki/HET_observability
#and the ligo skymaps tutorials
#https://github.com/gw-odw/odw-2018/tree/master/skymaps

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

    parser = argparse.ArgumentParser(description='FIND GALAXIES TO OBSERVE IN TWO CATALOGS')
    parser.add_argument('--http', dest='fits', default='https://dcc.ligo.org/public/0146/G1701985/001/LALInference_v2.fits.gz', action=GetLoc, help='HTTPS link to LIGO event localization. It will download the file if not cached.')
    parser.add_argument('-cat', dest='cat', default='2MASS', action=GetLoc, help='Specify which catalog to use: 2MASS or GLADE')
    args = parser.parse_args()

    return args

def write_catalog2MASS(fits):

        # Reading in the skymap prob and header
        locinfo, header = hp.read_map(fits, field=range(4), h=True)
        prob, distmu, distsigma, distnorm = locinfo
        #Getting healpix resolution and pixel area in deg^2
        npix = len(prob)
        nside = hp.npix2nside(npix)

        # Area per pixel in steradians
        pixarea = hp.nside2pixarea(nside)
        Vizier.ROW_LIMIT = -1 # This gets the complete catalog
        cat1, = Vizier.get_catalogs('J/ApJS/199/26/table3') # Downloading the 2MRS Galaxy Catalog
        completeness = 0.5
        alpha = -1.0
        MK_star = -23.55
        MK_max = MK_star + 2.5*np.log10(gammaincinv(alpha + 2, completeness))
        z = (u.Quantity(cat1['cz'])/c.c).to(u.dimensionless_unscaled)
        MK = cat1['Ktmag']-cosmo.distmod(z)
        keep = (z > 0) & (MK < MK_max) & (cat1['DEJ2000']>-12.7)&(cat1['DEJ2000']<74.1)
        cat1 = cat1[keep]
        z = z[keep]
        r = cosmo.luminosity_distance(z).to('Mpc').value
        theta = 0.5*np.pi - cat1['DEJ2000'].to('rad').value
        phi = cat1['RAJ2000'].to('rad').value
        ipix = hp.ang2pix(nside, theta, phi)
        logdp_dV = np.log(prob[ipix]) + np.log(distnorm[ipix])+np.log( norm(distmu[ipix], distsigma[ipix]).pdf(r))-np.log(pixarea)

        top99i = logdp_dV-np.max(logdp_dV) > np.log(1/100)
        #Now working only with event with probability 99% lower than the most probable
        logdp_dV = logdp_dV[top99i]
        cattop = cat1[top99i]
        isort = np.argsort(logdp_dV)[::-1]
        cattop = cattop[isort]
        logptop = logdp_dV[isort]


        index = Column(name='index',data=np.arange(len(cattop)))
        logprob = Column(name='LogProb',data=logptop)
        exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
        Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
        #Normalizing the probability of the most probable galaxies
        cattop.add_columns([index,logprob,exptime,Nvis])
        ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb'], 'galaxies2MASS.dat', overwrite=True)
        return cattop,logptop


def write_catalogGLADE(fits):
    
        # Reading in the skymap prob and header
        locinfo, header = hp.read_map(fits, field=range(4), h=True)
        prob, distmu, distsigma, distnorm = locinfo
        #Getting healpix resolution and pixel area in deg^2
        npix = len(prob)
        nside = hp.npix2nside(npix)

        # Area per pixel in steradians
        pixarea = hp.nside2pixarea(nside)
        cat1 = pd.read_csv("./GLADE2.3HETz.csv", sep=',',usecols = [1,2,3],names=['RAJ2000','DEJ2000','z'],header=0,dtype=pd.np.float64)

        z = cat1['z']
        r = cosmo.luminosity_distance(z).to('Mpc').value
        theta = 0.5*np.pi - cat1['DEJ2000']*np.pi/180
        phi = cat1['RAJ2000']*np.pi/180
        ipix = hp.ang2pix(nside, theta, phi)
        logdp_dV = np.log(prob[ipix]) + np.log(distnorm[ipix])+np.log( norm(distmu[ipix], distsigma[ipix]).pdf(r))-np.log(pixarea)

        top99i = logdp_dV-np.max(logdp_dV) > np.log(1/100)
        #Now working only with event with probability 99% lower than the most probable
        logdp_dV = logdp_dV[top99i]
        cattop = cat1[top99i]
        isort = np.argsort(logdp_dV)[::-1]
        cattop = Table.from_pandas(cattop.iloc[isort])
        logptop = logdp_dV[isort]
        index = Column(name='index',data=np.arange(len(cattop)))
        logprob = Column(name='LogProb',data=logptop)
        exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
        Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
        #Normalizing the probability of the most probable galaxies
        cattop.add_columns([index,logprob,exptime,Nvis])
        ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb'], 'galaxiesGLADE.dat', overwrite=True)
        return cattop,logptop

def main():

    args = parseargs()
    if args.cat == 'GLADE':
        write_catalogGLADE(args.fits)
    elif args.cat == '2MASS':
        write_catalog2MASS(args.fits)
    else:
        print('Must specify either GLADE or 2MASS as catalogs.')

if __name__== "__main__":
    main()

