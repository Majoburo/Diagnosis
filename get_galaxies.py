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
from ligo.skymap.distance import conditional_pdf

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
    parser.add_argument('-cat', dest='cat', default='MANGROVE', help='Specify which catalog to use: MANGROVE or GLADE')
    args = parser.parse_args()

    return args

def cdf(pdf):
    #Calculate contour in probability
    sortedpix = np.argsort(pdf)[::-1]
    cumsum = np.cumsum(pdf[sortedpix])
    cls = np.empty_like(pdf)
    cls[sortedpix] = cumsum*100
    return cls

def write_catalog(params,catalog):
    fits = params['skymap_fits']
    event = params['GraceID']
    probability = params['skymap_array']

    if catalog == 'MANGROVE':

        # Reading in the skymap prob and header
        locinfo, header = hp.read_map(fits, field=range(4), h=True)
        probb, distmu, distsigma, distnorm = locinfo

        # Getting healpix resolution and pixel area in deg^2
        npix = len(probb)
        nside = hp.npix2nside(npix)
        # Area per pixel in steradians
        pixarea = hp.nside2pixarea(nside)
        # Get the catalog
        cat1 = pd.read_csv("./mangroveHET.csv", sep=',',usecols = [1,2,3,4],names=['RAJ2000','DEJ2000','d','m'],header=0,dtype=pd.np.float64)
        #cat1 = pd.read_csv("./GLADE2.3.csv", sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','d','B_Abs','K_Abs'],header=0,dtype=pd.np.float64)
        theta = 0.5*np.pi - cat1['DEJ2000']*np.pi/180
        phi = cat1['RAJ2000']*np.pi/180
        cls = cdf(probb)
        ipix = hp.ang2pix(nside, theta, phi)
        cls = cls[ipix]

        dist = cat1['d']
        logdp_dV = np.log(probability[ipix]) + np.log(conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()) - np.log(pixarea)

        #cutting to select only 90 % confidence in position
        cattop = cat1[cls<90]
        logdp_dV= logdp_dV[cls<90]
        s_m = cat1['m'][cls<90]
        s_m = s_m/s_m.sum()
        cls = cls[cls<90]
        logdp_dV = np.log(s_m) + logdp_dV

        #Now working only with event with overall probability 99% lower than the most probable
        top99i = logdp_dV-np.max(logdp_dV) > np.log(1/100)

        cattop = cattop[top99i]
        logdp_dV = logdp_dV[top99i]
        cls = cls[top99i]

        #sorting by probability
        isort = np.argsort(logdp_dV)[::-1]
        cattop = Table.from_pandas(cattop.iloc[isort])
        logptop = logdp_dV.iloc[isort]
        cls = cls[isort]

        index = Column(name='index',data=np.arange(len(cattop)))
        logprob = Column(name='LogProb',data=logptop)
        exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
        contour = Column(name='contour',data = cls)
        Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
        cattop.add_columns([index,logprob,exptime,Nvis,contour])
        ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb','contour'], 'galaxiesMANGROVE_%s.dat'%event, overwrite=True)
        return cattop,logptop

    if catalog == 'GLADE':

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
        #cat1 = pd.read_csv("./GLADE2.3.csv", sep=',',usecols = [1,2,3,4,5],names=['RAJ2000','DEJ2000','d','B_Abs','K_Abs'],header=0,dtype=pd.np.float64)
        theta = 0.5*np.pi - cat1['DEJ2000']*np.pi/180
        phi = cat1['RAJ2000']*np.pi/180
        cls = cdf(probb)

        ipix = hp.ang2pix(nside, theta, phi)
        cls = cls[ipix]

        dist = cat1['d']
        logdp_dV = np.log(probability[ipix]) + np.log(conditional_pdf(dist,distmu[ipix],distsigma[ipix],distnorm[ipix]).tolist()) - np.log(pixarea)

        #cutting to select only 90 % confidence in position
        cattop = cat1[cls<90]
        logdp_dV= logdp_dV[cls<90]
        s_lumK = 10**(-0.4*cat1['K_Abs'][cls<90])
        s_lumK = s_lumK/s_lumK.sum()
        #s_lumB = 10**(-0.4*cat1['B_Abs'][cls>90])
        #s_lumB = s_lumB/s_lumB.sum()
        cls = cls[cls<90]
        #only using K for now
        logdp_dV = np.log(s_lumK) + logdp_dV

        #Now working only with event with overall probability 99% lower than the most probable
        top99i = logdp_dV-np.max(logdp_dV) > np.log(1/100)

        cattop = cattop[top99i]
        logdp_dV = logdp_dV[top99i]
        cls = cls[top99i]

        #sorting by probability
        isort = np.argsort(logdp_dV)[::-1]
        cattop = Table.from_pandas(cattop.iloc[isort])
        logptop = logdp_dV.iloc[isort]
        cls = cls[isort]

        index = Column(name='index',data=np.arange(len(cattop)))
        logprob = Column(name='LogProb',data=logptop)
        exptime = Column(name='exptime',data=60*20*np.ones(len(cattop)))
        contour = Column(name='contour',data = cls)
        Nvis = Column(name='Nvis',data=np.ones(len(cattop)))
        cattop.add_columns([index,logprob,exptime,Nvis,contour])
        ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb','contour'], 'galaxiesGLADE_%s.dat'%event, overwrite=True)
        return cattop,logptop

def main():

    args = parseargs()
    prob, header = hp.read_map(args.fits, h=True)
    header = dict(header)
    params = {'skymap_fits':args.fits,'skymap_array':prob,'GraceID':header['OBJECT']}
    if args.cat == 'MANGROVE' or args.cat == 'GLADE':
        write_catalog(params,args.cat)
    else:
        print('Must specify GLADE or MANGROVE as catalog.')

if __name__== "__main__":
    main()

