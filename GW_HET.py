#Mainly a simplified copy of HET_obs.py 
#https://github.com/sjanowiecki/HET_observability
#and the ligo skymaps tutorials
#https://github.com/gw-odw/odw-2018/tree/master/skymaps

import healpy as hp # for working with HEALPix files
import numpy as np # needed for vector operations
from matplotlib import pyplot as plt # plotting skymaps
from scipy.stats import norm # probability functions
import sys
from astropy.utils.data import download_file
from astropy.time import Time
from astropy.io import fits, ascii
import argparse
import datetime
import math
from astropy.table import Table
from astropy.table import Column
from astroquery.vizier import Vizier
from scipy.special import gammaincinv
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import astropy.constants as c
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
from astropy.coordinates import get_sun

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

    parser = argparse.ArgumentParser(description='FIND GALAXIES TO OBSERVE AND TIMES TO DO SO AFTER A LIGO EVENT')
    parser.add_argument('-mp', dest='moll_plot_true', action='store_true', help = 'Plot Mollweid projection of 90 %% probable region for event and HET')
    parser.add_argument('-p', dest='obs_plot', action='store_true', help = 'Plot OBS plan')
    parser.add_argument('--LIGO', dest='loc_ligo', default='https://dcc.ligo.org/public/0146/G1701985/001/LALInference_v2.fits.gz', action=GetLoc, help='HTTPS link to LIGO event localization')
    parser.add_argument('--HET', dest='HET_track',default='HET_opt_tracking.txt')
    args = parser.parse_args()

    return args



def main():
    args = parseargs()
    if True:

        # Reading in the skymap prob and header
        locinfo, header = hp.read_map(args.loc_ligo, field=range(4), h=True)
        prob, distmu, distsigma, distnorm = locinfo
        #Getting healpix resolution and pixel area in deg^2
        npix = len(prob)
        nside = hp.npix2nside(npix)
        pixarea_deg2 = hp.nside2pixarea(nside, degrees=True)
        probperdeg2 = prob/pixarea_deg2

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
        keep = (z > 0) & (MK < MK_max)
        cat1 = cat1[keep]
        z = z[keep]
        r = cosmo.luminosity_distance(z).to('Mpc').value
        theta = 0.5*np.pi - cat1['DEJ2000'].to('rad').value
        phi = cat1['RAJ2000'].to('rad').value
        ipix = hp.ang2pix(nside, theta, phi)
        dp_dV = prob[ipix] * distnorm[ipix] * norm(distmu[ipix], distsigma[ipix]).pdf(r)/pixarea
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
        ascii.write(cattop['index','RAJ2000','DEJ2000','exptime','Nvis','LogProb'], 'galaxies.dat', overwrite=True)

        #HETDEX WINDOW
        # Getting the local sideral time right now
        if args.moll_plot_true:
            HET_loc = ('-104.01472d','30.6814d','2025m')
            d = datetime.datetime.now(tz=datetime.timezone.utc)
            #t = Time(dict(header)['DATE-OBS'],scale='utc',location=HET_loc)
            t = Time(d,scale='utc',location=('-104.01472d','30.6814d','2025m'))
            LST = np.deg2rad(t.sidereal_time('mean').deg)
            #Getting HETDEX track
            het = ascii.read(args.HET_track)

            h_dec = het.columns[0].data
            h_tott = het.columns[1].data
            h_optaz = het.columns[2].data
            h_ha1 = het.columns[3].data
            h_ha2 = het.columns[4].data
            h_max = np.max(np.array(list(zip(h_ha1,h_ha2))),axis=1)
            h_min = np.min(np.array(list(zip(h_ha1,h_ha2))),axis=1)

            #Convert to radians
            theta = 0.5*np.pi - np.deg2rad(h_dec)
            phi2 = np.deg2rad(h_max*360/24)
            phi1 = np.deg2rad(h_min*360/24)
            phiss = [np.linspace(p1,p2,100) for p1,p2 in zip(phi1,phi2)]
            hetang = []
            for th,phis in zip(theta,phiss):
                for phi in phis:
                    hetang.append([th, +phi+LST]) #NOTE! converting to RA
                    hetang.append([th, -phi+LST]) #NOTE! converting to RA
            hetang = np.array(hetang)
            hetnow = hp.ang2pix(nside, hetang[:,0], hetang[:,1])

            hetfull = hp.query_strip(nside, (90-max(h_dec))*np.pi/180, (90-min(h_dec))*np.pi/180)
    
    
            #Getting 90% prob region
    
            i = np.flipud(np.argsort(prob)) # Array of pixel indices ordered from highest to lowest probability densities
            cumsum = np.cumsum(prob[i]) # Array of the cumulative sums of the ordered probability density pixels
            cls = np.empty_like(prob) # A new array with the same shape and type as prob
            cls[i] = cumsum*100 # The new array filled in with the value of the cumulative sums by the time we reach that pixel index
    
            for i in range(len(prob)):
                if cls[i] <= 90:
                    probperdeg2[i] = 1
    
            probperdeg2[hetfull]=np.max(probperdeg2)/10
            probperdeg2[hetnow]=np.max(probperdeg2)/2
            hp.mollview(probperdeg2, coord=['C'], title='GW170817 LALInference - HET', max=np.max(probperdeg2),cbar=False)
            hp.graticule(local=True)
            plt.savefig('MOLL_GWHET.pdf')
            plt.show()
    
    
    #trimester year and number
    c_y = 2019
    c_t = 1

    #####----------------------------#####
    #NOTE: only works for 2019-1 right now, with HETDEX data included.






    #targets file with headers: ID, RA, Dec, exptime, Nvisits, Prob (etc ok):
    targf = 'galaxies.dat' #change if needed

    setup_time = 300 #assume 300 seconds. 







    #set alias/wrap-around point based on trimester
    if (c_t==1):
        #hwrap = +18.0 #NOT USED
        offs= 20#-8 #limit offsets
    elif (c_t==2):
        ##hwrap = +8.0 #i.e., values below 8 get +24. lims shift by hwrap
        offs=+2#+8 #limit offsets
    elif (c_t==3):
        ##hwrap = 12.0 #NOT USED
        offs=+12 #0 #limits/save_LST offsets
    else:
        print("issue with trimester")
        asplodez





    plt.register_cmap(name='viridis', cmap=plt.cm.viridis)

    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)






    def dolab1(fig):#:,thisplot):
            #18, far over!
            thisplot = fig.add_axes([0.808,0.05,0.05,0.05],visible=True)
            thisplot.spines['bottom'].set_color('white')
            thisplot.spines['top'].set_color('white')
            thisplot.spines['left'].set_color('white')
            thisplot.spines['right'].set_color('white')
            thisplot.xaxis.set_visible(False)
            thisplot.yaxis.set_visible(False)
            thisplot.set_zorder(500)
            thisplot.fill('white')
            thisplot.text(-0.045,-0.015,'18',fontsize=17)
            #12, near middle!
            thisplot = fig.add_axes([0.623,0.05,0.05,0.05],visible=True)
            thisplot.spines['bottom'].set_color('white')
            thisplot.spines['top'].set_color('white')
            thisplot.spines['left'].set_color('white')
            thisplot.spines['right'].set_color('white')
            thisplot.xaxis.set_visible(False)
            thisplot.yaxis.set_visible(False)
            thisplot.set_zorder(500)
            thisplot.fill('white')
            thisplot.text(-0.045,-0.015,'12',fontsize=17)
            #6, at left
            thisplot = fig.add_axes([0.442,0.05,0.05,0.05],visible=True)
            thisplot.spines['bottom'].set_color('white')
            thisplot.spines['top'].set_color('white')
            thisplot.spines['left'].set_color('white')
            thisplot.spines['right'].set_color('white')
            thisplot.xaxis.set_visible(False)
            thisplot.yaxis.set_visible(False)
            thisplot.set_zorder(500)
            thisplot.fill('white')
            thisplot.text(-0.045,-0.015,'6',fontsize=17)
            return fig#,thisplot

    def dolab2(fig):
            a=0
            if (False):
                thisplot = fig.add_axes([0.80,0.05,0.05,0.05],visible=True)
                thisplot.spines['bottom'].set_color('white')
                thisplot.spines['top'].set_color('white')
                thisplot.spines['left'].set_color('white')
                thisplot.spines['right'].set_color('white')
                thisplot.xaxis.set_visible(False)
                thisplot.yaxis.set_visible(False)
                thisplot.set_zorder(500)
                thisplot.fill('grey')
                thisplot.set_alpha(0.5)
                thisplot.text(-0.02,0.0,'6',fontsize=17)
            return fig#thisplot

    def dolab3(fig):
            thisplot = fig.add_axes([0.680,0.045,0.05,0.05],visible=True)
            thisplot.spines['bottom'].set_color('white')
            thisplot.spines['top'].set_color('white')
            thisplot.spines['left'].set_color('white')
            thisplot.spines['right'].set_color('white')
            thisplot.xaxis.set_visible(False)
            thisplot.yaxis.set_visible(False)
            thisplot.set_zorder(500)
            thisplot.fill('grey')
            thisplot.set_alpha(0.5)
            thisplot.text(-0.02,0.0,'6',fontsize=15)

            thisplot2 = fig.add_axes([0.860,0.045,0.05,0.05],visible=True)
            thisplot2.spines['bottom'].set_color('white')
            thisplot2.spines['top'].set_color('white')
            thisplot2.spines['left'].set_color('white')
            thisplot2.spines['right'].set_color('white')
            thisplot2.xaxis.set_visible(False)
            thisplot2.yaxis.set_visible(False)
            thisplot2.set_zorder(500)
            thisplot2.fill('grey')
            thisplot2.set_alpha(0.5)
            thisplot2.text(-0.02,0.0,'12',fontsize=15)
            return fig#,thisplot





    #read in target ra/dec, exptime, nvis, moon
    targs = ascii.read(targf)

    #parse. IDs first, assume
    targ_id = targs.columns[0].data
    targ_ra = targs.columns[1].data
    targ_dec = targs.columns[2].data+90
    targ_exptime = targs.columns[3].data
    targ_nvis = targs.columns[4].data
    targ_logprob = targs.columns[5].data

    #read in HET observability data file
    hetf = args.HET_track
    het = ascii.read(hetf)
    #print(het)
    h_dec = het.columns[0].data
    h_tott = het.columns[1].data
    h_optaz = het.columns[2].data
    h_ha1 = het.columns[3].data
    h_ha2 = het.columns[4].data

#add ha3/ha4 for tracks of appropriate dec:
#  dec boundaries of double-valued?
#    -4.318553207530732 < dec < 65.6814360000000
    d2min = -4.318553207530732 
    d2max = 65.6814360000000
    h_ha3 = np.array([-h if ((d>d2min)&(d<d2max)) else -99 for h,d in zip(h_ha2,h_dec)])
    h_ha4 = np.array([-h if ((d>d2min)&(d<d2max)) else -99 for h,d in zip(h_ha1,h_dec)])

#determine LSTs of available tracks for these targets
# and prepare arrays for saving
    LST1_start = targ_ra*0.0-99.
    LST1_stop  = targ_ra*0.0-99.
    LST2_start = targ_ra*0.0-99.
    LST2_stop  = targ_ra*0.0-99.
    for i,r,d,e,n in zip(targ_id,targ_ra,targ_dec,targ_exptime,targ_nvis):
            print('TARGET: '+str(i))
            #ok, so for this ra,dec (and exptime/nvis), find time per visit
            epv = e/n
            #also, get RA in decimal hours
            ra_h = r/15.
            
            #consider only a single visit, since will be identical
            #assumed setup time is: setup_time (in seconds)
            # so, total time for a single visit is:
            tott = (epv+setup_time)*u.s
            
            
            
            #for this dec, find the optimal start/stop times in LST
            # so, find closest h_dec and get the HA values (1,2,3,4)
            dd = np.abs(h_dec - d)
            #Is it close enough?
            sh_dec = np.sort(h_dec)
            max_dh_dec = max(abs(sh_dec[1:]-sh_dec[:-1]))
            if dd.min() > max_dh_dec:
                print('99% confidence region not visible by HET, ¯\_(ツ)_/¯')
                continue
            ha1 = h_ha1[dd==dd.min()][0]
            ha2 = h_ha2[dd==dd.min()][0]
            ha3 = h_ha3[dd==dd.min()][0]
            ha4 = h_ha4[dd==dd.min()][0]

            #verify that ha1/ha2 are valid: i.e. the first (usu west) track)
            if ((np.abs(ha1)<5)&(np.abs(ha2)<5)):
                #good!  these are always set up so ha1 is before ha2?
                hami = np.min([ha1,ha2])
                hama = np.max([ha1,ha2])
                #find midpoint of track:
                hamid = hami + (hama-hami)/2.
                #print(ha1,ha2,hamid)
                
                #VERIFY! is there enough time for this requested exptime within each visit?
                ha_total = hama-hami #in hours
                req_h = tott/(3600*u.s)
                if (req_h > ha_total):
                    print(' ')
                    #print('target: '+str(i))
                    print(ha1,ha2,ha3,ha4)
                    print("{:.2f}h on first track available. you requested: {:.2f}h per ({:.0f}) visit".format(ha_total, req_h.value,n))
                    #how many visits would it need?
                    fix_nv = np.int(1.+(e/3600.)/ha_total)
                    print('   that is a PROBLEM -- need more visits. for now you get '+str(fix_nv)+' visits instead')
                    print('epv before '+str(epv))
                    epv = e/fix_nv
                    print('epv after '+str(epv))
                    n=fix_nv
                    targ_nvis[i] = n

                    tott = (epv+setup_time)*u.s
                
                #use total time to extend on either side of this, optimal.
                ha_start = hami*u.h#hamid*u.h - tott/2.
                ha_stop = hama*u.h#hamid*u.h + tott/2.
                #print(ha_start,ha_stop,tott)
                
                
                #combine with target RA to get LST start/stop
                LST1_start[targ_id==i] = ra_h + ha_start/u.h
                LST1_stop[targ_id==i] = ra_h + ha_stop/u.h

            #verify that ha3/ha4 are valid: i.e. the second (usu east) track)
            if ((np.abs(ha3)<5)&(np.abs(ha4)<5)):
                #good!  set up min/max
                hami = np.min([ha3,ha4])
                hama = np.max([ha3,ha4])
                #find midpoint of track:
                hamid = hami + (hama-hami)/2.
                #print(ha3,ha4,hamid)
                
                #VERIFY! is there enough time for this requested exptime within each visit?
                ha_total = hama-hami #in hours
                req_h = tott/(3600*u.s)
                if (req_h > ha_total):  #should never happen here!
                    print("{:.2f}h on second track available. you requested: {:.2f}h per visit".format(ha_total,req_h))
                    print("   error. still not enough visits....")
                    asplode

                #use total time to extend on either side of this, optimal.
                ha_start = hami*u.h#hamid*u.h - tott/2.
                ha_stop = hama*u.h#hamid*u.h + tott/2.
                #print(ha_start,ha_stop,tott)
                
                
                #combine with target RA to get LST start/stop
                LST2_start[targ_id==i] = ra_h + ha_start/u.h
                LST2_stop[targ_id==i] = ra_h + ha_stop/u.h
            else:
                print("No second track")
            #done:
            #print(' ')

#done:
    print(' ')
    print(' ')
        
        
#fix up rounding issues/etc. depends on trimester.
#print('checking')
    if (c_t==1):
            #for trimester 1, 
            tmp = np.array([-99 if (l==-99) else l if (l>22) else l+24.0 for l in LST1_start])
            LST1_start = tmp
            tmp = np.array([-99 if (l==-99) else l if (l>22) else l+24.0 for l in LST1_stop])
            LST1_stop = tmp
            tmp = np.array([-99 if (l==-99) else l if (l>22) else l+24.0 for l in LST2_start])
            LST2_start = tmp
            tmp = np.array([-99 if (l==-99) else l if (l>22) else l+24.0 for l in LST2_stop])
            LST2_stop = tmp
    if (c_t==2):
            #for trimester 2, we keep 0-24  ####OLD: allow up to 32h, and take negatives +24
            tmp = np.array([-99 if (l==-99) else l if ((l>=2)&(l<=25)) else l-24.0 if (l>25) else l+24.0 for l in LST1_start])
            LST1_start = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>=2)&(l<=25)) else l-24.0 if (l>25) else l+24.0 for l in LST1_stop])
            LST1_stop = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>=2)&(l<=25)) else l-24.0 if (l>25) else l+24.0 for l in LST2_start])
            LST2_start = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>=2)&(l<=25)) else l-24.0 if (l>25) else l+24.0 for l in LST2_stop])
            LST2_stop = tmp
    if (c_t==3):
            #for trimester 3, less than 12, send up to >24
            tmp = np.array([-99 if (l==-99) else l if ((l>12)&(l<24)) else l+24.0 if (l<12) else l-24.0 for l in LST1_start])
            LST1_start = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>12)&(l<24)) else l+24.0 if (l<12) else l-24.0 for l in LST1_stop])
            LST1_stop = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>12)&(l<24)) else l+24.0 if (l<12) else l-24.0 for l in LST2_start])
            LST2_start = tmp
            tmp = np.array([-99 if (l==-99) else l if ((l>12)&(l<24)) else l+24.0 if (l<12) else l-24.0 for l in LST2_stop])
            LST2_stop = tmp
        
        
        
        
    if args.obs_plot:
        fig = plt.figure(figsize=(13,6))
        a = fig.add_subplot(111)
            
        plt.subplots_adjust(wspace=0.25, hspace=0)

        minorLocator1 = MultipleLocator(1) #x ticks/etc
        majorLocator1 = MultipleLocator(6)
        minorLocator2 = MultipleLocator(0.25) #y ticks/etc
        majorLocator2 = MultipleLocator(1)
        a.xaxis.set_minor_locator(minorLocator1)
        a.xaxis.set_major_locator(majorLocator1)
        a.yaxis.set_minor_locator(minorLocator2)
        a.yaxis.set_major_locator(majorLocator2)


        a.tick_params(axis='both', which='major', labelsize=16)

        a.set_xlim(-0.5+offs,24.5+offs)
        ystep = (max(targ_logprob)-min(targ_logprob))/10
        a.set_ylim(min(targ_logprob)-ystep,max(targ_logprob)+ystep)

        a.set_xlabel(r'LST [h]', fontname='Arial', fontsize=22, fontweight='normal')
        a.set_ylabel(r'Probability (LOG SCALE)', fontname='Arial', fontsize=22, fontweight='normal')
            
#fix manually in two cases:
        if (c_t==1): dolab1(fig)
        if (c_t==2): dolab2(fig) #none needed
        if (c_t==3): dolab3(fig)



    with open('exptimes.out','w') as f:
        f.write("ID RA DEC LST1_start LST1_stop LST2_start LST2_stop Nvis Exptime LogProb \n")
        for ra,dec,t1,t2,t3,t4,i,nv,texp,tprob in zip(targ_ra, targ_dec, LST1_start, LST1_stop, LST2_start, LST2_stop, targ_id, targ_nvis,targ_exptime, targ_logprob):
            #if valid, add first trajectory for as many visits as needed
            if args.obs_plot:
                if ((t1>-30)&(t2>-30)&(t3>-30)&(t4>-30)):
                    #two tracks! are they both definitely in darkness w/ appropriate moon?
                    # for now, ALWAYS include, even if impossible....
                    # add first track to appropriate vector
                    a.plot(np.linspace(t1,t2,2), [tprob,tprob],color='red',lw=1)
                    a.plot(np.linspace(t3,t4,2), [tprob,tprob],color='green',lw=1)
                    a.text(t4+0.1, tprob,i,color='black', fontname='Arial',fontsize=10)
                    # add second track to appropriate vector
                    #print("  two tracks good")
                elif ((t1>-30)&(t2>-30)):
                    #only 1 is good. it is #1:
                    #add to singles
                    a.plot(np.linspace(t1,t2,2), [tprob,tprob],color='blue',lw=1)
                    a.text(t2+0.1, tprob,i,color='black', fontname='Arial',fontsize=10)
                    #print("  track 1 good")
                elif ((t3>-30)&(t4>-30)):
                    a.plot(np.linspace(t3,t4,2), [tprob,tprob],color='blue',lw=1)
                    a.text(t4+0.1, tprob,i,color='black', fontname='Arial',fontsize=10)
                    #print("  track 2 good")
            outp = "%i %6.f %.6f %.6f %.6f %.6f %.6f %i %i %.6f \n"%(i,ra,dec,t1,t2,t3,t4,nv,texp,tprob) 
            f.write(outp)


    if args.obs_plot:
    #labels:
        a.text(0+offs,np.max(targ_logprob),'Single tracks (N/S)',color='blue', fontname='Arial', fontsize=19, fontweight='bold',alpha=0.8)
        a.text(0+offs,np.max(targ_logprob)-ystep,'West tracks',color='red', fontname='Arial', fontsize=19, fontweight='normal',alpha=0.4)
        a.text(0+offs,np.max(targ_logprob)-2*ystep,'East tracks',color='green', fontname='Arial', fontsize=19, fontweight='normal',alpha=0.4)
    #note names are off.... but ok.

    #horizontal lines
        llines = np.arange(np.min(targ_logprob),np.max(targ_logprob),ystep)
        for l in llines:
                a.plot([-1+offs,25+offs],[l,l],color='grey',lw=1,ls=':',alpha=0.8)


        graout = 'LST_visits_PI_'+str(c_y)+'-'+str(c_t)+'.pdf'
        fig.savefig(graout, bbox_inches='tight')



if __name__== "__main__":
    main()

