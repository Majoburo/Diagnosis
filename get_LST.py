#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 13:48:25 2019

@author: majoburo
"""
import numpy as np # needed for vector operations
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt # plotting skymaps
from astropy.io import  ascii
import astropy.units as u
from matplotlib.ticker import MultipleLocator
from matplotlib import rc


HET_track = 'HET_opt_tracking.txt'
obs_plot = True

def get_LST(targf = 'galaxies2MASS.dat'):
    GraceID = targf.split('_')[1].split('.')[0]
    #trimester year and number
    c_y = 2019
    c_t = 1

    #####----------------------------#####
    #NOTE: only works for 2019-1 right now, with HETDEX data included.





    #targets file with headers: ID, RA, Dec, exptime, Nvisits, Prob (etc ok):


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





    plt.register_cmap(name='viridis', cmap=plt.cm.viridis)

    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)


    def dolab1(fig):#:,thisplot):
            #18, far over!
            thisplot = fig.add_axes([0.824,0.06,0.05,0.05],visible=True)
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
            thisplot = fig.add_axes([0.628,0.06,0.05,0.05],visible=True)
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
            thisplot = fig.add_axes([0.438,0.06,0.05,0.05],visible=True)
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
    #Only show 50 most probable galaxies
    if len(targs)>50:
        targs = targs[:50]
    #parse. IDs first, assume
    targ_id = targs.columns[0].data
    targ_ra = targs.columns[1].data
    targ_dec = targs.columns[2].data
    targ_exptime = targs.columns[3].data
    targ_nvis = targs.columns[4].data
    targ_logprob = targs.columns[5].data

    #read in HET observability data file
    hetf = HET_track
    het = ascii.read(hetf)
    #print(het)
    h_dec = het.columns[0].data
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
            dd = abs(h_dec - d)
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
                #hamid = hami + (hama-hami)/2.
                #print(ha1,ha2,hamid)
                
                #VERIFY! is nthere enough time for this requested exptime within each visit?
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
                
                #VERIFY! is there enough time for this requested exptime within each visit?
                ha_total = hama-hami #in hours
                req_h = tott/(3600*u.s)
                if (req_h > ha_total):  #should never happen here!
                    print("{:.2f}h on second track available. you requested: {:.2f}h per visit".format(ha_total,req_h))
                    print("   error. still not enough visits....")

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
        
        
        
        
    if obs_plot:
        fig = plt.figure(figsize=(13,6))
        a = fig.add_subplot(111)
            
        plt.subplots_adjust(wspace=0.25, hspace=0)

        minorLocator1 = MultipleLocator(1) #x ticks/etc
        majorLocator1 = MultipleLocator(6)
        minorLocator2 = MultipleLocator(0.25) #y ticks/etc
        majorLocator2 = MultipleLocator(1)
        a.xaxis.set_minor_locator(minorLocator1)
        a.xaxis.set_major_locator(majorLocator1)
        #a.yaxis.set_minor_locator(minorLocator2)
        #a.yaxis.set_major_locator(majorLocator2)


        a.tick_params(axis='both', which='major', labelsize=16)

        a.set_xlim(offs,24+offs)
        ystep = (max(targ_logprob)-min(targ_logprob))/10
        a.set_ylim(min(targ_logprob)-ystep,max(targ_logprob)+ystep)

        a.set_xlabel(r'LST [h]', fontname='Arial', fontsize=22, fontweight='normal')
        a.set_ylabel(r'Probability (LOG SCALE)', fontname='Arial', fontsize=22, fontweight='normal')
            
#fix manually in two cases:
        if (c_t==1): dolab1(fig)
        if (c_t==2): dolab2(fig) #none needed
        if (c_t==3): dolab3(fig)



    with open('LSTs_%s.out'%GraceID,'w') as f:
        f.write("ID RA DEC LST1_start LST1_stop LST2_start LST2_stop Nvis Exptime LogProb \n")
        for ra,dec,t1,t2,t3,t4,i,nv,texp,tprob in zip(targ_ra, targ_dec, LST1_start, LST1_stop, LST2_start, LST2_stop, targ_id, targ_nvis,targ_exptime, targ_logprob):
            #if valid, add first trajectory for as many visits as needed
            if obs_plot:
                if ((t1>-30)&(t2>-30)&(t3>-30)&(t4>-30)):
                    #two tracks! are they both definitely in darkness w/ appropriate moon?
                    # for now, ALWAYS include, even if impossible....
                    # add first track to appropriate vector
                    if t1 > 24+offs: #we wrapped around
                        t1 %= 24
                    if t2 > 24+offs:
                        t2 %= 24
                    if t3 > 24+offs:
                        t3 %= 24
                    if t4 > 24+offs:
                        t4 %= 24
                    if t2 < t1:
                        a.plot(np.linspace(offs,t2,2), [tprob,tprob],color='red',lw=1)
                        a.plot(np.linspace(t1,24+offs,2), [tprob,tprob],color='red',lw=1)
                    else:
                        a.plot(np.linspace(t1,t2,2), [tprob,tprob],color='red',lw=1)
                    if t4<t3:
                        a.plot(np.linspace(-0.5+offs,t4,2), [tprob,tprob],color='green',lw=1)
                        a.plot(np.linspace(t3,24.5+offs,2), [tprob,tprob],color='green',lw=1)
                    else:
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


    if obs_plot:
    #labels:
        a.text(0.5+offs,np.max(targ_logprob),'Single tracks (N/S)',color='blue', fontname='Arial', fontsize=19, fontweight='bold',alpha=0.8)
        a.text(0.5+offs,np.max(targ_logprob)-ystep,'West tracks',color='red', fontname='Arial', fontsize=19, fontweight='normal',alpha=0.4)
        a.text(0.5+offs,np.max(targ_logprob)-2*ystep,'East tracks',color='green', fontname='Arial', fontsize=19, fontweight='normal',alpha=0.4)
    #note names are off.... but ok.

    #horizontal lines
        llines = np.arange(np.min(targ_logprob),np.max(targ_logprob),ystep)
        for l in llines:
                a.plot([-1+offs,25+offs],[l,l],color='grey',lw=1,ls=':',alpha=0.8)


        graout = 'LSTs_'+GraceID+'.pdf'

        fig.savefig(graout, bbox_inches='tight')

def main():
    get_LST(targf = 'galaxies2MASS.dat')

if __name__== "__main__":
    main()
