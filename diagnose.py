import gcn
import healpy as hp
from prob_obs import prob_observable
from astropy.time import Time
from astropy.utils.data import download_file
import email_ip, get_galaxies, get_LST, make_phaseii
import os, urllib
import argparse
import urllib.request
import lxml
import numpy as np
plotting = True
params = {}
args = 0

def parseargs():

    parser = argparse.ArgumentParser(description='Receive and parse GCN alerts, alert observers and create observing tools.')
    parser.add_argument('recipients', help = 'Specify python file with list of recipients.')
    parser.add_argument('-cat', dest='cat', default = 'GLADE', help='Specify which catalog to use: 2MASS or GLADE')
    parser.add_argument('-fits', dest='fits', default = None, help='If given a skymap fits, override everything and use this fits file.')
    parser.add_argument('-t','--test', dest='test', action='count', default =0, help='DEFAULT:0. Run a test. With -tt will not query gcn but use a file stored in the git repo. With -t will query gcn for test alerts.')
    parser.add_argument('-e','--send_notification', dest = 'send_notification', action='store_false', help='DEFAULT:TRUE. Send emails and txt msjs to recipients.')
    args = parser.parse_args()
    args.cat = [x.replace(' ', '') for x in args.cat.split(',')]

    return args

# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE)

def process_gcn(payload, root):
    with open('time_last.txt', 'w') as f:
        f.write('%s'%Time.now().jd)
    if args.test > 0:
        if root.attrib['role'] != 'test':
            return
    elif root.attrib['role'] != 'observation' or root.attrib['role'] != 'test':
            return

    # Read all of the VOEvent parameters from the "What" section.
    global params
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    params['role'] = root.attrib['role']
    if params['role'] == 'test':
        recipients = 'recipients_hour_test.py'
    else:
        recipients = args.recipients
    # Respond only to 'CBC' events that have a change of EMBRIGHT.
    if params['Group'] != 'CBC' or (float(params['BNS']) + float(params['NSBH']) < 0.1):
        return

    # Print and save some parameters.
    interesting_parameters = ['GraceID', 'AlertType', 'EventPage',
    'Instruments', 'FAR', 'skymap_fits',
    'BNS', 'NSBH', 'BBH', 'MassGap',
    'Terrestrial', 'HasNS', 'HasRemnant', 'role']
    with open('./'+params['GraceID']+'.dat','w') as f:
        for key in interesting_parameters:
            print(key, '=', params[key])
            f.write('%s = %s\n'%(key,params[key]))

    https, skymapfits = os.path.split(params['skymap_fits'])
    skymap_local,_ = urllib.request.urlretrieve(params['skymap_fits'],'./'+params['GraceID']+'.fits')
    #skymap_local = skymapfits
    #skymap_local,_ = urllib.urlretrieve(params['skymap_fits'],'./'+params['GraceID']+'.fits')
    params['skymap_fits'] = skymap_local
    print(params['skymap_fits'])
    print("Skymap downloaded.")

    if 'skymap_fits' in params:
        process_fits()
    return
def process_fits():
        global params
        # Read the HEALPix sky map and the FITS header.
        skymap, header = hp.read_map(params['skymap_fits'],
                                     h=True, verbose=False)
        # Print some values from the FITS header.
        header = dict(header)
        header['GraceID'] = params['GraceID']
        with open('./'+params['GraceID']+'.dat','a') as f:
            print('Distance =', header['DISTMEAN'], '+/-', header['DISTSTD'])
            f.write('Distance = {} +/- {}\n'.format(header['DISTMEAN'], header['DISTSTD']))
        time = Time.now()
        prob, probfull, timetill90, m = prob_observable(skymap, header, time, plot=plotting)
        params['skymap_array'] = m
        if timetill90 ==-99:
            print("HET can't observe the source.")
            return
        else:
            print("Source has a {:.1f}% chance of being observable now.".format(
                 int(round(100 * prob))))
            print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(
                int(round(100 * probfull))))
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(
                timetill90))
            if args.send_notification:
                send_notifications(params,timetill90,text=True,email=False)
            for catalog in args.cat:
                get_galaxies.write_catalog(params,catalog)
                mincontour = get_LST.get_LST(targf = 'galaxies%s_%s.dat'%(catalog,params['GraceID']))
                make_phaseii.make_phaseii('LSTs_{}.out'.format(params['GraceID']))
            if args.send_notification:
                send_notifications(params,timetill90)

def send_notifications(params,timetill90,text=False,email=True):
    if text:
        msjstring=''
        if args.test > 0:
            msjstring = 'TEST '
        msjstring += 'GW ALERT!: Time till 90% prob region is {:.1f} hours  '.format(timetill90)
        email_ip.SendText(msjstring, emails=[], recipients = recipients)
    if email:
        with open('./'+params['GraceID']+'.dat','r') as f:
            emailcontent = '### {} GW ALERT ###\n'.format(params['AlertType'])
            if args.test > 0:
                emailcontent += 'THIS IS A TEST \n'
            emailcontent += 'Time until 90% probability region: {:.1f} hours\n\n'.format(timetill90)
            emailcontent += 'CALL MAJO: +15125763501 \n'
            emailcontent += 'JOIN THE ZOOM LINK TO TALK TO THE MOUNTAIN: \n'
            emailcontent += 'https://zoom.us/j/715671547 \n'
            emailcontent += f.read()
            emailcontent += '\n\n'
            emailcontent += "If you happen to find the location of the source, please submit coordinates to GraceDB using submit_gracedb.py. \n"
            emailcontent += '######\n'
            emailcontent += "Alert created with DIAGNOSIS, for more information on the software and data products, please refer to the wiki: \n"
            emailcontent += "https://github.com/Majoburo/Diagnosis \n"
            email_ip.SendText(emailcontent,
                    plotfiles = [x.format(params['GraceID']) for x in ['LSTs_{}.pdf','MOLL_GWHET_{}.pdf']],
                    datafiles = [x.format(params['GraceID']) for x in ['LSTs_{}.out','{}.tsl']] ,
                    numbers = [],
                    recipients = recipients)

def main():
    # Listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    global args
    global params
    args = parseargs()
    if args.fits:
        skymap, header = hp.read_map(args.fits, h=True, verbose=False)
        time = Time.now()
        header = dict(header)
        header['GraceID'] = header['OBJECT']
        params = {'AlertType':'fits','skymap_fits':args.fits,'skymap_array':np.copy(skymap),'GraceID':header['OBJECT']}
        with open('./'+params['GraceID']+'.dat','w') as f:
            f.write('THIS ALERT DOES NOT CONTAIN RELEVANT PARAMETER INFORMATION SINCE IT WAS CREATED FROM A FITS FILE (NOT A GCN NOTICE).')
        process_fits()
        return

    if args.test > 1:
        import lxml.etree
        payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
        root = lxml.etree.fromstring(payload)
        process_gcn(payload, root)
    else:
        gcn.listen(handler=process_gcn, port=8099)
if __name__== "__main__":
   main()
