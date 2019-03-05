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

plotting = True
params = {}
args = 0

def parseargs():

    parser = argparse.ArgumentParser(description='Receive and parse GCN alerts, alert observers and create observing tools.')
    parser.add_argument('recipients', help = 'Specify python file with list of recipients.')
    parser.add_argument('-cat', dest='cat', default = 'GLADE', help='Specify which catalog to use: 2MASS or GLADE')
    parser.add_argument('-t','--test', dest='test', action = 'store_true', help='DEFAULT:FALSE. Run a test. Will not query gcn but use a file stored in the git repo.')
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
    # Respond only to 'test' events.
    # VERY IMPORTANT! Replace with the following code
    # to respond to only real 'observation' events.
    if args.test:
        if root.attrib['role'] != 'test':
            return
    elif root.attrib['role'] != 'observation':
            return
    #print('GW ALERT!')

    # Read all of the VOEvent parameters from the "What" section.
    global params
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Respond only to 'CBC' events that have a change of EMBRIGHT.
    if params['Group'] != 'CBC' or (float(params['BNS']) + float(params['NSBH']) < 0.5):
        return

    # Print and save all parameters.
    with open('./'+params['GraceID']+'.dat','w') as f:
        for key, value in params.items():
            print(key, '=', value)
            f.write('%s = %s\n'%(key,value))

    https, skymapfits = os.path.split(params['skymap_fits'])
    skymap_local,_ = urllib.request.urlretrieve(params['skymap_fits'],'./'+params['GraceID']+'.fits')
    #skymap_local,_ = urllib.urlretrieve(params['skymap_fits'],'./'+params['GraceID']+'.fits')
    params['skymap_fits'] = skymap_local
    print("Skymap downloaded.")

    if 'skymap_fits' in params:
        # Read the HEALPix sky map and the FITS header.
        skymap, header = hp.read_map(params['skymap_fits'],
                                     h=True, verbose=False)
        header = dict(header)
        header['GraceID'] = params['GraceID']
        # Print some values from the FITS header.
        print('Distance =', header['DISTMEAN'], '+/-', header['DISTSTD'])
        with open('./'+params['GraceID']+'.dat','a') as f:
            f.write('Distance = {} +/- {}\n'.format(header['DISTMEAN'], header['DISTSTD']))
        time = Time.now()

        prob, probfull, timetill90 = prob_observable(skymap, header, time, plot=plotting)
        if timetill90 ==-99:
            print("HET can't observe the source.")
        else:
            print("Source has a {:.1f}% chance of being observable now.".format(
                    int(round(100 * prob))))
            print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(
                    int(round(100 * probfull))))
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(
                    timetill90))
            if args.send_notification:
                email_ip.SendText('GW ALERT! Time till 90% prob region is {:.1f} hours  '.format(timetill90),emails=[],recipients = args.recipients)
        for catalog in args.cat:
            get_galaxies.write_catalog(params,catalog)
            get_LST.get_LST(targf = 'galaxies%s_%s.dat'%(catalog,params['GraceID']))
            make_phaseii.make_phaseii('LSTs_{}.out'.format(params['GraceID']))
        with open('./'+params['GraceID']+'.dat','r') as f:
            emailcontent = '### {} GW ALERT ###\n'.format(params['AlertType'])
            emailcontent += 'Time until 90% probability region: {:.1f} hours\n\n'.format(timetill90)
            emailcontent += f.read()
            emailcontent += '\n\n'
            emailcontent += "If you happen to find the location of the source, please submit coordinates to GraceDB using submit_gracedb.py. \n"
            emailcontent += '######\n'
            emailcontent += "Alert created with DIAGNOSIS, for more information on the software and data products, please refer to the wiki: \n"
            emailcontent += "https://github.com/Majoburo/Diagnosis \n"
        if args.send_notification:
            email_ip.SendText(emailcontent,
                    plotfiles = [x.format(params['GraceID']) for x in ['LSTs_{}.pdf','MOLL_GWHET_{}.pdf']],
                    datafiles = [x.format(params['GraceID']) for x in ['LSTs_{}.out','{}.tsl']] ,
                    numbers = [],
                    recipients = args.recipients)
    return

def main():
    # Listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    global args
    args = parseargs()
    if args.test:
        import lxml.etree
        payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
        root = lxml.etree.fromstring(payload)
        process_gcn(payload, root)
    else:
        gcn.listen(handler=process_gcn, port=8099)
if __name__== "__main__":
    main()
