import gcn
from tabulate import tabulate as tb
import healpy as hp
from prob_obs import prob_observable
from astropy.time import Time
import astropy
from astropy.utils.data import download_file
import email_ip, get_galaxies, get_LST, make_phaseii
import os, urllib
import argparse
import urllib.request
import lxml
import numpy as np
import json
import matplotlib.pyplot as plt
import lxml.etree

plotting = True
params = {}
interesting_parameters = ['GraceID', 'AlertType','time', 'EventPage',\
    'Instruments', 'FAR', 'skymap_fits',\
    'BNS', 'NSBH', 'BBH', 'MassGap',\
    'Terrestrial', 'HasNS', 'HasRemnant', 'Distance','role']

def parseargs():

    parser = argparse.ArgumentParser(description='Receive and parse GCN alerts, alert observers and create observing tools.')
    parser.add_argument('recipients', help = 'Specify python file with list of recipients.')
    parser.add_argument('-g','--graceid', dest='graceid', default = None, help='graceID of event')
    args = parser.parse_args()

    return args

# Function to call every time a GCN is received.
# Run only for notices of type
# LVC_PRELIMINARY, LVC_INITIAL, or LVC_UPDATE.
@gcn.handlers.include_notice_types(
    gcn.notice_types.LVC_PRELIMINARY,
    gcn.notice_types.LVC_INITIAL,
    gcn.notice_types.LVC_UPDATE,
    gcn.notice_types.LVC_RETRACTION)

def process_gcn(payload, root):
    # Timekeeping of latest notice for connection checkup.
    with open('time_last.txt', 'w') as f:
        f.write('%s'%Time.now().jd)

    # Gather information on the type of event
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}
    params['role'] = root.attrib['role']

    if params['AlertType'] == 'Retraction':
        email_ip.SendText(params['GraceID']+' Retracted.', AlertType = params['role'].upper(), recipients = args.recipients )
        return

    if params['role'] == 'test' and not args.graceid:
        return

    # Filter non 'CBC' events or those have no chance of EMBRIGHT.
    if params['Group'] != 'CBC' or (float(params['BNS']) + float(params['NSBH']) + float(params['MassGap']) < 0.1):
        return


    print("Processing fits.")
    process_fits(args.recipients, params)
    print("Fits processed.")

    return

def process_fits(recipients, params = None):

        # Read the HEALPix sky map and the FITS header.
        skymap,header = hp.read_map(params['skymap_fits'], h=True, verbose=False)

        # Print and save some values from the FITS header.
        header = dict(header)
        params['time'] = Time(header['DATE-OBS'],format='isot',scale='utc')
        time = Time.now()
        params['Distance'] = str(header['DISTMEAN']) + ' +/- ' + str(header['DISTSTD'])
        header['GraceID'] = params['GraceID']
        with open('./'+params['GraceID']+'.dat','w') as f:
            data_p = []
            tableheaders = ['PARAMETER','VALUE']
            for prm in interesting_parameters:
                if prm in list(params.keys()):
                    data_p.append((prm,params[prm]))
            print(tb(data_p,headers=tableheaders,tablefmt='fancy_grid'))
            f.write(tb(data_p,headers=tableheaders,tablefmt='html'))

        # Making a pie chart of the type of event for the email
        labels = ['BNS', 'NSBH', 'BBH', 'MassGap', 'Terrestrial']
        sizes = [float(params[label])*100 for label in labels]
        labels = ['%s (%.1f %%)'%(lab,pct) for lab,pct in zip(labels,sizes)]
        fig1, ax1 = plt.subplots()
        patches,texts = ax1.pie(sizes, startangle=90)
        ax1.legend(patches, labels, loc="best")
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig('piechart_%s.png'%params['GraceID'])

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
            send_notifications(params,timetill90,text=True,email=False)
            get_galaxies.write_catalog(params,'GLADE')
            mincontour = get_LST.get_LST(targf = 'galaxies%s_%s.dat'%('GLADE',params['GraceID']))
            make_phaseii.make_phaseii('LSTs_{}.out'.format(params['GraceID']))
            send_notifications(params,timetill90)

def send_notifications(params,timetill90,text=False,email=True):
    if text:
        msjstring=''
        msjstring = params['role'].upper() + ' '
        if timetill90 == -99:
            msjstring = 'OUT-OF-HET-REACH '
        msjstring += 'GW ALERT: Time till 90% prob region is {:.1f} hours  '.format(timetill90)
        email_ip.SendText(msjstring, AlertType = params['role'].upper(), emails=[], recipients = args.recipients )
    if email:
        with open('./'+params['GraceID']+'.dat','r') as f:
            emailcontent = '<style>\
                            td {width: 50% ;}\
                            tr {width: 100% ;}\
                            </style>'
            if timetill90 == -99:
                params["AlertType"] = 'OUT-OF-HET-REACH' #Check if this introduces problems with another part of the code.
            title = "{} GW ALERT\n <br/>".format(params["AlertType"])
            emailcontent += '<h1 style="text-align:center">%s</h1>\
                            <table>\
                            <tr>\
                            <td rowspan="2">'%title
            if params['role']=='test':
                emailcontent += 'THIS IS A TEST \n <br/>'
            emailcontent += 'Time until 90% probability region: {:.1f} hours\n <br/>\n <br/>'.format(timetill90)
            emailcontent += 'CALL MAJO: +15125763501 \n <br/>'
            emailcontent += 'JOIN THE ZOOM LINK TO TALK TO THE MOUNTAIN: \n <br/>'
            emailcontent += 'https://zoom.us/j/715671547 \n <br/>'
            emailcontent += f.read()
            emailcontent += '\n <br/>\n <br/>'
            emailcontent += "If you happen to find the location of the source, please submit coordinates to GraceDB using submit_gracedb.py. \n <br/>"
            emailcontent += '######\n <br/>'
            emailcontent += "Alert created with DIAGNOSIS, for more information on the software and data products, please refer to the wiki: \n <br/>"
            emailcontent += "https://github.com/Majoburo/Diagnosis \n <br/>"
            emailcontent += '</td>\
                            <td> <img width="100%%" src="cid:%s"> </td>\
                            </tr>\
                            <tr>\
                            <td> <img width="100%%" src="cid:%s"> </td>\
                            </tr>\
                            </table>'%(*[x.format(params['GraceID']) for x in ['MOLL_GWHET_{}.png','piechart_{}.png']],)
            email_ip.SendText(emailcontent,
                    AlertType = params['role'].upper(),
                    plotfiles = [x.format(params['GraceID']) for x in ['LSTs_{}.pdf','MOLL_GWHET_{}.png','piechart_{}.png']],
                    datafiles = [x.format(params['GraceID']) for x in ['LSTs_{}.out','{}.tsl']] ,
                    numbers = [],
                    recipients = args.recipients )

def main():
    # Process a fits file or listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    global args
    global params
    args = parseargs()
    if args.graceid:
        urllib.request.urlretrieve('https://gracedb.ligo.org/apiweb/superevents/'+args.graceid+'/files/',"index.html")
        with open('index.html') as f: a = json.load(f)
        xmlfiles = [key for key in a.keys() if key.endswith('xml')]
        latestxml = sorted(xmlfiles)[-1]
        urllib.request.urlretrieve('https://gracedb.ligo.org/apiweb/superevents/'+args.graceid+'/files/'+latestxml,latestxml)
        payload = open(latestxml, 'rb').read()
        root = lxml.etree.fromstring(payload)
        process_gcn(payload, root)
        return
    else:
        gcn.listen(handler = process_gcn, port=8099)
if __name__== "__main__":
   main()
