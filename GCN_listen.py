import gcn
import healpy as hp
from prob_obs import prob_observable
from astropy.time import Time
from astropy.utils.data import download_file
import email_ip, get_galaxies, get_LST
import os


plotting = True
params = {}

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
    # if root.attrib['role'] != 'observation':
    #    return
    if root.attrib['role'] != 'test':
        return

    # Read all of the VOEvent parameters from the "What" section.
    global params
    params = {elem.attrib['name']:
              elem.attrib['value']
              for elem in root.iterfind('.//Param')}

    # Respond only to 'CBC' events. Change 'CBC' to "Burst'
    # to respond to only unmodeled burst events.
    if params['Group'] != 'CBC':
        return

    # Print all parameters.
    for key, value in params.items():
        print(key, '=', value)
    https, skymapfits = os.path.split(params['skymap_fits'])
    params['skymap'] = skymapfits
    download_file(params['skymap_fits'])
    print("Skymap downloaded.")

    if 'skymap_fits' in params:
        # Read the HEALPix sky map and the FITS header.
        skymap, header = hp.read_map(skymapfits,
                                     h=True, verbose=False)
        header = dict(header)
        # Print some values from the FITS header.
        print('Distance =', header['DISTMEAN'], '+/-', header['DISTSTD'])
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
            email_ip.SendText('GW ALERT! 90 % prob region accesible in {:.1f}'.format(
                    timetill90))

    return

def main():
    # Listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    #gcn.listen(handler=process_gcn, port=8099)
    import lxml.etree
    payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
    root = lxml.etree.fromstring(payload)
    process_gcn(payload, root)
    
    get_galaxies.write_catalog2MASS(params['skymap'])
    #get_galaxies.write_catalogGLADE(params['skymap'])
    get_LST.get_LST(targf = 'galaxies2MASS.dat')
if __name__== "__main__":
    main()
