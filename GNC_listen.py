import gcn
import healpy as hp
import numpy as np
from prob_obs import prob_observable
from astropy.time import Time, TimeDelta

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

    #skymap_url = root.find(".//Param[@name='skymap_fits']").attrib['value']

    if 'skymap_fits' in params:
        # Read the HEALPix sky map and the FITS header.
        skymap, header = hp.read_map(params['skymap_fits'],
                                     h=True, verbose=False)
        header = dict(header)
        # Print some values from the FITS header.
        print('Distance =', header['DISTMEAN'], '+/-', header['DISTSTD'])
        time = Time.now()
        prob, probfull = prob_observable(skymap, header, time, plot=plotting)
        print("Source has a {:.1f}% chance of being observable now.".format(
                int(round(100 * prob))))
        print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(
                int(round(100 * probfull))))

    return
plotting = False
def main():
    # Listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    gcn.listen(handler=process_gcn, port=8099)

if __name__== "__main__":
    main()
