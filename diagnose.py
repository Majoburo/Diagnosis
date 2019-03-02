import email_ip, get_galaxies, get_LST, GCN_listen

def main():
    # Listen for GCNs until the program is interrupted
    # (killed or interrupted with control-C).
    #gcn.listen(handler=process_gcn, port=8099)
    #Uncomment this lines when testing:
    import lxml.etree
    payload = open('MS181101ab-1-Preliminary.xml', 'rb').read()
    root = lxml.etree.fromstring(payload)
    GCN_listen.process_gcn(payload, root)
    params = GCN_listen.params
    get_galaxies.write_catalog2MASS(params['skymap_fits'])
    #get_galaxies.write_catalogGLADE(params['skymap_fits'])
    get_LST.get_LST(targf = 'galaxies2MASS.dat')
if __name__== "__main__":
    main()

