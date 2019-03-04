from ligo.gracedb.rest import GraceDbBasic, HTTPError
import argparse

def parseargs():

    parser = argparse.ArgumentParser(description='Submit event to GraceDB')
    parser.add_argument('group', default = 'test', help = 'MOU group responsible')
    parser.add_argument('grace_id',default = 'EVENTNAME', help ='Identifier for the GW event')
    parser.add_argument('ralist', help = 'List of ra of centers of footprints (degrees)')
    parser.add_argument('declist', help = 'List of ra of centers of footprints (degrees)')
    parser.add_argument('rawithlist', help = 'List (or one for all) of footprint widths in ra (degrees)')
    parser.add_argument('decwithlist', help = 'List (or one for all) of footprint widths in dec (degrees)')
    parser.add_argument('durationlist', default = '20.0', help = 'List (or one for all) of exposure times in sec')
    parser.add_argument('starttimelist', help = 'List of beginnings of exposures (utc)', default = [
        '2015-05-31t12:45:00',
        '2015-05-31t12:49:00',
        '2015-05-31t12:53:00'])
    parser.add_argument('-c', dest='comment', default = ' ', help = 'Comments')

    args = parser.parse_args()

    return args

def submit_gracedb(graceid, group, ralist, rawidthlist,
        declist, decwidthlist, starttimelist, durationlist, comment):
    service = 'https://gracedb.ligo.org/apibasic/'
    # Instantiate the GraceDB client
    g = GraceDbBasic(service)
    # Write the EMObservation record to GraceDB
    g.writeemobservation(graceid, group, ralist, rawidthlist,
        declist, decwidthlist, starttimelist, durationlist, comment)
    if g.status == 201:       # 201 means 'Created'
        print('Success!')
    return

def main():
    a = parseargs()
    submit_gracedb(a.graceid, a.group, a.ralist, a.rawidthlist,
        a.declist, a.decwidthlist, a.starttimelist, a.durationlist, a.comment)
    return

if __name__=='__main__':
    main()
