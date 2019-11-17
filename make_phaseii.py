import numpy as np
import sys
from astropy import units as u
from astropy.coordinates import Angle
def make_phaseii(lstfile):
    common = {
                'PROGRAM':'HET19-3-400',
                'VIFU':'047',
                'EXP':'360',
                'NUMEXP':'3',
                'EQUINOX':'2000.0',
                'INSTRUMENT':'VIRUS',
                'MAG':'22',
                'SKYBRIGHT':'18.0',
                'SEEING': '3.0',
                'SKYTRANS': 'S',
                'SKYCALS': 'Y',
                'FLUX': 'Y',
                'PRI':'0',
                'SETUPMETHOD':'DirectGuider',
                'COMMENT':'"Usual Dither, look for new object in target IFU"',
                }
    GraceID = lstfile.split('_')[1].split('.')[0]
    with open(GraceID+'.tsl','w') as f:
        f.write('COMMON\n')
        for key,value in common.items():
            f.write('\t{}\t{}\n'.format(key,value))
        f.write('TRACK_LIST\n')
        f.write(' OBJECT\tRA\tDEC\tPIPRI\n')
        targets = np.loadtxt(lstfile,skiprows=1,dtype=np.str)
        targets = np.atleast_2d(targets)
        dec = Angle(targets[:,2],u.degree).dms#.to_string(unit=u.degree,sep=':')
        ra = Angle(targets[:,1],u.degree).hms#.to_string(unit=u.hour,sep=':')
        dec = ["{:+03.0f}:{:02.0f}:{:05.2f}".format(dd[0],abs(float(dd[1])),abs(float(dd[2]))) for dd in zip(dec[0],dec[1],dec[2])]
        ra = ["{:02.0f}:{:02.0f}:{:06.3f}".format(rr[0],rr[1],rr[2]) for rr in zip(ra[0],ra[1],ra[2])]
        for i,target in enumerate(targets):
            if float(target[2])>0:
                target[2]='+'+target[2]
            f.write('GW{}\t{}\t{}\t{}\n'.format(target[1]+target[2],ra[i],dec[i],str(int(target[0]))))
def main():
    make_phaseii(sys.argv[1])
if __name__=='__main__':
    main()
