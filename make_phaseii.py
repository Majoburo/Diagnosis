import numpy as np
import sys

def make_phaseii(lstfile):
    common = {
                'PROGRAM':'UT19-2-009',
                'VIFU':'046',
                'EXP':'120',
                'NUMEXP':'3',
                'EQUINOX':'2000.0',
                'INSTRUMENT':'VIRUS',
                'MAG':'22',
                'SKYBRIGHT':'18.0',
                'SEEING': '3.0',
                'SKYTRANS': '1',
                'SKYCALS': 'Y',
                'FLUX': 'Y',
                'PRI':'0',
                'SETUPMETHOD':'DirectGuider',
                'COMMENT':'Usual Dither, look for new object in target IFU',
                }
    GraceID = lstfile.split('_')[1].split('.')[0]
    with open(GraceID+'.tsl','w') as f:
        f.write('COMMON\n')
        for key,value in common.items():
            f.write('\t{}\t{}\n'.format(key,value))
        f.write('TRACK_LIST\n')
        f.write(' OBJECT\tRA\tDEC\tPIPRI\n')
        targets = np.loadtxt(lstfile,skiprows=1,dtype=np.str)
        for target in targets:
            f.write('GW{}\t{}\t{}\t{}\n'.format(target[1]+'_'+target[2],target[1],target[2],str(int(target[0])+1)))
def main():
    make_phaseii(sys.argv[1])
if __name__=='__main__':
    main()
