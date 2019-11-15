import sys
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_from_healpix
from astropy.io import fits
from astropy.wcs import WCS
filename_ligo = get_pkg_data_filename(sys.argv[1])
target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                 1000
NAXIS2  =                  800
CTYPE1  = 'RA---MOL'
CRPIX1  =                  500
CRVAL1  =                180.0
CDELT1  =                 -0.4
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--MOL'
CRPIX2  =                  400
CRVAL2  =                  0.0
CDELT2  =                  0.4
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')
array, footprint = reproject_from_healpix(filename_ligo, target_header)
from astropy.wcs import WCS
import matplotlib.pyplot as plt

ax = plt.subplot(1,1,1, projection=WCS(target_header))
#ax.imshow(array, vmin=0, vmax=1.e-8)
#ax.coords.grid(color='white')
ax.coords.frame.set_color('none')
import numpy as np

np.random.seed(19680801)

x = 30*np.random.randn(10000)
mu = x.mean()
median = np.median(x)
sigma = x.std()
textstr = '\n'.join((
    r'$\mu=%.2f$' % (mu, ),
    r'$\mathrm{median}=%.2f$' % (median, ),
    r'$\sigma=%.2f$' % (sigma, )))

# these are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.cm import get_cmap
cmap = get_cmap('terrain')

legend_elements = [Patch(facecolor=cmap(0.2),
                         label='HETDEX current coverage'),
                   Patch(facecolor=cmap(0.4),
                         label='90% area')]

# Create the figure
ax.legend(handles=legend_elements, loc='center')


plt.show()
