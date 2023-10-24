# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 09:46:19 2023

@author: finnl
"""

from astropy.io import fits
from astropy import wcs
from astropy.wcs import WCS
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astroquery.sdss import SDSS
from astropy import constants as const
from astropy import units as u
from astroquery.simbad import Simbad

import matplotlib.pyplot as plt
import numpy as np

ra = fits.open('C:/Users/finnl/Downloads/red_antennae.fits')

dat = ra[0].data
#prints all header info
print(ra[0].header)

#prints date
print(ra[0].header['DATE'])

# xdata = ra[0].data[]
# ydata = ra[0].data[]

plt.imshow(np.log(ra[0].data))

#brightess pixel indices
print(np.argwhere(dat == np.amax(dat)))

plt.scatter(607,474,c='k')
plt.axhline(y=474,c='r')
plt.axvline(x=607,c='r')

#%%
m1 = fits.open('C:/Users/finnl/Downloads/m51_optical_R.fits')[0]
           
dat = m1.data
dat[:, :, 2] = np.minimum(data[:, :, 2] + 50, 255) 
wcs = WCS(m1.header)
ax = plt.subplot(projection = wcs)
ax.imshow(dat)


#sources within 10 arcmins of galactic centre
r = 2*u.arcmin
result = SDSS.query_region("m51",radius=r)
print(result)

#%%

def bol_lum(R,T):
    return (T**4*u.K)*4*np.pi*(R**2*u.m)*(const.sigma_sb)
    
mydat = np.loadtxt('my_data_file.dat')

def radius(m,g):
    return np.sqrt((const.G * (m*1.9e30*u.m)/(g/100)*u.m/u.s))

radii=[]
for i in range(len(mydat)):
    radii.append(radius(mydat[i,0],mydat[i,1]))

lums=[]
temp = []
for i in range(len(mydat)):
    lums.append(bol_lum(radii[i],mydat[i,2]))
    temp.append(mydat[i,2]*u.K)
    
plt.plot(temp*u.dimensionless_unscaled,lums*u.dimensionless_unscaled)

