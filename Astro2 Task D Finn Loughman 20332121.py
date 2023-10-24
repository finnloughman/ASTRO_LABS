# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 10:03:21 2023

@author: finnl
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.cosmology import FlatLambdaCDM, LambdaCDM



print('Finn Loughman 20332121')


#Definiitions
A = 0.141
Aerr = 0.006
B = 3.101
Berr = 0.075
MB = -19.05
MBerr = 0.02
delm = -0.07 
delmerr = 0.02 
H0 = 70

SN = pd.read_csv (r"C:\Users\finnl\Downloads\Data_for_download2\Data_for_download2\Allsupernovae.csv")
mb=SN.loc[:,"mb"]
mberr=SN.loc[:,"error_mb"]
x1=SN.loc[:,"X1"]
x1err=SN.loc[:,"error_X1"]
c=SN.loc[:,"c"]
cerr=SN.loc[:,"error_c"]
LogMst=SN.loc[:,"LogMst"]
LogMsterr=SN.loc[:,"error_Mst"]
zcmb=SN.loc[:,"zcmb"]

MB_new =[]
MBerr_new=[]

# Changes absolute magnitude depending if galaxy mass > 10**10 solar masses
for i in range(len(LogMst)):
    if LogMst[i] >= 10:
        MB_new.append(MB+delm)
        MBerr_new.append(MBerr+delmerr)
    else:
        MB_new.append(MB)
        MBerr_new.append(MBerr)
MB=np.array(MB)
MBerr=np.array(MBerr)

#Redshift vs distance modulus

distmodulus=(mb-MB_new)+(A*x1) -B*(c)

distmodulus_err = np.array(np.sqrt((1*mberr)**2+(-1*MBerr)**2+(x1*Aerr)**2+(A*x1err)**2+(-c*Berr)**2+(-B*cerr)**2))

# Average error in distance modulus
print(np.mean(distmodulus_err/distmodulus)*100)


plt.scatter(zcmb,distmodulus,s=6)
plt.xscale('log')
plt.title('Distance Modulus vs Redshift')
plt.ylabel('Distance Modulus')
plt.xlabel('Redshift')
plt.show()

#%%
distances=(3e5*(np.array(zcmb)))/(H0)*(1+(np.array(zcmb)/2)) #massless

cosmodark = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
cosmodist_dark=(cosmodark.luminosity_distance(np.array(zcmb)))
cosmo = FlatLambdaCDM(H0=70, Om0=1)
cosmodist=(cosmo.luminosity_distance(np.array(zcmb)))

dm = 5*np.log10(distances*(1e6)) - 5
cosmodark_dm = 5*np.log10(cosmodist_dark.value*(1e6)) -5
cosmo_dm = 5*np.log10(cosmodist.value*(1e6)) -5

plt.scatter(zcmb,dm,s=4,label='Massless')
plt.scatter(zcmb,cosmodark_dm,s=4,label='No Dark Energy')
plt.scatter(zcmb,cosmo_dm,s=4,label='70% Dark Energy')
plt.scatter(zcmb,distmodulus,label='Data',s=4,alpha=0.5)


plt.legend()
#plt.xscale('log')
plt.title('Distance Modulus vs Redshift')
plt.ylabel('Distance Modulus')
plt.xlabel('Redshift')
plt.show()

#%%
plt.scatter(zcmb,(dm-distmodulus),s=4,label='Massless model')
plt.scatter(zcmb,(cosmodark_dm -distmodulus),s=4,label='70% Dark Energy model')
plt.scatter(zcmb,(cosmo_dm- distmodulus),s=4,label='No Dark Energy model')
plt.legend()
plt.title('Residual vs Redshift')
plt.ylabel('Residuals')
plt.xlabel('Redshift')
plt.show()


# Chi squared statistics
chi2_massless=sum(((dm-distmodulus)**2)/distmodulus)
chi2_darkenergy=sum(((cosmodark_dm-distmodulus)**2)/distmodulus)
chi2_nodarkenergy=sum(((cosmo_dm-distmodulus)**2)/distmodulus)

print(chi2_massless)
print(chi2_darkenergy) # smallest, most accurate
print(chi2_nodarkenergy)

#Relative chi-squared values
print(chi2_massless/chi2_darkenergy)
print(chi2_nodarkenergy/chi2_darkenergy)
