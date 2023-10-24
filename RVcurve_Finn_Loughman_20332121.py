# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 10:33:11 2023

@author: finnl
"""
#ASTRO COMPUTATIONAL LAB ASSIGNMENT 1 - RV Curve

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import os 

print('Finn Loughman: 20332121')

#Changing working directory
os.chdir('C:/Users/finnl/Downloads')


#Loading data file
rv = np.load('JSAstroLab2023_rv_data_20332121.npz')


#Defining the spectral template used for cross-correlation.
spec_template = rv['spec_template']


#For loop that compiles all the time data into one array
times=[]
for i in range(40):
    time = rv['time_BJD_{}'.format(i+1)]
    times.append(time)

#For loop that complies the spectrum arrays into one array of arrays.
spectra =[]
for i in range(40):
    spectrum = rv['spectrum_{}'.format(i+1)]
    spectra.append(spectrum)


#Empty lists used for appending RV and error values
radial_vels = []
err = []


#Iterates over the total spectra array of length 40, cross correlates each spectrum with the spectral template,
#then fits a Gaussian to the cross-correlation with the velocity spectrum to find the radial velocity for the spectrum.
#The 'mean (x0)' parameter is taken as the radial velocity.
for i in spectra:
    cc = np.correlate(i - i.mean(),spec_template-spec_template.mean(),mode='same')
    def gauss(x, H, A, x0, sigma):
        return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

    parameters, covariance = opt.curve_fit(gauss, rv['velocity'], cc)
    radial_vels.append(parameters[2])
    err.append((np.sqrt(np.diag(covariance))[2]))
      




#%%

#Fitting a sinusoidal function to the RV curve, to find the parameter values:
#a=range, 2pi/b = period, c = phase, d = offset
def model(x,A,b,c,d):
    return A*np.sin((b*np.array(x) + c)) + d


#Curve Fit
popt, pcov = opt.curve_fit(model, times, radial_vels,p0=[25,5.5,10,7500])

#Uncertainty in fitted data

plt.plot(times, model(times, *popt),'--o',markersize = '4',label= 'Curve Fit',c='r')

#Actual Data plotted
plt.plot(times,radial_vels,label = 'Data',c='b')
plt.errorbar(times,radial_vels,err,markersize ='4',fmt='o',elinewidth=0.7,ecolor='k',capsize=4,c='b')


plt.title('Stellar Radial Velocity vs. Time')
plt.xlabel('Time(BJD days$_{exo}$)')
plt.ylabel('Radial Velocity')
plt.legend(loc=4)
plt.grid()
plt.tick_params(axis='x',direction = 'in')
plt.tick_params(axis='y',direction = 'in')
plt.show()



#Optimised parameters and respective uncertainties
print(*popt)
print(np.sqrt(np.diag(pcov)))

#%%
#Phase folded RV curve - plots the phase vs. RV instead of the time vs.566 RV
P = 1.14176084
def phase(t):
    return (t-times[0])/P

phase_1period = []
for i in times:
    if phase(i)<=1:
        phase_1period.append(phase(i))
    else: phase_1period.append(phase(i)-1)    



#Phase folded fit
#a=range, 2pi/b = period, c = phase, d = offset
def phase_model(x,a,b,c,d):
    return a*np.sin((b*np.array(x) + c)) + d


#Sin Fit for phase folded RV curve
popt_pf, pcov_pf = opt.curve_fit(phase_model, phase_1period, radial_vels,p0=[3,5,4,7520])

#Plotting both phase folded light curve and the sinusoidal fit
plt.scatter(phase_1period,phase_model(phase_1period,*popt_pf),label='Fitted Data',c='r')
plt.scatter(phase_1period,radial_vels,linestyle='None',label='Phase Data',c='b',s=100,marker='1')
plt.title('Phase Folded RV curve')
plt.xlabel('Phase')
plt.ylabel('Radial Velocity')
plt.legend(loc=2)
plt.show()



#%%
#Plot of RV curve without systemic velocity
rv_zeroed =[]
for i in radial_vels:
    rv_zeroed.append(i-popt[3])
    
plt.plot(times,rv_zeroed,label = 'Data',c='b')
