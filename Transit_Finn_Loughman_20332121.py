# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 14:02:53 2023

@author: finnl
"""
#ASTRO COMPUTATIONAL LAB ASSIGNMENT 1 - Transit Function
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


print('Finn Loughman: 20332121')

#Loading data file
tc = np.load('JSAstroLab2023_transit_data_20332121.npz')

#Defining data columns
time_transit = tc['time']
flux_transit = tc['flux']


#Period of planet
P = 1.14176084

#parameter guesses
P0 = [np.mean(time_transit),10,0.158,np.pi/2,1]


#Transit function
def transit_model(t,T0,a,rho,i,f_oot):
    
    
    #Orbital Phase
    phi = ((2*np.pi)/P) * (t-T0)

   
    #Z in terms of phi
    z = a * (np.sqrt(((np.sin(phi))**2) + ((np.cos(i) * (np.cos(phi)))**2)))
    
    
    #Initial flux array of 1
    flux = np.ones(len(time_transit))
    
    #Index for when planet is transitting
    index = (z > (1-rho)) * (z <= (1+rho))
 
    
    #Kappa definitions
    k0 = np.arccos(((rho**2) + (z[index]**2) - 1)/((2 * rho) * (z[index])))
    k1 = np.arccos((1 - (rho**2) + (z[index]**2))/(2 * (z[index])))
    
    #Defining the middle piecewise equation 
    f = 1 - ((1/np.pi) * ((((rho**2)*k0))+k1 - np.sqrt((4*(z[index]**2) - (1+(z[index]**2)-(rho**2))**2)/(4))))
   
    #Conditional indexing to form piecewise transit function
    flux[z > (1+rho)] = 1
    flux[index] = f
    flux[z <= (1-rho)] = 1 - rho**2
    
    

    return flux * f_oot



#Applying curve fit
popt_transit, pcov_transit = opt.curve_fit(transit_model, time_transit,flux_transit,P0,bounds=(0,np.inf))


#Priniting optimised parameters and uncertainties
print ('Optimised transit parameters',popt_transit)
print('Variance in the parameters',np.sqrt(np.diag(pcov_transit)))

#Uncertainties in parameters
uncertainty = np.sqrt(np.diag(pcov_transit))
error = transit_model(time_transit,*popt_transit + uncertainty) - transit_model(time_transit,*popt_transit - uncertainty)

#Plotting actual data 
plt.plot(time_transit,flux_transit,c='b',label='Transit Data',linewidth=0.7)
plt.title('Normalised Flux vs. Time')
plt.xlabel('Time(BJD days$_{exo}$)')
plt.ylabel('Flux$_{normalised}$')
plt.tick_params(axis='x',direction = 'in')
plt.tick_params(axis='y',direction = 'in')

#Plotting fitted data
plt.errorbar(time_transit, transit_model(time_transit,*popt_transit),yerr= error,fmt='o',markersize=1.5,ecolor='k',capsize =2,elinewidth=0.5,c='r',label='Fitted Data')
plt.grid()
plt.legend(loc=3)
plt.show()

#%%

