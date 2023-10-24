# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 15:42:18 2023

@author: finnl
"""
import numpy as np
from matplotlib import pyplot as plt 
from scipy import optimize
from astropy.cosmology import FlatLambdaCDM
import scipy.stats as stats

print('Finn Loughman 20332121')
#Loading spectra
spec1=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn1997bp_spec.ascii',unpack=True)
spec2=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn1997E_spec.ascii',unpack=True)
spec3=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn1998es_spec.ascii',unpack=True)
spec4=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn1999aa_spec.ascii',unpack=True)
spec5=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn1999dq_spec.ascii',unpack=True)
spec6=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn2000cn_spec.ascii',unpack=True)
spec7=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_spectra\sn2000dk_spec.ascii',unpack=True)


spectrumlist=(spec1,spec2,spec3,spec4,spec5,spec6,spec7)

#Empty lists for wavelengths and fluxes
waves=[]
fluxes =[]


for j in range(len(spectrumlist)):
    #to get the x and y axis from each spectrum:
    wave1=list(spectrumlist[j][0])
    flux1=list(spectrumlist[j][1])
    waves.append(wave1)
    fluxes.append(flux1)
    plt.plot(wave1,flux1)
    plt.show()
#%%

#Defining Gaussian function
def gauss(x, H, A, x0, sigma):
        return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

# Rest wavelengths of most prominent emission lines
H=6563
NaI=5896


# The emission line used in each spectrum, first index = spec1 etc...
transition=(H,H,H,H,NaI,H,H)
# Initial estimates for  mean and amplitude in gaussian fit
means=(6620,6655,6620,6650,5950,6715,6675)
amps=(2e-16,2e-16,2e-15,0.3,1,2e-16,2e-16)

# The max and min create the interval in which the peak is seen
mins=(6200,6200,6500,6400,5600,6400,6200)
maxes=(7200,7200,7400,7400,6400,7400,7200)

# Width parameter of the gaussian
width=20

# Fitting the Gaussian to the main emission peak in each spectru
for j in range(len(spectrumlist)):
    # Separating wavelength and flux columns
    wave1=list(spectrumlist[j][0])
    flux1=list(spectrumlist[j][1])
    waverange=[]
    
    # Initial estimates:
    mean=means[j]
    amp=amps[j]
    
    # Interval bounds for the specific peak in question
    bottom = mins[j]
    top = maxes[j]
    
    for i in range(len(wave1)):
        # Finds the wavelength range in the interval:
        if bottom<wave1[i]<top:
            waverange.append(wave1[i])
    # Finds the index of highest and lowest wavelngth in range
    wavemin=wave1.index(min(waverange))
    wavemax=wave1.index(max(waverange))
    
    # Flux range correspinding to wavelngth range:
    fluxrange=flux1[wavemin:wavemax+1]
    # Finds in the index of the peak:
    ind=fluxrange.index(max(fluxrange))
    
    # Creates a portion of the data to plot on either side of the peak:
    x_gauss=waverange[ind-width:ind+width]
    y_gauss=fluxrange[ind-width:ind+width]
    plt.plot(x_gauss,y_gauss,label='Data')
    
    # Initial parameters for Gaussian fit:
    p0=(0.3,amp,mean,4)
    parameters, covariance = optimize.curve_fit(gauss, x_gauss, y_gauss,p0)
    
    # Plots Gaussian onto peak for each spectrum
    plt.plot((x_gauss),gauss(x_gauss,parameters[0],parameters[1],parameters[2],parameters[3]),label='Gaussian Fit',c='r')
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (erg/s/cm$^2$)')
    plt.title('Spectrum %s'% (1+j))
    plt.legend()
    plt.show()
    
#%%

# Empty lists
vels=[]
vels_err=[]


actualmean=[]
errors=[]
zvalues=[]
zvalueserrors=[]
fit_errlist=[]

# Finding the error in the Gaussian fit 
for j in range(7):
    #all same as before
    wave1=list(spectrumlist[j][0])
    flux1=list(spectrumlist[j][1])
    waverange=[]
    mean=means[j]
    amp=amps[j]
    bottom=mins[j]
    top=maxes[j]
    meanlist=[]
    for k in range(100):
        # Plots the Gaussian in the same way as above, except for serval different widths:        width=random.randint(10,30)
        for i in range(len(wave1)):
            if bottom<wave1[i]<top:
                waverange.append(wave1[i])
        wavemin=wave1.index(min(waverange))
        wavemax=wave1.index(max(waverange))
        fluxrange=flux1[wavemin:wavemax+1]
        ind=fluxrange.index(max(fluxrange))
        gaussx=waverange[ind-width:ind+width]
        gaussy=fluxrange[ind-width:ind+width]
        p0=(0.3,amp,mean,4)
        parameters, covariance = optimize.curve_fit(gauss, gaussx, gaussy,p0)
        # Error for each peak using square root of diagonal of covariance matrix
        fit_err=np.sqrt(np.diag(covariance))[2]
        fit_errlist.append(fit_err)
        
        meanlist.append(parameters[2])
        
    
    # Adding standard deviation of the 100 gaussian means to the average fit error
    error=(np.std(meanlist)/len(meanlist)+(sum(fit_errlist)/len(fit_errlist)))

    errors.append(error)
    
    actualmean.append(sum(meanlist)/len(meanlist))
    
    
    meanshift=(np.array(actualmean)[j]-transition[j])
    
    # Calculating redshift and redshift error for each spectrum
    zvalues.append(meanshift/transition[j])
    zvalueserrors.append(np.array(errors[j])/transition[j])
    


print('Redshift:',zvalues)
print('Error in redshift:',zvalueserrors)



#%%

#PART B

spectrum1=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn1997bp_UBVRI.dat',unpack=True)
spectrum2=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn1997E_UBVRI.dat',unpack=True)
spectrum3=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn1998es_UBVRI.dat',unpack=True)
spectrum4=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn1999aa_UBVRI.dat',unpack=True)
spectrum5=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn1999dq_UBVRI.dat',unpack=True)
spectrum6=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn2000cn_UBVRI.dat',unpack=True)
spectrum7=np.loadtxt(r'C:\Users\finnl\Downloads\Data_for_download_2223\Data_for_upload\SN_lightcurves\sn2000dk_UBVRI.dat',unpack=True)


max_brightness=[]
max_brightness_err=[]

m15B_err =[]
mday15 =[]

z = zvalues

speclist=(spectrum1,spectrum2,spectrum3,spectrum4,spectrum5,spectrum6,spectrum7)


for j in range(7):
    # Separates dates and apparent magnitudes
    date=list(speclist[j][0]*(1+z[j])) # 
    Bband=list(speclist[j][3])
    Bband_err=list(speclist[j][4])
    daterange=[]
    Bbandrange=[]
    Bband_new=[]
    Date_new=[]
    Bbandnew_err = []
    Bbandrange_err =[]
    
    
    for i in range(len(Bband)):
        # Disregards fluxes of 99.9 as filtered measurement was not taken
        if Bband[i]<99.9:
            Bband_new.append(Bband[i])
            Date_new.append(date[i])
            Bbandnew_err.append(Bband_err[i])
            
    # Finds max brightness (minimum app mag)
    maxindex=Bband_new.index(min(Bband_new))
    for i in range(len(Bband_new)):
        peakdate=(Date_new[maxindex])
        # Interval with 20 days of the peak brightness:
        if peakdate<=Date_new[i]<(peakdate+20):
            daterange.append(Date_new[i]-peakdate)
            Bbandrange.append(Bband_new[i])  
            Bbandrange_err.append(Bbandnew_err[i])
            
      
   
    plt.ylim(max(Bbandrange)+1,min(Bbandrange)-1)
    plt.xlabel('Days from peak app. mag.')
    plt.ylabel('Apparent Magnitude')
    plt.title('Apparent Magnitude vs Days of Spectrum %s'% (1+j))
    
    # Cubic polyfit of the scatter plot
    x=(np.polyfit(daterange,Bbandrange,3,full='False',cov='True'))

    
    # Finding the 15th day after peak brightness
    day15=(x[0][0]*(15**3)+(x[0][1]*(15**2))+(x[0][2]*15)+(x[0][3]))
    day0 = x[0][3]
   
 
    # Plotting day 15 and day 0
    plt.scatter(15,day15,color='g',label='Day 15')
    plt.scatter(0,day0,c='k',label='Day 0: Peak Brightness')
    # Y values of the fitting
    y=np.polyval(x[0],daterange)
    
    # Plotting the polyfit
    plt.plot(daterange,y,'--o',label='PolyFit',alpha=0.7)
    plt.legend()
    plt.show()
    
    # Finds brightness index of peak:
    indmax=(list(y).index(min(y)))
    # Finds corresponding date index:
    datemax=(daterange[indmax])
    
    # Appending values
    max_brightness.append(day0)
    mday15.append(day15)
     
    # Error in polyfit
    day0_err = x[3][3]
    
    #Apparent magnitude error
    max_brightness_err.append(day0_err + Bbandrange_err[0])
    
    #Luminosity decline error
    m15B_err.append(Bbandrange_err[0])
# Finds the luminosity decline width
m15B=(-np.array(max_brightness)+np.array(mday15))

print('Peak Brightness',max_brightness)
print('Peak Brightness error',max_brightness_err)

print('Day0-15 luminosity decline',m15B)
print('Day0-15 luminosity decline',m15B_err)



#%%

#PART C


# Converting redshifts to parsec
cosmo = FlatLambdaCDM(H0=70, Om0=0)
distances = cosmo.luminosity_distance(zvalues)


# Absolute magnitude conversion
absmag=np.array(max_brightness)-5*np.log10(np.array(distances)*1e6/10)

# Plots absolute magnitude vs luminosity decline width
plt.scatter(m15B,absmag,label='Data')

# Y-axis inversion
plt.ylim(max(absmag+1),min(absmag-1))

# Calculating a linear fit of the relation to find correction slope
linregress = stats.linregress(m15B,absmag)
plt.plot(m15B,linregress.slope*m15B + linregress.intercept,c='r',label='Linear Fit')
print('Correction slope =', linregress.slope)

plt.legend()
plt.ylabel('Absolute Magnitude')
plt.xlabel('Luminosity-Decline rate ($\Delta m_{15}(B)$)')
plt.title('Absolute Magnitude vs Luminosity-Decline rate')



#%%



# Correction for magnitude using linear fit slope found above
corrected_mag=max_brightness-(linregress.slope*(m15B-1.1))

# Plotting corrected and uncorrected Hubble Diagrams
logcz = [] 
for i in zvalues:
    logcz.append(np.log10(3e8* i))


# Linear Fit
linreg = stats.linregress(logcz,max_brightness)
fit = np.array(linreg.slope*np.array(logcz) + linreg.intercept)

# Plotting uncorrected HD
plt.scatter(logcz,max_brightness)
plt.plot(logcz, fit,c='r')
plt.title('Uncorrected Hubble Diagram')
plt.xlabel('log(cz)')
plt.ylabel('Apparent Magnitude')
plt.show()


# Linear fit corrected
linreg_corr = stats.linregress(logcz,corrected_mag)
corrected_fit = np.array(linreg_corr.slope*np.array(logcz) + linreg_corr.intercept)

# Plotting corrected HD
plt.scatter(logcz,corrected_mag)
plt.plot(logcz, corrected_fit ,c='r')
plt.title('Corrected Hubble Diagram')
plt.xlabel('log(cz)')
plt.ylabel('Apparent Magnitude')



chi2=(sum((fit - max_brightness)**2/(fit)))
chi2_corrected = (sum((corrected_fit - corrected_mag)**2/(corrected_fit)))
print(chi2)
print(chi2_corrected) #smaller value = better fit