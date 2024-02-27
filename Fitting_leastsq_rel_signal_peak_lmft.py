# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:03:01 2023

@author: P.E. Trevisanutto
"""

import matplotlib.pyplot as plt
import numpy as np
#from astropy import units as u
#from astropy.modeling import models
#from specutils.spectra import Spectrum1D
#from specutils.fitting import fit_lines
#import scipy.integrate as it
#import scipy.special as special
#import scipy.io
import math as math
#from numpy import pi
#import pandas as pd
#import numpy as np
#import colorama
#from colorama import Fore, Style
# first part with least squares
from scipy.optimize import curve_fit, leastsq, least_squares
from scipy.signal import find_peaks,argrelextrema, peak_widths,find_peaks_cwt
from lmfit.models import LorentzianModel, QuadraticModel
#from scipy.signal import chirp, find_peaks
#from detecta import detect_peaks
#import cupy as cp

# second part about ODR
#from scipy.odr import ODR, Model, Data, RealData

# style and notebook integration of the plots
#import seaborn as sns
#%matplotlib inline
#sns.set(style="whitegrid")
Xeas = []
Yeas = []
values =[]
x_array=[]
y_array_2gauss=[]
parameter0 = []
temp_file= []
fhwm=[]
Inte=[]
gamma=[]

File_parameter='parameter.dat'
generalWidth = 1
restart= False
#restart= True
Nist =False
#Nist =True
###########---------------------------------------------------
def add_peak(prefix, center, amplitude=0.005, sigma=0.05):
    peak = LorentzianModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center)
    pars[prefix + 'amplitude'].set(amplitude)
    pars[prefix + 'sigma'].set(sigma, min=0)
    return peak, pars

#-------------------------------------------------------------
def remove_exponent(value):
    #valuecomma = value.replace(",","")

    decial = value.split('e')

    exp = decial[1].replace("0","")
    
    if(exp=="+"):

        ret_val=float(decial[0])

    else:
        ret_val = float(decial[0])*(10**(int(exp)))



    return ret_val
####--------------------------------------------
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    #####
def takey(x0,x,y):
    dx_1=abs(x[0]-x[len(x)-1])
    for i in range(0,len(x)):
        
        dx= abs(float(x0)-x[i])
        
        
        if dx<dx_1:
            dx_1=dx
            j_value=i
            
    return y[j_value]            
####--------------------------------------------
def save_all(x,y,name):

    
    output = '\n'.join('\t'.join(map(str,row)) for row in zip(x,y))
    
    with open(name+".dat", 'w') as f:
        f.write(output)
        f.write("\n")
    
    f.close()
#-----------------------------------------
def From_str2float(str):
    rvect_f =[]
    for item in str:
        rvect_f.append(float(item))
        
    rv = np.asarray(rvect_f)
    return rv
######-------------------------------------------
def lorentzian( x, x0, a, gam ):
    return abs(a* gam*x) / (( gam*x)**2 + ( x**2 - x0**2 )**2)
######-------------------------------------------
def lorentzian_Nist( x, x0, a, gam, B):
    if(x0 < 300):
        F =1
    else:
        F = 1.0/(1.0+np.exp(-(x-B)))
        
    y = abs(a * gam)*x*F / (( gam*x)**2 + ( x**2 - x0**2 )**2)
    return y

######-------------------------------------------
def multi_lorentz( x, params ):
    off = params[0]
    #print('off '+str(off))
    paramsRest = params[1:]
    assert not ( len( paramsRest ) % 3 )
    return off + sum( [ lorentzian( x, *paramsRest[ i : i+3 ] ) for i in range( 0, len( paramsRest ), 3 ) ] )
######-------------------------------------------
def multi_lorentz_Nist( x, params ):
    off = params[0]
    #print('off '+str(off))
    paramsRest = params[1:]
    assert not ( len( paramsRest ) % 4 )
    return off + sum( [ lorentzian_Nist( x, *paramsRest[ i : i+4 ] ) for i in range( 0, len( paramsRest ), 4 ) ] )


######-------------------------------------------
def res_multi_lorentz( params, xData, yData ):
    
    diff = [ multi_lorentz( x, params ) - y for x, y in zip( xData, yData ) ]
    return diff
######-------------------------------------------
def res_multi_lorentz_Nist( params, xData, yData ):
    
    diff = [ multi_lorentz_Nist( x, params ) - y for x, y in zip( xData, yData ) ]
    return diff

######-------------------------------------------
#######################
#print (Fore.MAGENTA +"Festina lente")
print ("Festina lente",flush=True)
#import os
#print(os.environ.get('CUDA_PATH'))
#for line in open('mediated_1000_tdlda.dat.brd', mode='r'):
#for line in open('Nist_ELF_CsPbBr3_4kev.dat', mode='r'):
for line in open('Nist_ELF_CsPbBr3_0_2_4kev.dat', mode='r'):
    values=line.split()
    if (is_number(values[0]) == False):
        pass
    else:
        Xeas.append(values[0])
        Yeas.append(values[1]) 
    values.clear()    

yData = From_str2float(Yeas)

#yData = yData / max(yData)
xData = From_str2float(Xeas)


yDataLoc = yData    

if (restart == True):
    
    for line in open(File_parameter, mode='r'):
        
       for s in line.split():
           parameter0.append(float(s))
           
    counter = len(parameter0)
    
    counter = int((counter - 1)/3) 
    
    print(parameter0,flush=True)
    guess= parameter0
   
    #ax.loglog( xData, testData,'-')
else:
    guess = [min(yData)]


testData = [ multi_lorentz(x, guess ) for x in xData ]
fig = plt.figure()
ax = fig.add_subplot( 1, 1, 1 )
ax.loglog( xData, yData )
ax.loglog( xData, testData,'-')
save_all( xData, testData, 'fitted_ELF.dat')
plt.title('Initial guess')
plt.show()



peaks, properties = find_peaks( yDataLoc,width=0.5,height=0.00000000001 )
print(peaks)
#fhwm =peak_widths(yDataLoc, peaks,rel_height=0.5 )
#print (len(peaks))
#print (peaks)
    #pk = argrelextrema(yDataLoc,np.greater)
#extract peak heights and fwhm 
#results_h = peak_widths(yDataLoc, peaks,rel_height=1 )
#print("here we are:")
#print(results_h[1])
#print (results_half[0])
#print("Number of the found peaks: " +str(len(peaks)))
#print (peaks)
#print("Properties: ")
#print (properties)
    #print (properties)
I = properties [ 'peak_heights' ]
fwhm =properties["width_heights"]

for i in range( len( peaks) ): 
   

   guess.append( xData[peaks[i]] )
   Inte.append(I[i])
   gamma.append(fwhm[i])
#   guess.append(I[i])
#   guess.append(fwhm[i]/xData[peaks[i]])

model = QuadraticModel(prefix='bkg_')
params = model.make_params(a=0, b=0, c=0)

rough_peak_positions = From_str2float(guess)
Amp = From_str2float(Inte)
gamm =From_str2float(gamma)
print (rough_peak_positions)
for i in range(len(peaks)):
    peak, pars = add_peak('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
    model = model + peak
    params.update(pars)

init = model.eval(params, x=xData)
result = model.fit(yData, params, x=xData)
comps = result.eval_components()

print(result.fit_report(min_correl=0.5))

plt.loglog(xData,yData, label='data')
plt.loglog(xData, result.best_fit, label='best fit')
plt.show()

for name, comp in comps.items():
    plt.loglog(Xeas, comp, '--')
plt.legend(loc='upper right')
plt.show()


with open("gnufile_fit", 'w') as f:

  f.write("f(x)=")
  for i in range( len( peaks) ):
      f.write("+abs("+ str(I[i])+"*"+str(fhwm[i])+"*x) / (( "+str(fhwm[i])+"*x)**2 + ( x**2 - "+ str(guess[i])+"**2 )**2)")
    
    
print ("Sutor ne ultra crepidam",flush=True)
