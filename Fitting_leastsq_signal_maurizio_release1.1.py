# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:03:01 2023

@author: tjs39879
"""

import matplotlib.pyplot as plt
from lmfit import Model
from lmfit.models import QuadraticModel
# first part with least squares
from scipy.signal import find_peaks,argrelextrema, peak_widths,find_peaks_cwt
from lmfit.models import LorentzianModel, QuadraticModel,SplitLorentzianModel,ThermalDistributionModel
import classutility as cu



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
Sep_value= 120
#Nist =True
#
a=0
b=0
c=0

myutil = cu.utilities()

class fitting_with_lorentzians():
###########---------------------------------------------------
    def add_peak_Mau_F(self,prefix,  center, amplitude, sigma,x0F,A=1,kt=-1):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(myutil.Mau_Lorentizian, prefix=prefix)*ThermalDistributionModel(prefix=prefix,form='fermi')
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        pars[prefix+'amplitude'].set(A,vary=False)
        pars[prefix+'kt'].set(kt,vary=False)
        pars[prefix+'center'].set(x0F)
        return gmodel, pars

###########---------------------------------------------------
    def add_peak_Mau_h_lor(self,prefix, center, amplitude, sigma,sigma_r):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(myutil.Mau_Lorentizian_h, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma,gam_r=sigma_r)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        pars[prefix + 'gam_r'].set(sigma_r, min=0)
    
        return gmodel, pars
###########---------------------------------------------------
###########---------------------------------------------------
    def add_peak_Mau(self,prefix, center, amplitude,sigma):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(myutil.Mau_Lorentizian, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        return gmodel, pars
###########---------------------------------------------------
    def add_peak(self,prefix, center, amplitude, sigma):
        peak = LorentzianModel(prefix=prefix)
        pars = peak.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma, min=0)
        return peak, pars    
###########---------------------------------------------------
    def add_peak_fermi(self,prefix_Lor,prefix_fermi, center,amplitude, sigma, x0,A=1,kt=-1):
        peak = LorentzianModel(prefix=prefix_Lor)*ThermalDistributionModel(prefix=prefix_fermi,form='fermi')

        pars = peak.make_params()
        pars[prefix_Lor + 'center'].set(center)
        pars[prefix_Lor + 'amplitude'].set(amplitude)
        
        pars[prefix_Lor + 'sigma'].set(sigma, min=0)
        pars[prefix_fermi+'amplitude'].set(A)
        pars[prefix_fermi+'kt'].set(kt)
        pars[prefix_fermi+'center'].set(x0)
        return peak, pars
####--------------------------------------------
###########---------------------------------------------------
    def add_peak_h_lor(self,prefix_Lor,prefix_fermi, center,amplitude, sigma, sigma_r):
 
        peak = SplitLorentzianModel(prefix=prefix_Lor)
        pars = peak.make_params()
        pars[prefix_Lor + 'center'].set(center)
        pars[prefix_Lor + 'amplitude'].set(amplitude, min=0)
    
        pars[prefix_Lor + 'sigma'].set(sigma, min=0)
        pars[prefix_Lor + 'sigma_r'].set(sigma_r, min=0)
 
        return peak, pars
###########---------------------------------------------------

fit= fitting_with_lorentzians()

myutil.letsstart()
#import os
#print(os.environ.get('CUDA_PATH'))
#for line in open('mediated_1000_tdlda.dat.brd', mode='r'):
#for line in open('Nist_ELF_CsPbBr3_4kev.dat', mode='r'):
for line in open('Nist_ELF_CsPbBr3_0_2_4kev.dat', mode='r'):
    values=line.split()
    if (myutil.is_number(values[0]) == False):
        pass
    else:
        Xeas.append(values[0])
        Yeas.append(values[1]) 
    values.clear()    

yData = myutil.From_str2float(Yeas)

#yData = yData / max(yData)
xData = myutil.From_str2float(Xeas)

#---------------------------------------------------------------------------
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
#----------------------------------------------------------------------
peaks, properties = find_peaks( yDataLoc,width=0.1,height=0.000001)
fig = plt.figure()
ax = fig.add_subplot( 1, 1, 1 )
ax.loglog( xData, yData )
ax.loglog(xData[peaks], yData[peaks], "x")
ax.set_xlabel('$\omega$ (eV)')
ax.set_ylabel('ELF')

plt.title('Initial guess')
plt.show()
print("Number of guessed peaks: " +str(len(peaks)))

I = properties [ 'peak_heights' ]
fwhm =properties["width_heights"]

for i in range( len( peaks) ): 

   guess.append( xData[peaks[i]] )
   Inte.append(I[i])

   if (xData[peaks[i]] <= Sep_value):
       gamma.append(fwhm[i])
   else:
       gamma.append(fwhm[i]*10)

gmodel = QuadraticModel(prefix='bkg_')

params = gmodel.make_params()
params['bkg_a'].set(a)
params['bkg_b'].set(b)
params['bkg_c'].set(c)

rough_peak_positions = myutil.From_str2float(guess)
Amp = myutil.From_str2float(Inte)
gamm =myutil.From_str2float(gamma)

print (rough_peak_positions)


                                                 

for i in range(len(peaks)):

    if (rough_peak_positions[i]< Sep_value):
         peak, pars = fit.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
    else:
         peak, pars = fit.add_peak_Mau_h_lor('flz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i], gamm[i]*10)

    # supply initial values, attributes for bounds, etcetera:
        
    
    gmodel = gmodel + peak
    params.update(pars)

init = gmodel.eval(params, x=xData)
result = gmodel.fit(yData, params, x=xData)
comps = result.eval_components()

print(result.fit_report(),flush=True)

result.params.pretty_print()
fig = plt.figure()
ax = fig.add_subplot( 1, 1, 1 )
ax.loglog(xData,yData, label='data')
ax.loglog(xData, result.best_fit, label='best fit')
ax.set_xlabel('$\omega$ (eV)')
ax.set_ylabel('ELF')
plt.show()

fig = plt.figure()
ax.set_xlabel('$\omega$ (eV)')
ax.set_ylabel('ELF')
ax.loglog(xData,yData, label='data')
ax.loglog(xData, result.best_fit, label='best fit')

for name, comp in comps.items():

    plt.loglog(xData, comp, '--')

plt.show()

##x0,a,gam,A=1,X0,kt=-1
with open("parameter_lmft_quadratic.dat", 'w') as fq:
    for name, par in result.params.items():
        if name[0]=='b':
         
         fq.write(str(par.value)+' ')
    fq.write("\n") 
with open("parameter_lmft_lor.dat", 'w') as f:
    for name, par in result.params.items():
        if name[0]=='l':
         
         f.write(str(par.value)+' ')
    f.write("\n") 
with open("parameter_lmft_fermi.dat", 'w') as ff:
    for name, par in result.params.items():
        if name[0]=='f':
         
         ff.write(str(par.value)+' ')
    ff.write("\n") 
fq.close()
f.close()
ff.close()    


myutil.sutor()