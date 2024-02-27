# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:42:14 2023
ALC_Sutor Project 
Main for the fitting procedure with Lorentians
Authors: PE Trevisanutto.
"""
import classutility as cu
import class_fit_spectrum_w_lorentzians as fit
from sys import argv
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from lmfit.models import QuadraticModel
###########---------------------------------------------------

##main ##main ##main ##main ##main ##main ##main

###########---------------------------------------------------
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
guess =[]
#File_parameter='parameter.dat'
generalWidth = 1

Sep_value= 120

a=0
b=0
c=0

myutil = cu.utilities()

lor = fit.fitting_with_lorentzians(myutil)

myutil.letsstart()
height =0.0000001


if (len(argv) == 2):
    filename =argv[1]
    
elif(len(argv) == 3):   
    filename =argv[1]
    height = float(argv[2])
    print ("Height: " + argv[2],flush=True)
else:
    print ("missing file/too many files")

#for line in open('Nist_ELF_CsPbBr3_4kev.dat', mode='r'):
#Nist_ELF_CsPbBr3_0_2_4kev.dat
print ("File: " + filename,flush=True)

    
Xeas,Yeas =myutil.taking_data(filename)


yData = myutil.From_str2float(Yeas)

#yData = yData / max(yData)
xData = myutil.From_str2float(Xeas)

#----------------------------------------------------------------------



peaks, properties = find_peaks( yData,width=0.1,height=height)

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
         peak, pars = lor.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
    else:
         peak, pars = lor.add_peak_Mau_h_lor('flz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i], gamm[i]*10)

    # supply initial values, attributes for bounds, etcetera:
        
    
    gmodel = gmodel + peak
    params.update(pars)



print("Calculating the best fit",flush=True)

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
