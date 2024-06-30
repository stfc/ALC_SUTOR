# -*- coding: utf-8 -*-

import classutility as cu
import class_fit_spectrum_w_lorentzians as fit
from sys import argv
import matplotlib.pyplot as plt

from lmfit.models import QuadraticModel,SplineModel,ExponentialModel,PolynomialModel,PowerLawModel
import numpy as np



###########---------------------------------------------------

##main ##main ##main ##main ##main ##main ##main

###########---------------------------------------------------
"""
Created on Wed Jul 26 11:42:14 2023
ALC_Sutor Project 
Main for the fitting procedure with Lorentians

Author: P.E. Trevisanutto


"""
PLOTTING = True
Guessing = False

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

cutoff1 =1.0*10**4
cutoff2 = 1*10**5
#cutoff2 = 4.0*10**5
cutoff3 = 4.0*10**5
Sep_value=100
offyougo=1

a=0
b=0
c=0
knot_xvals = np.array([1000, 3000, 5000, 7000, 8000, 9000, 10000, 11000, 12000,13000, 14000, 15000, 16000, 17000, 18000, 20000,21000,22000, 23000,24000, 25000, 27000, 28000,29000, 30000 ])
cApprox = 'L'
myutil = cu.utilities()

lor = fit.fitting_with_lorentzians(myutil)

myutil.letsstart()

height =10**(-9)

if (len(argv) == 2):
    filename =argv[1]

elif(len(argv) == 3):
    filename =argv[1]
    cApprox = argv[2]    
    
elif(len(argv) == 4):   
    filename =argv[1]
    cApprox = argv[2]
    height = float(argv[3])

elif(len(argv) == 5):   
    filename =argv[1]
    cApprox = argv[2]
    height = float(argv[3])
    Guessing = argv[4]
    
elif(len(argv) == 7):
    Sep_value=argv[5]
    cutoff1=argv[6]
    cutoff2=argv[7]
    
    print ("Height: " + argv[2],flush=True)
    
    print ("Approximation:", flush=True)
    if (cApprox == "F"):
        
        print ("Fermi", flush=True)
        
    elif(cApprox == "LN"):
        print ("Log-Normal", flush=True)
        
    elif(cApprox == "FA"):
        print ("Fano", flush=True)
    
    elif(cApprox == "L"):
         print ("Lorentzian", flush=True)
    else:
        print ("Not recognised: Fermi approximation", flush=True)        
else:
    print ("missing input/too many inputs", flush=True)
    print ("[ELF], [F,LN,FA], [Height]", flush=True)


print ("File: " + filename,flush=True)

    
Xeas,Yeas =myutil.taking_data(filename)


yData = myutil.From_str2float(Yeas)


xData = myutil.From_str2float(Xeas)


ix1 = myutil.index_of(xData, cutoff1)
ix11= myutil.index_of(xData, 5*cutoff1)
ix2 = myutil.index_of(xData, cutoff2)
ix3 = myutil.index_of(xData, cutoff3)


#----------------------------------------------------------------------

peaks, rough_peak_positions, Amp, gamm= myutil.taking_initialguess(xData,yData,height,Guessing)


if (PLOTTING==True):
     myutil.plotting_elf(xData,yData,peaks)



gmodel = QuadraticModel(prefix='bkg_')

#gmodel = PolynomialModel(degree=7,prefix='bkg_')
##Quadratic Model

#params = gmodel.make_params()
#params = gmodel.guess(yData, x=xData)

#params= gmodel1.guess(yData[ix1:], x=xData[ix1:])
print(f'parameter names: {gmodel.param_names}')
print(f'independent variables: {gmodel.independent_vars}')


params = gmodel.make_params()
params.update( gmodel.guess(yData[:ix1], x=xData[:ix1]))



exp1= PowerLawModel(prefix='pow3_')
par_exp= exp1.make_params()

par_exp.update(exp1.guess(yData[ix1:ix2], xData[ix1:ix2]))


exp2 = PowerLawModel(prefix='pow3_')
par_exp2 = exp2.make_params()

par_exp2.update(exp2.guess(yData[ix2:ix3], xData[ix2:ix3]))



for i in range(len(rough_peak_positions)):

    if (rough_peak_positions[i]< Sep_value):
         offyougo=1
         peak, pars = lor.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
         #peak, pars = lor.add_peak_Mau_h_lor('flz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i], gamm[i]*10)
    elif ((rough_peak_positions[i]>= Sep_value)and(rough_peak_positions[i]<= cutoff1)) : 
        offyougo=1
        if (cApprox=='FA'):
            print("Peak at: "+str(rough_peak_positions[i]))
            #peak, pars = lor.add_peak_FanoModel('fano%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i],1)
            peak, pars = lor.add_peak_FanoModel_Mau('fano%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i],1,1)
            
            #pars.update( peak.guess(yData, x=xData[ix1:]))
            
        elif (cApprox=='LN'):
            peak, pars = lor.add_peak_GaussianlogModel('gaussl%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i])
                    
        elif (cApprox=='F'):
            peak, pars = lor.add_peak_fermi('hlz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i],rough_peak_positions[i])
            #peak, pars = lor.add_peak_h_lor('hlz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i], gamm[i])
            #pars.update(peak.guess(yData[:ix1], xData[:ix1]))
        else:
            peak, pars = lor.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
            
    elif ((rough_peak_positions[i]>= cutoff1)and(rough_peak_positions[i]<= cutoff2)) : 
        offyougo=2
        if (cApprox=='FA'):
            print("Peak at zone2: "+str(rough_peak_positions[i]))
            #peak2, pars2 = lor.add_peak_FanoModel('fano%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i],1)
            peak2, pars2 = lor.add_peak_FanoModel_Mau('fano%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i],1,1)
            
        elif (cApprox=='LN'):
            peak2, pars2 = lor.add_peak_GaussianlogModel('gaussl%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i])
                    
        elif (cApprox=='F'):
            
            peak2, pars2 = lor.add_peak_fermi('hlz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i],rough_peak_positions[i])
            
        else:
            peak2, pars2 = lor.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
           
    else:
        offyougo=3
        if (cApprox=='FA'):
            print("Peak at zone2: "+str(rough_peak_positions[i]))
            
            peak3, pars3 = lor.add_peak_FanoModel('fano%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i],1)
            
        elif (cApprox=='LN'):
            peak3, pars3 = lor.add_peak_GaussianlogModel('gaussl%d_' % (i+1),rough_peak_positions[i], Amp[i], gamm[i])
        
        elif (cApprox=='F'):
           
            peak3, pars3 = lor.add_peak_fermi('hlz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i],rough_peak_positions[i])       
        else:
           
            peak3, pars3 = lor.add_peak_Mau('lz%d_' % (i+1), rough_peak_positions[i], Amp[i], gamm[i])
            
    # supply initial values, attributes for bounds, etcetera:
    if (offyougo==1):    
       
        gmodel = gmodel + peak
        params.update(pars)
    elif (offyougo==2):
        print (peak2)
        exp1 = exp1 + peak2
        par_exp.update(pars2)
        
    else:
        print (peak2)
        exp2 = exp2 + peak3
        par_exp2.update(pars3)
print ("\n")
print("Calculating the best fit",flush=True)

init = gmodel.eval(params, x=xData[:ix11])

result = gmodel.fit(yData[:ix1], params, x=xData[:ix1],nan_policy='omit' )
print ("done")

init1 =exp1.eval(par_exp, x=xData[ix1:ix2])

result1 = exp1.fit(yData[ix1:ix2], par_exp, x=xData[ix1:ix2],nan_policy='omit' )
print ("done")

init2 =exp2.eval(par_exp2, x=xData[ix2:ix3])

result2 = exp2.fit(yData[ix2:ix3], par_exp2, x=xData[ix2:ix3],nan_policy='omit' )
print ("done")

comps = result.eval_components()

print(result.fit_report(),flush=True)
print(result1.fit_report(),flush=True)
print(result2.fit_report(),flush=True)


result.params.pretty_print()

###da cambiare!!!
if(PLOTTING==True):
    myutil.resultbestfitplot(xData[:ix1],yData[:ix1], init, result.best_fit,"body")
    
    myutil.resultbestfitplot(xData[ix1:ix2],yData[ix1:ix2], init1, result1.best_fit,"tail")

   
    
    myutil.resultbestfitplot(xData[ix2:ix3],yData[ix2:ix3], init, result2.best_fit,"tail104")


#Printing converged data
with open("parameter_lmft_quadratic.dat", 'w') as fq:
    for name, par in result.params.items():
        if name[0]=='b':
         
         fq.write(str(par.value)+' ')
    fq.write("\n") 
fq.close()
with open("parameter_lmft_lor.dat", 'w') as f:
    for name, par in result.params.items():
        if name[0]=='l':
         
         f.write(str(par.value)+' ')
    f.write("\n") 
f.close()    



    
    
    
with open ("parameter_powlaw2sec.dat",'w') as fp:
    
    for name, par in result1.params.items():
        if name[0]=='p':
         
         fp.write(str(par.value)+' ')
    fp.write("\n")
fp.close()

with open ("parameter_powlaw3sec.dat",'w') as fp:
    
    for name, par in result2.params.items():
        if name[0]=='p':
         
         fp.write(str(par.value)+' ')
    fp.write("\n")
fp.close()

if (cApprox=='F'):
    
    with open("parameter_lmft_fermi_1sec.dat", 'w') as ff:
        for name, par in result.params.items():
            if name[0]=='h':
         
                ff.write(str(par.value)+' ')
        ff.write("\n")
    ff.close()
    
    with open("parameter_lmft_fermi_2sec.dat", 'w') as ff:
        for name, par in result1.params.items():
            if name[0]=='h':
         
                ff.write(str(par.value)+' ')
        ff.write("\n")
    ff.close()
    
    with open("parameter_lmft_fermi_3sec.dat", 'w') as ff:
        for name, par in result2.params.items():
            if name[0]=='h':
         
                ff.write(str(par.value)+' ')
        ff.write("\n")
    ff.close()
    
    
elif (cApprox=='LN'):    
    with open("parameter_lmft_loggaussian_1sec.dat", 'w') as fln:
        for name, par in result.params.items():
            if name[0]=='g':
         
                fln.write(str(par.value)+' ')
         
        fln.write("\n")
    fln.close()
    
    with open("parameter_lmft_loggaussian_2sec.dat", 'w') as fln2:
        for name, par in result1.params.items():
            if name[0]=='g':
         
                fln2.write(str(par.value)+' ')
         
        fln2.write("\n")
    fln2.close()
    
elif (cApprox=='FA'): 
    ##first sector up to cutoff1
    with open("parameter_lmft_fano_1sec.dat", 'w') as ffa1:
        for name, par in result.params.items():
            if name[0]=='f':
         
                ffa1.write(str(par.value)+' ')
        ffa1.write("\n")
    ffa1.close()
    ##secondsector cutoff1< >cutoff2
    
    with open("parameter_lmft_fano_2sec.dat", 'w') as ffa2:
        for name, par in result1.params.items():
            if name[0]=='f':
         
                ffa2.write(str(par.value)+' ')
        ffa2.write("\n")
    ffa2.close()
    
    with open("parameter_lmft_fano_3sec.dat", 'w') as ffa2:
        for name, par in result2.params.items():
            if name[0]=='f':
         
                ffa2.write(str(par.value)+' ')
        ffa2.write("\n")
    ffa2.close()
else:
    ##secondsector cutoff1< >cutoff2
    
    with open("parameter_lmft_lor_2sec.dat", 'w') as ffa2:
        for name, par in result1.params.items():
            if name[0]=='l':
         
                ffa2.write(str(par.value)+' ')
        ffa2.write("\n")
    ffa2.close()
    
    with open("parameter_lmft_lor_3sec.dat", 'w') as ffa2:
        for name, par in result2.params.items():
            if name[0]=='l':
         
                ffa2.write(str(par.value)+' ')
        ffa2.write("\n")
    ffa2.close()
    
fq.close()

   



myutil.sutor()