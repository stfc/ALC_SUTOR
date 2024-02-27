# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:04:47 2023

Main for the Inelastic Mean Free Path calculations.
Command lines parameters passed: ELF file, Fermi energy EF, Band gap Eb, and Max Kinetic energy Tmax

Authors: P. E. Trevisanutto
"""
import numpy as np
import classutility as cu
from tqdm import tqdm
import sys
import class_IMFP as impf
#-----------------------------------------------------
###main ###main ###main ###main ###main ###main ###main
#-----------------------------------------------------
EF=0.5
#Band gap (eV)
Eb = 1.24 
Tmax =2000
ComuProb =[]
list_W_step_ComProb =[]
myutil =cu.utilities()

myutil.letsstart()


#'Nist_ELF_CsPbBr3_0_2_4kev.dat'
try:
    filename = sys.argv[1]
    EF = float (sys.argv[2])
    Eb = float (sys.argv[3])
    Tmax = float (sys.argv[4])
except:
    print("Command line parameter wrongly inserted: ELF file, Fermi Energy EF, Band gap Eb, Max Kinetic Energy Tmax.")
else:
    print ("ELF file: "+ filename, flush=True)
    print ("Fermi Energy [eV]: "+ sys.argv[2], flush=True)
    print ("Band gap Energy [eV]: "+sys.argv[3], flush=True)
    print ("Max Kinetic Energy [eV]: "+sys.argv[4], flush=True)
    
imfp_cl =  impf.class_imfp(EF,Eb,Tmax,myutil)
#   f.write("File starting!")
par_q,par_l,par_f = imfp_cl.reading_parameters(filename)

Integral_tot = 0.0

W = imfp_cl.get_Wmin()

imfp_cl.inel_mean_free_path(par_q,par_l,par_f)


list_of_T_energies= imfp_cl.generate_array(imfp_cl.get_Wmin(), Tmax,step1=15)

for T in list_of_T_energies:
    
      if ((T+Eb)/2 > T-EF):
              
          Wmax=T-EF 
                  
      else:
         
          Wmax=(T+Eb)/2 
    
      print ("Kinetic energy: "+ str(T) +" eV \n")

    
      Integral_tot = 0.0
    
      Wmax= round(Wmax,imfp_cl.get_digit())
    
      print ("W_min ="+str(imfp_cl.get_Wmin()))
    
      print("W_max="+str(Wmax))
    
      
      listW =imfp_cl.generate_array(imfp_cl.get_Wmin(), Wmax,step2=20)
        
      W_in = imfp_cl.get_Wmin()
      
      for W in tqdm( listW):
          
          Integral_step =imfp_cl.diimfp(W_in,W, T, par_q,par_l,par_f)
       
          Integral_tot = Integral_tot + Integral_step
        
          ComuProb.append(Integral_tot)
          
          list_W_step_ComProb.append(W)
          
          W_in = W
    
    
      Cumu_Prob_renorm = np.divide(ComuProb,Integral_tot)
    
    
      myutil.save_all(list_W_step_ComProb, Cumu_Prob_renorm, "P_inel_T="+str(T)+"_eV.dat")
      
      myutil.plotty(list_W_step_ComProb, Cumu_Prob_renorm, "P_inel_T= "+ str(T)+" eV", color="cyan",scale="normal",ylabel="Com. Prob.")
    
      list_W_step_ComProb.clear()
      ComuProb.clear()
      Cumu_Prob_renorm = 0      

myutil.sutor()
