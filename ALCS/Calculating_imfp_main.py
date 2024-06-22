 
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:04:47 2023

Main for the Inelastic Mean Free Path calculations.
Command lines parameters passed: ELF file, Fermi energy EF, Band gap Eb, and Max Kinetic energy Tmax

BSD 3-Clause License

 

Copyright (c) 2024, Ada Lovelace Centre, Science and Technology Facilities Council, part of the UKRI (UK).

Author:             Paolo Emilio Trevisanutto.

 

All rights reserved.

 

Redistribution and use in source and binary forms, with or without

modification, are permitted provided that the following conditions are met:

 

1. Redistributions of source code must retain the above copyright notice, this

   list of conditions and the following disclaimer.

 

2. Redistributions in binary form must reproduce the above copyright notice,

   this list of conditions and the following disclaimer in the documentation

   and/or other materials provided with the distribution.

 

3. Neither the name of the copyright holder nor the names of its

   contributors may be used to endorse or promote products derived from

   this software without specific prior written permission.

 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"

AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE

IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE

DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE

FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL

DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR

SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER

CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,

OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE

OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import colorama
from colorama import Fore, Style
import numpy as np
import classutility as cu
from tqdm import tqdm
import sys
import class_IMFP as impf
#-----------------------------------------------------
###main ###main ###main ###main ###main ###main ###main
#-----------------------------------------------------
EF=0
#Band gap (eV)
Eb=0
Tmin=10
Tmax =2*10**3
cutoff1 = 1.5*10**3
#cutoff1 = 1*10**4
cutoff2 = 2*10**3


Integral_tot = 0.0

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
    sys.exit()
else:
    print ("ELF file: "+ filename+"\n", flush=True)
    print ("Fermi Energy [eV]: "+ sys.argv[2]+"\n", flush=True)
    print ("Band gap Energy [eV]: "+sys.argv[3]+"\n", flush=True)
    print ("Max Kinetic Energy [eV]: "+sys.argv[4]+"\n", flush=True)
    
if (len(sys.argv))> 5:
   cutoff1= float(sys.argv[5])
   cutoff2= float(sys.argv[6])

print ("Cutoff1 [eV]: "+str(cutoff1)+"\n", flush=True)
print ("Cutoff2 [eV]: "+str(cutoff2)+"\n", flush=True)
    
imfp_cl =  impf.class_imfp(EF,Eb,Tmax, myutil,cutoff1,cutoff2)
#   f.write("File starting!")
#self.__params_q,self.__params_p,self.__params_l,self.__params_Nist1,self.__params_pl1 , self.__params_Nist2,self.__params_pl2, final_elf, self.__cApprox 
print(f"{Fore.RED}Approximation----number of parameters \n{Style.RESET_ALL}",flush=True)


par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3,par_f3,tot_elf, Approximation = imfp_cl.reading_parameters(filename)


W = imfp_cl.get_Wmin()
##ora modificare qui.
print(f"{Fore.YELLOW}Inelastic mean free path\n",flush=True)

imfp_cl.inel_mean_free_path(par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3,par_f3, Approximation)

print (imfp_cl.get_Wmin())
list_of_T_energies= imfp_cl.generate_array(Tmin, Tmax)

print (list_of_T_energies)
for T in list_of_T_energies:
    
      if ((T+Eb)/2 > T-EF):
              
          Wmax=T-EF 
                  
      else:
         
          Wmax=(T+Eb)/2 
    
      print (f"{Fore.GREEN}\n Kinetic energy: "+ str(T) +" eV \n")

    
      Integral_tot = 0.0
    
      Wmax= round(Wmax,imfp_cl.get_digit())
    
      print ("W_min ="+str(imfp_cl.get_Wmin()))
    
      print("W_max="+str(Wmax))
    
      
      listW =imfp_cl.generate_array(imfp_cl.get_Wmin(), Wmax,step2=20)
        
      
      
      for W in tqdm( listW):
          
          Integral_step = imfp_cl.integral_step (imfp_cl.get_Wmin(),W, T, par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3, par_f3, Approximation)
       
          Integral_tot = Integral_tot + abs(Integral_step)
        
          ComuProb.append(Integral_tot)
          
          list_W_step_ComProb.append(W)
          
          W_in = W
    
    
      Cumu_Prob_renorm = np.divide(ComuProb,Integral_tot)
    
    
      myutil.save_pinel(list_W_step_ComProb, Cumu_Prob_renorm, "PINEL"+str(T))
      
      myutil.plotty(list_W_step_ComProb, Cumu_Prob_renorm, "PINEL"+ str(T), color="cyan",scale="normal",ylabel="Com. Prob.")
    
      list_W_step_ComProb.clear()
      ComuProb.clear()
      Cumu_Prob_renorm = 0      

myutil.sutor()
