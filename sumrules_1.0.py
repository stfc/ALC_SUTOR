# -*- coding: utf-8 -*-
"""
Sutor Project
Script for the Sum rule calculations.
Authors: P.E. Trevisanutto, G. Teobaldi.
"""

import sys
from tqdm import tqdm
import classutility as cu
import matplotlib as plt        

##############################MAIN####################################################
#MAIN_MAIN_ MAIN_ MAIN_ MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_MAIN_
##===================================================================================
myutil = cu.utilities()
##Definizioni
SR_P=[]
SR_Z=[]
om_jump=[]
omega =[]
Imeps = []

myutil.letsstart()

name_of_script = sys.argv[1]
#name_of_script ='Nist_ELF_CsPbBr3_0_2_4kev.dat'
#om,imeps = myutil.taking_data(name_of_script)
      
    
   
for line in open(name_of_script, mode='r'):
           
           values=line.split()
           print (values)
           if (myutil.is_number(values[0]) == False):
               pass
           else:
               omega.append(float(values[0])/27.2116)
               Imeps.append(values[1])
           values.clear()   
               
           
       
#myutil.plotty(omega,Imeps,'blue') 

omega=myutil.From_str2float(omega)
imels=myutil.From_str2float(Imeps)

SR_Z.append(0)
om_jump.append(omega[0])



for t in tqdm(range(1,len(omega),50)):


        sumrulef = myutil.Integration(omega[:t],imels[:t])
        
        om_jump.append(omega[t])
        
        SR_Z.append(sumrulef)

#plotty(omega, SR_Z, 'sumrule', 'orange','lin')

print(sumrulef,flush=True)

myutil.plotty(om_jump,SR_Z, "sum_rule",ylabel="Sum rule")

myutil.save_all(om_jump, SR_Z,"sumrulenew")

myutil.sutor()

###END END END END END
#####