# -*- coding: utf-8 -*-
"""
ALC_SUTOR Project

Script for the Sum rule calculations.

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
Authors: P.E. Trevisanutto.
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
