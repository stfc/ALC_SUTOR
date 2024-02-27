# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 12:10:50 2023

@author: PE Trevisanutto
"""
import sys

import numpy as np

import matplotlib.pyplot as plt

from scipy import interpolate

from scipy.interpolate import UnivariateSpline

#_________________________________________________________________________
#-----------------------------------------
def From_str2float(str):
    rvect_f =[]
    for item in str:
        rvect_f.append(float(item))
        
    rv = np.asarray(rvect_f)
    return rv
#_________________________________________________________________________

def save_all(x,y,name):

    
    output = '\n'.join('\t'.join(map(str,row)) for row in zip(x,y))
    
    with open(name+".dat", 'w') as f:
        f.write(output)
        f.write("\n")
    
    f.close()
#_________________________________________________________________________
#_________________________________________________________________________

def save_all2(x,y1,y2,name):

    
    output = '\n'.join('\t'.join(map(str,row)) for row in zip(x,y1,y2))
    
    with open(name+".dat", 'w') as f:
        f.write(output)
        f.write("\n")
    
    f.close()
 
#-------------------------------------------------------------------

om=[]
f1=[] 
f2=[]
print("Festina lente!")
#name= 'Br_f1_f2_Nist.txt'
#name= 'Cs_f1_f2_Nist.txt'
#name= 'Pb_f1_f2_Nist.txt'
name =sys.argv[1]
with open(name, 'r') as f:
  next(f)
  for line in f:
      values = [float(s) for s in line.split()]
      om.append(values[0]*10**3)
      f1.append(values[1])
      f2.append(values[2])
      #pip.append(values[2])
  
om_ar=From_str2float(om)
f1_ar=From_str2float(f1)
f2_ar=From_str2float(f2)

xnew = np.arange(20, 430000, 1)


f1_int=interpolate.interp1d(om_ar,f1_ar)



ynew_f1 = f1_int(xnew)   # use interpolation function returned by `interp1d`





f2_int=interpolate.interp1d(om_ar,f2_ar)


ynew_f2 = f2_int(xnew)   # use interpolation function returned by `interp1d`


save_all2(xnew,ynew_f1,ynew_f2,name+'_inter')

plt.loglog(om,f1, 'o')
plt.loglog( xnew, ynew_f1, '-',label='f1')
plt.loglog(om,f2, 'o')
plt.loglog(xnew, ynew_f2, '-',label='f2')
plt.xlabel("$\omega (eV)$")
plt.title(name)
plt.legend(loc="lower left")
plt.grid(True)
plt.show()

#save_all2(xnew,ynew_f1,ynew_f2,'Br_f1_f2_Nist_inter.txt')
#save_all2(xnew,ynew_f1,ynew_f2,'Cs_f1_f2_Nist_inter.txt')

#timestamp = (0,5,10,15,30,35,40,50,55,60)

#distance = (100,90,65,85,70,30,40,45,20,0)


#plt.plot(Xeas, Yeas, 'o')

#plt.show()
print("Sutor ne ultra crepidam.")
