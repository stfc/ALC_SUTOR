# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:07:59 2023

ALC_SUTOR Project

Class for calculations of Inelastic Mean Free Path and Cumulative Probabilities for the IMFP 
at different initial kinetic energies.

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

import matplotlib.pyplot as plt
import numpy as np
import math as math
import scipy.integrate as it
from tqdm import tqdm
import classutility as cu
import sys
##eq 4.64 pag46 Transport of Energetic Electrons in Solids.




class class_imfp():
        

    __EF=0.5
#Band gap (eV)
    __Eb = 1.24 



    __hb2m=7.615 ##(hb*c)^2/0.5 MeV=[eV*Angs]^2
    __ao=0.529 #Bohr radius
##Here to generate not starting from beginning
    __T  = 10
    __dT = 5
    __Tmax=2000.
    __off = 0
    __dens_mix =5.01
    __N_mix= 1


###Atomic Number
#    A_Br = 79.90400
#    A_Cs = 132.9054
#    A_Pb =  207.2000
#    N_mix= dens_mix*N_avo/(((12*A_Br+ 4*A_Pb+ 4*A_Cs))*10)
#    N_mix=1

##Input that I have to insert in an input file.
#Fermi Energy

##these are just for plots
    om_min = 0.3
    om_max = __Tmax

    __File_quadratic  = 'parameter_lmft_quadratic.dat'
    __File_lorentzian = 'parameter_lmft_lor.dat'
    __File_fermi      = 'parameter_lmft_fermi.dat'

##other parameters

    __digit=0
    
    __DeltaW = 10**(__digit)
    
    __Wmin=round(__Eb,__digit)



    __Xeas=[]
    __Yeas=[]

    __dW =[]
    __W_cp = []
    __inte = []
    
    __params_q = []
    __params_l = []
    __params_f = []
####--------------------------------------------
    def __init__(self, EF, Eb,Tmax,cu):
        
        self.__EF   = EF
        self.__Eb   = Eb
        self.__Tmax = Tmax
        self.myut = cu

####--------------------------------------------
    def get_Wmin(self):
        return self.__Wmin
####--------------------------------------------
    def get_digit(self):
        return self.__digit
####--------------------------------------------
    def generate_array(self,W_min, W_max, step1 =1, step2=100,step3=1000):
        list_w= []
    
        W=W_min+1
    
   
    
        while W<W_max:
    
            if W<=100.0 : 
    
                W_step = step1
    
            elif (W>100) and (W<1000):
    
            
                W_step = step2
            
            elif (W>100) and (W<1000):
            
                W_step = step3
            
            
            W = W + W_step       
        
            list_w.append(W)
        
        
        return list_w
###-------------------------------------------------------------------------
    def integrand (self,kappa, E,T, par_q,par_l,par_f):
    
        parati_l=[]
    
        parati_h = []
    
    
        q = kappa*self.__hb2m**(0.5)
        
        for i in range (0, len( par_l), 3):
            omega= par_l[i]
            Ampl = par_l[i+1]
            gam = par_l[i+2]
        
        #if (Ampl > 0):
            if (omega > 0):  
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            else:
                omega = - (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            
            gam = (gam**2 + q**2/2 + q**4/4)**(0.5)

        #print (par_l[i])
            parati_l.append(omega)
            parati_l.append(Ampl)
            parati_l.append(gam)
        #print(par_l[i])
        
        
        for i in range (0, len( par_f), 4):  
            omega_h = par_f[i]
            Ampl  = par_f[i+1]
            gam_h = par_f[i+2]
            gam_r = par_f[i+3]
        
            if (omega_h > 0):   
                omega_h = (omega_h**2+ 12*self.EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            else:
                omega_h = (omega_h**2+ 12*self.EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            
            gam_h   = (gam_h**2 + q**2/2 + q**4/4)**(0.5)
            gam_r = (gam_r**2 + q**2/2 + q**4/4)**(0.5)
            
            parati_h.append( omega)
            parati_h.append(Ampl)
            parati_h.append( gam_h)
            parati_h.append( gam_r)
        
    ##Here E comes in.
    
        ELF = self.myut.multi_lorentz(E, par_q,parati_l,parati_h)    
        parati_l.clear()
        parati_h.clear()
        return ELF/(math.pi*self.__ao*T*kappa)
######-------------------------------------------
    def reading_parameters(self,filename):
        try:
            for line in open(self.__File_quadratic , mode='r'):
        
                for s in line.split():
                    self.__params_q.append(float(s))
    
            for line in open(self.__File_lorentzian, mode='r'):
        
                for s in line.split():
                    self.__params_l.append(float(s))

            for line in open(self.__File_fermi, mode='r'):
        
                for s in line.split():
                    self.__params_f.append(float(s))  
        except:
            print("Wrong files")
        counter =len(self.__params_q)+len(self.__params_f)+len(self.__params_l)
    
        print ("Parameter total length:" +str(counter)+"\n")
        print ("Quadratic: "+str(len(self.__params_q))+"\n")
        print ("Lorentzians: "+str(len(self.__params_l))+"\n")
        print ("SplitLorentzianModel: "+str(len(self.__params_f))+"\n")
    
    
        for line in open(filename, mode='r'):
            values=line.split()
            if (self.myut.is_number(values[0]) == False):
                pass
            else:
                self.__Xeas.append(values[0])
                self.__Yeas.append(values[1]) 
            values.clear()    

        yData = self.myut.From_str2float(self.__Yeas)

        xData = self.myut.From_str2float(self.__Xeas)    

    


        test = [ self.myut.multi_lorentz(x, self.__params_q,self.__params_l,self.__params_f ) for x in xData ]
        fig = plt.figure()
        ax = fig.add_subplot( 1, 1, 1 )
        ax.loglog( xData, yData )
        ax.loglog( xData, test )

        plt.title('Reconstructed spectrum')
        plt.show()
    
        return self.__params_q,self.__params_l,self.__params_f 
#####

#-----------------------------------------------------
    def diimfp(self,Emin,Emax,T, par_q,par_l,par_f):
        
#        kmin=(2./hb2m)**(0.5)*((T)**(0.5) - (abs(T-Emin))**(0.5))
        #print(kmin)
#        kmax=(2./hb2m)**(0.5)*((T)**(0.5) + (abs(T-Emax))**(0.5)) 
        #print(kmax)
            if Emax < self.__Wmin:
                result = 0
            else:
                da_integrare = lambda kappa, E ,kin,q,l,h: self.integrand(kappa, E,kin,q,l,h )
                imfp_integral = it.dblquad(da_integrare,Emin,Emax,lambda E:(2./self.__hb2m)**(0.5)*(T**(0.5) - (abs(T-E))**(0.5)),\
                                   lambda E:(2./self.__hb2m)**(0.5)*((T)**(0.5) + (abs(T-E))**(0.5)),args=(T,par_q,par_l,par_f))
                
            result = imfp_integral[0]

            return  result
##====================================================
    def inel_mean_free_path(self,par_q,par_l,par_f):
    
    
    
        imfp_integral=[]
        T_variation_list = self.generate_array(self.__Wmin, self.__Tmax,step1=5)
        print(T_variation_list)
    
        for T in tqdm(T_variation_list):
        #print (T)
            if ((T+self.__Eb)/2 > T-self.__EF):
            
                    Wmax=T-self.__EF 
                
            else:
            
                    Wmax=(T+self.__Eb)/2
                
            Wmax= round(Wmax,self.__digit)
        
            Integral_tot =self.diimfp(self.__Wmin,Wmax, T, par_q,par_l,par_f)
            imfp_integral.append(1/(self.__N_mix*Integral_tot))
        
        print("\n")
        print ('min= '+ str(min(imfp_integral))+"\n")
        i = np.argmin(imfp_integral)
        print ('Minimum energy: ' + str(T_variation_list[i]))
    
        self.myut.save_all(T_variation_list, imfp_integral, "IMFP.dat")
        self.myut.plotty(T_variation_list,imfp_integral, "Inelastic_mean_free_path",scale='loglog')
    
    
#-----------------------------------------------------
###main ###main ###main ###main ###main ###main ###main
#-----------------------------------------------------
"""
EF=0.5
#Band gap (eV)
Eb = 1.24 
Tmax =2000
ComuProb =[]
list_W_step_ComProb =[]
myutil =cu.utilities()

myutil.letsstart()

filename = sys.argv[1]

#'Nist_ELF_CsPbBr3_0_2_4kev.dat'
imfp_cl =  class_imfp(EF,Eb,Tmax,myutil)
#   f.write("File starting!")
par_q,par_l,par_f = imfp_cl.reading_parameters(filename, imfp_cl.File_quadratic, imfp_cl.File_lorentzian, imfp_cl.File_fermi)

Integral_tot = 0.0

W = imfp_cl.Wmin

imfp_cl.inel_mean_free_path(par_q,par_l,par_f)


list_of_T_energies= imfp_cl.generate_array(imfp_cl.Wmin, imfp_cl.Tmax,step1=15)

for T in list_of_T_energies:
    
      if ((T+imfp_cl.Eb)/2 > T-imfp_cl.EF):
              
          Wmax=T-imfp_cl.EF 
                  
      else:
         
          Wmax=(T+imfp_cl.Eb)/2 
    
      print ("Kinetic energy: "+ str(T) +" eV \n")

    
      Integral_tot = 0.0
    
      Wmax= round(Wmax,imfp_cl.digit)
    
      print ("W_min ="+str(imfp_cl.Wmin))
    
      print("W_max="+str(Wmax))
    
      
      listW =imfp_cl.generate_array(imfp_cl.Wmin, Wmax,step2=20)
        
      W_in = imfp_cl.Wmin
      
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
"""
