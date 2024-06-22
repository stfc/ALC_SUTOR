# -*- coding: utf-8 -*-
"""
Created on Tue May 30 11:07:59 2023

Sutor Project
Class for calculations of Inelastic Mean Free Path and Comulative Probabilities for the IMFP 
at different initial kinetic energies.

Authors: P.E. Trevisanutto, 
"""
#from multipledispatch import dispatch
import matplotlib.pyplot as plt
import numpy as np
import math as math
import scipy.integrate as it
from tqdm import tqdm
import classutility as cu
import sys
import os.path
from typing import overload
##eq 4.64 pag46 Transport of Energetic Electrons in Solids.




class class_imfp():
        

    __EF=0.5
#Band gap (eV)
    __Eb = 1.24 
    __xcutoff = 0.0



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
    om_min            = 0.3
    om_max            = __Tmax
    __File_powerlaw   = 'parameter_powlaw.dat'
    __File_powerlaw1  = 'parameter_powlaw2sec.dat'
    __File_powerlaw2  = 'parameter_powlaw3sec.dat'
    __File_quadratic  = 'parameter_lmft_quadratic.dat'
    __File_lorentzian = 'parameter_lmft_lor.dat'
    __File_lorentzian1 = 'parameter_lmft_lor_1sec.dat'
    __File_lorentzian2 = 'parameter_lmft_lor_2sec.dat'
    __File_lorentzian3 = 'parameter_lmft_lor_3sec.dat'
    __File_fermi      = 'parameter_lmft_fermi.dat'
    __File_fermi1     = 'parameter_lmft_fermi_1sec.dat'
    __File_fermi2     = 'parameter_lmft_fermi_2sec.dat'
    __File_fermi3     = 'parameter_lmft_fermi_3sec.dat'
    __File_fano       = 'parameter_lmft_fano.dat'
    __File_fano1      = 'parameter_lmft_fano_1sec.dat'
    __File_fano2      = 'parameter_lmft_fano_2sec.dat'
    __File_fano3      = 'parameter_lmft_fano_3sec.dat'
    __File_gaussian   = 'parameter_lmft_loggaussian.dat'
    __File_gaussian1  = 'parameter_lmft_loggaussian_1sec.dat'
    __File_gaussian2  = 'parameter_lmft_loggaussian_2sec.dat'
    __File_polynomial = 'parameter_lmft_poly.dat'

 
    __cApprox         = "No_approx"
    
    final_elf= []

##other parameters

    __digit           = 0
    
    __DeltaW          = 10**(__digit)
    
    __Wmin            = round(__Eb,__digit)


    
    __Xeas=[]
    __Yeas=[]

    __dW =[]
    __W_cp = []
    __inte = []
    
    
    __params_Nist1 = []
    __params_Nist2 = []
    __params_Nist3 = []
    __params_q     = []
    __params_l     = []
    __params_p     = []
    __params_pl1   = []
    __params_pl2   = []
#    __params_f    = []
#    __params_fa   = []
#    __params_g    = []
    
####--------------------------------------------
    def __init__(self, EF, Eb,Tmax, cu,cutoff1, cutoff2):
        
        self.__EF      = EF
        self.__Eb      = Eb
        self.__Tmax    = Tmax
        self.myut      = cu
        self.cutoff1   = cutoff1
        self.cutoff2   = cutoff2
        self.__Wmin    = self.__Eb
        #self.__Wmin    = round(self.__Eb,self.__digit)
        
####--------------------------------------------
####--------------------------------------------
    def set_approx(self,appr):
        self.__cApprox=appr
###----------

    def get_Wmin(self):
        return self.__Wmin
####--------------------------------------------
    def get_digit(self):
        return self.__digit
####--------------------------------------------
    def generate_array(self,W_min, W_max, step1 =5, step2=10,step3=10**2,step4=10**3,step5=10**4):
        list_w= []
        
        W_st1=100
        W_st2=200
        W_st3=10**3
        W_st4=10**4
        
        
        W=W_min
        print (step1)
        print (step2)
   
    
        while W<W_max:
            #100
            if W<W_st1 : 
    
                W_step = step1
                #100-200
            elif (W>=W_st1) and (W<W_st2):
    
                W_step = step2
                
            #200-1000
            elif (W>=W_st2) and (W<W_st3):
            
                W_step = step3
            #1000-10000
            elif (W >=W_st3) and (W < W_st4):
            
                W_step = step4
             #10000   
            elif (W >=W_st4):
             
                 W_step = step5
                        
            W = W + W_step       
        
            list_w.append(W)
        
        
        return list_w
######-------------------------------------------
    def reading_background(self):
        try:
            ##Reading power law1
            for line in open(self.__File_powerlaw1 , mode='r'):
        
                for s in line.split():
                    self.__params_pl1.append(float(s))
                    
            print ("Power law 1st sector: "+str(len(self.__params_pl1))+"\n")
            ##Reading power law 2
            for line in open(self.__File_powerlaw2 , mode='r'):
        
                for s in line.split():
                    self.__params_pl2.append(float(s))
                    
            print ("Power law 2nd sector: "+str(len(self.__params_pl2))+"\n")

            ##Reading polynomial
            
            for line in open(self.__File_quadratic , mode='r'):
        
                for s in line.split():
                    self.__params_q.append(float(s))
                    
            print ("Quadratic: "+str(len(self.__params_q))+"\n")
            
            ##Reading polynomial
            #
            #for line in open(self.__File_polynomial , mode='r'):
        
            #    for s in line.split():
            #        self.__params_p.append(float(s))
                    
            #print ("Polynomial: "+str(len(self.__params_p))+"\n")
            ##Reading Lorentzian
            for line in open(self.__File_lorentzian, mode='r'):
              ##until omega<co1
                for s in line.split():
                    self.__params_l.append(float(s))
    
            print ("Lorentzians: "+str(len(self.__params_l))+"\n")
            
            counter =len(self.__params_q) + len(self.__params_l)
        except:
            print("No Lorentzian and quadratic files",flush=True)
            self.myut.sutor_issue()
            
        return counter
######-------------------------------------------
    def reading_first_sector(self):
        
        if (os.path.exists(self.__File_fermi1)):
            
            self.__cApprox="Fermi"
            
            print("Sector1 Fermi: "+self.__File_fermi1, flush=True)
            
            for line in open(self.__File_fermi1, mode='r'):
        
                for s in line.split():
                    self.__params_Nist1.append(float(s)) 
                    
            print ("Sector1: SplitLorentzianModel: "+str(len(self.__params_Nist1))+"\n")
            
            
                    
        elif os.path.exists(self.__File_fano1):
            
            self.__cApprox="Fano"
            
            for line in open(self.__File_fano1, mode='r'):
        
                for s in line.split():
                    self.__params_Nist1.append(float(s))
                    
            print ("Sector 1: FanoModel: "+str(len(self.__params_Nist1))+"\n")

            
        elif os.path.exists(self.__File_gaussian1):
           
           
           
           for line in open(self.__File_gaussian1, mode='r'):
       
               for s in line.split():
                   self.__params_Nist1.append(float(s))
                   
           print ("Sector 1: LogNormal: "+str(len(self.__params_Nist1))+"\n")
           
           
        elif os.path.exists(self.__File_lorentzian1):
            
          
           
           for line in open(self.__File_lorentzian1, mode='r'):
       
               for s in line.split():
                   self.__params_Nist1.append(float(s))
                   
           print ("Sector 1: Lorentzian "+str(len(self.__params_Nist1))+"\n")
 
        
        #else:
        #   print ("Not found parameters for Nist ELF first sector\n")            
            
        counter = len(self.__params_Nist1)       
            
        return counter
########------------------------------------------------------------------------------------        
    def reading_second_sector(self):
            
            if (os.path.exists(self.__File_fermi2)):
                
                self.set_approx('Fermi')
                
                print("Sector2 Fermi: "+self.__File_fermi2, flush=True)
                
                for line in open(self.__File_fermi2, mode='r'):
            
                    for s in line.split():
                        self.__params_Nist2.append(float(s)) 
                        
                print ("Sector2: SplitLorentzianModel: "+str(len(self.__params_Nist2))+"\n")
                
                
                        
            elif os.path.exists(self.__File_fano2):
                               
                self.set_approx('Fano')
                
                for line in open(self.__File_fano2, mode='r'):
            
                    for s in line.split():
                        self.__params_Nist2.append(float(s))
                        
                print ("Sector 2: FanoModel: "+str(len(self.__params_Nist2))+"\n")
            
            
            elif os.path.exists(self.__File_gaussian2):
                
                self.set_approx("LN")
                
                for line in open(self.__File_gaussian2, mode='r'):
            
                    for s in line.split():
                        self.__params_Nist2.append(float(s))
                        
                print ("Sector 2: LogNormal: "+str(len(self.__params_Nist2))+"\n")
            
            elif os.path.exists(self.__File_lorentzian2):
                   
                  
                   self.set_approx('Lorentz')
                   
                   for line in open(self.__File_lorentzian2, mode='r'):
               
                       for s in line.split():
                           self.__params_Nist2.append(float(s))
                           
                   print ("Sector 2: Lorentzian "+str(len(self.__params_Nist2))+"\n")    
            else:
                print ("Not found parameters for Nist ELF sector 2 \n")

                
            counter = len(self.__params_Nist2)       
                
            return counter
########------------------------------------------------------------------------------------        
    def reading_third_sector(self):
            
        
            counter = 0
            if (os.path.exists(self.__File_fermi3)):
                
 
                
                print("Sector3 Fermi: "+self.__File_fermi3, flush=True)
                
                for line in open(self.__File_fermi3, mode='r'):
            
                    for s in line.split():
                        self.__params_Nist3.append(float(s)) 
                        
                print ("Sector2: SplitLorentzianModel: "+str(len(self.__params_Nist3))+"\n")
                
                        
            elif os.path.exists(self.__File_fano3):
                                
                for line in open(self.__File_fano3, mode='r'):
            
                    for s in line.split():
                        self.__params_Nist3.append(float(s))
                        
                    print ("Sector 3: FanoModel: "+str(len(self.__params_Nist3))+"\n")
            
            elif os.path.exists(self.__File_lorentzian3):
                   
                   
                   
                   for line in open(self.__File_lorentzian3, mode='r'):
               
                       for s in line.split():
                           self.__params_Nist3.append(float(s))
                           
                   print ("Sector 3: Lorentzian "+str(len(self.__params_Nist3))+"\n")    
                
            else:
                 print ("Not found parameters for Nist ELF sector 3 \n")
                
            counter = len(self.__params_Nist3)       
                
            return counter
#########################---------------------------------
    def reading_parameters(self,filename):
        
        
        
        counter = self.reading_background()
 
# NisT            
        counter1 =self.reading_first_sector()
        
       
        counter2 =self.reading_second_sector()
        
        counter3 =self.reading_third_sector()
    
        counter = counter + counter1 + counter2 + counter3
        
        print ("Parameter total length:" +str(counter) + "\n")
        
        
        final_elf = self.Plotting_analytical_vs_numerical(filename)

        return self.__params_q,self.__params_p,self.__params_l,self.__params_Nist1,self.__params_pl1 , self.__params_Nist2,self.__params_pl2, self.__params_Nist3, final_elf, self.__cApprox

#########################---------------------------------
    def Plotting_analytical_vs_numerical(self,filename):
        final_elf= []
        
        for line in open(filename, mode='r'):
            if (line != "\n"):
                values=line.split()
            
                if (self.myut.is_number(values[0]) == False):
                    pass
                else:
                    self.__Xeas.append(values[0])
                    self.__Yeas.append(values[1]) 
                values.clear()    

        yData = self.myut.From_str2float(self.__Yeas)

        xData = self.myut.From_str2float(self.__Xeas)    

        ind_co1 = self.myut.index_of(xData,self.cutoff1)
        
        ind_co2 = self.myut.index_of(xData,self.cutoff2)
        
        print ("Cutoff first sector:"+str(self.cutoff1), flush=True)
        print ("Cutoff second sector "+str(self.cutoff2), flush=True)
        #Here reproducing the quiiii domani
        
        test1 = [ self.myut.multi_lorentz(x,self.__params_q,self.__params_p,self.__params_l,self.__params_Nist1,self.__cApprox ) for x in xData[:ind_co1 ]]
        
        if (self.__cApprox=="Fermi"):
            
            test2 = [ self.myut.powerlaw_fermi( x,self.__params_pl1 , self.__params_Nist2) for x in xData[ ind_co1 : ind_co2] ]
        
            test3 = [ self.myut.powerlaw_fermi( x, self.__params_pl2, self.__params_Nist3) for x in xData[ind_co2: ] ]
        
        elif (self.__cApprox=="Lorentz"):
            
            test2 = [ self.myut.powerlaw_lorentz( x,self.__params_pl1 , self.__params_Nist2) for x in xData[ ind_co1 : ind_co2] ]
        
            test3 = [ self.myut.powerlaw_lorentz( x, self.__params_pl2, self.__params_Nist3) for x in xData[ind_co2: ] ]
            
        elif (self.__cApprox=="Fano"):
            
            test2 = [ self.myut.powerlaw_fano( x,self.__params_pl1 , self.__params_Nist2) for x in xData[ ind_co1 : ind_co2] ]
           
            test3 = [ self.myut.powerlaw_fano( x, self.__params_pl2, self.__params_Nist3) for x in xData[ind_co2: ] ] 
            
        elif (self.__cApprox=="LN"):
            
            test2 = [ self.myut.powerlaw_ln( x,self.__params_pl1 , self.__params_Nist2) for x in xData[ ind_co1 : ind_co2] ]
        else:
            test2 = [ self.myut.powerlaw_lorentz( x,self.__params_pl1 , self.__params_Nist2) for x in xData[ ind_co1 : ind_co2] ]
           
            test3 = [ self.myut.powerlaw_lorentz( x, self.__params_pl2, self.__params_Nist3) for x in xData[ind_co2: ] ] 
    
        

        fig = plt.figure()
        ax = fig.add_subplot( 1, 1, 1 )
        ax.loglog( xData, yData )
        ax.loglog( xData[: ind_co1], test1,color='red' )
        ax.loglog( xData[ ind_co1:ind_co2], test2, color='red')
        ax.loglog( xData[ind_co2:], test3, color='red' )
        plt.title('Reconstructed spectrum', color='red')
        plt.show()
        
        final_elf.append(test1)
        final_elf.append(test2)
        final_elf.append(test3)
        
         
#####----------------------------------------------------------------------------
    def evolution_Lor(self,q,par_f):
        
        parati_l=[]
        
        for i in range (0, len( par_f), 3):
            omega = par_f[i]
            Ampl  = par_f[i+1]
            gam   = par_f[i+2]
        
        #if (Ampl > 0):
            if (omega > 0):  
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
                gam = (gam**2 + q**2/2 + q**4/4)**(0.5)
            else:
                omega = omega
                gam = gam
            
           #if (gam>0):
            #    gam = (gam**2 + q**2/2 + q**4/4)**(0.5)
            #else:
            #    gam = -(gam**2 + q**2/2 + q**4/4)**(0.5)
        #print (par_l[i])
            parati_l.append(omega)
            parati_l.append(Ampl)
            parati_l.append(gam)
        
        return parati_l
#####----------------------------------------------------------------------------
    def evolution_Maur(self,q,par_f):
        parati_h = []
        
        for i in range (0, len( par_f), 4):  
            omega_h = par_f[i]
            Ampl    = par_f[i+1]
            gam_h   = par_f[i+2]
            gam_r   = par_f[i+3]
    
            if (omega_h > 0):   
                omega_h = (omega_h**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            else:
                omega_h = (omega_h**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
        
            gam_h   = (gam_h**2 + q**2/2 + q**4/4)**(0.5)
            gam_r = (gam_r**2 + q**2/2 + q**4/4)**(0.5)
        
            parati_h.append( omega_h)
            parati_h.append(Ampl)
            parati_h.append( gam_h)
            parati_h.append( gam_r)
        
        return parati_h
#####----------------------------------------------------------------------------
    def evolution_Maur_F(self,q,par_f):
        parati_h = []
        
        for i in range (0, len( par_f), 4):  
            omega = par_f[i]
            Ampl  = par_f[i+1]
            gam   = par_f[i+2]
            B     = par_f[i+3]
    
            if (omega > 0):   
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
                gam   = (gam**2 + q**2/2 + q**4/4)**(0.5)
            else:
                omega = -(omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
        
            
            
        
            parati_h.append( omega)
            parati_h.append(Ampl)
            parati_h.append(gam)
            parati_h.append(B)
      
        return parati_h
#####----------------------------------------------------------------------------
    def evolution_Fano(self,q,par_f):
        parati_h = []
        
        
        for i in range (0, len( par_f), 5):  
            Ampl = par_f[i]
            omega = par_f[i+1]
            gam   = par_f[i+2]
            q_par = par_f[i+3]
            
            
            
            if (omega > 0):   
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
                
            else:
                omega = -(omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
        
            gam   = (gam**2 + q**2/2 + q**4/4)**(0.5)
            
            
            parati_h.append(Ampl)
            parati_h.append(omega)
            parati_h.append( gam)
            parati_h.append(q_par)
            parati_h.append(1)
            
        return parati_h
#####----------------------------------------------------------------------------
    def evolution_LogNorm(self,q,par_f):
        
        parati_h = []
        
        for i in range (0, len( par_f), 3):  
            Ampl = par_f[i]
            omega  = par_f[i+1]
            gam   = par_f[i+2]
    
            if (omega > 0):   
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
            else:
                omega = (omega**2+ 12*self.__EF*q**2/10 + (0.5*q**2)**2)**(0.5)
        
            gam   = (gam**2 + q**2/2 + q**4/4)**(0.5)
            
        
            parati_h.append( Ampl)
            parati_h.append(omega)
            parati_h.append( gam)
            
        
        return parati_h
###-------------------------------------------------------------------------
    
    def integrand2(self,kappa, E,T, par_pl,par_f,Approximation):
        
        
  ### second sector power law +fano      
        parati_Nist = []
        
        q = kappa*self.__hb2m**(0.5)
        
        if (Approximation =="Fermi"):
            
           #parati_Nist = self.evolution_Maur(q,par_f)
           parati_Nist = self.evolution_Maur_F(q, par_f)
           
           ELF = self.myut.powerlaw_fermi( E,par_pl , parati_Nist)
           
        elif (Approximation =="Fano"):
            
            parati_Nist = self.evolution_Fano(q,par_f)
            
            ELF = self.myut.powerlaw_fano( E,par_pl , parati_Nist)
            
        elif (Approximation =="LogNormal"):
            
            parati_Nist = self.evolution_LogNorm(q,par_f)
            ELF = self.myut.powerlaw_ln( E,par_pl , parati_Nist)
        elif (Approximation =="Lorentz"):
            print("ttr")
            parati_Nist = self.evolution_Lor(q,par_f)
            ELF = self.myut.powerlaw_lorentz( E,par_pl , parati_Nist)
        else:
            
            #self.myut.sutor_issue("Integrand did not find any approximation")
            ELF=self.myut.Power_law(E,par_pl)
            
        
        
        #
        
        parati_Nist.clear()
        return ELF/(math.pi*self.__ao*T*kappa)
###-------------------------------------------------------------------------
    
    def integrand3 (self,kappa, E,T, par_pl,par_f, Approximation):
        ### third sector power law   
        parati_Nist = []
        
        q = kappa*self.__hb2m**(0.5)
        
        if (Approximation =="Fermi"):
            
           #parati_Nist = self.evolution_Maur(q,par_f)
           parati_Nist = self.evolution_Maur_F(q, par_f)
           
           ELF = self.myut.powerlaw_fermi( E,par_pl , parati_Nist)
           
           
        elif (Approximation =="Fano"):
            
            parati_Nist = self.evolution_Fano(q,par_f)
            
            ELF = self.myut.powerlaw_fano( E,par_pl , parati_Nist)
            
        elif (Approximation =="LogNormal"):
            
            parati_Nist = self.evolution_LogNorm(q,par_f)
        
        else:
            ELF=self.myut.power_law(par_pl)
            #self.myut.sutor_issue("Integrand did not find any approximation")
        
        
        parati_Nist.clear()
        return ELF/(math.pi*self.__ao*T*kappa)
            
###-------------------------------------------------------------------------
    
    def integrand1 (self,kappa, E,T, par_q,par_p,par_l,par_f,par_pl1,par_f1,par_pl2, par_f2,Approximation):
    
        parati_l=[]
    
        parati_Nist = []
        parati_Nist1 = []
        parati_Nist2 = []
        
        co1=self.cutoff1
        co2=self.cutoff2
        
        q = kappa*self.__hb2m**(0.5)
        
        parati_l = self.evolution_Lor(q, par_l)

        #print(par_l[i])
        
        if (Approximation =="Fermi"):
            
           #parati_Nist = self.evolution_Maur(q,par_f)
           parati_Nist = self.evolution_Maur_F(q, par_f)
           parati_Nist1 = self.evolution_Maur_F(q, par_f1)
           parati_Nist2 = self.evolution_Maur_F(q, par_f2)
           
           ELF1 = self.myut.powerlaw_fermi( E,par_pl1 , parati_Nist1)
           
           ELF2 = self.myut.powerlaw_fermi( E,par_pl2 , parati_Nist2)
           
        elif (Approximation =="Fano"):
            
            parati_Nist = self.evolution_Fano(q,par_f)
            parati_Nist1 = self.evolution_Fano(q,par_f1)
            parati_Nist2 = self.evolution_Fano(q,par_f2)
            
            ELF1 = self.myut.powerlaw_fano( E,par_pl1 , parati_Nist1)
            
            ELF2 = self.myut.powerlaw_fano( E,par_pl2 , parati_Nist2)
            
        elif (Approximation =="LogNormal"):
            
            parati_Nist = self.evolution_LogNorm(q,par_f)
        
        elif (Approximation =="Lorentz"):
            
            parati_Nist = self.evolution_Lor(q,par_f)
            parati_Nist1 = self.evolution_Lor(q,par_f1)
            parati_Nist2 = self.evolution_Lor(q,par_f2)
            
            ELF1 = self.myut.powerlaw_lorentz( E,par_pl1 , parati_Nist1)
            ELF2 = self.myut.powerlaw_lorentz( E,par_pl2 , parati_Nist2)
        
       # else:
            
        #    self.myut.sutor_issue("Integrand did not find any approximation")
            
    ##Here E comes in.
        ELF = self.myut.multi_lorentz(E, par_q,par_p,parati_l,parati_Nist,Approximation)
       
        
        ELF = ELF*np.heaviside(co1-E,1/2)+ ELF1*np.heaviside(E-co1,1/2)*np.heaviside(co2-E,1/2) + ELF2*np.heaviside(E-co2,1/2)
        
        
        
        y = ELF/(math.pi*self.__ao*T*kappa)
        
        parati_l.clear()
        parati_Nist.clear()
        parati_Nist1.clear()
        parati_Nist2.clear()
        #print (y)
        return y
#-----------------------------------------------------
    #parameter pl: power law, f: fano.
 #-----------------------------------------------------   
    def diimfp3(self,Emin,Emax,T, par_pl, par_f,Approximation):
    ##2nd sector    
#        kmin=(2./hb2m)**(0.5)*((T)**(0.5) - (abs(T-Emin))**(0.5))
        #print(kmin)
#        kmax=(2./hb2m)**(0.5)*((T)**(0.5) + (abs(T-Emax))**(0.5)) 
        #print(kmax)
            
            if Emax < self.__Wmin:
                result = 0
            else:
                da_integrare = lambda kappa, E ,kin,pl,Approximation: self.integrand3(kappa, E,kin,pl, par_f, Approximation)
                imfp_integral = it.dblquad(da_integrare,Emin,Emax,lambda E:(2./self.__hb2m)**(0.5)*(T**(0.5) - (abs(T-E))**(0.5)),\
                                   lambda E:(2./self.__hb2m)**(0.5)*((T)**(0.5) + (abs(T-E))**(0.5)),args=(T,par_pl))
                
            result = imfp_integral[0]

            return  result
#-----------------------------------------------------
    #parameter pl: power law, f: fano.
    
    def diimfp2(self,Emin,Emax,T, par_pl,par_f,Approximation):
    ##2nd sector    
#        kmin=(2./hb2m)**(0.5)*((T)**(0.5) - (abs(T-Emin))**(0.5))
        #print(kmin)
#        kmax=(2./hb2m)**(0.5)*((T)**(0.5) + (abs(T-Emax))**(0.5)) 
        #print(kmax)
            
            if Emax < self.__Wmin:
                result = 0
            else:
                da_integrare = lambda kappa, E ,kin,pl,f,Approximation: self.integrand2(kappa, E,kin,pl,f,Approximation )
                imfp_integral = it.dblquad(da_integrare,Emin,Emax,lambda E:(2./self.__hb2m)**(0.5)*(T**(0.5) - (abs(T-E))**(0.5)),\
                                   lambda E:(2./self.__hb2m)**(0.5)*((T)**(0.5) + (abs(T-E))**(0.5)),args=(T,par_pl,par_f,Approximation))
                
            result = imfp_integral[0]

            return  result
#-----------------------------------------------------
    
    def diimfp1(self,Emin,Emax,T, par_q,par_p,par_l,par_f,par_pl1,par_f1,par_pl2,par_f2 ,Approximation): 
    
            
        if Emax < self.get_Wmin():
            result = 0
        else:
                
            da_integrare = lambda kappa, E ,kin,q,p,l,f,pl1,f1,pl2,f2,Approximation: self.integrand1(kappa, E,kin,q,p,l,f,pl1,f1,pl2,f2,Approximation )
            imfp_integral = it.dblquad(da_integrare,Emin,Emax,lambda E:(2./self.__hb2m)**(0.5)*(T**(0.5) - (abs(T-E))**(0.5)),\
                                   lambda E:(2./self.__hb2m)**(0.5)*((T)**(0.5) + (abs(T-E))**(0.5)),args=(T,par_q,par_p,par_l,par_f,par_pl1, par_f1,par_pl2,par_f2,Approximation))
                
        result = imfp_integral[0]
            
            #print(result)
        return  result
##====================================================

    def inel_mean_free_path(self,par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3,par_f3,cApprox):
    
    
    
        imfp_integral=[]
        T_variation_list = self.generate_array(self.__Wmin, self.__Tmax,step1 =3, step2=5)
        
    
        #for T in tqdm(T_variation_list):
        for T in (T_variation_list):   
            if ((T+self.__Eb)/2 > T-self.__EF):
            
                    Wmax=T-self.__EF 
                
            else:
            
                    Wmax=(T+self.__Eb)/2
                
            Wmax= round(Wmax,self.__digit)
            print (self.get_Wmin())
            print ("Wmax")
            print (Wmax)
           
            Integral_tot =self.diimfp1( self.get_Wmin(), Wmax, T, par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3,par_f3,cApprox)
           
            print (Integral_tot)
            imfp_integral.append((1/(self.__N_mix*Integral_tot)))
        
        imfp_np =np.asarray(imfp_integral)
        
        print("\n")
        print ('min IMFP= '+ str(min(imfp_np[imfp_np > 0]))+"$\AA$ \n")
        
        i = np.argmin(imfp_np)
        print ('Minimum Kinetic energy: ' + str(T_variation_list[i])+" eV \n")
    
        self.myut.save_all(T_variation_list, imfp_integral, "IMFP")
        self.myut.plotty(T_variation_list,imfp_integral, "Inelastic Mean Free Path",scale='loglog')
##====================================================
    def integral_step (self,W_in,Wmax, T,par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3, par_f3,cApprox):


        Integral_tot =self.diimfp1(self.__Wmin,Wmax, T, par_q,par_p,par_l,par_f1, par_pl2, par_f2,par_pl3,par_f3,cApprox)
                

        return Integral_tot    
##====================================================