# -*- coding: utf-8 -*-
"""

Class for interpreting the Nist f1 and f2 files to calculate the ELF after 100 eV.

P. E. Trevisanutto, G. Teobaldi STFC-UKRI
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi

import sys
from scipy import interpolate

#import csv


Xeas, Yeas = [], []
omeg_final =[]

ELF_semiempirical =[]
###https://www.rsc.org/periodic-table
###https://www.periodic-table.org/what-is-atomic-number-density-definition/
###https://www.periodic-table.org/what-is-atomic-number-density-definition/

class Nist_interpolation():
    N_mix = 0
    File_TDDFT=""
####--------------------------------------------
    def __init__(self,utilities):

        self.ut = utilities
####--------------------------------------------    
    def reading_input(self):

        atom =[]
        N_atom = []
        Number = []
        files = []
        dens = []
        
        denominator = 0 
        iteration = 0
        filename =sys.argv[1]
        
        print ("Reading the file: " + filename)
        
        with open(filename, 'r') as f:
            
            for line in f:
                
                iteration = iteration + 1
                    
                values = [s for s in line.split()]
                
                try: 
                    
                    if iteration == 1 :
                        N_species = int(values[0])
                        print ("Species: " + str(N_species),flush = True)
                        dens_mix = float(values[1])
                        if (values[2] != 0):
                            self.File_TDDFT = values[2]
                            print ("TD-DFT file: "+ self.File_TDDFT)
                        
                    else:
                        
                        atom.append(values[0])
                        Number.append(float(values[1]))
                        N_atom.append(float(values[2]))
                        dens.append(float(values[3]))
                        files.append(values[4])
                except:
                    print ("Wrong input file!", flush=True)
                    print("input file must have the following structure:", flush=True)
                    print ("Number of species, density of system [g/cm^3], TD-DFT file ", flush=True)
                    print ("Atom, Atomic_number, Atomic density[g/cm^3], file Nist with f1 and f2", flush=True)
                    exit()
                    
                else:
                    if (iteration > 1):
                    
                        print ("Number of atoms in the supercell: " +str(Number),flush = True)
                        print ("Atomic number: "+str(N_atom),flush = True)
                        print (" Densities: " +str(dens),flush = True)
                        print ("Files: " + str(files),flush = True)
        
        
        for i in range(N_species):
            denominator = N_atom[i]*Number[i] +denominator
            
        self.N_mix= dens_mix*self.ut.get_Avogadro()/denominator

        return N_species, atom, Number, N_atom, dens, files   
            
####--------------------------------------------    
    def Interpolating_Nist(self,filename):
        om=[]
        f1=[] 
        f2=[]
#name= 'Br_f1_f2_Nist.txt'
#name= 'Cs_f1_f2_Nist.txt'
#name= 'Pb_f1_f2_Nist.txt'
       
        with open(filename, 'r') as f:
#            next(f)
            for line in f:
                values = [float(s) for s in line.split()]
                om.append(values[0]*10**3)
                f1.append(values[1])
                f2.append(values[2]) 
      #pip.append(values[2])
  
        om_ar=self.ut.From_str2float(om)
        f1_ar=self.ut.From_str2float(f1)
        f2_ar=self.ut.From_str2float(f2)

        xnew = np.arange(20, 430000, 1)


        f1_int=interpolate.interp1d(om_ar,f1_ar)

        ynew_f1 = f1_int(xnew)   # use interpolation function returned by `interp1d`


        f2_int=interpolate.interp1d(om_ar,f2_ar)


        ynew_f2 = f2_int(xnew)   # use interpolation function returned by `interp1d`


        self.ut.save_all2(xnew,ynew_f1,ynew_f2,filename+'_inter')

        plt.loglog(om,f1, 'o')
        plt.loglog( xnew, ynew_f1, '-',label='f1')
        plt.loglog(om,f2, 'o')
        plt.loglog(xnew, ynew_f2, '-',label='f2')
        plt.xlabel("$\omega (eV)$")
        plt.title(filename)
        plt.legend(loc="lower left")
        plt.grid(True)
        plt.show()
        f.close()
#=====================================================
    def Get_Core_elfi_Nist(self,file):
        omega =[]
        f1 = []
#        f2 =[]
#        rho =[]
    
        for line in open(file, mode='r'):
        #print (line)
            values = [float(s) for s in line.split()]
            if (values[0] >=0):
                omega.append(values[0]*10**(3))
                f1.append(values[1])
            #f2.append(values[2])
            #rho.append(values[7]*10**(-7))
        
        self.ut.plotty(omega, f1,"f1",'green','loglog')
    
        return omega, f1
####-------------------------------------------- 
    def Get_Core_f1_f2_Nist(self,file):
        omega =[]
        f1 = []
        f2 =[]
        rho =[]
    
        for line in open(file, mode='r'):
        #print (line)
            values = [float(s) for s in line.split()]
            if (values[0] >=0):
                omega.append(values[0]*10**(3))
                f1.append(values[1])
                f2.append(values[2])
                rho.append(values[7]*10**(-7))
        
#        self.ut.plotty(omega, f1,"f1",'green','loglog')
#        self.ut.plotty(omega, f2,"f2",'green','loglog')
    
        print ((len(omega)))
        return omega, f1, f2, rho
####-------------------------------------------- 
    def Get_Core_f1_f2_Nist_interpolated(self,file):
        omega =[]
        f1 = []
        f2 =[]

    
        for line in open(file, mode='r'):
        #print (line)
            values = [float(s) for s in line.split()]
            if (values[0] >=0):
                omega.append(values[0])
                f1.append(values[1])
                f2.append(values[2])

        #self.ut.plotty(omega, f1,"f1",'green','loglog')
        #self.ut.plotty(omega, f2,"f2",'green','loglog')
   
        #print ((len(omega)))
        return omega, f1, f2
####-------------------------------------------- 
    def Get_Core_f1_f2_Henke(self,file):
        omega =[]
        f1 = []
        f2 =[]
    
        for line in open(file, mode='r'):
        #print (line)
            values = [float(s) for s in line.split()]
            if (values[0] >=0):
                omega.append(values[0])
                f1.append(values[1])
                f2.append(values[2])
            
            
            self.ut.plotty(omega, f1,"f1",'green','loglog')
    
        self.ut.plotty(omega, f2,"f2",'green','loglog')
  
        print ((len(omega)))
        return omega, f1, f2
####--------------------------------------------
    def transforming_f1__n(self,omega,rho,f1,N_i,atom):
        n_refr=[]
    
    
    #inum = len(omega)
        f1_np  = self.ut.From_str2float(f1)
        rho_np = self.ut.From_str2float(rho)
        omg_np = self.ut.From_str2float(omega)
    
    #rho =eV2Ang(omg_np)
    
        for i in range(len(omega)):
    
            n_refr.append(1.0-((rho_np[i])*(N_i/(2*pi))*(1/(self.ut.C_henke*pi*omg_np[i]))*f1_np[i]))
    
        self.ut.save_all(omg_np,n_refr, "refractive_"+atom)
    
        self.ut.plotty(omg_np,n_refr,"n refractive index "+ atom,'blue','semi')
    
        n_np = self.ut.From_str2float(n_refr)
 
        return  n_np


####--------------------------------------------
    def transforming_f2__k(self,omega,rho,f2,N_i,atom):
        k_ext=[]
    #k_ext_henke=[]
    
        f2_np = self.ut.From_str2float(f2)
        rho_np =self.ut.From_str2float(rho)
        omg_np = self.ut.From_str2float(omega)
    
    #print ("rho-----------------")
    #print (rho)
        rho_calc = self.ut.eV2cm(omg_np)
    #print ("rho_calc---------------")
    #print (rho_calc)
    #print ("omega-----------------")
    #print (omg_np)
    
        for i in range(len(rho_np)):
        
            k_ext.append((rho_np[i])*(N_i/(2*pi))*(1/(self.ut.C_henke*pi*omg_np[i]))*f2_np[i])
        
        #k_ext_henke.append((rho_np[i]/(2*pi))*(1/(C_henke*pi*omg_np[i]))*f2_np[i])
        #print(-1.0*N_i*r_e*f2_np[i]*rho[i]**2/(2.0*math.pi))
              
        k_np = self.ut.From_str2float(k_ext)
    #k_np_h = From_str2float(k_ext_henke)
        self.ut.save_all(omg_np, k_np, "extinction_" + atom)
  
        self.ut.plotty(omega, k_ext,"extinction "+atom, 'magenta','semi')
   
    
        return k_np
####--------------------------------------------
    def  n_k2elfi(self,n,k, omg, name):
        elf=[]
    
        for i in range (len(omg)):
            elf.append(abs(2.0*n[i]*k[i])/((n[i]**2-k[i]**2)**2+(2.0*n[i]*k[i])**2))
    
        elf_np=self.ut.From_str2float(elf)
        self.ut.save_all(omg, elf_np, "elf_"+ name)
    
        self.ut.plotty(omg, elf_np, "Elf_"+ name,'red','semi')
   
        return elf_np 

####--------------------------------------------
    def transforming_f12__nk_henke_formula(self,omega,rho,f1,f2):
        n_refr=[]
        k_ext =[]
        elf=[]
    #inum = len(omega)
        f1_tot = self.ut.From_str2float(f1)
        f2_tot = self.ut.From_str2float(f2)
    
    
        rho_np = self.ut.From_str2float(rho)
        omg_np = self.ut.From_str2float(omega)
    
        for i in range(len(omega)):
        
        
            n_refr.append(1.0-((rho_np[i]*self.N_mix/(2*pi))*(1/(self.ut.get_henke()*pi*omg_np[i]))*f1_tot[i]))
            k_ext.append((rho_np[i]*self.N_mix/(2*pi))*(1/(self.ut.get_henke()*pi*omg_np[i]))*f2_tot[i])
        
        n = self.ut.From_str2float(n_refr)
        k =self.ut.From_str2float(k_ext)
    
    
        for i in range (len(omega)):
            elf.append(abs(2.0*n[i]*k[i])/((n[i]**2-k[i]**2)**2+(2.0*n[i]*k[i])**2))
    
        elf_final = self.ut.From_str2float(elf)
    
        return elf_final
####--------------------------------------------
    def calc_ELF_compound_Nist(self):
           
            N_species, atom, Number_in_cell, N_atom, dens, files=self.reading_input()
            
            
            
            for i in range(N_species):
                
                print (str(files[i])+ " to be interpolated",flush= True)
                self.Interpolating_Nist(files[i])
                file_inter = files[i]+"_inter.dat"
                print ("file interpolated: " +file_inter ,flush= True)
                om, f1_atom, f2_atom = self.Get_Core_f1_f2_Nist_interpolated(file_inter)
                
                
                if (i==0):
                    
                    f1_henke = np.zeros(len(f1_atom))
                    f2_henke = np.zeros(len(f2_atom))
                    
                    rho_np    = self.ut.eV2cm(om)
                    
                    om_np = self.ut.From_str2float(om)

                f2_np = self.ut.From_str2float(f2_atom)
                
                f2_np = Number_in_cell[i]*f2_np
                
                f1_np = self.ut.From_str2float(f1_atom)
                
                f1_np = Number_in_cell[i]*f1_np
                
                ##according to formula 4.47,8 "Electron Beam interanction with solids Dapor" (2003)
                f1_henke = np.add(f1_np,f1_henke)
                f2_henke = np.add(f2_np,f2_henke)
                
            Elf_henke = self.transforming_f12__nk_henke_formula(om_np,rho_np,f1_henke,f2_henke)
            
                

            self.ut.save_all(om,Elf_henke, "elf_henke_f_Nist" )

            return om, Elf_henke

