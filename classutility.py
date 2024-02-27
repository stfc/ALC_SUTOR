# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:32:35 2023

Class utility
Authors: P.E. Trevisanutto, G. Teobaldi
"""
import matplotlib.pyplot as plt
import numpy as np
import math as math
import scipy.integrate as it
import sys
class utilities():

    __hb2m = 7.615 ##(hb*c)^2/0.5 MeV=[eV*Angs]^2
    __ao   = 0.529  #Borh radius in Angs
    #Avogadro number
    __N_avo=6.022*10**(23)
    #electronic radius
    ##BransdenAngs
    __r_e =2.8179410**(-13)

    #atomic density for atom
    __C_henke=0.911*10**(16)

    ##hc[eV*cm]
    __hc=12.398*10**(-5)
    #Vol= 5681.2801

###----------------------------------------------------------------------------    
    def __init__(self, Volume= 5681.2801):
        
            self.Vol=Volume
###----------------------------------------------------------------------------  
        
    def plotty(self,x, y, name, color='blue',scale='linear',ylabel="IMFP ($\AA$)"):
        if (scale=='loglog'):
            
            plt.loglog(x, y,color= color)
            
        elif (scale=='semi'):
            plt.semilogx(x, y,color= color)
            
        elif('linear'):
            plt.plot(x, y,color= color)
            
        plt.xlabel("$\omega (eV)$")
        plt.ylabel(ylabel)
        plt.title(name)
        plt.grid(True)
        plt.show()
#-----------------------------------------
    def taking_data (self,filename):
        Xeas = []
        Yeas = []
        values =[]
        ##
        for line in open(filename, mode='r'):
            values=line.split()
            if (self.is_number(values[0]) == False):
                pass
            else:
                Xeas.append(values[0])
                Yeas.append(values[1]) 
            values.clear()    
    
        return Xeas,Yeas  
###-------------------------------------------------------------------------
    def Integration(self,omega, imels):
     
        da_integrare=[]
        
        da_integrare = (self.Vol/(2*math.pi**2))*omega*imels
     
        sum_peff = it.simpson(da_integrare,omega,'avg')
     
        return sum_peff
###--------------------------------------------
    def plot_data2(self,Ysym_data,Yasym_data, X_data, save):

            plt.figure(figsize=(8, 6))
            plt.plot(X_data,Ysym_data ,linestyle='solid' ,marker='o', label= 'Z')
            plt.plot(X_data,Yasym_data ,linestyle='solid' ,marker='o', label= 'P')


            plt.tick_params(labelsize= 22)
            plt.grid(True)  
            plt.ylabel(str,size= 'x-large')
            plt.xlabel('$\omega$(eV)', size= 'x-large')
            plt.legend(ncol=1)
            plt.title(str)

            if save == True: plt.savefig(str+".pdf")

            plt.show()
###--------------------------------------------
    def save_all(self,x,y,name):

    
        output = '\n'.join('\t'.join(map(str,row)) for row in zip(x,y))
    
        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")
    
        f.close()
###--------------------------------------------
    def save_all2(self,x,y1,y2,name):

    
        output = '\n'.join('\t'.join(map(str,row)) for row in zip(x,y1,y2))
    
        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")
    
        f.close()
###--------------------------------------------
###hc=12398.5[eV*Angs] E=hc/lambda
    def eV2cm(self,omg):
    
        #rho=[]
    
        
            
            #rho.append(self.hc/omg[i])
        omg_np=self.From_str2float(omg)
        rho_np = np.reciprocal(omg_np)
        rho_np = self.hc*rho_np
    
        return(rho_np)
####--------------------------------------------
    def is_number(self,s):
        
        ret=True
        try:
            float(s)
        
        except ValueError:
           ret=False
           
        return ret
######-------------------------------------------
    def lorentzian(self, x,a, x0, gam ):
   # print (str(a) +" "+str(x0) +" "+str(gam) +" ")
        
        return a* gam/math.pi / (( gam)**2 + ( x - x0 )**2)
###########---------------------------------------------------
    def Mau_Lorentizian_h(self, x,x0, a, gam, gam_r ):
        Const =2*a/(gam+gam_r)
    
        f_lor_h= Const*gam*x*np.heaviside(x0,x) / (( gam*x)**2 + ( x**2 - x0**2 )**2) +\
            Const*gam_r*x*np.heaviside(x,x0) / (( gam_r*x)**2 + ( x**2 - x0**2 )**2)
        
        return f_lor_h
##############################################################
    def Mau_Lorentizian(self, x,x0, a, gam ):
        return a* gam*x / (( gam*x)**2 + ( x**2 - x0**2 )**2)
###########---------------------------------------------------
    def quadratic(self,x,a,b,c):
        return a*x**2 +b*x + c
######-------------------------------------------
    def multi_lorentz(self, x, par_q,par_l,par_f):
    
        off= self.quadratic(x,par_q[0],par_q[1],par_q[2])
    #print('off '+str(off))
    #off =0
    
        assert not ( len( par_l ) %  3)
    
    #lorentzian x0,a,gam
        multi_l = off + sum( [ self.Mau_Lorentizian( x, par_l[ i ],par_l[i+1],par_l[i+2] ) for i in range( 0, len( par_l ), 3 ) ] )
    
        assert not ( len( par_f ) %  4)
    
    #lorentzian x0,a,gam,gam_r
        multi_lf = multi_l + sum( [ self.Mau_Lorentizian_h( x, par_f[ i ],par_f[i+1],par_f[i+2],par_f[i+3] ) for i in range( 0, len( par_f ), 4 ) ] )
    #print("passed")
        return multi_lf
#-----------------------------------------
    def reading_inputfile(self):
        
        name_of_script = sys.argv[1]
        return name_of_script
#-----------------------------------------
    def From_str2float(self,str):
        rvect_f =[]
        for item in str:
            rvect_f.append(float(item))
        
        rv = np.asarray(rvect_f)
        return rv
#####----------------------------------------
    def letsstart(self):
        print ("Festina lente",flush=True)
#####----------------------------------------
    def sutor(self):
        print ("Sutor ne ultra crepidam",flush=True)