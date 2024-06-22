# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:32:35 2023
Sutor Project
Class utility

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
import sys
from scipy.special import expit
from scipy.signal import find_peaks

class utilities():

    __hb2m = 7.615  # (hb*c)^2/0.5 MeV=[eV*Angs]^2
    __ao = 0.529  # Borh radius in Angs
    # Avogadro number
    __N_avo = 6.022*10**(23)
    # electronic radius
    # BransdenAngs
    __r_e = 2.8179410**(-13)

    # atomic density for atom
    __C_henke = 0.911*10**(16)

    # hc[eV*cm]
    __hc = 12.398*10**(-5)
    #Vol= 5681.2801

    def plotting_elf(self, xData, yData, peaks):

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.loglog(xData, yData)
        ax.loglog(xData[peaks], yData[peaks], "x")
        ax.set_xlabel('$\omega$ (eV)')
        ax.set_ylabel('ELF')

        plt.title('Initial guess')
        plt.show()
# ----------------------------------------------------------------------------

    def resultbestfitplot(self, x, y, init, result,name):

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        #ax.loglog(x, init, '--', label='initial fit')
        ax.loglog(x, y, label='data')
        ax.loglog(x, result, label='best fit')
        ax.set_xlabel('$\omega$ (eV)')
        ax.set_ylabel('ELF')
        plt.show()
        plt.savefig(name+'_figure.pdf')
        #fig = plt.figure()
        #ax = fig.add_subplot( 1, 1, 1 )
        #ax.set_xlabel('$\omega$ (eV)')
        # ax.set_ylabel('ELF')
        #ax.loglog(x,y, linewidth=3,color='red',label='data')
        #ax.loglog(xData, result.best_fit, label='best fit')

# ----------------------------------------------------------------------------
    def __init__(self, Volume=5681.2801):

        self.Vol = Volume
# ----------------------------------------------------------------------------

    def get_Avogadro(self):

        return self.__N_avo
# ----------------------------------------------------------------------------
    def get_hc(self):
        return self.__hc
    
    # ----------------------------------------------------------------------------
    def get_henke(self):
            return self.__C_henke
# ----------------------------------------------------------------------------

    def index_of(self, arrval, value):
        """Return index of array *at or below* value."""

        
        if value < min(arrval):
            return 0
        return max(np.where(arrval <= value)[0])

# ----------------------------------------------------------------------------

    def plotty(self, x, y, name, color='blue', scale='linear', ylabel="IMFP ($\AA$)",save=True):
        if (scale == 'loglog'):

            plt.loglog(x, y, color=color)

        elif (scale == 'semi'):
            plt.semilogx(x, y, color=color)

        elif('linear'):
            plt.plot(x, y, color=color)

        plt.xlabel("$\omega (eV)$")
        plt.ylabel(ylabel)
        plt.title(name)
        plt.grid(True)
        plt.show()
        if save == True:
            plt.savefig(name + ".pdf")
        
# -----------------------------------------
# ----------------------------------------- 
    def taking_initialguess(self,xData,yData,height,Guessing):
        
        guess=[]
        Amp =[]
        gamma =[]
        Inte = []
        print (Guessing)
        if (Guessing==True):
            print ("Looking for the peaks")
            peaks, properties = find_peaks( yData,width=0.000000001,height=height)

       
            I = properties [ 'peak_heights' ]
            fwhm =properties["width_heights"]

            for i in range( len( peaks) ): 

                guess.append( xData[peaks[i]] )

                Inte.append(I[i])
     
                gamma.append(fwhm[i]/2)
           


            rough_peak_positions = self.From_str2float(guess)
            Amp = self.From_str2float(Inte)
            gamm =self.From_str2float(gamma)

            print (rough_peak_positions)
            with open("peaks.dat", 'w') as ffa2:
                for i in range (len(guess)):
                
                    ffa2.write(str(guess[i])+' ')
                    ffa2.write(str(Amp[i])+' ')
                    ffa2.write(str(gamm[i])+' ')
                    ffa2.write(str(peaks[i])+ ' ')
                    ffa2.write("\n")

            ffa2.close()
            peas = self.From_str2Integer(peaks)
        else:
            print ("Reading file peaks.dat")
            peaks =[]
            for line in open("peaks.dat",  mode='r'):
        
          
                if (line != "\n"):
                    values = line.split()

                    if (self.is_number(values[0]) == False):
                        pass
                    else:
                        guess.append(values[0])
                        Inte.append(values[1])
                        gamma.append(values[2])
                        peaks.append(values[3])   
                    values.clear()
                    
            rough_peak_positions = self.From_str2float(guess)
            Amp = self.From_str2float(Inte)
            gamm =self.From_str2float(gamma)
            peas = self.From_str2Integer(peaks)
            print( peas)
        return peas, rough_peak_positions, Amp, gamm
            
# -----------------------------------------

    def taking_data(self, filename):
        Xeas = []
        Yeas = []
        values = []
        ##
        for line in open(filename, mode='r'):
            if (line != "\n"):
                values = line.split()

                if (self.is_number(values[0]) == False):
                    pass
                else:
                    Xeas.append(values[0])
                    Yeas.append(values[1])
                values.clear()

        return Xeas, Yeas
# -------------------------------------------------------------------------

    def Integration(self, omega, imels):

        da_integrare = []

        da_integrare = (self.Vol/(2*math.pi**2))*omega*imels

        sum_peff = it.simpson(da_integrare, omega, 'avg')

        return sum_peff
# --------------------------------------------

    def plot_data2(self, Ysym_data, Yasym_data, X_data, save):

        plt.figure(figsize=(8, 6))
        plt.plot(X_data, Ysym_data, linestyle='solid', marker='o', label='Z')
        plt.plot(X_data, Yasym_data, linestyle='solid', marker='o', label='P')

        plt.tick_params(labelsize=22)
        plt.grid(True)
        plt.ylabel(str, size='x-large')
        plt.xlabel('$\omega$(eV)', size='x-large')
        plt.legend(ncol=1)
        plt.title(str)

        if save == True:
            plt.savefig(str+".pdf")

        plt.show()
# --------------------------------------------


    def save_pinel(self, x, y, name):

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y))

        with open(name+".txt", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------
    def save_all(self, x, y, name):

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y))

        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------

    def save_all2(self, x, y1, y2, name):

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y1, y2))

        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------
# hc=12398.5[eV*Angs] E=hc/lambda

    def eV2cm(self, omg):

        # rho=[]

        # rho.append(self.hc/omg[i])
        omg_np = self.From_str2float(omg)
        rho_np = np.reciprocal(omg_np)
        rho_np = self.get_hc()*rho_np

        return(rho_np)
# --------------------------------------------

    def is_number(self, s):

        ret = True
        try:
            float(s)

        except ValueError:
            ret = False

        return ret
# -------------------------------------------
    def Polynomial(self, x, coeff):

        y = 0

        y = sum(x**(i)*coeff[i] for i in range(0, len(coeff)))

        return y

# ---------------------------------------------------    
    def lorentzian(self, x, a, x0, gam):
       # print (str(a) +" "+str(x0) +" "+str(gam) +" ")

        return a * gam/math.pi / ((gam)**2 + (x - x0)**2)


# ---------------------------------------------------
    def Power_law(self,x, par):
        
        y =par[0]*x**par[1] 

        return y

# ---------------------------------------------------
    def Power_law_h(self,x, x_in,x_fin,Ampl,coeff):
        
        y =np.heaviside(x-x_in, 0.5)*np.heaviside(x_fin-x, 0.5)*Ampl*x**coeff 

        return y      
####------------------------------------------------------------------

    def Mau_Lorentizian_h(self, x, x0, a, gam, gam_r):
        Const = 2*a/(gam+gam_r)

        f_lor_h = Const*gam*x*np.heaviside( x-x0,0.5) / ((gam*x)**2 + (x**2 - x0**2)**2) +\
            Const*gam_r*x*np.heaviside(x0-x, 0.5) / \
            ((gam_r*x)**2 + (x**2 - x0**2)**2)

        return f_lor_h
##############################################################

    def Mau_Lorentizian(self, x, mu, a, gam):
        return a * gam*x / ((gam*x)**2 + (x**2 - mu**2)**2)
##############################################################

    def Mau_Lorentizian_F(self, x, mu, a, gam,B):
         
        #s = -1.*(x-B)
        #F = 1.0/(expit(s) + 1.)
        F = np.heaviside(x-B,0.5)
        y =  a * gam*x / ((gam*x)**2 + (x**2 - mu**2)**2)
        #y*np.heaviside(mu-x, 0.5)
        return y*F
##############################################################

    def Fano_function(self, x,amplitude,center, sigma,q):

        f = amplitude*(q*sigma/2 + x - center)**2/((sigma/2)**2+(x-center)**2)

        return f
#############################################################

    def Fano_Mau(self, x,amplitude,center, sigma,q,alpha):
        #alpha =0
        if (alpha == 1):
            beta = 2
        elif (alpha == 0):
            beta =1

        f = amplitude*(q*sigma*x/2 +alpha*(x - center))**beta/((sigma*x/2)**2+(x-center)**2)

        return f
##############################################################

    def lognormal_function(self, x, A, mu, sigma):

        f = (A/(x*sigma*(2*math.pi)**(0.5))) * \
            math.exp((-1.0/(2*sigma**2))*(math.log(x) - mu)**2)

        return f
##############################################################
# ---------------------------------------------------

    def quadratic(self, x, a, b, c):
        return a*x**2 + b*x + c
##############################################################
# -------------------------------------------

    def multi_lorentz(self, x, par_q, par_pl, par_l, par_f, approx):

        off = self.quadratic(x, par_q[0], par_q[1], par_q[2])
 
        
    # lorentzian x0,a,gam
        multi_l = off + sum([self.Mau_Lorentizian(x, par_l[i],
                            par_l[i+1], par_l[i+2]) for i in range(0, len(par_l), 3)])

        if (approx == "Fermi"):
            
            assert not (len(par_f) % 4)

    # lorentzian x0,a,gam,gam_r
#            multi_lf = multi_l + sum([self.Mau_Lorentizian_h(
#                x, par_f[i], par_f[i+1], par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])
            multi_l = multi_l + sum([self.Mau_Lorentizian_F(x, par_f[i], par_f[i+1], par_f[i+2], par_f[i+3])
                                          for i in range(0, len(par_f), 4)])
        elif (approx == "Fano"):

            assert not (len(par_f) % 5)

            multi_l = multi_l + sum([self.Fano_Mau(x, par_f[i], par_f[i+1],
                                     par_f[i+2], par_f[i+3],par_f[i+4]) for i in range(0, len(par_f), 5)])

        elif (approx == "LN"):

            assert not (len(par_f) % 3)
            
            multi_l = multi_l + sum([self.lognormal_function(
                x, par_f[i], par_f[i+1], par_f[i+2]) for i in range(0, len(par_f), 3)])
        
            
        
            #self.sutor_issue()
    # print("passed")
        return multi_l
# -------------------------------------------   
    def powerlaw_fermi(self, x, par_pl, par_f):

       off = self.Power_law(x, par_pl)
  

       assert not (len(par_f) % 4)

       multi_pf = off + sum([self.Mau_Lorentizian_F(x, par_f[i], par_f[i+1],
                                    par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])

       return multi_pf
# -------------------------------------------   
    def powerlaw_lorentz(self, x, par_pl, par_f):

       off = self.Power_law(x, par_pl)
  

       assert not (len(par_f) % 3)

       multi_pf = off + sum([self.Mau_Lorentizian(x, par_f[i], par_f[i+1],
                                    par_f[i+2]) for i in range(0, len(par_f), 3)])

       return multi_pf
 # -------------------------------------------   
    def powerlaw_fano(self, x, par_pl, par_f):

        off = self.Power_law(x, par_pl)
   

        assert not (len(par_f) % 5)

        #multi_pf = off + sum([self.Fano_function_Mau(x, par_f[i], par_f[i+1],
        #                             par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])
        multi_pf = off + sum([self.Fano_Mau(x, par_f[i], par_f[i+1],
                                     par_f[i+2], par_f[i+3],par_f[i+4]) for i in range(0, len(par_f), 5)])
        return multi_pf
# -------------------------------------------   
    def powerlaw_ln(self, x, par_pl, par_f):

       off = self.Power_law(x, par_pl)
  
       print (len(par_f))
       assert not (len(par_f) % 3)

       multi_pf = off + sum([self.lognormal_function(x, par_f[i], par_f[i+1],
                                    par_f[i+2]) for i in range(0, len(par_f), 3)])

       return multi_pf
# -------------------------------------------
# Not used.

    def multi_lorentz_h(self, x, par_f):

        assert not (len(par_f) % 4)

    # lorentzian x0,a,gam,gam_r
        multi_lf = sum([self.Mau_Lorentizian_h(x, par_f[i], par_f[i+1],
                       par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])
    # print("passed")
        return multi_lf
# -------------------------------------------
# Not used.

    def multi_Fano(self, x, par_fano):

        assert not (len(par_fano) % 4)

    # lorentzian x0,a,gam,gam_r
        multi_lf = sum([self.Fano_function(x, par_fano[i], par_fano[i+1],
                       par_fano[i+2], par_fano[i+3]) for i in range(0, len(par_fano), 4)])
    # print("passed")
        return multi_lf
# -----------------------------------------

    def reading_inputfile(self):

        name_of_script = sys.argv[1]
        return name_of_script
# -----------------------------------------

    def From_str2float(self, str):
        rvect_f = []
        for item in str:
            rvect_f.append(float(item))

        rv = np.asarray(rvect_f)
        return rv
# -----------------------------------------

    def From_str2Integer(self, str):
        rvect_f = []
        for item in str:
            rvect_f.append(int(item))

        rv = np.asarray(rvect_f)
        return rv
# ----------------------------------------

    def letsstart(self):
        print("Festina lente.\n", flush=True)
# ----------------------------------------

    def sutor(self):
        print("Sutor ne ultra crepidam.", flush=True)
# ----------------------------------------

    def sutor_issue(self, cissue):
        print(cissue, flush=True)
        print("Morituro te salutat", flush=True)
        sys.exit("failure")
