# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 11:32:35 2023
ALC Sutor Project
Class utility
 class with several plotting , maths and lmfit methods

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
# ----------------------------------------------------------------------------
    def plotting_elf(self, xData, yData, peaks):
        """
          method that produces the ELF plot 
           Args:
              xData (float): Energy (eV)
              yData (float): ELF 

           Returns:
              plot
       """

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.loglog(xData, yData)
        ax.loglog(xData[peaks], yData[peaks], "x")
        ax.set_xlabel('$\omega$ (eV)')
        ax.set_ylabel('ELF')

        plt.title('Initial guess')
        plt.show()
# ----------------------------------------------------------------------------
#not used
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
        """
          method get Avogadro value 
           
       """

        return self.__N_avo
# ----------------------------------------------------------------------------
    def get_hc(self):
        """
        method get Planck constant X Light speed value
        """
        return self.__hc
    
    # ----------------------------------------------------------------------------
    def get_henke(self):
        """
            method get  constant for Henke calculations
        """
        return self.__C_henke
# ----------------------------------------------------------------------------

    def index_of(self, arrval, value):
        """
           Method that Returns the index of an array *at or below* value.
        
         Args:
              arrval (numpy array): 
              value (float): chosen value 

           Returns:
              index of array *at or below* value
         
         
         """
      

        
        if value < min(arrval):
            return 0
        return max(np.where(arrval <= value)[0])

# ----------------------------------------------------------------------------

    def plotty(self, x, y, name, color='blue', scale='linear', ylabel="IMFP ($\AA$)",save=True):
        """
          General method for plotting 
           Args:
              x (float): abscissa
              y (float): ordinate 
              color (string): color of line 
              scale (string):  linear, semilog, loglog
              ylabel (string): label for ordinate
              save (boolean): true if file .dat saved

           Returns:
              plot and file .dat
       """
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
        """
        Method for the determination of the initial guess for the lmfit (Lorentzian) parameter
        
          Args:
             xData (float): energy (eV)
             yData (float): ELF
             height (float): height sensitivity to individuate peaks
             Guessing (boolean): True lmfit find_peaks looks for parameter peaks, False it takes the parameter from file peaks.dat

          Returns:
             peas (integer): position in the ELF file
             rough_peak_positions (float) : guessed peak position
             Amp (float) : guessed peak heights
             gamm (float): guessed peak widths
      """
        
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
        """
         Method for storing abscissa and ordinate in two lists
         
           Args:
              filename (string): file 

           Returns:
              Xeas (list)
              Yeas  (list)
       """
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
        """
         Method for ELF integration 
         
           Args:
              filename (string): file 

           Returns:
              Xeas (list)
              Yeas  (list)
       """

        da_integrare = []

        da_integrare = (self.Vol/(2*math.pi**2))*omega*imels

        sum_peff = it.simpson(da_integrare, omega, 'avg')

        return sum_peff
# --------------------------------------------

    def plot_data2(self, Ysym_data, Yasym_data, X_data, save):
        """
          General method for plotting three columns
           Args:
              Ysym_data  (float): First ordinate 
              Yasym_data (float): Second ordinate 
              X_data     (float): abscissa
              
              save (boolean): true if file .dat saved

         
       """

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
        """
          method for saving Cumulative Probabilities files with name (two columns)
           Args:
              x  (float):  abscissa
              y (float):   ordinate
              name (string): Name of saved file
              
         
       """

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y))

        with open(name+".txt", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------
    def save_all(self, x, y, name):
        """
          method for saving files with name (two columns)
           Args:
              x  (float):  abscissa
              y (float):   ordinate
              name (string): Name of saved file
              
         
       """

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y))

        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------

    def save_all2(self, x, y1, y2, name):
        """
         method for saving files with name (three columns)
          Args:
             x  (float):  abscissa
             y1 (float):  first ordinate
             y2 (float):  second ordinate
             name (string): Name of saved file           
        
      """

        output = '\n'.join('\t'.join(map(str, row)) for row in zip(x, y1, y2))

        with open(name+".dat", 'w') as f:
            f.write(output)
            f.write("\n")

        f.close()
# --------------------------------------------
# hc=12398.5[eV*Angs] E=hc/lambda

    def eV2cm(self, omg):
        
        """
          method transforming eV to cm-1
          
          Args:
           omg (float): frequency (eV)
           
          Returns:
              rho_np (float)= frequency (cm-1)
       """

        # rho=[]

        # rho.append(self.hc/omg[i])
        omg_np = self.From_str2float(omg)
        rho_np = np.reciprocal(omg_np)
        rho_np = self.get_hc()*rho_np

        return(rho_np)
# --------------------------------------------

    def is_number(self, s):
        """
          method to detect if the variable is a number
          Args:
           omg (float): frequency (eV)
           
          Returns:
              rho_np (float)= frequency (cm-1)
       """

        ret = True
        try:
            float(s)

        except ValueError:
            ret = False

        return ret
# -------------------------------------------
    def Polynomial(self, x, coeff):
        """
          method to Reconstruct Lmfit Polynomial function
          Args:
           x (float):  x-value
           coeff (list): list of coefficients
           
          Returns:
              y (float)= Polynomial y value
       """

        y = 0

        y = sum(x**(i)*coeff[i] for i in range(0, len(coeff)))

        return y

# ---------------------------------------------------    
    def lorentzian(self, x, a, x0, gam):
        
        """
          method to Reconstruct Lmfit Lorentzian function
          Args:
           x  (float):  x-value
           a  (float):  Amplitude
           x0  (float): Position
           gam (float): width
           
          Returns:
              y (float)= Lorentzian y value
       """
       # print (str(a) +" "+str(x0) +" "+str(gam) +" ")

        return a * gam/math.pi / ((gam)**2 + (x - x0)**2)


# ---------------------------------------------------
    def Power_law(self,x, par):
        
        """
          method to Reconstruct Lmfit Power law function
          Args:
           x  (float):  x-value
           par (list):  list of coefficients
           
           
          Returns:
              y (float)= Lorentzian y value
       """
        y =par[0]*x**par[1] 

        return y

# ---------------------------------------------------
    def Power_law_h(self,x, x_in,x_fin,Ampl,coeff):
        
        ## not used 
        y =np.heaviside(x-x_in, 0.5)*np.heaviside(x_fin-x, 0.5)*Ampl*x**coeff 

        return y      
#------------------------------------------------------------------

    def Mau_Lorentizian_h(self, x, x0, a, gam, gam_r):
        """
          method to Reconstruct Asymmetric Lorentzian
          Args:
           x  (float):  x-value
           a  (float):  Amplitude
           x0  (float): Position
           gam (float): width
           gam_r (float): asymmetric width 
           
          Returns:
              y (float)= Lorentzian y value
       """ 
        Const = 2*a/(gam+gam_r)

        f_lor_h = Const*gam*x*np.heaviside( x-x0,0.5) / ((gam*x)**2 + (x**2 - x0**2)**2) +\
            Const*gam_r*x*np.heaviside(x0-x, 0.5) / \
            ((gam_r*x)**2 + (x**2 - x0**2)**2)

        return f_lor_h
# ---------------------------------------------------

    def Mau_Lorentizian(self, x, mu, a, gam):
        """
          method to Reconstruct Asymmetric Lorentzian
          Args:
           x  (float):  x-value
           mu (float): Position
           a  (float):  Amplitude
           gam (float): width
                   
          Returns:
              y (float)= Lorentzian y value
       """ 
        return a * gam*x / ((gam*x)**2 + (x**2 - mu**2)**2)
# ---------------------------------------------------

    def Mau_Lorentizian_F(self, x, mu, a, gam,B):
         
        #s = -1.*(x-B)
        #F = 1.0/(expit(s) + 1.)
        F = np.heaviside(x-B,0.5)
        y =  a * gam*x / ((gam*x)**2 + (x**2 - mu**2)**2)
        #y*np.heaviside(mu-x, 0.5)
        return y*F
# ---------------------------------------------------
    def Fano_function(self, x,amplitude,center, sigma,q):
        """
          method to Reconstruct Lmfit Fano function
          
          Args:
           x  (float):  x-value
           center (float):  position
           sigma (float):   width
           q (float):    q asymmetric parameter
           
          Returns:
              y (float)= Fano function y value
       """

        f = amplitude*(q*sigma/2 + x - center)**2/((sigma/2)**2+(x-center)**2)

        return f
# ---------------------------------------------------

    def Fano_Mau(self, x,amplitude,center, sigma,q,alpha):
        
        #alpha =0
        if (alpha == 1):
            beta = 2
        elif (alpha == 0):
            beta =1

        f = amplitude*(q*sigma*x/2 +alpha*(x - center))**beta/((sigma*x/2)**2+(x-center)**2)

        return f
# ---------------------------------------------------

    def lognormal_function(self, x, A, mu, sigma):

        f = (A/(x*sigma*(2*math.pi)**(0.5))) * \
            math.exp((-1.0/(2*sigma**2))*(math.log(x) - mu)**2)

        return f
# ---------------------------------------------------

    def quadratic(self, x, a, b, c):
        """
          method to Reconstruct quadratic lmfit function
          Args:
           x  (float):  x-value
           a (float):  x^2 coefficient
           b  (float):  x coefficient
           c (float): Constant
                   
           
          Returns:
              y (float)= quadratic y value
       """
        return a*x**2 + b*x + c
# -------------------------------------------

    def multi_lorentz(self, x, par_q, par_pl, par_l, par_f, approx):
        """
          method to sum the Lorentzians functions
          Args:
           x  (float):  x-value
           par_q (list): parameters for quadratic function
           par_pl  (list):parameters for polynomial function (not used)
           par_l (list): parameters for Lorentzians
           par_f (list): parameters for several different functions (Fermi_lorentzian, Fano, LogNormal)
           approx (string): approximation ("Fermi_lorentz", "Fano", "LN")
           
           
          Returns:
              y (float)= reconstructed spectrum
       """
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
       """
          method to sum the fermi_lorentz with lmfit powerlaw function
          Args:
           x  (float):  x-value  
           par_pl  (list): powerlaw  parameters
           par_f (list): Fermi_lorentzian parameters

          Returns:
              y (float)= quadratic y value
       """

       off = self.Power_law(x, par_pl)
  

       assert not (len(par_f) % 4)

       multi_pf = off + sum([self.Mau_Lorentizian_F(x, par_f[i], par_f[i+1],
                                    par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])

       return multi_pf
# -------------------------------------------   
    def powerlaw_lorentz(self, x, par_pl, par_f):
       """
          method to sum the lorentzian with lmfit powerlaw function
          Args:
           x  (float):  x-value  
           par_pl  (list): powerlaw  parameters
           par_f (list): lorentzian parameters

          Returns:
              y (float)= quadratic y value
       """
        
       off = self.Power_law(x, par_pl)
  

       assert not (len(par_f) % 3)

       multi_pf = off + sum([self.Mau_Lorentizian(x, par_f[i], par_f[i+1],
                                    par_f[i+2]) for i in range(0, len(par_f), 3)])

       return multi_pf
 # -------------------------------------------   
    def powerlaw_fano(self, x, par_pl, par_f):
        """
          method to sum the Fano with lmfit powerlaw function
          Args:
           x  (float):  x-value  
           par_pl  (list): powerlaw  parameters
           par_f (list): Fano parameters

          Returns:
              y (float)= quadratic y value
       """
        off = self.Power_law(x, par_pl)
   

        assert not (len(par_f) % 5)

        #multi_pf = off + sum([self.Fano_function_Mau(x, par_f[i], par_f[i+1],
        #                             par_f[i+2], par_f[i+3]) for i in range(0, len(par_f), 4)])
        multi_pf = off + sum([self.Fano_Mau(x, par_f[i], par_f[i+1],
                                     par_f[i+2], par_f[i+3],par_f[i+4]) for i in range(0, len(par_f), 5)])
        return multi_pf
# -------------------------------------------   
    def powerlaw_ln(self, x, par_pl, par_f):
       """
          method to sum the lognorm with lmfit powerlaw function
          Args:
           x  (float):  x-value  
           par_pl  (list): powerlaw  parameters
           par_f (list): lognorm parameters

          Returns:
              y (float)= quadratic y value
       """

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
        ##not used
        

        name_of_script = sys.argv[1]
        return name_of_script
# -----------------------------------------

    def From_str2float(self, str):
        """
           method transforms string to floats
           Args:
            str  (string):  string to transform

           Returns:
               rv (float)= float value
        """
        
        rvect_f = []
        for item in str:
            rvect_f.append(float(item))

        rv = np.asarray(rvect_f)
        return rv
# -----------------------------------------

    def From_str2Integer(self, str):
        """
           method transforms string to integer
           Args:
            str  (string):  string to transform

           Returns:
               rv (float)= integer value
        """
        rvect_f = []
        for item in str:
            rvect_f.append(int(item))

        rv = np.asarray(rvect_f)
        return rv
# ----------------------------------------

    def letsstart(self):
        """
           method let's start
    """
        print("Festina lente.\n", flush=True)
# ----------------------------------------

    def sutor(self):
        """
           method It is over
    """
        print("Sutor ne ultra crepidam.", flush=True)
# ----------------------------------------

    def sutor_issue(self, cissue):
        """
           method to show an error
    """
        print(cissue, flush=True)
        print("Morituro te salutat", flush=True)
        sys.exit("failure")
