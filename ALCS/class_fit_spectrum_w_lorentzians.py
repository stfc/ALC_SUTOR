# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:03:01 2023
ALC_Sutor Project
Class for fitting the ELF spectrum with lorentians and other lmfit funcitions
Import lmfit model classes



Author:             Paolo Emilio Trevisanutto (STFC-UKRI)

"""

from lmfit import Model
from lmfit.models import LorentzianModel, SplitLorentzianModel,ThermalDistributionModel,ExponentialGaussianModel,GaussianModel,LognormalModel,BreitWignerModel,PowerLawModel
import math


class fitting_with_lorentzians():
    """
    Created on Tue May  2 17:03:01 2023
    ALC_Sutor Project
    Class for fitting the ELF spectrum with lorentians and other lmfit funcitions
    Import lmfit model classes



    Author:             Paolo Emilio Trevisanutto (STFC-UKRI)

    """
    
#--------------------------------------------
    def __init__(self,utilities):

        self.ut = utilities

#---------------------------------------------------
    def add_peak_FanoModel_Mau(self,prefix,  center, amplitude, sigma,q,al):
        """
           Method that creates parameters for the lmfit Fano Model modified 
        
         Args:
             prefix (string): Name of the Fano model
             center (float): position function
             amplitude (float): height for Fano function
             sigma(float): width
             q(float): q for asymmetry of the Fano function

           Returns:
              peaks (ndarray): lmfit peak
              params (ndarray):related paramaters for the lmfit peak
    """
        #peak = BreitWignerModel(prefix=prefix)
        peak = Model(self.ut.Fano_Mau,prefix=prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        pars = peak.make_params()
        
        pars[prefix + 'alpha'].set(1,vary=False)
        
        if (center>=10**4):
            pars[prefix + 'center'].set(center)
            pars[prefix + 'amplitude'].set(amplitude)
            pars[prefix + 'sigma'].set(sigma,min=0)
            pars[prefix + 'q'].set(q,min=0,vary=True)
            
        else:
            pars[prefix + 'center'].set(center)
            pars[prefix + 'amplitude'].set(amplitude)
            pars[prefix + 'sigma'].set(sigma,min=0)
            pars[prefix + 'q'].set(q,min=0,vary=True)
        
        
            
           
    
        
        return peak, pars
###########---------------------------------------------------
    def add_peak_FanoModel(self,prefix,  center, amplitude, sigma,q=1):
        """
       Method that creates parameters for the lmfit Fano Model modified 
    
     Args:
         prefix (string): Name of the Fano model
         center (float): position function
         amplitude (float): height for Fano function
         sigma(float): width
         q(float): q for asymmetry of the Fano function

       Returns:
          peaks (ndarray): lmfit peak
          params (ndarray):related paramaters for the lmfit peak
"""
 
        peak = BreitWignerModel(prefix=prefix)
        
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        pars = peak.make_params()
        
        
        if (center>=10**4):
            pars[prefix + 'center'].set(center)
            pars[prefix + 'amplitude'].set(amplitude)
            pars[prefix + 'sigma'].set(sigma,min=0)
            pars[prefix + 'q'].set(q,min=0,vary=True)
            
        else:
            pars[prefix + 'center'].set(center)
            pars[prefix + 'amplitude'].set(amplitude)
            pars[prefix + 'sigma'].set(sigma,min=0)
            pars[prefix + 'q'].set(q,min=0,vary=True)
            
            
           
    
        
        return peak, pars

###########---------------------------------------------------
    def add_peak_GaussianlogModel(self,prefix, center, amplitude,sigma):
        """
           Method that creates parameters for the lmfit lognorm Model
        
         Args:
             prefix (string): Name of the lmfit model
             center (float): position function
             amplitude (float): height 
             sigma(float): width
            

           Returns:
              peaks (ndarray): lmfit peak
              params (ndarray):related paramaters for the lmfit peak
    """
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        peak = LognormalModel(prefix=prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        pars = peak.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma)
    
        
        return peak, pars
###########---------------------------------------------------
    def add_peak_GaussianModel(self,prefix, amplitude, center,sigma):
        """
           Method that creates parameters for the lmfit gaussian Model
        
         Args:
             prefix (string): Name of the lmfit model
             center (float): position function
             amplitude (float): height 
             sigma(float): width
            

           Returns:
              peaks (ndarray): lmfit peak
              params (ndarray):related paramaters for the lmfit peak
    """
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        peak = GaussianModel(prefix=prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        pars = peak.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma)
    
        
        return peak, pars
    
###########---------------------------------------------------
    def add_peak_ExponentialGaussianModel(self,prefix, amplitude, center, sigma,gamma=0.1):
        ##not used
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        peak = ExponentialGaussianModel(prefix=prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        pars = peak.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma)
        pars[prefix + 'gamma'].set(gamma)
        
        return peak, pars
    
###########---------------------------------------------------    
    
#    def add_peak_Mau_F(self,prefix,  center, amplitude, sigma,A=1,kt=-1):
#    #Mau_LorentzianModel= model(Mau_Lorentizian)     
#        gmodel = Model(self.ut.Mau_Lorentizian, prefix=prefix)*ThermalDistributionModel(prefix=prefix,form='fermi')
#        print(f'parameter names: {gmodel.param_names}')
#        print(f'independent variables: {gmodel.independent_vars}')  
#        
#        gmodel.set_param_hint(prefix + 'x0', min=0)
#        #mod.set_param_hint('fwhm', expr='2.3548*wid')
#        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma)
#        pars[prefix + 'x0'].set(center)
#        pars[prefix + 'a'].set(amplitude)
#        pars[prefix + 'gam'].set(sigma, min=0)
#        pars[prefix+'amplitude'].set(A,vary=False)
#        pars[prefix+'kt'].set(kt,vary=False)
#        pars[prefix+'center'].set(center)
#        return gmodel, pars

#---------------------------------------------------
    def add_peak_Mau_h_lor(self,prefix, center, amplitude, sigma,sigma_r):
        """
           Method that creates parameters for the lmfit Asymmetric Lorentzian Model
        
         Args:
             prefix (string): Name of the lmfit model
             center (float): position function
             amplitude (float): height 
             sigma(float): width left from the center
             sigma_r(float): width right from the center
            

           Returns:
              peaks (ndarray): lmfit peak
              params (ndarray):related paramaters for the lmfit peak
    """
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Mau_Lorentizian_h, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma,gam_r=sigma_r)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma)
        pars[prefix + 'gam_r'].set(sigma_r)
    
        return gmodel, pars
#---------------------------------------------------
##not used
    def add_peak_PowerL_h(self,prefix, xin,xfin,A,c=1):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Power_law_h,prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x_in=xin, x_fin=xfin,Ampl=A,coeff=c)
        pars[prefix + 'x_in'].set(xin)
        pars[prefix + 'x_fin'].set(xfin)
        pars[prefix + 'Ampl'].set(A)
        pars[prefix + 'coeff'].set(c)
        return gmodel, pars
#---------------------------------------------------
    def add_peak_Mau(self,prefix, center, amplitude,sigma):
        """
           Method that creates parameters for a Lorentzian Model
        
         Args:
             prefix (string): Name of the lmfit model
             center (float): position function
             amplitude (float): height 
             sigma(float): width 
             
            

           Returns:
              peaks (ndarray): lmfit peak
              params (ndarray):related paramaters for the lmfit peak
    """
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Mau_Lorentizian, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
        
        pars = gmodel.make_params(prefix=prefix, mu=center, a=amplitude,gam=sigma)
        if (center > 10**3):
            pars[prefix + 'mu'].set(center, vary=True)
            pars[prefix + 'a'].set(amplitude, vary=True)
            pars[prefix + 'gam'].set(sigma, min=0)
        else:
            pars[prefix + 'mu'].set(center, vary=True)
            pars[prefix + 'a'].set(amplitude, vary=True)
            pars[prefix + 'gam'].set(sigma, min=0)
        return gmodel, pars
###########---------------------------------------------------
    def add_peak(self,prefix, center, amplitude, sigma):
        gmodel = LorentzianModel(prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
        
        pars = gmodel.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma, min=0)
        return gmodel, pars    
###########---------------------------------------------------
    def add_peak_fermi(self,prefix, center,amplitude, sigma, B):
        
        gmodel =  Model(self.ut.Mau_Lorentizian_F, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}') 
        
        #gmodel.set_param_hint(prefix + 'mu',expr=center)
        #gmodel.set_param_hint(prefix + 'a',expr=amplitude)
        #gmodel.set_param_hint(prefix + 'gam',expr=sigma)
        pars = gmodel.make_params()
        
        pars[prefix+ 'mu'].set(center, vary=True)
        pars[prefix + 'a'].set(amplitude, vary=True)
        pars[prefix + 'gam'].set(sigma, min=0)
        pars[prefix + 'B'].set(B, vary=True)
        
        return gmodel, pars
###########---------------------------------------------------
    def add_peak_h_lor(self,prefix, center,amplitude, sigma, sigma_r):
        ##Split Lorentzian!!!
        peak = SplitLorentzianModel(prefix=prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}') 
        pars = peak.make_params()
        
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
    
        pars[prefix + 'sigma'].set(sigma)
        pars[prefix + 'sigma_r'].set(sigma_r)
 
        return peak, pars
