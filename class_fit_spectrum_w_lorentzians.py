# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:03:01 2023
Sutor Project
Class for fitting the ELF spectrum with lorentians.
Import lmfit model classes

Authors: P.E. Trevisanutto, G. Teobaldi.

"""

from lmfit import Model
from lmfit.models import LorentzianModel, SplitLorentzianModel,ThermalDistributionModel


class fitting_with_lorentzians():
    
####--------------------------------------------
    def __init__(self,utilities):

        self.ut = utilities
###########---------------------------------------------------
    def add_peak_Mau_F(self,prefix,  center, amplitude, sigma,x0F,A=1,kt=-1):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Mau_Lorentizian, prefix=prefix)*ThermalDistributionModel(prefix=prefix,form='fermi')
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        pars[prefix+'amplitude'].set(A,vary=False)
        pars[prefix+'kt'].set(kt,vary=False)
        pars[prefix+'center'].set(x0F)
        return gmodel, pars

###########---------------------------------------------------
    def add_peak_Mau_h_lor(self,prefix, center, amplitude, sigma,sigma_r):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Mau_Lorentizian_h, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma,gam_r=sigma_r)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        pars[prefix + 'gam_r'].set(sigma_r, min=0)
    
        return gmodel, pars
###########---------------------------------------------------
###########---------------------------------------------------
    def add_peak_Mau(self,prefix, center, amplitude,sigma):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        gmodel = Model(self.ut.Mau_Lorentizian, prefix=prefix)
        print(f'parameter names: {gmodel.param_names}')
        print(f'independent variables: {gmodel.independent_vars}')  
    
        pars = gmodel.make_params(prefix=prefix, x0=center, a=amplitude,gam=sigma)
        pars[prefix + 'x0'].set(center)
        pars[prefix + 'a'].set(amplitude)
        pars[prefix + 'gam'].set(sigma, min=0)
        return gmodel, pars
###########---------------------------------------------------
    def add_peak(self,prefix, center, amplitude, sigma):
        peak = LorentzianModel(prefix=prefix)
        pars = peak.make_params()
        pars[prefix + 'center'].set(center)
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'sigma'].set(sigma, min=0)
        return peak, pars    
###########---------------------------------------------------
    def add_peak_fermi(self,prefix_Lor,prefix_fermi, center,amplitude, sigma, x0,A=1,kt=-1):
        peak = LorentzianModel(prefix=prefix_Lor)*ThermalDistributionModel(prefix=prefix_fermi,form='fermi')

        pars = peak.make_params()
        pars[prefix_Lor + 'center'].set(center)
        pars[prefix_Lor + 'amplitude'].set(amplitude)
        
        pars[prefix_Lor + 'sigma'].set(sigma, min=0)
        pars[prefix_fermi+'amplitude'].set(A)
        pars[prefix_fermi+'kt'].set(kt)
        pars[prefix_fermi+'center'].set(x0)
        return peak, pars
####--------------------------------------------
###########---------------------------------------------------
    def add_peak_h_lor(self,prefix_Lor,prefix_fermi, center,amplitude, sigma, sigma_r):
 
        peak = SplitLorentzianModel(prefix=prefix_Lor)
        pars = peak.make_params()
        pars[prefix_Lor + 'center'].set(center)
        pars[prefix_Lor + 'amplitude'].set(amplitude, min=0)
    
        pars[prefix_Lor + 'sigma'].set(sigma, min=0)
        pars[prefix_Lor + 'sigma_r'].set(sigma_r, min=0)
 
        return peak, pars
