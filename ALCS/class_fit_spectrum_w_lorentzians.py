# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:03:01 2023
Sutor Project
Class for fitting the ELF spectrum with lorentians.
Import lmfit model classes

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

from lmfit import Model
from lmfit.models import LorentzianModel, SplitLorentzianModel,ThermalDistributionModel,ExponentialGaussianModel,GaussianModel,LognormalModel,BreitWignerModel,PowerLawModel
import math


class fitting_with_lorentzians():
    
####--------------------------------------------
    def __init__(self,utilities):

        self.ut = utilities
###########---------------------------------------------------
    def add_PowerLawModel(self,xData,yData,prefix,amplitude,exponent):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
        peak = Model(self.ut.PowerLawModel,prefix)
        print(f'parameter names: {peak.param_names}')
        print(f'independent variables: {peak.independent_vars}')  
    
        #pars = peak.guess(yData, xData)
        pars = peak.make_params()
        pars[prefix + 'amplitude'].set(amplitude)
        pars[prefix + 'exponent'].set(exponent) 
        
        return peak, pars
###########---------------------------------------------------

    def add_peak_FanoModel_Mau(self,prefix,  center, amplitude, sigma,q,al):
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
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
    #Mau_LorentzianModel= model(Mau_Lorentizian)     
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

###########---------------------------------------------------
    def add_peak_Mau_h_lor(self,prefix, center, amplitude, sigma,sigma_r):
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
###########---------------------------------------------------
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
###########---------------------------------------------------
    def add_peak_Mau(self,prefix, center, amplitude,sigma):
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
