# -*- coding: utf-8 -*-
"""
Main for the Nist-ELF

Authors: PE Trevisanutto, G. Teobaldi.
"""
import classutility as cu
import class_nist as Nist
import matplotlib.pyplot as plt

Xeas, Yeas = [], []

ut = cu.utilities()

ut.letsstart()

Nist_class =Nist.Nist_interpolation(ut)

om ,elf_release = Nist_class.calc_ELF_compound_Nist()


for line in open(Nist_class.File_TDDFT,'r'):

    values = [float(s) for s in line.split()]
    
    Xeas.append(values[0]) 
    Yeas.append(values[2])

plt.loglog(om, elf_release,linestyle='dashed',color='blue',label='Nist w henke for')
plt.loglog(Xeas,Yeas,color='red',linestyle='dashed', label='TDDFT')
plt.xlabel("$\omega (eV)$")
plt.legend(loc="lower left")
plt.title("ELF")
plt.grid(True)
plt.show()

ut.sutor()