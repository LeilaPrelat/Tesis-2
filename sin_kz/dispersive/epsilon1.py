#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 10:30:50 2020

@author: leila

permeabilidad electrica del medio activo + nanocristal en el interior del cilindro
"""

import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/dispersive','')

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

def epsilon1(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    epsi_ci: Im(epsilon1) del colorante (dye)
    
    Returns
        epsilon1, permeabilidad electrica del medio 1 
            con nanocristal + colorante
    -------
    """
    # omega = omegac*c
    energy = omegac*c*hb

    num = Ep**2
    den = energy**2 + 1j*gamma_DL*energy
    epsi1 = epsiinf_DL - num/den + 1j*epsi_ci
    
    return epsi1

#%%
