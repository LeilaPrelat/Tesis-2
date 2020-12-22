#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:50:07 2020

@author: leila
"""

import numpy as np
import sys 
import os

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

def omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,epsi_ci,R,hbaramu):
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    epsi_ci: Im(epsilon1) del colorante (dye)
    modo : modo
    R : radio del cilindro en micrometros
    hbaramu : hbar*mu potencial quimico del grafeno in eV

    Returns
    -------
    complex frequency (pole) omega/c QE approximation
    dispersive case

    """


    TK = 300               	# temperature in K
    akb = 8.6173324E-5           ### Boltzman constant in eV/K 

    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02 = -4*pi*1j*intra*modo/R #omega0n**2

    gammaDL_ev = gamma_DL/hb
    gamma_c_ev = hbargama/hb

    omegap = Ep/hb
    num_real = omegap**2 + omega02
    den_real = epsiinf_DL + 1j*epsi_ci + epsi2
    omega_real = (num_real/den_real)**(1/2)
    num_imag = gammaDL_ev*(omegap**2) + gamma_c_ev*omega02
    den_imag = 2*(omegap**2 + omega02)
    omega_imag = num_imag/den_imag
    rta = omega_real - 1j*omega_imag    
    
    return rta/c

#%%