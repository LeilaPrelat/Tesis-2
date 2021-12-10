#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 10:58:12 2020

@author: leila
"""

import numpy as np
import sys
import os

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/con_kz_real/complex_freq','') 

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

def omegac_QE(im_epsi1,modo,re_epsi1,R,hbaramu):
    """
    Parameters
    ----------
    im_epsi1 : imaginary part of electric permeability of medium 1
    modo : mode
    re_epsi1 : real part of electric permeability of medium 1
    R : radius of cylinder in micrometer
    hbaramu : hbar*mu chemical potential of graphene in eV

    Returns
    -------
    complex frequency (pole) omega/c QE approximation
    non-dispersive case

    """


    TK = 300               	# temperature in K
    akb = 8.6173324E-5           ### Boltzman constant in eV/K 

    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02n = -4*pi*1j*intra*modo/R #omega0n**2
    gamma = hbargama/hb
    
    # term1 = omega02n/(epsi1 + epsi2)
    # term = (term1 - (gamma/2)**2)**(1/2) - 1j*gamma/2
    epsi1 = re_epsi1 + 1j*im_epsi1
    cte = epsi1 + epsi2
    term = (omega02n/cte)**(1/2) - 1j*gamma/2
    
    return term/c