#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con roots_of_gn.py:
misma formula que gn pero 
adimensionalizada
    
(determinante de 4x4 sin kz)

"""

# import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/dispersive','')

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene)
    path_graphene2 = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene2)
    from graphene_sigma import sigma

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

def determinante(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,mu_c):
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
    mu_c : potencial quimico del grafeno

    Returns
    -------
    1 bloque del determinante de la matriz de 4x4
	---> caso sin kz 

    """
    energy = omegac*(c*hb)
    
    sigmatot, inter, intra = sigma(energy,mu_c,hbargama) 

    num = Ep**2
    den = energy**2 + 1j*gamma_DL*energy
    epsi1 = epsiinf_DL - num/den + 1j*epsi_ci
    
    k0 = omegac
    Rbarra = R*k0 #adimensional
    
    x1,x2 = (epsi1*mu1)**(1/2), (epsi2*mu2)**(1/2)
    x1t,x2t = x1, x2
    
    # if (x1t*Rbarra).real>=0:
 	  #   x1t = x1t
    # else:
 	  #   x1t = -x1t
    
    if (x2t*Rbarra).real>=0:
	    x2t = x2t
    else:
	    x2t = -x2t
    
    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    Sigma_ad = 4*pi*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    
    def Bessel(mode):
        J = special.jn(mode,x1t*Rbarra)+0j
        derJ = special.jvp(mode,x1t*Rbarra)+0j
        H =  special.hankel1(mode,x2t*Rbarra)+0j
        derH =  special.h1vp(mode,x2t*Rbarra)+0j
        return [J,derJ,H,derH]
    
    [J,derJ,H,derH] = Bessel(modo)    

    det_ad = -epsi1*x2*J*derH + epsi2*x1*H*derJ - Sigma_ad*1j*x1*x2*derH*derJ
    
    # print((x1t*Rbarra/1j).real>=0,(x2t*Rbarra/1j).real>=0)
    
    return det_ad

#%%
