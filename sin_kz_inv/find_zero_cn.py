#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con roots_of_gn.py:
misma formula que gn pero 
adimensionalizada
    
(determinante de 4x4 sin kz)

en lugar de minimizar el denominador de los coef
an y cn (caso sin kz no dispersivo) al que llamamos
gn, ahora vamos a minimizar al numerador de an
(cn es el coef de Hz1, medio interior del cilindro)
 

"""

# import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz_inv','') 

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

def determinante(omegac,epsi1,mode,R,mu_c):
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio1
    mode : modo
    R : radio del cilindro en micrometros
    mu_c : potencial quimico del grafeno
    Ao : coef del campo incidente

    Returns
    -------

    en lugar de minimizar el denominador de los coef
    an y cn (caso sin kz no dispersivo) al que llamamos
    gn, ahora vamos a minimizar al numerador de cn (medio 1)

    sol no trivial: Ao =! 0

    """
    Ao = 1
    mode = int(mode)
    
    k0 = omegac
    Rbarra = R*k0 #adimensional
    
    x1,x2 = (epsi1*mu1)**(1/2), (epsi2*mu2)**(1/2)
    x1t,x2t = x1, x2
    
    if (x1t*Rbarra).real>=0:
 	    x1t = x1t
    else:
 	    x1t = -x1t
    
    if (x2t*Rbarra).real>=0:
	    x2t = x2t
    else:
	    x2t = -x2t
    
    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    # Sigma_ad = 4*pi*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    
    def Bessel(modo):
        J2 =  special.jn(modo,x2t*Rbarra)+0j
        derJ2 =  special.jvp(modo,x2t*Rbarra)+0j
        H2 =  special.hankel1(modo,x2t*Rbarra)+0j
        derH2 =  special.h1vp(modo,x2t*Rbarra)+0j
        return [J2,derJ2,H2,derH2]
    
    [J2,derJ2,H2,derH2] = Bessel(mode)    

    det_ad = derJ2*H2 - derH2*J2
    
    # print((x1t*Rbarra/1j).real>=0,(x2t*Rbarra/1j).real>=0)
    
    cte_aux = x2*epsi1*Ao*(1j**mode)
    return det_ad*cte_aux

#%%
