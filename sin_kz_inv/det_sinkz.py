#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

en lugar de minimizar el denominador de los coef
an y cn (caso sin kz no dispersivo) al que llamamos
gn, ahora vamos a minimizar al numerador de an
(an es el coef de Hz1, medio interior del cilindro)

"""

# import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz_inv/non-dispersive','') 

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene2 = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene2)
    from graphene_sigma import sigma

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

def determinante(omegac,epsi1,mode,R,mu_c,Ao):
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
    gn, ahora vamos a minimizar al numerador de an
    (an es el coef de Hz1, medio interior del cilindro)

    sol no trivial: Ao =! 0

    """
    mode = int(mode)
    energy = omegac*c*hb
    
    sigmatot, inter, intra = sigma(energy,mu_c,hbargama) 
    
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
    Sigma_ad = 4*pi*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    
    def Bessel(modo):
        J = special.jn(modo,x1t*Rbarra)+0j
        derJ = special.jvp(modo,x1t*Rbarra)+0j
        J2 =  special.jn(modo,x2t*Rbarra)+0j
        derJ2 =  special.jvp(modo,x2t*Rbarra)+0j
        return [J1,derJ1,J2,derJ2]
    
    [J1,derJ1,J2,derJ2] = Bessel(mode)    

    det_ad = epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*1j*x1*x2*derJ1*derJ2
    
    # print((x1t*Rbarra/1j).real>=0,(x2t*Rbarra/1j).real>=0)
    
    cte_aux = -Ao*(1j**modo)
    return det_ad*cte_aux

#%%
