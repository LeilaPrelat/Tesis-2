#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

k transversal cuadrado = k**2 - kz**2 ---> determina el tipo de onda
"""

import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real/extra','') 


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


def kt(kz_var,omegac,epsi1,mode,R,mu_c):
    """
    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio1
    mode : modo
    R : radio del cilindro en micrometros
    mu_c : potencial quimico del grafeno

    Returns
    -------
    k transversal del medio 1---> el signo determina si es onda guiada u onda localizada
	(el medio 1 porque es el interior del cilindro en donde se comienza a propagar la onda)
    
    """
    
    k0 = omegac #=omega/c
    Rbarra = R*k0
    xz = kz_var/k0
    xt2 = -xz**2 + mu1*epsi1          

    inside = xt2 + 0j
    xtmedio = (inside)**(1/2)   
    argument = xtmedio*Rbarra

    if argument.real >=0:
	    xtmedio = xtmedio
    else:
	    xtmedio = -xtmedio
    
    if argument <= 1e-50:
        raise TypeError('El argumento de Hankel es muy chico y Hankel diverge en el 0')

    return xtmedio*k0 #k transversal

#%%
