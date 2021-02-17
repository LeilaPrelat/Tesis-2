#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
constantes pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 que se mantienen 

Dentro del graphene_sigma.py hay otras constantes que se mantienen
a lo largo de la carpeta Tesis-2 y Tesis-1 (aca estaba mal el grafeno)

"""
import numpy as np

#%%

def constantes():
    """
    Returns
    -------
    pi : cte universal pi
    hb : cte universal hbar en eV*s
    c : cte universal vel de la luz en micron/seg
    alfac : cte universal de estructura fina
    hbargama : freq de colision del grafeno en eV
    mu1 : permeabilidad magnetica del medio 1 (interior)
    mu2 : permeabilidad magnetica del medio 2 (exterior)
    epsi2 : permeabilidad electrica del medio 2 (exterior)
    """

    global mu1,mu2,epsi2,hbargama ### variables globales ---> no pueden cambiarse
    pi = np.pi
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    c = 3*10**(14)                ### light velocity in micron/seg
    alfac = 1/137.0359            ### fine structure
    hbargama = 0.0001             ### collision frequency in eV 
    
    ###medio 2

    mu2, epsi2 = 1,1
    mu1 = 1

    # if hbargama!=0.0001:
    #     raise TypeError('Create new repository')
    # if mu1!=1:
    #     raise TypeError('Create new repository')
    # if mu2!=1:
    #     raise TypeError('Create new repository')
    # if epsi2!=1:
    #     raise TypeError('Create new repository')

    return pi,hb,c,alfac,hbargama,mu1,mu2,epsi2

#%%
