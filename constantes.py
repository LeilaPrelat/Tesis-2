#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

ctes universales ---> IMPORTANTE. Importar en todas las funciones

constantes que se mantienen a lo largo del trabajo (espero jajaja)

Dentro del graphene_sigma.py hay otras constantes que se mantienen
a lo largo de la tesis
"""
import numpy as np

#%%

def constantes():
    pi = np.pi
    hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s
    c = 3*10**(14)                ### light velocity in micron/seg
    alfac = 1/137.0359            ### fine structure
    hbargama = 0.0001             ### collision frequency in eV 
    
    ###medio 2
    mu2, epsi2 = 1,1
    mu1 = 1

    if hbargama!=0.0001:
        raise TypeError('Create new repository')
    if mu1!=1:
        raise TypeError('Create new repository')
    if mu2!=1:
        raise TypeError('Create new repository')
    if epsi2!=1:
        raise TypeError('Create new repository')

    return pi,hb,c,alfac,hbargama,mu1,mu2,epsi2

#%%
