#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 20:49:34 2020

@author: leila
electromagnetic energy density of the mode definido en el paper
    de Dionne
"""

import numpy as np
import sys
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/con_kz_real','') 
name_this_py = 'Ver ' + name_this_py

try:
    sys.path.insert(1, path_basic)
    from Etot_conkz import Etot
except ModuleNotFoundError:
    print('Etot_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra Etot_conkz.py')
    sys.path.insert(1, path_basic)
    from Etot_conkz import Etot
    
try:
    sys.path.insert(1, path_basic)
    from Htot_conkz import Htot
except ModuleNotFoundError:
    print('Htot_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra Htot_conkz.py')
    sys.path.insert(1, path_basic)
    from Htot_conkz import Htot

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene)
    from constantes import constantes

pi,hb,c,alfac,mu1,mu2,epsi2 = constantes()

#%%
#print('Definir la funcion densidad de energia')

def energy_density(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z):
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez
    rho : coordenada radial en micrometros
    phi : coordenada azimutal
    z : coordenada longitudinal en micrometros
    
    Returns
    -------
    omega_em = electromagnetic energy density of the mode definido en el paper
    de Dionne

    """
    
    rho = np.abs(rho) #radio positivo
    
    H1_tot,H2_tot = Htot(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z) #|Htot|^2
    E1_tot,E2_tot = Etot(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z) #|Etot|^2

    if rho <=R:
        E_tot,H_tot = E1_tot,H1_tot
    else:
        E_tot,H_tot = E2_tot,H2_tot
    
    rta = epsi1*E_tot + mu1*H_tot
    
    return 0.5*(rta.real)

#%%

