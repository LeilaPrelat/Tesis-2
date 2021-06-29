#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos longitudinales para el caso con kz:
    
    seccion 3.2 del cuaderno corto

"""

import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real','') 
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from coef_matrix_AoBo import coef
except ModuleNotFoundError:
    print('coef_matrix_AoBo.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra coef_matrix_AoBo.py')
    sys.path.insert(1, path_basic)
    from coef_matrix_AoBo import coef

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

#print('Definir los campos longitudinales')

#Ojo: aca solo van frecuencias reales

def fieldsZ(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
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
        Ez1,Ez2,Hz1,Hz2
    -------
    """

    def coeficientes(modo):
        
        # cte1 = (-1)**modo
        # cte2 = 1j*pi/2
        # cte3 = (-1j)**modo
        # cte41 = (1j)**modo
        # cte4 = cte41*cte2
        
        # coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu)
        coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu,Ao,Bo)
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3
    
        # coef_A1,coef_B2 = coeff[0][0],coeff[1][0]
        # coef_C1,coef_D2 = coeff[2][0],coeff[3][0]
        
        # coef_A1 = coef_A1/cte1
        # coef_C1 = coef_C1/cte1
        
        # coef_B2 = coef_B2/cte2
        # coef_D2 = coef_D2/cte2    
        
        coef_A1 = complex(coef_A1)
        coef_C1 = complex(coef_C1)
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
        # if pol == 's':
        #     coef_D2 = 0
        # elif pol == 'p':
        #     coef_B2 = 0
        
        # print(coef_D2,coef_B2)
        return coef_A1,coef_C1,coef_B2,coef_D2

    k0 = omegac #=omega/c
    xz = kz_var2/k0 #adimensional
    rhobarra = rho*k0 #adimensional
    Rbarra = R*k0

    def epsi(medio):
        if medio ==1:
            return epsi1
        elif medio ==2:
            return epsi2
    
    def mu(medio):
        if medio ==1:
            return mu1
        elif medio ==2:
            return mu2
            
    def xt2(medio):
        return -xz**2 + mu(medio)*epsi(medio)
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   
        if (xtmedio*Rbarra).real>=0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        return xtmedio

    
    Ez1_tot, Ez2_tot = 0,0
    Hz1_tot, Hz2_tot = 0,0
    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    exp2 = np.e**(1j*kz_var2*z) 
    for nu in list_modos: 
        
        # cte1 = (-1)**nu
        # cte2 = 1j*pi/2
        # cte3 = (-1j)**nu
        # cte41 = (1j)**nu
        # cte4 = cte41*cte2
        
        coef_A1,coef_C1,coef_B2,coef_D2 = coeficientes(nu)
        exp1 = np.e**(1j*nu*phi)
        
        J1 = special.jn(nu,xt(1)*rhobarra)+0j
        J2 = special.jn(nu,xt(2)*rhobarra)+0j
        H =  special.hankel1(nu,xt(2)*rhobarra)+0j
    
        # Ez1 = coef_A1*J*cte3*exp1*exp2
        # Ez2 = coef_B2*H*cte4*exp1*exp2
        
        # Hz1 = coef_C1*J*cte3*exp1*exp2
        # Hz2 = coef_D2*H*cte4*exp1*exp2
        
        aux = (1j)**nu
        Ez1 = coef_A1*J1*exp1*exp2
        Ez2 = (coef_B2*H + Bo*aux*J2)*exp1*exp2
        
        Hz1 = coef_C1*J1*exp1*exp2
        Hz2 = (coef_D2*H + Ao*aux*J2)*exp1*exp2
        
        Ez1_tot = Ez1_tot + Ez1
        Ez2_tot = Ez2_tot + Ez2
        
        Hz1_tot = Hz1_tot + Hz1
        Hz2_tot = Hz2_tot + Hz2
    
    return Ez1_tot,Ez2_tot,Hz1_tot,Hz2_tot

#%%
