#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

determinante de la matriz de 4x4
    (caso con kz) seccion 3.5 del cuaderno corto 

(reescribir matriz para recuperar el caso sin kz cuando kz = 0)
    
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

pi,hb,c,alfac,mu1,mu2,epsi2 = constantes()

#%%


def determinante(kz_var,omegac,epsi1,mode,R,mu_c):
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
    determinante de la matriz de 4x4
    
    """
    
    ω = omegac*c
    E = ω*hb 
    
    k0 = omegac #=omega/c
    if kz_var == 0 :
        xz = 0
    else:
        xz = kz_var/k0
    Rbarra = R*k0 #adimensional
 
    # def epsi(medio):
    #     if medio ==1:
    #         epsirta = epsi1
    #     elif medio ==2:
    #         epsirta = epsi2
    #     else:
    #         raise TypeError('Wrong value for medium in epsilon')
    #     return epsirta 
        
    # def mu(medio):
    #     if medio ==1:
    #         murta = mu1
    #     elif medio ==2:
    #         murta = mu2
    #     else:
    #         raise TypeError('Wrong value for medium in mu')
    #     return murta 
            
    def xt2(medio):
        if medio == 1:
            xt2rta = -xz**2 + mu1*epsi1            
        elif medio ==2:   
            xt2rta = -xz**2 + mu2*epsi2            
        return xt2rta
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   
        argument = xtmedio*Rbarra

        if argument.real>=0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        if medio ==2:
            if np.abs(argument) <= 1e-50:
            	raise TypeError('El argumento de Hankel es muy chico y Hankel diverge en el 0')
        return xtmedio
    
    def a(j,l):
        if mode!=0:
            rta1 = 1j*xz*mode*xt(j)
            rta2 = Rbarra*xt(l)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def Bessel(mode):
        J = special.jn(mode,xt(1)*Rbarra)+0j
        derJ = special.jvp(mode,xt(1)*Rbarra)+0j
        H =  special.hankel1(mode,xt(2)*Rbarra)+0j
        derH =  special.h1vp(mode,xt(2)*Rbarra)+0j
        return [J,derJ,H,derH]
    
    [J,derJ,H,derH] = Bessel(mode)   
    
    cte_b = xt(1)*xt(2)
    
    cond = sigma(E,mu_c)[0]
    cond3 = cond*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    
    #INTERCAMBIAR LA COLUMNA 2 Y 3 Y LUEGO FILAS 2 CON 4
    
    A = [J*cte_b,-H*cte_b,0,0]

    B = [cond3*1j*J*cte_b - derJ*epsi1*xt(2),derH*epsi2*xt(1),-a(2,1)*J,a(1,2)*H]
    
    C = [cond3*a(2,1)*J,0,J*cte_b - cond3*1j*mu1*xt(2)*derJ,-H*cte_b]
    
    D = [-a(2,1)*J, a(1,2)*H, derJ*mu1*xt(2), -derH*mu2*xt(1)]
    

    M = np.matrix([A,B,C,D], dtype='complex')
    # det = np.linalg.det(M)
    # print(np.abs(det))
    
    cte_det = (mu1*mu2)**2
    
    try:
        rta = (np.linalg.det(M))
    except np.linalg.LinAlgError as Error:
        raise Error 
    return rta/cte_det


# coef(kz,omegac0,crit)

#%%
