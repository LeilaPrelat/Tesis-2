#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

coeficientes de la matriz de 4x4
      (caso con kz)
caso homogeneo (sin campos incidentes)  
"""

import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel
from numpy import linalg as LA

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real','') 
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_basic)
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


def coef(kz_var,omegac,epsi1,mode,R,mu_c):
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
    coef de la matriz de 4x4, 2 medios, con grafeno
        (A1,B2,C1,D2) del min{|autovalor1|,|autovalor2|,|autovalor3|,|autovalor4|}
    
    """
    
    ω = omegac*c
    E = ω*hb 
    
    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
 
    def epsi(medio):
        if medio ==1:
            epsirta = epsi1
        elif medio ==2:
            epsirta = epsi2
        else:
            raise TypeError('Wrong value for medium in epsilon')
        return epsirta 
        
    def mu(medio):
        if medio ==1:
            murta = mu1
        elif medio ==2:
            murta = mu2
        else:
            raise TypeError('Wrong value for medium in mu')
        return murta 
            
    def xt2(medio):
        xt2rta = -xz**2 + mu(medio)*epsi(medio)
        if medio == 2 and xt2rta <= 1e-23:
            raise TypeError('Se anulo el argumento de Hankel y Hankel diverge en el 0')
        return xt2rta
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   

        if (xtmedio*Rbarra).real>=0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        return xtmedio
    
    def a(j):
        if mode!=0:
            rta1 = xz*mode
            rta2 = Rbarra*xt2(j)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(j):
        return 1j*mu(j)/xt(j)

    def d(j):
        return 1j*epsi(j)/xt(j)
    
    def Bessel(mode):
        J = special.jn(mode,xt(1)*Rbarra)+0j
        derJ = special.jvp(mode,xt(1)*Rbarra)+0j
        H =  special.hankel1(mode,xt(2)*Rbarra)+0j
        derH =  special.h1vp(mode,xt(2)*Rbarra)+0j
        return [J,derJ,H,derH]
    
    [J,derJ,H,derH] = Bessel(mode)   
    
    cond = sigma(E,mu_c)[0]
    cond3 = cond*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
     
    #INTERCAMBIAR LA COLUMNA 2 Y 3 Y LUEGO FILAS 2 CON 4
    
    A = [J,-H,0,0]

    B = [cond3*J + d(1)*derJ,-d(2)*derH,-a(1)*J,a(2)*H]
    
    C = [cond3*a(1)*J,0,J+cond3*b(1)*derJ,-H]
    
    D = [-a(1)*J, a(2)*H,-b(1)*derJ,b(2)*derH]
    

    M = np.matrix([A,B,C,D], dtype='complex')
    
    # Madj = M.getH()
    # Mherm = np.dot(M,Madj)
    # try:
    #     λ, v = LA.eig(Mherm)
    #     index_min = int(np.argmin(λ))        
    # except np.linalg.LinAlgError as Error:
    #     raise Error 
    
    # det = np.linalg.det(M)
    # print(np.abs(det))
    
    try:
        λ, v = LA.eig(M)
        mod_λ = []
        for autovalor in λ:
            mod_λ.append(np.abs(autovalor))
        minλ = np.min(mod_λ)
        index_min = int(np.argmin(mod_λ))        
    except np.linalg.LinAlgError as Error:
        raise Error 
    # print(minλ)
    # print(v[:,index_min])
    return v[:,index_min]


# coef(kz,omegac0,crit)

#%%
