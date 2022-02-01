#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

determinante de la matriz de 4x4
    (caso con kz)
tal como sale de plantear las condiciones de borde
"""

import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real/real_freq_lorentziana','') 

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene)
    path_graphene2 = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene2)
    from graphene_sigma import sigma
    
try:
    sys.path.insert(1, path_basic)
    from epsilon_lorentziana import epsilon 
except ModuleNotFoundError:
    print('epsilon_lorentziana.py no se encuentra en ' + path_basic)
    path_basic2 = input('path de la carpeta donde se encuentra epsilon_lorentziana.py')
    sys.path.insert(1, path_basic)
    from epsilon_lorentziana import epsilon

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()
aux_cte = c*1e-12 # freq en THz para epsilon

#%%

def determinante(kz_var,omegac,eta,mode,R,hbaramu): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    eta : inversion de poblacion
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV

    Returns
    -------
    determinante de la matriz de 4x4
    
    """
    E = omegac*(c*hb)
    omegga = omegac*aux_cte
    freqq = omegga/(2*np.pi)
    k0 = omegac #=omega/c
    xz = kz_var/k0
    
    Rbarra = R*k0 #adimensional
 
    def epsi(medio):
        if medio == 1:
            epsirta = complex(epsilon(freqq, eta))  # habia problemas con listas en solve_det_conkz.pys
        elif medio == 2:
            epsirta = epsi2
        else:
            raise TypeError('Wrong value for medium in epsilon')
        return epsirta 
        
    def mu(medio):
        if medio == 1:
            murta = mu1
        elif medio == 2:
            murta = mu2
        else:
            raise TypeError('Wrong value for medium in mu')
        return murta 
            
    def xt2(medio):
        xt2rta = -xz**2 + mu(medio)*epsi(medio)      
        return xt2rta
            
    def xt(medio):
        inside = xt2(medio)+0j
        xtmedio = (inside)**(1/2)   
        argument = xtmedio*Rbarra

        if argument.real>=0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
        if medio == 2: 
            if np.abs(argument) <= 1e-50:
            	raise TypeError('El argumento de Hankel es muy chico y Hankel diverge en el 0')
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
    
    def Bessel(modo):
        J = special.jn(modo,xt(1)*Rbarra)+0j
        derJ = special.jvp(modo,xt(1)*Rbarra)+0j
        H =  special.hankel1(modo,xt(2)*Rbarra)+0j
        derH =  special.h1vp(modo,xt(2)*Rbarra)+0j
        return [J,derJ,H,derH]
    
    [J,derJ,H,derH] = Bessel(mode)   
    
    cond = sigma(E,hbaramu,hbargama)[0]
    cond3 = cond*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
    
    #INTERCAMBIAR LA COLUMNA 2 Y 3 Y LUEGO FILAS 2 CON 4
    
    A = [J,-H,0,0]

    B = [cond3*J + d(1)*derJ,-d(2)*derH,-a(1)*J,a(2)*H]
    
    C = [cond3*a(1)*J,0,J+cond3*b(1)*derJ,-H]
    
    D = [-a(1)*J, a(2)*H,-b(1)*derJ,b(2)*derH]
    

    M = np.matrix([A,B,C,D], dtype='complex')
    # det = np.linalg.det(M)
    # print(np.abs(det))
    
    try:
        rta = (np.linalg.det(M))
    except np.linalg.LinAlgError as Error:
        raise Error 
    return rta


# coef(kz,omegac0,crit)

#%%
