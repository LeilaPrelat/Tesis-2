#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

coeficientes de la matriz de 4x4
      (caso con kz)
caso inhomogeneo (con campos incidentes) + dipolo dentro del cilindro 
"""

import numpy as np
import sys
import os 
from scipy import special #funciones de Bessel
#from numpy import linalg as LA

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real','') 

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
    from dipolo_functions import hzDIP, ezDIP
except ModuleNotFoundError:
    print('dipolo_functions.py no se encuentra en ' + path_basic)
    path_graphene2 = input('path de la carpeta donde se encuentra dipolo_functions.py')
    sys.path.insert(1, path_basic)
    from dipolo_functions import hzDIP, ezDIP

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


def coef(kz,omegac,epsi1,modo,R,mu_c,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
    """

    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    modo: modo entero
    
    mu_c: potencial quimico del grafeno en eV
    
    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez

    Returns
    -------
    coef de la matriz de 4x4, 2 medios, con grafeno
        (A1,B2,C1,D2) del min{|autovalor1|,|autovalor2|,|autovalor3|,|autovalor4|}
    
    """
    
    E = omegac*(c*hb) 
    k0 = omegac #=omega/c
    xz = kz/k0
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
        if modo!=0:
            rta1 = xz*modo
            rta2 = Rbarra*xt2(j)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(j):
        return 1j*mu(j)/xt(j)

    def d(j):
        return 1j*epsi(j)/xt(j)
    
    ez,der_ez = ezDIP(kz,omegac,epsi1,modo,mu_c,R,p1,p2,pz,rho_D,theta_D,z_D)
    hz,der_hz = hzDIP(kz,omegac,epsi1,modo,mu_c,R,p1,p2,pz,rho_D,theta_D,z_D)
    
    J = special.jn(modo,xt(1)*Rbarra) + 0j
    J2 = special.jn(modo,xt(2)*Rbarra) + 0j
    derJ = special.jvp(modo,xt(1)*Rbarra) + 0j
    derJ2 = special.jvp(modo,xt(2)*Rbarra) + 0j
    H =  special.hankel1(modo,xt(2)*Rbarra) + 0j
    derH =  special.h1vp(modo,xt(2)*Rbarra) + 0j

    cond = sigma(E,mu_c,hbargama)[0]
    cond3 = cond*alfac*4*pi  #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
  
    cte = (1j)**modo
    cte_z = xz*modo/Rbarra
    #INTERCAMBIAR LA COLUMNA 2 Y 3 Y LUEGO FILAS 2 CON 4
    
    A = [J,-H,0,0]

    B = [cond3*J + d(1)*derJ,-d(2)*derH,-a(1)*J,a(2)*H]
    
    C = [cond3*a(1)*J,0,J+cond3*b(1)*derJ,-H]
    
    D = [-a(1)*J, a(2)*H,-b(1)*derJ,b(2)*derH]
    
    col1 = Bo*J2 - ez
    col2 = cte_z*(hz/xt2(1) - cte*Ao*J2/xt2(2)) + 1j*(epsi2*cte*Bo*derJ2/xt(2) - epsi1*der_ez/xt2(1)) - cond3*ez
    col3 = Ao*cte*J2 - hz - (cond3/xt2(1))*(ez*cte_z + 1j*mu1*der_hz)
    col4 = 1j*(mu1*der_hz -mu2*cte*Ao*derJ2/xt(2)) + cte_z*(ez/xt2(1) - cte*Bo*J2/xt2(2))
    
    b = np.array([col1, col2, col3, col4])
    M = np.matrix([A,B,C,D], dtype='complex')
    
    try:
        x = np.linalg.solve(M, b)      
    except np.linalg.LinAlgError as Error:
        raise Error 
    return x


# coef(kz,omegac0,crit)

#%%
