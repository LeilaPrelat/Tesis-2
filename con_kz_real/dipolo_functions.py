#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

funciones del dipolo : h_{zm}, h'_{zm}, e_{zm}, e'_{zm}

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

#print('Definir hz del dipolo: hz y hz'/k0')

def hzDIP(kz,omegac,epsi1,modo,hbaramu,rho,p1,p2,pz,rho_D,theta_D,z_D): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    modo: modo entero
    
    hbaramu: potencial quimico del grafeno en eV
    rho : coordenada radial en micrometros
    
    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Returns
        h_{zm}, h'_{zm}/k0
        son diferentes funciones si rho<rho_D o rho>rho_D
    -------
    """
    k0 = omegac #=omega/c
    
    if kz != 0 :
        xz = kz/k0 #adimensional
    else:
        xz = 0
        
    if z_D != 0 :
        x_D = z_D*k0 #adimensional
    else:
        x_D = 0
        
    rhoDbarra = rho_D*k0 #adimensional
    rhobarra = rho*k0        
    
    xt1_2 = -xz**2 + mu1*epsi1
    inside = xt1_2 + 0j
    xt1 = (inside)**(1/2)   
    if (xt1*rhobarra).real>=0:
	    xt1 = xt1
    else:
	    xt1 = -xt1

    def J(nu,x):
        """
        Parameters
        ----------
        nu : modo entero nu
        x : rho*k0 (rhobarra o rhoDbarra)

        Returns
        -------
        J : Bessel function de modo nu 
        y argumento k_{t,1}*rho = x_{t,1}*rho*k0
        """
        J = special.jn(nu,xt1*x) + 0j
        return J
    
    def H(nu,x):
        H =  special.hankel1(nu,xt1*x)+0j
        return H
 
    def derJ(nu,x):
        derJ = special.jvp(nu,xt1*x) + 0j
        return derJ
    
    def derH(nu,x):
        derH =  special.h1vp(nu,xt1*x)+0j
        return derH   
 
    exp1 = np.e**(1j*theta_D)
    exp2 = np.e**(-1j*theta_D)
    
    if rhobarra < rhoDbarra:
        h = -J(modo,rhobarra)*(p1*exp2*H(modo+1,rhoDbarra) + p2*exp1*H(modo-1,rhoDbarra))
        derh = -xt1*derJ(modo,rhobarra)*(p1*exp2*H(modo+1,rhoDbarra) + p2*exp1*H(modo-1,rhoDbarra)) 
        # se devuelve h'/k0 por eso aparece xt1 en lugar de kt1
    else:
        h = -H(modo,rhobarra)*(p1*exp2*J(modo+1,rhoDbarra) + p2*exp1*J(modo-1,rhoDbarra))
        derh = -xt1*derH(modo,rhobarra)*(p1*exp2*J(modo+1,rhoDbarra) + p2*exp1*J(modo-1,rhoDbarra))
         # se devuelve h'/k0 por eso aparece xt1 en lugar de kt1

    exp = np.e**(-1j*(modo*theta_D + xz*x_D))
    cte = 0.25*1j*k0*(xt1*k0)

    return h*exp*cte, derh*exp*cte

#%%

#print('Definir ez del dipolo: ez y ez'/k0')

def ezDIP(kz,omegac,epsi1,modo,hbaramu,rho,p1,p2,pz,rho_D,theta_D,z_D): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    modo: modo entero
    
    hbaramu: potencial quimico del grafeno en eV
    rho : coordenada radial en micrometros
    
    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Returns
        e_{zm}, e'_{zm}/k0
        son diferentes funciones si rho<rho_D o rho>rho_D
    -------
    """
    k0 = omegac #=omega/c
    
    if kz != 0 :
        xz = kz/k0 #adimensional
    else:
        xz = 0
        
    if z_D != 0 :
        x_D = z_D*k0 #adimensional
    else:
        x_D = 0
        
    rhoDbarra = rho_D*k0 #adimensional
    rhobarra = rho*k0        
    
    xt1_2 = -xz**2 + mu1*epsi1
    inside = xt1_2 + 0j
    xt1 = (inside)**(1/2)   
    if (xt1*rhobarra).real>=0:
	    xt1 = xt1
    else:
	    xt1 = -xt1

    def J(nu,x):
        """
        Parameters
        ----------
        nu : modo entero nu
        x : rho*k0 (rhobarra o rhoDbarra)

        Returns
        -------
        J : Bessel function de modo nu 
        y argumento k_{t,1}*rho = x_{t,1}*rho*k0
        """
        J = special.jn(nu,xt1*x) + 0j
        return J
    
    def H(nu,x):
        H =  special.hankel1(nu,xt1*x)+0j
        return H
 
    def derJ(nu,x):
        derJ = special.jvp(nu,xt1*x) + 0j
        return derJ
    
    def derH(nu,x):
        derH =  special.h1vp(nu,xt1*x)+0j
        return derH   
 
    exp1 = np.e**(1j*theta_D)
    exp2 = np.e**(-1j*theta_D)
    aux = modo*(modo+1)*(1/(rhobarra*xt1*rho)-xt1*k0)
    
    if rhobarra < rhoDbarra:
        h = J(modo,rhobarra)*(-p1*kz*exp2*H(modo+1,rhoDbarra) + p2*kz*exp1*H(modo-1,rhoDbarra) + 2*1j*pz*H(modo,rhoDbarra)*aux)
        derh = xt1*derJ(modo,rhobarra)*(-p1*kz*exp2*H(modo+1,rhoDbarra) + p2*kz*exp1*H(modo-1,rhoDbarra) + 2*1j*pz*H(modo,rhoDbarra)*aux)
        # se devuelve h'/k0 por eso aparece xt1 en lugar de kt1
    else:
        h = H(modo,rhobarra)*(-p1*kz*exp2*J(modo+1,rhoDbarra) + p2*kz*exp1*J(modo-1,rhoDbarra) + 2*1j*pz*J(modo,rhoDbarra)*aux)
        derh = xt1*derH(modo,rhobarra)*(-p1*kz*exp2*J(modo+1,rhoDbarra) + p2*kz*exp1*J(modo-1,rhoDbarra) + 2*1j*pz*J(modo,rhoDbarra)*aux)
         # se devuelve h'/k0 por eso aparece xt1 en lugar de kt1

    exp = np.e**(-1j*(modo*theta_D + xz*x_D))
    cte = 0.25*epsi1*(xt1*k0)

    return h*exp*cte, derh*exp*cte

#%%