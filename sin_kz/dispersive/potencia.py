#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 10:30:50 2020

@author: leila

Antes de integrar en theta para obtener Qscat (despreciando campos incidentes) hecho SOLAMENTE para la polarizacion p
para el medio dispersivo (medio activo + nanocristal en el interior del cilindro)

la idea es ver la potencia es los casos degenerados 

"""

import numpy as np
import sys
import os 
from scipy import special as sp #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/dispersive','')

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

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir el coeficiente Qscat y Qabs')

#Ojo: aca solo van frecuencias reales

def potencia(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,nmax,hbaramu,Ao,rho,phi): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    epsi_ci: Im(epsilon1) del colorante (dye)
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    hbaramu: potencial quimico del grafeno en eV
    Ao: amplitud del Hz incidente (solo pol p)
    rho: coordenada radial en micrometros
    phi : variable angular de las coordenadas polares (de 0 a 2pi)
    
    Returns
        potencia disipada caso DISPERSIVO
        despreciando campos incidentes  (pol p UNICAMENTE)
    -------
    """
    nmax = int(nmax)
    # omega = omegac*c
    energy = omegac*c*hb

    num = Ep**2
    den = energy**2 + 1j*gamma_DL*energy
    epsi1 = epsiinf_DL - num/den + 1j*epsi_ci

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac #=omega/c
    Rbarra = rho*k0 #adimensional
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
    
    if (x1*Rbarra).real>=0:
 	    x1 = x1
    else:
 	    x1 = -x1
    
    if (x2*Rbarra).real>=0:
 	    x2 = x2
    else:
 	    x2 = -x2
         
    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    Sigma_ad = 4*pi*1j*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
  
           
    def coef_an(modo):
        
        def bessel1(modo):
            CJ1 = sp.jn(modo,x1*Rbarra)+0j
            DCJ1 = sp.jvp(modo,x1*Rbarra)+0j 
            return CJ1,DCJ1
    
        def bessel2(modo):
            CJ2 = sp.jn(modo,x2*Rbarra)+0j
            DCJ2 = sp.jvp(modo,x2*Rbarra)+0j 
            CH2 = sp.hankel1(modo,x2*Rbarra)+0j
            DCH2 = sp.h1vp(modo,x2*Rbarra)+0j  
            return CJ2,DCJ2,CH2,DCH2
        
        J1,derJ1 = bessel1(modo)
        J2,derJ2,H2,derH2 = bessel2(modo)
        
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        arg = x2*Rbarra - modo*0.5*pi - 0.25*pi
        # print(an,bn)
        return an,arg
        
    potencia = 0
    cte = (2/(pi*x2*Rbarra))**(1/2) #devuelve Qscat/c ---> adimensional
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    Ao2 = np.conjugate(Ao)
    
    for nu2 in list_modos:
        for nu in list_modos: 
            cte_modo = (1j)**nu
            cte_modo2 = (-1j)**nu2
            
            an,arg = coef_an(nu)
            an2,arg2 = coef_an(nu2)
            
            term1 = cte*(Ao*cte_modo*np.sin(arg) - 1j*an*(np.e**arg))
            term2 = cte*(Ao2*cte_modo2*np.cos(arg2) + an2*(np.e**(-1j*arg2)))
    
            potencia = potencia + 1j*(mu2/x2)*term1*np.e**(1j*nu*phi) + term2
    
    return potencia.real

#%%
