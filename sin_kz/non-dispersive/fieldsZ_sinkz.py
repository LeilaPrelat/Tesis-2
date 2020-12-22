#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

obtuve gn (minimizandolo) ---> obtuve
los denominadores de los coef an y cn
(entonces obtuve an y cn) 
---> podemos graficar los campos

"""

import scipy.special as sp #bessel functions
import numpy as np
import os 
import sys

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/non-dispersive','') 

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma

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
    
#print('Definir la funcion Hz para el medio 1 y el 2')

def Hz(omegac,epsi1,nmax,R,hbaramu,Ao,rho,phi): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao: pol p ---> Hz (Amplitud del campo incidente en Hz del medio 2)
    rho: coordenada radial en micrometros
    phi: coordenada azimutal

    Returns
        Hz del medio 1 si rho<R y del medio 2 si rho>R
    -------
    """
    energy = omegac*c*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac
    Rbarra = R*k0 #adimensional
    rhobarra = rho*k0
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
    
    def coef_an_cn(modo):
        
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
        
        Sigma_ad = 4*pi*1j*sigmatot*alfac #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        cn_num = Ao*aux*x2*epsi1*(derJ2*H2-derH2*J2) 
        cn_den = -an_den
        cn = cn_num/cn_den
        
        # print(an,bn)
        return an,cn
    

    Hz1tot, Hz2tot = 0,0 
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for modo in list_modos:
        j1,j2,h2 = sp.jn(modo,x1*rhobarra), sp.jn(modo,x2*rhobarra), sp.hankel1(modo,x2*rhobarra)
        an,cn = coef_an_cn(modo)
        aux = (1j)**(modo)
        
        Hz1 = cn*j1*np.e**(1j*modo*phi)
        Hz2 = (Ao*aux*j2 + an*h2)*np.e**(1j*modo*phi)
        
        Hz1tot = Hz1 + Hz1tot
        Hz2tot = Hz2 + Hz2tot
    
    return [Hz1tot,Hz2tot]

#%%        

#print('Definir la funcion Ez para el medio 1 y el 2')

def Ez(omegac,epsi1,nmax,R,hbaramu,Bo,rho,phi): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Bo: pol s ---> Ez (Amplitud del campo incidente en Ez del medio 2)
    rho: coordenada radial en micrometros
    phi: coordenada azimutal

    Returns
        Ez del medio 1 si rho<R y del medio 2 si rho>R
    -------
    """
    energy = omegac*c*hb

    sigmatot, inter, intra = sigma(energy,hbaramu) 
        
    k0 = omegac
    Rbarra = R*k0 #adimensional
    rhobarra = rho*k0
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
 
    def coef_bn_dn(modo):
        
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
        
        Sigma_ad = 4*pi*1j*sigmatot*alfac #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
        aux =   (1j)**(modo)
        bn_num = Bo*aux*(epsi1*x2*J2*derJ1 - epsi2*x1*J1*derJ2 - Sigma_ad*x1*x2*J2*J1)
        bn_den = epsi2*x1*J1*derH2 - epsi1*x2*H2*derJ1 + Sigma_ad*x1*x2*H2*J1
        bn = bn_num/bn_den
        
        dn_num = Bo*aux*x1*epsi2*(derJ2*H2-derH2*J2) 
        dn_den = -bn_den
        dn = dn_num/dn_den
        
        # print(an,bn)
        return bn,dn   


    Ez1tot, Ez2tot = 0,0 
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for modo in list_modos:
        
        j1,j2,h2 = sp.jn(modo,x1*rhobarra), sp.jn(modo,x2*rhobarra), sp.hankel1(modo,x2*rhobarra)
        bn,dn = coef_bn_dn(modo)
        aux = (1j)**(modo)   
        
        Ez1 = dn*j1*np.e**(1j*modo*phi)
        Ez2 = (Bo*aux*j2 + bn*h2)*np.e**(1j*modo*phi)

        Ez1tot = Ez1 + Ez1tot
        Ez2tot = Ez2 + Ez2tot 
    
    return [Ez1tot,Ez2tot]

#%%

