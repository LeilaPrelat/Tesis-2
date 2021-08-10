#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Usamos las nuevas formulas (simplificacion de Qabs ---> para Qabs2, ver cuaderno) 
para calcular las secciones eficaces, hecho SOLAMENTE para la polarizacion p
para el medio NO dispersivo (solamente medio activo en el interior del cilindro)

"""
import numpy as np
import os 
import sys
from scipy import special as sp #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/non-dispersive','') 

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene)
    path_graphene2 = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene2)
    from graphene_sigma import sigma

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

#print('Definir el coeficiente Qscat')

#Ojo: aca solo van frecuencias reales

def Qscat(omegac,epsi1,nmax,R,hbaramu,Ao): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    Ao: amplitud del campo Hz incidente (pol p)
    
    Returns
        Qscat ADIMENSIONAL (dividido por c) caso NO-DISPERSIVO
        despreciando campos incidentes  (pol p UNICAMENTE)
    -------
    """
    nmax = int(nmax)
    # epsi1 = re_epsi1 + 1j*im_epsi1
    # omega_c = 2*pi/lambbda
    omega = omegac*c
    energy = omega*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac #=omega/c
    Rbarra = R*k0 #adimensional
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
    x1t,x2t = x1, x2
    
    if (x1t*Rbarra).real>=0:
 	    x1t = x1t
    else:
 	    x1t = -x1t
    
    if (x2t*Rbarra).real>=0:
 	    x2t = x2t
    else:
 	    x2t = -x2t
        
    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    Sigma_ad = 4*pi*1j*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
           
    def coef_an_bn(modo):
        
        def bessel1(modo):
            CJ1 = sp.jn(modo,x1t*Rbarra)+0j
            DCJ1 = sp.jvp(modo,x1t*Rbarra)+0j 
            return CJ1,DCJ1
    
        def bessel2(modo):
            CJ2 = sp.jn(modo,x2t*Rbarra)+0j
            DCJ2 = sp.jvp(modo,x2t*Rbarra)+0j 
            CH2 = sp.hankel1(modo,x2t*Rbarra)+0j
            DCH2 = sp.h1vp(modo,x2t*Rbarra)+0j  
            return CJ2,DCJ2,CH2,DCH2
        
        J1,derJ1 = bessel1(modo)
        J2,derJ2,H2,derH2 = bessel2(modo)
        
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        bn = 0
        
        # print(an,bn)
        return an,bn
        
    Qscat = 0
    cte = 4*pi*epsi2*Rbarra #devuelve Qscat/c ---> adimensional
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for nu in list_modos: 
        nu = int(nu)
        an,bn = coef_an_bn(nu)

        term2 = np.abs(an)**2  #pol p 

        Qscat = Qscat + term2 
        
    Qscatf = Qscat/cte
    
    return Qscatf

#%%

def Qabs(omegac,epsi1,nmax,R,hbaramu,Ao): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    Ao: amplitud del campo Hz incidente (pol p)
    
    Returns
        Qabs ADIMENSIONAL (dividido por c) caso NO-DISPERSIVO
        despreciando campos incidentes  (pol p UNICAMENTE)
    -------
    """
    nmax = int(nmax)
    # epsi1 = re_epsi1 + 1j*im_epsi1
    # omega_c = 2*pi/lambbda
    omega = omegac*c
    energy = omega*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac #=omega/c
    Rbarra = R*k0 #adimensional
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
    x1t,x2t = x1, x2
    
    if (x1t*Rbarra).real>=0:
 	    x1t = x1t
    else:
 	    x1t = -x1t
         

    if (x2t*Rbarra).real>=0:
 	    x2t = x2t
    else:
 	    x2t = -x2t

    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    Sigma_ad =  4*pi*1j*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
           
    def coef_an_bn(modo):
        
        def bessel1(modo):
            CJ1 = sp.jn(modo,x1t*Rbarra)+0j
            DCJ1 = sp.jvp(modo,x1t*Rbarra)+0j 
            return CJ1,DCJ1
    
        def bessel2(modo):
            CJ2 = sp.jn(modo,x2t*Rbarra)+0j
            DCJ2 = sp.jvp(modo,x2t*Rbarra)+0j 
            CH2 = sp.hankel1(modo,x2t*Rbarra)+0j
            DCH2 = sp.h1vp(modo,x2t*Rbarra)+0j  
            return CJ2,DCJ2,CH2,DCH2
        
        J1,derJ1 = bessel1(modo)
        J2,derJ2,H2,derH2 = bessel2(modo)
        
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        bn = 0  
        
        # print(an,bn)
        return an,bn

    Qemi_SR2 = 0    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for nu in list_modos: 
        nu = int(nu)
        aux = (-1j)**nu
        
        an,bn = coef_an_bn(nu)
        
        an2= np.abs(an)**2  
        AoC = Ao.conjugate() 
   
        term3 = ((AoC*aux*an).real) + an2 #S(R+) #el conjugado de la parte real qeu aparece en la eq 2.81 de la tesis
        Qemi_SR2 = Qemi_SR2 + term3  # S(R+) se calculo con normal exterior +rho ---> en realidad es Qemision
    
    cte_f = 4*pi*epsi2*Rbarra
    
    Qabs_SR2f = -Qemi_SR2/cte_f
    
    return Qabs_SR2f


#%%

# Qext = Qabs + Qscat, tengo miedo de que suceda cancelacion catastrofica, asi que lo voy a escribir a mano


def Qext(omegac,epsi1,nmax,R,hbaramu,Ao): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    Ao: amplitud del campo Hz incidente (pol p)
    
    Returns
        Qext ADIMENSIONAL (dividido por c) caso NO-DISPERSIVO
        despreciando campos incidentes  (pol p UNICAMENTE)
    -------
    """
    nmax = int(nmax)
    # epsi1 = re_epsi1 + 1j*im_epsi1
    # omega_c = 2*pi/lambbda
    omega = omegac*c
    energy = omega*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac #=omega/c
    Rbarra = R*k0 #adimensional
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
    x1t,x2t = x1, x2
    
    if (x1t*Rbarra).real>=0:
 	    x1t = x1t
    else:
 	    x1t = -x1t
    
    if (x2t*Rbarra).real>=0:
 	    x2t = x2t
    else:
 	    x2t = -x2t
        
    # cte_misterio = alfac*c
    # Sigma_ad = 4*pi*sigmatot*cte_misterio/c
    Sigma_ad = 4*pi*1j*alfac*sigmatot #sigma devuelve la conductividad teniendo que multiplicar por alfac*c ---> no hay que dividir por c
           
    def coef_an_bn(modo):
        
        def bessel1(modo):
            CJ1 = sp.jn(modo,x1t*Rbarra)+0j
            DCJ1 = sp.jvp(modo,x1t*Rbarra)+0j 
            return CJ1,DCJ1
    
        def bessel2(modo):
            CJ2 = sp.jn(modo,x2t*Rbarra)+0j
            DCJ2 = sp.jvp(modo,x2t*Rbarra)+0j 
            CH2 = sp.hankel1(modo,x2t*Rbarra)+0j
            DCH2 = sp.h1vp(modo,x2t*Rbarra)+0j  
            return CJ2,DCJ2,CH2,DCH2
        
        J1,derJ1 = bessel1(modo)
        J2,derJ2,H2,derH2 = bessel2(modo)
        
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        bn = 0
        
        # print(an,bn)
        return an,bn
        

    Qemi_SR2 = 0    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for nu in list_modos: 
        nu = int(nu)
        aux = (-1j)**nu
        
        an,bn = coef_an_bn(nu)
    
        AoC = Ao.conjugate() 
   
        term3 = ((AoC*aux*an).real) 
        Qemi_SR2 = Qemi_SR2 + term3  # S(R+) se calculo con normal exterior +rho ---> en realidad es Qemision
    
    cte_f = 4*pi*epsi2*Rbarra
    
    Qabs_SR2f = -Qemi_SR2/cte_f
    
    return Qabs_SR2f


