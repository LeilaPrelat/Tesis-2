#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 10:30:50 2020

@author: leila

Usamos las nuevas formulas (simplificacion de Qabs ---> para Qabs2, ver cuaderno) 
para calcular las secciones eficaces, hecho SOLAMENTE para la polarizacion p
para el medio NO dispersivo (solamente medio activo en el interior del cilindro)

"""

import numpy as np
import sys
import os 
from scipy import special as sp #funciones de Bessel

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz_inv/non-dispersive','') 

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

#print('Definir parametros del problema')

#Ao = 1 #el resultado depende de este parametro. parece que ya no: chequear
#Bo = 0
#info_pol = 'p'

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
           
    def coeficientes(modo):
    
        J1,derJ1 = bessel1(modo)
        J2,derJ2,H2,derH2 = bessel2(modo)
        
        aux =   (1j)**(modo)
        an_num = -Ao*aux*(epsi1*x2*J1*derJ2 - epsi2*x1*J2*derJ1 + Sigma_ad*x1*x2*derJ2*derJ1)
        an_den = epsi1*x2*J1*derH2 - epsi2*x1*H2*derJ1 + Sigma_ad*x1*x2*derH2*derJ1
        an = an_num/an_den
        
        bn = 0 #polarizacion p 
        
        cn_num = Ao*aux*x2*epsi1*(derJ2*H2-derH2*J2) 
        cn_den = -an_den
        cn = cn_num/cn_den
        
        dn = 0 #polarizacion p 
        
        return an,bn,cn,dn
    
    
    cte11 = mu1/epsi1

    Qemi_SR1 = 0
    Qemi_SR2 = 0    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    for nu in list_modos: 
        nu = int(nu)
        aux = (-1j)**nu
        
        J1,derJ1 = bessel1(nu)        
        J2,derJ2,H2,derH2 = bessel2(nu)
        an,bn,cn,dn = coeficientes(nu)
        
        an2= np.abs(an)**2  
        cn2 = np.abs(cn)**2     
        J1C = J1.conjugate()
   
        term1 = 1j*cn2*J1C*derJ1*((cte11)**(1/2)) #viene de <S(R-)>
        SR1 = (-term1).real # = S(R-)
   
        term3 = Ao*((aux*an).real) + an2 #S(R+)
              
        # Qabs = Qabs + SR2 - SR1 #es con normal exterior rho ---> cambiar el signo 
                                        # SR2 - SR1 ---> -SR2 + SR1
        # Qabs = Qabs - SR2 + SR1

        Qemi_SR1 = Qemi_SR1 + SR1  # S(R-) se calculo con normal exterior +rho ---> dar vuelta el signo
        Qemi_SR2 = Qemi_SR2 + term3  # S(R+) se calculo con normal exterior +rho ---> en realidad es Qemision
    
    cte_f = 4*pi*epsi2*Rbarra
    
    Qabs_SR2f, Qabs_SR1f = -Qemi_SR2/cte_f, -Qemi_SR1*0.125
    
    #Qemision tot = Qemi_SR2 - Qemi_SR1
    
    #Qabsorcion tot = -Qemision tot = Qemi_SR1 - Qemi_SR2 = Qabs_SR2 - Qabs_SR1
    
    return Qabs_SR2f

#%%
