#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos transversales + longitudinales caso con kz:
    
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

#print('Definir los campos electricos transversales: Erho y Ephi')

#Ojo: aca solo van frecuencias reales

def Et(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        Erho1,Erho2,Ephi1,Ephi2
    -------
    """
    
    def coeficientes(modo):
        
        coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu,Ao,Bo)
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3   
        
        coef_A1 = complex(coef_A1)
        coef_C1 = complex(coef_C1)
        coef_B2 = complex(coef_B2)
        coef_D2 = complex(coef_D2)
        
        return coef_A1,coef_C1,coef_B2,coef_D2

    k0 = omegac #=omega/c
    if kz_var2 != 0 :
        xz = kz_var2/k0 #adimensional
    else:
        xz = 0
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
    
    def Bessel(modo):
        J = special.jn(modo,xt(1)*rhobarra)+0j
        derJ = special.jvp(modo,xt(1)*rhobarra)+0j
        H =  special.hankel1(modo,xt(2)*rhobarra)+0j
        derH =  special.h1vp(modo,xt(2)*rhobarra)+0j
        J2 = special.jn(modo,xt(2)*rhobarra)+0j
        derJ2 = special.jvp(modo,xt(2)*rhobarra)+0j
        return [J,J2,derJ,derJ2,H,derH]  
    
    
    def a(j,modo):  #no esta incluido el xz
        if modo!=0:
            rta1 = modo
            rta2 = rhobarra*xt2(j)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(j):
        return 1j*mu(j)/xt(j)

    def d(j):
        return 1j*epsi(j)/xt(j)
    
    def f(j):
        if kz_var2 != 0:
            return 1j*xz/xt(j)
        else:
            return 0
    
    Erho1_tot, Erho2_tot = 0,0
    Ephi1_tot, Ephi2_tot = 0,0
    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    exp2 = np.e**(1j*kz_var2*z)   
    for nu in list_modos: 
        coef_A1,coef_C1,coef_B2,coef_D2 = coeficientes(nu)
        
        # cte1 = (-1)**nu
        # cte2 = 1j*pi/2
        # cte3 = (-1j)**nu
        # cte41 = (1j)**nu
        # cte4 = cte41*cte2
        aux = (1j)**nu        
        exp1 = np.e**(1j*nu*phi)
         
        [J,J2,derJ,derJ2,H,derH]  = Bessel(nu)
        
        # Erho1 = (-coef_C1*J*a(1,nu)*mu(1) + coef_A1*derJ*f(1))*exp1*exp2*cte3
        # Erho2 =  (-coef_D2*H*a(2,nu)*mu(2) + coef_B2*derH*f(2))*exp1*exp2*cte4
        
        # Ephi1 = (-coef_A1*J*a(1,nu)*xz - coef_C1*derJ*b(1))*exp1*exp2*cte3
        # Ephi2 = (-coef_B2*H*a(2,nu)*xz - coef_D2*derH*b(2))*exp1*exp2*cte4
        
        Erho1 = (-coef_C1*J*a(1,nu)*mu(1) + coef_A1*derJ*f(1))*exp1*exp2        
        Ephi1 = (-coef_A1*J*a(1,nu)*xz - coef_C1*derJ*b(1))*exp1*exp2

        Erho2 =  (-a(2,nu)*mu(2)*(coef_D2*H + Ao*aux*J2) + f(2)*(coef_B2*derH + Bo*aux*derJ2))*exp1*exp2
        Ephi2 = (-a(2,nu)*xz*(coef_B2*H + Bo*aux*J2) - b(2)*(coef_D2*derH + Ao*aux*derJ2))*exp1*exp2
        
        Erho1_tot = Erho1_tot + Erho1
        Erho2_tot = Erho2_tot + Erho2
        
        Ephi1_tot = Ephi1_tot + Ephi1
        Ephi2_tot = Ephi2_tot + Ephi2

    return Erho1_tot,Erho2_tot,Ephi1_tot,Ephi2_tot

#print('Definir los campos longitudinales')

#Ojo: aca solo van frecuencias reales

def Ez(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        Ez1,Ez2
    -------
    """
    
    def coeficientes(modo):
        
        coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu,Ao,Bo)
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3   
        
        coef_A1 = complex(coef_A1)
        coef_C1 = complex(coef_C1)
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
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
    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    exp2 = np.e**(1j*kz_var2*z) 
    for nu in list_modos: 
        coef_A1,coef_C1,coef_B2,coef_D2 = coeficientes(nu)
        exp1 = np.e**(1j*nu*phi)
        
        # cte1 = (-1)**nu
        # cte2 = 1j*pi/2
        # cte3 = (-1j)**nu
        # cte41 = (1j)**nu
        # cte4 = cte41*cte2
        
        J1 = special.jn(nu,xt(1)*rhobarra)+0j
        J2 = special.jn(nu,xt(2)*rhobarra)+0j
        H =  special.hankel1(nu,xt(2)*rhobarra)+0j
    
        # Ez1 = coef_A1*J*cte3*exp1*exp2
        # Ez2 = coef_B2*H*cte4*exp1*exp2
        
        aux = (1j)**nu
        Ez1 = coef_A1*J1*exp1*exp2
        Ez2 = (coef_B2*H + Bo*aux*J2)*exp1*exp2
        
        Ez1_tot = Ez1_tot + Ez1
        Ez2_tot = Ez2_tot + Ez2
    
    return Ez1_tot,Ez2_tot

def Etot(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        |Etot|^2 =  |Erho|^2 + |Ephi|^2 + |Ez|^2
    -------
    """
    
    Ez1, Ez2 = Ez(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    Erho1,Erho2,Ephi1,Ephi2 = Et(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    
    Erho1_tot, Erho2_tot = np.abs(Erho1)**2, np.abs(Erho2)**2
    Ephi1_tot, Ephi2_tot = np.abs(Ephi1)**2, np.abs(Ephi2)**2
    Ez1_tot, Ez2_tot = np.abs(Ez1)**2, np.abs(Ez2)**2
    
    E1_tot = Ez1_tot + Erho1_tot + Ephi1_tot
    E2_tot = Ez2_tot + Erho2_tot + Ephi2_tot
    return E1_tot,E2_tot 

#%%

#print('Definir los campos magneticos transversales: Hrho y Hphi')

#Ojo: aca solo van frecuencias reales

def Ht(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        Hrho1,Hrho2,Hphi1,Hphi2
    -------
    """
    
    def coeficientes(modo):
        
        coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu,Ao,Bo)
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3   
        
        coef_A1 = complex(coef_A1)
        coef_C1 = complex(coef_C1)
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
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
    
    def Bessel(mode):
        J = special.jn(mode,xt(1)*rhobarra)+0j
        derJ = special.jvp(mode,xt(1)*rhobarra)+0j
        H =  special.hankel1(mode,xt(2)*rhobarra)+0j
        derH =  special.h1vp(mode,xt(2)*rhobarra)+0j
        J2 = special.jn(mode,xt(2)*rhobarra)+0j
        derJ2 = special.jvp(mode,xt(2)*rhobarra)+0j
        return [J,J2,derJ,derJ2,H,derH]  
    
    def a(j,mode):
        if mode!=0:
            rta1 = mode
            rta2 = rhobarra*xt2(j)
        else: 
            rta1 = 0
            rta2 = 1
        return rta1/rta2
    
    def b(j):
        return 1j*mu(j)/xt(j)

    def d(j):
        return 1j*epsi(j)/xt(j)
    
    def f(j):
        return 1j*xz/xt(j)

    Hrho1_tot, Hrho2_tot = 0,0
    Hphi1_tot, Hphi2_tot = 0,0
    
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    exp2 = np.e**(1j*kz_var2*z)   
    for nu in list_modos: 
        coef_A1,coef_C1,coef_B2,coef_D2 = coeficientes(nu)
        
        # cte1 = (-1)**nu
        # cte2 = 1j*pi/2
        # cte3 = (-1j)**nu
        # cte41 = (1j)**nu
        # cte4 = cte41*cte2
        aux = (1j)**nu           
        exp1 = np.e**(1j*nu*phi)
         
        [J,J2,derJ,derJ2,H,derH]  = Bessel(nu)
        
        # Hrho1 = (coef_A1*J*a(1,nu)*epsi(1) + coef_C1*derJ*f(1))*exp1*exp2*cte3
        # Hrho2 =  (coef_B2*H*a(2,nu)*epsi(2) + coef_D2*derH*f(2))*exp1*exp2*cte4
        
        # Hphi1 = (-coef_C1*J*a(1,nu)*xz + coef_A1*derJ*d(1))*exp1*exp2*cte3
        # Hphi2 = (-coef_D2*H*a(2,nu)*xz + coef_B2*derH*d(2))*exp1*exp2*cte4
        
        Hrho1 = (coef_A1*J*a(1,nu)*epsi(1) + coef_C1*derJ*f(1))*exp1*exp2
        Hphi1 = (-coef_C1*J*a(1,nu)*xz + coef_A1*derJ*d(1))*exp1*exp2
        
        Hrho2 =  (a(2,nu)*epsi(2)*(coef_B2*H + Bo*aux*J2) + f(2)*(coef_D2*derH + Ao*aux*derJ2))*exp1*exp2
        Hphi2 = (-a(2,nu)*xz*(coef_D2*H + Ao*aux*J2) + d(2)*(coef_B2*derH + Bo*aux*derJ2))*exp1*exp2
        
        Hrho1_tot = Hrho1_tot + Hrho1
        Hrho2_tot = Hrho2_tot + Hrho2
        
        Hphi1_tot = Hphi1_tot + Hphi1
        Hphi2_tot = Hphi2_tot + Hphi2
    
    return Hrho1_tot,Hrho2_tot,Hphi1_tot,Hphi2_tot

#print('Definir los campos longitudinales')

#Ojo: aca solo van frecuencias reales

def Hz(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        Hz1,Hz2
    -------
    """

    def coeficientes(modo):
        
        coeff = coef(kz_var2,omegac,epsi1,modo,R,hbaramu,Ao,Bo)
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3   
        
        coef_A1 = complex(coef_A1)
        coef_C1 = complex(coef_C1)
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
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
        
        # Hz1 = coef_C1*J*cte3*exp1*exp2
        # Hz2 = coef_D2*H*cte4*exp1*exp2
        aux = (1j)**nu
        Hz1 = coef_C1*J1*exp1*exp2
        Hz2 = (coef_D2*H + Ao*aux*J2)*exp1*exp2
        
        Hz1_tot = Hz1_tot + Hz1
        Hz2_tot = Hz2_tot + Hz2
    
    return Hz1_tot,Hz2_tot


def Htot(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z): 
    """
    Parameters
    ----------
    kz_var2 : kz en 1/micrometros
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
        |Htot|^2 =  |Hrho|^2 + |Hphi|^2 + |Hz|^2
    -------
    """
    
    Hz1, Hz2 = Hz(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    Hrho1,Hrho2,Hphi1,Hphi2 = Ht(kz_var2,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    
    Hrho1_tot, Hrho2_tot = np.abs(Hrho1)**2, np.abs(Hrho2)**2
    Hphi1_tot, Hphi2_tot = np.abs(Hphi1)**2, np.abs(Hphi2)**2
    Hz1_tot, Hz2_tot = np.abs(Hz1)**2, np.abs(Hz2)**2
    
    H1_tot = Hz1_tot + Hrho1_tot + Hphi1_tot
    H2_tot = Hz2_tot + Hrho2_tot + Hphi2_tot
    return H1_tot,H2_tot 

#%%