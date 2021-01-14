#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

coeficiente an

"""
import scipy.special as sp #bessel functions
import sys
import os

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_graphene = path.replace('/' + name_this_py,'')

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_graphene)
    path_graphene = input('path de la carpeta donde se encuentra graphene_sigma.py')
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

def an_nondisp(omegac,epsi1,modo,R,hbaramu,Ao): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    modo: modo (entero)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao: pol p ---> Hz (Amplitud del campo incidente en Hz del medio 2)

    Returns
        coeficiente complejo an (el del campo Hz)
        para el caso no dispersivo 
    -------
    """
    energy = omegac*c*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    k0 = omegac
    Rbarra = R*k0 #adimensional
    x1,x2 = (epsi1*mu1)**(1/2)+0j,(epsi2*mu2)**(1/2)+0j
        
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
    
    return an

#%%

def an_disp(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,hbaramu,Ao): 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrometros
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    epsi_ci: Im(epsilon1) del colorante (dye)
    modo : modo
    R : radio del cilindro en micrometros
    hbaramu : potencial quimico del grafeno en eV
    
    Ao: pol p ---> Hz (Amplitud del campo incidente en Hz del medio 2)

    Returns
        coeficiente complejo an (el del campo Hz)
        para el caso dispersivo
    -------
    """
    energy = omegac*c*hb

    sigmatot, inter, intra = sigma(energy,hbaramu,hbargama) 
        
    num = Ep**2
    den = energy**2 + 1j*gamma_DL*energy
    epsi1 = epsiinf_DL - num/den + 1j*epsi_ci
    
    k0 = omegac
    Rbarra = R*k0 #adimensional
    
    x1,x2 = (epsi1*mu1)**(1/2), (epsi2*mu2)**(1/2)
        
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
    
    return an

#%%
