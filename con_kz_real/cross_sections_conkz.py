#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar las formulas del Qscat y Qabs con kz
(que aparece en el overleaf
           en la seccion 4.8)

ESTABAN MAL ESAS FORMULAS DEL CUADERNO (movidas a extra)
VER TESIS
(revisar cuentas)

"""
import numpy as np
import os 
import sys

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
    print('coef_matrix_AoBo.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra coef_matrix_AoBo.py')
    sys.path.insert(1, path_basic)
    from coef_matrix_AoBo import coef
    
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

#print('Definir el coeficiente Qscat')

#Ojo: aca solo van frecuencias reales

def Qscat(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao: pol p dentro de Hz
    Bo: pol s dentro de Ez

    Returns
        Qscat/c adimensional
    -------
    """
    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,Ao,Bo)
    
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3
        
        # cte1 = (-1)**mode
        # cte2 = 1j*pi/2
        
        # coef_A1 = coef_A1/cte1
        # coef_C1 = coef_C1/cte1
        # coef_B2 = coef_B2/cte2
        # coef_D2 = coef_D2/cte2    
        
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
        
        # if Ao == 0:
        #     coef_D2 = 0
        # if Bo== 0:
        #     coef_B2 = 0
                 
        # print(coef_D2,coef_B2)
        return coef_D2,coef_B2

    #xt2 = -xz**2 + mu2*epsi2 #xt**2

    Qscat = 0
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    # list_modos = np.linspace(0,nmax,nmax+1)
    # list_modos =  [nmax]
    
    for nu in list_modos: 
        nu = int(nu)
        
        Dn,Bn = coef_an_bn(nu)
        Dn2 = np.abs(Dn)**2
        Bn2 = np.abs(Bn)**2

        Qscat = Qscat + mu2*Dn2 - epsi2*Bn2
     
    ctef = 1/(4*pi*Rbarra*xt2)
    Qscatf = Qscat*ctef
    
    return Qscatf # ---> esta dividido por c el resultado

#%%

def Qabs(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao: pol p dentro de Hz
    Bo: pol s dentro de Ez

    Returns
        Qabs/c adimensional
    -------
    """

    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,Ao,Bo)
    
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3
        
        # cte1 = (-1)**mode
        # cte2 = 1j*pi/2
        
        # coef_A1 = coef_A1/cte1
        # coef_C1 = coef_C1/cte1
        # coef_B2 = coef_B2/cte2
        # coef_D2 = coef_D2/cte2    
        
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
        
        # if Ao == 0:
        #     coef_D2 = 0
        # if Bo== 0:
        #     coef_B2 = 0
                 
        # print(coef_D2,coef_B2)
        return coef_D2,coef_B2

    #xt2 = -xz**2 + mu2*epsi2 #xt**2

    Qabs = 0
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    # list_modos = np.linspace(0,nmax,nmax+1)
    # list_modos =  [nmax]
    
    for nu in list_modos: 
        nu = int(nu)
        aux = (-1j)**nu
        
        Dn,Bn = coef_an_bn(nu)
        Aoconj,Boconj = Ao.conjugate(), Bo.conjugate()
        Dn2 = np.abs(Dn)**2
        Bn2 = np.abs(Bn)**2
        term1 = mu2*(Dn2 + (aux*Aoconj*Dn).real)
        term2 =  epsi2*(Bn2 + (aux*Boconj*Bn).real)

        Qabs = Qabs - term1 + term2
     
    ctef = 1/(4*pi*Rbarra*xt2)
    Qabstf = Qabs*ctef
    
    return Qabstf # ---> esta dividido por c el resultado

#%%

# Qext = Qabs + Qscat, tengo miedo de que suceda cancelacion catastrofica, asi que lo voy a escribir a mano

def Qext(kz_var,omegac,epsi1,nmax,R,hbaramu,Ao,Bo): 
    """
    Parameters
    ----------
    kz : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    R: radio del cilindro en micrometros
    hbaramu: potencial quimico del grafeno en eV
    
    Ao: pol p dentro de Hz
    Bo: pol s dentro de Ez

    Returns
        Qext/c adimensional
    -------
    """

    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,Ao,Bo)
    
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3
        
        # cte1 = (-1)**mode
        # cte2 = 1j*pi/2
        
        # coef_A1 = coef_A1/cte1
        # coef_C1 = coef_C1/cte1
        # coef_B2 = coef_B2/cte2
        # coef_D2 = coef_D2/cte2    
        
        coef_D2 = complex(coef_D2)
        coef_B2 = complex(coef_B2)
        
        
        # if Ao == 0:
        #     coef_D2 = 0
        # if Bo== 0:
        #     coef_B2 = 0
                 
        # print(coef_D2,coef_B2)
        return coef_D2,coef_B2

    #xt2 = -xz**2 + mu2*epsi2 #xt**2

    Qext = 0
    list_modos = np.linspace(-nmax,nmax,2*nmax+1)
    # list_modos = np.linspace(0,nmax,nmax+1)
    # list_modos =  [nmax]
    
    for nu in list_modos: 
        nu = int(nu)
        aux = (-1j)**nu
        
        Dn,Bn = coef_an_bn(nu)
        Aoconj,Boconj = Ao.conjugate(), Bo.conjugate()

        term1 = mu2*((aux*Aoconj*Dn).real)
        term2 =  epsi2*((aux*Boconj*Bn).real)

        Qext = Qext - term1 + term2
     
    ctef = 1/(4*pi*Rbarra*xt2)
    Qexttf = Qext*ctef
    
    return Qexttf # ---> esta dividido por c el resultado



#%%
