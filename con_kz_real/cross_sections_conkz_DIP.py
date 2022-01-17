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
    con un dipolo interior en el cilindro
    
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
    from coef_matrix_AoBo_DIP import coef
except ModuleNotFoundError:
    print('coef_matrix_AoBo.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra coef_matrix_AoBo_DIP.py')
    sys.path.insert(1, path_basic)
    from coef_matrix_AoBo_DIP import coef
    
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

def Qscat(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo): 
    """
    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    
    hbaramu: potencial quimico del grafeno en eV

    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez

    Returns
        Qscat/c adimensional
    -------
    """
    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
    
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

def Qabs(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo): 
    """
    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    
    hbaramu: potencial quimico del grafeno en eV
    
    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez

    Returns
        Qabs/c adimensional
    -------
    """

    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
    
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

def Qext(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo): 
    """
    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    
    hbaramu: potencial quimico del grafeno en eV

    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez

    Returns
        Qext/c adimensional
    -------
    """

    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    xt2 = -xz**2 + mu2*epsi2

    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
    
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

def Qscat2(kz_var,omegac,epsi1,nmax,R,hbaramu,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo): 
    """
    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    epsi1 : permeabilidad electrica del medio 1
    
    nmax: sumatoria en modos desde -nmax hasta +nmax (2*nmax+1 modos)
    
    hbaramu: potencial quimico del grafeno en eV
    rho : coordenada radial en micrometros
    phi : coordenada azimutal
    z : coordenada longitudinal en micrometros
    
    p1 = p+ = px + 1j*py
    p2 = p- = px - 1j*py
    pz 
    
    rho_D : posicion rho del dipolo en micrometros    
    theta_D : posicion theta del dipolo en radianes    
    z_D : posicion z del dipolo en micrometros    

    Ao : Amplitud del campo incidente en Hz
    Bo : Amplitud del campo incidente en Ez

    Returns
        formula del Qscat antes de integrar en theta (pregunta de Claudio
                                                      en la defensa)
    -------
    """
    k0 = omegac #=omega/c
    xz = kz_var/k0
    Rbarra = R*k0 #adimensional
    rhobarra = rho*k0
    zbar = z*k0
    xt2 = -xz**2 + mu2*epsi2

    if k0.imag != 0:
        raise TypeError('input: Freq reales, no complejas')
        
    if kz_var.imag != 0:
        raise TypeError('input: kz real, no complejo')
        
    def coef_an_bn(mode):
        
        coeff = coef(kz_var,omegac,epsi1,mode,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
    
        coef_A1,coef_B2,coef_C1,coef_D2 = coeff  #intercambio de la columna 2 con 3
        
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
    
    cte = (2/np.pi)**(1/2)
    
    for nu in list_modos: 
        nu = int(nu)
        for nu2 in list_modos: #nu prima
            nu2 = int(nu2)

            den1 = xt2*((xt2*rhobarra)**(1/2))
            den2 = (xt2*rhobarra)**(1/2)
            xt3 = xt2.conjugate()
            den3 = xt3*((xt3*rhobarra)**(1/2))
                
            Dn,Bn = coef_an_bn(nu)
            Dn2,Bn2 = coef_an_bn(nu2)
            Dn2A = np.abs(Dn)**2
            Bn2A = np.abs(Bn)**2
            
            Dn2B = np.abs(Dn2)**2 #nu prima
            Bn2B = np.abs(Bn2)**2 #nu prima
            Dn2B = Dn2B.conjugate()
            Bn2B = Bn2B.conjugate()
            
            exp0 = np.e**(1j*nu*phi)
            exp1 = np.e**(1j*zbar*xz)
            exp2 = np.e**(1j*(xt2*rhobarra -nu*np.pi/2 - np.pi/4))            
            expA = exp0*exp1*exp2

            exp3 = np.e**(-1j*nu2*phi)
            exp4 = np.e**(-1j*zbar*xz)
            exp5 = np.e**(-1j*(xt3*rhobarra -nu2*np.pi/2 - np.pi/4))
            expB = exp3*exp4*exp5            
            

            term1A = Dn2A*expA*mu2/den1
            term1B = Dn2B*expB/den2
            
            term1 = term1A*term1B
            
            term2A = Bn2A*expA/den2
            term2B = Bn2B*expB*epsi2/den3
    
            term2 = term2A*term2B
            
            Qscat = Qscat + term1 - term2
     
    ctef = cte/(4*pi*Rbarra*xt2)
    Qscatf = Qscat*ctef
    
    return Qscatf # ---> esta dividido por c el resultado