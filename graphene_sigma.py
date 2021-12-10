#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

conductividad del grafeno  ---> IMPORTANTE. Importar en todas las funciones

"""
import numpy as np

#%%

def sigma(hbw,hbmu,hbgama):
    """
    Parameters
    ----------
    hbw : energia = hbar*omega en eV
    hbmu : potencial quimico del grafeno mu (en eV)
    hbgama : collision frequency in eV 

    Returns
    -------
    devuelve la conductividad del grafeno
    en unidades de e^2/hbar 
	---> multiplicar despues por fine_structure*c (unidades gaussianas)

    """
    global Tk  ### variables globales ---> no pueden cambiarse
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  
    Tk = 300 			### temperature in K ---> no cambia
    
    TKmu = Tk*akb/hbmu
    if TKmu==0.:
        # intra = 1j*(1/np.pi)*abs(hbmu)/(hbw+1j*hbgama)
        # if hbw-2*abs(hbmu)>0:
        #     inter1 = 0.25
        # else: 
        #     inter1 = 0
         
        # inter = inter1 + 1j*np.log((hbw-2*abs(hbmu))**2/(hbw+2*abs(hbmu))**2)/(4*np.pi)
        raise ValueError('no se la formula para este caso')
    elif TKmu!=0:

        aux2 = hbw + 1j*hbgama
        aux3 = (aux2 - 2*hbmu)**2 + (2*Tk*akb)**2
        aux4 = (aux2 + 2*hbmu)**2 
        aux5 = (aux2 - 2*hbmu)/(2*akb*Tk)   
        inter = 0.25*(0.5 + np.arctan(aux5)/np.pi - 1j*np.log(aux4/aux3)/(2*np.pi))   

        auxTkmu = hbmu/(2*Tk*akb)
        aux = np.log(np.exp(auxTkmu) + np.exp(-auxTkmu))
        intra_aux = aux/aux2
        intra = 2j*akb*Tk*intra_aux/np.pi  

    sigmatot = inter + intra
    return sigmatot, inter, intra  

#%%

