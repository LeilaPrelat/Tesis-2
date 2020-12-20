#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 21:00:49 2020

@author: leila
"""
import numpy as np
import matplotlib.pyplot as plt
import os

#%%

save_graphs = 1


path = os.path.abspath(__file__) #path absoluto del .py actual
name_this_py = os.path.basename(__file__)
path = path.replace(name_this_py,'')
os.chdir(path) 

name_this_py = '. Ver ' + name_this_py

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

print('Constantes universales')

pi = np.pi
c = 3*10**(14)             ### light velocity in micron/seg
alfac = 1/137.0359            ### fine structure
akb = 8.6173324E-5           ### Boltzman constant in eV/K  
hb = 6.58211899*10**(-16)     ### Planck constant hbar in eV*s

def sigma(x,hbmu,hbgama): #x = hbw/mu_c  #sigma del libro de depine: comparar con el del paper de tunable
    hbw = x*hbmu
    Tk = 300               # temperature in K
    TKmu = Tk*akb/hbmu
    if TKmu==0.:
        intra = 1j*(1/np.pi)*abs(hbmu)/(hbw+1j*hbgama)
        if hbw-2*abs(hbmu)>0:
            inter = 0.25
        else:  
            inter = 0
         
        inter = inter + 1j*np.log((hbw-2*abs(hbmu))**2/(hbw+2*abs(hbmu))**2)/(4*np.pi)
        
    elif TKmu!=0:

        aux2 = hbw + 1j*hbgama
        aux3 = (aux2 - 2*hbmu)**2 + (2*Tk*akb)**2
        aux4 = (aux2 + 2*hbmu)**2 
        aux5 = (aux2 - 2*hbmu)/(2*akb*Tk)   
        inter = 0.25*(0.5 + np.arctan(aux5)/pi - 1j*np.log(aux4/aux3)/(2*pi))   

        auxTkmu = hbmu/(2*Tk*akb)
        aux = np.log(np.exp(auxTkmu) + np.exp(-auxTkmu))
        intra_aux = aux/aux2
        intra = 2j*akb*Tk*intra_aux/pi        


    sigmatot = inter + intra
    
    return sigmatot, inter, intra

def sigma2(omega_c,hbmu,hbgama): #sigma del libro de depine: comparar con el del paper de tunable
    
    hbw = omega_c*c*hb
    Tk = 300               # temperature in K
    TKmu = Tk*akb/hbmu
    if TKmu==0.:
        intra = 1j*(1/np.pi)*abs(hbmu)/(hbw+1j*hbgama)
        if hbw-2*abs(hbmu)>0:
            inter = 0.25
        else:  
            inter = 0
         
        inter = inter + 1j*np.log((hbw-2*abs(hbmu))**2/(hbw+2*abs(hbmu))**2)/(4*np.pi)
        
    elif TKmu!=0:
        
        aux2 = hbw + 1j*hbgama
        aux3 = (aux2 - 2*hbmu)**2 + (2*Tk*akb)**2
        aux4 = (aux2 + 2*hbmu)**2 
        aux5 = (aux2 - 2*hbmu)/(2*akb*Tk)   
        inter = 0.25*(0.5 + np.arctan(aux5)/pi - 1j*np.log(aux4/aux3)/(2*pi))   

        auxTkmu = hbmu/(2*Tk*akb)
        aux = np.log(np.exp(auxTkmu) + np.exp(-auxTkmu))
        intra_aux = aux/aux2
        intra = 2j*akb*Tk*intra_aux/pi   
        
    sigmatot = inter + intra
    ### Total (Inter+Intra) conductivity, units of e**2/hb 
    ###     
    return sigmatot, inter, intra

def sigmaMAL(omega_c,hbmu,hbgama): #sigma del paper de tunable: comparar con el del libro de depine
    
    hbw = omega_c*c*hb
    Tk = 300               # temperature in K
    TKmu = Tk*akb/hbmu
    if TKmu==0.:
        intra = 1j*(1/np.pi)*abs(hbmu)/(hbw+1j*hbgama)
        if hbw-2*abs(hbmu)>0:
            inter=0.25
        else: 
            inter=0
         
        inter = inter + 1j*np.log((hbw-2*abs(hbmu))**2/(hbw+2*abs(hbmu))**2)/(4*np.pi)
        
    elif TKmu!=0:
        TK = TKmu*hbmu/akb
        aux = np.exp(1/(2*TKmu))+np.exp(-1/(2*TKmu))
        intra = 2j*(1/np.pi)*akb*TK*np.log(aux)/(hbw+1j*hbgama)
        aux2 = (hbw-2*hbmu)/(2*TK*akb)
        aux3 = (hbw-2*hbmu)**2 
        aux4 = (2*TK*akb)**2
        aux5 = aux3 + aux4
        inter = 0.25*(0.5+np.arctan(aux2)/np.pi-1j*np.log(aux3/aux5)/(2*np.pi))   
        
    sigmatot = inter+intra
    
    return sigmatot, inter, intra

#%%

hbargama = 0.0001      # collision frequency in eV
hbaramu = 0.3           #eV mu_c

n = int(1e3)
list_E = np.linspace(0.1,5,n)
list_x = np.linspace(0.15,6,n)
list_x2 = np.linspace(0.15,20,n)

sigma_grafeno = sigma(list_x,hbaramu,hbargama)[0]
sigma_grafeno2 = sigma2(list_x2,hbaramu,hbargama)[0]
sigma_grafenoMAL2 = sigmaMAL(list_x2,hbaramu,hbargama)[0]
#sigma_grafeno = np.array(sigma_grafeno)*alfac*c

plt.figure(figsize=tamfig)
plt.plot(list_x,sigma_grafeno.real,'.r',ms = 8,label = 'Re($\sigma$)')
plt.plot(list_x,sigma_grafeno.imag,'.m',ms = 8,label = 'Im($\sigma$)')
plt.title('Grafeno con $\mu_c$ = %.1f eV, $\gamma$ = %.4f eV' %(hbaramu,hbargama)+name_this_py,fontsize=tamtitle)
plt.ylabel('Conductividad',fontsize=tamletra)
plt.xlabel('$\hbar \omega /\mu_c$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('grafeno-depine')

plt.figure(figsize=tamfig)
plt.plot(list_x2,sigma_grafeno2.real,'.r',ms = 8,label = 'libro Re($\sigma$)')
plt.plot(list_x2,sigma_grafenoMAL2.real,'.g',ms = 8,label = 'tunable Re($\sigma$)')
plt.plot(list_x2,sigma_grafeno2.imag,'.m',ms = 8,label = 'libro Im($\sigma$)')
plt.plot(list_x2,sigma_grafenoMAL2.imag,'.b',ms = 8,label = 'tunable Im($\sigma$)')
plt.title('Grafeno con $\mu_c$ = %.1f eV, $\gamma$ = %.4f eV' %(hbaramu,hbargama)+name_this_py,fontsize=tamtitle)
plt.ylabel('Conductividad',fontsize=tamletra)
plt.xlabel('$\omega/c$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=3,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('grafeno-depine2')

    



#%%
    
    
    
    
    
    
