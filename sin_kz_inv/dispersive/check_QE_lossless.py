#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

barrido en mu viendo si el QE sin perdidas (QE_lossless) da lo mismo que 
minimizar Im(omega/c) de QE
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

#%% 

save_data_opt = 1 #guardar data de la minimizacion
save_graphs = 1 #guardar los graficos

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
del path

try:
    sys.path.insert(1, path_basic)
    from complex_omegac_QE import omegac_QE
except ModuleNotFoundError:
    print('complex_omegac_QE.py no se encuentra en ' + path_basic)
    path_basic1 = input('path de la carpeta donde se encuentra complex_omegac_QE.py')
    sys.path.insert(1, path_basic1)
    from complex_omegac_QE import omegac_QE

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic + '/real_freq')
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en ' + path_basic + '/real_freq')
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

#%%

print('Definir parametros del problema')

R = 5              #micrones
# modo = 4

Ep = 0.6
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

list_mu =  np.linspace(0.3,0.9,6001)  
list_modes = [1,2,3,4]
# list_mu = [0.3]

#%%

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')
    
# if modo not in [1,2,3,4]:
#     raise TypeError('Wrong value for mode')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1 or save_graphs ==1:

    path_det = r'/complex_freq/epsiinf_DL_%.2f_vs_mu/R_%.2f' %(epsiinf_DL,R)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)

#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(hbaramu,modo):
    """
    Parameters
    ----------
    hbaramu : potencial quimico del grafeno en eV

    Returns
    -------
    [im(epsilon1)]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu)   
    b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu) 
    return [b]    

#%%    

print('Minimizacion de im(omega/c) de la aprox QE usando como parametro im(epsilon1)')

mu0 = list_mu[0]
tol_NM = 1e-13
ite_NM = 1150

for modo in list_modes:
    
    info1 = 'R = %.2f $\mu$m, $E_p$ = %.3f eV, modo = %i' %(R,Ep,modo)
    info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
    info =  ', ' + info1 + ', ' + info2  + ', ' + name_this_py
    title = info1 +'\n' + info2  + ', ' + name_this_py
    
    cond_inicial = fcond_inicial(mu0,modo)
    
    epsi1_imag_opt = []
    re_omegac_opt = []
    im_omegac_opt = []
    list_mu_opt = []
    
    for mu in list_mu:
        mu = np.round(mu,4)
        print('')
        print(mu)
    
               
        def im_omegac_QE(x):
            rta = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,x,R,mu)
            return rta.imag
            
        
        sol = fsolve(im_omegac_QE,cond_inicial,xtol=tol_NM,maxfev=ite_NM)[0]
    
        omegac = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,sol,R,mu)
        re_omegac_opt.append(omegac.real)
        im_omegac_opt.append(omegac.imag)
        epsi1_imag_opt.append(sol)
        list_mu_opt.append(mu)
            
        cond_inicial = [sol]
    
       
    if save_data_opt==1:
        os.chdir(path)
        print('Guardar data de minimizacion en .txt')
    
        tabla = np.array([list_mu_opt,epsi1_imag_opt,re_omegac_opt,im_omegac_opt])
        tabla = np.transpose(tabla)
        header1 = 'mu [eV]     Im(epsi1)     Re(omega/c) [1/micrones]     Im(omega/c) [1/micrones]' + info
        np.savetxt('opt_QE_vs_mu_modo%i.txt' %(modo), tabla, fmt='%1.11e', delimiter='\t', header = header1)
    
    
    im_epsi1_QE_lossless = []
    omegac_QE_lossless = []
    for mu in list_mu:
        a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,mu)
        if a.imag == 0:
            a = a.real 
        b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,mu)
        if b.imag == 0:
            b = b.real     
        im_epsi1_QE_lossless.append(b)
        omegac_QE_lossless.append(a)
    
    label_QE1 = 'QE sin perdidas'
    label_QE2 = 'QE minimizando im(omega/c)'
    labelx = '$\mu_c$ [eV]'
    
    plt.figure(figsize=tamfig)
    plt.plot(list_mu,im_epsi1_QE_lossless,'.-m',ms=15,label=label_QE1)
    plt.plot(list_mu_opt,epsi1_imag_opt,'.-b',ms=10,label=label_QE2)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'$\epsilon_{ci}$',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        os.chdir(path)
        plt.savefig('Im_epsi1_vs_mu_QE_%i'%(modo))
    
    plt.figure(figsize=tamfig)
    plt.plot(list_mu,omegac_QE_lossless,'.-m',ms=15,label=label_QE1)
    plt.plot(list_mu_opt,re_omegac_opt,'.-b',ms=10,label=label_QE2)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'Re($\omega/c$) [1/$\mu$m]',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        plt.savefig('re_omegac_vs_mu_QE_%i'%(modo))
        
    plt.figure(figsize=tamfig)
    plt.plot(list_mu_opt,im_omegac_opt,'.m',ms=10,label=label_QE2)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'Im($\omega/c$) [1/$\mu$m]',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        plt.savefig('im_omegac_vs_mu_QE_%i'%(modo))

#%%
