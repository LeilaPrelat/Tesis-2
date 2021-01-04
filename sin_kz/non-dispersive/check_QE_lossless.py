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

re_epsi1 = 4.9
R = 0.5              #micrones
modo = 4

list_mu =  np.linspace(0.3,0.9,6001)  
# list_mu = [0.3]

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1 or save_graphs ==1:

    path_det = r'/complex_freq/re_epsi1_%.2f_vs_mu' %(re_epsi1)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)

#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(hbaramu):
    """
    Parameters
    ----------
    hbaramu : potencial quimico del grafeno en eV

    Returns
    -------
    [im(epsilon1)]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_cuasi(modo,R,re_epsi1,hbaramu)
    if a.imag == 0:
        a = a.real    
    b = im_epsi1_cuasi(a,modo,R,hbaramu) 
    if b.imag == 0:
        b = b.real
    return [b]    

#%%    

print('Minimizacion de im(omega/c) de la aprox QE usando como parametro im(epsilon1)')

mu0 = list_mu[0]
cond_inicial = fcond_inicial(mu0)

epsi1_imag_opt = []
re_omegac_opt = []
im_omegac_opt = []
list_mu_opt = []

tol_NM = 1e-13
ite_NM = 1150

for mu in list_mu:
    mu = np.round(mu,4)
    print('')
    print(mu)

           
    def im_omegac_QE(x):
        rta = omegac_QE(x,modo,re_epsi1,R,mu)
        return rta.imag
        
    
    sol = fsolve(im_omegac_QE,cond_inicial,xtol=tol_NM,maxfev=ite_NM)[0]

    epsi1 = re_epsi1 + 1j*sol
    omegac = omegac_QE(sol,modo,re_epsi1,R,mu)
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
    info = '.Opt det SIN kz, R=%.1f \mum, Re(epsi1)=%.2f' %(R,re_epsi1) 
    header1 = 'mu [eV]     Im(epsi1)     Re(omega/c) [1/micrones]     Im(omega/c) [1/micrones]' + info + ', ' + name_this_py
    np.savetxt('opt_QE_vs_mu_modo%i.txt' %(modo), tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

im_epsi1_QE_lossless = []
omegac_QE_lossless = []
for mu in list_mu:
    a = omegac_cuasi(modo,R,re_epsi1,mu)
    b = im_epsi1_cuasi(a,modo,R,mu) 
    im_epsi1_QE_lossless.append(b)
    omegac_QE_lossless.append(a)

label_QE1 = 'QE sin perdidas'
label_QE2 = 'QE minimizando im(omega/c)'
title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) +  ', ' + name_this_py
labelx = '$\mu_c$ [eV]'

plt.figure(figsize=tamfig)
plt.plot(list_mu,im_epsi1_QE_lossless,'.-m',ms=15,label=label_QE1)
plt.plot(list_mu_opt,epsi1_imag_opt,'.-b',ms=10,label=label_QE2)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\epsilon_1$)',fontsize=tamletra)
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
