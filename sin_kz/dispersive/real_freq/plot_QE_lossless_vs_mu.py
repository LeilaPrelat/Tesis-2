#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

comparar la solucion QE sin perdidas con la misma formula
	despreciando gamma_c**2 (que queda como -mu^{-1/2})
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt

#%% 

save_graphs = 1 #guardar los graficos
save_data_opt = 1

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

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi,im_epsi1_cuasi_aprox,omegac_cuasi_aprox
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi,im_epsi1_cuasi_aprox,omegac_cuasi_aprox

#%%

print('Definir parametros del problema')

R = 0.05              #micrones
modo = 4

Ep = 0.6
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

list_mu =  np.linspace(0.3,0.9,6001)  
path_save = path_basic + '/R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)

info1 = 'R = %.2f $\mu$m, $E_p$ = %.3f eV, modo = %i' %(R,Ep,modo)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
title = info1 +'\n' + info2  + ', ' + name_this_py
info = info1 + ', ' + info2  + ', ' + name_this_py

#%%

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

label_graph = 'Opt det sin kz'
label_QE = 'QE approx sin perdidas'
label_QE2 = 'QE approx sin perdidas sin ' + r'$\gamma^2_c$' 
labelx = '$\mu_c$ [eV]'

im_epsi1_QE = []
omegac_QE = []
im_epsi1_QE_aprox = []
omegac_QE_aprox = []
for mu in list_mu:
    a1 = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,mu)
    b1 = im_epsi1_cuasi(a1,Ep,epsiinf_DL,gamma_DL,modo,R,mu) 
    
    if a1.imag == 0:
        a1 = a1.real
    if b1.imag == 0:
        b1 = b1.real         
    
    a2 = omegac_cuasi_aprox(modo,epsiinf_DL,gamma_DL,R,mu)
    b2 = im_epsi1_cuasi_aprox(a2,epsiinf_DL,gamma_DL,modo,R,mu) 

    if a2.imag == 0:
        a2 = a2.real
    if b2.imag == 0:
        b2 = b2.real        

    im_epsi1_QE.append(b1)
    omegac_QE.append(a1)
    im_epsi1_QE_aprox.append(b2)
    omegac_QE_aprox.append(a2)

#%%

if save_data_opt==1:
    print('Guardar data de minimizacion en .txt')
    os.chdir(path_save)
    tabla = np.array([list_mu,omegac_QE,im_epsi1_QE])
    tabla = np.transpose(tabla)
    header1 = 'mu [eV]     Im(epsi1)     omega/c' + info 
    np.savetxt('QE_lossless_vs_mu_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)

#%%

plt.figure(figsize=tamfig)
plt.plot(list_mu,im_epsi1_QE,'.m',ms=10,label=label_QE)
plt.plot(list_mu,im_epsi1_QE_aprox,'.b',ms=10,label=label_QE2)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\epsilon_1$)',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('Im_epsi1_vs_mu_%i_QE'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_mu,omegac_QE,'.m',ms=10,label=label_QE)
plt.plot(list_mu,omegac_QE_aprox,'.b',ms=10,label=label_QE2)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omegac_vs_mu_%i_QE'%(modo))

#%%
