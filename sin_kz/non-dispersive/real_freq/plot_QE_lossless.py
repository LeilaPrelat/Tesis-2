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

re_epsi1 = 4.9
R = 0.5              #micrones
modo = 4

list_mu =  np.linspace(0.3,0.9,6001)  
path_save = path_basic + '/re_epsi1_%.2f_vs_mu' %(re_epsi1)

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')
        
#%%

label_QE = 'QE approx sin perdidas'
label_QE2 = 'QE approx sin perdidas sin ' + r'$\gamma^2_c$' 
title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) +  ', ' + name_this_py
labelx = '$\mu_c$ [eV]'

im_epsi1_QE = []
omegac_QE = []
im_epsi1_QE_aprox = []
omegac_QE_aprox = []
for mu in list_mu:
    a1 = omegac_cuasi(modo,R,re_epsi1,mu)
    b1 = im_epsi1_cuasi(a1,modo,R,mu) 
    
    if a1.imag == 0:
        a1 = a1.real
    if b1.imag == 0:
        b1 = b1.real         
    
    a2 = omegac_cuasi_aprox(modo,R,re_epsi1,mu)
    b2 = im_epsi1_cuasi_aprox(a2,modo,R,mu) 

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
    info = '.QE lossless sin despreciar gamma_c**2, R=%.1f \mum, Re(epsi1)=%.2f' %(R,re_epsi1) 
    header1 = 'mu [eV]     Im(epsi1)     omega/c' + info + name_this_py
    np.savetxt('QE_lossless_vs_mu_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)
  

plt.figure(figsize=tamfig)
plt.plot(list_mu,im_epsi1_QE,'.-m',ms=10,label=label_QE)
plt.plot(list_mu,im_epsi1_QE_aprox,'-b',lw=5,label=label_QE2)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\epsilon_1$)',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Im_epsi1_vs_mu_%i_QE'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_mu,omegac_QE,'.-m',ms=10,label=label_QE)
plt.plot(list_mu,omegac_QE_aprox,'-b',lw=5,label=label_QE2)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omegac_vs_mu_%i_QE'%(modo))

#%%
