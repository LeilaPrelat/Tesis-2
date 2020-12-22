#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

gain factor definido en el paper
    de Dionne
    
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/non-dispersive','') 
path_save = path_basic + '/' + 'gain_factor'
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene)
    from constantes import constantes

pi,hb,c,alfac,mu1,mu2,epsi2 = constantes()

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
modo = 4
R = 0.5         #micrones 
index = 0
    
#%%

if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1)
os.chdir(path_load)
name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en el path_load')

data_load = np.transpose(data_load)
[barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load

hbaramu = barrido_mu[index]
omegac0 = omegac_opt[index]
im_epsi1c = epsi1_imag_opt[index]
epsi1 = re_epsi1 + 1j*im_epsi1c

info1 = 'R = %.1f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_1$ = %.1f - i%.2e' %(R,hbaramu,re_epsi1,-im_epsi1c)
info2 = ' y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(omegac0,modo)

title2 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f, modo = %i' %(R,re_epsi1,modo)
title2 = title2 + ', ' + name_this_py

#%%

print('Calcular gain factor para diferentes Im(epsi1)')

tol = 1e-3
list_im_epsi1_fino = [0,im_epsi1c+3*tol,im_epsi1c+2*tol,im_epsi1c+tol,im_epsi1c,im_epsi1c-tol,im_epsi1c-2*tol]

N = int(1e3)               
omegac1,omegac2 = omegac0*0.95,omegac0*1.05
# lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
list_omegac = np.linspace(omegac1,omegac2,N)

list_gain_factor_tot = []
for im_epsi1 in list_im_epsi1_fino:
    im_epsi1 = np.round(im_epsi1,7)
    print(im_epsi1)
    epsi1 = re_epsi1 + 1j*im_epsi1
    list_gain_factor = []
    for omeggac in list_omegac:
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omeggac*im_epsi1/n1
        list_gain_factor.append(gain_factor)  
    list_gain_factor_tot.append(list_gain_factor)  

#%%

colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']
labely = 'g gain factor [1/$\mu$m]'
title = info1 +'\n' + info2  + ', ' + name_this_py

plt.figure(figsize=tamfig)
for k in range(len(list_gain_factor_tot)):
    list_gain_factor = list_gain_factor_tot[k]
    im_epsi1 = list_im_epsi1_fino[k]
    if im_epsi1 == im_epsi1c:
        labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
    elif im_epsi1 == 0:
        labell = 'Im($\epsilon_1$) = 0'  
    else:
        labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)

    plt.plot(list_omegac,list_gain_factor,'o',color = colores[k],ms = 4,alpha = 0.8,label = labell)
plt.title(title,fontsize=tamtitle)
mini,maxi = np.min(list_gain_factor_tot), np.max(list_gain_factor_tot)
plt.plot(omegac0*np.ones(10),np.linspace(mini,maxi,10),'k-')
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1)  
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor_modo%i_mu%.4f' %(modo,hbaramu), format='png') 

#%%


gain_factor_tot = []
for j in range(len(barrido_mu)):
    crit = epsi1_imag_opt[j]
    omegac = omegac_opt[j] 
    epsi1 = re_epsi1 + 1j*crit

    n1 = (epsi1*mu1)**(1/2)
    n1 = n1.real
    gain_factor = -omegac*crit/n1

    gain_factor_tot.append(gain_factor)
    
plt.figure(figsize=tamfig)
plt.plot(omegac_opt,gain_factor_tot,'.m',ms=10)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\omega/c$ $[1/\mu m]$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1) 
if save_graphs==1:
    plt.savefig('gain_factor_vs_omegac_modo%i' %(modo)) 

plt.figure(figsize=tamfig)
plt.plot(epsi1_imag_opt,gain_factor_tot,'.m',ms=10)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('Im($\epsilon_1$)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1) 
if save_graphs==1:
    plt.savefig('gain_factor_vs_im_epsi1_modo%i' %(modo)) 
    
plt.figure(figsize=tamfig)
plt.plot(barrido_mu,gain_factor_tot,'.m',ms=10)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1) 
if save_graphs==1:
    plt.savefig('gain_factor_vs_mu_modo%i' %(modo)) 
    
    
#%%
