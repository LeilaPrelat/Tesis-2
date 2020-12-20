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
path_graphene = path_basic.replace('/sin_kz/non-dispersive/complex_freq','') 

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene2 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene2)
    from constantes import constantes

pi,hb,c,alfac,mu1,mu2,epsi2 = constantes()

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Usar los datos de Minimizo_perdidas.txt para graficar los campos en la condicion de spaser')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5 #micrones
hbaramu = 0.3        #eV mu_c
modo = 1

name_this_py = os.path.basename(__file__)
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    path_g = path_basic + '/' + 're_epsi1_%.2f' %(re_epsi1)  
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)

try:
    path_load = path_basic + '/' + 're_epsi1_%.2f' %(re_epsi1)  + '/'  + 'find_im_epsi1'
    os.chdir(path_load)
except FileNotFoundError: 
    print('No existe el path_load = lugar en donde se encuentra Minimizo_perdidas.txt')
try:
    data = np.loadtxt('Minimizo_perdidas_modo_%i.txt' %(modo),delimiter = '\t', skiprows=1)
except OSError or IOError:
    print('El archivo Minimizo_perdidas.txt no se encuentra en el path_load')

#%%

data = np.transpose(data)
barrido_mu,rta_Im_epsi1,rta_Im_omega_c,rta_Re_omega_c = data

info = 'R = %.1f$\mu$m, Re($\epsilon_1$) = %.1f, modo = %i, $\mu_c$ = %.1f eV' %(R,re_epsi1,modo,hbaramu)

gain_factor_tot = []
for j in range(len(barrido_mu)):
    crit = rta_Im_epsi1[j]
    omegac = rta_Re_omega_c[j] + 1j*rta_Im_omega_c[j]
    omegac = omegac.real
    epsi1 = re_epsi1 + 1j*crit

    n1 = (epsi1*mu1)**(1/2)
    n1 = n1.real
    gain_factor = -omegac*crit/n1

    gain_factor_tot.append(gain_factor)
    
#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

labely = 'g gain factor [1/$\mu$m]'
title = info + '\n' + name_this_py

plt.figure(figsize=tamfig)
plt.plot(rta_Re_omega_c,gain_factor_tot,'.m',ms=10)
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1)  
if save_graphs==1:
    os.chdir(path_g)
    plt.savefig('gain_factor_vs_omegac_modo%i' %(modo))

plt.figure(figsize=tamfig)
plt.plot(rta_Im_epsi1,gain_factor_tot,'.m',ms=10)
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('Im($\epsilon_1$)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.grid(1) 
if save_graphs==1:
    plt.savefig('gain_factor_vs_im_epsi1_modo%i' %(modo)) 

#%%
