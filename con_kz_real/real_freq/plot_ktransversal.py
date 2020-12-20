#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)

determinante de 4x4

barrido en kz para mu = 0.3.  Hallar omega/c real y im(epsilon1) que minimizan
el determinante con kz para diferentes valores de kz
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt

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
path_basic2 = path_basic.replace('/' + 'real_freq','')

try:
    sys.path.insert(1, path_basic2)
    from def_kt import kt
except ModuleNotFoundError:
    print('def_kt.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra def_kt.py')
    sys.path.insert(1, path_basic2)
    from def_kt import kt

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
kz_real = 1       #eV mu_c
modo = 4
    
kzlim = 0.14

path_load = path_basic + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1)
path_save = path_load

#%%

if R != 0.5:
    raise TypeError('Los .txt fueron hechos para R = 0.5 micrones')

#%%

os.chdir(path_load)
name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz_real,modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data_load = np.transpose(data_load)
[list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
#%%

label_graph = 'Opt det con kz = %.4f 1/$\mu$m' %(kz_real)
labelsinkz = 'Opt sin kz'
labelx = '$\mu_c$ [eV]'
title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) 
title2 = title + '\n' + name_this_py

kt_imag = [] # Imag(kt)
kt_real = [] # Real(kt)
kt_abs = [] #modulo del k transversal |kt|

for j in range(len(list_mu_opt)):
    hbaramu = list_mu_opt[j]
    omegac = omegac_opt[j]
    epsi1_imag = epsi1_imag_opt[j]
    epsi1 = re_epsi1 + 1j*epsi1_imag
    kt_value = kt(kz_real,omegac,epsi1,modo,R,hbaramu)    
    kt_abs.append(np.abs(kt_value))
    kt_imag.append(kt_value.imag)
    kt_real.append(kt_value.real)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_abs,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|kt|_vs_mu_kz%.4f_%i'%(kz_real,modo), format='png')
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_real,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel('Re($k_t$) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('kt_real_vs_mu_kz%.4f_%i'%(kz_real,modo), format='png')
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_imag,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('kt_imag_vs_mu_kz%.4f_%i'%(kz_real,modo), format='png')
    
#%%
