#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

plotear la diferencia entre los valores de Im(epsilon1)
    del caso kz y del caso sin kz
    
"""

import numpy as np
import os 
import matplotlib.pyplot as plt

save_graphs = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'resta_conkz_sinkz'
path_sinkz = path_basic.replace('/' + 'con_kz_real','') + '/' + 'sin_kz' + '/' + 'non-dispersive/real_freq'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5              #micrones
kz_real = 0    
modo = 4

#%%

path_sinkz = path_sinkz + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1) 
path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1) 

#%%
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

print('')
print('Importar los valores de SPASER')

os.chdir(path_load)
name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz_real,modo)
try:
    data_conkz = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)
print('')

data_conkz = np.transpose(data_conkz)
[list_mu_kz,omegac_opt,epsi1_imag_opt,eq_det] = data_conkz

os.chdir(path_sinkz)
name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_sinkz)
data_load = np.transpose(data_load)
[barrido_mu_sinkz,omegac_opt_sinkz,epsi1_imag_opt_sinkz,eq_det_sinkz] = data_load

if list_mu_kz[0]  != barrido_mu_sinkz[0]:
    raise TypeError('Los barridos no empiezan del mismo mu')

#%%

title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f, kz = %.4f 1/ $\mu$m' %(modo,R,re_epsi1,kz_real) 
title2 = title + '\n' + name_this_py
labelx = '$\mu_c$ [eV]'
labely = '(Im($\epsilon_1$) con kz - Im($\epsilon_1$) sin kz)'
#print('Graficar el campo |Etot|^2 para el medio 1 y 2')

resta_values = []
for j in range(len(list_mu_kz)):
    im_epsi1_conkz = epsi1_imag_opt[j]
    im_epsi1_sinkz = epsi1_imag_opt_sinkz[j]
    resta_values.append((im_epsi1_conkz-im_epsi1_sinkz))

plt.figure(figsize=tamfig)
plt.plot(list_mu_kz,resta_values,'.',ms=10,label = 'modo = %i' %(modo))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('resta_kz%.4f_modo%i' %(kz_real,modo), format='png')

#%%

