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
path_det = path_basic.replace('/real_freq','') 
path_save = path_basic + '/' + 'plot_det_sinkz_2D'
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_det)
    from det_sinkz import determinante
except ModuleNotFoundError:
    print('det_sinkz.py no se encuentra en ' + path_det)
    path_basic2 = input('path de la carpeta donde se encuentra det_sinkz.py')
    sys.path.insert(1, path_det)
    from det_sinkz import determinante

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

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
modo = 2
R = 0.5         #micrones 

#%%
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1)
os.chdir(path_load)
name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data_load = np.transpose(data_load)
[barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load

index = 0
hbaramu = barrido_mu[index]
im_epsi1c = epsi1_imag_opt[index]
omegac = omegac_opt[index] 
epsi1 = re_epsi1 + 1j*im_epsi1c

info1 = 'R = %.1f $\mu$m, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f - i%.2e' %(R,hbaramu,re_epsi1,-im_epsi1c)
info2 = ' y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(omegac,modo)
title = info1 +'\n' + info2 + name_this_py

#%%

def fcond_inicial(hbaramu):
    """
    Parameters
    ----------
    hbaramu : potencial quimico del grafeno en eV

    Returns
    -------
    [re(omega/c), im(epsilon1)]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_cuasi(modo,R,re_epsi1,hbaramu)
    b = im_epsi1_cuasi(a,modo,R,hbaramu) 
    return [a,b]  

#%%

print('Calcular el |det sin kz|')

tol = 1e-1
n = 450
labelx = '$\omega/c$'
labely = 'Im($\epsilon_1$)'
cond_init = fcond_inicial(hbaramu)

list_omegac0 = np.linspace(omegac - tol,omegac + tol,n)
list_im_epsi10 = np.linspace(im_epsi1c - tol,im_epsi1c + tol,n)
X, Y = np.meshgrid(list_omegac0, list_im_epsi10)

def det_2variables2(x0,y0):    
    epsi1 = re_epsi1 + 1j*y0 
    rta = determinante(x0,epsi1,modo,R,hbaramu)
    return np.log10(np.abs(rta))

f2 = np.vectorize(det_2variables2)
Z = f2(X, Y)
plt.figure(figsize=tamfig)
limits = [min(list_omegac0) , max(list_omegac0), min(list_im_epsi10) , max(list_im_epsi10)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.plot([omegac],[im_epsi1c],'.m',ms = 10,label = 'Nelder-Mead')
plt.plot([cond_init[0].real],[cond_init[1].real],'.b',ms = 10,label = 'QE aprox')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.title(title,fontsize=int(tamtitle*0.8)) 
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label('log|det sin kz|',size=tamlegend)

if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('log_det_sinkz_modo%i'%(modo))

def det_2variables3(x0,y0):    
    epsi1 = re_epsi1 + 1j*y0 
    rta = determinante(x0,epsi1,modo,R,hbaramu)
    return np.abs(rta)

f2 = np.vectorize(det_2variables3)
Z = f2(X, Y)
plt.figure(figsize=tamfig)
limits = [min(list_omegac0) , max(list_omegac0), min(list_im_epsi10) , max(list_im_epsi10)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.plot([omegac],[im_epsi1c],'.m',ms = 10,label = 'Nelder-Mead')
plt.plot([cond_init[0]],[cond_init[1]],'.b',ms = 10,label = 'QE aprox')
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.title(title,fontsize=int(tamtitle*0.8)) 
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label('|det sin kz|',size=tamlegend)

del Z
del X
del Y

if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('det_sinkz_modo%i'%(modo))
