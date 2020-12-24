#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
    
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
name_this_py = 'Ver ' + name_this_py

try:
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante
except ModuleNotFoundError:
    print('det_sinkz_nano.py no se encuentra en el path_det definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra det_sinkz_nano.py')
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en el path_basic definido/carpeta de trabajo')
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
modo = 3
R = 0.5         #micrones 
#2) Mezcla: DL
epsiinf_DL = 4.9
gamma_DL = 0.01 #unidades de omega
Ep = 0.3

if epsiinf_DL == 4.9:
    Ep = 0.3

#%%
    
if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Importar los valores de SPASER')

try:
    path_load = path_basic + '/R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
    os.chdir(path_load)
except OSError or IOError:
    path_load = path_basic + '/R_%.2f/epsiinf_DL_%.2f_vs_mu' %(R,epsiinf_DL)
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

m = len(barrido_mu)
index = -1
index = int(index)
hbaramu = barrido_mu[index]
delta_ci = epsi1_imag_opt[index]
omegac = omegac_opt[index] 

if np.abs(index) > 0: #index != 0 
    cond_init2 = [omegac_opt[index-1],epsi1_imag_opt[index-1]]

info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.3f eV, $\epsilon_\infty$ = %.1f, Ep = %.1f eV, $\gamma_{DL}$ = %.2f eV' %(R,hbaramu,epsiinf_DL,Ep,gamma_DL)
info2 = '$\Delta_{ci}$ = %.3e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(delta_ci,omegac,modo)
title = info1 + '\n' + info2 + '\n' + name_this_py

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
    a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu)
    b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu) 
    return [a,b]  

#%%

print('Calcular el |det sin kz|')

cond_init = fcond_inicial(hbaramu)

tol1 = np.abs(omegac-cond_init[0].real)
tol2 = np.abs(delta_ci-cond_init[1].real)
tol = np.max([tol1,tol2])
tol = tol*1.2
tol = np.round(tol,3)

n = 250
labelx = '$\omega/c$'
labely = 'Im($\epsilon_1$)'

labelQE = 'QE aprox'
labelci = '$\mu$ anterior'

list_omegac0 = np.linspace(omegac - tol,omegac + tol,n)
list_im_epsi10 = np.linspace(delta_ci - tol,delta_ci + tol,n)
X, Y = np.meshgrid(list_omegac0, list_im_epsi10)

def det_2variables2(x0,y0):    
    rta = determinante(x0,Ep,epsiinf_DL,gamma_DL,y0,modo,R,hbaramu)
    return np.log10(np.abs(rta))

f2 = np.vectorize(det_2variables2)
Z = f2(X, Y)
plt.figure(figsize=tamfig)
limits = [min(list_omegac0) , max(list_omegac0), min(list_im_epsi10) , max(list_im_epsi10)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.plot([omegac],[delta_ci],'.m',ms = 10,label = 'Nelder-Mead')
plt.plot([cond_init[0].real],[cond_init[1].real],'.b',ms = 10,label = labelQE)
if np.abs(index) > 0: #index != 0 
    plt.plot([cond_init2[0].real],[cond_init2[1].real],'.g',ms = 10,label = labelci)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.title(title,fontsize=int(tamtitle*0.8)) 
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label('log|det sin kz|',size=tamlegend)

if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('log_det_sinkz_mu%.3f_%i'%(hbaramu,modo), format='png')

def det_2variables3(x0,y0):    
    rta = determinante(x0,Ep,epsiinf_DL,gamma_DL,y0,modo,R,hbaramu)
    return np.abs(rta)

f2 = np.vectorize(det_2variables3)
Z = f2(X, Y)
plt.figure(figsize=tamfig)
limits = [min(list_omegac0) , max(list_omegac0), min(list_im_epsi10) , max(list_im_epsi10)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.plot([omegac],[delta_ci],'.m',ms = 10,label = 'Nelder-Mead')
plt.plot([cond_init[0].real],[cond_init[1].real],'.b',ms = 10,label = labelQE)
if np.abs(index) > 0: #index != 0 
    plt.plot([cond_init2[0].real],[cond_init2[1].real],'.g',ms = 10,label = labelci)
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
    plt.savefig('det_sinkz_mu%.3f_%i'%(hbaramu,modo), format='png')

#%%