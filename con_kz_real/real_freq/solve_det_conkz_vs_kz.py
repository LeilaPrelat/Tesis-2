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
from scipy.optimize import minimize   
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
    from det_conkz import determinante
except ModuleNotFoundError:
    print('det_conkz.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra det_conkz.py')
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante

try:
    sys.path.insert(1, path_basic2)
    from def_kt import kt
except ModuleNotFoundError:
    print('def_kt.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra def_kt.py')
    sys.path.insert(1, path_basic2)
    from def_kt import kt

try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
hbaramu = 0.3        #eV mu_c
modo = 1
    
list_kz = np.linspace(0,4,8001) 
# list_kz = np.linspace(0.105,1.5,1395+1)
# list_kz = [0]

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

if list_kz[0]!= 0:
    raise TypeError('El barrido en kz debe empezar de 0 para usar el QE_approx del caso sin kz')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/re_epsi1_%.2f_vs_kz' %(re_epsi1)
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
    re_epsi1 : parte real de la permeabilidad electrica del medio 1

    Returns
    -------
    [re(omega/c), im(epsilon1)]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_cuasi(modo,R,re_epsi1,hbaramu)
    b = im_epsi1_cuasi(a,modo,R,hbaramu) 
    return [a,b]    

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

cond_inicial = fcond_inicial(hbaramu)

epsi1_imag_opt = []
omegac_opt = []
eq_det = []
list_kz_opt = []

tol_NM = 1e-13
ite_NM = 1150

for kz in list_kz:
    kz = np.round(kz,4)
    print('')
    print(kz)

           
    def det_2variables(x):
        [omegac,im_epsi1] = x
        epsi1 = re_epsi1 + 1j*im_epsi1
        rta = determinante(kz,omegac,epsi1,modo,R,hbaramu)
        return np.abs(rta)
        
    res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
    if res.message == 'Optimization terminated successfully.':
        omegac_opt.append(res.x[0])
        epsi1_imag_opt.append(res.x[1])
        eq_det.append(det_2variables([res.x[0],res.x[1]]))
        list_kz_opt.append(kz)
        
    cond_inicial = [res.x[0],res.x[1]]
    
        
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det])
    tabla = np.transpose(tabla)
    info = '.Opt det con kz, R=%.1f \mum, Re(epsi1)=%.2f, $\mu_c$ = %.3f eV' %(R,re_epsi1,hbaramu) 
    header1 = 'kz [eV]     Omega/c [1/micrones]    Im(epsi1)     Eq(det)' + info + ', ' + name_this_py
    np.savetxt('opt_det_conkz_vs_kz_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)

#%%

label_graph = 'Opt det con kz'
labelx = '$k_z$ [1/$\mu$m]'
labelpng = '_vs_kz_%i' %(modo)
title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f, $\mu_c$ = %.3f eV' %(modo,R,re_epsi1,hbaramu) 
title2 = title + '\n' + name_this_py
a = omegac_cuasi(modo,R,re_epsi1,hbaramu)
b = im_epsi1_cuasi(a,modo,R,hbaramu) 

kt_imag1 = [] # Imag(kt)
kt_real1 = [] # Real(kt)
kt_abs1 = [] #modulo del k transversal |kt|
list_kz_opt1 = []
epsi1_imag_opt1 = []
omegac_opt1 = []

kt_imag2 = [] # Imag(kt)
kt_real2 = [] # Real(kt)
kt_abs2 = [] #modulo del k transversal |kt|
list_kz_opt2 = []
epsi1_imag_opt2 = []
omegac_opt2 = []
for j in range(len(list_kz_opt)):
    kz_var = list_kz_opt[j]
    omegac = omegac_opt[j]
    epsi1_imag = epsi1_imag_opt[j]
    epsi1 = re_epsi1 + 1j*epsi1_imag
    kt_value = kt(kz_var,omegac,epsi1,modo,R,hbaramu)
    if kt_value < 1: 
        kt_abs1.append(np.abs(kt_value))
        kt_real1.append(kt_value.real)
        kt_imag1.append(kt_value.imag)
        list_kz_opt1.append(kz_var)
        epsi1_imag_opt1.append(epsi1_imag)
        omegac_opt1.append(omegac)
    else:
        kt_abs2.append(np.abs(kt_value))
        kt_real2.append(kt_value.real)
        kt_imag2.append(kt_value.imag)
        list_kz_opt2.append(kz_var)
        epsi1_imag_opt2.append(epsi1_imag)
        omegac_opt2.append(omegac)        

hspace = 0.15
fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(hspace = hspace)
ax1.set_title(title2,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_kz_opt1,epsi1_imag_opt1,'.m',ms=10)
# ax1.plot([0],[b],'.b', ms=10, label = 'QE approx (sin kz)')
ax2.plot(list_kz_opt2,epsi1_imag_opt2,'.m',ms=10)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(),fontsize = tamnum)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)
ax1.set_ylabel('Im($\epsilon_1$)',fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel('Im($\epsilon_1$)',fontsize=tamletra)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    os.chdir(path)
    plt.savefig('Im_epsi1' + labelpng)
    
fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(hspace = hspace)
ax1.set_title(title2,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_kz_opt1,omegac_opt1,'.m',ms=10)
# ax1.plot([0],[a],'.b', ms=10, label = 'QE approx (sin kz)')
ax2.plot(list_kz_opt2,omegac_opt2,'.m',ms=10)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(), fontsize = tamnum)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)
ax1.set_ylabel('$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel('$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    plt.savefig('Omegac' + labelpng)

plt.figure(figsize=tamfig)
plt.plot(list_kz_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'|det con kz|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('det' + labelpng)

fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(hspace = hspace)
ax1.set_title(title2,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_kz_opt1,kt_abs1,'.m',ms=10)
ax2.plot(list_kz_opt2,kt_abs2,'.m',ms=10)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(), fontsize = tamnum)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)
ax1.set_ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    plt.savefig('|kt|' + labelpng)
    
fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(hspace = hspace)
ax1.set_title(title2,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_kz_opt1,kt_real1,'.m',ms=10)
ax2.plot(list_kz_opt2,kt_real2,'.m',ms=10)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(), fontsize = tamnum)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)
ax1.set_ylabel('Re($k_t$) [1/$\mu$m]',fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel('Re($k_t$) [1/$\mu$m]',fontsize=tamletra)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    plt.savefig('kt_real' + labelpng)
    
fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
fig.subplots_adjust(hspace = hspace)
ax1.set_title(title2,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_kz_opt1,kt_imag1,'.m',ms=10)
ax2.plot(list_kz_opt2,kt_imag2,'.m',ms=10)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(), fontsize = tamnum)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)
ax1.set_ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    plt.savefig('kt_imag' + labelpng)
    
#%%
