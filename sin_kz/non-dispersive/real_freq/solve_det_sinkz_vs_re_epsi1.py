#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

barrido en mu. Hallar omega/c real y im(epsilon1) que minimizan
el determinante sin kz para diferentes valores de mu
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from scipy.optimize import minimize   

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
    from det_sinkz import determinante
except ModuleNotFoundError:
    print('det_sinkz.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra det_sinkz.py')
    sys.path.insert(1, path_basic2)
    from det_sinkz import determinante

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

#%%

print('Definir parametros del problema')

R = 0.5
hbaramu = 0.9           #eV
modo = 4

list_re_epsi1 =  np.linspace(1,9,801)  
# list_mu = [0.3]

#%%

"""
if hbaramu != 0.3:
    raise TypeError('Wrong value for mu')
"""  
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/R_%.1f_vs_re_epsi1/mu_%.1f' %(R,hbaramu)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)
        
#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(re_epsi1):
    """
    Parameters
    ----------
    re_epsi1 : la parte real de la permeabilidad electrica del medio 1

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

re_epsi1_0 = list_re_epsi1[0]
cond_inicial = fcond_inicial(re_epsi1_0)

epsi1_imag_opt = []
omegac_opt = []
eq_det = []
list_re_epsi1_opt = []

tol_NM = 1e-13
ite_NM = 1150
for re_epsi1 in list_re_epsi1:
    re_epsi1 = np.round(re_epsi1,4)
    print('')
    print(re_epsi1)

           
    def det_2variables(x):
        [omegac,im_epsi1] = x
        epsi1 = re_epsi1 + 1j*im_epsi1
        rta = determinante(omegac,epsi1,modo,R,hbaramu)
        return np.abs(rta)
        
    res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
    if res.message == 'Optimization terminated successfully.':
        omegac_opt.append(res.x[0])
        epsi1_imag_opt.append(res.x[1])
        eq_det.append(det_2variables([res.x[0],res.x[1]]))
        list_re_epsi1_opt.append(re_epsi1)
        
    cond_inicial = [res.x[0],res.x[1]]

        
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([list_re_epsi1_opt,omegac_opt,epsi1_imag_opt,eq_det])
    tabla = np.transpose(tabla)
    info = '.Opt det SIN kz, mu = %.1f eV, R = %.1f micrones' %(hbaramu,R) 
    header1 = 'Re(epsi1)     Omega/c [1/micrones]    Im(epsi1)     Eq(det)' + info + ', ' + name_this_py
    np.savetxt('opt_det_sinkz_vs_re_epsi1_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)

label_graph = 'Opt det sin kz'
label_QE = 'QE approx sin perdidas'
title = 'Modo = %i, $\mu$ = %.1f eV, R = %.1f $\mu$m' %(modo,hbaramu,R) +  ', ' + name_this_py
labelx = 'Re($\epsilon_1$)'

im_epsi1_QE = []
omegac_QE = []
for re_epsi1 in list_re_epsi1:
    a = omegac_cuasi(modo,R,re_epsi1,hbaramu)
    b = im_epsi1_cuasi(a,modo,R,hbaramu) 
    im_epsi1_QE.append(b)
    omegac_QE.append(a)

plt.figure(figsize=tamfig)
plt.plot(list_re_epsi1,im_epsi1_QE,'.m',ms=10,label=label_QE)
plt.plot(list_re_epsi1_opt,epsi1_imag_opt,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\epsilon_1$)',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path)
    plt.savefig('Im_epsi1_vs_re_epsi1_%i'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_re_epsi1,omegac_QE,'.m',ms=10,label=label_QE)
plt.plot(list_re_epsi1_opt,omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path)
    plt.savefig('Omegac_vs_re_epsi1_%i'%(modo))

#%%
