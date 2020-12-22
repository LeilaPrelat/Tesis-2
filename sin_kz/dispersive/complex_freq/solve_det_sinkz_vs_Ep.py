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
path_det = path_basic.replace('/' + 'complex_freq','')
path_ctes = path_basic.replace('/sin_kz/dispersive/complex_freq','') 


try:
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante
except ModuleNotFoundError:
    print('det_sinkz_nano.py no se encuentra en ' + path_det)
    path_basic2 = input('path de la carpeta donde se encuentra det_sinkz_nano.py')
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante

#Para condiciones iniciales
try:
    sys.path.insert(1, path_det)
    from complex_omegac_QE import omegac_QE
except ModuleNotFoundError:
    print('complex_omegac_QE.py no se encuentra en ' + path_det)
    path_det = input('path de la carpeta donde se encuentra complex_omegac_QE.py')
    sys.path.insert(1, path_det)
    from complex_omegac_QE import omegac_QE

try:
    sys.path.insert(1, path_ctes)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_ctes)
    path_ctes = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_ctes)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 1

hbaramu = 0.3
epsi_ci = 0
epsiinf_DL = 4.9
gamma_DL = 0.01 #unidades de energia

list_Ep =  np.linspace(0,0.9,901)  

info1 = 'R = %.1f $\mu$m, $\mu_c$ = %.1f eV, modo = %i, $\epsilon_{ci}$ = %i' %(R,hbaramu,modo,epsi_ci)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
title = info1 +'\n' + info2  + ', ' + name_this_py
info = info1 + ', ' + info2 + ', ' + name_this_py

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

if epsi_ci != 0:
    raise TypeError('Wrong value for epsilon_{ci}: must be zero')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/epsiinf_DL_%.2f_vs_Ep' %(epsiinf_DL)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)
        
#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(Ep_var):
    """
    Parameters
    ----------
    Ep_var : hbar*omega_p siendo omega_p el de la permeabilidad electrica de DL

    Returns
    -------
    [re(omega/c), im(omega/c),]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_QE(modo,Ep_var,epsiinf_DL,gamma_DL,epsi_ci,R,hbaramu)
    return [a.real,a.imag]    

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

Ep0 = list_Ep[0]
cond_inicial = fcond_inicial(Ep0)

re_omegac_opt = []
im_omegac_opt = []
eq_det = []
list_Ep_opt = []

tol_NM = 1e-13
ite_NM = 1150
for Ep in list_Ep:
    Ep = np.round(Ep,4)
    print('')
    print(Ep)
       
    def det_2variables(x):
        [re_omegac,im_omegac] = x
        omegac = re_omegac + 1j*im_omegac
        rta = determinante(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,hbaramu)
        return np.abs(rta)
        
    res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
    value = det_2variables([res.x[0],res.x[1]])
    # if res.message == 'Optimization terminated successfully.' :
    if value <= 1e-10:
    # if res.x[1] <= fcond_inicial(Ep)[1]: #QE tiene que requerir menor medio activo que la sol numerica
        re_omegac_opt.append(res.x[0])
        im_omegac_opt.append(res.x[1])
        eq_det.append(value)
        list_Ep_opt.append(Ep)
        
        cond_inicial = [res.x[0],res.x[1]]
        
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([list_Ep_opt,re_omegac_opt,im_omegac_opt,eq_det])
    tabla = np.transpose(tabla)
    header1 = 'Ep [eV]     Re(omega/c) [1/micrones]    Im(omega/c) [1/micrones]     Eq(det)' + info
    np.savetxt('opt_det_nano_sinkz_vs_Ep_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)

#%%

label_graph = 'Opt det sin kz'
label_QE = 'QE approx sin perdidas'
labelx = '$E_p$ [eV]'

re_omegac_QE = []
im_omegac_QE = []
for Ep in list_Ep:
    a = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,epsi_ci,R,hbaramu)

    re_omegac_QE.append(a.real)
    im_omegac_QE.append(a.imag)
    

plt.figure(figsize=tamfig)
plt.plot(list_Ep,re_omegac_QE,'.m',ms=10,label=label_QE)
plt.plot(list_Ep_opt,re_omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Re($\omega/c$) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path)
    plt.savefig('re_omegac_vs_Ep_%i'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_Ep,im_omegac_QE,'.m',ms=10,label=label_QE)
plt.plot(list_Ep_opt,im_omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\omega/c$) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('im_omegac_vs_Ep_%i'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_Ep_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'|det sin kz|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('detsinkz_vs_Ep_%i'%(modo))

#%%
