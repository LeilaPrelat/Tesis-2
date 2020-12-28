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
path_det = path_basic.replace('/' + 'real_freq','')
path_ctes = path_basic.replace('/sin_kz/dispersive/real_freq','') 


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
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

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

R = 0.01              #micrones
hbaramu = 0.3
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

list_Ep =  np.linspace(0,0.9,9001)  

list_modos = [1,2,3,4]

#%%

# if R != 0.5:
#     raise TypeError('Wrong value for radium')

if list_Ep[0] != 0:
    raise TypeError('El barrido en Ep tiene que empezar en 0 (facil convergencia)') 

if hbaramu != 0.3:
    raise TypeError('Wrong value for mu_c') 
    #estos valores se van a usar para los barridos en mu que parten de mu = 0.3

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/R_%.2f/epsiinf_DL_%.2f_vs_Ep' %(R,epsiinf_DL)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)
        
#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(Ep_var,modo):
    """
    Parameters
    ----------
    Ep_var : hbar*omega_p siendo omega_p el de la permeabilidad electrica de DL

    Returns
    -------
    [re(omega/c), im(epsilon1)]
    
    Uso como condiciones iniciales las funciones de QE_lossless.py (ver seccion 1.7 del cuaderno corto)
    """
    a = omegac_cuasi(modo,Ep_var,epsiinf_DL,gamma_DL,R,hbaramu)
    b = im_epsi1_cuasi(a,Ep_var,epsiinf_DL,gamma_DL,modo,R,hbaramu) 
    return [a.real,b.real]    

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

Ep0 = list_Ep[0]

for modo in list_modos:

    info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.1f eV, modo = %i' %(R,hbaramu,modo)
    info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
    title = info1 +'\n' + info2  + ', ' + name_this_py
    info = info1 + ', ' + info2    

    cond_inicial = fcond_inicial(Ep0,modo)
    
    epsi1_imag_opt = []
    omegac_opt = []
    eq_det = []
    list_Ep_opt = []
    
    tol_NM = 1e-13
    ite_NM = 1150
    for Ep in list_Ep:
        Ep = np.round(Ep,4)
        print('')
        print(Ep)
           
        def det_2variables(x):
            [omegac,epsi_ci] = x
            rta = determinante(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,hbaramu)
            return np.abs(rta)
            
        res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                       options={'maxiter':ite_NM})
    #        print(res.message)
        # if res.message == 'Optimization terminated successfully.' :
        if res.x[1] <= fcond_inicial(Ep,modo)[1]: #QE tiene que requerir menor medio activo que la sol numerica
            omegac_opt.append(res.x[0])
            epsi1_imag_opt.append(res.x[1])
            eq_det.append(det_2variables([res.x[0],res.x[1]]))
            list_Ep_opt.append(Ep)
            
            cond_inicial = [res.x[0],res.x[1]]
        
            
    if save_data_opt==1:
        os.chdir(path)
        print('Guardar data de minimizacion en .txt')
    
        tabla = np.array([list_Ep_opt,omegac_opt,epsi1_imag_opt,eq_det])
        tabla = np.transpose(tabla)
        info2 = '.Opt det SIN kz nano,' + info
        header1 = 'Ep [eV]     Omega/c [1/micrones]    Im(epsi1)     Eq(det)' + info2 + ', ' + name_this_py
        np.savetxt('opt_det_sinkz_vs_Ep_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)
    
    
    label_graph = 'Opt det sin kz'
    label_QE = 'QE approx sin perdidas'
    labelx = '$E_p$ [eV]'
    labely = r'$\epsilon_{ci}$'
    
    im_epsi1_QE = []
    omegac_QE = []
    for Ep in list_Ep:
        a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu)
        b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu) 
        im_epsi1_QE.append(b)
        omegac_QE.append(a)
    
    plt.figure(figsize=tamfig)
    plt.plot(list_Ep,im_epsi1_QE,'.m',ms=10,label=label_QE)
    plt.plot(list_Ep_opt,epsi1_imag_opt,'.-r',ms=10,label=label_graph)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(labely,fontsize=int(tamletra*1.2))
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        os.chdir(path)
        plt.savefig('Im_epsi1_vs_Ep_%i'%(modo))
    
    plt.figure(figsize=tamfig)
    plt.plot(list_Ep,omegac_QE,'.m',ms=10,label=label_QE)
    plt.plot(list_Ep_opt,omegac_opt,'.-r',ms=10,label=label_graph)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        os.chdir(path)
        plt.savefig('Omegac_vs_Ep_%i'%(modo))
    
    plt.figure(figsize=tamfig)
    plt.plot(list_Ep_opt,eq_det,'.-r',ms=10,label=label_graph)
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'|det sin kz|',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        os.chdir(path)
        plt.savefig('detsinkz_vs_Ep_%i'%(modo))
        
    list_omegap = (np.array(list_Ep)/hb)*1e-12
    list_omegap_opt = (np.array(list_Ep_opt)/hb)*1e-12
    omega_QE = (np.array(omegac_QE)*c)*1e-12
    omega_opt = (np.array(omegac_opt)*c)*1e-12
    
    plt.figure(figsize=tamfig)
    plt.plot(list_omegap,omega_QE,'.m',ms=10,label=label_QE)
    plt.plot(list_omegap_opt,omega_opt,'.-r',ms=10,label=label_graph)
    # plt.plot(omega_opt,omega_opt,'-k',label ='y = x')
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(r'$\omega$ [THz]',fontsize=tamletra)
    plt.xlabel(r'$\omega_p$ [THz]',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
    if save_graphs==1:
        os.chdir(path)
        plt.savefig('Omega_vs_omegap_%i'%(modo))

#%%
