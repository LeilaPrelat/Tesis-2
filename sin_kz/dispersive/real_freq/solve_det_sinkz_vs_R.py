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
from scipy.interpolate import interp1d

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

try:
    sys.path.insert(1, path_ctes)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_ctes = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_ctes)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

modo = 2
hbaramu = 0.3

Ep = 0.3
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

list_R =  np.linspace(0.5,5.5,5001)
# list_mu = [0.3]

if modo == 4:
    list_R = np.linspace(0.5,5,4501)

info1 = '$\mu$ = %.1f eV, $E_p$ = %.3f eV, modo = %i' %(hbaramu,Ep,modo)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
title = info1 +'\n' + info2  + ', ' + name_this_py
info = info1 + ', ' + info2

#%%

if list_R[0] != 0.5:
    raise TypeError('Las condiciones iniciales (barrido en Ep) fueron hechas para R = 0.5')

if gamma_DL != 0.01:
    raise TypeError('Las condiciones iniciales (barrido en Ep) fueron hechas para gamma_DL = 0.01')
    
if hbaramu != 0.3:
    raise TypeError('Las condiciones iniciales (barrido en Ep) fueron hechas para mu = 0.3')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_save = r'/epsiinf_DL_%.2f_vs_R' %(epsiinf_DL)
    path = path_basic + path_save

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)

path = path + '/' + 'Ep_%.2f' %(Ep)
        
#%%

print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(Ep):
    """
    Parameters
    ----------
    Ep : hbar*omega_p (unidades eV)

    Returns
    -------
    [re(omega/c), im(epsilon1)]
    
    Uso como condiciones iniciales el barrido en Ep hecho
    para R = 0.5 $\mu$m, $\mu_c$ = 0.3 eV, $\gamma_{DL}$ = 0.01 eV
    """
    os.chdir(path_basic + '/' + 'epsiinf_DL_%.2f_vs_Ep' %(epsiinf_DL))
    cond_init = np.loadtxt('opt_det_sinkz_vs_Ep_modo%i.txt' %(modo),delimiter='\t', skiprows = 1)
    cond_init = np.transpose(cond_init)
    [Ep_opt,omegac_opt,epsi1_imag_opt,eq_det] = cond_init
    f1 = interp1d(Ep_opt,omegac_opt)
    f2 = interp1d(Ep_opt,epsi1_imag_opt)
    a = float(f1(Ep))
    b = float(f2(Ep))
    return [a,b]      

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

R0 = list_R[0]
cond_inicial = fcond_inicial(Ep)

epsi1_imag_opt = []
omegac_opt = []
eq_det = []
list_R_opt = []

tol_NM = 1e-13
ite_NM = 1150

for R in list_R:
    R = np.round(R,3)
    print('')
    print(R)

    def det_2variables(x):
        [omegac,epsi_ci] = x
        rta = determinante(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,hbaramu)
        return np.abs(rta)
        
    res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
    # if res.message == 'Optimization terminated successfully.':
    a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu)
    value = det_2variables([res.x[0],res.x[1] ])
    if res.x[1] <=  im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu) and value < tol_NM: #QE tiene que requerir menor medio activo que la sol numerica

        omegac_opt.append(res.x[0])
        epsi1_imag_opt.append(res.x[1])
        eq_det.append(value)
        list_R_opt.append(R)
        
        cond_inicial = [res.x[0],res.x[1]]
    
        
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([list_R_opt,omegac_opt,epsi1_imag_opt,eq_det])
    tabla = np.transpose(tabla)
    info2 = '.Opt det SIN kz nano,' + info
    header1 = 'R [micrones]     Omega/c [1/micrones]    Im(epsi1)     Eq(det)' + info + ', ' + name_this_py
    np.savetxt('opt_det_sinkz_vs_R_modo%i.txt' %(modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)

#%%

label_graph = 'Opt det sin kz'
label_QE = 'QE approx sin perdidas'
labelx = 'R [$\mu$m]'
labely = r'$\epsilon_{ci}$'

im_epsi1_QE = []
omegac_QE = []
for R in list_R:
    a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu)
    b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu) 
    im_epsi1_QE.append(b)
    omegac_QE.append(a)

fig = plt.figure(figsize=tamfig)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex = ax1)
fig.subplots_adjust(hspace=0.05)
ax1.set_title(title,loc = 'center')
ax1.title.set_size(tamtitle)
ax1.plot(list_R,im_epsi1_QE,'.m',ms=10,label=label_QE)
ax2.plot(list_R_opt,epsi1_imag_opt,'.-r',ms=10,label=label_graph)
plt.setp(ax2.get_xticklabels(), fontsize = tamnum)
plt.setp(ax2.get_yticklabels(), fontsize = tamnum)
plt.setp(ax1.get_xticklabels(), visible = False)
plt.setp(ax1.get_yticklabels(), fontsize = tamnum)

ax1.set_ylabel(labely,fontsize=tamletra)
ax2.set_xlabel(labelx,fontsize=tamletra)
ax2.set_ylabel(labely,fontsize=tamletra)

ax1.legend(loc='best',markerscale=2,fontsize=tamlegend)
ax2.legend(loc='best',markerscale=2,fontsize=tamlegend)
ax1.grid(1)      
ax2.grid(1)  
if save_graphs==1:
    os.chdir(path)
    plt.savefig('Im_epsi1_vs_R_%i'%(modo))


# plt.figure(figsize=tamfig)
# plt.plot(list_R,im_epsi1_QE,'.m',ms=10,label=label_QE)
# plt.plot(list_R_opt,epsi1_imag_opt,'.-r',ms=10,label=label_graph)
# plt.title(title,fontsize=tamtitle)
# plt.ylabel(labely,fontsize=tamletra)
# plt.xlabel(labelx,fontsize=tamletra)
# plt.tick_params(labelsize = tamnum)
# plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
# plt.grid(1)
# if save_graphs==1:
#     os.chdir(path)
#     plt.savefig('Im_epsi1_vs_R_%i'%(modo))
    
plt.figure(figsize=tamfig)
plt.plot(list_R,omegac_QE,'.m',ms=10,label=label_QE)
plt.plot(list_R_opt,omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omegac_vs_R_%i'%(modo))

plt.figure(figsize=tamfig)
plt.plot(list_R_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'|det|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('det_vs_R_%i'%(modo))

omega_QE = (np.array(omegac_QE)*c)*1e-12
omega_opt = (np.array(omegac_opt)*c)*1e-12
omegap = (Ep/hb)*1e-12

plt.figure(figsize=tamfig)
plt.plot(list_R,omega_QE,'.m',ms=10,label=label_QE)
# plt.plot(list_mu, np.ones(len(list_mu))*omegap,'.b',label = '$\omega_p$')
plt.plot(list_R_opt,omega_opt,'.-r',ms=10,label=label_graph)
plt.plot([],[],'.w',label = '$\omega_p$ = %.2f THz' %(omegap))
# plt.plot(omega_opt,omega_opt,'-k',label ='y = x')
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\omega$ [THz]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omega_vs_R_%i'%(modo))

#%%
