#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)

determinante de 4x4

barrido en mu para modelando el epsilon(omega) como una lorentziana.  
Hallar omega/c real y re(kz) que minimizan
el determinante con kz para diferentes valores de mu, 
el valor de im(epsi1) esta fijo por la funcion lorentziana de epsilon(omega)

"""

import numpy as np
import os
import sys
from scipy.optimize import minimize   
import matplotlib.pyplot as plt
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
path_basic2 = path_basic.replace('/' + 'real_freq_lorentziana','')
path_graphene = path_basic2.replace('/' + 'con_kz_real','') 
path_sinkz = path_basic2.replace('/' + 'con_kz_real','') + '/' + 'sin_kz' + '/' + 'non-dispersive/real_freq'

try:
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante
except ModuleNotFoundError:
    print('det_conkz.py no se encuentra en ' + path_basic2)
    path_basic2 = input('path de la carpeta donde se encuentra det_conkz.py')
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante

try:
    sys.path.insert(1, path_basic2 + '/' + 'extra')
    from def_kt import kt
except ModuleNotFoundError:
    print('def_kt.py no se encuentra en ' + path_basic2 + '/' + 'extra')
    path_basic2 = input('path de la carpeta donde se encuentra def_kt.py')
    sys.path.insert(1, path_basic2)
    from def_kt import kt

#freq en THz
try:
    sys.path.insert(1, path_basic2)
    from epsilon_lorentziana import epsilon 
except ModuleNotFoundError:
    print('epsilon_lorentziana.py no se encuentra en ' + path_basic2)
    path_basic2 = input('path de la carpeta donde se encuentra epsilon_lorentziana.py')
    sys.path.insert(1, path_basic2)
    from epsilon_lorentziana import epsilon

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

re_epsi1 = 16
R = 0.25              #micrones
   
modo = 2
eta = 0.9
aux_cte = c*1e-12 # freq en THz para epsilon
   
list_mu = np.linspace(0.3,0.9,601) 
kzlim = 0.132

#%%

if R not in [0.5, 0.25]:
    raise TypeError('Mal el valor de R: El barrido en re(epsi1) se hizo para R = 0.5 o 0.25, y se usa como cond inicial')

if list_mu[0]!= 0.3:
    raise TypeError('El barrido en re(epsi1) se hizo para mu = 0.3, para usarlo como cond inicial list_mu debe empezar en 0.3')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/re_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)
    path = path_basic + path_det

    if not os.path.exists(path):
        print('Creating folder to save data')
        os.mkdir(path)

#%%

#print('Definir las condiciones iniciales para el metodo de minimizacion: usar las funciones de QE_lossless.py')

def fcond_inicial(R,mu0):
    """
    Parameters
    ----------
    R (unidad micrones)
    mu0 : potencial quimico (unidad eV)

    Returns
    -------
    [re(omega/c), kz=0]
    
    Uso como condiciones iniciales el barrido en kz hecho
    para R = 0.50 $\mu$m o R = 0.25 $\mu$m, $\mu_c$ = 0.3 eV o 0.6 eV
    """
  
    os.chdir(path_sinkz + '/' + 'R_%.2f_vs_re_epsi1/mu_%.1f' %(R,mu0))
    cond_init = np.loadtxt('opt_det_sinkz_vs_re_epsi1_modo%i.txt' %(modo),delimiter='\t', skiprows = 1)
    cond_init = np.transpose(cond_init)
    [list_re_epsi1, omegac_opt, epsi1_imag_opt, eq_det] = cond_init
    f1 = interp1d(list_re_epsi1,omegac_opt)
    f2 = interp1d(list_re_epsi1,epsi1_imag_opt)
    a = float(f1(re_epsi1))
    b = 0
    return [a,b] 

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

labeltxt = '_vs_mu_modo%i' %(modo)
mu0 = list_mu[0]
cond_inicial = fcond_inicial(R,mu0)

kz_real_opt = []
omegac_opt = []

eq_det = []
list_mu_opt = []

tol_NM = 1e-13
ite_NM = 1150

for mu in list_mu:
    mu = np.round(mu,4)
    print('')
    print(mu)

           
    def det_2variables(x):
        [omegac,kz_real] = x
        freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon
        epsi1 = epsilon(freq, eta)
        rta = determinante(kz_real,omegac,epsi1,modo,R,mu)
        return np.abs(rta)
        
    res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
    if res.message == 'Optimization terminated successfully.':
        omegac_opt.append(res.x[0])
        kz_real_opt.append(res.x[1])
        eq_det.append(det_2variables([res.x[0],res.x[1]]))
        list_mu_opt.append(mu)
        
    cond_inicial = [res.x[0],res.x[1]]
    
        
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([list_mu_opt,omegac_opt,kz_real_opt,eq_det])
    tabla = np.transpose(tabla)
    info = '.Opt det con kz, R=%.2f \mum, Re(epsi1)=%.2f, $\eta$ = %.2f' %(R,re_epsi1,eta) 
    header1 = 'mu [eV]     Omega/c [1/micrones]    kz [1/micrones]    Eq(det)' + info + ', ' + name_this_py
    np.savetxt('opt_det_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

label_graph = 'Opt det con $\eta$ = %.2f' %(eta)
labelsinkz = 'Opt sin kz'
labelx = '$\mu_c$ [eV]'
labelpng = '_vs_mu_%i' %(modo)
title = 'Modo = %i, R = %.2f$\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) 
title2 = title + '\n' + name_this_py

kt_imag = [] # Imag(kt)
kt_real = [] # Real(kt)
kt_abs = [] #modulo del k transversal |kt|

for j in range(len(list_mu_opt)):
    hbaramu = list_mu_opt[j]
    omegac = omegac_opt[j]
    kz_real = kz_real_opt[j]
    
    freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon

    epsi1 = epsilon(freq, eta)
    kt_value = kt(kz_real,omegac,epsi1,modo,R,hbaramu)    
    kt_abs.append(np.abs(kt_value))
    kt_imag.append(kt_value.imag)
    kt_real.append(kt_value.real)
    

plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kz_real_opt,'.-r',ms=10, label = label_graph)
# if kz_real<kzlim:
#     plt.plot(barrido_mu_sinkz,epsi1_imag_opt_sinkz,'.-m',ms=10, label = labelsinkz)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'kz [1/$mu$]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    os.chdir(path)
    plt.savefig('kz' + labelpng, format='png')

plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omegac' + labelpng, format='png')
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'|det|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('det' + labelpng, format='png')
 
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_abs,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
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
    plt.savefig('|kt|' + labelpng, format='png')
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_imag,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('kt_imag' + labelpng, format='png')
    
#%%
