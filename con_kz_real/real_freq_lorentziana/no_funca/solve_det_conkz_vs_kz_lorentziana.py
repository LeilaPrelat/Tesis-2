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
Hallar omega/c real y mu que minimizan
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

R = 0.25              #micrones

modo = 1
eta = 0.2
aux_cte = c*1e-12 # freq en THz para epsilon
   
list_mu = np.linspace(0.3,0.9,601) 

#%%

if R not in [0.5, 0.25]:
    raise TypeError('Mal el valor de R: El barrido en re(epsi1) se hizo para R = 0.5 o 0.25, y se usa como cond inicial')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_data_opt==1:

    path_det = r'/eta_%.2f_vs_mu/R_%.2f' %(eta,R)
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
    R (unidad 1/micrones)
    mu : potencial quimico (unidades eV)

    Returns
    -------
    [re(omega/c)]
    
    Uso como condiciones iniciales el barrido en kz hecho
    para R = 0.50 $\mu$m o R = 0.25 $\mu$m, $\mu_c$ = 0.3 eV o 0.6 eV
    """
    os.chdir(path_sinkz + '/' + 'R_%.2f_vs_re_epsi1/mu_%.1f' %(R,mu0))
    cond_init = np.loadtxt('opt_det_sinkz_vs_re_epsi1_modo%i.txt' %(modo),delimiter='\t', skiprows = 1)
    cond_init = np.transpose(cond_init)
    [list_re_epsi1, omegac_opt, epsi1_imag_opt, eq_det] = cond_init
    f1 = interp1d(list_re_epsi1,omegac_opt)
    # f2 = interp1d(list_re_epsi1,epsi1_imag_opt)
    re_epsi1 = 15
    a = float(f1(re_epsi1))
    return a 

#%%        
"""
print('Minimizacion del determinante de 4x4 para un barrido en kz')

labeltxt = '_vs_kz_modo%i' %(modo)
cond_inicial = [fcond_inicial(0.25,0.3), 0.25, 0.5]

mu_opt = []

omegac_opt = []
R_opt = []
kz_opt = []

eq_det = []

omegac_min, omegac_max = 10/aux_cte, fcond_inicial(0.25,0.3) 
R_min, R_max = 0.1, 0.7
#mu_max, mu_min = 0.9, 0.3
kz_min, kz_max = 0, 0.5 

#if modo == 1:
#    kz_max = 0.13
    
bnds = ((omegac_min, omegac_max), (R_min, R_max), (kz_min, kz_max))

tol_NM = 1e-6
ite_NM = 1150

for mu in list_mu:
    mu = np.round(mu,4)
    print('')
    print(mu)

    def det_3variables(x):
        [omegac, R, kz_real] = x
        freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon
        epsi1 = epsilon(freq, eta)
        
        rta = determinante(kz_real,omegac,epsi1,modo,R,mu)
        return np.abs(rta)
        
    res = minimize(det_3variables, cond_inicial, method='TNC', bounds=bnds, tol=tol_NM, 
                   options={'maxiter':ite_NM})

    if res.success == True:
#        res.x = float(res.x)
        value = det_3variables([res.x[0],res.x[1],res.x[2]])

        omegac_opt.append(res.x[0])
        R_opt.append(res.x[1])
        kz_opt.append(res.x[2])
        
        eq_det.append(value)
        
        mu_opt.append(mu)
        
        cond_inicial = [res.x[0],res.x[1],res.x[2]] 

"""
#%%
"""
from scipy.optimize import minimize_scalar
cond_inicial = fcond_inicial(0.25,0.3)

R = 0.25
list_kz = [0,0.05,0.1,0.15,0.2]
mu_opt_tot = []
omegaTHz_opt_tot = []
cond_init_tot = []

eq_det = []

omegac_min, omegac_max = 10/aux_cte, fcond_inicial(0.25,0.3) 
bnds = (omegac_min, omegac_max)

tol_NM = 1e-6
ite_NM = 1150

j = 1
for kz in list_kz:
 
    mu_opt = []
    omegaTHz_opt = []
    
    for mu in list_mu:
        mu = np.round(mu,4)
        print('')
        print(mu)
    
        def det_1variable(x):
            omegac = x
            freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon
            epsi1 = epsilon(freq, eta)
            
            rta = determinante(kz,omegac,epsi1,modo,R,mu)
            return np.abs(rta)
            
        res = minimize_scalar(det_1variable, cond_inicial, tol=tol_NM)
    
        if res.success == True:
            res.x = float(res.x)
            value = det_1variable(res.x)
    
            omegaTHz_opt.append((res.x)*aux_cte)

            eq_det.append(value)
            mu_opt.append(mu)
            
            cond_inicial = res.x
        if mu == 0.3:
            cond_init_tot.append(res.x)
    
    cond_inicial = cond_init_tot[j-1]
#    cond_inicial = fcond_inicial(0.25,0.3)
    mu_opt_tot.append(mu_opt)
    omegaTHz_opt_tot.append(omegaTHz_opt)

#%%
plt.figure(figsize=tamfig)

for j in range(len(list_kz)):
    kz = list_kz[j]
    omegaTHz_opt = omegaTHz_opt_tot[j]
    mu_opt = mu_opt_tot[j]

    plt.plot(mu_opt,omegac_opt,'.-',ms=10,label='kz = %.2f 1/$\mu$m'%(kz))
    plt.ylabel(r'$\omega$ [THz]',fontsize=tamletra)
    plt.xlabel(r'$\mu$ [eV]',fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    plt.grid(1)
"""
#%%

from scipy.optimize import minimize_scalar
cond_inicial = [5/aux_cte,0.05]

bounded = 1

omegac_min, omegac_max = 1e-3/aux_cte, 90/aux_cte
bnds = ((omegac_min, omegac_max), (0,4))

tol_NM = 1e-6
ite_NM = 1150

mu_opt = []
eq_det = []
kz_opt = []

im_epsi1_opt = []
re_epsi1_opt = []
omegaTHz_opt = []

for mu in list_mu:
    mu = np.round(mu,4)
    print('')
    print(mu)

    def det_2variables(x):
        [omegac, kz] = x
        freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon
        epsi1 = epsilon(freq, eta)
        
        rta = determinante(kz,omegac,epsi1,modo,R,mu)
        return np.abs(rta)

    if bounded == 0:        
        res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                       options={'maxiter':ite_NM})
    
        if res.message == 'Optimization terminated successfully.':
            # res.x = float(res.x)
            value = det_2variables([res.x[0],res.x[1]])
            kz_opt.append(res.x[1])
            omegaTHz_opt.append((res.x[0])*aux_cte)
    
            freq = (res.x[0])*aux_cte # freq tiene que estar en THz para la funcion epsilon
            epsi1 = epsilon(freq, eta)
            im_epsi1_opt.append(epsi1.imag)
            re_epsi1_opt.append(epsi1.real)
    
            eq_det.append(value)
            mu_opt.append(mu)
            
            cond_inicial = res.x
    else:
        res = minimize(det_2variables, cond_inicial, method='TNC', bounds=bnds, tol=tol_NM, 
                   options={'maxiter':ite_NM})
    
        if res.success == True:
            # res.x = float(res.x)
            value = det_2variables([res.x[0],res.x[1]])
            kz_opt.append(res.x[1])
            omegaTHz_opt.append((res.x[0])*aux_cte)
    
            freq = (res.x[0])*aux_cte # freq tiene que estar en THz para la funcion epsilon
            epsi1 = epsilon(freq, eta)
            im_epsi1_opt.append(epsi1.imag)
            re_epsi1_opt.append(epsi1.real)
    
            eq_det.append(value)
            mu_opt.append(mu)
            
            cond_inicial = res.x        

#%%    
labeltxt = '_vs_mu_modo%i' %(modo)
    
if save_data_opt==1:
    os.chdir(path)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([mu_opt,omegaTHz_opt,kz_opt,eq_det])
    tabla = np.transpose(tabla)
    info = '.Opt det con kz, R=%.2f \mum, $\eta$ = %.2f' %(R,eta) 
    header1 = 'mu [eV]     Omega [THz]    kz [1/mu]    Eq(det)' + info + ', ' + name_this_py
    np.savetxt('opt_det_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

label_graph = 'Opt det con $\eta$ = %.2f' %(eta)
labelsinkz = 'Opt sin kz'
labelx = 'kz [$1/mu$]'
labelpng = '_vs_mu_%i.png' %(modo)
title = 'Modo = %i, R = %.2f$\mu$m, $\eta$ = %.2f' %(modo,R,eta) 
title2 = title + '\n' + name_this_py

plt.figure(figsize=tamfig)
plt.plot(mu_opt,omegaTHz_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'$\omega$ [THz]',fontsize=tamletra)
plt.xlabel(r'$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('OmegaTHz' + labelpng, format='png')


plt.figure(figsize=tamfig)
plt.plot(mu_opt,im_epsi1_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'Im($\epsilon_1$)',fontsize=tamletra)
plt.xlabel(r'$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('im_epsi1' + labelpng, format='png')
    
plt.figure(figsize=tamfig)
plt.plot(mu_opt,re_epsi1_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'Re($\epsilon_1$)',fontsize=tamletra)
plt.xlabel(r'$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('re_epsi1' + labelpng, format='png')

plt.figure(figsize=tamfig)
plt.plot(mu_opt,kz_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'kz [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(r'$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('kz' + labelpng, format='png')
    
plt.figure(figsize=tamfig)
plt.plot(mu_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'|det|',fontsize=tamletra)
plt.xlabel(r'$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('det' + labelpng, format='png')
       

#%%
"""
label_graph = 'Opt det con $\eta$ = %.2f' %(eta)
labelsinkz = 'Opt sin kz'
labelx = 'kz [$1/mu$]'
labelpng = '_vs_mu_%i' %(modo)
title = 'Modo = %i, R = %.2f$\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) 
title2 = title + '\n' + name_this_py

kt_imag = [] # Imag(kt)
kt_real = [] # Real(kt)
kt_abs = [] #modulo del k transversal |kt|

hbaramu = mu0
for j in range(len(list_kz_opt)):

    omegac = omegac_opt[j]
    kz_real = list_kz_opt[j]
    
    freq = omegac*aux_cte # freq tiene que estar en THz para la funcion epsilon

    epsi1 = epsilon(freq, eta)
    kt_value = kt(kz_real,omegac,epsi1,modo,R,hbaramu)    
    kt_abs.append(np.abs(kt_value))
    kt_imag.append(kt_value.imag)
    kt_real.append(kt_value.real)
    

# plt.figure(figsize=tamfig)
# plt.plot(list_kz_opt,mu_opt,'.-r',ms=10, label = label_graph)
# # if kz_real<kzlim:
# #     plt.plot(barrido_mu_sinkz,epsi1_imag_opt_sinkz,'.-m',ms=10, label = labelsinkz)
# plt.title(title2,fontsize=tamtitle)
# plt.ylabel(r'$\mu_c$ [eV]',fontsize=tamletra)
# plt.xlabel(labelx,fontsize=tamletra)
# plt.tick_params(labelsize = tamnum)
# plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
# plt.grid(1)
# if save_graphs==1:
#     os.chdir(path)
#     plt.savefig('mu' + labelpng, format='png')

plt.figure(figsize=tamfig)
plt.plot(list_kz_opt,omegac_opt,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'$\omega/c$ [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('Omegac' + labelpng, format='png')
    
plt.figure(figsize=tamfig)
plt.plot(list_kz_opt,eq_det,'.-r',ms=10,label=label_graph)
plt.title(title2,fontsize=tamtitle)
plt.ylabel(r'|det|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)
if save_graphs==1:
    plt.savefig('det' + labelpng, format='png')
 """   
#%%
