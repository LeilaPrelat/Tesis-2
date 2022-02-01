#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)

determinante de 4x4

graficar el log|det(kz)| en mapa de color en funcion de 
omega en THz y del potencial quimico mu
y minimizar

"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import minimize  
# from scipy.optimize import minimize_scalar

#%% 

avisar_python_termino = 1
save_data_opt = 1 #guardar data de la minimizacion
save_graphs = 1 #guardar los graficos

## 1 variable #########################

minimizar_omega = 0 # data para un mu0 y kz0 fijos
graficar_omega = 0 # graficar plot del omega 

## 2 variables ############################

minimizar_omega_kz = 1 # data para un mu0,eta0 fijos

minimizar_omega_mu = 0 # data para un kz0,eta0 fijos

minimizar_omega_eta = 0 # data para un mu0,kz0 fijos
    
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real/real_freq_lorentziana','') 

try:
    sys.path.insert(1, path_basic)
    from det_conkz_lorentziana import determinante
except ModuleNotFoundError:
    print('det_conkz_lorentziana.py no se encuentra en ' + path_basic)
    path_basic2 = input('path de la carpeta donde se encuentra det_conkz_lorentziana.py')
    sys.path.insert(1, path_basic)
    from det_conkz_lorentziana import determinante

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()
aux_cte = c*1e-12

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 1

#%%

eta0 = 0.95
kz0 = 0
mu0 = 0.3

list_kz0 = np.linspace(0.1,1,10) # para mu0 fijo 
list_mu0 = np.linspace(0.1,0.9,9) # para un kz0 fijo


n = 251
list_kz = np.linspace(0, 2.5, n)
list_mu = np.linspace(0.1,0.9,n)
list_eta = np.linspace(0,1,n)

list_omegaTHz = np.linspace(5,100,n)
freq1 = 3.15 #THz
freq2 = 5.09 #THz
omegaTHz1 = freq1*2*np.pi
omegaTHz2 = freq2*2*np.pi

cond_inicial = [20]    

#%%        

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16
labelpady = 0
labelpadx = 0
pad = 0
lw = 1

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%

path_save = path_basic + '/' +  r'R_%.2f' %(R)
info = 'modo = %i, R = %.2f$\mu$m' %(modo,R)  

if save_data_opt==1:
    if not os.path.exists(path_save):
        print('Creating folder to save data')
        os.mkdir(path_save)

if minimizar_omega_kz == 1:
    inf = '$\mu_c$ = %.3f eV, $\eta$ = %.2f, ' %(mu0,eta0) + info + ', ' + name_this_py
    title = '$\mu_c$ = %.3f eV, $\eta$ = %.2f, ' %(mu0,eta0) + info + '\n' + name_this_py
    labely2 = '$\mu_c$ [eV]'
    inf_fig = 'det_modo%i_mu%.4f_eta%.2f' %(modo,mu0,eta0)
    
elif minimizar_omega_mu == 1:
    inf = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f, ' %(kz0,eta0) + info + ', ' + name_this_py
    title = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f, ' %(kz0,eta0) + info + '\n' + name_this_py
    labely2 = '$k_z$ [$\mu m^{-1}$]'
    inf_fig = 'det_modo%i_kz%.4f_eta%.2f' %(modo,kz0,eta0)
        
elif minimizar_omega_eta == 1:
    inf = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$, ' %(mu0,kz0) + info + ', ' + name_this_py
    title = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$, ' %(mu0,kz0) + info + '\n' + name_this_py
    labely2 = '$\eta$'
    inf_fig = 'det_modo%i_kz%.4f_mu%.2f' %(modo,kz0,mu0)

labelx = '$\omega$ [THz]'
labely1 = '|det|'

#%% 1 
   
print('Graficar y minimizar |det|')

tol_NM = 1e-10
ite_NM = 4000

if minimizar_omega_kz == 1:
    omegaTHz_opt = []
    kz_opt = []
    eq_det = []
    for kz_var in list_kz:
    # cond_inicial = list_cond_init[j]
        kz_var = float(kz_var)
        def det_1var(omegaTHz):
            
            omegac = omegaTHz/aux_cte    
            rta = determinante(kz_var,omegac,eta0,modo,R,mu0)
            return np.abs(rta)
        
        res = minimize(det_1var, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                          options={'maxiter':ite_NM})
     
        
        if res.message == 'Optimization terminated successfully.':

            # res.x = float(res.x)
            value = det_1var(res.x[0])
            eq_det.append(value)
            kz_opt.append(kz_var)
            omegaTHz_opt.append(res.x[0])
            
        cond_inicial = res.x[0]
  
    if save_data_opt == 1:
        os.chdir(path_save)
        print('Guardar data de minimizacion en .txt')
    
        tabla = np.array([kz_opt,omegaTHz_opt,eq_det])
        tabla = np.transpose(tabla)
        labeltxt = inf_fig + '.txt'
        header1 = inf
        header2 = 'kz [1/micrones]     Omega [THz]    |det|' + header1
        np.savetxt('opt_lorentz_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header2)

    del omegaTHz,kz,omegac
    def det_2var(omegaTHz,kz):
        omegac = omegaTHz/aux_cte    
        rta = determinante(kz,omegac,eta0,modo,R,mu0)
        
        return np.abs(rta)
    
    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')

    X, Y = np.meshgrid(list_omegaTHz,list_kz, sparse=True)
    f = np.vectorize(det_2var)
    Z = f(X, Y)
    
    vmin, vmax = np.min(Z), np.max(Z)
    
    # plt.plot(list_x,np.ones(N)*mu0,'--',lw = lw,color = 'green')
    # plt.plot(np.ones(N)*omega0,list_y,'--', lw = lw,color = 'green')
    
    pcm = plt.pcolormesh(X, Y, Z,
                      norm=colors.SymLogNorm(linthresh=0.5, linscale=0.5,
                                          vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    
    cbar = plt.colorbar(pcm, extend='both')
    plt.plot(omegaTHz_opt,kz_opt,'g.',ms =15)  
    # cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    cbar.set_label('|det|',fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig(inf_fig + '.png', format='png') 

#%% 2

elif minimizar_omega_mu == 1:
    omegaTHz_opt = []
    mu_opt = []
    eq_det = []
    for mu_var in list_mu:
    # cond_inicial = list_cond_init[j]
        def det_1var(omegaTHz):
            omegac = omegaTHz/aux_cte    
            rta = determinante(kz0,omegac,eta0,modo,R,mu_var)
            return np.abs(rta)
        
        res = minimize(det_1var, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                          options={'maxiter':ite_NM})
     
        
        if res.message == 'Optimization terminated successfully.':

            res.x = float(res.x)
            value = det_1var(res.x)
            eq_det.append(value)
            mu_opt.append(mu_var)
            omegaTHz_opt.append(res.x)
            
        cond_inicial = res.x
  
    if save_data_opt == 1:
        os.chdir(path_save)
        print('Guardar data de minimizacion en .txt')
    
        tabla = np.array([mu_opt,omegaTHz_opt,eq_det])
        tabla = np.transpose(tabla)
        labeltxt = inf_fig + '.txt'
        header1 = inf
        header2 = 'mu [eV]     Omega [THz]    |det|' + header1
        np.savetxt('opt_lorentz_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header2)

    def det_2var(omegaTHz,mu):
        omegac = omegaTHz/aux_cte    
        rta = determinante(kz0,omegac,eta0,modo,R,mu)
        return np.abs(rta)
    
    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')

    X, Y = np.meshgrid(list_omegaTHz,list_mu, sparse=True)
    f = np.vectorize(det_2var)
    Z = f(X, Y)
    
    vmin, vmax = np.min(Z), np.max(Z)
    
    # plt.plot(list_x,np.ones(N)*mu0,'--',lw = lw,color = 'green')
    # plt.plot(np.ones(N)*omega0,list_y,'--', lw = lw,color = 'green')
    
    pcm = plt.pcolormesh(X, Y, Z,
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=vmin, vmax=vmax),   
    				cmap='RdBu')
    
    
    cbar = plt.colorbar(pcm, extend='both')
    plt.plot(omegaTHz_opt,mu_opt,'g.',ms =15)  
    # cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    cbar.set_label('|det|',fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig(inf_fig + '.png', format='png') 
        
#%% 3

elif minimizar_omega_eta == 1:
    omegaTHz_opt = []
    eta_opt = []
    eq_det = []
    for eta_var in list_eta:
    # cond_inicial = list_cond_init[j]
        def det_1var(omegaTHz):
            omegac = omegaTHz/aux_cte    
            rta = determinante(kz0,omegac,eta_var,modo,R,mu0)
            return np.abs(rta)
        
        res = minimize(det_1var, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                          options={'maxiter':ite_NM})
     
        
        if res.message == 'Optimization terminated successfully.':

            res.x = float(res.x)
            value = det_1var(res.x)
            eq_det.append(value)
            eta_opt.append(eta_var)
            omegaTHz_opt.append(res.x)
            
        cond_inicial = res.x
  
    if save_data_opt == 1:
        os.chdir(path_save)
        print('Guardar data de minimizacion en .txt')
    
        tabla = np.array([eta_opt,omegaTHz_opt,eq_det])
        tabla = np.transpose(tabla)
        labeltxt = inf_fig + '.txt'
        header1 = inf
        header2 = 'eta     Omega [THz]    |det|' + header1
        np.savetxt('opt_lorentz_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header2)


    def det_2var(omegaTHz,eta):
        omegac = omegaTHz/aux_cte    
        rta = determinante(kz0,omegac,eta,modo,R,mu0)
        return np.abs(rta)

    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')

    X, Y = np.meshgrid(list_omegaTHz,list_eta, sparse=True)
    f = np.vectorize(det_2var)
    Z = f(X, Y)
    
    vmin, vmax = np.min(Z), np.max(Z)
    
    # plt.plot(list_x,np.ones(N)*mu0,'--',lw = lw,color = 'green')
    # plt.plot(np.ones(N)*omega0,list_y,'--', lw = lw,color = 'green')
    
    pcm = plt.pcolormesh(X, Y, Z,
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=vmin, vmax=vmax),   
    				cmap='RdBu')
    
    
    cbar = plt.colorbar(pcm, extend='both')
    plt.plot(omegaTHz_opt,eta_opt,'g.',ms =15)  
    # cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    cbar.set_label('|det|',fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
            
#%%