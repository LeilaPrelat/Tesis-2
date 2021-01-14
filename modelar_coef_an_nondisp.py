#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

-modelar_coef_an.py : 
    modelar al coeficiente an como (omega-omega_inv)/(omega-omega_res)
 caso no dispersivo 
     - hacer graficos 2D
     - hacer un ajuste (1D) del modelo

"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sys
from scipy.interpolate import interp1d
# from scipy.optimize import curve_fit

#%% 

graph_3D = 0 #graficar 2D el coef an y el modelo para el coef an (graficos de color)
log_scale = 0
save_graphs = 1 #guardar los graficos

tamfig = (12,9)
tamlegend = 16
tamletra = 16
tamtitle = 16
tamnum = 14

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_load_inv = path_basic + '/sin_kz_inv/non-dispersive/complex_freq'
path_load_res = path_basic + '/sin_kz/non-dispersive/complex_freq'

try:
    sys.path.insert(1, path_basic)
    from coef_an import an_nondisp
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra coef_an.py')
    sys.path.insert(1, path_basic)
    from coef_an import an_nondisp

#%%

list_color = ['purple','darkblue','darkgreen']

print('Definir parametros del problema')

modo = 2
R = 0.5 #micrones
re_epsi1 = 4.9
mu = 0.3
Ao = 1

info0 = 'Modo = %i, R = %.2f $\mu$m, Re($\epsilon_1$) = %.2f, $\mu_c$ = %.2f eV, Ao = %i' %(modo,R,re_epsi1,mu,Ao)
info = info0  +  ', ' + name_this_py
path_load = '/re_epsi1_%.2f_vs_mu/find_Lambda/modo_%i' %(re_epsi1,modo)
path_save = '/coef_an_nondisp'

title = info0 +'\n' + name_this_py
labelx = 'Re($\omega/c$)'
labely = 'Im($\omega/c$)'

#%%
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')

#%%

def list_omega(typpe,index,tol,N):

    if typpe == 'res':
        # print('Importar los omega complejos de resonancia')
        os.chdir(path_load_res + path_load) 

    elif typpe == 'inv':
        # print('Importar los omega complejos de invisibilidad')
        os.chdir(path_load_inv + path_load)
    else:
        raise TypeError('Wrong type for frequency: res or inv')    
        
    tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
    tabla = np.transpose(tabla)       
    try: 
        [list_epsi1_imag,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
    except ValueError: 
        [list_epsi1_imag,lambda_real_opt,lambda_imag_opt] = tabla
    
    list_omegac_real = []
    list_omegac_imag = []
    for k in range(len(list_epsi1_imag)):
        lambda_real = lambda_real_opt[k]
        lambda_imag = lambda_imag_opt[k]
        lambbda = lambda_real + 1j*lambda_imag
        omega_c = 2*np.pi/lambbda
        list_omegac_real.append(omega_c.real)
        list_omegac_imag.append(omega_c.imag)    

    omegac_real = list_omegac_real[index]
    omegac_imag = list_omegac_imag[index]    
    
    omegac12_real,omegac22_real = omegac_real*(1-tol), omegac_real*(1+tol)
    omegac12_imag,omegac22_imag = omegac_imag-tol, omegac_imag + tol
    list_omegac_real = np.linspace(omegac12_real, omegac22_real,N)
    list_omegac_imag = np.linspace(omegac12_imag, omegac22_imag,N)
    
    return list_omegac_real,list_omegac_imag,list_epsi1_imag, omegac_real,omegac_imag

#%%

index = 0
tol = 1e-3
N = 200

listx,listy,list_epsi1_imag, omegac_real,omegac_imag = list_omega('inv',index,tol,N)
omegac_inv = omegac_real + 1j*omegac_imag 

listx,listy,list_epsi1_imag, omegac_real,omegac_imag = list_omega('res',index,tol,N)
omegac_res = omegac_real + 1j*omegac_imag 

epsi1_imag_opt = list_epsi1_imag[index]
epsi1 = re_epsi1 + 1j*epsi1_imag_opt

def an_modelo1D(x,y):
    omegac = x + 1j*y
    num = omegac - omegac_inv
    den = omegac - omegac_res
    return np.abs(num/den)

def coef_an2D(x,y):
    omegac = x + 1j*y
    an_value = an_nondisp(omegac,epsi1,modo,R,mu,Ao)
    return np.abs(an_value)    


X, Y = np.meshgrid(listx, listy, sparse=True)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

#%%

if graph_3D == 1:
    print('Grafico 2D del coeficiente an cerca de los valores criticos')
    f = np.vectorize(coef_an2D)
    Z = f(X, Y)
    vmin,vmax = np.min(Z), np.max(Z)
    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
    
    if vmin < 0 :
        tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                          + [0] 
                          + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
    else:
        tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = tamtitle)
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel(labely,fontsize=int(tamletra*1.2))
    plt.tick_params(labelsize = tamnum)
    plt.tick_params(labelsize = tamnum)
    plt.plot(listx,np.ones(N)*omegac_imag,'--',color = 'green')
    plt.plot(np.ones(N)*omegac_real,listy,'--',color = 'green')
    if log_scale == 1:
        pcm = plt.pcolormesh(X, Y, Z,norm = SymLogNorm(linthresh=0.03, linscale=0.03,
                                                    vmin=vmin, vmax=vmax), cmap='RdBu_r')
    else:
        pcm = plt.pcolormesh(X, Y, Z, cmap='RdBu_r')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_ticks(tick_locations)
    cbar.set_label('|an|',fontsize=tamletra)
    
    if save_graphs==1:
        #plt.tight_layout(pad=wspace)
        os.chdir(path_basic + path_save)
        plt.savefig('coef_an%i'%(modo))
    
    print('Grafico 2D del MODELO del coeficiente an cerca de los valores criticos')
    
    f = np.vectorize(an_modelo1D)
    Z = f(X, Y)
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = tamtitle)
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel(labely,fontsize=int(tamletra*1.2))
    plt.tick_params(labelsize = tamnum)
    plt.tick_params(labelsize = tamnum)
    plt.plot(listx,np.ones(N)*omegac_imag,'--',color = 'green')
    plt.plot(np.ones(N)*omegac_real,listy,'--',color = 'green')
    pcm = plt.pcolormesh(X, Y, Z, cmap='RdBu_r')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label('$|\dfrac{\omega/c - \omega_{zero}/c}{\omega/c - \omega_{polo}/c}|$',fontsize=tamletra)
    if save_graphs==1:
        #plt.tight_layout(pad=wspace)
        plt.savefig('modelo_coef_an%i'%(modo))

#%%

index = 100
N = 400
graficar = 0

list_cte_real_crit = []
list_cte_imag_crit = []
    
for index in range(len(list_epsi1_imag)):

    listx,listy,list_epsi1_imag, omegac_real,omegac_imag = list_omega('res',index,tol,N)
    omegac_res = omegac_real + 1j*omegac_imag 
    
    epsi1 = re_epsi1 + 1j*epsi1_imag_opt
    
    def an_modelo_1D(x):
        # x es la parte real de omega/c
        omegac = x + 1j*omegac_imag
        num = omegac - omegac_inv
        den = omegac - omegac_res
        return num/den
    
    def coef_an1D(x):
        # x es la parte real de omega/c
        omegac = x + 1j*omegac_imag
        an_value = an_nondisp(omegac,epsi1,modo,R,mu,Ao)
        return an_value 
    
    
    list_an_ajuste = [] # lo que queda al multiplicar por la cte compleja
    list_an = []
    list_cte_real = []
    list_cte_imag = []
    
    for value in listx:
        cte = coef_an1D(value)/an_modelo_1D(value)
        list_cte_real.append(cte.real)
        list_cte_imag.append(cte.imag)
        
        list_an_ajuste.append(np.abs(an_modelo_1D(value)*cte))
        list_an.append(np.abs(coef_an1D(value)))
    
    f_cte_real = interp1d(listx,list_cte_real)
    f_cte_imag = interp1d(listx,list_cte_imag)
    cte_real_crit = f_cte_real(omegac_real)
    cte_imag_crit = f_cte_imag(omegac_real)
    list_cte_real_crit.append(cte_real_crit)
    list_cte_imag_crit.append(cte_imag_crit)
    
    if graficar == 1:
    
        plt.figure(figsize=tamfig)
        plt.title(title, fontsize = tamtitle)
        plt.xlabel(labelx,fontsize=int(tamletra*1.2))
        plt.plot(np.ones(len(listx))*omegac_real,list_cte_imag,'--',color = 'green')
        plt.ylabel('Im(cte multiplicativa)',fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.plot(listx,list_cte_imag,'b.',ms = 10)
        plt.grid(1)
        
        plt.figure(figsize=tamfig)
        plt.title(title, fontsize = tamtitle)
        plt.xlabel(labelx,fontsize=int(tamletra*1.2))
        plt.plot(np.ones(len(listx))*omegac_real,list_cte_real,'--',color = 'green')
        plt.ylabel('Re(cte multiplicativa)',fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.plot(listx,list_cte_real,'m.',ms = 10)
        plt.grid(1)
        
        plt.figure(figsize=tamfig)
        plt.title(title, fontsize = tamtitle)
        plt.xlabel(labelx,fontsize=int(tamletra*1.2))
        plt.tick_params(labelsize = tamnum)
        plt.plot(np.ones(len(listx))*omegac_real,list_an,'--',color = 'green')
        plt.plot(listx,list_an_ajuste,'m.',ms = 10, label = '|(modelo de an)*cte|')
        plt.plot(listx,list_an,'b.',ms = 10,label = '|an|')
        plt.legend(loc = 'best',markerscale = 2,fontsize=tamlegend)
        plt.grid(1)

#%%

plt.figure(figsize=tamfig)
plt.title(title, fontsize = tamtitle)
plt.xlabel('Im($\epsilon_1$)',fontsize=int(tamletra*1.2))
plt.ylabel('Re(cte multiplicativa)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
if modo > 2:
    plt.yscale('log')
plt.plot(list_epsi1_imag,list_cte_real_crit,'m.',ms = 10)
plt.grid(1)

plt.figure(figsize=tamfig)
plt.title(title, fontsize = tamtitle)
plt.xlabel('Im($\epsilon_1$)',fontsize=int(tamletra*1.2))
plt.ylabel('Im(cte multiplicativa)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
if modo > 2:
    plt.yscale('log')
plt.plot(list_epsi1_imag,list_cte_imag_crit,'b.',ms = 10)
plt.grid(1)


#%%




