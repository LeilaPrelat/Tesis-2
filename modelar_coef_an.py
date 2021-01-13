#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

-modelar_coef_an.py : modelar al coeficiente an como (omega-omega_inv)/(omega-omega_res)
 caso no dispersivo

"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import sys

#%% 

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
    from coef_an import an
except ModuleNotFoundError:
    print('graphene_sigma.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra coef_an.py')
    sys.path.insert(1, path_basic)
    from coef_an import an

#%%

list_color = ['purple','darkblue','darkgreen']

print('Definir parametros del problema')

modo = 1
R = 0.5 #micrones
re_epsi1 = 4.9
mu = 0.3
Ao = 1

info = 'Modo = %i, R = %.2f $\mu$m, Re($\epsilon_1$) = %.2f, $\mu_c$ = %.2f eV, Ao = %i' %(modo,R,re_epsi1,mu,Ao)
info = info  +  ', ' + name_this_py
path_load = '/re_epsi1_%.2f_vs_mu/find_Lambda/modo_%i' %(re_epsi1,modo)
path_save = '/coef_an'

title = info +'\n' + name_this_py
labelx = 'Re($\omega/c$)'
labely = 'Im($\omega/c$)'

#%%
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')
    
#%%

print('res: Sacar los valores positivos de Im(omega/c)')

list_im_epsi1_res = []
index_crit_res = []

os.chdir(path_load_res + path_load)
tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
tabla = np.transpose(tabla)       
try: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
except ValueError: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
    
omega_c_real = []
omega_c_imag = []
for k in range(len(lambda_real_opt)):
    lambda_real = lambda_real_opt[k]
    lambda_imag = lambda_imag_opt[k]
    lambbda = lambda_real + 1j*lambda_imag
    omega_c = 2*np.pi/lambbda
    omega_c_real.append(omega_c.real)
    omega_c_imag.append(omega_c.imag)

#hallar en donde se anula la parte imaginaria de omega/c (valor critico)

omega_c_imag2 = np.abs(omega_c_imag)
ind3 = np.argmin(omega_c_imag2)

list_im_epsi1_res.append(epsi1_imag_opt[ind3])
index_crit_res.append(ind3)

#%%

print('inv: Sacar los valores positivos de Im(omega/c)')

list_im_epsi1_inv = []
index_crit_inv = []

os.chdir(path_load_inv + path_load)
tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
tabla = np.transpose(tabla)       
try: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
except ValueError: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
    
omega_c_real = []
omega_c_imag = []
for k in range(len(lambda_real_opt)):
    lambda_real = lambda_real_opt[k]
    lambda_imag = lambda_imag_opt[k]
    lambbda = lambda_real + 1j*lambda_imag
    omega_c = 2*np.pi/lambbda
    omega_c_real.append(omega_c.real)
    omega_c_imag.append(omega_c.imag)

#hallar en donde se anula la parte imaginaria de omega/c (valor critico)

omega_c_imag2 = np.abs(omega_c_imag)
ind3 = np.argmin(omega_c_imag2)

list_im_epsi1_inv.append(epsi1_imag_opt[ind3])
index_crit_inv.append(ind3)
    
#%%        

print('Importar los omega complejos de resonancia y de invisibilidad')

os.chdir(path_load_res + path_load) 
tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
tabla = np.transpose(tabla)       
try: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
except ValueError: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla

omega_c_real = []
omega_c_imag = []
epsi1_imag_opt_res = []
for k in range(index_crit_res[0]):
    lambda_real = lambda_real_opt[k]
    lambda_imag = lambda_imag_opt[k]
    lambbda = lambda_real + 1j*lambda_imag
    omega_c = 2*np.pi/lambbda
    omega_c_real.append(omega_c.real)
    omega_c_imag.append(omega_c.imag)    
    epsi1_imag_opt_res.append(epsi1_imag_opt[k])


os.chdir(path_load_inv + path_load) 
tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
tabla = np.transpose(tabla)       
try: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
except ValueError: 
    [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla

omega_c_real2 = []
omega_c_imag2 = []
epsi1_imag_opt_inv = []
for k in range(index_crit_inv[0]):
    lambda_real = lambda_real_opt[k]
    lambda_imag = lambda_imag_opt[k]
    lambbda = lambda_real + 1j*lambda_imag
    omega_c = 2*np.pi/lambbda
    omega_c_real2.append(omega_c.real)
    omega_c_imag2.append(omega_c.imag)    
    epsi1_imag_opt_inv.append(epsi1_imag_opt[k])

#%%

index = 0
tol = 1e-3
N = 200


def list_omega(typpe,index,tol,N):
    if typpe == 'res':
        omegac_real = omega_c_real[index]
        omegac_imag = omega_c_imag[index]
    elif typpe == 'inv':
        omegac_real = omega_c_real2[index]
        omegac_imag = omega_c_imag2[index]
    else:
        raise TypeError('Wrong type for frequency: res or inv')
       
    omegac12_real,omegac22_real = omegac_real*(1-tol), omegac_real*(1+tol)
    omegac12_imag,omegac22_imag = omegac_imag-tol, omegac_imag + tol
    list_omegac_real = np.linspace(omegac12_real, omegac22_real,N)
    list_omegac_imag = np.linspace(omegac12_imag, omegac22_imag,N)
    return list_omegac_real,list_omegac_imag,omegac_real,omegac_imag

listx,listy,omegac_real,omegac_imag = list_omega('inv',index,tol,N)
omegac_inv = omegac_real + 1j*omegac_imag 
listx,listy,omegac_real,omegac_imag = list_omega('res',index,tol,N)
omegac_res = omegac_real + 1j*omegac_imag 

epsi1 = re_epsi1 + 1j*epsi1_imag_opt_res[index]

def an_modelo(x,y):
    omegac = x + 1j*y
    num = omegac - omegac_inv
    den = omegac - omegac_res
    return np.abs(num/den)

def coef_an2D(x,y):
    omegac = x + 1j*y
    an_value = an(omegac,epsi1,modo,R,mu,Ao)
    return np.abs(an_value)    


X, Y = np.meshgrid(listx, listy, sparse=True)
limits = [min(listx) , max(listx), min(listy) , max(listy)]

#%%

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

#%%

f = np.vectorize(an_modelo)
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

