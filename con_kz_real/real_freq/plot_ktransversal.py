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
import matplotlib.pyplot as plt

#%% 

save_data_opt = 1 #guardar data de la minimizacion
save_graphs = 1 #guardar los graficos
paper = 1

if paper == 0:
    tamfig = (11,9)
    tamlegend = 18
    tamletra = 18
    tamtitle = 18
    tamnum = 16
    ms,ms2 = 1,10
    pad = 0
    labelpady = 0.2
    labelpadx = 0.2
else:
    tamfig = (3, 2)
    tamletra = 9
    tamnum = 8.5
    tamlegend = 8.5
    ms,ms2 = 1,2
    pad = -1.5
    import seaborn as sns
    sns.set()
# ticks_x = [0.4,0.6,0.8]
# ticks_y1_kz1 = [-0.036,-0.034,-0.032]
# ticks_y2_kz1 = [-0.012,-0.020,-0.028]

# ticks_y1_kz2 = [-0.761, -0.770, -0.779]
# ticks_y2_kz2 = [-0.012,-0.020,-0.028]
    
# columnspace = 0.5
# markerscale = 0.7
# loc1 = [0.275,0.9]  
# length_marker = 1
    labelpady = 1.1
    labelpadx = 0.7

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')

try:
    sys.path.insert(1, path_basic2 + '/extra')
    from def_kt import kt
except ModuleNotFoundError:
    print('def_kt.py no se encuentra en '+ path_basic2 + '/extra')
    path_basic2 = input('path de la carpeta donde se encuentra def_kt.py')
    sys.path.insert(1, path_basic2)
    from def_kt import kt

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
kz_real = 0.13      #eV mu_c
modo = 1
    
# kzlim = 0.14

path_load = path_basic + '/' + 're_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)
path_save = path_load

#%%

if R != 0.5:
    raise TypeError('Los .txt fueron hechos para R = 0.5 micrones')

#%%

os.chdir(path_load)
name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz_real,modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data_load = np.transpose(data_load)
[list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
#%%

label_graph = 'Opt det con kz = %.4f 1/$\mu$m' %(kz_real)
labelsinkz = 'Opt sin kz'
labelx = '$\mu_c$ [eV]'
title = 'Modo = %i, R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) 
title2 = title + '\n' + name_this_py
dpi = 300

kt_imag = [] # Imag(kt)
kt_real = [] # Real(kt)
kt_abs = [] #modulo del k transversal |kt|

kt_imag2 = [] # Imag(kt)
kt_real2 = [] # Real(kt)
kt_abs2 = [] #modulo del k transversal |kt|

lambbda_R_re = []
lambbda_R_im = []

mu2,epsi2 = 1,1
for j in range(len(list_mu_opt)):
    hbaramu = list_mu_opt[j]
    omegac = omegac_opt[j]
    lambbda = 2*np.pi/omegac
    
    factor_re = np.real(R/lambbda)
    factor_im = np.imag(R/lambbda)
    
    epsi1_imag = epsi1_imag_opt[j]
    epsi1 = re_epsi1 + 1j*epsi1_imag
    kt_value = kt(kz_real,omegac,epsi1,modo,R,hbaramu)    
    k0 = omegac
    xz = kz_real/k0
    Rbarra = R*k0
    xt_value2 = -xz**2 + mu2*epsi2
    xtmedio = (xt_value2)**(1/2)   
    if (xtmedio*Rbarra).real>=0:
	    xtmedio = xtmedio
    else:
	    xtmedio = -xtmedio
    
    kt_value2 = xtmedio*k0
    kt_abs2.append(np.abs(kt_value2))
    kt_imag2.append(kt_value2.imag)
    kt_real2.append(kt_value2.real)
    
    kt_abs.append(np.abs(kt_value))
    kt_imag.append(kt_value.imag)
    kt_real.append(kt_value.real)
    
    lambbda_R_re.append(factor_re)
    lambbda_R_im.append(factor_im)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_abs,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|kt|_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_real,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Re($k_t$) [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('kt_real_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_imag,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('kt_imag_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
 
    ########
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_abs2,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('|$k_t$| [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|kt|_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_real2,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Re($k_t$) [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('kt_real_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
    
plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,kt_imag2,'-r',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Im($k_t$) [1/$\mu$m]',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('kt_imag_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)


plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,lambbda_R_re,'-m',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Re($R/\lambda$)',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('lambdaR_re_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)

plt.figure(figsize=tamfig)
plt.plot(list_mu_opt,lambbda_R_im,'-m',lw=ms2,label=label_graph)
if paper == 0:
    plt.title(title2,fontsize=tamtitle)
    plt.legend(loc='best',markerscale=ms,fontsize=tamlegend)
    plt.grid(1)
plt.ylabel('Im($R/\lambda$)',fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if paper == 1:
    plt.tight_layout()
if save_graphs==1:
    plt.savefig('lambdaR_imag_vs_mu_kz%.4f_%i.png'%(kz_real,modo), format='png', dpi = dpi)
    
#%%
