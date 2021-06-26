#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

gain factor definido en el paper
    de Dionne
    
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/con_kz_real','') 
path_save = path_basic + '/' + 'gain_factor'
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene)
    from constantes import constantes

pi,hb,c,alfac,mu1,mu2,epsi2 = constantes()

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5 #micrones
hbaramu = 0.3        #eV mu_c
modo = 4

if modo != 1:
    kz_chicos = [0,0.05,0.1,0.13,0.135,0.14]
    kz_grandes = [0.5,0.75,1]
else:
    kz_chicos = [0,0.05,0.1,0.13]
    kz_grandes = [0.135,0.14,0.5,0.75,1]    

labely = 'g gain factor [1/$\mu$m]'
title = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f, modo = %i' %(R,re_epsi1,modo)
title = title + ', ' + name_this_py

path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1) 

#%%

if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('Importar los valores de SPASER')

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_chicos: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    
    plt.plot(omegac_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\omega/c$ $[1/\mu m]$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor1_vs_omegac_modo%i' %(modo)) 

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_chicos: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
        
    plt.plot(epsi1_imag_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('Im($\epsilon_1$)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor1_vs_im_epsi1_modo%i' %(modo)) 

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_chicos: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    plt.plot(list_mu_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\mu$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor1_vs_mu_modo%i' %(modo)) 
    
os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_chicos: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    plt.plot(list_mu_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title + '\n' + 'con perdidas',fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\mu$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor1_vs_mu_loss_modo%i' %(modo)) 

#%%

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_grandes: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    
    plt.plot(omegac_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\omega/c$ $[1/\mu m]$',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor2_vs_omegac_modo%i' %(modo)) 

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_grandes: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
        
    plt.plot(epsi1_imag_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('Im($\epsilon_1$)',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor2_vs_im_epsi1_modo%i' %(modo)) 

os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_grandes: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 + 1j*crit
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    plt.plot(list_mu_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\mu$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor2_vs_mu_modo%i' %(modo)) 


os.chdir(path_load)
plt.figure(figsize=tamfig)
for kz in kz_grandes: 
    
    name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data_load = np.transpose(data_load)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    gain_factor_tot = []
    for j in range(len(list_mu_opt)):
        crit = epsi1_imag_opt[j]
        omegac = omegac_opt[j] 
        epsi1 = re_epsi1 
        n1 = (epsi1*mu1)**(1/2)
        n1 = n1.real
        gain_factor = -omegac*crit/n1
        gain_factor_tot.append(gain_factor)
    plt.plot(list_mu_opt,gain_factor_tot,'.',ms=10, label = 'kz = %.4f 1/$\mu$m' %(kz))
plt.title(title + '\n' + 'con perdidas',fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel('$\mu$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.8))
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('gain_factor2_vs_mu_loss_modo%i' %(modo)) 

#%%
