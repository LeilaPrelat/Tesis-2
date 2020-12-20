#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos transversales + longitudinales:
    
    seccion 3.2 del cuaderno corto
    
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo
non_active_medium = 1 #plotear campos con im(epsilon1) = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'fields_2D'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
try:
    sys.path.insert(1, path_basic)
    from Htot_conkz import Htot
except ModuleNotFoundError:
    print('Htot_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra Htot_conkz.py')
    sys.path.insert(1, path_basic)
    from Htot_conkz import Htot

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
R = 0.5              # micrones
hbaramu = 0.3        # eV mu_c

list_modos = [1,2,3,4]

z = 0
phi = 0
nmax = 10
ind = 50

Ao,Bo = 1,1

#%%

if Ao*Bo != 0:
    path_save0 = path_save + '/' + '2pol'
elif Bo == 0:
    path_save0 = path_save + '/' + 'polAo'
elif Ao == 0:
    path_save0 = path_save + '/' + 'polBo'
    
if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

n = 500
list_rho = np.linspace(-2*R,2*R,n)
labelx = r'$\rho$ [$\mu$m]'
labely = '|H|$^2$'

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

print('')
print('Importar los valores de SPASER')

modo = 1
plt.figure(figsize=tamfig)
os.chdir(path_load)
name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)
try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)
print('')

data = np.transpose(data)
[list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
m = len(list_kz_opt)
ind = int(ind)
kz = list_kz_opt[ind] #micrones
print('kz = ', kz)
print('modo = ', modo)
print('')

im_epsi1 = epsi1_imag_opt[ind]
omegac = omegac_opt[ind] 
epsi1 = re_epsi1 + 1j*im_epsi1

info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m' %(kz,z,phi,R)
info2 = 'nmax = %i, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m' %(nmax,hbaramu,re_epsi1,-im_epsi1,omegac)
title = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + '\n' + name_this_py

if kz < 0.13:
    path_save = path_save0 + '/' + 'kz_chico'
else:
    path_save = path_save0 + '/' + 'kz_grande'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#print('Graficar el campo |Etot|^2 para el medio 1 y 2')

Htot_values = []
for rho in list_rho:
    rho = np.abs(rho)
    H1_tot,H2_tot = Htot(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    if rho <=R:
        Htot_values.append(H1_tot)
    else:
        Htot_values.append(H2_tot)

plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i' %(modo))
plt.title(title,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|Htot|2_modo%i_index%i.png' %(modo,ind))

#%%

plt.figure(figsize=tamfig)
os.chdir(path_load)
name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)
try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)
print('')

data = np.transpose(data)
[list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
m = len(list_kz_opt)
ind = int(ind)
kz = list_kz_opt[ind] #micrones
print('kz = ', kz)
print('modo = ', modo, 'con perdidas')
print('')

im_epsi1 = epsi1_imag_opt[ind]
omegac = omegac_opt[ind] 
epsi1 = re_epsi1 + 1j*im_epsi1        

info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m' %(kz,z,phi,R)
info2_loss = 'nmax = %i, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f, $\omega/c$ = %.2e 1/$\mu$m' %(nmax,hbaramu,re_epsi1,omegac)
title_loss = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

if kz < 0.13:
    path_save = path_save0 + '/' + 'kz_chico'
else:
    path_save = path_save0 + '/' + 'kz_grande'

Htot_values = []
for rho in list_rho:
    rho = np.abs(rho)
    H1_tot,H2_tot = Htot(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    if rho <=R:
        Htot_values.append(H1_tot)
    else:
        Htot_values.append(H2_tot)

plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i' %(modo))
plt.title(title_loss,fontsize=tamtitle)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 

if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|Htot|2_loss_modo%i_index%i.png' %(modo,ind))

#%%

print('')
print('Importar los valores de SPASER')

plt.figure(figsize=tamfig)
for modo in list_modos:
    os.chdir(path_load)
    name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)
    try:
        data = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('values de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    print('')
    
    data = np.transpose(data)
    [list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
    m = len(list_kz_opt)
    ind = int(ind)
    kz = list_kz_opt[ind] #micrones
    print('kz = ', kz)
    print('modo = ', modo)
    print('')
    
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = re_epsi1 + 1j*im_epsi1
    
    info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,phi,R,nmax,hbaramu)
    info2 = '$\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)
    title = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + '\n' + name_this_py
    
    if kz < 0.13:
        path_save = path_save0 + '/' + 'kz_chico'
    else:
        path_save = path_save0 + '/' + 'kz_grande'
    
    if save_graphs==1:    
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)
    
    #print('Graficar el campo |Etot|^2 para el medio 1 y 2')
    
    Htot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        H1_tot,H2_tot = Htot(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if rho <=R:
            Htot_values.append(H1_tot)
        else:
            Htot_values.append(H2_tot)
    
    plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i' %(modo))
    plt.title(title,fontsize=tamtitle)
    plt.ylabel(labely,fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
    plt.grid(1) 
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig('|Htot|2_index%i.png' %(ind))
    
if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    plt.figure(figsize=tamfig)
    for modo in list_modos:   
        
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)
        try:
            data = np.loadtxt(name,delimiter = '\t', skiprows=1)
            for line in (l.strip() for l in open(name) if l.startswith('#')):
                print('values de ', name, ':', line)
        except OSError or IOError:
            print('El archivo ' + name + ' no se encuentra en ' + path_load)
        print('')
        
        data = np.transpose(data)
        [list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
        m = len(list_kz_opt)
        ind = int(ind)
        kz = list_kz_opt[ind] #micrones
        print('kz = ', kz)
        print('modo = ', modo)
        print('')
        
        im_epsi1 = epsi1_imag_opt[ind]
        omegac = omegac_opt[ind] 
        epsi1 = re_epsi1 + 1j*im_epsi1        
        
        info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,phi,R,nmax,hbaramu)
        info2_loss = '$\epsilon_1$ = %.1f, $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,omegac,modo)
        title_loss = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

        if kz < 0.13:
            path_save = path_save0 + '/' + 'kz_chico'
        else:
            path_save = path_save0 + '/' + 'kz_grande'
        
        Htot_values = []
        for rho in list_rho:
            rho = np.abs(rho)
            H1_tot,H2_tot = Htot(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            if rho <=R:
                Htot_values.append(H1_tot)
            else:
                Htot_values.append(H2_tot)
        
        plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i' %(modo))
        plt.title(title_loss,fontsize=tamtitle)
        plt.ylabel(labely,fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
        plt.grid(1) 
        
        if save_graphs==1:
            os.chdir(path_save)
            plt.savefig('|Htot|2_loss_index%i.png' %(ind))

#%%
