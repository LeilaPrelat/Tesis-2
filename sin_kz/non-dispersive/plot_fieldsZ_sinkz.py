#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

obtuve gn (minimizandolo) ---> obtuve
los denominadores de los coef an y cn
(entonces obtuve an y cn) 
---> podemos graficar los campos

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt

#%%

save_graphs = 1
modulo = 0 #si modulo == 1 ---> |Hz| (si modulo == 0 ---> Re(Hz))

non_active_medium = 1 #plotear campos con im(epsilon1) = 0

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    path_save = path_basic + '/' + 'fields'
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from fieldsZ_sinkz import Ez,Hz
except ModuleNotFoundError:
    print('fieldsZ_sinkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_fields = input('path de la carpeta donde se encuentra fieldsZ_sinkz.py')
    sys.path.insert(1, path_fields)
    from fieldsZ_sinkz import Ez,Hz

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
Ao,Bo = 1,1
modo = 1
R = 0.5 #micrones 
nmax = 10

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1)
os.chdir(path_load)
name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en el path_load')

data_load = np.transpose(data_load)
[barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load

index = 0
hbaramu = barrido_mu[index]
omegac = omegac_opt[index]
im_epsi1 = epsi1_imag_opt[index]
epsi1 = re_epsi1 + 1j*im_epsi1

info1 = 'R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(R,nmax,hbaramu)
info2 = '$\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)
title1 = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2  + '\n' + name_this_py
title2 = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2  + '\n' + name_this_py

if non_active_medium == 1:
    info2_loss = '$\epsilon_1$ = %.1f, $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,omegac,modo)
    title1_loss = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2_loss  + '\n' + name_this_py
    title2_loss = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

#%%

if R != 0.5:
    raise TypeError('Wrong value for R')
    
#%%

print('Graficar el campo Hz para el medio 1 y 2')

def Hz_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Hz1,Hz2] = Hz(omegac,epsi1,nmax,R,hbaramu,Ao,rho,phi)   
    if modulo==1:
        Hz1_tot, Hz2_tot = np.abs(Hz1), np.abs(Hz2)
    else:
        Hz1_tot, Hz2_tot = Hz1.real, Hz2.real
    if np.abs(rho)<=R: #medio1
        return Hz1_tot
    else: #medio2
        return Hz2_tot
    
n2 = 100
cota = 2*R
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
f1 = np.vectorize(Hz_2variable)
Z = f1(X, Y)
# Z = Z/np.max(Z)
labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'

plt.figure(figsize=tamfig)
limits = [min(x) , max(x), min(y) , max(y)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.title(title1,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
if modulo==1:
    cbar.set_label('|Hz|',size=tamlegend)
else:
    cbar.set_label('Re(Hz)',size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modHz_modo%i_mu%.4f' %(modo,hbaramu), format='png')  
    else:
        plt.savefig('reHz_modo%i_mu%.4f' %(modo,hbaramu), format='png')  

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    def Hz_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Hz1,Hz2] = Hz(omegac,re_epsi1,nmax,R,hbaramu,Ao,rho,phi)   
        if modulo==1:
            Hz1_tot, Hz2_tot = np.abs(Hz1), np.abs(Hz2)
        else:
            Hz1_tot, Hz2_tot = Hz1.real, Hz2.real
        if np.abs(rho)<=R: #medio1
            return Hz1_tot
        else: #medio2
            return Hz2_tot
        
    n2 = 100
    cota = 2*R
    x = np.linspace(-cota,cota,n2)
    y = np.linspace(-cota,cota,n2)
    X, Y = np.meshgrid(x, y)
    f1 = np.vectorize(Hz_2variable)
    Z = f1(X, Y)
    # Z = Z/np.max(Z)
    labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
    
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.title(title1_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    if modulo==1:
        cbar.set_label('|Hz|',size=tamlegend)
    else:
        cbar.set_label('Re(Hz)',size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modHz_loss_modo%i_mu%.4f' %(modo,hbaramu), format='png')  
        else:
            plt.savefig('reHz_loss_modo%i_mu%.4f' %(modo,hbaramu), format='png')  

#%%

print('Graficar el campo Ez para el medio 1 y 2')

def Ez_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Ez1,Ez2] = Ez(omegac,epsi1,nmax,R,hbaramu,Bo,rho,phi)   
    if modulo==1:
        Ez1_tot, Ez2_tot = np.abs(Ez1), np.abs(Ez2)
    else:
        Ez1_tot, Ez2_tot = Ez1.real, Ez2.real
    if np.abs(rho)<=R: #medio1
        return Ez1_tot
    else: #medio2
        return Ez2_tot
    
# x = np.linspace(-cota,cota,n2)
# y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
f2 = np.vectorize(Ez_2variable)
Z = f2(X, Y)
# Z = Z/np.max(Z)

plt.figure(figsize=tamfig)
limits = [min(x) , max(x), min(y) , max(y)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.title(title2,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
if modulo==1:
    cbar.set_label('|Ez|',size=tamlegend)
else:
    cbar.set_label('Re(Ez)',size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modEz_modo%i_mu%.4f' %(modo,hbaramu), format='png')  
    else:
        plt.savefig('reEz_modo%i_mu%.4f' %(modo,hbaramu), format='png')  

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
       
    def Ez_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Ez1,Ez2] = Ez(omegac,re_epsi1,nmax,R,hbaramu,Bo,rho,phi)   
        if modulo==1:
            Ez1_tot, Ez2_tot = np.abs(Ez1), np.abs(Ez2)
        else:
            Ez1_tot, Ez2_tot = Ez1.real, Ez2.real
        if np.abs(rho)<=R: #medio1
            return Ez1_tot
        else: #medio2
            return Ez2_tot
        
    # x = np.linspace(-cota,cota,n2)
    # y = np.linspace(-cota,cota,n2)
    X, Y = np.meshgrid(x, y)
    f2 = np.vectorize(Ez_2variable)
    Z = f2(X, Y)
    # Z = Z/np.max(Z)
    
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.title(title2_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    if modulo==1:
        cbar.set_label('|Ez|',size=tamlegend)
    else:
        cbar.set_label('Re(Ez)',size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modEz_loss_modo%i_mu%.4f' %(modo,hbaramu), format='png')  
        else:
            plt.savefig('reEz_loss_modo%i_mu%.4f' %(modo,hbaramu), format='png')  
    
#%%
