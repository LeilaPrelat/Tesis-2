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
modulo = 1 #si modulo == 1 ---> |Hz| (si modulo == 0 ---> Re(Hz))

non_active_medium = 1 #plotear campos con im(epsilon1) = 0
paper = 1  # sacar el titulo y guardar la info (formato paper)

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
    path_save0 = 'fields'  + '/' + 'barrido_R'
    if paper == 0:
        path_save = path_basic + '/' + path_save0
    else:
        path_save = path_basic + '/' + path_save0 + '/' + 'paper'
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from fieldsZ_sinkz import Ez,Hz
except ModuleNotFoundError:
    print('fieldsZ_sinkz.py no se encuentra en ' + path_basic)
    path_fields = input('path de la carpeta donde se encuentra fieldsZ_sinkz.py')
    sys.path.insert(1, path_fields)
    from fieldsZ_sinkz import Ez,Hz

#%%

print('Definir parametros del problema')

Ao,Bo = 1,1
hbaramu = 0.3
modo = 3

Ep = 0.3
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia
nmax = 10

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + 'epsiinf_DL_%.2f_vs_R' %(epsiinf_DL) + '/' + 'Ep_%.2f' %(Ep)
os.chdir(path_load)
name = 'opt_det_sinkz_vs_R_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)
    try:
        os.chdir(path_load + '/' + 'Ep_%.1f' %(Ep))
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError as error:
        print(error)

data_load = np.transpose(data_load)
[barrido_R,omegac_opt,epsi1_imag_opt,eq_det] = data_load

m = len(barrido_R)
index = 2860
index = int(index)
R = barrido_R[index]
omegac = omegac_opt[index]
delta_ci = epsi1_imag_opt[index]
print(R)

info1 = 'R = %.3f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, Ep = %.1f, $\gamma_{DL}$ = %.2f eV' %(R,hbaramu,epsiinf_DL,Ep,gamma_DL)
info2 = '$\Delta_{ci}$ = %.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i' %(delta_ci,omegac,modo,nmax)
title1 = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2  + '\n' + name_this_py
title2 = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2  + '\n' + name_this_py

labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
infotot = 'Ao = %i, Bo = %i' %(Ao,Bo) + ', ' + info1 +  ', ' + info2 +  ', ' + name_this_py

if non_active_medium == 1:
    info2_loss = '$\omega/c$ = %.2e 1/$\mu$m del modo = %i, nmax = %i' %(omegac,modo,nmax)
    title1_loss = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2_loss  + '\n' + name_this_py
    title2_loss = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

#%%

if hbaramu != 0.3:
    raise TypeError('Wrong value for mu_c')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')
    
#%%

print('Graficar el campo Hz para el medio 1 y 2')

if modulo==1:
    label = '$|H_z|$'
else:
    label = 'Re($H_z$)'

n2 = 250
cota = 2*R
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
limits = [min(x) , max(x), min(y) , max(y)]

if non_active_medium == 1: #epsi_ci = 0 ---> no hay medio activo
    def Hz_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Hz1,Hz2] = Hz(omegac,Ep,epsiinf_DL,gamma_DL,0,nmax,R,hbaramu,Ao,rho,phi)   
        if modulo==1:
            Hz1_tot, Hz2_tot = np.abs(Hz1), np.abs(Hz2)
        else:
            Hz1_tot, Hz2_tot = Hz1.real, Hz2.real
        if np.abs(rho)<=R: #medio1
            return Hz1_tot
        else: #medio2
            return Hz2_tot
        
    f1 = np.vectorize(Hz_2variable)
    Z = f1(X, Y)
    maxH = np.max(Z) 
    Z = Z/maxH
    
    plt.figure(figsize=tamfig)  
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    if paper == 0:
        plt.title(title1_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(label,fontsize=tamletra)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modHz_loss_modo%i_R%.3f.png' %(modo,R), format='png') 
        else:
            plt.savefig('reHz_loss_modo%i_R%.3f.png' %(modo,R), format='png') 


def Hz_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Hz1,Hz2] = Hz(omegac,Ep,epsiinf_DL,gamma_DL,delta_ci,nmax,R,hbaramu,Ao,rho,phi)   
    if modulo==1:
        Hz1_tot, Hz2_tot = np.abs(Hz1), np.abs(Hz2)
    else:
        Hz1_tot, Hz2_tot = Hz1.real, Hz2.real
    if np.abs(rho)<=R: #medio1
        return Hz1_tot
    else: #medio2
        return Hz2_tot
    
f1 = np.vectorize(Hz_2variable)
Z = f1(X, Y)
if non_active_medium == 1:
    Z = Z/maxH

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
if paper == 0:
    plt.title(title1,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(label,fontsize=tamletra)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modHz_modo%i_R%.3f.png' %(modo,R), format='png') 
    else:
        plt.savefig('reHz_modo%i_R%.3f.png' %(modo,R), format='png') 

#%%

print('Graficar el campo Ez para el medio 1 y 2')

if modulo==1:
    label = '$|E_z|$'
else:
    label ='Re($E_z$)'
    
if non_active_medium == 1: #epsi_ci = 0 ---> no hay medio activos
       
    def Ez_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Ez1,Ez2] = Ez(omegac,Ep,epsiinf_DL,gamma_DL,0,nmax,R,hbaramu,Bo,rho,phi)   
        if modulo==1:
            Ez1_tot, Ez2_tot = np.abs(Ez1), np.abs(Ez2)
        else:
            Ez1_tot, Ez2_tot = Ez1.real, Ez2.real
        if np.abs(rho)<=R: #medio1
            return Ez1_tot
        else: #medio2
            return Ez2_tot
    
    f2 = np.vectorize(Ez_2variable)
    Z = f2(X, Y)
    maxE = np.max(Z) 
    Z = Z/maxE
    
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    if paper == 0:
        plt.title(title2_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(label,fontsize=tamletra)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modEz_loss_modo%i_R%.3f.png' %(modo,R), format='png') 
        else:
            plt.savefig('reEz_loss_modo%i_R%.3f.png' %(modo,R), format='png') 

def Ez_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Ez1,Ez2] = Ez(omegac,Ep,epsiinf_DL,gamma_DL,delta_ci,nmax,R,hbaramu,Bo,rho,phi)   
    if modulo==1:
        Ez1_tot, Ez2_tot = np.abs(Ez1), np.abs(Ez2)
    else:
        Ez1_tot, Ez2_tot = Ez1.real, Ez2.real
    if np.abs(rho)<=R: #medio1
        return Ez1_tot
    else: #medio2
        return Ez2_tot
    
f2 = np.vectorize(Ez_2variable)
Z = f2(X, Y)
# Z = Z/np.max(Z)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
if paper == 0:
    plt.title(title2,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(label,size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modEz_modo%i_R%.3f.png' %(modo,R), format='png') 
    else:
        plt.savefig('reEz_modo%i_R%.3f.png' %(modo,R), format='png') 

if paper == 1:
    np.savetxt('info_fields_disp_modo%i_R%.3f.txt' %(modo,R), [infotot],fmt='%s')   
    
#%%
