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
path_save = path_basic + '/' + 'fields_3D'

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
R = 0.5 #micrones
hbaramu = 0.3        #eV mu_c
modo = 4

z = 0
nmax = 8
Ao,Bo = 1,1

#%%

print('Importar los valores de SPASER')

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 
os.chdir(path_load)
name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)

try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)


data = np.transpose(data)
[list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
m = len(list_kz_opt)
ind = -1
ind = int(ind)
kz = list_kz_opt[ind] #micrones
print('kz = ', kz)
print('modo = ', modo)
print('')

im_epsi1 = epsi1_imag_opt[ind]
omegac = omegac_opt[ind] 
epsi1 = re_epsi1 + 1j*im_epsi1

info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,R,nmax,hbaramu)
info2 = '$\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)
title = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + '\n' + name_this_py

if non_active_medium == 1:
    info2_loss = '$\epsilon_1$ = %.1f, $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,omegac,modo)
    title_loss = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

if Ao*Bo != 0:
    path_save = path_save + '/' + '2pol'
elif Bo == 0:
    path_save = path_save + '/' + 'polAo'
elif Ao == 0:
    path_save = path_save + '/' + 'polBo'
    
if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

if kz < 0.13:
    path_save = path_save + '/' + 'kz_chico'
else:
    path_save = path_save + '/' + 'kz_grande'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('Graficar el campo |Htot|^2 para el medio 1 y 2')

def Htot_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [H1_tot,H2_tot] = Htot(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    if np.abs(rho) <= R: #medio1
        return H1_tot
    else: #medio2
        return H2_tot

n2 = 100
cota = 2*R
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
f1 = np.vectorize(Htot_2variable)
Z1 = f1(X, Y)
# Z = Z/np.max(Z)
labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
labelz = '|Htot|$^2$'

plt.figure(figsize=tamfig)
limits = [min(x) , max(x), min(y) , max(y)]
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.title(title,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z1, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.set_label(labelz,size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|Htot|2_modo%i_kz%.4f' %(modo,kz), format='png') 

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    def Htot_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [H1_tot,H2_tot] = Htot(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if np.abs(rho) <= R: #medio1
            return H1_tot
        else: #medio2
            return H2_tot
        

    # x = np.linspace(-cota,cota,n2)
    # y = np.linspace(-cota,cota,n2)
    # X, Y = np.meshgrid(x, y)
    f1 = np.vectorize(Htot_2variable)
    Z2 = f1(X, Y)
    # Z = Z/np.max(Z)
    
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.title(title_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z2, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.set_label(labelz,size=tamlegend)
    
    if save_graphs==1:
        plt.savefig('|Htot|2_loss_modo%i_kz%.4f' %(modo,kz), format='png') 


del X,Y,Z1,Z2
           
#%%
