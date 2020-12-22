#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos longitudinales para el caso con kz:
    
    seccion 3.2 del cuaderno corto

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo
modulo = 1 #si modulo = 0 ---> grafica partes reales, si modulo = 1 grafica el modulo de los campos

non_active_medium = 1 #plotear campos con im(epsilon1) = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
name_this_py = 'Ver ' + name_this_py
path_graphene = path_basic.replace('/' + 'con_kz_real','') 

if save_graphs==1:
    path_save = path_basic + '/' + 'fields_3D'
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from fieldsZ_conkz import fieldsZ
except ModuleNotFoundError:
    print('fieldsZ_conkz.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra fieldsZ_conkz.py')
    sys.path.insert(1, path_basic)
    from fieldsZ_conkz import fieldsZ

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
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
modo = 1

z = 0
nmax = 8
Ao,Bo = 0,1

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

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
ind = 50
ind = int(ind)
kz = list_kz_opt[ind] #micrones
print('kz = ', kz)
print('modo = ', modo)
print('')

im_epsi1 = epsi1_imag_opt[ind]
omegac = omegac_opt[ind] 
epsi1 = re_epsi1 + 1j*im_epsi1

# ktot = omegac*((mu2*epsi2)**(1/2))
# ang_theta = np.arcsin(kz/np.abs(ktot)) + np.pi/2
# ang2 = np.cos(np.pi/2 - ang_theta) 

# print(360*ang_theta/(2*np.pi))
# print(360*ang2/(2*np.pi))


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

print('Graficar el campo Hz para el medio 1 y 2')

def Hz_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
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
plt.title(title,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
if modulo==1:
    cbar.set_label('|Hz|',size=tamlegend)
else:
    cbar.set_label('Re(Hz)',size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modHz_modo%i_kz%.4f' %(modo,kz), format='png') 
    else:
        plt.savefig('reHz_modo%i_kz%.4f' %(modo,kz), format='png') 

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    def Hz_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
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
    plt.title(title_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    if modulo==1:
        cbar.set_label('|Hz|',size=tamlegend)
    else:
        cbar.set_label('Re(Hz)',size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modHz_loss_modo%i_kz%.4f' %(modo,kz), format='png') 
        else:
            plt.savefig('reHz_loss_modo%i_kz%.4f' %(modo,kz), format='png') 
            
#%%

print('Graficar el campo Ez para el medio 1 y 2')

def Ez_2variable(x,y):   
    phi = np.arctan2(y,x)
    rho = (x**2+y**2)**(1/2)
    [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
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
plt.title(title,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
if modulo==1:
    cbar.set_label('|Ez|',size=tamlegend)
else:
    cbar.set_label('Re(Ez)',size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modEz_modo%i_kz%.4f' %(modo,kz), format='png') 
    else:
        plt.savefig('reEz_modo%i_kz%.4f' %(modo,kz), format='png') 

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
       
    def Ez_2variable(x,y):   
        phi = np.arctan2(y,x)
        rho = (x**2+y**2)**(1/2)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
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
    plt.title(title_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    if modulo==1:
        cbar.set_label('|Ez|',size=tamlegend)
    else:
        cbar.set_label('Re(Ez)',size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modEz_loss_modo%i_kz%.4f' %(modo,kz), format='png') 
        else:
            plt.savefig('reEz_loss_modo%i_kz%.4f' %(modo,kz), format='png') 
    

#%%
