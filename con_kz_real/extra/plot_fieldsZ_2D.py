#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos longitudinales:
    
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
path_save = path_basic + '/' + 'fields_2D'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
try:
    sys.path.insert(1, path_basic)
    from fieldsZ_conkz import fieldsZ
except ModuleNotFoundError:
    print('fieldsZ_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra fieldsZ_conkz.py')
    sys.path.insert(1, path_basic)
    from fieldsZ_conkz import fieldsZ

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

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

def names(modo):
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
    kz = list_kz_opt[ind] #micrones
    print('kz = ', kz)
    print('modo = ', modo)
    print('')
    
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = re_epsi1 + 1j*im_epsi1
    
    if kz < 0.13:
        path_save = path_save0 + '/' + 'kz_chico'
    else:
        path_save = path_save0 + '/' + 'kz_grande'
    
    if save_graphs==1:    
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)

    info1 = 'z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m' %(z,phi,R)
#    info2 = 'nmax = %i, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m' %(nmax,hbaramu,re_epsi1,-im_epsi1,omegac)
    info2 = 'nmax = %i, $\mu_c$ = %.1f eV, Re($\epsilon_1$) = %.1f' %(nmax,hbaramu,re_epsi1)
    title = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + ', ' + name_this_py    

    info2_loss = 'nmax = %i, $\mu_c$ = %.1f eV, Re($\epsilon_1$) = %.1f' %(nmax,hbaramu,re_epsi1)
    title_loss = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2_loss  + ', ' + name_this_py

    return kz,epsi1,omegac,path_save,title,title_loss

#%%

print('')
print('Importar los valores de SPASER')

plt.figure(figsize=tamfig)
for modo in list_modos:
    
    kz,epsi1,omegac,path_save,title,title_loss = names(modo)

    #print('Graficar el campo |Etot|^2 para el medio 1 y 2')
    
    Etot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if np.abs(rho)<=R:
            if modulo == 1:
                Etot_values.append(np.abs(Ez1))
            else:
                Etot_values.append(Ez1.real)
        else:
            if modulo == 1:
                Etot_values.append(np.abs(Ez2))
            else:
                Etot_values.append(Ez2.real)
    
    plt.plot(list_rho,Etot_values,'.',ms=10,label = 'modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
    plt.title(title,fontsize=tamtitle)
    
    if modulo == 1:
        plt.ylabel('|Ez|',fontsize=tamletra)
    else:
        plt.ylabel('Re(Ez)',fontsize=tamletra)
        
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
    plt.grid(1) 
    if save_graphs==1:
        os.chdir(path_save)
        if modulo == 1: 
            plt.savefig('|Ez|_index%i.png' %(ind))
        else:
            plt.savefig('ReEz_index%i.png' %(ind))

if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    plt.figure(figsize=tamfig)
    for modo in list_modos:   
        
        kz,epsi1,omegac,path_save,title,title_loss = names(modo)     
        
        Etot_values = []
        for rho in list_rho:
            rho = np.abs(rho)
            [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            if np.abs(rho)<=R:
                if modulo == 1:
                    Etot_values.append(np.abs(Ez1))
                else:
                    Etot_values.append(Ez1.real)
            else:
                if modulo == 1:
                    Etot_values.append(np.abs(Ez2))
                else:
                    Etot_values.append(Ez2.real)
        
        plt.plot(list_rho,Etot_values,'.',ms=10,label = 'modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
        plt.title(title_loss,fontsize=tamtitle)
        if modulo == 1:
            plt.ylabel('|Ez|',fontsize=tamletra)
        else:
            plt.ylabel('Re(Ez)',fontsize=tamletra)
        
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
        plt.grid(1) 
        
        if save_graphs==1:
            os.chdir(path_save)
            if modulo == 1: 
                plt.savefig('|Ez|_loss_index%i.png' %(ind))
            else:
                plt.savefig('ReEz_loss_index%i.png' %(ind))

#%%

print('')
print('Importar los valores de SPASER')

plt.figure(figsize=tamfig)
for modo in list_modos:
    
    kz,epsi1,omegac,path_save,title,title_loss = names(modo)

    #print('Graficar el campo |Etot|^2 para el medio 1 y 2')
    
    Htot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if np.abs(rho)<=R:
            if modulo == 1:
                Htot_values.append(np.abs(Hz1))
            else:
                Htot_values.append(Hz1.real)
        else:
            if modulo == 1:
                Htot_values.append(np.abs(Hz2))
            else:
                Htot_values.append(Hz2.real)
    
    plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
    plt.title(title,fontsize=tamtitle)
    
    if modulo == 1:
        plt.ylabel('|Hz|',fontsize=tamletra)
    else:
        plt.ylabel('Re(Hz)',fontsize=tamletra)
        
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
    plt.grid(1) 
    if save_graphs==1:
        os.chdir(path_save)
        if modulo == 1: 
            plt.savefig('|Hz|_index%i.png' %(ind))
        else:
            plt.savefig('ReHz_index%i.png' %(ind))
    
if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
    plt.figure(figsize=tamfig)
    for modo in list_modos:   
        
        kz,epsi1,omegac,path_save,title,title_loss = names(modo)     

        
        Htot_values = []
        for rho in list_rho:
            rho = np.abs(rho)
            [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            if np.abs(rho)<=R:
                if modulo == 1:
                    Htot_values.append(np.abs(Hz1))
                else:
                    Htot_values.append(Hz1.real)
            else:
                if modulo == 1:
                    Htot_values.append(np.abs(Hz2))
                else:
                    Htot_values.append(Hz2.real)
        
        plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
        plt.title(title_loss,fontsize=tamtitle)
        if modulo == 1:
            plt.ylabel('|Hz|',fontsize=tamletra)
        else:
            plt.ylabel('Re(Hz)',fontsize=tamletra)
        
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
        plt.grid(1) 
        
        if save_graphs==1:
            os.chdir(path_save)
            if modulo == 1: 
                plt.savefig('|Hz|_loss_index%i.png' %(ind))
            else:
                plt.savefig('ReHz_losss_index%i.png' %(ind))

#%%
