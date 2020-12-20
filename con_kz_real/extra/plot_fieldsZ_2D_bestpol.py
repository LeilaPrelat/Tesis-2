#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos longitudinales:
    
    seccion 3.2 del cuaderno corto

para la mejor polarizacion (autoestados del caso inhomogeneo)    

"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
import math

save_graphs = 1 #guardar los graficos 2D del campo
modulo = 1 #si modulo = 0 ---> grafica partes reales, si modulo = 1 grafica el modulo de los campos

normal_pol = 1 #plotear campos con polarizacion normal (la del principio)

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'fields_2D_bestpol'

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
ind = -1 #0 50 -20 -1

Ao,Bo = 0,1

#%%

if Ao*Bo != 0:
    path_save0 = path_save + '/' + '2pol'
elif Bo == 0:
    path_save0 = path_save + '/' + 'polAo'
elif Ao == 0:
    path_save0 = path_save + '/' + 'polBo'
    
if save_graphs==1:    
    if not os.path.exists(path_save0):
        print('Creating folder to save graphs')
        os.mkdir(path_save0)

n = 500
list_rho = np.linspace(-2*R,2*R,n)
labelx = r'$\rho$ [$\mu$m]'

#path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 
path_loadAoBO = path_basic  + '/' + 'best_pol'

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

print('')
print('Importar los valores de SPASER')

def names(modo):
    os.chdir(path_loadAoBO)
    if Ao*Bo != 0:
        nameAo1 = 'C1best_polAoBo_vs_kz_modo%i.txt' %(modo)
        nameBo1 = 'A1best_polAoBo_vs_kz_modo%i.txt' %(modo)
    elif Bo == 0:
        nameAo1 = 'C1best_polAo_vs_kz_modo%i.txt' %(modo)
        nameBo1 = 'A1best_polAo_vs_kz_modo%i.txt' %(modo)
    elif Ao == 0:
        nameAo1 = 'C1best_polBo_vs_kz_modo%i.txt' %(modo)
        nameBo1 = 'A1best_polBo_vs_kz_modo%i.txt' %(modo)
        
        
    try:
        dataAo1 = np.loadtxt(nameAo1,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(nameAo1) if l.startswith('#')):
            print('values de ', nameAo1, ':', line)
    except OSError or IOError:
        print('El archivo ' + nameAo1 + ' no se encuentra en ' + path_loadAoBO)
    print('')
    
    dataAo1 = np.transpose(dataAo1)
    [list_re_Ao1,list_im_Ao1,list_kz_optAo1,omegac_opt,epsi1_imag_opt] = dataAo1
    
    dataBo1 = np.loadtxt(nameBo1,delimiter = '\t', skiprows=1)
    dataBo1 = np.transpose(dataBo1)
    [list_re_Bo1,list_im_Bo1,list_kz_optBo1,omegac_opt,epsi1_imag_opt] = dataBo1    
    
    if list_kz_optAo1[ind] != list_kz_optBo1[ind]:
        raise TypeError('Differents values of kz')
    
    kz = list_kz_optAo1[ind] #micrones
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = re_epsi1 + 1j*im_epsi1
    Ao1 = list_re_Ao1[ind] + 1j*list_im_Ao1[ind]
    Bo1 = list_re_Bo1[ind] + 1j*list_im_Bo1[ind]
    print('kz = ', kz)
    
    # Ao = coef_D2        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
    # Bo = coef_B2  
    
    # Ao = coef_C1        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
    # Bo = coef_A1  
    
    order1 = math.floor(math.log(np.abs(Ao1), 10))
    order2 = math.floor(math.log(np.abs(Bo1), 10))
    
    exp1,exp2 = 10**order1,10**order2
    
    Ao1 = Ao1/exp1       #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
    Bo1 = Bo1/exp2

    if kz < 0.13:
        path_save = path_save0 + '/' + 'kz_chico'
    else:
        path_save = path_save0 + '/' + 'kz_grande'
    
    if save_graphs==1:    
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)    

    if np.sign(Ao1.imag) == 1 and np.sign(Bo1.imag) == 1:
        labelAo1Bo1 = 'Ao = %.3f + i%.3f, Bo = %.3f + i%.3f' %(Ao1.real,Ao1.imag,Bo1.real,Bo1.imag)
    elif np.sign(Ao1.imag) == 1 and np.sign(Bo1.imag) == -1:
        labelAo1Bo1 = 'Ao = %.3f + i%.3f, Bo = %.3f - i%.3f' %(Ao1.real,Ao1.imag,Bo1.real,-Bo1.imag)
    elif np.sign(Ao1.imag) == -1 and np.sign(Bo1.imag) == 1:
        labelAo1Bo1 = 'Ao = %.3f - i%.3f, Bo = %.3f + i%.3f' %(Ao1.real,-Ao1.imag,Bo1.real,Bo1.imag)
    elif np.sign(Ao1.imag) == -1 and np.sign(Bo1.imag) == -1:
        labelAo1Bo1 = 'Ao = %.3f - i%.3f, Bo = %.3f - i%.3f' %(Ao1.real,-Ao1.imag,Bo1.real,-Bo1.imag)

    im_epsi1 = epsi1.imag

    info1 = 'z = %i $\mu$m, $\phi$ = %i, R = %.1f $\mu$m' %(z,phi,R)
    info2 = 'nmax = %i, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m' %(nmax,hbaramu,re_epsi1,-im_epsi1,omegac)

    title1 = labelAo1Bo1 + info1 + '\n' + info2  + '\n'  + name_this_py
    title2 =  'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + '\n'  + name_this_py
    return kz,epsi1,omegac,Ao1,Bo1,labelAo1Bo1,path_save,title1,title2

#%%


for modo in list_modos:
    
    kz,epsi1,omegac,Ao1,Bo1,labelAo1Bo1,path_save,title1,title2 = names(modo)

    title = title1
    
    #print('Graficar el campo |Etot|^2 para el medio 1 y 2')
    
    Etot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao1,Bo1,rho,phi,z)
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
    
    plt.figure(figsize=tamfig)
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
            plt.savefig('|Ez|_index%i_modo%i.png' %(ind,modo))
        else:
            plt.savefig('ReEz_index%i_modo%i.png' %(ind,modo))

if normal_pol == 1: 
    
    for modo in list_modos:   
        
        kz,epsi1,omegac,Ao1,Bo1,labelAo1Bo1,path_save,title1,title2 = names(modo)
    
        title = title2
        
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
        
        plt.figure(figsize=tamfig)
        plt.plot(list_rho,Etot_values,'.',ms=10,label ='modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
        plt.title(title,fontsize=tamtitle)
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
                plt.savefig('|Ez|_normalpol_index%i_modo%i.png' %(ind,modo))
            else:
                plt.savefig('ReEz_normalpol_index%i_modo%i.png' %(ind,modo))

#%%

print('')
print('Importar los valores de SPASER')

for modo in list_modos:
    
    kz,epsi1,omegac,Ao1,Bo1,labelAo1Bo1,path_save,title1,title2 = names(modo)

    title = title1
    
    #print('Graficar el campo |Etot|^2 para el medio 1 y 2')
    
    Htot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao1,Bo1,rho,phi,z)
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
                
    plt.figure(figsize=tamfig)
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
            plt.savefig('|Hz|_index%i_modo%i.png' %(ind,modo))
        else:
            plt.savefig('ReHz_index%i_modo%i.png' %(ind,modo))
    
if normal_pol == 1: 
    
    for modo in list_modos:   
        
        kz,epsi1,omegac,Ao1,Bo1,labelAo1Bo1,path_save,title1,title2 = names(modo)
    
        title = title2
    
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
        
        plt.figure(figsize=tamfig)
        plt.plot(list_rho,Htot_values,'.',ms=10,label = 'modo = %i, kz = %.4f 1/$\mu$m' %(modo,kz))
        plt.title(title,fontsize=tamtitle)
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
                plt.savefig('|Hz|_normalpol_index%i_modo%i.png' %(ind,modo))
            else:
                plt.savefig('ReHz_normalpol_index%i_modo%i.png' %(ind,modo))

#%%
