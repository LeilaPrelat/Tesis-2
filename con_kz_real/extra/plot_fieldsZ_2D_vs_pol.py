#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos longitudinales:
    
    seccion 3.2 del cuaderno corto

plotear el maximo de los campos (para los valores de spaser)
para un barrido en polarizacion (Ao, Bo)  y comparar con los 
autovectores de la matriz del caso homogeneo

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
path_save = path_basic + '/' + 'fields_2D_vs_pol'

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

try:
    sys.path.insert(1, path_basic)
    from coef_matrix import coef
except ModuleNotFoundError:
    print('coef_matrix.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra coef_matrix.py')
    sys.path.insert(1, path_basic)
    from coef_matrix import coef

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

modo = 1

z = 0
phi = 0
nmax = 10
ind = 50

Ao,Bo = 1,1

#%%

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

n = 500
list_rho = np.linspace(1e-5,2*R,n)

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')
    
#%%

print('')
print('Importar los valores de SPASER')


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

info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m,$\phi$ = %i, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,phi,R,nmax,hbaramu)
info2 = '$\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)

if kz < 0.13:
    path_save = path_save + '/' + 'kz_chico'
else:
    path_save = path_save + '/' + 'kz_grande'

if save_graphs==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#print('Graficar el campo |Etot|^2 para el medio 1 y 2')

coeff = coef(kz,omegac,epsi1,modo,R,hbaramu)
coef_A1,coef_B2,coef_C1,coef_D2 = coeff
coef_A1 = complex(coef_A1)
coef_C1 = complex(coef_C1)
coef_D2 = complex(coef_D2)
coef_B2 = complex(coef_B2)

Ao = coef_D2        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
Bo = coef_B2
n = 100
list_Bo = np.linspace(Bo*0.7,Bo*1.5,n)

Etot_values_max = []
for Bo in list_Bo:
    Etot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if np.abs(rho)<=R:
            Etot_values.append(np.abs(Ez1))
        else:
            Etot_values.append(np.abs(Ez2))
    Etot_values_max.append(np.max(Etot_values)/np.min(Etot_values))

#%%

plt.figure(figsize=tamfig)         
title = 'Ao = %.3f + i%.3f, ' %(Ao.real,Ao.imag) + info1 + '\n' + info2  + '\n' + name_this_py
plt.plot(list_Bo,Etot_values_max,'.',ms=10,label = 'modo = %i' %(modo))
plt.title(title,fontsize=tamtitle)
plt.ylabel('max|Ez|/min|E_z|',fontsize=tamletra)
plt.xlabel('Bo',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('max|Ez|_index%i' %(ind), format='png')

#%%

print('')
print('Importar los valores de SPASER')


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

info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m,$\phi$ = %i, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,phi,R,nmax,hbaramu)
info2 = '$\epsilon_1$ = %.1f - i%.2e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)

#print('Graficar el campo |Etot|^2 para el medio 1 y 2')

coeff = coef(kz,omegac,epsi1,modo,R,hbaramu)
coef_A1,coef_B2,coef_C1,coef_D2 = coeff
coef_A1 = complex(coef_A1)
coef_C1 = complex(coef_C1)
coef_D2 = complex(coef_D2)
coef_B2 = complex(coef_B2)

Ao = coef_D2        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
Bo = coef_B2   
n = 100
list_Ao = np.linspace(Ao*0.7,Ao*1.5,n)

Htot_values_max = []
for Ao in list_Ao:
    Htot_values = []
    for rho in list_rho:
        rho = np.abs(rho)
        [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if np.abs(rho)<=R:
            Htot_values.append(np.abs(Hz1))
        else:
            Htot_values.append(np.abs(Hz2))
    Htot_values_max.append(np.max(Htot_values)/np.min(Htot_values))

#%%

plt.figure(figsize=tamfig)          
title = 'Bo = %.3f + i%.3f, ' %(Bo.real,Bo.imag) + info1 + '\n' + info2  + '\n' + name_this_py
plt.plot(list_Ao,Htot_values_max,'.',ms=10,label = 'modo = %i' %(modo))
plt.title(title,fontsize=tamtitle)
plt.ylabel('max|Hz|/min|Hz|',fontsize=tamletra)
plt.xlabel('Ao',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('max|Hz|_index%i' %(ind), format='png')

#%%

coeff = coef(kz,omegac,epsi1,modo,R,hbaramu)
coef_A1,coef_B2,coef_C1,coef_D2 = coeff
coef_A1 = complex(coef_A1)
coef_C1 = complex(coef_C1)
coef_D2 = complex(coef_D2)
coef_B2 = complex(coef_B2)

Ao = coef_D2        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
Bo = coef_B2  

Ao = coef_C1        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
Bo = coef_A1  

Htot_values1 = []
Etot_values1 = []
for rho in list_rho:
    rho = np.abs(rho)
    [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    if np.abs(rho)<=R:
        Htot_values1.append(np.abs(Hz1))
        Etot_values1.append(np.abs(Ez1))
    else:
        Htot_values1.append(np.abs(Hz2))
        Etot_values1.append(np.abs(Ez2))

Ao,Bo = 1,1

Htot_values2 = []
Etot_values2 = []
for rho in list_rho:
    rho = np.abs(rho)
    [Ez1,Ez2,Hz1,Hz2] = fieldsZ(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
    if np.abs(rho)<=R:
        Htot_values2.append(np.abs(Hz1))
        Etot_values2.append(np.abs(Ez1))
    else:
        Htot_values2.append(np.abs(Hz2))
        Etot_values2.append(np.abs(Ez2))
        
#%%

labelAoBo = 'Ao = %.3f + i%.3f, Bo = %.3f + i%.3f' %(coef_D2.real,coef_D2.imag,coef_B2.real,coef_B2.imag)
labelx = r'$\rho$ [$\mu$m]'

plt.figure(figsize=tamfig)          
title = info1 + '\n' + info2  + '\n' + name_this_py
plt.plot(list_rho,Htot_values1,'.m',ms=10,label = labelAoBo)
plt.plot(list_rho,Htot_values2,'.b',ms=10,label = 'Ao = %i, Bo = %i' %(Ao,Bo))
plt.title(title,fontsize=tamtitle)
plt.ylabel('|Hz|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
# plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|Hz|_index%i' %(ind), format='png')

plt.figure(figsize=tamfig)          
title = info1 + '\n' + info2  + '\n' + name_this_py
plt.plot(list_rho,Etot_values1,'.m',ms=10,label = labelAoBo)
plt.plot(list_rho,Etot_values2,'.b',ms=10,label = 'Ao = %i, Bo = %i' %(Ao,Bo))
plt.title(title,fontsize=tamtitle)
plt.ylabel('|Ez|',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
# plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
plt.grid(1) 
if save_graphs==1:
    os.chdir(path_save)
    plt.savefig('|Ez|_index%i' %(ind), format='png')
    
    
#%%