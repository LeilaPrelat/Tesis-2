#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

electromagnetic energy density of the mode definido en el paper
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
path_save = path_basic + '/' + 'energy_density'
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from energy_density import energy_density
except ModuleNotFoundError:
    print('energy_density.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra energy_density.py')
    sys.path.insert(1, path_basic2)
    from energy_density import energy_density

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
R = 0.5         #micrones

list_modos = [1,2,3,4]

z = 0
phi = 0
nmax = 8
index = 0

Ao,Bo = 1,1

#%%

n = 500
list_rho = np.linspace(-2*R,2*R,n)

list_kz = [0,0.05,0.1,0.13,0.135,0.14,0.5,0.75,1]

labely = '$\omega_{em}$ energy density'
labelx = r'$\rho$ [$\mu$m]'
info1 = 'z = %i, $\phi$ = %i, nmax = %i, R = %.1f $\mu$m' %(z,phi,nmax,R)

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

path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1) 

#%%

if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('')
print('Importar los valores de SPASER')

for modo in list_modos:
    info2 = 'Re($\epsilon_1$) = %.2f, modo = %i' %(re_epsi1,modo)
    title0 = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + ', ' + info2 
    
    for kz in list_kz: 
        
        os.chdir(path_load)
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
        
        hbaramu = list_mu_opt[index] 
        crit = epsi1_imag_opt[index]
        omegac = omegac_opt[index] 
        epsi1 = re_epsi1 + 1j*crit
        info3 = 'kz = %.4f 1/$\mu$m, $\mu$ = %.4f eV, Im($\epsilon_1$) = %.3f, $\omega/c$ = %.3f 1/$\mu$m' %(kz,hbaramu,crit,omegac)
        title = title0 + '\n' + info3  + '\n' + name_this_py    
        
        omega_em_tot1 = []
        for rho in list_rho:
            omega_em1 = energy_density(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            omega_em_tot1.append(omega_em1)
        
        plt.figure(figsize=tamfig)
        plt.plot(list_rho,omega_em_tot1,'.m',ms=10)
        plt.title(title,fontsize=tamtitle)
        plt.ylabel(labely,fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        # plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
        plt.grid(1) 
        if save_graphs==1:
            os.chdir(path_save)
            plt.savefig('energy_density_modo%i_kz%.4f' %(modo,kz), format='png')
      
        info3 = 'kz = %.4f 1/$\mu$m, $\mu$ = %.3f eV, $\omega/c$ = %.3f 1/$\mu$m' %(kz,hbaramu,omegac)
        title = title0 + '\n' + info3  + '\n' + name_this_py    
        
        omega_em_tot2 = []
        for rho in list_rho:
            omega_em2 = energy_density(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            omega_em_tot2.append(omega_em2)
        
        plt.figure(figsize=tamfig)
        plt.plot(list_rho,omega_em_tot2,'.b',ms=10,label = 'Im($\epsilon_1$) = 0')
        
        plt.title(title,fontsize=tamtitle)
        plt.ylabel(labely,fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='best',markerscale=2,fontsize = tamlegend)
        plt.grid(1) 
        if save_graphs==1:
            os.chdir(path_save)
            plt.savefig('energy_density_loss_modo%i_kz%.4f' %(modo,kz), format='png')        

#%%
