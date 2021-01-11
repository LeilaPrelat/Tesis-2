#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

Antes de integrar en theta para obtener Qscat (despreciando campos incidentes) hecho SOLAMENTE para la polarizacion p
para el medio dispersivo (medio activo + nanocristal en el interior del cilindro)

la idea es ver la potencia es los casos degenerados 

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import ast

#%%

save_graphs = 1

non_active_medium = 1 #plotear campos con im(epsilon1) = 0
paper = 0  # sacar el titulo y guardar la info (formato paper)

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 15

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
name_this_py = 'Ver ' + name_this_py


try:
    sys.path.insert(1, path_basic)
    from potencia import potencia
except ModuleNotFoundError:
    print('Qscat_nano.py no se encuentra  en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra Qscat_nano.py')
    sys.path.insert(1, path_basic)
    from potencia import potencia


#%%

print('Definir parametros del problema')

Ao,Bo = 1,1
R = 0.35              #micrones

Ep = 0.8
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

nmax = 10

if save_graphs==1:
    path_save0 = 'potencia_disipada'  + '/' + 'barrido_mu'
    if paper == 0:
        path_save = path_basic + '/' + path_save0
    else:
        path_save = path_basic + '/' + path_save0 + '/' + 'paper'
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
        
#%%

deg_values = [[0.3,0.9],[0.35,0.8],[0.35,0.9],[0.4,0.7]]

if [R,Ep] not in deg_values:
    raise TypeError('R y Ep no son valores en los que hay degeneracion')

if epsiinf_DL != 3.9:
    raise TypeError('Wrong value for epsilon infinito')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Importar los valores de SPASER de deg')

path_load = path_basic + '/' + 'real_freq' + '/' + 'R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
os.chdir(path_load)
deg_txt = 'info_degenerations_dispR_%.2f.txt' %(R)
with open(deg_txt, 'r') as file: 
    list_mu = []
    list_omegac = []
    list_epsi1_imag_tot = []
    list_modes_tot = []
    
    for line in file:
        line.replace('\n','')
        mydict = ast.literal_eval(line)
        hbaramu,omegac = mydict['mu_c deg'],mydict['omega/c deg']
        list_mu.append(hbaramu)
        list_omegac.append(omegac)
        
        list_modes = []
        list_epsi1_imag = []
        for nu in mydict['modes']:
            list_modes.append(nu)
            delta_ci = mydict['im_epsi1 del modo %i' %(nu)]
            list_epsi1_imag.append(delta_ci)
            
        list_modes_tot.append(list_modes)
        list_epsi1_imag_tot.append(list_epsi1_imag)

index1 = 0 # len(barrido_mu) == len(omega_opt)
hbaramu = list_mu[index1]
omegac = list_omegac[index1]

index2 = 0
 # len(list_modes) == len(delta_ci)
delta_ci = list_epsi1_imag_tot[index1][index2]
modo = list_modes_tot[index1][index2]
list_modes = list_modes_tot[index1]
print('Modos degenerados:', list_modes_tot)
print('Degeneracion de los modos:',list_modes,'cerca del modo ', modo)
name_graph = '_modos%i%i_modo%i.png' %(list_modes[0],list_modes[1],modo) 

info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV, Ep = %.1f eV' %(R,hbaramu,epsiinf_DL,gamma_DL,Ep)
info2 = '$\Delta_{ci}$ = %.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(delta_ci,omegac,modo,nmax,Ao)
info3 = '$\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(omegac,modo,nmax,Ao)

labelx = '$\phi$'
labely = 'Re(potencia)'
title = info1 +'\n' + info2  +'\n' + name_this_py
title_loss = info1 +'\n' + info3  +'\n' + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py
name_graph = '_modos%i%i_modo%i.png' %(list_modes[0],list_modes[1],modo) 

#%%

print('Graficar 1D el campo potencia para el medio 1 y 2')

N = 200
list_phi = np.linspace(0,2*np.pi,N)
rho = R
list_rho = [R/4,R/2,R,1.5*R,2*R]
tick_label = ['R/4','R/2','R','1.5*R','2*R']
i = 0
if non_active_medium == 1: #epsi_ci = 0 ---> no hay medio activo
    
    plt.figure(figsize=tamfig)
    for rho in list_rho:
        list_pot = []
        for phi in list_phi:
            pot = potencia(omegac,Ep,epsiinf_DL,gamma_DL,0,nmax,hbaramu,Ao,rho,phi) 
            if modulo==1:
                pot = np.abs(pot)
            else:
                pot = pot.real
            list_pot.append(pot)
        
        plt.plot(list_phi,list_pot,'.',ms = 10,label = r'$\rho$ = %s' %(tick_label[i]))
        i = i + 1
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel(labely,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.legend(loc='best',markerscale=2, fontsize = tamlegend)
    if paper == 0:
        plt.title(title_loss,fontsize=int(tamtitle*0.9))
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modpot_loss' + name_graph, format='png') 
        else:
            plt.savefig('repot_loss' + name_graph, format='png') 

rho = R
plt.figure(figsize=tamfig)
list_pot = []
for phi in list_phi:
    pot = potencia(omegac,Ep,epsiinf_DL,gamma_DL,delta_ci,nmax,hbaramu,Ao,rho,phi) 
    if modulo==1:
        pot = np.abs(pot)
    else:
        pot = pot.real
    list_pot.append(pot)

plt.plot(list_phi,list_pot,'.',ms = 10)
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.yscale('log')
plt.tick_params(labelsize = tamnum)
if paper == 0:
    plt.title(title_loss,fontsize=int(tamtitle*0.9))
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modpot' + name_graph, format='png') 
    else:
        plt.savefig('repot' + name_graph, format='png') 

#%%

# def circle(radius):
#     listx = []
#     listy = []
#     listphi = np.linspace(0,2*np.pi,100)
#     for phi in listphi:
#         x = radius*np.cos(phi)
#         y = radius*np.sin(phi)
#         listx.append(x)
#         listy.append(y)
#     return listx,listy

# circlex,circley = circle(R)

# print('Graficar 2D el campo potencia para el medio 1 y 2')

# n2 = 200
# cota = 2*R
# x = np.linspace(-cota,cota,n2)
# y = np.linspace(-cota,cota,n2)
# X, Y = np.meshgrid(x, y)
# limits = [min(x) , max(x), min(y) , max(y)]

# if non_active_medium == 1: #epsi_ci = 0 ---> no hay medio activo
#     def Hz_2variable(x,y):   
#         phi = np.arctan2(y,x)
#         rho = (x**2+y**2)**(1/2)
        
#         pot = potencia(omegac,Ep,epsiinf_DL,gamma_DL,delta_ci,nmax,hbaramu,Ao,rho,phi) 
#         if modulo==1:
#             pot = np.abs(pot)
#         else:
#             pot = pot.real

#             return pot
        
#     f1 = np.vectorize(Hz_2variable)
#     Z = f1(X, Y)
#     maxH = np.max(Z) 
#     Z = Z/maxH
    
#     plt.figure(figsize=tamfig)
#     plt.xlabel(labelx,fontsize=tamletra)
#     plt.ylabel(labely,fontsize=tamletra)
#     plt.tick_params(labelsize = tamnum)
#     if paper == 0:
#         plt.title(title_loss,fontsize=int(tamtitle*0.9))
#     im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
#     cbar = plt.colorbar(im)
#     cbar.ax.tick_params(labelsize = tamnum)
#     cbar.set_label(label,size=tamlegend)
#     plt.plot(circlex,circley)
#     if save_graphs==1:
#         os.chdir(path_save)
#         if modulo==1:
#             plt.savefig('modpot_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png') 
#         else:
#             plt.savefig('repot_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  

# def Hz_2variable(x,y):   
#     phi = np.arctan2(y,x)
#     rho = (x**2+y**2)**(1/2)
#     [Hz1,Hz2] =potencia(omegac,Ep,epsiinf_DL,gamma_DL,delta_ci,nmax,hbaramu,Ao,rho,phi)   
#     if modulo==1:
#         Hz1_tot, Hz2_tot = np.abs(Hz1), np.abs(Hz2)
#     else:
#         Hz1_tot, Hz2_tot = Hz1.real, Hz2.real
#     if np.abs(rho)<=R: #medio1
#         return Hz1_tot
#     else: #medio2
#         return Hz2_tot
    
# f1 = np.vectorize(Hz_2variable)
# Z = f1(X, Y)
# if non_active_medium == 1: 
#     Z = Z/maxH
# plt.figure(figsize=tamfig)
# plt.xlabel(labelx,fontsize=tamletra)
# plt.ylabel(labely,fontsize=tamletra)
# plt.tick_params(labelsize = tamnum)
# if paper == 0:
#     plt.title(title,fontsize=int(tamtitle*0.9))
# im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
# cbar = plt.colorbar(im)
# cbar.ax.tick_params(labelsize = tamnum)
# cbar.set_label(label,size=tamlegend)
# plt.plot(circlex,circley)
# if save_graphs==1:
#     os.chdir(path_save)
#     if modulo==1:
#         plt.savefig('modpot_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
#     else:
#         plt.savefig('repot_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  

if paper == 1:
    np.savetxt('info_potencia_modo%i_mu%.4f.txt' %(modo,hbaramu), [inf_tot],fmt='%s')   

#%%
