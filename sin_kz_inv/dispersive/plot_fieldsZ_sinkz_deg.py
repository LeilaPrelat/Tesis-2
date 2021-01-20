#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

obtuve gn (minimizandolo) ---> obtuve
los denominadores de los coef an y cn
(entonces obtuve an y cn) 
---> podemos graficar los campos
    en el caso donde hay degeneracion
    (misma frecuencia de spaser para diferentes modos, peero
     diferente valor im(epsi1) critico)

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import ast

#%%

save_graphs = 1
modulo = 1 #si modulo == 1 ---> |Hz| (si modulo == 0 ---> Re(Hz))

non_active_medium = 1 #plotear campos con im(epsilon1) = 0
paper = 0  # sacar el titulo y guardar la info (formato paper)

plot_degenerations = 1 #encontrar degeneraciones en frecuencia: modos con el mismo omega/c

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 15

if modulo==1:
    label = '$|H_z|$'
else:
    label = 'Re($H_z$)'
    
if modulo==1:
    label2 = '$|E_z|$'
else:
    label2 ='Re($E_z$)'

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
name_this_py = 'Ver ' + name_this_py


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
R = 0.1              #micrones

Ep = 0.9
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

nmax = 10

if save_graphs==1:
    path_save0 = 'fields'  + '/' + 'degenerations' + '/' + 'R_%.2f' %(R) + '/' + 'Ep_%.1f' %(Ep)
    if paper == 0:
        path_save = path_basic + '/' + path_save0
    else:
        path_save = path_basic + '/' + path_save0 + '/' + 'paper'
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
        
#%%

deg_values =  [[0.1, 3.9, 0.9], [0.3, 3.9, 0.9], [0.4, 3.9, 0.9], [0.5, 3.9, 0.9]] 


if [R,epsiinf_DL,Ep] not in deg_values:
    raise TypeError('R, epsiinf_DL, Ep no son valores en los que hay degeneracion')

if epsiinf_DL != 3.9:
    raise TypeError('Wrong value for epsilon infinito')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + 'R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
os.chdir(path_load)
deg_txt1 = 'info_degenerations_dispR_%.2f.txt' %(R)
deg_txt2 = 'info_degenerations_inv_dispR_%.2f.txt' %(R)

def open_txt(path):
    
    with open(path, 'r') as file: 
    
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
            for [nu,nu2] in mydict['modes']:
                list_modes.append([nu,nu2])
                delta_ci = mydict['im_epsi1 del modo %i con %i' %(nu,nu2)]
                list_epsi1_imag.append(delta_ci)
                
            list_modes_tot.append(list_modes)
            list_epsi1_imag_tot.append(list_epsi1_imag)
        
    return list_mu,list_omegac,list_epsi1_imag_tot,list_modes_tot

try:
    list_mu,list_omegac,list_epsi1_imag_tot,list_modes_tot = open_txt(deg_txt1)
    print('deg_txt1 not found')
except OSError or IOError:
    list_mu,list_omegac,list_epsi1_imag_tot,list_modes_tot = open_txt(deg_txt2)
    
index1 = 1 #0 ---> saltos de a 1 modo, 1 ---> saltos de a 2 modos
index2 = 1
index3 = 1 # 1 o 0 tiene que ser

try: 
    hbaramu = list_mu[index1][index2]
    omegac = list_omegac[index1][index2]
    delta_ci = list_epsi1_imag_tot[index1][index2]
    list_modes = list_modes_tot[index1][index2]
    modo = list_modes_tot[index1][index2][index3] #---> solo si hay 2 diccionarios
    
except IndexError: #tenemos solo un diccionario en el deg_txt
    index2 = 0 
    hbaramu = list_mu[index1][index2]
    omegac = list_omegac[index1][index2]
    delta_ci = list_epsi1_imag_tot[index1][index2]
    list_modes = list_modes_tot[index1]    
    modo = list_modes_tot[index1][index2]

print('Degeneracion de los modos:',list_modes,'cerca del modo ', modo)
name_graph = '_modos%i%i_modo%i.png' %(list_modes[0],list_modes[1],modo) 

info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, Ep = %.1f eV, $\gamma_{DL}$ = %.2f eV' %(R,hbaramu,epsiinf_DL,Ep,gamma_DL)
info2 = '$\Delta_{ci}$ = %.5e del modo = %i y $\omega/c$ = %.5e 1/$\mu$m , nmax = %i' %(delta_ci,modo,omegac,nmax)
title1 = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2  + '\n' + name_this_py
title2 = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2  + '\n' + name_this_py

labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
infotot = 'Ao = %i, Bo = %i' %(Ao,Bo) + ', ' + info1 +  ', ' + info2 +  ', ' + name_this_py

if non_active_medium == 1:
    info2_loss = '$\omega/c$ = %.2e 1/$\mu$m del modo = %i, nmax = %i' %(omegac,modo,nmax)
    title1_loss = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2_loss  + '\n' + name_this_py
    title2_loss = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py

#%%

def circle(radius):
    listx = []
    listy = []
    listphi = np.linspace(0,2*np.pi,100)
    for phi in listphi:
        x = radius*np.cos(phi)
        y = radius*np.sin(phi)
        listx.append(x)
        listy.append(y)
    return listx,listy

circlex,circley = circle(R)
    
print('Graficar el campo Hz para el medio 1 y 2')

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
    plt.plot(circlex,circley)
    plt.tick_params(labelsize = tamnum)
    if paper == 0:
        plt.title(title1_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(label,size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modHz_loss' + name_graph, format='png') 
        else:
            plt.savefig('reHz_loss' + name_graph, format='png') 

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
plt.plot(circlex,circley)
plt.tick_params(labelsize = tamnum)
if paper == 0:
    plt.title(title1,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(label,size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modHz' + name_graph, format='png') 
    else:
        plt.savefig('reHz' + name_graph, format='png') 

#%%

print('Graficar el campo Ez para el medio 1 y 2')

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
    plt.plot(circlex,circley)
    plt.tick_params(labelsize = tamnum)
    if paper == 0:
        plt.title(title2_loss,fontsize=int(tamtitle*0.9))
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(label2,size=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modEz_loss' + name_graph, format='png') 
        else:
            plt.savefig('reEz_loss' + name_graph, format='png') 

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
if non_active_medium == 1: 
    Z = Z/maxE
# Z = Z/np.max(Z)

plt.figure(figsize=tamfig)
plt.xlabel(labelx,fontsize=tamletra)
plt.ylabel(labely,fontsize=tamletra)
plt.plot(circlex,circley)
plt.tick_params(labelsize = tamnum)
if paper == 0:
    plt.title(title2,fontsize=int(tamtitle*0.9))
im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
cbar = plt.colorbar(im)
cbar.ax.tick_params(labelsize = tamnum)
cbar.set_label(label2,size=tamlegend)
if save_graphs==1:
    os.chdir(path_save)
    if modulo==1:
        plt.savefig('modEz' + name_graph, format='png') 
    else:
        plt.savefig('reEz' + name_graph, format='png') 

if paper == 1:
    np.savetxt('info_fields_disp_modo%i_mu%.4f.txt' %(modo,hbaramu), [infotot],fmt='%s')   

#%%
