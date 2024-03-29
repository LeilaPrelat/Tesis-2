#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

obtuve gn (minimizandolo) ---> obtuve
los denominadores de los coef an y cn
(entonces obtuve an y cn) 
---> podemos graficar los campos
        en la condicion de spaser

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt

#%%

save_graphs = 1
modulo = 1      #si modulo == 1 ---> |Hz| (si modulo == 0 ---> Re(Hz))

non_active_medium = 1    #plotear campos con im(epsilon1) = 0
paper = 1   # sacar el titulo y guardar la info (formato paper)

#%%

print('Definir parametros para graficos')

if paper == 1: 
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    labelpady = -1.5
    labelpadx = -0.5
    pad = 0.5
else:
    tamfig = (10,8)
    tamlegend = 18
    tamletra = 18
    tamtitle = 18
    tamnum = 15
    labelpady = 0   

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

if save_graphs==1:
    if paper == 0:
        path_save = path_basic + '/' + 'fields'
    else:
        path_save = path_basic + '/' + 'fields' + '/' + 'paper'
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

list_modos = [1,2,3,4]
R = 0.5      #micrones 
nmax = 10
index = 0

#%%

if R != 0.5:
    raise TypeError('Wrong value for R')

labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
n2 = 250
cota = 2*R
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    # if paper == 1:
    #     plt.tight_layout(1)
    return 

#%%

for modo in list_modos:

    print('Importar los valores de SPASER')
    
    path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu' %(re_epsi1)
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
    
    data_load = np.transpose(data_load)
    [barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    hbaramu = barrido_mu[index]
    omegac = omegac_opt[index]
    im_epsi1 = epsi1_imag_opt[index]
    epsi1 = re_epsi1 + 1j*im_epsi1
    
    info1 = 'R = %.2f $\mu$m, nmax = %i, $\mu_c$ = %.4f eV' %(R,nmax,hbaramu)
    info2 = '$\epsilon_1$ = %.1f - i%.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)
    title1 = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2  + '\n' + name_this_py
    title2 = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2  + '\n' + name_this_py
    
    infotot = 'Ao = %i, Bo = %i' %(Ao,Bo) + ', ' + info1 +  ', ' + info2 +  ', ' + name_this_py
    
    if non_active_medium == 1:
        info2_loss = '$\epsilon_1$ = %.1f, $\omega/c$ = %.5e 1/$\mu$m del modo = %i' %(re_epsi1,omegac,modo)
        title1_loss = 'Ao = %i, ' %(Ao) + info1 + '\n' + info2_loss  + '\n' + name_this_py
        title2_loss = 'Bo = %i, ' %(Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py
    
    print('Graficar el campo Hz para el medio 1 y 2')


    limits = [min(x) , max(x), min(y) , max(y)]
    X, Y = np.meshgrid(x, y)
    
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
            
        f1 = np.vectorize(Hz_2variable)
        Z = f1(X, Y)
        maxH1 = np.max(Z) 
        Z = Z/maxH1
        
        graph(title1_loss,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        if paper == 0:
            cbar.set_label(label,fontsize=tamlegend,labelpad = 1)
        if save_graphs==1:
            os.chdir(path_save)
            if modulo==1:
                plt.savefig('modHz_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
            else:
                plt.savefig('reHz_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
    
        del Z
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
    
    f1 = np.vectorize(Hz_2variable)
    Z = f1(X, Y)
    if non_active_medium == 1 and paper == 0 :
        Z = Z/maxH1
    if paper == 1: #normarlizar TODOS los campos : con y sin medio activo
        maxH2 = np.max(Z) 
        Z = Z/maxH2
    
    graph(title1,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if paper == 0:
        cbar.set_label(label,fontsize=tamlegend,labelpad = 1)
    if save_graphs==1:  
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modHz_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
        else:
            plt.savefig('reHz_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
    del Z
  
    print('Graficar el campo Ez para el medio 1 y 2')
        
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
            
        X, Y = np.meshgrid(x, y)
        f2 = np.vectorize(Ez_2variable)
        
        Z = f2(X, Y)
        maxE1 = np.max(Z) 
        Z = Z/maxE1
        
        graph(title2_loss,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        if paper == 0:
            cbar.set_label(label2,fontsize=tamlegend,labelpad = 1)

        if save_graphs == 1:
            if modulo==1:
                plt.savefig('modEz_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
            else:
                plt.savefig('reEz_loss_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
        del Z
        
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
        
    f2 = np.vectorize(Ez_2variable)
    Z = f2(X, Y)
    if non_active_medium == 1 and paper == 0:
        Z = Z/maxE1
    elif paper == 1: #normarlizar TODOS los campos : con y sin medio activo
        maxE2 = np.max(Z) 
        Z = Z/maxE2
        
    graph(title2,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    if paper == 0:
        cbar.set_label(label2,fontsize=tamlegend,labelpad = 1)
    if save_graphs==1:
        os.chdir(path_save)
        if modulo==1:
            plt.savefig('modEz_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
        else:
            plt.savefig('reEz_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
                
    if paper == 1:
        np.savetxt('info_fields_nondisp_modo%i_mu%.4f.txt' %(modo,hbaramu), [infotot],fmt='%s')        
    del Z

#%%
