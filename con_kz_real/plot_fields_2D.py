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

#%%

save_graphs = 1 #guardar los graficos 2D del campo
non_active_medium = 1 #plotear campos con im(epsilon1) = 0
paper = 1

list_field = ['Ez', 'Hz', 'Etot', 'Htot'] 
type_field = list_field[3]
if type_field == 'Ez' or type_field == 'Hz':
    modulo = 1 #if modulo == 1 graficar |Hz| o |Ez| y si vale 0 grafica la parte real

#%%

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
    labelpadx = 0
    pad = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save0 = path_basic + '/' + 'fields_3D'

if save_graphs==1:    
    if not os.path.exists(path_save0):
        print('Creating folder to save graphs')
        os.mkdir(path_save0)

#print('Importar modulos necesarios para este codigo')
err = 'fields_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo'
err2 = 'path de la carpeta donde se encuentra fields_conkz.py'
if type_field == 'Ez':
    try:
        sys.path.insert(1, path_basic)
        from fields_conkz import Ez
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_conkz import Ez

    def fields(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z):
        rta = Ez(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if modulo == 1:
            rta2 = np.abs(rta)
        else:
            rta2 = rta.real
        return rta2
    
elif type_field == 'Hz':
    try:
        sys.path.insert(1, path_basic)
        from fields_conkz import Hz
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_conkz import Hz

    def fields(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z):
        rta = Hz(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
        if modulo == 1:
            rta2 = np.abs(rta)
        else:
            rta2 = rta.real
        return rta2


elif type_field == 'Htot':
    try:
        sys.path.insert(1, path_basic)
        from fields_conkz import Htot
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_conkz import Htot

    def fields(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z):
        return Htot(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)

else: 

    try:
        sys.path.insert(1, path_basic)
        from fields_conkz import Etot
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_conkz import Etot

    def fields(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z):
        return Etot(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)

#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5              # micrones
hbaramu = 0.3        # eV mu_c
list_modos = [1]
list_ind = [0,50,500,1000,5000,-1]
list_ind = [976]
z = 0
nmax = 10
Ao,Bo = 1,1

#%%
    
n2 = 100
cota = 2*R
x = np.linspace(-cota,cota,n2)
y = np.linspace(-cota,cota,n2)
X, Y = np.meshgrid(x, y)
limits = [min(x) , max(x), min(y) , max(y)]

labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
if type_field == 'Htot':
    labelz = '|Htot|$^2$'
    campo = 'Htot'
elif type_field == 'Etot':
    labelz = '|Etot|$^2$'
    campo = 'Etot'
else: 
    if modulo == 1:
        labelz = '|' + type_field + '|'
        campo = labelz
    else:
        labelz = 'Re(' + type_field + ')'
        campo  =  'Re' + type_field 
#    name = type_field


def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    if paper == 1:
        plt.tight_layout()
    return   

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

for modo in list_modos: 
    for ind in list_ind: 
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
    
        ind = int(ind)
        kz = list_kz_opt[ind] #micrones
        print('kz = ', kz)
        print('modo = ', modo)
        print('')
        
        im_epsi1 = epsi1_imag_opt[ind]
        omegac = omegac_opt[ind] 
        epsi1 = re_epsi1 + 1j*im_epsi1
        
        # ktot = omegac*((mu1*epsi1)**(1/2))
        
        # ang_theta = np.arcsin(kz/np.abs(ktot))
        # ang2 = np.cos(np.pi/2 - ang_theta) 
        # print(360*ang2/(2*np.pi))
        
        
        info1 = 'kz = %.4f 1/$\mu$m, z = %i $\mu$m, R = %.1f $\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,z,R,nmax,hbaramu)
        info2 = '$\epsilon_1$ = %.1f - i%.4e y $\omega/c$ = %.4e 1/$\mu$m del modo = %i' %(re_epsi1,-im_epsi1,omegac,modo)
        title = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2  + '\n' + name_this_py
        
        if non_active_medium == 1:
            info2_loss = '$\epsilon_1$ = %.1f, $\omega/c$ = %.2e 1/$\mu$m del modo = %i' %(re_epsi1,omegac,modo)
            title_loss = 'Ao = %i, Bo = %i, ' %(Ao,Bo) + info1 + '\n' + info2_loss  + '\n' + name_this_py
        
        if Ao*Bo != 0:
            path_save1 = path_save0 + '/' + '2pol'
        elif Bo == 0:
            path_save1 = path_save0 + '/' + 'polAo'
        elif Ao == 0:
            path_save1 = path_save0 + '/' + 'polBo'
            
        if save_graphs==1:    
            if not os.path.exists(path_save1):
                print('Creating folder to save graphs')
                os.mkdir(path_save1)
            
        if kz < 0.13:
            path_save2 = path_save1 + '/' + 'kz_chico'
        else:
            path_save2 = path_save1 + '/' + 'kz_grande'
        
        if save_graphs==1:    
            if not os.path.exists(path_save2):
                print('Creating folder to save graphs')
                os.mkdir(path_save2)
    
  
        print('Graficar el campo ' + labelz + ' para el medio 1 y 2')
          
        if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
               
            def Etot_2variable(x,y):   
                phi = np.arctan2(y,x)
                rho = (x**2+y**2)**(1/2)
                [E1_tot,E2_tot] = fields(kz,omegac,re_epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
                if np.abs(rho) <= R: #medio1
                    return E1_tot
                else: #medio2
                    return E2_tot
            
            f2 = np.vectorize(Etot_2variable)
            Z2 = f2(X, Y)
            maxi = np.max(Z2)
            Z2 = Z2/maxi
            
            graph(title_loss,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
            im = plt.imshow(Z2, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
            cbar = plt.colorbar(im)
            cbar.ax.tick_params(labelsize = tamnum)
            if paper == 0:
                cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
            
            if save_graphs==1:
                plt.tight_layout(1)
                os.chdir(path_save2)
                plt.savefig(campo + '_loss_modo%i_kz%.4f.png' %(modo,kz), format='png')
        
        def Etot_2variable(x,y):   
            phi = np.arctan2(y,x)
            rho = (x**2+y**2)**(1/2)
            [E1_tot,E2_tot] = fields(kz,omegac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            if np.abs(rho) <= R: #medio1
                return E1_tot
            else: #medio2
                return E2_tot
            
        
        f1 = np.vectorize(Etot_2variable)
        Z1 = f1(X, Y)
        if non_active_medium == 1 and paper == 0:
            Z1 = Z1/maxi
        elif paper == 1: #normarlizar TODOS los campos : con y sin medio activo
            maxi2 = np.max(Z1) 
            Z1 = Z1/maxi2
        
            
        graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        im = plt.imshow(Z1, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        if paper == 0:
            cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
        
        if save_graphs==1:
            plt.tight_layout(1)
            os.chdir(path_save2)
            plt.savefig(campo + '_modo%i_kz%.4f.png' %(modo,kz), format='png')
        
        del Z1,Z2
        if paper == 1:
            np.savetxt('info_'+ campo + '_modo%i_kz%.4f.txt' %(modo,kz), [title], fmt='%s')

#%%
