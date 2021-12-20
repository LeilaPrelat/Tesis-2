#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

campos transversales + longitudinales:

    con un dipolo interior en el cilindro
    
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 2D del campo
non_active_medium = 1 #plotear campos con im(epsilon1) = 0
paper = 0

list_field = ['Ez', 'Hz', 'Etot', 'Htot'] 
type_field = list_field[1]
if type_field == 'Ez' or type_field == 'Hz':
    modulo = 0 #if modulo == 1 graficar |Hz| o |Ez| y si vale 0 grafica la parte real

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
    tamfig = (5,5)
#    tamlegend = 18
#    tamletra = 18
#    tamtitle = 18
#    tamnum = 15
#    labelpady = 0
#    labelpadx = 0
#    pad = 0
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 7
    labelpady = -1.5
    labelpadx = -0.5
    pad = 0.5


#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save0 = path_basic + '/' + 'fields_DIP_3D'

if save_graphs==1:    
    if not os.path.exists(path_save0):
        print('Creating folder to save graphs')
        os.mkdir(path_save0)

#print('Importar modulos necesarios para este codigo')
err = 'fields_DIP_conkz.py no se encuentra en ' + path_basic
err2 = 'path de la carpeta donde se encuentra fields_DIP_conkz.py'
if type_field == 'Ez':
    try:
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Ez
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Ez

    def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        [rta1,rta2] = Ez(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
        if modulo == 1:
            [rta3,rta4] = [np.abs(rta1),np.abs(rta2)] 
        else:
            [rta3,rta4] = [rta1.real,rta2.real] 
        return [rta3,rta4]
    
elif type_field == 'Hz':
    try:
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Hz
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Hz

    def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        [rta1,rta2] = Hz(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
        if modulo == 1:
            [rta3,rta4] = [np.abs(rta1),np.abs(rta2)] 
        else:
            [rta3,rta4] = [rta1.real,rta2.real] 
        return [rta3,rta4]


elif type_field == 'Htot':
    try:
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Htot
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Htot

    def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        return Htot(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)

else: 

    try:
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Etot
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from fields_DIP_conkz import Etot

    def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        return Etot(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)

#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.05              # micrones
hbaramu = 0.3        # eV mu_c
list_modos = [1,2,3,4]
list_ind = [0,50,500,1000,5000,-1]
list_ind = [500, -1]
z = 0
nmax = 10
Ao,Bo = 1,1


### datos del dipolo ####
theta_D = 0
z_D = 0
rho_D = R/2

px = 1
py = 0
pz = 0
p1 = px + 1j*py # p+
p2 = px - 1j*py # p-

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

def circle(radio):
    listx = []
    listy = []
    for theta in np.linspace(0,2*np.pi):
        x = radio*np.cos(theta)
        y = radio*np.sin(theta)
        listx.append(x)
        listy.append(y)
    return listx,listy

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R not in [0.5,0.05]:
    raise TypeError('Wrong value for radium')

#%%

for modo in list_modos: 
    for ind in list_ind: 
        print('Importar los valores de SPASER')
        
        path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/R_%.2f/mu_%.1f' %(re_epsi1,R,hbaramu) 
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
        
        info1 = r'Ao = %i, Bo = %i, kz = %.4f1/$\mu$m, z = %i$\mu$m' %(Ao,Bo,kz,z)
        info2 = r'R = %.2f$\mu$m, nmax = %i, $\mu_c$ = %.1f eV, $\nu$ = %i' %(R,nmax,hbaramu,modo)
        info3 = r'$\epsilon_1$ = %.1f - i%.4e, $\omega/c$ = %.4e1/$\mu$m, $\rho_D$ = %.3f$\mu$m' %(re_epsi1,-im_epsi1,omegac,rho_D)
        infoDIP = r'px = %i, py = %i, pz = %i, $\theta_D$ = %i, $z_D$ = %i$\mu$m' %(px,py,pz,theta_D,z_D)
        title =  info1 + '\n' + info2 + '\n' + info3 + '\n' + infoDIP
        
        if non_active_medium == 1:
            info3_loss =  r'$\epsilon_1$ = %.1f, $\omega/c$ = %.4e 1/$\mu$m, $\rho_D$ = %.3f$\mu$m' %(re_epsi1,omegac,rho_D)
            title_loss =  info1 + '\n' + info2 + '\n' + info3_loss + '\n' + infoDIP
        
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
         
        listx1,listy1 = circle(rho_D)
        listx2,listy2 = circle(R)
        
        if non_active_medium == 1: #epsi1 = re_epsi1 ---> no hay medio activo
               
            def Etot_2variable(x,y):   
                phi = np.arctan2(y,x)
                rho = (x**2+y**2)**(1/2)
                [E1_tot,E2_tot] = fields(kz,omegac,re_epsi1,nmax,R,hbaramu,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
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
            
            plt.plot(listx1,listy1,'-',color = 'blue')
            plt.plot(listx2,listy2,'-',color = 'blue')
        
            if paper == 0:
                cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
            
            if save_graphs==1:
                plt.tight_layout(1)
                os.chdir(path_save2)
                plt.savefig(campo + '_loss_modo%i_kz%.4f.png' %(modo,kz), format='png')
     
        def Etot_2variable(x,y):   
            phi = np.arctan2(y,x)
            rho = (x**2+y**2)**(1/2)
            [E1_tot,E2_tot] = fields(kz,omegac,epsi1,nmax,R,hbaramu,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
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
        
        plt.plot(listx1,listy1,'-',color = 'blue')
        plt.plot(listx2,listy2,'-',color = 'blue')
        
        if paper == 0:
            cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
        
        if save_graphs==1:
            plt.tight_layout(1)
            os.chdir(path_save2)
            plt.savefig(campo + '_modo%i_kz%.4f.png' %(modo,kz), format='png')
        
        del Z1,Z2

        if paper == 1:
            np.savetxt('info_'+ campo + '_modo%i_kz%.4f.txt' %(modo,kz), [title + ', ' + name_this_py], fmt='%s')

#%%
