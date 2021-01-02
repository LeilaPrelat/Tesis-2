#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

graficar epsilon1(omega) cerca de los valores criticos
de omega/c real y de epsilon_ci cerca de los dos polos que aparecieron
---> queremos ver si epsilon se anula en alguno de esos polos

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D/1D del campo
graph_2D = 1    #graficos 2D
graph_1D = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_epsilon1 = path_basic.replace('/' + 'real_freq','')

try:
    sys.path.insert(1, path_epsilon1)
    from epsilon1 import epsilon1
except ModuleNotFoundError:
    print('epsilon1.py no se encuentra  en ' + path_epsilon1)
    path_epsilon1 = input('path de la carpeta donde se encuentra epsilon1.py')
    sys.path.insert(1, path_epsilon1)
    from epsilon1 import epsilon1

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 15

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 1
hbaramu = 0.3005

Ep = 0.6
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

if R != 0.5:
    raise ValueError('El otro polo aparece para R = 0.5')

#%%

print('Importar los valores de SPASER')

path_save = path_basic + '/epsilon1/R_%.2f' %(R)
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

path_load = path_epsilon1 + '/seccion_eficaz/' + 'R_%.2f' %(R)
os.chdir(path_load)
name1 = 'minimus_omegac_modo%i_mu%.4f.txt' %(modo,hbaramu)
name2 = 'minimus_epsici_modo%i_mu%.4f.txt' %(modo,hbaramu)

try:
    data_load1 = np.loadtxt(name1,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name1) if l.startswith('#')):
        print('valores de ', name1, ':', line)
except OSError or IOError:
    print('El archivo ' + name1 + ' no se encuentra en ' + path_load)
    
try:
    data_load2 = np.loadtxt(name2,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name1) if l.startswith('#')):
        print('valores de ', name2, ':', line)
except OSError or IOError:
    print('El archivo ' + name2 + ' no se encuentra en ' + path_load)
    
list_omegac = np.transpose(data_load1)
list_im_epsi1 = np.transpose(data_load2)

N = 500
list_omegac = np.linspace(1.3,1.6,N)
list_im_epsi1 = np.linspace(0,-1.5,N)


info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV, Ep = %.1f eV' %(R,hbaramu,epsiinf_DL,gamma_DL,Ep)
info2 = 'modo = %i' %(modo)

labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely1 = 'Re($\epsilon_1(\omega)$)'
labely2 = 'Im($\epsilon_1(\omega)$)'
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py
title = info1 +'\n' + info2  +'\n' + name_this_py

#%%

# if modo in [2,3,4]:
#     print('Ojo: Modo ' + modo + ' no excitado')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

if graph_2D==1:
    print('Graficar Qabs en 2D') 
    
    def epsi1_2D(Omegac,Im_epsi1): 
        epsi1 = epsilon1(Omegac,Ep,epsiinf_DL,gamma_DL,Im_epsi1)
        return epsi1.real
    
    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(epsi1_2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel('$\epsilon_{ci}$',fontsize=int(tamletra*1.5))
    plt.tick_params(labelsize = tamnum)
    # plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')  
    
    pcm = plt.pcolormesh(X, Y, Z,cmap='RdBu_r')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label('|$\epsilon_1$|',fontsize=int(tamletra*1.5))
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save + '/2D')
        plt.savefig('epsilon1_modo%i_mu%.4f.png' %(modo,hbaramu), format='png') 


if graph_1D==1:
    n = 100
    list_omegac = np.linspace(1.3,1.6,n)
    ejey = []
    for omegac in list_omegac:
        
        epsi1 = epsilon1(omegac,Ep,epsiinf_DL,gamma_DL,0)
        ejey.append(epsi1.real)
        
    plt.figure(figsize=tamfig)
    plt.plot(list_omegac,ejey,lw = 1)
    plt.plot(list_omegac,np.ones(n)*0,'k-')
    plt.ylabel('Re($\epsilon_1$)',fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    

#%%
