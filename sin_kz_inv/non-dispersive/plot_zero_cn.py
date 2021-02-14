#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

graficar el numerador del coeficiente cn (medio 1)
para ver que no se anule

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D/1D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/non-dispersive','') 
path_sinkz = path_basic.replace('/sin_kz_inv/non-dispersive','')
path_sinkz = path_sinkz + '/sin_kz/non-dispersive/real_freq'

path_save = path_basic + '/' + 'zero_cn'
name_this_py = ' .Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
        
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_graphene)
    from find_zero_cn import determinante
except ModuleNotFoundError:
    print('find_zero_cn.py no se encuentra en ' + path_graphene)
    path_basic = input('path de la carpeta donde se encuentra find_zero_cn.py')
    sys.path.insert(1, path_graphene)
    from find_zero_cn import determinante

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
Ao = 1          #pol p unicamente
modo = 4
R = 0.5         #micrones 
nmax = 10        #sumatoria desde -nmax hasta +nmax (se suman 2*nmax + 1 modos)

zoom = 0
index = 0       #value of chemical potential

#%%

print('Importar los valores de SPASER del caso sin kz INVISIBILIDAD')

path_load = path_basic + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)
os.chdir(path_load)
name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data_load = np.transpose(data_load)
[barrido_mu_inv,omegac_opt_inv,epsi1_imag_opt_inv,eq_det_inv] = data_load

m = len(barrido_mu_inv)
hbaramu_inv = barrido_mu_inv[index]
omegac_inv = omegac_opt_inv[index]
im_epsi1_inv = epsi1_imag_opt_inv[index]
epsi1_inv = re_epsi1 + 1j*im_epsi1_inv

del barrido_mu_inv,omegac_opt_inv,epsi1_imag_opt_inv,eq_det_inv

if zoom == 0:
    print('Importar los valores de SPASER del caso sin kz RESONANTE')
    
    path_sinkz = path_sinkz + '/re_epsi1_3.90_vs_mu'
    os.chdir(path_sinkz)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    
    try:
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_sinkz)
    
    data_load = np.transpose(data_load)
    [barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    
    hbaramu = barrido_mu[index]
    
    if hbaramu_inv != hbaramu:
        raise TypeError('el mu sin kz y el mu sin kz invisible no son iguales')
    
    omegac = omegac_opt[index]
    im_epsi1 = epsi1_imag_opt[index]
    epsi1 = re_epsi1 + 1j*im_epsi1
    
    if modo == 1:
        name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(2)
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        data_load = np.transpose(data_load)
        [barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load
            
        hbaramu2 = barrido_mu[index]
        if hbaramu2 != hbaramu :
            raise TypeError('mu del modo 1 y mu del modo 2 no coinciden')
        omegac_mod2 = omegac_opt[index]
        im_epsi12 = epsi1_imag_opt[index]


#%%

info1 = 'R = %.2f $\mu$m, nmax = %i, $\mu_c$ = %.4f eV, $\epsilon_1$ = %.1f - i%.5e' %(R,nmax,hbaramu_inv,re_epsi1,-im_epsi1_inv)
info2 = ' y $\omega/c$ = %.5e 1/$\mu$m del modo = %i' %(omegac_inv,modo)

labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely = '|numerador cn adimensional|'
title = info1 +'\n' + info2 + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py 

#%%

if R != 0.5:
    raise TypeError('Wrong value for R')
    
if modo in [2,3,4]:
    print('Ojo: Modo ', modo,' no excitado')
    
#%%

print('Calcular Qscat para diferentes Im(epsi1)')

tol = 1e-1
list_im_epsi1_fino = [0,im_epsi1_inv + 3*tol,im_epsi1_inv + 2*tol,im_epsi1_inv + tol,im_epsi1_inv,im_epsi1_inv - tol,im_epsi1_inv - 2*tol]
if zoom == 0:
    list_im_epsi1_fino.append(im_epsi1)

N = int(3*1e3)               
omegac1,omegac2 = omegac_inv*0.5,omegac_inv*1.4
list_omegac = np.linspace(omegac1,omegac2,N)

list_zero_cn_tot1 = []
for im_epsi1 in list_im_epsi1_fino:
    im_epsi1 = np.round(im_epsi1,7)
    print(im_epsi1)
    list_zero_cn = []
    for omeggac in list_omegac:
        epsi1 = re_epsi1 + 1j*im_epsi1
        zero = determinante(omeggac,epsi1,modo,R,hbaramu_inv)
        list_zero_cn.append(np.abs(zero))  
    list_zero_cn_tot1.append(list_zero_cn)  
            
colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']

print('Graficar Qscat para diferentes Im(epsi1)')

plt.figure(figsize=tamfig)
plt.title(title, fontsize = int(tamtitle*0.9))
labelomegac = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac_inv)
if zoom == 0:
    labelomegac2 = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac)
    if modo == 1:
        labelomegac3 = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac_mod2)


for j in range(len(list_zero_cn_tot1)):
    list_zero_cn = list_zero_cn_tot1[j]
    im_epsi1 = list_im_epsi1_fino[j]
        
    if im_epsi1 == im_epsi1_inv:
        labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
    elif im_epsi1 == 0:
        labell = 'Im($\epsilon_1$) = 0'  
    else:
        labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)

    plt.plot(list_omegac,list_zero_cn,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)

n = 10
mini,maxi = np.min(list_zero_cn_tot1),np.max(list_zero_cn_tot1)
eje_Lambda2 = np.linspace(mini,maxi,n)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))

if save_graphs==1:
    os.chdir(path_save)
    # plt.tight_layout(1)
    plt.savefig('zero_cn_modo%i_mu%.4f.png' %(modo,hbaramu_inv), format='png')  

#%%

