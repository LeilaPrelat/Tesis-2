#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar la formula que obtuve
del Qabs (que aparece en el overleaf
           en la seccion 4.8)

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo
graph_2D = 0    #graficos 2D

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'seccion_eficaz'
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
        
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Qabs_nano import Qabs
except ModuleNotFoundError:
    print('Qabs_nano.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra Qabs_nano.py')
    sys.path.insert(1, path_basic)
    from Qabs_nano import Qabs

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 3
Ep = 0.3
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia
nmax = 10
Ao = 1

#%%

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + r'/epsiinf_DL_%.2f_vs_mu' %(epsiinf_DL)
os.chdir(path_load)
name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)

try:
    data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('valores de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en el path_load')
    try:
        os.chdir(path_load + '/' + 'Ep_%.1f' %(Ep))
        data_load = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError as error:
        print(error)
    
data_load = np.transpose(data_load)
[barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load

index = 0
hbaramu = barrido_mu[index]
omegac = omegac_opt[index]
delta_ci = epsi1_imag_opt[index]
    
info1 = 'R = %.1f $\mu$m, $\mu_c$ = %.1f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(R,hbaramu,epsiinf_DL,gamma_DL)
info2 = '$\Delta_{ci}$ = %.3e y $\omega/c$ = %.2e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(delta_ci,omegac,modo,nmax,Ao)
title = info1 +'\n' + info2  +'\n' + name_this_py

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('')
print('Calcular Qabs para diferentes Im(epsi1)')

tol = 1e-3
list_im_epsi1_fino = [0,delta_ci + 3*tol,delta_ci + 2*tol,delta_ci + tol,delta_ci,delta_ci - tol,delta_ci - 2*tol]
list_im_epsi1_grueso = [0,-delta_ci,0.5,-0.001,-0.01,delta_ci,-0.5]

N = int(1e3)               
omegac1,omegac2 = omegac*0.97,omegac*1.03
if modo ==3 or modo ==4: 
    omegac1,omegac2 = omegac*0.9999,omegac*1.0001
# lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
list_omegac = np.linspace(omegac1,omegac2,N)

list_Qabs_tot1 = []
for im_epsi1 in list_im_epsi1_fino:
    im_epsi1 = np.round(im_epsi1,7)
    print(im_epsi1)
    list_Qabs1 = []
    for omeggac in list_omegac:
        Qabss = Qabs(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
        list_Qabs1.append(Qabss)  
    list_Qabs_tot1.append(list_Qabs1)  
    

del list_Qabs1
            
list_Qabs_tot2 = []
for im_epsi1 in list_im_epsi1_grueso:
    im_epsi1 = np.round(im_epsi1,7)
    print(im_epsi1)
    list_Qabs2 = []
    for omeggac in list_omegac:
        Qabss = Qabs(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
        list_Qabs2.append(Qabss)  
    list_Qabs_tot2.append(list_Qabs2)  
    
del list_Qabs2

#%%

colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']


print('Graficar Qabs para diferentes Im(epsi1)')

plt.figure(figsize=tamfig)

plt.title(title, fontsize = int(tamtitle*0.9))
labelomegac = '$\omega/c$ = %.2f$\mu m^{-1}$' %(omegac)
labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely = '|Qabs$_{ad}|$'
 
for j in range(len(list_Qabs_tot1)):
    list_Qabs2 = np.abs(list_Qabs_tot1[j])
    im_epsi1 = list_im_epsi1_fino[j]
        
    if im_epsi1 == delta_ci:
        labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
    elif im_epsi1 == 0:
        labell = '$\epsilon_{ci}$ = 0'  
    else:
        labell = '$\epsilon_{ci}$ = %.7f' %(im_epsi1)

    plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)


n = 10
mini,maxi = np.min(np.abs(list_Qabs_tot1)),np.max(np.abs(list_Qabs_tot1))
eje_Lambda2 = np.linspace(mini,maxi,n)
plt.plot(omegac*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))

if save_graphs==1:
    os.chdir(path_save)
    # plt.tight_layout(1)
    plt.savefig('Qabs_fino_modo%i' %(modo)) 

#%%

plt.figure(figsize=tamfig)
plt.title(title, fontsize = int(tamtitle*0.9))
 
for j in range(len(list_Qabs_tot2)):
    list_Qabs2 = np.abs(list_Qabs_tot2[j])
    im_epsi1 = list_im_epsi1_grueso[j]
        
    if im_epsi1 == delta_ci:
        labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
    elif im_epsi1 == -delta_ci:
        labell = '$\epsilon_{ci}$ = -$\epsilon_{ci}$ crit'    
    elif im_epsi1 == 0:
        labell = '$\epsilon_{ci}$ = 0'  
    else:
        labell = '$\epsilon_{ci}$ = %.3f'%(im_epsi1)
        
    plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)

n = 10
mini,maxi = np.min(np.abs(list_Qabs_tot2)),np.max(np.abs(list_Qabs_tot2))
eje_Lambda2 = np.linspace(mini,maxi,n)
plt.plot(omegac*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
if save_graphs==1:
    os.chdir(path_save)
    # plt.tight_layout(1)
    plt.savefig('Qabs_grueso_modo%i' %(modo)) 

#%%

from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors

zoom = 0

if graph_2D==1:
    print('Graficar Qabs en 2D') 
    
    def Qscat2D(Lambbda,Im_epsi1): 
        Qscatt = Qabs(Lambbda,Im_epsi1)   
        return Qscatt
       
    N = 300
    
    if zoom==1:
        tol = 1e-3
    else:
        tol = 5*1e-1
    
    omegac12,omegac22 = omegac*(1-tol),omegac*(1+tol)
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5
    list_im_epsi1 = np.linspace(delta_ci - delta,delta_ci + delta,N)
    
    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=tamletra)
    plt.ylabel('$\epsilon_{ci}$',fontsize=tamletra)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    pcm = plt.pcolormesh(X, Y, Z,
                        norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=-1.0, vmax=1.0),cmap='RdBu_r')
    
    plt.plot(x,np.ones(N)*delta_ci,'k-',label= '$\epsilon_{ci}$ crit')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.set_label(labely,fontsize=int(tamletra*0.8))
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig('Qabs2D_modo%i' %(modo)) 
    
#%%
