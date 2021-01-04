#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar la formula que obtuve
del Qscat (que aparece en el overleaf
           en la seccion 4.8 del cuaderno corto)

usar datos de frecuencia real

Qscat 1D y 2D

caso INVISIBILIDAD

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

save_graphs = 1 #guardar los graficos 2D/1D del campo
graph_2D = 1    #graficos 2D
graph_1D = 0   #graficos 1D

if graph_2D == 1:
    save_data = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'seccion_eficaz'
name_this_py = ' .Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)
        
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Qscat_sinkz import Qscat
except ModuleNotFoundError:
    print('Qscat_sinkz.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra Qscat_sinkz.py')
    sys.path.insert(1, path_basic)
    from Qscat_sinkz import Qscat

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
modo = 1
R = 0.5         #micrones 
nmax = 10        #sumatoria desde -nmax hasta +nmax (se suman 2*nmax + 1 modos)

zoom = 0

#%%

print('Importar los valores de SPASER')

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
[barrido_mu,omegac_opt,epsi1_imag_opt,eq_det] = data_load

m = len(barrido_mu)
index = 0
hbaramu = barrido_mu[index]
omegac = omegac_opt[index]
im_epsi1c = epsi1_imag_opt[index]
epsi1 = re_epsi1 + 1j*im_epsi1c
    
info1 = 'R = %.2f $\mu$m, nmax = %i, $\mu_c$ = %.4f eV, $\epsilon_1$ = %.1f - i%.5e' %(R,nmax,hbaramu,re_epsi1,-im_epsi1c)
info2 = ' y $\omega/c$ = %.5e 1/$\mu$m del modo = %i' %(omegac,modo)

labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely = 'Qscat$_{ad}$'
title = info1 +'\n' + info2 + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py

#%%

if R != 0.5:
    raise TypeError('Wrong value for R')
    
if modo in [2,3,4]:
    print('Ojo: Modo ' + modo + ' no excitado')
    
#%%

if graph_1D==1:
    print('Calcular Qscat para diferentes Im(epsi1)')
    
    tol = 1e-2
    list_im_epsi1_fino = [0,im_epsi1c + 3*tol,im_epsi1c + 2*tol,im_epsi1c + tol,im_epsi1c,im_epsi1c - tol,im_epsi1c - 2*tol]
    list_im_epsi1_grueso = [0,-im_epsi1c,0.5,-0.001,-0.01,im_epsi1c,-0.5]
    
    N = int(3*1e3)               
    omegac1,omegac2 = omegac*0.5,omegac*1.4
    if modo ==3 or modo ==4: 
        omegac1,omegac2 = omegac*0.9999,omegac*1.0001
    # lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
    list_omegac = np.linspace(omegac1,omegac2,N)
    
    list_Qscat_tot1 = []
    for im_epsi1 in list_im_epsi1_fino:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qscat1 = []
        for omeggac in list_omegac:
            epsi1 = re_epsi1 + 1j*im_epsi1
            Qscatt = Qscat(omeggac,epsi1,nmax,R,hbaramu,Ao)
            list_Qscat1.append(Qscatt)  
        list_Qscat_tot1.append(list_Qscat1)  
    
    del list_Qscat1
                
    list_Qscat_tot2 = []
    for im_epsi1 in list_im_epsi1_grueso:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qscat2 = []
        for omeggac in list_omegac:
            epsi1 = re_epsi1 + 1j*im_epsi1
            Qscatt = Qscat(omeggac,epsi1,nmax,R,hbaramu,Ao)
            list_Qscat2.append(Qscatt)  
        list_Qscat_tot2.append(list_Qscat2)  
        
    del list_Qscat2
    
    colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']
    
    
    print('Graficar Qscat para diferentes Im(epsi1)')
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = int(tamtitle*0.9))
    labelomegac = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac)

    
    for j in range(len(list_Qscat_tot1)):
        list_Qscat2 = list_Qscat_tot1[j]
        im_epsi1 = list_im_epsi1_fino[j]
            
        if im_epsi1 == im_epsi1c:
            labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
        elif im_epsi1 == 0:
            labell = 'Im($\epsilon_1$) = 0'  
        else:
            labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)
    
        plt.plot(list_omegac,list_Qscat2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qscat_tot1),np.max(list_Qscat_tot1)
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
        plt.savefig('Qscat_fino_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = int(tamtitle*0.9))
     
    for j in range(len(list_Qscat_tot2)):
        list_Qscat2 = list_Qscat_tot2[j]
        im_epsi1 = list_im_epsi1_grueso[j]
            
        if im_epsi1 == im_epsi1c:
            labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
        elif im_epsi1 == -im_epsi1c:
            labell = 'Im($\epsilon_1$) = -Im($\epsilon_1$)$_c$'    
        elif im_epsi1 == 0:
            labell = 'Im($\epsilon_1$) = 0'  
        else:
            labell = 'Im($\epsilon_1$) = %.3f'%(im_epsi1)
            
        plt.plot(list_omegac,list_Qscat2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qscat_tot2),np.max(list_Qscat_tot2)
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
        plt.savefig('Qscat_grueso_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  

#%%

if graph_2D==1:
    print('Graficar Qscat en 2D') 
    
    def Qscat2D(Omegac,Im_epsi1): 
        Epsi1 = re_epsi1 + 1j*Im_epsi1
        Qscatt = Qscat(Omegac,Epsi1,nmax,R,hbaramu,Ao)
        return Qscatt
       
    N = 400
    
    if zoom==1:
        tol = 1e-3
    else:
        tol = 1e-1
        # tol = 5*1e-1 #para ver el otro polo
    
    omegac12,omegac22 = omegac*(1-tol),omegac*(1+tol)
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    list_im_epsi1 = np.linspace(im_epsi1c - 5*tol,im_epsi1c + 5*tol,N)

    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    limits = [np.min(x) , np.max(x), np.min(y) , np.max(y)]
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel('Im($\epsilon_1$) ',fontsize=int(tamletra*1.2))
    plt.tick_params(labelsize = tamnum)
    # plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    vmin, vmax = np.min(Z), np.max(Z)
    maxlog = int(np.ceil( np.log10( np.abs(vmax) )))
    minlog = int(np.ceil( np.log10( np.abs(vmin) )))
    
    if np.abs(maxlog-minlog)>3:
    
        if vmin < 0:
            tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                              + [0] 
                              + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
        
        else:
            tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
        
        pcm = plt.pcolormesh(X, Y, Z,
                            norm = SymLogNorm(linthresh=0.03, linscale=0.03,
                                                vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    else:
        pcm = plt.pcolormesh(X, Y, Z,cmap='RdBu_r')
        
    plt.plot(x,np.ones(N)*im_epsi1c,'--',color = 'green')
    plt.plot(np.ones(N)*omegac,y,'--',color = 'green')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(labely,fontsize=tamletra)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save + '/2D')
        plt.savefig('Qscat2D_modo%i_mu%.4f.png' %(modo,hbaramu), format='png')  

    if save_data == 1:
        np.savetxt('info_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu), [inf_tot + ', se desprecian los campos incidentes en Qscat'],fmt='%s')
        np.savetxt('x_lambda_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),x, fmt='%1.5e', delimiter='\t', header = inf_tot)
        np.savetxt('y_im_epsi1_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),y, fmt='%1.5e', delimiter='\t', header = inf_tot)
        np.savetxt('Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),Z,fmt='%1.5e', delimiter='\t', header = inf_tot)
    
#%%
