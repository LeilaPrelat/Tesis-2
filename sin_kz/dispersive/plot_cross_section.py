#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

graficar Qscat Qabs o Qext (elegir en la variable cross_section)

ESTABAN MAL ESAS FORMULAS DEL CUADERNO (movidas a extra)
VER TESIS
(revisar cuentas)

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 
graph_2D = 1    #graficos 2D
zoom = 1        #graficos 2D con o sin zoom
graph_1D = 0    #graficos 1D

if graph_2D == 1:
    paper = 1

list_cross_section = ['Qscat', 'Qabs', 'Qext'] 
cross_section = list_cross_section[2]
    
#%%

print('Definir parametros del problema')

R = 0.05              #micrones
modo = 1

Ep = 0.6
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

nmax = 5
Ao = 1
index = 0

#%% 

#print('importar funciones para graficarlas')

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

if save_graphs==1:
    path_g = path_basic + '/' + 'seccion_eficaz/R_%.2f' %(R)
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)

#print('Importar modulos necesarios para este codigo')
err = 'cross_sections_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo'
err2 = 'path de la carpeta donde se encuentra cross_sections_conkz.py'
if cross_section == 'Qscat':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections import Qscat
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections import Qscat
        
    def Cross_Section(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao):
        return Qscat(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)

elif cross_section == 'Qabs':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections import Qabs
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections import Qabs
        
    def Cross_Section(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao):
        return Qabs(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)

else: 

    try:
        sys.path.insert(1, path_basic)
        from cross_sections import Qext
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections import Qext
        
    def Cross_Section(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao):
        return Qext(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)

#%%

#print('Definir parametros para graficos')

if paper == 1: 
    tamfig = (3.5,3.5)
    tamlegend = 7
    tamletra = 6
    tamtitle = 6
    tamnum = 6
    labelpady = -2.5
    labelpadx = 0.7
    pad = 0.25
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

print('Importar los valores de SPASER')

path_load = path_basic + '/' + 'real_freq' + '/' + 'R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
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

m = len(barrido_mu)
hbaramu = barrido_mu[index]
omegac0 = omegac_opt[index]
crit = epsi1_imag_opt[index]
    
info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV, Ep = %.1f eV' %(R,hbaramu,epsiinf_DL,gamma_DL,Ep)
info2 = ', $\Delta_{ci}$ = %.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(crit,omegac0,modo,nmax,Ao)

title = info1 +'\n' + info2 + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py

del barrido_mu,omegac_opt,epsi1_imag_opt,eq_det

labelomegac = '$\omega/c$ = %.2f$\mu m^{-1}$' %(omegac0)
labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely2 = '$\epsilon_{ci}$'
labely1 = cross_section + '$_{ad}$'
name = cross_section
inf_fig = '_modo%i_mu%.4f.png' %(modo,hbaramu)
if zoom == 1: 
    inf_fig = '_modo%i_mu%.4f_zoom.png' %(modo,hbaramu)
    
#%%
    
if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')
    
if modo in [2,3,4]:
    print('Ojo: Modo ', modo, ' no excitado')
    
#%%

if graph_1D==1:
    print('')
    print('Calcular '+ name + ' para diferentes Im(epsi1)')
    
    tol = 1e-3
    list_im_epsi1_fino = [0,crit+3*tol,crit+2*tol,crit+tol,crit,crit-tol,crit-2*tol]
    list_im_epsi1_grueso = [0,-crit,0.5,-0.001,-0.01,crit,-0.5]
    
    N = int(5*1e3)               
    omegac1,omegac2 = omegac0*0.98,omegac0*1.02
    # lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
    list_omegac = np.linspace(omegac1,omegac2,N)
    
    list_Qabs_tot1 = []
    for im_epsi1 in list_im_epsi1_fino:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qabs1 = []
        for omeggac in list_omegac:
            Qabss = Cross_Section(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
            list_Qabs1.append(Qabss)  
        list_Qabs_tot1.append(list_Qabs1)  
    
    del list_Qabs1
                
    list_Qabs_tot2 = []
    for im_epsi1 in list_im_epsi1_grueso:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qabs2 = []
        for omeggac in list_omegac:
            Qabss = Cross_Section(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
            list_Qabs2.append(Qabss)  
        list_Qabs_tot2.append(list_Qabs2)  
    
    del list_Qabs2
    
    colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    
    print('Graficar ' + name  + ' para diferentes Im(epsi1)')
    
    plt.figure(figsize=tamfig)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    
    for j in range(len(list_Qabs_tot1)):
        list_Qabs2 = np.abs(list_Qabs_tot1[j])
        im_epsi1 = list_im_epsi1_fino[j]
            
        if im_epsi1 == crit:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0' 
        else:
            labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)
    
        plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    
    n = 10
    mini,maxi = np.min(list_Qabs_tot1),np.max(list_Qabs_tot1)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    plt.ylabel(labely1,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    
    if save_graphs==1:
        os.chdir(path_g)
        # plt.tight_layout(1)
        plt.savefig(name + '_fino' + inf_fig, format='png')  
        if paper == 0:
            np.savetxt('info_' + name + '_modo%i_mu%.4f_zoom1D.txt' %(modo,hbaramu), [inf_tot],fmt='%s')
    
    plt.figure(figsize=tamfig)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
     
    for j in range(len(list_Qabs_tot2)):
        list_Qabs2 = np.abs(list_Qabs_tot2[j])
        im_epsi1 = list_im_epsi1_grueso[j]
            
        if im_epsi1 == crit:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == -crit:
            labell = '$\epsilon_{ci}$ = -$\epsilon_{ci}$ crit'    
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0'  
        else:
            labell = '$\epsilon_{ci}$ = %.3f'%(im_epsi1)
            
        plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qabs_tot2),np.max(list_Qabs_tot2)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    plt.ylabel(labely1,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    if save_graphs==1:
        os.chdir(path_g)
        # plt.tight_layout(1)
        plt.savefig(name + '_grueso' + inf_fig, format='png') 

#%%

from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors

if graph_2D==1:
    print('Graficar '+ name + ' en 2D') 
    
    def Qscat2D(Omegac,Im_epsi1): 
        Qscatt = Cross_Section(Omegac,Ep,epsiinf_DL,gamma_DL,Im_epsi1,nmax,R,hbaramu,Ao)
#        return np.log10(Qscatt)
        return Qscatt
       
    N = 500
    
    if zoom==1:
        tol = 1e-3
    else:
        tol = 1e-2
        # tol = 5*1e-1 #para ver el otro polo
    
    omegac12,omegac22 = omegac0*(1-tol),omegac0*(1+tol)
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    list_im_epsi1 = np.linspace(crit - 5*tol,crit + 5*tol,N)

    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    limits = [np.min(x) , np.max(x), np.min(y) , np.max(y)]
    plt.ylabel(labely2,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    vmin,vmax = np.min(Z), np.max(Z)
    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
    
    
    # if np.abs(maxlog-minlog)>3:
    
    if vmin < 0:
        tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                          + [0] 
                          + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
    
    else:
        tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
    
    pcm = plt.pcolormesh(X, Y, Z,
                        norm = SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    # else:
    #     pcm = plt.pcolormesh(X, Y, Z,cmap='RdBu_r')
    
    plt.plot(x,np.ones(N)*crit,'--',lw = 0.5,color = 'green')
    plt.plot(np.ones(N)*omegac0,y,'--',lw = 0.5,color = 'green')
    
#    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    
    cbar = plt.colorbar(pcm, extend='both')
    cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize = tamnum)

    if save_graphs==1:
        os.chdir(path_g + '/2D')
        plt.savefig(name + '2D' + inf_fig, format='png') 
        if paper == 1:
            np.savetxt('info_' + name + '_modo%i_mu%.4f_zoom2D.txt' %(modo,hbaramu), [inf_tot],fmt='%s')

#%%
