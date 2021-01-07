#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar la formula que obtuve
del Qscat (que aparece en el overleaf
           en la seccion 4.8 del cuaderno corto)

Qscat 1D y 2D
"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import ast

save_graphs = 1 #guardar los graficos 2D/1D del campo
graph_2D = 0    #graficos 2D
graph_1D = 1    #graficos 1D

if graph_2D == 1:
    save_data = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
        
#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Qscat_nano import Qscat
except ModuleNotFoundError:
    print('Qscat_nano.py no se encuentra  en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra Qscat_nano.py')
    sys.path.insert(1, path_basic)
    from Qscat_nano import Qscat

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

R = 0.4              #micrones

Ep = 0.7
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

nmax = 10
Ao = 1

zoom = 0

path_save = path_basic + '/seccion_eficaz'  + '/' + 'R_%.2f' %(R)  + '/' + 'degenerations' + '/' + 'Ep_%.1f' %(Ep)
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
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

print('Importar los valores de SPASER')

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

index2 = 1
 # len(list_modes) == len(delta_ci)
delta_ci = list_epsi1_imag_tot[index1][index2]
modo = list_modes_tot[index1][index2]
list_modes = list_modes_tot[index1]
print('Modos degenerados:', list_modes_tot)
print('Degeneracion de los modos:',list_modes,'cerca del modo ', modo)
name_graph = '_modos%i%i_modo%i.png' %(list_modes[0],list_modes[1],modo) 

info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV, Ep = %.1f eV' %(R,hbaramu,epsiinf_DL,gamma_DL,Ep)
info2 = ', $\Delta_{ci}$ = %.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(delta_ci,omegac,modo,nmax,Ao)

labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely = 'Qscat$_{ad}$'
title = info1 +'\n' + info2  +'\n' + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py

#%%

if graph_1D==1:
    print('Calcular Qscat para diferentes Im(epsi1)')
    
    tol = 1e-3
    if modo == 4:
        tol = 1e-4
    list_im_epsi1_fino = [0,delta_ci + 3*tol,delta_ci + 2*tol,delta_ci + tol,delta_ci,delta_ci - tol,delta_ci - 2*tol]
    list_im_epsi1_grueso = [0,-delta_ci,0.5,-0.001,-0.01,delta_ci,-0.5]
    
    N = int(1e3)               
    omegac1,omegac2 = omegac*0.5,omegac*1.4
    list_omegac = np.linspace(omegac1,omegac2,N)
    
    list_Qscat_tot1 = []
    for im_epsi1 in list_im_epsi1_fino:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qscat1 = []
        for omeggac in list_omegac:
            Qscatt = Qscat(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
            list_Qscat1.append(Qscatt)  

        list_Qscat_tot1.append(list_Qscat1)  
        if im_epsi1 == 0:
            norm = np.max(list_Qscat1)
    
    del list_Qscat1
                
    list_Qscat_tot2 = []
    for im_epsi1 in list_im_epsi1_grueso:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qscat2 = []
        for omeggac in list_omegac:
            Qscatt = Qscat(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu,Ao)
            list_Qscat2.append(Qscatt)
        list_Qscat_tot2.append(list_Qscat2)  
        
    del list_Qscat2
    
    list_Qscat_tot2 = np.array(list_Qscat_tot2)/norm
    list_Qscat_tot1 = np.array(list_Qscat_tot1)/norm
    
    
    colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']
    
    
    print('Graficar Qscat para diferentes Im(epsi1)')
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = int(tamtitle*0.9))
    labelomegac = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac)
    
    for j in range(len(list_Qscat_tot1)):
        list_Qscat2 = list_Qscat_tot1[j]
        im_epsi1 = list_im_epsi1_fino[j]
            
        if im_epsi1 == delta_ci:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0'  
        else:
            labell = '$\epsilon_{ci}$ = %.7f'%(im_epsi1)
    
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
            
        if im_epsi1 == delta_ci:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == -delta_ci:
            labell = '$\epsilon_{ci}$ = -$\epsilon_{ci}$ crit'    
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0'  
        else:
            labell = '$\epsilon_{ci}$ = %.3f'%(im_epsi1)
            
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
        Qscatt = Qscat(Omegac,Ep,epsiinf_DL,gamma_DL,Im_epsi1,nmax,R,hbaramu,Ao)
        return Qscatt
       
    N = 450
    
    if zoom==1:
        tol = 1e-3
        list_im_epsi1 = np.linspace(delta_ci - 50*tol,delta_ci + 50*tol,N)
    else:
        if R == 0.5:
            #hay otro polo para R = 0.5
            tol = 0.9*1e-1
            list_im_epsi1 = np.linspace(delta_ci - 20*tol,delta_ci + 20*tol,N)
            
        else:
            #no hay otro polo para R = 0.05
            tol = 1e-2
            list_im_epsi1 = np.linspace(delta_ci - 30*tol,delta_ci + 30*tol,N)
        
    omegac12,omegac22 = omegac*(1-tol),omegac*(1+tol)
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5
    
    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    limits = [min(x) , max(x), min(y) , max(y)]
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel('$\epsilon_{ci}$',fontsize=int(tamletra*1.5))
    plt.tick_params(labelsize = tamnum)
    # plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    vmin, vmax = np.min(Z), np.max(Z)
    maxlog = int(np.ceil( np.log10( np.abs(vmax) )))
    minlog = int(np.ceil( np.log10( np.abs(vmin) )))
    
    if vmin < 0 :
        tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                          + [0] 
                          + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
    else:
        tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    

    
    pcm = plt.pcolormesh(X, Y, Z,
                        norm = SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=vmin, vmax=vmax),cmap='RdBu_r')
    
    plt.plot(x,np.ones(N)*delta_ci,'--',color = 'green')
    plt.plot(np.ones(N)*omegac,y,'--',color = 'green')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(labely,fontsize=tamletra)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save + '/2D')
        if zoom == 1:
            plt.savefig('Qscat2D_modo%i_mu%.4f_zoom.png' %(modo,hbaramu), format='png') 
        else:
            plt.savefig('Qscat2D_modo%i_mu%.4f.png' %(modo,hbaramu), format='png') 

    if save_data == 1:
        if zoom == 1:
            np.savetxt('info_Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu), [inf_tot + ', se desprecian los campos incidentes en Qscat'],fmt='%s')
            np.savetxt('x_lambda_Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu),x, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('y_im_epsi1_Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu),y, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu),Z,fmt='%1.5e', delimiter='\t', header = inf_tot)
        else:
            np.savetxt('info_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu), [inf_tot + ', se desprecian los campos incidentes en Qscat'],fmt='%s')
            np.savetxt('x_lambda_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),x, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('y_im_epsi1_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),y, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu),Z,fmt='%1.5e', delimiter='\t', header = inf_tot)             
    
#%%
