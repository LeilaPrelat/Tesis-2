#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar la formula que obtuve
del Qscat (que aparece en el overleaf
           en la seccion 4.8 del cuaderno corto)

Qscat 1D y 2D para los casos donde hay degeneracion
"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import ast

save_graphs = 1 #guardar los graficos 2D/1D del campo
graph_2D = 1    #graficos 2D
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

R = 0.1              #micrones

Ep = 0.9
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

nmax = 10
Ao = 1

zoom = 0
path_save = path_basic + '/seccion_eficaz/degenerations/R_%.2f' %(R)
name_this_py = 'Ver ' + name_this_py

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#%%

deg_values =  [[0.1, 3.9, 0.9], [0.3, 3.9, 0.9], [0.4, 3.9, 0.9], [0.5, 3.9, 0.9]] 

if [R,epsiinf_DL,Ep] not in deg_values:
    raise TypeError('R, epsiinf_DL, Ep no son valores en los que hay degeneracion')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Importar los valores de SPASER')

def open_txt(path,index1,index2,index3):
    
    with open(path, 'r') as file: 
    
        list_mu = []
        list_omegac = []
        list_epsi1_imag_tot = []
        list_modes_tot = []
        
        for line in file:
            line.replace('\n','')
            mydict = ast.literal_eval(line)
            hbaramu,omegac = mydict['mu_c deg'], mydict['omega/c deg']
            list_mu.append(hbaramu)
            list_omegac.append(omegac)
            
            list_modes = mydict['modes']
            list_modes_tot.append(list_modes)
            list_epsi1_imag = []
            for [nu1,nu2] in list_modes:
                delta_ci1 = mydict['im_epsi1 del modo %i con %i' %(nu1,nu2)]
                delta_ci2 = mydict['im_epsi1 del modo %i con %i' %(nu2,nu1)]
                list_epsi1_imag.append([delta_ci1,delta_ci2])
            list_epsi1_imag_tot.append(list_epsi1_imag)    
        
    try: 
        hbaramu = list_mu[index1][index2]
        omegac = list_omegac[index1][index2]
        list_modes = list_modes_tot[index1][index2]
        modo = list_modes_tot[index1][index2][index3]
        delta_ci = list_epsi1_imag_tot[index1][index2][index3]
        
        #---> solo si hay 2 diccionarios
        
    except IndexError: #tenemos solo un diccionario en el deg_txt
        index2 = 0 
        hbaramu = list_mu[index1][index2]
        omegac = list_omegac[index1][index2]
        list_modes = list_modes_tot[index1]    
        modo = list_modes_tot[index1][index2]
        delta_ci = list_epsi1_imag_tot[index1][index2]
        
    return hbaramu,omegac,delta_ci,list_modes,modo

#%%

path_load = path_basic + '/' + 'real_freq' + '/' + 'R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
os.chdir(path_load)
deg_txt1 = 'info_degenerations_dispR_%.2f.txt' %(R)
deg_txt2 = 'info_degenerations_inv_dispR_%.2f.txt' %(R)

index1 = 0 #0 ---> saltos de a 1 modo, 1 ---> saltos de a 2 modos
index2 = 3
index3 = 0 # 1 o 0 tiene que ser (longitud de list_modes)

os.chdir(path_load)
try:
    np.loadtxt(deg_txt1)
    hbaramu_inv,omegac_inv,delta_ci_inv,list_modes,modo = open_txt(deg_txt1,index1,index2,index3)
    print('deg_txt1 not found')
except OSError or IOError:
    hbaramu_inv,omegac_inv,delta_ci_inv,list_modes,modo = open_txt(deg_txt2,index1,index2,index3)

print('Degeneracion de los modos:',list_modes,'cerca del modo ', modo)
print('Im(epsilon1): ', delta_ci_inv, 'mu:', hbaramu_inv, 'omega/c:', omegac_inv)
name_graph = '_modos%i%i_modo%i.png' %(list_modes[0],list_modes[1],modo) 

info1 = 'R = %.2f $\mu$m, $\mu_c$ = %.4f eV, $\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV, Ep = %.1f eV' %(R,hbaramu_inv,epsiinf_DL,gamma_DL,Ep)
info2 = '$\Delta_{ci}$ = %.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i, nmax = %i, Ao = %i' %(delta_ci_inv,omegac_inv,modo,nmax,Ao)

labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely = 'Qscat$_{ad}$'
title = info1 +'\n' + info2  +'\n' + name_this_py
inf_tot = info1 + ', ' + info2  + ', ' + name_this_py

#%%

if graph_1D==1:
    print('Calcular Qscat para diferentes Im(epsi1)')
    
    tol = 1e-3
    list_im_epsi1_fino = [0,delta_ci_inv + 3*tol,delta_ci_inv + 2*tol,delta_ci_inv + tol,delta_ci_inv,delta_ci_inv - tol,delta_ci_inv - 2*tol]
    list_im_epsi1_grueso = [0,-delta_ci_inv,0.5,-0.001,-0.01,delta_ci_inv,-0.5]
    # list_im_epsi1_fino.append(delta_ci)
    # list_im_epsi1_grueso.append(delta_ci)
    N = int(1e3)               
    omegac1,omegac2 = omegac_inv*0.5,omegac_inv*1.4
    list_omegac = np.linspace(omegac1,omegac2,N)
    
    list_Qscat_tot1 = []
    for im_epsi1 in list_im_epsi1_fino:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qscat1 = []
        for omeggac in list_omegac:
            Qscatt = Qscat(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu_inv,Ao)
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
            Qscatt = Qscat(omeggac,Ep,epsiinf_DL,gamma_DL,im_epsi1,nmax,R,hbaramu_inv,Ao)
            list_Qscat2.append(Qscatt)
        list_Qscat_tot2.append(list_Qscat2)  
        
    del list_Qscat2

    list_Qscat_tot2 = np.array(list_Qscat_tot2)/norm
    list_Qscat_tot1 = np.array(list_Qscat_tot1)/norm    
    
    colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']
    
    
    print('Graficar Qscat para diferentes Im(epsi1)')
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = int(tamtitle*0.9))
    labelomegac = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac_inv)
    # if zoom == 0:
    #     labelomegac2 = '$\omega/c$ = %.4f$\mu m^{-1}$' %(omegac)
    
    for j in range(len(list_Qscat_tot1)):
        list_Qscat2 = list_Qscat_tot1[j]
        im_epsi1 = list_im_epsi1_fino[j]
            
        if im_epsi1 == delta_ci_inv:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0'  
        else:
            labell = '$\epsilon_{ci}$ = %.7f'%(im_epsi1)
    
        plt.plot(list_omegac,list_Qscat2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qscat_tot1),np.max(list_Qscat_tot1)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac_inv*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    # if zoom == 0:
    #     plt.plot(omegac*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac2)
    plt.ylabel(labely,fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    
    if save_graphs==1:
        os.chdir(path_save)
        # plt.tight_layout(1)
        plt.savefig('Qscat_fino' + name_graph, format='png')  
    
    plt.figure(figsize=tamfig)
    plt.title(title, fontsize = int(tamtitle*0.9))
     
    for j in range(len(list_Qscat_tot2)):
        list_Qscat2 = list_Qscat_tot2[j]
        im_epsi1 = list_im_epsi1_grueso[j]
            
        if im_epsi1 == delta_ci_inv:
            labell = '$\epsilon_{ci}$ = $\epsilon_{ci}$ crit'
        elif im_epsi1 == -delta_ci_inv:
            labell = '$\epsilon_{ci}$ = -$\epsilon_{ci}$ crit'    
        elif im_epsi1 == 0:
            labell = '$\epsilon_{ci}$ = 0'  
        else:
            labell = '$\epsilon_{ci}$ = %.3f'%(im_epsi1)
            
        plt.plot(list_omegac,list_Qscat2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qscat_tot2),np.max(list_Qscat_tot2)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac_inv*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    # if zoom == 0:
    #     plt.plot(omegac*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac2)
    plt.ylabel(labely,fontsize=tamletra)
    plt.xlabel(labelx,fontsize=tamletra)
    plt.tick_params(labelsize = tamnum)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    if save_graphs==1:
        os.chdir(path_save)
        # plt.tight_layout(1)
        plt.savefig('Qscat_grueso' + name_graph, format='png')  

#%%

if graph_2D==1:
    print('Graficar Qscat en 2D') 
    
    def Qscat2D(Omegac,Im_epsi1): 
        Qscatt = Qscat(Omegac,Ep,epsiinf_DL,gamma_DL,Im_epsi1,nmax,R,hbaramu_inv,Ao)
        return Qscatt
     
    N = 400
    if zoom == 1:
        tol1 = 1e-1
        tol2 = 0.5
        omegac12, omegac22 = omegac_inv*(1-tol1), omegac_inv*(1+tol1)
        im_epsi1_1, im_epsi1_2 = delta_ci_inv - tol2, delta_ci_inv + tol2
    else:
        # omegac0 = np.mean([omegac_inv,omegac])
        # im_epsi10 = np.mean([delta_ci_inv,delta_ci])
        
        tol1 = 0.3
        tol2 = 0.8
        omegac12,omegac22 = omegac_inv*(1-tol1), omegac_inv*(1+tol1)
        im_epsi1_1, im_epsi1_2 = delta_ci_inv - tol2, delta_ci_inv + tol2
        # tol = 5*1e-1 #para ver el otro polo
    
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    list_im_epsi1 = np.linspace(im_epsi1_1,im_epsi1_2,N)

    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    plt.figure(figsize=tamfig)
    plt.title(title,fontsize=int(tamtitle*0.9))
    limits = [np.min(x) , np.max(x), np.min(y) , np.max(y)]
    plt.xlabel(labelx,fontsize=int(tamletra*1.2))
    plt.ylabel('Im($\epsilon_1$) ',fontsize=int(tamletra*1.2))
    plt.tick_params(labelsize = tamnum)
    # plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    vmin, vmax = np.min(Z), np.max(Z)
    maxlog = int(np.ceil( np.log10( np.abs(vmax) )))
    minlog = int(np.ceil( np.log10( np.abs(vmin) )))
    
    if np.abs(maxlog-minlog)>5:
    
        if vmin < 0:
            tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                              + [0] 
                              + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+2)] )
        
        else:
            tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
        
        pcm = plt.pcolormesh(X, Y, Z,
                            norm = SymLogNorm(linthresh=vmin*1e-3,#desde donde hasta donde la escala es lineal (no diverge en el 0)
                                                vmin=vmin, vmax=vmax),cmap='RdBu_r')
    else:
        pcm = plt.pcolormesh(X, Y, Z, norm = SymLogNorm(linthresh=vmin*1e-3,#desde donde hasta donde la escala es lineal (no diverge en el 0)
                                                vmin=vmin, vmax=vmax),cmap='RdBu_r')
        
    plt.plot(x,np.ones(N)*delta_ci_inv,'--',color = 'green')
    plt.plot(np.ones(N)*omegac_inv,y,'--',color = 'green')
    cbar = plt.colorbar(pcm, extend='both')
    # if zoom == 0:
    #     plt.plot(x,np.ones(N)*delta_ci,'--',color = 'magenta')
    #     plt.plot(np.ones(N)*omegac,y,'--',color = 'magenta')
    cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize = tamnum)
    cbar.set_label(labely,fontsize=tamletra)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save + '/2D')
        if zoom == 1:
            plt.savefig('Qscat2D'  + name_graph, format='png')  
        elif zoom == 0:
            plt.savefig('Qscat2D_zoom'  + name_graph, format='png')  

    if save_data == 1:
        np.savetxt('info_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu_inv), [inf_tot + ', se desprecian los campos incidentes en Qscat'],fmt='%s')
        if zoom == 1:
            np.savetxt('x_lambda_Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu_inv),x, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('y_im_epsi1_Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu_inv),y, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('Qscat_modo%i_mu%.4f_zoom.txt' %(modo,hbaramu_inv),Z,fmt='%1.5e', delimiter='\t', header = inf_tot)
        else:
            np.savetxt('x_lambda_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu_inv),x, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('y_im_epsi1_Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu_inv),y, fmt='%1.5e', delimiter='\t', header = inf_tot)
            np.savetxt('Qscat_modo%i_mu%.4f.txt' %(modo,hbaramu_inv),Z,fmt='%1.5e', delimiter='\t', header = inf_tot)   
            
#%%
