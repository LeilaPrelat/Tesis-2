# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values.py: graficar optical gain y omega/c obtenidos formato paper

dif con critical_values_allmodes: aca se grafica para cada modo (un grafico por modo)

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

#%%

save_graphs = 1
paper = 1

if paper == 1: 
    tamfig = (6,3.5)
    tamlegend = 10
    tamletra = 10
    tamtitle = 10
    tamnum = 9
    
    labelpady = -0.45
    labelpadx = -0.5
    lw = 1.5
    pad = -2
    loc1 = [0.21,0.88]
    loc2 = [0.24,1]
else:
    tamfig = (11,7)
    tamlegend = 18
    tamletra = 20
    tamtitle = 20
    tamnum = 16
    labelpady = 0
    labelpadx = 2
    lw = 3.5
    pad = 0
    loc1 = [0.15,0.88]
    loc2 = [0.025,1]

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/non-dispersive/real_freq','') 


try:
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en ' + path_basic)
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
 
path_load = path_basic + '/re_epsi1_%.2f_vs_mu' %(re_epsi1)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py
labelx = '$\mu_c$ [eV]'
labely = '[Im($\epsilon_1$)]$_c$'
label1 = 'Rigurosa'
label2 = 'Aproximada'

os.chdir(path_load)
np.savetxt('info_critical_values_nondisp.txt', [inf_tot], fmt='%s')

colors = ['darkred','steelblue','coral','yellowgreen']
symbols = ['-','--','-.','-']

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')
    
#%% 

labely = '[Im($\epsilon_1$)]$_c$'

k = 0
for modo in [1,2,3,4]:
    plt.figure(figsize=tamfig)
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    

    im_epsi1_QE = []
    omegac_QE = []
    for mu in list_mu_opt:
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        b = im_epsi1_cuasi(a,modo,R,mu) 
        im_epsi1_QE.append(b)
        
    plt.plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw = lw,color = colors[k],label = label1)
    plt.plot(list_mu_opt,im_epsi1_QE,symbols[k+1],lw = lw,color = colors[k+1],label = label2) 

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 2,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('Im_epsi1_vs_mu%i' %(modo))    

#%%

labely = '$\omega/c$ $[\mu m]^{-1}$'


k = 0
for modo in [1,2,3,4]:
    plt.figure(figsize=tamfig)
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    

    im_epsi1_QE = []
    omegac_QE = []
    for mu in list_mu_opt:
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        omegac_QE.append(a)
        
    plt.plot(list_mu_opt,omegac_opt,symbols[k],lw = lw,color = colors[k],label = label1)
    plt.plot(list_mu_opt,omegac_QE,symbols[k+1],lw = lw,color = colors[k+1],label = label2 ) 

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 2,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('omegac_vs_mu%i' %(modo))  

#%%

labely = '$\omega$ [THz]'
k = 0
for modo in [1,2,3,4]:
    plt.figure(figsize=tamfig)
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    

    im_epsi1_QE = []
    omega_QE = []
    omega_opt = []
    for j in range(len(list_mu_opt)):
        mu = list_mu_opt[j]
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        omega_QE.append(a*c*1e-12)
        omega_opt.append(omegac_opt[j]*c*1e-12)
        
    plt.plot(list_mu_opt,omega_opt,symbols[k],lw = lw,color = colors[k],label = label1)
    plt.plot(list_mu_opt,omega_QE,symbols[k+1],lw = lw,color = colors[k+1],label = label2 ) 

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 2,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('omega_vs_mu%i' %(modo))  

#%%

labely = 'Energía [eV]'


k = 0
for modo in [1,2,3,4]:
    plt.figure(figsize=tamfig)
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    

    im_epsi1_QE = []
    energy_QE = []
    for mu in list_mu_opt:
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        energy_QE.append(a*c*hb)
        
    plt.plot(list_mu_opt,np.array(omegac_opt)*hb*c,symbols[k],lw = lw,color = colors[k],label = label1)
    plt.plot(list_mu_opt,energy_QE,symbols[k+1],lw = lw,color = colors[k+1],label =label2) 

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 2,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('energy_vs_mu%i' %(modo))  

#%%


