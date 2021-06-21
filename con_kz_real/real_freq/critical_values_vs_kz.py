# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values.py: graficar optical gain y omega/c obtenidos formato paper

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import seaborn as sns

save_graphs = 0
paper = 1

#%%

if paper == 1: 
    tamfig = (5.5,3.5)
    tamlegend = 7
    tamletra = 6
    tamtitle = 6
    tamnum = 6
    
    labelpady = 0.5
    labelpadx = 0.5
    lw = 1.5
    loc1 = [0.03,1]
    loc2 = [0.03,1]
    
else:
    tamfig = (11,7)
    tamlegend = 18
    tamletra = 20
    tamtitle = 20
    tamnum = 16
    
    labelpady = 0
    labelpadx = 0
    lw = 3.5
    loc1 = [-0.1,1]
    loc2 = [-0.1,1]

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/con_kz_real/real_freq','') 


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
list_kz = [0,0.0001,0.05,0.1,0.135,0.14,0.5] 
list_kz_chicos = [0,0.05,0.1,0.13]

path_load = path_basic + '/re_epsi1_%.2f_vs_mu' %(re_epsi1)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_conkz.txt', [inf_tot], fmt='%s')

colors = ['darkred','yellowgreen','steelblue','coral']
symbols = ['-','-','--','-.']
sns.set()

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

#%%


for modo in [1,2,3,4]: 

    plt.figure(figsize=tamfig)

    labelx = '$\mu_c$ [eV]'
    labely = '[Im($\epsilon_1$)]$_c$'
    
    k = 0
    for kz in list_kz_chicos:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        line1, = plt.plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw = lw,color = colors[k],label = 'kz %.2f$\mu m^{-1}$'%(kz))
        if k == 0:
            line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
        k = k + 1
        
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum)
                
    #fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('Im_epsi1_vs_mu_tot_modo%i.png'%(modo) ,format = 'png')
        
    labely = '$\omega/c$'
    
    plt.figure(figsize=tamfig)
    
    k = 0
    for kz in list_kz_chicos:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        line1, = plt.plot(list_mu_opt,omegac_opt,symbols[k],lw = lw,color = colors[k],label = 'kz %.2f$\mu m^{-1}$'%(kz))
        if k == 0:
            line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
            
        k = k + 1
    
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 4,markerscale=1,fontsize=tamlegend,frameon=0,handletextpad=0.2)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('omegac_vs_mu_tot_modo%i.png'%(modo) ,format = 'png')
        
#%%
 
variacion_im_epsi1 = []

kz1, kz2 = [np.min(list_kz_chicos),np.max(list_kz_chicos)]

for modo in [1,2,3,4]: 
    os.chdir(path_load)
    name1 = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz1,modo)
    tabla1 = np.loadtxt(name1, delimiter='\t', skiprows=1)
    tabla1 = np.transpose(tabla1)
    [list_mu_opt1,omegac_opt1,epsi1_imag_opt1,eq_det1] = tabla1
    
    
    name2 = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz2,modo)
    tabla2 = np.loadtxt(name2, delimiter='\t', skiprows=1)
    tabla2 = np.transpose(tabla2)
    [list_mu_opt2,omegac_opt2,epsi1_imag_opt2,eq_det2] = tabla2
    
    dif = []
    indmin = np.min([len(epsi1_imag_opt1), len(epsi1_imag_opt2)])
    for ind in range(indmin):
        delta = (epsi1_imag_opt1[ind] - epsi1_imag_opt2[ind])
        dif.append(np.abs(delta/epsi1_imag_opt1[ind]))
    
    variacion_im_epsi1.append(np.mean(dif))

#%%
    













