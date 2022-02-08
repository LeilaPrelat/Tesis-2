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

save_graphs = 1
paper = 1

#%%

if paper == 1: 
    tamfig = (6,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    
    labelpady = -0.5
    labelpadx = -0.5
    lw = 1.5
    pad = -2.5
    loc1 = [0.16,0.88]
    loc2 = [0.05,1]
    columnspace = 0.5
    
else:
    tamfig = (11,7)
    tamlegend = 18
    tamletra = 20
    tamtitle = 20
    tamnum = 16
    
    labelpady = 0
    labelpadx = 4
    lw = 3.5
    pad = 0
    loc1 = [0.15,0.88]
    loc2 = [0.025,1]

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

path_load = path_basic + '/re_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)

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

hspace = 0.1
wspace = 0.08

for kz in list_kz: 
    fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
    plt.subplots_adjust(hspace =hspace,wspace = wspace)
    
    labelx = '$\mu_c$ [eV]'
    labely = '[Im($\epsilon_1$)]$_c$'
    
    if kz == 0.05:
        ticks_y1 = np.linspace(-0.036,-0.032,5)
        ticks_y2 = np.linspace(-0.03,-0.01,5)
        
    if kz == 0.14:
        ticks_y1 = np.linspace(-0.78,-0.76,5)
        ticks_y2 = np.linspace(-0.03,-0.01,5)
    
    k = 0
    for modo in [1,2,3,4]:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
    
        if modo == 1:
            line1, = axs[0].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = 'mode %i'%(k+1))
            line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
            if kz in [0.05,0.14]:
                axs[0].set_yticks(ticks_y1)
        else:
            axs[1].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = 'mode %i'%(k+1))
            if kz in [0.05,0.14]:
                axs[1].set_yticks(ticks_y2)
        k = k + 1
    
    
    for i in [0,1]:
        axs[i].minorticks_on()
        axs[i].tick_params(labelsize = tamnum,pad = pad)

    
    axs[0].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
    axs[1].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
    axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad =labelpadx) 
    
    #fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
    
    # s1,s2 = 'a)', 'b)'
    # xs = 0.133
    # fig.text(xs,0.83,s1,fontsize = tamletra) 
    # fig.text(xs,0.422,s2,fontsize = tamletra) 
    
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('Im_epsi1_vs_mu_kz%.4f.png'%(kz) ,format = 'png')
        
    labely = '$\omega/c$ [$\mu$m$^{-1}$]'
    
    plt.figure(figsize=tamfig)
    
    k = 0
    for modo in [1,2,3,4]:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        line1, = plt.plot(list_mu_opt,omegac_opt,symbols[k],lw =lw,color = colors[k],label = 'modo %i'%(k+1))
        if k == 0:
            line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
            
        k = k + 1
    
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('omegac_vs_mu_kz%.4f.png'%(kz) ,format = 'png')
        
#%%
 
