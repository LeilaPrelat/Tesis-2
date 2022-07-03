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
opcion1 = 1

#%%

if paper == 1: 

    if opcion1 == 1:
    
        tamfig = (3, 2)
        tamletra = 9
        tamnum = 8.5
        tamlegend = 8.5
        ticks_x = [0.4,0.6,0.8]
        ticks_y1_kz1 = [-0.036,-0.034,-0.032]
        ticks_y2_kz1 = [-0.012,-0.020,-0.028]
    
        ticks_y1_kz2 = [-0.761, -0.770, -0.779]
        ticks_y2_kz2 = [-0.012,-0.020,-0.028]
            
        columnspace = 0.5
        markerscale = 0.7
        loc1 = [0.275,0.9]  
        length_marker = 1
        labelpady = 1.1
        labelpadx = 0.7
        
    else:    
        tamfig = (3, 2)
        tamletra = 7
        tamnum = 6
        tamlegend = 6
        ticks_x = [0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        ticks_y1_kz1 = [-0.036,-0.035,-0.034,-0.033,-0.032]
        ticks_y2_kz1 = [-0.010,-0.015,-0.020,-0.025,-0.030]
        
        ticks_y1_kz2 = [-0.76 , -0.765, -0.77 , -0.775, -0.78 ]
        ticks_y2_kz2 = [-0.010,-0.015,-0.020,-0.025,-0.030]

        columnspace = 3.5
        markerscale = 0.7
        loc1 = [0.052,1]
        length_marker = 3
        labelpady = 0.2
        labelpadx = 0.2

    lw = 1
    pad = -4        
    loc2 = [0.231,0.9]
    columnspace2 = 0.9

    dpi = 500
    
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
path_save = path_basic + '/check_plasmon_behavior'

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)

try:
    sys.path.insert(1, path_graphene)
    from graphene_sigma import sigma_DL
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)


pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

def eq_check_plasmon(mu,omegac,epsi1,kz,R,modo):
    
    energy = omegac*(hb*c)
    
#    sigma_graphene, inter, intra  = sigma(energy,mu,hbargama)
    
    sigma_graphene  = sigma_DL(energy,mu,hbargama)
    
    den = 4*np.pi*alfac*sigma_graphene
    
    tot = omegac*(epsi1 + epsi2)/den
    
    return tot**2 + kz**2 + (modo/R)**2 

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
list_kz = [0,0.0001,0.05,0.1,0.135,0.14,0.5] 
list_kz = [0,0.0001,0.05,0.1,0.135,0.14,0.5] 
list_modos = [1,2,3,4]

path_load = path_basic + '/re_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_conkz.txt', [inf_tot], fmt='%s')

colors = ['darkred','yellowgreen','coral','steelblue']
symbols = ['-','-','--','-.']
sns.set()

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

#%%

hspace = -0.1
wspace = 0.03

for kz in list_kz: 
    
    k = 0
    for modo in list_modos: 
        
        labelx = '$\mu_c$ [eV]'
        labely1 = 'Re(Eq = 0)'
        labely2 = 'Im(Eq = 0)'
        labely3 = '|Eq = 0|'
        
        label_eje_mu = r'$\nu$ = %i, kz = %.2f 1/$\mu$m'%(modo,kz)
      

        
        os.chdir(path_load)
        
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        
        list_y_re = []
        list_y_im = []
        list_y_abs = []
        for j in range(len(list_mu_opt)):
            mu = list_mu_opt[j]
            omegac = omegac_opt[j]
            epsi1_imag = epsi1_imag_opt[j]
            epsi1 = re_epsi1 + 1j*epsi1_imag
            ejey = eq_check_plasmon(mu,omegac,epsi1,kz,R,modo)
            list_y_re.append(np.real(ejey))
            list_y_im.append(np.imag(ejey))
            list_y_abs.append(np.abs(ejey))

        plt.figure(figsize=tamfig)        
        plt.plot(list_mu_opt,list_y_re,lw =lw,color = colors[k])
        plt.ylabel(labely1,fontsize=tamletra,labelpad =labelpady)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.tick_params(labelsize = tamnum,pad = pad)
        plt.xticks(ticks_x)  
                
        #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        plt.legend(loc = 'best',frameon=0,handletextpad=0.2, handlelength=length_marker)     
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_save)
            plt.savefig('re_Eq_vs_mu_kz%.4f_modo%i.png'%(kz,modo) ,format = 'png',dpi = dpi)


        plt.figure(figsize=tamfig)   
        plt.plot(list_mu_opt,list_y_im,lw =lw,color = colors[k])
        plt.ylabel(labely2,fontsize=tamletra,labelpad =labelpady)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.tick_params(labelsize = tamnum,pad = pad)
        plt.xticks(ticks_x)  
                
        #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        plt.legend(loc = 'best',frameon=0,handletextpad=0.2, handlelength=length_marker)     
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_save)
            plt.savefig('im_Eq_vs_mu_kz%.4f_modo%i.png'%(kz,modo) ,format = 'png',dpi = dpi)
            


        plt.figure(figsize=tamfig)   
        plt.plot(list_mu_opt,list_y_abs,lw =lw,color = colors[k])
        plt.ylabel(labely3,fontsize=tamletra,labelpad =labelpady)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.tick_params(labelsize = tamnum,pad = pad)
        plt.xticks(ticks_x)  
                
        #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        plt.legend(loc = 'best',frameon=0,handletextpad=0.2, handlelength=length_marker)     
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_save)
            plt.savefig('abs_Eq_vs_mu_kz%.4f_modo%i.png'%(kz,modo) ,format = 'png',dpi = dpi)            
    
        k = k + 1
            
#%%
