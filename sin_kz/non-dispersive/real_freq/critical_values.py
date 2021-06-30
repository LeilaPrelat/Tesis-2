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

import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns

#%%

save_graphs = 1
paper = 1

if paper == 1: 
    tamfig = (5.5,3.5)
    tamlegend = 7
    tamletra = 7
    tamtitle = 6
    tamnum = 6
    labelpady = 0
    labelpadx = 1
    lw = 1.5
    pad = -1
    loc1 = [0.22,0.88]
    loc2 = [0.14,1]
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

os.chdir(path_load)
np.savetxt('info_critical_values_nondisp.txt', [inf_tot], fmt='%s')

colors = ['darkred','yellowgreen','steelblue','coral']
symbols = ['-','-','--','-.']

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

#%%

sns.set()

hspace = 0.1
wspace = 0.08

fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace =hspace,wspace = wspace)

labelx = '$\mu_c$ [eV]'
labely = '[Im($\epsilon_1$)]$_c$'

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla

    if modo == 1:
        line1, = axs[0].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw = lw,color = colors[k],label = 'mode %i'%(k+1))
        line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
    else:
        axs[1].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw = lw,color = colors[k],label = 'mode %i'%(k+1))
        
    k = k + 1


ticks1 = [-0.033,-0.035,-0.037,-0.039,-0.041]
ticks2 = [-0.011,-0.017,-0.023,-0.029,-0.035]

ticks = [ticks1,ticks2]

for i in [0,1]:
    axs[i].minorticks_on()
    axs[i].tick_params(labelsize = tamnum,pad = pad)
    if re_epsi1 == 4.9:
        axs[i].set_yticks(ticks[i])
        axs[i].set_yticklabels(ticks[i])

axs[0].set_ylabel(labely,fontsize = tamletra,labelpad=labelpady) 
axs[1].set_ylabel(labely,fontsize = tamletra,labelpad=labelpady) 
axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad=labelpadx)

legend_mode2 = {'mode 1' : colors[0], 'mode 2' : colors[1],'mode 3' : colors[2], 'mode 4' : colors[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)

#Create custom artists
custom_lines = [Line2D([0], [0], color= colors[0], lw=4, linestyle = symbols[0]),
                Line2D([0], [1], color= colors[1], lw=4, linestyle = symbols[1]),
                Line2D([0], [2], color= colors[2], lw=4, linestyle = symbols[2]), 
                Line2D([0], [3], color= colors[3], lw=4, linestyle = symbols[3])]     

#fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
fig.legend(loc = loc1, ncol = 4,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

s1,s2 = 'a)', 'b)'
xs = 0.133
fig.text(xs,0.83,s1,fontsize = tamletra) 
fig.text(xs,0.422,s2,fontsize = tamletra) 


if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_mu')
    
#%% 

labely = '[Im($\epsilon_1$)]$_c$ QE'
plt.figure(figsize=tamfig)
k = 0
for modo in [1,2,3,4]:
    im_epsi1_QE = []
    omegac_QE = []
    for mu in list_mu_opt:
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        b = im_epsi1_cuasi(a,modo,R,mu) 
        im_epsi1_QE.append(b)
        omegac_QE.append(a)
        
    line1, = plt.plot(list_mu_opt,im_epsi1_QE,symbols[k],lw = lw,color = colors[k],label = 'mode %i'%(k+1))
    if k == 0:
        line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
        
    k = k + 1

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 4,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('Im_epsi1_vs_muQE')    

#%%

labely = '$\omega/c$ $[\mu m]^{-1}$'

plt.figure(figsize=tamfig)

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
    line1, = plt.plot(list_mu_opt,omegac_opt,symbols[k],lw = lw,color = colors[k],label = 'mode %i'%(k+1))
    if k == 0:
        line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
        
    k = k + 1

plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
        
#plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
plt.legend(loc = loc2, ncol = 4,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('omegac_vs_mu')
    
#%%

plt.figure(figsize=tamfig)
k = 0
for modo in [1,2,3,4]:
    im_epsi1_QE = []
    omegac_QE = []
    for mu in list_mu_opt:
        a = omegac_cuasi(modo,R,re_epsi1,mu)
        b = im_epsi1_cuasi(a,modo,R,mu) 
        im_epsi1_QE.append(b)
        omegac_QE.append(a)
        
    line1, = plt.plot(list_mu_opt,omegac_QE,symbols[k],lw = lw,color = colors[k],label = 'mode %i'%(k+1))
    if k == 0:
        line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
        
    k = k + 1

    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.xlabel(labelx,fontsize=tamletra,labelpad=labelpadx)
    plt.tick_params(labelsize = tamnum,pad = pad)
            
    #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    plt.legend(loc = loc2, ncol = 4,markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)
    
    if save_graphs==1:
        os.chdir(path_load)
        plt.savefig('omegac_vs_muQE')    

#%%