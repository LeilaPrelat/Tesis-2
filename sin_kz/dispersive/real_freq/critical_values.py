# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre paper.py: graficar optical gain y omega/c obtenidos

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

#%%

save_graphs = 1

tamfig = (12,7)
tamlegend = 20
tamletra = 20
tamtitle = 20
tamnum = 18

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')
path_graphene = path_basic.replace('/sin_kz/dispersive/real_freq','') 


try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

epsiinf_DL = 3.9
R = 0.35              #micrones
Ep = 0.7
gamma_DL = 0.01 #unidades de energia

path_load = path_basic + '/R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)

info1 = 'R = %.2f $\mu$m, $\epsilon_\infty$ = %.2f, Ep = %.1f, $\gamma_{DL}$ = %.2f' %(R,epsiinf_DL,Ep,gamma_DL)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_dispR_%.2f.txt' %(R), [inf_tot], fmt='%s')

list_color = ['darkred','yellowgreen','steelblue','coral']

#%%

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%


# colors = sns.set_palette("RdBu_r")

# current_palette = sns.color_palette("RdBu_r", 4)
# list_color = current_palette

# sns.palplot(current_palette)
sns.set()

hspace = 0.12
wspace = 0.1

fig,axs = plt.subplots(3,1, sharex=True, facecolor='w', figsize = (12,10))
plt.subplots_adjust(hspace =hspace,wspace = wspace)

labelx = '$\mu_c$ [eV]'
labely = '$\epsilon_{ci}$'

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla

    if modo == 1:
        axs[0].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
    elif modo in [3,4]:
        axs[2].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
    else:
        axs[1].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
    k = k + 1

if R == 0.05 and epsiinf_DL == 3.9 and Ep == 0.6:
    ticks1 = [-0.14,-0.152,-0.164,-0.176,-0.188]
    ticks2 = [-0.100,-0.12,-0.14,-0.16,-0.18]
    ticks3 = [-0.07,-0.10,-0.13,-0.16,-0.19]
elif R == 0.5 and epsiinf_DL == 3.9 and Ep == 0.6:
    ticks1 = [-0.929,-0.932,-0.935,-0.938,-0.941,-0.944]
    ticks2 = [-0.226,-0.228,-0.230,-0.232,-0.234,-0.236]
    ticks3 = [-0.157,-0.164,-0.171,-0.178,-0.185]    



for i in [0,1,2]:
    axs[i].minorticks_on()
    axs[i].tick_params(labelsize = tamnum)
    if R == 0.05 and epsiinf_DL == 3.9 and Ep == 0.6 or R == 0.5 and epsiinf_DL == 3.9 and Ep == 0.6 :
        ticks = [ticks1,ticks2,ticks3]
        axs[i].set_yticks(ticks[i])
        axs[i].set_yticklabels(ticks[i])
    axs[i].set_ylabel(labely,fontsize = int(tamletra*1.2),labelpad = 6) 
    
# axs[2].set_ylabel(labely,fontsize =  int(tamletra*1.2),labelpad=6) 
# axs[1].set_ylabel(labely,fontsize =  int(tamletra*1.2),labelpad=6) 

axs[2].set_xlabel(labelx,fontsize = tamletra,labelpad=4)

legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)
        
fig.legend(handles=patchList2,loc=[0.145,0.878],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 

if R == 0.5:
    s1,s2,s3 = 'a)', 'b)', 'c)'
else:
    s1,s2,s3 = 'd)', 'e)', 'f)'
xs = 0.133
fig.text(xs,0.84,s1,fontsize = tamletra) 
fig.text(xs,0.575,s2,fontsize = tamletra) 
fig.text(xs,0.307,s3,fontsize = tamletra) 

if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_muR_%.2f.png' %(R),format = 'png')
    
#%% 

labely = '$\omega/c$'

plt.figure(figsize=tamfig)

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
    plt.plot(list_mu_opt,omegac_opt,'-',lw = 5,color = list_color[k])
        
    k = k + 1

plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra,labelpad=2)
plt.tick_params(labelsize = tamnum)

legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)
        
plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 

if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('omegac_vs_muR_%.2f.png' %(R),format = 'png')
    
#%%
 
