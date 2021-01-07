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

tamfig = (11,7)
tamlegend = 18
tamletra = 20
tamtitle = 20
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')
path_graphene = path_basic.replace('/sin_kz_inv/non-dispersive/real_freq','') 


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

re_epsi1 = 4.9
R = 0.5              #micrones
 
path_load = path_basic + '/re_epsi1_%.2f_vs_mu' %(re_epsi1)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_nondisp.txt', [inf_tot], fmt='%s')

list_color = ['darkred','yellowgreen','steelblue','coral']

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

#%%


# colors = sns.set_palette("RdBu_r")

# current_palette = sns.color_palette("RdBu_r", 4)
# list_color = current_palette

# sns.palplot(current_palette)
sns.set()

hspace = 0.1
wspace = 0.08

fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace =hspace,wspace = wspace)

labelx = '$\mu_c$ [eV]'
labely = 'Im($\epsilon_1$)'

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla

    if modo == 1:
        axs[0].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
    else:
        axs[1].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
        
    k = k + 1


ticks1 = [-0.033,-0.035,-0.037,-0.039,-0.041]
ticks2 = [-0.011,-0.017,-0.023,-0.029,-0.035]

ticks = [ticks1,ticks2]

for i in [0,1]:
    axs[i].minorticks_on()
    axs[i].tick_params(labelsize = tamnum)
    axs[i].set_yticks(ticks[i])
    axs[i].set_yticklabels(ticks[i])

axs[0].set_ylabel(labely,fontsize = tamletra,labelpad=0) 
axs[1].set_ylabel(labely,fontsize = tamletra,labelpad=6) 
axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad=2)

legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)
        
fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 

s1,s2 = 'a)', 'b)'
xs = 0.133
fig.text(xs,0.83,s1,fontsize = tamletra) 
fig.text(xs,0.422,s2,fontsize = tamletra) 


if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_mu')
    
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
    plt.savefig('omegac_vs_mu')
    
#%%
 
