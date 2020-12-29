#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

polos complejos del caso no dispersivo
sin kz 

"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

#%% 

save_graphs = 1 #guardar los graficos
break_axes = 1

tamfig = (12,9)
tamlegend = 16
tamletra = 16
tamtitle = 16
tamnum = 14

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/sin_kz/non-dispersive','') 

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

modo = 4
R = 0.5 #micrones
re_epsi1 = 4.9

list_mu = [0.3,0.5,0.7]

info = 'Modo = %i, R = %.2f $\mu$m, Re($\epsilon_1$) = %.2f' %(modo,R,re_epsi1) +  ', ' + name_this_py

path_g = path_basic + r'/re_epsi1_%.2f_vs_mu' %(re_epsi1)
path_load = path_g  + '/' + 'find_Lambda' + '/' + 'modo_%i' %(modo)
os.chdir(path_load)

#%%

if re_epsi1 not in [3.9,4.9,0,-0.5,-0.05]:
    raise TypeError('Wrong value for re(epsi1)')

if R!=0.5:
    raise TypeError('Wrong value for R')
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')

#%%

#hallar minimos y maximos, y valores criticos para guardar en un .txt

list_im_epsi1 = []
index_crit = []

maxis_y = []
minis_y = []
maxis_x = []
minis_x = []
for mu in list_mu:
    
    tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
    tabla = np.transpose(tabla)       
    try: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
    except ValueError: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
        
    omega_c_real = []
    omega_c_imag = []
    for k in range(len(lambda_real_opt)):
        lambda_real = lambda_real_opt[k]
        lambda_imag = lambda_imag_opt[k]
        lambbda = lambda_real + 1j*lambda_imag
        omega_c = 2*pi/lambbda
        omega_c_real.append(omega_c.real)
        omega_c_imag.append(omega_c.imag)

    maxis_x.append(np.max(omega_c_real))
    minis_x.append(np.min(omega_c_real))
    
    maxis_y.append(np.max(omega_c_imag))
    minis_y.append(np.min(omega_c_imag))
    
    #hallar en donde se anula la parte imaginaria de omega/c (valor critico)

    omega_c_imag2 = np.abs(omega_c_imag)
    ind3 = np.argmin(omega_c_imag2)
    
    list_im_epsi1.append(epsi1_imag_opt[ind3])
    index_crit.append(ind3)

#%%        

print('Graficar')

list_color = ['purple','darkblue','darkgreen']


#fig,ax = plt.subplots(1,len(list_mu), sharey=True, gridspec_kw={'hspace': 1, 'wspace': 0.01}, figsize = (12,14))
wspace = 0.0005
if break_axes==1:
    wspace = 0.05

fig,ax = plt.subplots(1,len(list_mu), sharey=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace = 0.001,wspace=wspace)



i = 0
for mu in list_mu:
    
    tabla = np.loadtxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
    tabla = np.transpose(tabla)       
    try: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
    except ValueError: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
    
    print(len(epsi1_imag_opt))
    
    omega_c_real = []
    omega_c_imag = []
    for k in range(index_crit[i]):
        lambda_real = lambda_real_opt[k]
        lambda_imag = lambda_imag_opt[k]
        lambbda = lambda_real + 1j*lambda_imag
        omega_c = 2*pi/lambbda
        omega_c_real.append(omega_c.real)
        omega_c_imag.append(omega_c.imag)    
    
    omega_c_imag2 = np.abs(omega_c_imag)
    ind3 = np.argmin(omega_c_imag2)
    
#    list_im_epsi1.append(epsi1_imag_opt[ind3])
    label_graph_im_eps1_3 = 'Im($\epsilon_1)_c$ = %.7f' %(epsi1_imag_opt[ind3])
    
    ax[i].plot(omega_c_real,omega_c_imag,'-',lw = 5,color = list_color[i],label = '$\mu_c$ = %.1f' %(mu))
    #ax[i].plot(omega_c_real[ind3],omega_c_imag[ind3],'.', color = 'darkred',ms=10) 
    #label = label_graph_im_eps1_3
    #ax[i].set_xscale('log')
    i = i + 1

ax[0].set_ylabel('Im($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamletra,labelpad=5)
ax[1].set_xlabel('Re($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamtitle,labelpad=17)


ticks1_mod1 = [0.1221677,0.1221681]
ticks2_mod1 = [0.1574907,0.157491]
ticks3_mod1 = [0.1860559,0.1860565]

ticks1_mod2 = [0.1728754,0.1728758]
ticks2_mod2 = [0.2231034,0.2231037]
ticks3_mod2 = [0.2637985,0.2637987]

ticks1_mod3 = [0.2116442,0.2116446]
ticks2_mod3 = [0.2732946,0.2732949]
ticks3_mod3 = [0.3232751,0.3232753]

ticks1_mod4 = [0.244244,0.2442445]
ticks2_mod4 = [0.3155353,0.3155356]
ticks3_mod4 = [0.3733461,0.3733463]

ticks_mod1 = [ticks1_mod1,ticks2_mod1,ticks3_mod1]
ticks_mod2 = [ticks1_mod2,ticks2_mod2,ticks3_mod2]
ticks_mod3 = [ticks1_mod3,ticks2_mod3,ticks3_mod3]
ticks_mod4 = [ticks1_mod4,ticks2_mod4,ticks3_mod4]

ax[0].tick_params(axis='y',left=True, right=False,labelsize=tamnum)
ax[0].tick_params(axis='x',bottom=True, top=False,labelsize=tamnum)

ax[1].tick_params(axis='y', right=False, left = False)
ax[1].tick_params(axis='x', bottom=True, top=False,labelsize=tamnum)

ax[2].tick_params(axis='y', right=False, left = False)
ax[2].tick_params(axis='x', bottom=True, top=False,labelsize=tamnum)


if break_axes==1:

    # hide the spines between ax and ax2
    ax[0].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    #ax[0].yaxis.tick_left()
    #ax[0].tick_params(labelright='off')
    #ax[1].yaxis.tick_right()

    
    d = 0.015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform = ax[0].transAxes, color='k', clip_on=False)
    ax[0].plot((1-d,1+d), (-d,+d), **kwargs)
    ax[0].plot((1-d,1+d),(1-d,1+d), **kwargs)
    
    kwargs.update(transform=ax[1].transAxes)  # switch to the bottom axes
    ax[1].plot((-d,+d), (1-d,1+d), **kwargs)
    ax[1].plot((-d,+d), (-d,+d), **kwargs)
    
    # hide the spines between ax and ax2
    ax[1].spines['right'].set_visible(False)
    ax[2].spines['left'].set_visible(False)
    # ax[1].yaxis.tick_left()
    # ax[1].tick_params(labelright='off')
    # ax[2].yaxis.tick_right()
    
    #d = 0.015 # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform = ax[1].transAxes, color='k', clip_on=False)
    ax[1].plot((1-d,1+d), (-d,+d), **kwargs)
    ax[1].plot((1-d,1+d),(1-d,1+d), **kwargs)
    
    kwargs.update(transform=ax[2].transAxes)  # switch to the bottom axes
    ax[2].plot((-d,+d), (1-d,1+d), **kwargs)
    ax[2].plot((-d,+d), (-d,+d), **kwargs)


for i in [0,1,2]:

    if modo==1:
        ticks = ticks_mod1[i]
    elif modo==2:
        ticks = ticks_mod2[i]
    elif modo==3:
        ticks = ticks_mod3[i]     
    elif modo==4:
        ticks = ticks_mod4[i]
        
    if re_epsi1 == 4.9:
        ax[i].set_xticks(ticks)
        ax[i].set_xticklabels(ticks)
    
    ax[i].legend(loc=[0.25,1],markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

if save_graphs==1:
    #plt.tight_layout(pad=wspace)
    os.chdir(path_g)
    plt.savefig('complex_poles%i'%(modo))

tabla = np.array([list_mu, list_im_epsi1],dtype=object)
tabla = np.transpose(tabla)

header1 = 'mu[eV]     Im(epsi1)' + ' inf: ' + info
np.savetxt('info_nondisp_complex_poles_modo%i.txt' %(modo), tabla, fmt='%s',header = header1)

#%%

