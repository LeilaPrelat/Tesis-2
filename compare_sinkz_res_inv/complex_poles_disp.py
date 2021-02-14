#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

comparar los polos complejos del caso dispersivo
( sin kz ) del caso de resonancia y de invisibilidad

"""

import numpy as np
import os
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
tamnum = 15

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('compare_sinkz_res_inv/' + name_this_py,'')
path_load_inv = path_basic + 'sin_kz_inv/dispersive/complex_freq'
path_load_res = path_basic + 'sin_kz/dispersive/complex_freq'

#%%

print('Definir parametros del problema')

R = 0.1             #micrones
modo = 4

Ep = 0.9
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

list_mu = [0.3,0.5,0.7]

info1 = 'R = %.2f $\mu$m, $E_p$ = %.3f eV, modo = %i' %(R,Ep,modo)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
info =  ', ' + info1 + ', ' + info2  + ', ' + name_this_py

path_load = r'/epsiinf_DL_%.2f_vs_mu/R_%.2f/Ep_%.1f/find_Lambda/modo_%i' %(epsiinf_DL,R,Ep,modo)
name_txt = 'opt_det_nano_'
path_save = '/complex_poles_disp'

#%%
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')

#%%

print('resonancia: hallar minimos y maximos, y valores criticos para guardar en un .txt')
print('resonancia: separar lo que es modo (im(omega/c)) <= 0 de lo que no es')

maxis_x = []
minis_x = []

omegac_real_tot11 = []
omegac_imag_tot11 = []

omegac_real_tot12 = []
omegac_imag_tot12 = []

#hallar en donde se anula la parte imaginaria de omega/c (valor critico)

list_im_epsi1_res = []
index_crit_res = []

for mu in list_mu:
    os.chdir(path_load_res + path_load)
    tabla = np.loadtxt('opt_det_nano_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
    tabla = np.transpose(tabla)       
    try: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
    except ValueError: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
        
    omega_c_real1 = []
    omega_c_imag1 = []
    omega_c_real2 = []
    omega_c_imag2 = []
    
    for k in range(len(lambda_real_opt)):
        lambda_real = lambda_real_opt[k]
        lambda_imag = lambda_imag_opt[k]
        lambbda = lambda_real + 1j*lambda_imag
        omega_c = 2*np.pi/lambbda
        if omega_c.imag <=0 :  
            omega_c_real1.append(omega_c.real)
            omega_c_imag1.append(omega_c.imag)
        else:
            omega_c_real2.append(omega_c.real)
            omega_c_imag2.append(omega_c.imag)            

    omegac_real_tot11.append(omega_c_real1)
    omegac_imag_tot11.append(omega_c_imag1)   
    omegac_real_tot12.append(omega_c_real2)
    omegac_imag_tot12.append(omega_c_imag2)  

    maxis_x.append(np.max(omega_c_real1))
    try:
        minis_x.append(np.min(omega_c_real2))
    except ValueError:
        minis_x.append(np.min(omega_c_real1))
    
    #hallar en donde se anula la parte imaginaria de omega/c (valor critico)

    omega_c_imag_crit = np.abs(omega_c_imag1)
    ind3 = np.argmin(omega_c_imag_crit)
    
    list_im_epsi1_res.append(epsi1_imag_opt[ind3])
    index_crit_res.append(ind3)

#%%

print('invisibilidad: hallar minimos y maximos, y valores criticos para guardar en un .txt')

maxis2_x = []
minis2_x = []

omegac_real_tot21 = []
omegac_imag_tot21 = []

omegac_real_tot22 = []
omegac_imag_tot22 = []

#hallar en donde se anula la parte imaginaria de omega/c (valor critico)

list_im_epsi1_inv = []
index_crit_inv = []

i = 0 
for mu in list_mu:
    os.chdir(path_load_inv + path_load)
    tabla = np.loadtxt('opt_det_nano_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows=1)     
    tabla = np.transpose(tabla)       
    try: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_gn] = tabla
    except ValueError: 
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
        
    omega_c_real1 = []
    omega_c_imag1 = []
    omega_c_real2 = []
    omega_c_imag2 = []
    im_epsi1_crit_res = list_im_epsi1_res[i] 
    for k in range(len(lambda_real_opt)):
        lambda_real = lambda_real_opt[k]
        lambda_imag = lambda_imag_opt[k]
        im_epsi1_crit_inv = epsi1_imag_opt[k]
        lambbda = lambda_real + 1j*lambda_imag
        omega_c = 2*np.pi/lambbda
        
        
        if np.abs(im_epsi1_crit_inv) <= np.abs(im_epsi1_crit_res) :  
            omega_c_real1.append(omega_c.real)
            omega_c_imag1.append(omega_c.imag)
        else:
            omega_c_real2.append(omega_c.real)
            omega_c_imag2.append(omega_c.imag)            

    omegac_real_tot21.append(omega_c_real1)
    omegac_imag_tot21.append(omega_c_imag1)   
    omegac_real_tot22.append(omega_c_real2)
    omegac_imag_tot22.append(omega_c_imag2)  
    
    maxis2_x.append(np.max(omega_c_real1))
    try:
        minis2_x.append(np.min(omega_c_real2))
    except ValueError:
        minis2_x.append(np.min(omega_c_real1))
        
    #hallar en donde se anula la parte imaginaria de omega/c (valor critico)

    omega_c_imag_crit = np.abs(omega_c_imag1)
    ind3 = np.argmin(omega_c_imag_crit)
    
    list_im_epsi1_inv.append(epsi1_imag_opt[ind3])
    index_crit_inv.append(ind3)
    i = i + 1

#%%


def ticksx_res(modo):
    tick0 = []
    for j in range(len(list_mu)):
        maxi = maxis_x[j]
        mini = minis_x[j] 
        x = []
        delta = maxi - mini
        # print(delta/mini)
        for k in [0,2,4,6]:
            if delta/mini < 1e-3:
                value = np.round(mini + k*1e-4,4)
            else:
                value = np.round(mini + k*1.7*1e-4,4)
            if value < maxi :
                x.append(value)
        tick0.append(x)
    print('res:')
    print(tick0)
    return tick0

def ticksx_inv(modo):
    tick0 = []
    for j in range(len(list_mu)):
        maxi = maxis2_x[j]
        mini = minis2_x[j] 
        
        x = []
        for k in [0,2,4,6]:
            value = np.round(mini + k*1e-4,4)
            if value < maxi:
                x.append(value)
        tick0.append(x)
        
    print('inv:')
    print(tick0)    
    return tick0

#%%        

print('Graficar resonancia')

list_color = ['purple','darkblue','darkgreen']


#fig,ax = plt.subplots(1,len(list_mu), sharey=True, gridspec_kw={'hspace': 1, 'wspace': 0.01}, figsize = (12,14))
wspace = 0.0005
if break_axes==1:
    wspace = 0.05

fig,ax = plt.subplots(1,len(list_mu), sharey=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace = 0.001,wspace=wspace)

i = 0
for mu in list_mu:
    labelmu = '$\mu_c$ = %.1f' %(mu)
    omegac_real1 = omegac_real_tot11[i]
    omegac_imag1 = omegac_imag_tot11[i]
    
    omegac_real2 = omegac_real_tot12[i]
    omegac_imag2 = omegac_imag_tot12[i]
    
    ax[i].plot(omegac_real1,omegac_imag1,'-',lw = 5,color = list_color[i],label = labelmu)
    ax[i].plot(omegac_real2,omegac_imag2,'--',lw = 5,color = list_color[i],alpha = 0.5)
    
    i = i + 1


ax[0].set_ylabel('Im($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamletra,labelpad=5)
ax[1].set_xlabel('Re($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamtitle,labelpad=17)

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
    # if modo == 1:
    #     ticks = ticksx_res(modo)[i]
    #     ax[i].set_xticks(ticks)
    #     ax[i].set_xticklabels(ticks)
    ax[i].legend(loc=[0.25,1],markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

if save_graphs==1:
    #plt.tight_layout(pad=wspace)
    os.chdir(path_basic + path_save)
    plt.savefig('complex_poles_res_disp%i'%(modo))

#%%

print('Graficar invisbilidad')

#fig,ax = plt.subplots(1,len(list_mu), sharey=True, gridspec_kw={'hspace': 1, 'wspace': 0.01}, figsize = (12,14))
wspace = 0.0005
if break_axes==1:
    wspace = 0.05

fig,ax = plt.subplots(1,len(list_mu), sharey=True, facecolor='w', figsize = tamfig)
plt.subplots_adjust(hspace = 0.001,wspace=wspace)

i = 0
for mu in list_mu:
    labelmu = '$\mu_c$ = %.1f' %(mu)
    omegac_real1 = omegac_real_tot21[i]
    omegac_imag1 = omegac_imag_tot21[i]
    
    omegac_real2 = omegac_real_tot22[i]
    omegac_imag2 = omegac_imag_tot22[i]
    
    
    ax[i].plot(omegac_real1,omegac_imag1,'-',lw = 5,color = list_color[i],label = labelmu)
    ax[i].plot(omegac_real2,omegac_imag2,'--',lw = 5,color = list_color[i],alpha = 0.5)
    
    i = i + 1

ax[0].set_ylabel('Im($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamletra,labelpad=5)
ax[1].set_xlabel('Re($\omega/c$) [$\mu$m]$^{-1}$',fontsize = tamtitle,labelpad=17)

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
    # if modo == 1:
    #     ticks = ticksx_inv(modo)[i]
    #     ax[i].set_xticks(ticks)
    #     ax[i].set_xticklabels(ticks)
    ax[i].legend(loc=[0.25,1],markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)

if save_graphs==1:
    #plt.tight_layout(pad=wspace)
    os.chdir(path_basic + path_save)
    plt.savefig('complex_poles_inv_disp%i'%(modo))

#%%

tabla = np.array([list_mu, list_im_epsi1_res],dtype=object)
tabla = np.transpose(tabla)

header1 = 'mu[eV]     Im(epsi1)' + ' inf: caso resonancia' + info
np.savetxt('info_disp_res_complex_poles_modo%i.txt' %(modo), tabla, fmt='%s',header = header1)

tabla_inv = np.array([list_mu, list_im_epsi1_inv],dtype=object)
tabla_inv = np.transpose(tabla_inv)

header1 = 'mu[eV]     Im(epsi1)' + ' inf: caso invisibilidad' + info
np.savetxt('info_disp_inv_complex_poles_modo%i.txt' %(modo), tabla_inv, fmt='%s',header = header1)

#%%

