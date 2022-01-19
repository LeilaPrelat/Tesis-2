# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:42:33 2020

@author: Usuario
"""

"""
Sobre compare_re_epsi1.py: este .py toma los .txt de la carpeta find_im_epsi1
y grafica las curvas im_epsi1 vs R comparando con los distintos valores
de Re(epsi1). 

hecho para determinante sin kz

"""

import numpy as np
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

#%%

print('Definir parametros del problema')

barrido_mu = [0.3,0.6,0.9] 
modo = 4
R = 0.25 #micrones

# if R!= 0.5:
#     raise TypeError('Wrong value for radium')

if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for modo')

#%%

print('Graficar Im(epsi1) vs R para diferentes valores de Re(epsi1)') 
title = 'Modo = %i, R = %.2f $\mu$m' %(modo,R) +  ', ' + name_this_py
labelx = 'Re($\epsilon_1$)'

fig1 = plt.figure(figsize=tamfig)
for mu in barrido_mu:
    path_load = path_basic + r'/R_%.2f_vs_re_epsi1/mu_%.1f' %(R,mu)
    os.chdir(path_load)
    data_load = np.loadtxt('opt_det_sinkz_vs_re_epsi1_modo%i.txt' %(modo), delimiter='\t', skiprows = 1)
    data_load = np.transpose(data_load)
    [barrido_re_epsi1,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    plt.plot(barrido_re_epsi1,epsi1_imag_opt,'.',ms=10,label = r'$\mu$ = %.1f' %(mu))
#    plt.plot(barrido_R,rta_Im_epsi1_teo2,'*',ms=5,label='Modo = %i, aprox cuasi-estatica'%(modo))

plt.title(title,fontsize=tamtitle)
plt.ylabel('Im($\epsilon_1$)',fontsize=tamletra)
# plt.ylim([-0.4,0])
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)      
if save_graphs==1:
    os.chdir(path_basic + '/R_%.2f_vs_re_epsi1' %(R))
    fig1.savefig('Im_epsi1_re_epsi1_modo%i'%(modo))

fig2 = plt.figure(figsize=tamfig)
for mu in barrido_mu:
    path_load = path_basic + r'/R_%.2f_vs_re_epsi1/mu_%.1f' %(R,mu)
    os.chdir(path_load)
    data_load = np.loadtxt('opt_det_sinkz_vs_re_epsi1_modo%i.txt' %(modo), delimiter='\t', skiprows = 1)
    data_load = np.transpose(data_load)
    [barrido_re_epsi1,omegac_opt,epsi1_imag_opt,eq_det] = data_load
    plt.plot(barrido_re_epsi1,omegac_opt,'.',ms=10,label = r'$\mu$ = %.1f' %(mu))

plt.title(title,fontsize=tamtitle)
plt.ylabel('Re($\omega$/c) [1/$\mu$m]',fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)      
if save_graphs==1:
    os.chdir(path_basic + '/R_%.2f_vs_re_epsi1' %(R))
    fig2.savefig('omegac_vs_re_epsi1_modo%i'%(modo))
    
#%%
