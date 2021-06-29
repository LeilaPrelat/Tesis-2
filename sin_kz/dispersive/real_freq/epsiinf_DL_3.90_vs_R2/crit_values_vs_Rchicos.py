# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values.py: graficar optical gain y omega/c obtenidos
en formato paper (sin titulo y prolijo)


"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import seaborn as sns

#%%

save_graphs = 1 #guardar los graficos
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
    pad = -1
    loc1 = [0.22,0.88]
    loc2 = [0.1,1]
    
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
path_basic2 =  path_basic.replace('/epsiinf_DL_3.90_vs_R2','')

#Para condiciones iniciales
try:
    sys.path.insert(1, path_basic2)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi
except ModuleNotFoundError:
    print('QE_lossless.py no se encuentra en ' + path_basic2)
    path_basic = input('path de la carpeta donde se encuentra QE_lossless.py')
    sys.path.insert(1, path_basic)
    from QE_lossless import im_epsi1_cuasi,omegac_cuasi

#%%

print('Definir parametros del problema')

epsiinf_DL = 3.9
R = 0.05             #micrones
gamma_DL = 0.01 #unidades de energia
Ep = 0.1
mu = 0.3 
modo = 4

ind = 50

#%%

if epsiinf_DL!= 3.9:
    raise TypeError('Wrong value for Re(epsilon1)')

if mu != 0.3:
    raise TypeError('Wrong value for mu')

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')
    
colors = ['darkred','yellowgreen','steelblue','coral']
symbols = ['-','-','--','-.']
sns.set()

#%%

path_load = path_basic + r'/Ep_%.2f' %(Ep)
os.chdir(path_load)
data_load = np.loadtxt('opt_det_sinkz_vs_R_modo%i.txt' %(modo), delimiter='\t', skiprows = 1)
data_load = np.transpose(data_load)
[barrido_R,omegac_opt,epsi1_imag_opt,eq_det] = data_load
[barrido_R,omegac_opt,epsi1_imag_opt,eq_det]= [barrido_R[ind:-1],omegac_opt[ind:-1],epsi1_imag_opt[ind:-1],eq_det[ind:-1]]

im_epsi1_QE_lossless = []
omegac_QE_lossless = []
for R in barrido_R:
    a = omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,mu)
    b = im_epsi1_cuasi(a,Ep,epsiinf_DL,gamma_DL,modo,R,mu) 
    im_epsi1_QE_lossless.append(b)
    omegac_QE_lossless.append(a)

#%%

labelx = 'R [$\mu$m]'
labely1 = '$[\epsilon_{di}]_c$'
labely2 = '$\omega/c$ [$\mu$m$^{-1}$]'

print('Graficar Im(epsi1) vs R') 

plt.figure(figsize=tamfig)
plt.plot(barrido_R,epsi1_imag_opt,symbols[0],lw =lw,color = colors[0],label='sol num')
plt.plot(barrido_R,im_epsi1_QE_lossless,symbols[1],lw =lw,color = colors[1],label='Aprox QE')
plt.ylabel(labely1,fontsize=int(tamletra*1.3),labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)  
if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_Rchico_%i.png'%(modo) ,format = 'png')
 

plt.figure(figsize=tamfig)
plt.plot(barrido_R,omegac_opt,symbols[2],lw =lw,color = colors[2],label='sol num')
plt.plot(barrido_R,omegac_QE_lossless,symbols[3],lw =lw,color = colors[3],label='Aprox QE')
plt.ylabel(labely2,fontsize=int(tamletra*1.1),labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend,frameon=0,handletextpad=0.5)  
if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('omegac_vs_Rchico_%i.png'%(modo) ,format = 'png')

#%%

info1 = '$\mu_c$ = %.2f eV, $\epsilon_\infty$ = %.2f, Ep = %.1f, $\gamma_{DL}$ = %.2f' %(mu,epsiinf_DL,Ep,gamma_DL)
inf_tot = info1 + '. Ver ' + name_this_py
info_txt = 'info_critical_values_disp_vsR.txt' 

os.chdir(path_load)
np.savetxt(info_txt, [inf_tot], fmt='%s')

#%%
