#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

barrido en mu. Find omega/c complejo que minimizan
el determinante sin kz para diferentes valores de im(epsilon1) 
para cada valor de mu
"""

import numpy as np
import os
import sys
from scipy.optimize import minimize   
import matplotlib.pyplot as plt

#%% 

save_data_opt = 1 #guardar data de la minimizacion

graficar = 0
save_graphs = 0
if graficar ==1:
    close_graphs = 0

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_det = path_basic.replace('/' + 'complex_freq','')

try:
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante
except ModuleNotFoundError:
    print('det_sinkz_nano.py no se encuentra en ' + path_det)
    path_det = input('path de la carpeta donde se encuentra det_sinkz_nano.py')
    sys.path.insert(1, path_det)
    from det_sinkz_nano import determinante

try:
    sys.path.insert(1, path_det)
    from complex_omegac_QE import omegac_QE
except ModuleNotFoundError:
    print('complex_omegac_QE.py no se encuentra en ' + path_det)
    path_det = input('path de la carpeta donde se encuentra complex_omegac_QE.py')
    sys.path.insert(1, path_det)
    from complex_omegac_QE import omegac_QE

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 1

Ep = 0.3
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

if modo==1:
    list_im_epsi1 = np.linspace(0,-1,1001) 
elif modo==2:
    list_im_epsi1 = np.linspace(0,-0.3,301) 
elif modo==3 or modo == 4:
    list_im_epsi1 = np.linspace(0,-0.2,201) 
    
list_mu =  np.linspace(0.3,0.9,601)  
# list_mu = [0.3]

info1 = 'R = %.1f $\mu$m, $E_p$ = %.3f eV, modo = %i' %(R,Ep,modo)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
info =  ', ' + info1 + ', ' + info2  + ', ' + name_this_py
title = info1 +'\n' + info2  + ', ' + name_this_py

#%%

if list_im_epsi1[0]!= 0 or list_mu[0]!= 0.3:
    raise TypeError('Se necesita empezar im(epsi1)=0 y mu = 0.3 para usar las cond iniciales')

if R!=0.5:
    raise TypeError('Wrong value for radium')
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_graphs==1 or save_data_opt==1:
    try:
        path_data = r'/epsiinf_DL_%.2f_vs_mu/Ep_%.1f/find_Lambda' %(epsiinf_DL,Ep)
        path = path_basic + path_data
        os.chdir(path) 
    except OSError or IOError as error:
        print('Error en el path:', error)
        path = input('Copie y pegue la direc en la cual quiere guardar los graficos-archivos')
        os.chdir(path)

    path_g = path + '/' + 'graficos'
    path_d = path + '/' + 'modo_%i' %(modo)
    
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)
        
    if not os.path.exists(path_d):
        print('Creating folder to save data')
        os.mkdir(path_d)

#%%        

print('Minimizacion del determinante sin kz')
print('Usar QE como condiciones iniciales para el metodo de minimizacion')

mu0 = list_mu[0]
omegac_inicial = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,0,R,mu0)
cond_inicial = 2*np.pi/omegac_inicial

if mu0 == 0.3:
    cond_inicial0 = [cond_inicial.real,cond_inicial.imag]
else:
    cond_inicial0 = [] #mu anterior

if graficar==0:
    os.chdir(path_d)

cond_inicial_mu = []    

tol_NM = 1e-13
ite_NM = 1150
j = 0
for mu in list_mu:
    mu = np.round(mu,4)
    print('')
    print(mu)
    if mu == list_mu[0]:
        cond_inicial = cond_inicial0
    else:
        cond_inicial = cond_inicial_mu[j-1]
           
    epsi1_imag_opt = []
    lambda_real_opt = []
    lambda_imag_opt = []
    eq_det = []
    
    for epsi_ci in list_im_epsi1:
        
        epsi_ci = np.round(epsi_ci,4)
        
    
        def det_2variables(lmbd):
            [Re_lmbd,Im_lmbd] = lmbd
            lambdda = Re_lmbd + 1j*Im_lmbd      
            omegac = 2*np.pi/lambdda
            rta = determinante(omegac,Ep,epsiinf_DL,gamma_DL,epsi_ci,modo,R,mu)
            return np.abs(rta)
        
        res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                   options={'maxiter':ite_NM})
#        print(res.message)
        if res.message == 'Optimization terminated successfully.':
            lambda_real_opt.append(res.x[0])
            lambda_imag_opt.append(res.x[1])
            epsi1_imag_opt.append(epsi_ci)
            eq_det.append(det_2variables([res.x[0],res.x[1]]))
            
        cond_inicial = [res.x[0],res.x[1]]
        if epsi_ci==0:
            cond_inicial_mu.append([res.x[0],res.x[1]])
        
    if save_data_opt==1:
        print('Guardar data de minimizacion en .txt')
        if graficar==1:
            os.chdir(path_d)
        tabla = np.array([epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_det])
        tabla = np.transpose(tabla) 
        header1 = 'epsi_ci     Re(Lambda)      Im(Lambda)      Eq(det)' + info
        np.savetxt('opt_det_nano_mu%.4f_modo%i.txt' %(mu,modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)
  
    if graficar ==1:
        
        label_graph = 'Opt det, mu = %.4f' %(mu)
        labelx = r'$\epsilon_{ci}$'
        tamtitle = int(tamtitle*0.9)
        
        plt.figure(figsize=tamfig)
        plt.plot(epsi1_imag_opt,lambda_real_opt,'.g',ms = 10,label = label_graph)
        plt.plot([],[],'w',label = 'R = %.1f nm' %(R))
        plt.title(title+name_this_py,fontsize=tamtitle)
        plt.ylabel('Re($\Lambda$)',fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='lower right',markerscale=2,fontsize=tamlegend)
        plt.grid(1)
        if save_graphs==1:
            os.chdir(path_g)
            plt.savefig('Re_Lambda_mu%i_modo%i'%(mu*1e4,modo))
            
        if close_graphs==1:    
            plt.close()
            
        plt.figure(figsize=tamfig)
        plt.plot(epsi1_imag_opt,lambda_imag_opt,'.g',ms = 10,label = label_graph)
        plt.plot([],[],'w',label = 'R = %.1f nm' %(R))
        plt.title(title+name_this_py,fontsize=tamtitle)
        plt.ylabel('Im($\Lambda$)',fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.legend(loc='lower right',markerscale=2,fontsize=tamlegend)
        plt.grid(1)
        if save_graphs==1:
            plt.savefig('Im_Lambda_mu%i_modo%i'%(mu*1e4,modo))
    
        if close_graphs==1:    
            plt.close()
    
    j = j + 1

    del cond_inicial
    
#%%
