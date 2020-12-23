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
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#%% 

save_data_opt = 1 #guardar data de la minimizacion

graficar = 0
save_graphs = 0
if graficar ==1:
    close_graphs = 0

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_det = path_basic.replace('/' + 'complex_freq','')

try:
    sys.path.insert(1, path_det)
    from det_sinkz import determinante
except ModuleNotFoundError:
    print('det_sinkz.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_det = input('path de la carpeta donde se encuentra det_sinkz.py')
    sys.path.insert(1, path_det)
    from det_sinkz import determinante

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

re_epsi1 = 4.9
R = 0.5              #micrones
modo = 4

if modo==1:
    list_im_epsi1 = np.linspace(0,-0.05,501) 
else:
    list_im_epsi1 = np.linspace(0,-0.04,401) 

list_mu =  np.linspace(0.3,0.9,601)  
# list_mu = [0.3]

#%%

# if list_im_epsi1[0]!= 0 or list_mu[0]!= 0.3:
#     raise TypeError('Se necesita empezar im(epsi1)=0 y mu = 0.3 para usar las cond iniciales')

if R!=0.5:
    raise TypeError('Wrong value for radium')
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_graphs==1 or save_data_opt==1:
    try:
        path_data = r'/re_epsi1_%.2f_vs_mu/find_Lambda' %(re_epsi1)
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

#%%

print('Definir las condiciones iniciales para el metodo de minimizacion')

def fcond_inicial(re_epsi1):
    """
    Parameters
    ----------
    re_epsi1 : parte real de la permeabilidad electrica del medio 1

    Returns
    -------
    [re(lambda), im(lambda)]
    
    Uso las condiciones iniciales halladas en la seccion1 sin_kz, barrido_re_epsi1. Estas condiciones
    iniciales son validas cuando R = 0.5 micrones, Im(epsi1)=0, kz=0 y mu = 0.3
    """
    os.chdir(path_basic)
    cond_init = np.loadtxt('cond_init_modo%i.txt' %(modo),delimiter='\t', skiprows = 1)
    cond_init = np.transpose(cond_init)
    [list_re_epsi1,list_re_lambda,list_im_lambda] = cond_init
    f1 = interp1d(list_re_epsi1,list_re_lambda)
    f2 = interp1d(list_re_epsi1,list_im_lambda)
    a = float(f1(re_epsi1))
    b = float(f2(re_epsi1))
    return [a,b]    

#%%        

print('Minimizacion del determinante de 4x4 para un barrido en kz')

mu0 = list_mu[0]
if mu0 == 0.3:
    cond_inicial0 = fcond_inicial(re_epsi1)
else:
    cond_inicial0 = [] #mu anterior

if graficar==0:
    os.chdir(path_d)

cond_inicial_mu = []    

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
    
    for im_epsi1 in list_im_epsi1:
        
        im_epsi1 = np.round(im_epsi1,4)
        
        epsi1 = re_epsi1 + 1j*im_epsi1
    
        def det_2variables(lmbd):
            [Re_lmbd,Im_lmbd] = lmbd
            lambdda = Re_lmbd + 1j*Im_lmbd      
            omegac = 2*np.pi/lambdda
            rta = determinante(omegac,epsi1,modo,R,mu)
            return np.abs(rta)
        
        res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=1e-11, 
                       options={'maxiter':1050})
#        print(res.message)
        if res.message == 'Optimization terminated successfully.':
            lambda_real_opt.append(res.x[0])
            lambda_imag_opt.append(res.x[1])
            epsi1_imag_opt.append(im_epsi1)
            eq_det.append(det_2variables([res.x[0],res.x[1]]))
            
        cond_inicial = [res.x[0],res.x[1]]
        if im_epsi1==0:
            cond_inicial_mu.append([res.x[0],res.x[1]])
        
    if save_data_opt==1:
        print('Guardar data de minimizacion en .txt')
        if graficar==1:
            os.chdir(path_d)
        tabla = np.array([epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_det])
        tabla = np.transpose(tabla)
        info = '.Opt det SIN kz, R=%.1f \mum, Re(epsi1)=%.2f' %(R,re_epsi1) 
        header1 = 'Im(epsi1)     Re(Lambda)      Im(Lambda)      Eq(det)' + info + name_this_py
        np.savetxt('opt_det_mu%.4f_modo%i.txt' %(mu,modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)
  
    if graficar ==1:
        
        label_graph = 'Opt det, mu = %.4f' %(mu)
        title = 'Re($\epsilon_1$) = %.2f, modo = %i' %(re_epsi1,modo)
        tamtitle = int(tamtitle*0.9)
        
        plt.figure(figsize=tamfig)
        plt.plot(epsi1_imag_opt,lambda_real_opt,'.g',ms = 10,label = label_graph)
        plt.plot([],[],'w',label = 'R = %.1f nm' %(R))
        plt.title(title+name_this_py,fontsize=tamtitle)
        plt.ylabel('Re($\Lambda$)',fontsize=tamletra)
        plt.xlabel(r'Im($\epsilon_1$)',fontsize=tamletra)
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
        plt.xlabel(r'Im($\epsilon_1$)',fontsize=tamletra)
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
