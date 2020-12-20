#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)

determinante de 4x4

barrido NO en radio sino en 
kz. Find lambda que minimiza 
el determinante de 4x4 para diferentes 
valores de kz
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
path_basic2 = path_basic.replace('/' + 'complex_freq','')

try:
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante
except ModuleNotFoundError:
    print('det_conkz.py no se encuentra en el path_basic2 definido/carpeta de trabajo')
    path_basic2 = input('path de la carpeta donde se encuentra det_conkz.py')
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante

#%%

print('Definir parametros del problema')

re_epsi1 = 3.9
R = 0.5              #micrones
hbaramu = 0.3        #eV mu_c
modo = 1

if modo==1:
    list_im_epsi1 = np.linspace(0,-0.04,401) 
    list_im_epsi1 = np.linspace(0,-0.06,601) 
elif modo==2 or modo==3:
    list_im_epsi1 = np.linspace(0,-0.03,301) 
elif modo==4:
    list_im_epsi1 = np.linspace(0,-0.03,301) 
    
list_kz = np.linspace(0.1277,0.134,64) 
# list_kz = np.linspace(0.105,1.5,1395+1)
# list_kz = [0]

#%%

if list_im_epsi1[0]!=0:
    raise TypeError('Se necesita empezar con kz=0 para usar los datos de sin_kz, lo mismo para im(epsi1)')

if R!=0.5 or hbaramu!=0.3:
    raise TypeError('Wrong value for radium or for mu_c')
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')
    
#%%

print('Definir en donde vamos a guardar los datos de la minimizacion')

if save_graphs==1 or save_data_opt==1:
    try:
        path_det = r'/re_epsi1_%.2f/find_Lambda_vs_kz/modo_%i' %(re_epsi1,modo)
        path = path_basic + path_det
        os.chdir(path) 
    except OSError or IOError as error:
        print('Error en el path:', error)
        path = input('Copie y pegue la direc en la cual quiere guardar los graficos-archivos')
        os.chdir(path)

    path_g = path + '/' + 'graficos'
    path_d = path + '/' + 'archivos_txt'
    
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)
    if not os.path.exists(path_d):
        print('Creating folder to save data')
        os.mkdir(path_d)
        
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
    
    El barrido empieza con kz = 0 entonces puedo usar las condiciones iniciales
    halladas en la seccion1 sin_kz, barrido_re_epsi1. Estas condiciones
    iniciales son validas cuando R = 0.5 micrones, Im(epsi1)=0, kz=0
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

kz0 = list_kz[0]
if kz0 == 0:
    cond_inicial0 = fcond_inicial(re_epsi1)
else:
    cond_inicial0 = [4.686774571e+01, 2.004216225e-01] #kz anterior

if graficar==0:
    os.chdir(path_d)

cond_inicial_kz = []    

j = 0
for kz in list_kz:
    kz = np.round(kz,4)
    print('')
    print(kz)
    if kz == list_kz[0]:
        cond_inicial = cond_inicial0
    else:
        cond_inicial = cond_inicial_kz[j-1]
           
    epsi1_imag_opt = []
    lambda_real_opt = []
    lambda_imag_opt = []
    eq_det = []
    
    for im_epsi1 in list_im_epsi1:
        
        im_epsi1 = np.round(im_epsi1,4)
        
        epsi1 = re_epsi1 + 1j*im_epsi1
    
        def det_2variables(lmbd):
            [Re_lmbd,Im_lmbd] = lmbd
            lambdda = Re_lmbd+1j*Im_lmbd      
            omegac = 2*np.pi/lambdda
            rta = determinante(kz,omegac,epsi1,modo,R,hbaramu)
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
            cond_inicial_kz.append([res.x[0],res.x[1]])
        
    if save_data_opt==1:
        print('Guardar data de minimizacion en .txt')
        if graficar==1:
            os.chdir(path_d)
        tabla = np.array([epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_det])
        tabla = np.transpose(tabla)
        info = '.Opt det CON kz, R=%.1f \mum, mu_c=%.1f, Re(epsi1)=%.2f' %(R,hbaramu,re_epsi1) 
        header1 = 'Im(epsi1)     Re(Lambda)      Im(Lambda)      Eq(det)' + info + name_this_py
        np.savetxt('opt_det_kz%.4f_modo%i.txt' %(kz,modo), tabla, fmt='%1.9e', delimiter='\t', header = header1)
  
    if graficar ==1:
        
        label_graph = 'Opt det, kz = %.4f' %(kz)
        title = 'Re($\epsilon_1$) = %.2f, modo = %i, $\mu_c$ = %.1f eV' %(re_epsi1,modo,hbaramu)
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
            plt.savefig('Re_Lambda_kz%i_modo%i'%(kz*1e4,modo))
            
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
            plt.savefig('Im_Lambda_kz%i_modo%i'%(kz*1e4,modo))
    
        if close_graphs==1:    
            plt.close()
    
    j = j + 1

    del cond_inicial
    
#%%
