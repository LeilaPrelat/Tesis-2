#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)

determinante de 4x4

---> cython (no se grafica aca)

"""

import numpy as np
import os
import sys
from scipy.optimize import minimize   

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

try:
    sys.path.insert(1, path_basic)
    from det_conkz import determinante
except ModuleNotFoundError:
    print('det_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra det_conkz.py')
    sys.path.insert(1, path_basic)
    from det_conkz import determinante
#%%

def solve_det(kz,omegac,re_epsi1,list_im_epsi1,modo,R,hbaramu,cond_inicial):   
    """

    Parameters
    ----------
    kz_var : kz en 1/micrometros
    omegac : omega/c en 1/micrometros
    re_epsi1 : parte real de permeabilidad electrica del medio1
    list_im_epsi1 : barrido en im(epsi1) para el kz dado
    modo : modo
    R : radio del cilindro en micrometros
    mu_c : potencial quimico del grafeno
    cond_inicial: condicion inicial del barrido en im_epsi1 (que arranca de 0)

    Returns
    -------
    array de arrays: resultado de la minimizacion del |det| para un barrido en im(epsi1)
        devuelve : list(im_epsi1), lamba real, lambda imaginario, valor del |det(lambda)| 
    
    """   
    if list_im_epsi1[0]!=0:
        raise TypeError('Se necesita empezar con im(epsi1) = 0 para usar los datos de sin_kz y asi')
    
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
            rta = determinante(kz,omegac,epsi1,modo,R,hbaramu)
            return np.abs(rta)
        
        res = minimize(det_2variables, cond_inicial, method='Nelder-Mead', tol=1e-14, 
                       options={'maxiter':1050})
    #        print(res.message)
        if res.message == 'Optimization terminated successfully.':
            lambda_real_opt.append(res.x[0])
            lambda_imag_opt.append(res.x[1])
            epsi1_imag_opt.append(im_epsi1)
            eq_det.append(det_2variables([res.x[0],res.x[1]]))
            
        cond_inicial = [res.x[0],res.x[1]]

    return [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_det]

    
#%%
