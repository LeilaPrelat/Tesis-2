#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 20:56:58 2020

@author: leila

-find_values_R_Ep_deg.py : encontrar todos los [R,Ep] para los que exista
 degeneracion, buscar en las carpetas y ver si hay un.txt que empiezan con
 info_degenerations_inv_dispR_%.2f.txt o info_degenerations_dispR_%.2f.txt 

"""
import os 

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

#%%

print('Definir parametros del problema')

listR = [0.1,0.2,0.3,0.4,0.5]
list_Ep = [0.7,0.4,0.5,0.3,0.8,0.6,0.9]
list_epsiinf_DL = [3.9,4.9]

gamma_DL = 0.01 #unidades de energia

#%%

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')

#%%

print('Ver para que valores hay degeneracion')

deg_values = []
for epsiinf_DL in list_epsiinf_DL:
    for R in listR:
        for Ep in list_Ep:
            path_load = path_basic + '/' + 'real_freq' + '/' + 'R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
            os.chdir(path_load)
            name1 = 'info_degenerations_dispR_%.2f.txt' %(R)
            name2 = 'info_degenerations_inv_dispR_%.2f.txt' %(R)
            deg = [R,epsiinf_DL,Ep]
            print('R, epsiinf_DL, Ep', deg)
            if os.path.isfile(path_load + '/' + name1) or os.path.isfile(path_load + '/' + name2): 
                deg_values.append(deg)
                print('deg found')

#%%
