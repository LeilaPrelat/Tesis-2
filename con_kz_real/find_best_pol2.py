#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

hallar la mejor polarizacion (la que maximiza los campos)
para el caso homogeneo vs barrido en kz

"""

import numpy as np
import sys
import os 

save_data_opt = 1 #guardar los graficos 2D del campo

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'best_pol2'

if save_data_opt==1:    
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

try:
    sys.path.insert(1, path_basic)
    from coef_matrix import coef
except ModuleNotFoundError:
    print('coef_matrix.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra coef_matrix.py')
    sys.path.insert(1, path_basic)
    from coef_matrix import coef

#print('Definir parametros para graficos')

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5              # micrones
hbaramu = 0.3        # eV mu_c

modo = 4
    
path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 

#%%

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')
    
if R!= 0.5:
    raise TypeError('Wrong value for radium')

#%%

print('')
print('Importar los valores de SPASER')

os.chdir(path_load)
name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)
try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)
print('')

data = np.transpose(data)
[list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data

coef_A1_re_tot = []
coef_C1_re_tot = []
coef_B2_re_tot = []
coef_D2_re_tot = []

coef_A1_im_tot = []
coef_C1_im_tot = []
coef_B2_im_tot = []
coef_D2_im_tot = []

for ind in range(len(list_kz_opt)):

    ind = int(ind)
    kz = list_kz_opt[ind] #micrones
    print('kz = ', kz)
    
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = re_epsi1 + 1j*im_epsi1
    
    coeff = coef(kz,omegac,epsi1,modo,R,hbaramu)
    coef_A1,coef_B2,coef_C1,coef_D2 = coeff
    coef_A1 = complex(coef_A1)
    coef_C1 = complex(coef_C1)
    coef_D2 = complex(coef_D2)
    coef_B2 = complex(coef_B2)

    coef_A1_re_tot.append(coef_A1.real)
    coef_C1_re_tot.append(coef_C1.real)
    coef_D2_re_tot.append(coef_D2.real)
    coef_B2_re_tot.append(coef_B2.real)
    
    coef_A1_im_tot.append(coef_A1.imag)
    coef_C1_im_tot.append(coef_C1.imag)
    coef_D2_im_tot.append(coef_D2.imag)
    coef_B2_im_tot.append(coef_B2.imag)

#%%

info = '.Opt det con kz, R=%.1f \mum, Re(epsi1)=%.2f, $\mu_c$ = %.3f eV' %(R,re_epsi1,hbaramu) 
allstr = 'vs_kz_modo%i.txt' %(modo)

if save_data_opt==1:
    os.chdir(path_save)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([coef_A1_re_tot,coef_A1_im_tot,list_kz_opt,omegac_opt,epsi1_imag_opt])
    tabla = np.transpose(tabla)
    header1 = 'Re(A1)     Im(A1)     kz [eV]     Omega/c [1/micrones]    Im(epsi1)' + info + ', ' + name_this_py
    np.savetxt('A1best_' + allstr, tabla, fmt='%1.9e', delimiter='\t', header = header1)

    tabla = np.array([coef_C1_re_tot,coef_C1_im_tot,list_kz_opt,omegac_opt,epsi1_imag_opt])
    tabla = np.transpose(tabla)
    header1 = 'Re(C1)     Im(C1)     kz [eV]     Omega/c [1/micrones]    Im(epsi1)' + info + ', ' + name_this_py
    np.savetxt('C1best_' + allstr, tabla, fmt='%1.9e', delimiter='\t', header = header1)

    tabla = np.array([coef_D2_re_tot,coef_D2_im_tot,list_kz_opt,omegac_opt,epsi1_imag_opt])
    tabla = np.transpose(tabla)
    header1 = 'Re(D2)     Im(D2)     kz [eV]     Omega/c [1/micrones]    Im(epsi1)' + info + ', ' + name_this_py
    np.savetxt('D2best_' + allstr, tabla, fmt='%1.9e', delimiter='\t', header = header1)

    tabla = np.array([coef_B2_re_tot,coef_B2_im_tot,list_kz_opt,omegac_opt,epsi1_imag_opt])
    tabla = np.transpose(tabla)
    header1 = 'Re(B2)     Im(B2)     kz [eV]     Omega/c [1/micrones]    Im(epsi1)' + info + ', ' + name_this_py
    np.savetxt('B2best_' + allstr, tabla, fmt='%1.9e', delimiter='\t', header = header1)
    
#%%
