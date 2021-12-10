#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 14:48:39 2021

@author: leila

Purcell factor proporcional a lambda^3 * (Q/Veff) con Veff "Effective Volume" y Q "mode quality factor"

Q/Veff: cuantifica la fuerza de acople entre dipolo emisor y cavidad (polariton)

Problema con Qmodo = omega_k/(2*gamma_k) el gamma_k es el factor de damping del modo k-esimo (el medio 
                                                                                              interior es dielectrico, 
                                                                                              no tiene damping en principio)
Veff: integral de epsilon*|E(r)|^2 d^3r normalizada 

"""

import numpy as np
import sys
import os 
import scipy.integrate as integrate
import matplotlib.pyplot as plt

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real','') 

paper = 0
save_graphs = 1
normalizar = 0 #si normalizar = 1 entonces se normaliza Veff (tarda MUCHO)

#%%

err = 'fields_conkz.py no se encuentra en ' + path_basic
err2 = 'path de la carpeta donde se encuentra fields_conkz.py'
try:
    sys.path.insert(1, path_basic)
    from fields_conkz import Etot
except ModuleNotFoundError:
    print(err)
    path_basic = input(err2)
    sys.path.insert(1, path_basic)
    from fields_conkz import Etot
    
try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()


if save_graphs==1:
    path_g = path_basic + '/' + 'Purcell'
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)

#%%    

if paper == 1: 
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    labelpady = -2.5
    labelpadx = -0.5
    pad = 0.5
    lw = 1
else:
    tamfig = (10,8)
    tamlegend = 18
    tamletra = 18
    tamtitle = 18
    tamnum = 15
    labelpady = 0
    labelpadx = 0
    pad = 0
    lw = 1
    
#%%

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.05 #micrones
hbaramu = 0.3        #eV mu_c
modo = 1
Ao,Bo = 1,1

nmax = 5
ind = 99
#ind = 197

phi = 0
z = 0
gamma = 1 #factor de dumping

#%%

print('Importar datos de SPASER')

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/R_%.2f/mu_%.1f' %(re_epsi1,R,hbaramu) 
os.chdir(path_load)
name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)

try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data = np.transpose(data)
[list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data
kz = list_kz_opt[ind] #micrones
print('kz = ', kz)
crit = epsi1_imag_opt[ind]
omegac0 = omegac_opt[ind] 
epsi1 = re_epsi1 + 1j*crit

if Ao.imag == 0 and Bo.imag == 0:
    infAoBo = ', Ao = %i, Bo = %i' %(Ao,Bo)
elif Ao.imag != 0 and Bo.imag == 0:
    infAoBo = ', Ao = %i + i%i, Bo = %i' %(Ao.real,Ao.imag,Bo)
elif Ao.imag == 0 and Bo.imag != 0:
    infAoBo = ', Ao = %i, Bo = %i + i%i' %(Ao,Bo.real,Bo.imag)    
else:
    infAoBo = ', Ao = %i + i%i, Bo = %i + i%i' %(Ao.real,Ao.imag,Bo.real,Bo.imag)    
    
if Ao*Bo != 0:
    path_g = path_g + '/' + '2pol'
elif Bo == 0:
    path_g = path_g + '/' + 'polAo'
elif Ao == 0:
    path_g = path_g + '/' + 'polBo'
    
#%%

del list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det

#%%

info1 = 'kz = %.4f $\mu m^{-1}$, R = %.2f$\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,R,nmax,hbaramu) +infAoBo 
info2 = '$\epsilon_1$ = %.1f - i%.5e, $\omega/c$ = %.5e 1/$\mu$m, modo = %i, $\gamma$ = %.2f' %(re_epsi1,-crit,omegac0,modo,gamma)

inf_tot = info1 + ', ' + info2  + ', ' + name_this_py
title = info1 +'\n' + info2
labelomegac = '$\omega/c$ = %.2f$\mu m^{-1}$' %(omegac0)
labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely2 = 'Im($\epsilon_1$)'
labely1 = 'Purcell Factor'
inf_fig = '_modo%i_kz%.4f.png' %(modo,kz)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%

print('')
print('Calcular Purcell factor para diferentes Im(epsi1)')

tol = 1e-3
list_im_epsi1_fino = [0,crit+tol,crit,crit-tol,crit-2*tol]

N = int(5*1e2)               
omegac1,omegac2 = omegac0*0.98,omegac0*1.02
# lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
list_omegac = np.linspace(omegac1,omegac2,N)

list_Pfactor_tot = []
for im_epsi1 in list_im_epsi1_fino:
    im_epsi1 = np.round(im_epsi1,7)
    print(im_epsi1)
    list_P = []
    for omeggac in list_omegac:
        epsi1 = re_epsi1 + 1j*im_epsi1

        def Etotal(rho):
            [E1,E2] = Etot(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo,rho,phi,z)
            return E1
        
        result,err = integrate.quad(lambda rho: Etotal(rho), 0, R)
        
        Veff =  result*epsi1
        Q = omeggac*c/(2*gamma)
        cte = 0.75/(pi**2)
        lambdda = 2*pi/omeggac
        
        if normalizar == 1: 
            list_E = []
            list_r = np.linspace(1e-3,R,300)
            for r in list_r:
                list_E.append(Etotal(r))
        
            maxVeff = np.max(list_E)
            
        else:
            maxVeff = 1
            
        Veff = Veff/maxVeff
        
        P = cte*(lambdda**3)*Q/Veff  #purcell factor

        list_P.append(P)  
    list_Pfactor_tot.append(list_P)  

#%%

colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    
print('Graficar Purcell factor para diferentes Im(epsi1)')

graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
for j in range(len(list_Pfactor_tot)):
    list_P = np.abs(list_Pfactor_tot[j])
    im_epsi1 = list_im_epsi1_fino[j]
        
    if im_epsi1 == crit:
        labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
    elif im_epsi1 == 0:
        labell = 'Im($\epsilon_1$) = 0'  
    else:
        labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)

    plt.plot(list_omegac,list_P,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
n = 10
mini,maxi = np.min(list_Pfactor_tot),np.max(list_Pfactor_tot)
eje_Lambda2 = np.linspace(mini,maxi,n)
plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))

if save_graphs==1:
    os.chdir(path_g)
    # plt.tight_layout(1)
    namefig = 'Purcell_modo%i_kz%.4f' %(modo,kz)
    if normalizar == 1:
        namefig = namefig + '_norm.png'
    else:
        namefig = namefig + '.png'
    plt.savefig(namefig, format='png') 
    if paper == 0:
        np.savetxt('info_Purcell_modo%i_kz%.4f_zoom1D.txt' %(modo,kz), [inf_tot],fmt='%s')
    

#%%