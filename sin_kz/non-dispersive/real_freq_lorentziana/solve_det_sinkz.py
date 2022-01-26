#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 08:58:35 2020

@author: leila

graficar el log|det(sin kz)| en mapa de color en funcion de 
omega en THz y del potencial quimico mu

"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import minimize  
# from scipy.optimize import minimize_scalar

#%% 

save_data_opt = 1 #guardar data de la minimizacion
save_graphs = 1 #guardar los graficos

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16
labelpady = 0
labelpadx = 0
pad = 0
lw = 1
    
#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq_lorentziana','')
path_graphene = path_basic2.replace('/' + 'non-dispersive','') 

try:
    sys.path.insert(1, path_basic2)
    from det_sinkz import determinante
except ModuleNotFoundError:
    print('det_sinkz.py no se encuentra en ' + path_basic2)
    path_basic2 = input('path de la carpeta donde se encuentra det_sinkz.py')
    sys.path.insert(1, path_basic2)
    from det_sinkz import determinante

#freq en THz
try:
    sys.path.insert(1, path_basic2)
    from epsilon_lorentziana import epsilon 
except ModuleNotFoundError:
    print('epsilon_lorentziana.py no se encuentra en ' + path_basic2)
    path_basic2 = input('path de la carpeta donde se encuentra epsilon_lorentziana.py')
    sys.path.insert(1, path_basic2)
    from epsilon_lorentziana import epsilon

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

print('Definir parametros del problema')

R = 0.65              #micrones
eta = 0.95

#%%

# mu0 = 0.9

n = 250
aux_cte = c*1e-12 # freq en THz para epsilon
list_omegaTHz = np.linspace(1e-3,100,n)
list_mu = np.linspace(0.3,0.9,n)

list_modos = [1,2,3,4]
list_mu0 = [0.4,0.5,0.6]
cond_inicial = [25]   

# list_kz0 = np.linspace(2.5, 0.3, 221) # para mu0 fijo, iterar al reves por la bifurcacion
# cond_inicial = [175]   

path_save = path_basic + '/' +  r'eta_%.2f/R_%.2f' %(eta,R)

if save_data_opt==1:
    if not os.path.exists(path_save):
        print('Creating folder to save data')
        os.mkdir(path_save)

#%%

# if R not in [0.5, 0.25]:
#     raise TypeError('Mal el valor de R: El barrido en re(epsi1) se hizo para R = 0.5 o 0.25, y se usa como cond inicial')

#%%        

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%

for modo in list_modos:
    info1 = 'modo = %i, R = %.2f$\mu$m, $\eta$ = %.2f' %(modo,R,eta)  

#
    print('')
    print('Graficar |det| vs omega y mu para kz fijo para el modo = %i' %(modo))
           
    def det_2variables(omegaTHz,mu):
        fTHz = omegaTHz/(2*np.pi)
        epsi1 = epsilon(fTHz, eta)
        omegac = omegaTHz/aux_cte
        rta = determinante(omegac,epsi1,modo,R,mu)
        return np.log10(np.abs(rta))
    
    
    inf_tot =  info1 +  '\n' + name_this_py
    title = inf_tot
    labelx = '$\omega$ [THz]'
    labely2 = '$\mu_c$ [eV]'
    labely1 = '|det|'
    inf_fig = 'det_modo%i.png' %(modo)

#
    
    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    X, Y = np.meshgrid(list_omegaTHz, list_mu, sparse=True)
    f = np.vectorize(det_2variables)
    Z = f(X, Y)
    
    vmin,vmax = np.min(Z), np.max(Z)
    maxlog=int(np.ceil( np.log10( np.abs(vmax))))
    minlog=int(np.ceil( np.log10( np.abs(vmin))))
    
    if vmin < 0 :
          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                            + [0] 
                            + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+minlog+3)] )
    else:
          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
        
    pcm = plt.pcolormesh(X, Y, Z,
                          norm=colors.SymLogNorm(linthresh=0.5, linscale=0.5,
                                              vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    
    #    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(pcm, extend='both')
    # cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig(inf_fig, format='png') 
        # np.savetxt('info_det_modo%i_kz%.4f.txt' %(modo,kz0), [inf_tot],fmt='%s')
    
#
    for mu0 in list_mu0:
    
    #
        print('Minimizar |det| con mu = %.1f eV' %(mu0))
        
        tol_NM = 1e-10
        ite_NM = 4000
        
        # list_cond_init  = [43, 59, 73]
        omegaTHz_opt = []
        eq_det = []


        def det_1variable(x):
            omegaTHz = x
            
            fTHz = omegaTHz/(2*np.pi)
            epsi1 = epsilon(fTHz, eta)
            omegac = omegaTHz/aux_cte
            
            rta = determinante(omegac,epsi1,modo,R,mu0)
            return np.log10(np.abs(rta))
        
        res = minimize(det_1variable, cond_inicial, method='Nelder-Mead', tol=tol_NM, 
                          options={'maxiter':ite_NM})
     
        
        if res.message == 'Optimization terminated successfully.':
    
        # res = minimize_scalar(det_1variable, method='bounded', bracket=None,  bounds= bnds, tol=tol_NM)
    
        # if res.success == True:
    
            res.x = float(res.x)
            value = det_1variable(res.x)
            eq_det.append(value)
            omegaTHz_opt.append(res.x)
        cond_inicial = res.x
    
    #
        if save_data_opt == 1:
            os.chdir(path_save)
            print('Guardar data de minimizacion en .txt')
        
            tabla = np.array([omegaTHz_opt,eq_det])
            tabla = np.transpose(tabla)
            labeltxt = '_modo%i_mu%.4f.txt' %(modo,mu0)
            header1 = '     , $\mu_c$ = %.4f eV, ' %(mu0) + info1 +  ', ' + name_this_py
            header2 = 'Omega [THz]    log10(|det|)' + header1
            np.savetxt('opt_det_conkz_lorentziana' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header2)
    
    #
        graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
        plt.title(title,fontsize=int(tamtitle*0.9))
        # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
        
        plt.plot(list_mu0,omegaTHz_opt)
    
#%%
