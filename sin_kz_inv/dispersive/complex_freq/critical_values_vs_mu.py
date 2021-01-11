# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values_vs_mu.py: este .py toma los .txt de la carpeta hallar_Lambda
y obtiene los graficos que se guardaran en la carpeta find_critical_values (omega/c real y im_epsi1)

(ir de det_sinkz_vs_mu --> find_critical_values mediante critical_values_vs_mu.py)

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import minimize,fsolve,minimize_scalar,root_scalar

#%%

save_graphs = 1
save_data = 1

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_omegacQE = path_basic.replace('/' + 'complex_freq','')

try:
    sys.path.insert(1, path_omegacQE)
    from complex_omegac_QE import omegac_QE
except ModuleNotFoundError:
    print('complex_omegac_QE.py no se encuentra en ' + path_omegacQE)
    path_omegacQE = input('path de la carpeta donde se encuentra complex_omegac_QE.py')
    sys.path.insert(1, path_omegacQE)
    from complex_omegac_QE import omegac_QE 

#%%

print('Definir parametros para graficos')

tamfig = (10,8)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

#%%

print('Definir parametros del problema')

R = 0.4              #micrones
modo = 1

Ep = 0.9
epsiinf_DL = 3.9
gamma_DL = 0.01 #unidades de energia

barrido_mu = np.linspace(0.9,0.3,601)

info1 = 'R = %.1f $\mu$m, $E_p$ = %.3f eV, modo = %i' %(R,Ep,modo)
info2 = '$\epsilon_\infty$ = %.1f, $\gamma_{DL}$ = %.2f eV' %(epsiinf_DL,gamma_DL)
info =  ', ' + info1 + ', ' + info2  + ', ' + name_this_py
title = info1 +'\n' + info2  + ', ' + name_this_py

#%%

try:
    path0 = r'/epsiinf_DL_%.2f_vs_mu/R_%.2f/Ep_%.1f'  %(epsiinf_DL,R,Ep)
    path_load = path_basic + path0 + '/find_Lambda/modo_%i' %(modo)
    os.chdir(path_load)
except OSError or IOError:
    print('Los archivos .txt no se encuentran en ' + path_load)
    path_load = input('Copie y pegue la direc en donde se encuentran los archivos .txt del det_sinkz_vs_mu.py')    

if save_graphs == 1 or save_data == 1:
    try:
        path = path_basic + path0 + r'/find_critical_values' 
    except OSError or IOError:    
        print('No exist el path' + path + 'para guardar data/graficos')
        path = input('Copie y pegue la direc en la cual quiere guardar los graficos y los archivos .txt')

os.chdir(path_load)
        
#%%

if gamma_DL != 0.01:
    raise TypeError('Wrong value for gamma_DL')
    
if modo not in [1,2,3,4]:
    raise TypeError('Wrong value for mode')

#%%

tol_NM = 1e-13
ite_NM = 1150
N_inter = int(1e5) #puntos de la interpolacion

if R >= 0.5 :
    bounds = [-0.199,-0.13]
else:
    bounds = [-0.199,-0.05]

# ((min_first_var, min_second_var), (max_first_var, max_second_var))

rta_Im_epsi1 = []
rta_Im_epsi1_cuasi = [] #eq 16 del paper 2, aprox cuasi

rta_Im_omega_c = []
rta_Im_omega_c_cuasi = [] #eq 16 del paper 2, aprox cuasi

rta_Re_omega_c = []
rta_Re_omega_c_cuasi = [] #eq 16 del paper 2, aprox cuasi

if modo == 1:
    cond_inicial1 = -0.05
else:
    cond_inicial1 = -0.05
    
    
cond_inicial3 = cond_inicial1

for mu in barrido_mu:
    mu = np.round(mu,4)
    
    
    def Im_omegac_QE(epsi_ci): #eq 16 del paper 2, aprox cuasi
        rta = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,epsi_ci,R,mu)
        return rta.imag

    def Re_omegac_QE(epsi_ci): #eq 16 del paper 2, aprox cuasi
        rta = omegac_QE(modo,Ep,epsiinf_DL,gamma_DL,epsi_ci,R,mu)
        return rta.real


    tabla = np.loadtxt('opt_det_nano_mu%.4f_modo%i.txt' %(mu,modo), delimiter='\t', skiprows = 1)     
    tabla = np.transpose(tabla)
    try:
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt,eq_det] = tabla
    except ValueError:
        [epsi1_imag_opt,lambda_real_opt,lambda_imag_opt] = tabla
    
    omega_c_re = []
    omega_c_im = []
    for j in range(len(lambda_imag_opt)):
        lambbda = lambda_real_opt[j] + 1j*lambda_imag_opt[j]
        value = 2*np.pi/lambbda
        omega_c_re.append(value.real)
        omega_c_im.append(value.imag)
    
#    print('Interpolacion de los datos de minimizacion de gn para Im(omega/c) para hallar su cero')
    
    f2 = interp1d(epsi1_imag_opt,omega_c_im)
    ejex2 = np.linspace(np.min(epsi1_imag_opt),np.max(epsi1_imag_opt),N_inter)
    ejey2 = f2(ejex2)
    
    f3 = interp1d(epsi1_imag_opt,omega_c_re)
    ejex3 = ejex2
    ejey3 = f3(ejex3)
    

    # sol1 = fsolve(f2,cond_inicial1, xtol=1e-19,maxfev=1700)
    
    # sol1 = minimize(f2, cond_inicial1,method='Nelder-Mead', tol=tol_NM, 
    #         options={'maxiter':ite_NM})
    
    sol1 = fsolve(f2,cond_inicial3,xtol=1e-19,maxfev=1700)
    
    # sol1 = minimize_scalar(f2, bounds=bounds, method='bounded',tol = 1e-18)
        
    sol3 = fsolve(Im_omegac_QE,cond_inicial3,xtol=1e-19,maxfev=1700)
    
    print('mu = ', mu, 'ev')
    # print('Miminizar Im(omega/c):', sol1,f2(sol1[0]))   
    # rta_Im_epsi1.append(sol1[0])
    # rta_Im_omega_c.append(float(f2(sol1[0])))
    # rta_Re_omega_c.append(float(f3(sol1[0])))

    print('Miminizar Im(omega/c):', sol1[0],f2(sol1[0]))   
    rta_Im_epsi1.append(sol1[0])
    rta_Im_omega_c.append(f2(sol1[0]))
    rta_Re_omega_c.append(f3(sol1[0]))
    
    
    rta_Im_epsi1_cuasi.append(sol3[0])
    rta_Im_omega_c_cuasi.append(Im_omegac_QE(sol3[0]))   
    rta_Re_omega_c_cuasi.append(Re_omegac_QE(sol3[0]))   
    
    cond_inicial1 = sol1
    cond_inicial3 = sol3
    
del cond_inicial1
del cond_inicial3

#%%

print('Guardar data')

if save_data==1:
    os.chdir(path)
    tabla = np.array([barrido_mu,rta_Im_epsi1,rta_Im_omega_c,rta_Re_omega_c])
    tabla = np.transpose(tabla)
    header1 = 'mu [eV]     Im(epsi1)     Im(omega/c)      Re(omega/c)' + info
    np.savetxt('Minimizo_perdidas_nano_modo_%i.txt' %(modo), tabla, fmt='%1.11e', delimiter='\t', header = header1)

    tabla = np.array([barrido_mu,rta_Im_epsi1_cuasi,rta_Im_omega_c_cuasi,rta_Re_omega_c_cuasi])
    tabla = np.transpose(tabla)
    header1 = 'mu [eV]     Im(epsi1)     Im(omega/c)      Re(omega/c)' + info
    np.savetxt('Minimizo_perdidas_nano_modo_%i_QE.txt' %(modo), tabla, fmt='%1.11e', delimiter='\t', header = header1)

#%%

print('Importar datos obtenidos planteando frecuencia real')

path_load = path_omegacQE + '/real_freq/R_%.2f/epsiinf_DL_%.2f_vs_mu/Ep_%.1f' %(R,epsiinf_DL,Ep)
os.chdir(path_load)
opt_complex_freq = np.loadtxt('opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo),delimiter = '\t',skiprows=1)
opt_complex_freq = np.transpose(opt_complex_freq)
[barrido_mu2, Re_omegac_real_freq, Im_epsi1_real_freq, Eq] = opt_complex_freq


# os.chdir(path)
opt_complex_freq_QE = np.loadtxt('QE_lossless_vs_mu_modo%i.txt' %(modo),delimiter = '\t',skiprows=1)
opt_complex_freq_QE = np.transpose(opt_complex_freq_QE)
[barrido_mu3, Re_omegac_real_freq_QE, Im_epsi1_real_freq_QE] = opt_complex_freq_QE


#%%

label_graph = 'Opt complex freq'
label_graph2 = 'Opt real freq'
label_QE  = 'QE complex freq'
label_QE2  = 'QE real freq'

print('Graficar Im(epsi1) vs mu') 

plt.figure(figsize=tamfig)
plt.plot(barrido_mu2,Im_epsi1_real_freq,'-r',lw=5,alpha = 0.4,label=label_graph2)
plt.plot(barrido_mu,rta_Im_epsi1,'.-m',ms=10,label=label_graph)

plt.plot(barrido_mu3,Im_epsi1_real_freq_QE,'-b',lw=5,alpha = 0.4,label=label_QE2)
plt.plot(barrido_mu,rta_Im_epsi1_cuasi,'.-g',ms=10,label=label_QE)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'$\epsilon_{ci}$',fontsize=int(tamletra*1.2))
plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum) 
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)      
if save_graphs ==1:
    os.chdir(path)
    plt.savefig('Im_epsi1_vs_mu_modo%i'%(modo))
        
plt.figure(figsize=tamfig)
plt.plot(barrido_mu2,Re_omegac_real_freq,'-r',lw=5,alpha = 0.4,label=label_graph2)
plt.plot(barrido_mu,rta_Re_omega_c,'.-m',ms=10,label=label_graph)

plt.plot(barrido_mu3,Re_omegac_real_freq_QE,'-b',lw=5,alpha = 0.4,label=label_QE2)
plt.plot(barrido_mu,rta_Re_omega_c_cuasi,'.-g',ms=10,label=label_QE)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Re($\omega$/c)',fontsize=tamletra)
plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)      
if save_graphs ==1:
    plt.savefig('Re_omegac_vs_mu_modo%i'%(modo))

    
plt.figure(figsize=tamfig)
plt.plot(barrido_mu,rta_Im_omega_c,'.m',ms=10,label=label_graph)
plt.plot(barrido_mu,rta_Im_omega_c_cuasi,'.g',ms=10,label=label_QE)
plt.title(title,fontsize=tamtitle)
plt.ylabel(r'Im($\omega$/c)',fontsize=tamletra)
plt.xlabel('$\mu_c$ [eV]',fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
plt.grid(1)      
if save_graphs ==1:
    plt.savefig('Im_omegac_vs_mu_modo%i'%(modo))    
        
#%%
