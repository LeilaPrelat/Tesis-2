#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  11 22:07:37 2021

@author: leila

graficar Qscat Qabs o Qext (elegir en la variable cross_section)

ESTABAN MAL ESAS FORMULAS DEL CUADERNO (movidas a extra)
VER TESIS
(revisar cuentas)

epsilon como funcion de omega: eq 23 del paper 
passarelli2019
modelo lorentziana para el epsilon pero con los resultados
del modelo simple no dispersivo (real_freq)

plot_cross_section_lorentziana.py : usa las soluciones del epsilon1(omega) de la carpeta 
real_freq_lorentziana
"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 
graph_2D_freq_mu = 0   #graficos mapa de color en el plano omega-mu
graph_2D_freq_eta = 1    #graficos mapa de color en el plano eta- freq (omega/2*np.pi)
zoom = 0        #graficos 2D con o sin zoom
graph_1D = 0  #graficos 1D

if graph_2D_freq_eta == 1 or graph_2D_freq_mu == 1:
    paper = 1
else:
    paper = 0

list_cross_section = ['Qscat', 'Qabs', 'Qext'] 
cross_section = list_cross_section[0]

#%% 

#print('importar funciones para graficarlas')

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/' + 'con_kz_real/real_freq_lorentziana','') 


#print('Importar modulos necesarios para este codigo')
err = 'cross_sections_conkz_lorentziana.py no se encuentra en el path_basic definido/carpeta de trabajo'
err2 = 'path de la carpeta donde se encuentra cross_sections_conkz_lorentziana.py'

try:
    sys.path.insert(1, path_basic)
    if cross_section == 'Qscat':
        from cross_sections_conkz_lorentziana import Qscat as cross_section_f
    elif cross_section == 'Qabs':
        from cross_sections_conkz_lorentziana import Qabs as cross_section_f
    else:
        from cross_sections_conkz_lorentziana import Qext as cross_section_f
except ModuleNotFoundError:
    print(err)

def Cross_Section(kz,omeggac,eta,nmax,R,hbaramu,Ao,Bo):
    return cross_section_f(kz,omeggac,eta,nmax,R,hbaramu,Ao,Bo)

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()
aux_cte = c*1e-12

#%%   

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
eta0 = 0.95  # population invertion para el epsilon(omega) lorentziano
R = 0.5    #micrones
mu0 = 0.3       #eV mu_c
modo = 2

Ao, Bo = 1,1
nmax = 5

ind = 0

if save_graphs == 1:
    path_g = path_basic + '/' + 'seccion_eficaz_lorentziana/R_%.2f' %(R)
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)
        
#%%

#print('Definir parametros para graficos')

if paper == 1: 
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    labelpady = -1.5
    labelpadx = 0.5
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

print('Importar datos de SPASER')

path_load = path_basic  + '/' +  'eta_%.2f/R_%.2f' %(eta0,R) 
os.chdir(path_load)
name = 'opt_det_conkz_lorentziana_modo%i_mu%.4f.txt' %(modo,mu0)

try:
    data = np.loadtxt(name,delimiter = '\t', skiprows=1)
    for line in (l.strip() for l in open(name) if l.startswith('#')):
        print('values de ', name, ':', line)
except OSError or IOError:
    print('El archivo ' + name + ' no se encuentra en ' + path_load)

data = np.transpose(data)
[list_kz_opt,omegaTHz_opt,eq_det] = data
kz0 = list_kz_opt[ind] #micrones
omega0 = omegaTHz_opt[ind]
print('kz = ', kz0)
#epsi1 = re_epsi1 + 1j*crit

tol = 1e-1
list_eta = [eta0 - tol, eta0, eta0 + tol]
list_eta = [0.9,0.95,1]

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

del list_kz_opt,omegaTHz_opt,eq_det
 
info1 = r'kz = %.4f$\mu m^{-1}$, R = %.2f$\mu$m, $\mu_c$ = %.4feV, $\nu$ = %i' %(kz0,R,mu0,modo) + infAoBo 
info2 = r'$\eta$ = %.4f, nmax = %i' %(eta0,nmax) + ', ' + name_this_py

inf_tot = info1 + ', ' + info2  
title = info1 + '\n' + info2 
labelomegac = '$\omega$ = %.2f THz' %(omega0)
labelOmega = '$\omega$ [THz]'
labelkz = '$k_z$ [1/$\mu$m]'
labelmu = '$\mu_c$ [eV]'
labelEta = '$\eta$'
labelfreq = 'f = $\omega/2\pi$'

label_cross_section = cross_section + '$_{ad}$'
name = cross_section
inf_fig = '_modo%i_kz%.4f_mu%.2f_eta%.2f' %(modo,kz0,mu0,eta0)
if zoom == 1: 
    inf_fig = '_modo%i_kz%.4f_zoom.png' %(modo,kz0)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%

"""
if re_epsi1 != 3.9:
    raise TypeError('Wrong value for Re(epsi1): Re(epsi1) = 3.9')
    
if R != 0.5:
    raise TypeError('Wrong value for R: R = 0.5 micrones')
    
if hbaramu != 0.3:
    raise TypeError('Wrong value for chemical potential: mu = 0.3')
"""
#%%

if graph_1D==1:
    print('')
    print('Calcular '+ name + ' para diferentes $\eta$')
    
    N = int(5*1e3)               
    omega1,omega2 = omega0-3, omega0+3
    # lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
    list_omega = np.linspace(omega1,omega2,N)
    
    list_Qabs_tot1 = []
    for eta in list_eta:
        eta = np.round(eta,3)
        print(eta)
        list_Qabs1 = []
        for omegga in list_omega:
            omeggac = omegga/aux_cte
            Qabss = Cross_Section(kz0,omeggac,eta,nmax,R,mu0,Ao,Bo)
            list_Qabs1.append(Qabss)  
        list_Qabs_tot1.append(list_Qabs1)  
        if eta == eta0:
            arg = np.argmax(list_Qabs1)
            omega_maximum = list_omega[arg]           
    del list_Qabs1
                    
    colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    
    # print('Graficar ' + name  + ' para diferentes Im(epsi1)')
    
    graph(title,labelOmega,label_cross_section,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
    for j in range(len(list_Qabs_tot1)):
        list_Qabs2 = list_Qabs_tot1[j]  # no tomar valor absoluto porque hay valores de Qscat negativos !!! 
        eta = list_eta[j]
            
        labell = '$\eta$ = %.3f' %(eta)
    
        plt.plot(list_omega,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qabs_tot1),np.max(list_Qabs_tot1)
    print(np.abs(omega_maximum-omega0))
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omega0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    
    # plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    
    if save_graphs==1:
        os.chdir(path_g)
        # plt.tight_layout(1)
        plt.savefig(name + '_fino_modo%i_kz%.4f.png' %(modo,kz0), format='png') 
        if paper == 0:
            np.savetxt('info_' + name + '_modo%i_kz%.4f_zoom1D.txt' %(modo,kz0), [inf_tot],fmt='%s')
    
#%%

#from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors

if graph_2D_freq_mu ==1:
    print('Graficar '+ name + ' en 2D') 
    
    # def Qscat2D(Mu,Kz): 
    #     omegac0 = omega0/aux_cte         
    #     Qscatt = Cross_Section(Kz,omegac0,eta0,nmax,R,Mu,Ao,Bo)
    #     return Qscatt

    # kz0 = 0.5325
    def Qscat2D(Freq,Mu): 
        Omega = 2*np.pi*Freq
        Omegac = Omega/aux_cte   
        
        Qscatt = Cross_Section(kz0,Omegac,eta0,nmax,R,Mu,Ao,Bo)
        return Qscatt
        

    N = 250
    
    if zoom==1:
        tol = 1e-3
    else:
        tol = 1e-2
        # tol = 5*1e-1 #para ver el otro polo
    
    omega_min,omega_max = omega0*(1-tol),omega0*(1+tol)
    list_omega = np.linspace(omega_min,omega_max,N)
    delta = (omega_max-omega_min)*0.5
    # print(omega_min,omega_max)

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    list_omega = np.linspace(omega0 - 0.5, omega0 + 0.5,N)
    
    freq0 = omega0/(2*np.pi)
    list_x = np.linspace(freq0 - 1.5, freq0 + 3.5,N) # mu
    
    mu11,mu22 = mu0 - 5*1e-2, mu0 + 5*1e-2
    list_y = np.linspace(np.round(mu11,2), np.round(mu22,2) , N) # kz
    
#    list_y = np.linspace(0.3,0.9,N)

    # omega0 = list_x[195]
    # mu0 =  list_y[195]
    # list_x = np.linspace(omega0 - 5*1e-1, omega0 + 5*1e-1,N) # mu
    # list_y = np.linspace(mu0 - 5*1e-2, mu0 + 5*1e-2,N) # kz

    print(list_x[0],list_x[-1])

    X, Y = np.meshgrid(list_x, list_y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    limits = [np.min(list_x) , np.max(list_x), np.min(list_y) , np.max(list_y)]
    graph(title,labelfreq,labelmu,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    vmin, vmax = np.min(Z), np.max(Z)
    
    # plt.plot(list_x,np.ones(N)*mu0,'--',lw = lw,color = 'green')
    # plt.plot(np.ones(N)*omega0,list_y,'--', lw = lw,color = 'green')
    
    pcm = plt.pcolormesh(X, Y, Z,
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=vmin, vmax=vmax),   
    				cmap='RdBu')
    
    
    cbar = plt.colorbar(pcm, extend='both')
    cbar.ax.tick_params(labelsize=tamnum)
    if paper == 0:
         cbar.set_label(name,fontsize=tamlegend)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_g)
        if paper == 1:
            plt.tight_layout()
        plt.savefig(name + '2D_v1' + inf_fig + '.png', format='png') 
        if paper == 1:
            np.savetxt('info_' + name + inf_fig + '.txt', [inf_tot],fmt='%s')

#%%

if graph_2D_freq_eta ==1:
    print('Graficar '+ name + ' en 2D') 

    # kz0 = 0.5325
    def Qscat2D(Freq,Eta): 
        Omega = 2*np.pi*Freq
        Omegac = Omega/aux_cte       
        Qscatt = Cross_Section(kz0,Omegac,Eta,nmax,R,mu0,Ao,Bo) 
        return Qscatt
        

    N = 250
    
    if zoom==1:
        tol = 1e-3
    else:
        tol = 1e-2
        # tol = 5*1e-1 #para ver el otro polo
    
    omega_min,omega_max = omega0*(1-tol),omega0*(1+tol)
    list_omega = np.linspace(omega_min,omega_max,N)
    delta = (omega_max-omega_min)*0.5
    # print(omega_min,omega_max)

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    list_omega = np.linspace(omega0 - 0.5, omega0 + 0.5,N)
    
    freq0 = omega0/(2*np.pi)
    list_x = np.linspace(freq0 - 1, freq0 + 3.5,N) # mu

    list_y = np.linspace(eta0 - 5*1e-1, eta0 + 5*1e-2,N) # kz
#    list_y = np.linspace(0.3,0.9,N)

    # omega0 = list_x[195]
    # mu0 =  list_y[195]
    # list_x = np.linspace(omega0 - 5*1e-1, omega0 + 5*1e-1,N) # mu
    # list_y = np.linspace(mu0 - 5*1e-2, mu0 + 5*1e-2,N) # kz

    print(list_x[0],list_x[-1])

    X, Y = np.meshgrid(list_x, list_y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    limits = [np.min(list_x) , np.max(list_x), np.min(list_y) , np.max(list_y)]
    graph(title,labelfreq,labelEta,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    vmin, vmax = np.min(Z), np.max(Z)
    
    # plt.plot(list_x,np.ones(N)*mu0,'--',lw = lw,color = 'green')
    # plt.plot(np.ones(N)*omega0,list_y,'--', lw = lw,color = 'green')
    
    pcm = plt.pcolormesh(X, Y, Z,
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=vmin, vmax=vmax),   
    				cmap='RdBu')
    
    
    cbar = plt.colorbar(pcm, extend='both')
    cbar.ax.tick_params(labelsize=tamnum)
    if paper == 0:
         cbar.set_label(name,fontsize=tamlegend)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_g)
        if paper == 1:
            plt.tight_layout()
        plt.savefig(name + '2D_v2' + inf_fig + '.png', format='png') 
        if paper == 1:
            np.savetxt('info_' + name + inf_fig + '.txt', [inf_tot],fmt='%s')
            
#%%            
            
            
