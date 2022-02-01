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

data_freq_kz = 1   #importar datos a mu0,eta0 fijos.
data_freq_mu = 0   #importar datos a kz0,eta0 fijos. No funciona la minimizacion (no hay data)
data_freq_eta = 0  #importar datos a kz0,mu0 fijos. No funciona la minimizacion (no hay data)

graficar_2D = 0
graficar_3D = 1
paper = 0
save_graphs = 1 #guardar los graficos 

if data_freq_kz == 1 and graficar_3D == 1:
    graficar3D_freq_vs_kz = 0
    graficar3D_freq_vs_mu = 0
    graficar3D_freq_vs_eta = 1

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
    elif cross_section == 'Qext':
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
R = 0.5              #micrones
modo = 1
freq1 = 3.15 #en THz freq1 del medio activo
freq2 = 5.09 #en THz freq2 del medio activo
omegaTHz1 = freq1*2*np.pi
omegaTHz2 = freq2*2*np.pi
Ao, Bo = 1,1
nmax = 5
N = 250

infAoBo = r'$A_o$ = %i, $B_o$ = %i' %(Ao,Bo)
info = r'$\nu$ = %i, R = %.2f$\mu$m' %(modo,R) + ', ' + infAoBo
labelfreq = 'f = $\omega/2\pi$ [THz]'
tolkz = 0.1
tolmu = 0.005
tol_eta = 0.05
tolfreq = 0.1

#%% 1 
print('Importar datos de SPASER')

path_load = path_basic  + '/' +  'R_%.2f' %(R) 
os.chdir(path_load)

if data_freq_kz == 1:
    
    ################## cerca de freq2 ###########################
    [mu0,eta0,ind] = 0.6,0.8,60  # no funca muy bien cerca de freq2, nunca termina siendo el maximo de Qscat
    [mu0,eta0,ind] = 0.5,0.8,135 # no funca muy bien cerca de freq2, nunca termina siendo el maximo de Qscat
    [mu0,eta0,ind] = 0.4,0.8,207 # no funca muy bien cerca de freq2, nunca termina siendo el maximo de Qscat. (***)
    # [mu0,eta0,ind] = 0.6,0.85,60 
    # [mu0,eta0,ind] = 0.5,0.85,135 # (*)
    # [mu0,eta0,ind] = 0.6,0.9,61 
    # [mu0,eta0,ind] = 0.5,0.9,135
    # [mu0,eta0,ind] = 0.4,0.9,206
    # [mu0,eta0,ind] = 0.6,0.95,60 
    # [mu0,eta0,ind] = 0.5,0.95,134
    # [mu0,eta0,ind] = 0.4,0.95,205
    
    
    [mu0,eta0,ind] = 0.4,0.85,206  # este funca bien cambiando eta por 0.75 para el barrido en eta: (***)
#    [mu0,eta0,ind] = 0.4,0.66,206  # (***). cambiar por eta = 0.66
    
#    [mu0,eta0,ind] = 0.53,0.8,117 # match con (**), no poner eta0 = 0.74. Cambiar el freq0 a 5.39888 y el ind a 153 para el barrido en kz: (iv)
#    [mu0,eta0,ind] = 0.53,0.8,153 # (iv) cambiar el kz para que llegue al maximo de scattering
    
    
    ############################################################
    ################## cerca de freq1 ###########################
 #   [mu0,eta0,ind] = 0.2,0.8,117 # (**) # poner desp eta0 = 0.74 (estaba corrido el maximo de Qscat en la minimizacion) para el graf omega vs eta
#    [mu0,eta0,ind] = 0.2,0.85,117
#    [mu0,eta0,ind] = 0.2,0.9,90
#    [mu0,eta0,ind] = 0.2,0.95,117
#    [mu0,eta0,ind] = 0.1,0.85,135 # extra para que halla match con el de mu0 = 0.5 (*)
    
    
#    tolfreq = 0.05
    
    ############################################################
    
    labely = '$k_z$ [$\mu m^{-1}$]'
    inf = '$\mu_c$ = %.3f eV, $\eta$ = %.2f, ' %(mu0,eta0) + info + ', ' + name_this_py
    title = r'$\mu_c$ = %.3f eV, $\eta$ = %.2f' %(mu0,eta0) + '\n' + info 
    nametxt = 'opt_lorentz_conkzdet_modo%i_mu%.4f_eta%.2f' %(modo,mu0,eta0) + '.txt'
    
    
    try:
        data = np.loadtxt(nametxt,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(nametxt) if l.startswith('#')):
            print('values de ', nametxt, ':', line)
    except OSError or IOError:
        print('El archivo ' + nametxt + ' no se encuentra en ' + path_load)
    data = np.transpose(data)
    [row1,row2,row3] = data  
    kz0 = row1[ind] #micrones
    print('kz = ', kz0)
    omega0 = row2[ind]     
    freq0 = omega0/(2*np.pi)

    
    eta0 = 0.66  # estaba corrido el maximo de Qscat en la minimizacion para el caso [mu0,eta0,ind] = 0.2,0.8,117 
    # freq0 = 5.39888
    inf_fig = '_modo%i_mu%.4f_eta%.2f'  %(modo,mu0,eta0)  + '.png'
    
#%% 2

elif data_freq_mu == 1:
    eta0 = 0.9
    kz0 = 0.1   
    labely = '$\mu_c$ [eV]'
    inf = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f, ' %(kz0,eta0) + info + ', ' + name_this_py
    title = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f' %(kz0,eta0) + '\n' + info 
    nametxt = 'opt_lorentz_conkzdet_modo%i_kz%.4f_eta%.2f' %(modo,kz0,eta0) + '.txt'
    inf_fig = '_modo%i_kz%.4f_eta%.2f' %(modo,kz0,eta0) + '.png'
    
    try:
        data = np.loadtxt(nametxt,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(nametxt) if l.startswith('#')):
            print('values de ', nametxt, ':', line)
    except OSError or IOError:
        print('El archivo ' + nametxt + ' no se encuentra en ' + path_load)
    data = np.transpose(data)
    [row1,row2,row3] = data    
    mu0 = row1[ind] #micrones
    omega0 = row2[ind]
    freq0 = omega0/(2*np.pi)
    ejey1 = np.linspace(np.min(row1),np.max(row1),10)

    def Qscat3D(Freq,Mu): 
        Omega = 2*np.pi*Freq
        Omegac = Omega/aux_cte   
        
        Qscatt = Cross_Section(kz0,Omegac,eta0,nmax,R,Mu,Ao,Bo)
        return Qscatt   

#%% 3
    
elif data_freq_eta == 1:    
    mu0 = 0.6
    kz0 = 0.1
    labely = '$\eta$'
    inf = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$, ' %(mu0,kz0) + info + ', ' + name_this_py
    title = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$' %(mu0,kz0) + '\n' + info 
    nametxt = 'opt_lorentz_conkzdet_modo%i_kz%.4f_mu%.2f' %(modo,kz0,mu0) + '.txt'
    inf_fig = '_modo%i_kz%.4f_mu%.2f' %(modo,kz0,mu0) + '.png'
    
    try:
        data = np.loadtxt(nametxt,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(nametxt) if l.startswith('#')):
            print('values de ', nametxt, ':', line)
    except OSError or IOError:
        print('El archivo ' + nametxt + ' no se encuentra en ' + path_load)
    data = np.transpose(data)
    [row1,row2,row3] = data
    eta0 = row1[ind] #micrones
    omega0 = row2[ind]
    freq0 = omega0/(2*np.pi)
    ejey1 = np.linspace(np.min(row1),np.max(row1),10)
    
    def Qscat3D(Freq,Eta): 
        Omega = 2*np.pi*Freq
        Omegac = Omega/aux_cte   
        
        Qscatt = Cross_Section(kz0,Omegac,Eta,nmax,R,mu0,Ao,Bo)
        return Qscatt   
    
label_cross_section = cross_section + '$_{ad}$'
name = cross_section
ejey1 = np.linspace(np.min(row1),np.max(row1),10) # para marcar el valor critico con una recta
ejex1 = np.linspace(np.min(row2),np.max(row2),10) # para marcar el valor critico con una recta

#%%

del row1,row2,row3
 
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

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    return 

if save_graphs == 1:
    path_g = path_basic + '/' + 'seccion_eficaz_lorentziana/R_%.2f' %(R)
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)
    
#%%
if data_freq_kz == 1:
    if graficar_3D == 1:
    #from matplotlib.colors import SymLogNorm
        import matplotlib.colors as colors
        
        list_x = np.linspace(freq0 - tolfreq, freq0 + tolfreq, N) # eje x para el mapa de color
        if graficar3D_freq_vs_kz == 1: 
            print('')
            print('Graficar '+ name + ' en 3D (freq-kz)') 
            list_y = np.linspace(kz0 - tolkz, kz0 + tolkz, N) # kz    
            def Qscat3D(Freq,Kz): 
                Omega = 2*np.pi*Freq
                Omegac = Omega/aux_cte   
                
                Qscatt = Cross_Section(Kz,Omegac,eta0,nmax,R,mu0,Ao,Bo)
                return Qscatt
            
            # print(list_x[0],list_x[-1])
            
            X, Y = np.meshgrid(list_x, list_y, sparse=True)
            f = np.vectorize(Qscat3D)
            Z = f(X, Y)
                
            limits = [np.min(list_x) , np.max(list_x), np.min(list_y) , np.max(list_y)]
            graph(title,labelfreq,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
            plt.title(title,fontsize=int(tamtitle*0.9))
            vmin, vmax = np.min(Z), np.max(Z)
            
            
            pcm = plt.pcolormesh(X, Y, Z,
                                 norm=colors.SymLogNorm(linthresh=0.5, linscale=0.5,
                                                        vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
                    
            cbar = plt.colorbar(pcm, extend='both')
            cbar.ax.tick_params(labelsize=tamnum)
            cbar.set_label(name,fontsize=tamlegend)
            # plt.plot(freq0*np.ones(N),list_y,'-k',lw = 2.5)
            # plt.plot(list_x,kz0*np.ones(N),'-k',lw = 2.5)
            if save_graphs==1:
                os.chdir(path_g)
                plt.tight_layout()
                plt.savefig(name + '3D' + inf_fig, format='png') 
                
        if graficar3D_freq_vs_mu == 1:
            
            labely = '$\mu_c$ [eV]'
            inf = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f, ' %(kz0,eta0) + info + ', ' + name_this_py
            title = 'kz = %.4f $\mu m^{-1}$, $\eta$ = %.2f' %(kz0,eta0) + '\n' + info 
            nametxt = 'opt_lorentz_conkzdet_modo%i_kz%.4f_eta%.2f' %(modo,kz0,eta0) + '.txt'
            inf_fig = '_modo%i_kz%.4f_eta%.2f' %(modo,kz0,eta0) + '.png'
            print('')
            print('Graficar '+ name + ' en 3D (freq-mu)') 
            list_y = np.linspace(mu0 - tolmu, mu0 + tolmu, N) # kz    
            def Qscat3D(Freq,Mu): 
                Omega = 2*np.pi*Freq
                Omegac = Omega/aux_cte   
                
                Qscatt = Cross_Section(kz0,Omegac,eta0,nmax,R,Mu,Ao,Bo)
                return Qscatt
            
            # print(list_x[0],list_x[-1])
            
            X, Y = np.meshgrid(list_x, list_y, sparse=True)
            f = np.vectorize(Qscat3D)
            Z = f(X, Y)
                
            limits = [np.min(list_x) , np.max(list_x), np.min(list_y) , np.max(list_y)]
            graph(title,labelfreq,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
            plt.title(title,fontsize=int(tamtitle*0.9))
            vmin, vmax = np.min(Z), np.max(Z)
            
            
            pcm = plt.pcolormesh(X, Y, Z,
                                 norm=colors.SymLogNorm(linthresh=0.5, linscale=0.5,
                                                        vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
                    
            cbar = plt.colorbar(pcm, extend='both')
            cbar.ax.tick_params(labelsize=tamnum)
            cbar.set_label(name,fontsize=tamlegend)
            # plt.plot(freq0*np.ones(N),list_y,'-k',lw = 2.5)
            # plt.plot(list_x,mu0*np.ones(N),'-k',lw = 2.5)
            if save_graphs==1:
                os.chdir(path_g)
                plt.tight_layout()
                plt.savefig(name + '3D' + inf_fig, format='png')             

        if graficar3D_freq_vs_eta == 1:
            

            labely = '$\eta$'
            inf = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$, ' %(mu0,kz0) + info + ', ' + name_this_py
            title = '$\mu_c$ = %.3f eV, kz = %.4f $\mu m^{-1}$' %(mu0,kz0) + '\n' + info 
            nametxt = 'opt_lorentz_conkzdet_modo%i_kz%.4f_mu%.2f' %(modo,kz0,mu0) + '.txt'
            inf_fig = '_modo%i_kz%.4f_mu%.2f' %(modo,kz0,mu0) + '.png'
            print('')
            print('Graficar '+ name + ' en 3D (freq-eta)') 
            list_y = np.linspace(eta0 - tol_eta, eta0 + tol_eta, N) # kz    
            def Qscat3D(Freq,Eta): 
                Omega = 2*np.pi*Freq
                Omegac = Omega/aux_cte   
                
                Qscatt = Cross_Section(kz0,Omegac,Eta,nmax,R,mu0,Ao,Bo)
                return Qscatt
            
            # print(list_x[0],list_x[-1])
            
            X, Y = np.meshgrid(list_x, list_y, sparse=True)
            f = np.vectorize(Qscat3D)
            Z = f(X, Y)
                
            limits = [np.min(list_x) , np.max(list_x), np.min(list_y) , np.max(list_y)]
            graph(title,labelfreq,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
            plt.title(title,fontsize=int(tamtitle*0.9))
            vmin, vmax = np.min(Z), np.max(Z)
            
            
            pcm = plt.pcolormesh(X, Y, Z,
                                 norm=colors.SymLogNorm(linthresh=0.5, linscale=0.5,
                                                        vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
                    
            cbar = plt.colorbar(pcm, extend='both')
            cbar.ax.tick_params(labelsize=tamnum)
            cbar.set_label(name,fontsize=tamlegend)
            # plt.plot(freq0*np.ones(N),list_y,'-k',lw = 2.5)
            # plt.plot(list_x,eta0*np.ones(N),'-k',lw = 2.5)
            if save_graphs==1:
                os.chdir(path_g)
                plt.tight_layout()
                plt.savefig(name + '3D' + inf_fig, format='png')      
        
#%%
if data_freq_kz == 1:
    if graficar_2D == 1:
        print('')
        print('Calcular '+ name + ' para diferentes kz')
        
        N = int(5*1e3)               
        # lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
        list_freq = np.linspace(freq1-1,freq2+1,N)
        list_kz = [kz0 - 0.2, kz0 - 0.1, kz0, kz0 + 0.1, kz0 + 0.2]
        list_Qabs_tot1 = []
        for kz in list_kz:
            kz = np.round(kz,4)
            list_Qabs1 = []
            for freqq in list_freq:
                omegga = freqq*2*np.pi
                omeggac = omegga/aux_cte
                Qabss = Cross_Section(kz,omeggac,eta0,nmax,R,mu0,Ao,Bo)
                list_Qabs1.append(Qabss)  
            list_Qabs_tot1.append(list_Qabs1)  
            if kz == kz0:
                arg = np.argmax(list_Qabs1)
                omega_maximum = list_freq[arg]           
        del list_Qabs1
                        
        colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    
        # print('Graficar ' + name  + ' para diferentes Im(epsi1)')
        
        graph(title,labelfreq,label_cross_section,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
        for j in range(len(list_Qabs_tot1)):
            list_Qabs2 = list_Qabs_tot1[j]  # no tomar valor absoluto porque hay valores de Qscat negativos !!! 
            kz = list_kz[j]
                
            labell = '$k_z$ = %.3f 1/$\mu$m' %(kz)
        
            plt.plot(list_freq,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
        
        n = 10
        mini,maxi = np.min(list_Qabs_tot1),np.max(list_Qabs_tot1)
        # print(np.abs(omega_maximum-omega0))
        eje_Lambda2 = np.linspace(mini,maxi,n)
        plt.plot(freq0*np.ones(n),eje_Lambda2,'-r',lw = 1)
        plt.plot(freq1*np.ones(n),eje_Lambda2,'-k',lw = 1)
        plt.plot(freq2*np.ones(n),eje_Lambda2,'-k',lw = 1)
        # plt.yscale('log')
        plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
        
        if save_graphs==1:
            os.chdir(path_g)
            # plt.tight_layout(1)
            plt.savefig(name + '2D' + inf_fig, format='png') 
        
#%%
