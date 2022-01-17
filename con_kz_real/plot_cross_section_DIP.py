#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

graficar Qscat Qabs o Qext (elegir en la variable cross_section)

ESTABAN MAL ESAS FORMULAS DEL CUADERNO (movidas a extra)
VER TESIS
(revisar cuentas)

con un dipolo en el interior

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 
graph_2D = 1    #graficos 2D
zoom = 0        #graficos 2D con o sin zoom
graph_1D = 0   #graficos 1D

if graph_2D == 1:
    paper = 1
else:
    paper = 0

list_cross_section = ['Qscat', 'Qabs', 'Qext'] 
cross_section = list_cross_section[0]

#%%

#print('Definir parametros para graficos')

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

#print('importar funciones para graficarlas')

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

if save_graphs==1:
    path_g = path_basic + '/' + 'seccion_eficaz_DIP'
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)

#print('Importar modulos necesarios para este codigo')
err = 'cross_sections_DIP_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo'
err2 = 'path de la carpeta donde se encuentra cross_sections_DIP_conkz.py'
if cross_section == 'Qscat':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qscat
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qscat
        
    def Cross_Section(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        return Qscat(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)

elif cross_section == 'Qabs':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qabs
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qabs
        
    def Cross_Section(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        return Qabs(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)

else: 

    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qext
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz_DIP import Qext
        
    def Cross_Section(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
        return Qext(kz_var,omegac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)

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
ind = -1

### datos del dipolo ####
theta_D = 0
z_D = 0
rho_D = R/2

px = 1/2
py = -0.5*1j

# px = 0
# py = 0
pz = 0
p1 = px + 1j*py # p+
p2 = px - 1j*py # p-

#%%

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    return 

def labelp(pvalue):
    if pvalue.imag > 0:
        if pvalue.real > 0:
            rta = ' = %.2f + i%.2f' %(pvalue.real, pvalue.imag)
        elif pvalue.real < 0:
            rta = ' = -%.2f + i%.2f' %(-pvalue.real, pvalue.imag)
        else:
            rta = ' = i%i' %(pvalue.imag)
    elif pvalue.imag < 0:
        if pvalue.real > 0:
            rta = ' = %.2f - i%.2f' %(pvalue.real, -pvalue.imag)
        elif pvalue.real < 0:
            rta = ' = -%.2f - i%.2f' %(-pvalue.real, -pvalue.imag)
        else:
            rta = ' = - i%.2f' %(-pvalue.imag)
    else:
        if pvalue.real>0:
            rta = ' = %.2f' %(pvalue.real)
        elif pvalue.real<0:
            rta = ' = -%.2f' %(-pvalue.real)
        else:
            rta = ' = 0'
    return rta

#%%

print('Importar datos de SPASER')

path_load = path_basic  + '/' + 'real_freq' + '/' +  're_epsi1_%.2f_vs_kz/R_%.2f/mu_%.1f' %(re_epsi1,R,hbaramu) 
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

labelpx = 'px' + labelp(px)
labelpy = 'py' + labelp(py)

labelp1 = 'p+' + labelp(p1)
labelp2 = 'p-' + labelp(p2)
    
#%%

del list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det

info1 = r'Ao = %i, Bo = %i, kz = %.4f 1/$\mu$m' %(Ao,Bo,kz)
info2 = r'R = %.2f$\mu$m, nmax = %i, $\mu_c$ = %.1f eV, $\nu$ = %i' %(R,nmax,hbaramu,modo)
info3 = r'$\epsilon_1$ = %.1f - i%.4e, $\omega/c$ = %.4e 1/$\mu$m, $\rho_D$ = %.3f$\mu$m' %(re_epsi1,-crit,omegac0,rho_D)
infoDIP = labelpx + ', ' + labelpy + r', pz = %i, $\theta_D$ = %i, $z_D$ = %i$\mu$m' %(pz,theta_D,z_D)
title =  info1 + '\n' + info2 + '\n' + info3 + '\n' + infoDIP
inf_tot = info1 + ', ' + info2 + ', ' + info3 + ', ' + infoDIP

labelomegac = '$\omega/c$ = %.2f$\mu m^{-1}$' %(omegac0)
labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely2 = 'Im($\epsilon_1$)'
labely1 = cross_section + '$_{ad}$'
name = cross_section
inf_fig = '_modo%i_kz%.4f.png' %(modo,kz)
if zoom == 1: 
    inf_fig = '_modo%i_kz%.4f_zoom.png' %(modo,kz)
colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    

#%%
    
if re_epsi1 != 3.9:
    raise TypeError('Wrong value for Re(epsi1): Re(epsi1) = 3.9')
    
if R not in [0.5,0.05]:
    raise TypeError('Wrong value for radium')
    
if hbaramu != 0.3:
    raise TypeError('Wrong value for chemical potential: mu = 0.3')

#%%

if graph_1D==1:
    print('')
    print('Calcular '+ name + ' para diferentes Im(epsi1)')
    
    tol = 1e-3
    list_im_epsi1_fino = [0,crit+3*tol,crit+2*tol,crit+tol,crit,crit-tol,crit-2*tol]
    list_im_epsi1_grueso = [0,-crit,0.5,-0.001,-0.01,crit,-0.5]
    
    N = int(5*1e3)               
    omegac1,omegac2 = omegac0*0.98,omegac0*1.02
    # lambda1,lambda2 = lambbda_real*0.999998,lambbda_real*1.000002
    list_omegac = np.linspace(omegac1,omegac2,N)
    
    list_Qabs_tot1 = []
    for im_epsi1 in list_im_epsi1_fino:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qabs1 = []
        for omeggac in list_omegac:
            epsi1 = re_epsi1 + 1j*im_epsi1
            Qabss = Cross_Section(kz,omeggac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
            list_Qabs1.append(Qabss)  
        list_Qabs_tot1.append(list_Qabs1)  
    
    del list_Qabs1
                
    list_Qabs_tot2 = []
    for im_epsi1 in list_im_epsi1_grueso:
        im_epsi1 = np.round(im_epsi1,7)
        print(im_epsi1)
        list_Qabs2 = []
        for omeggac in list_omegac:
            epsi1 = re_epsi1 + 1j*im_epsi1
            Qabss = Cross_Section(kz,omeggac,epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
            list_Qabs2.append(Qabss)  
        list_Qabs_tot2.append(list_Qabs2)  
    
    del list_Qabs2
    
    print('Graficar ' + name  + ' para diferentes Im(epsi1)')
    
    graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
    for j in range(len(list_Qabs_tot1)):
        list_Qabs2 = np.abs(list_Qabs_tot1[j])
        im_epsi1 = list_im_epsi1_fino[j]
            
        if im_epsi1 == crit:
            labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
        elif im_epsi1 == 0:
            labell = 'Im($\epsilon_1$) = 0'  
        else:
            labell = 'Im($\epsilon_1$) = %.7f'%(im_epsi1)
    
        plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    
    n = 10
    mini,maxi = np.min(list_Qabs_tot1),np.max(list_Qabs_tot1)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    
    if save_graphs==1:
        os.chdir(path_g)
        plt.tight_layout(1)
        plt.savefig(name + '_fino_modo%i_kz%.4f.png' %(modo,kz), format='png') 
        if paper == 0:
            np.savetxt('info_' + name + '_modo%i_kz%.4f_zoom1D.txt' %(modo,kz), [inf_tot],fmt='%s')
    
    graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)   
    for j in range(len(list_Qabs_tot2)):
        list_Qabs2 = np.abs(list_Qabs_tot2[j])
        im_epsi1 = list_im_epsi1_grueso[j]
            
        if im_epsi1 == crit:
            labell = 'Im($\epsilon_1$) = Im($\epsilon_1$)$_c$'
        elif im_epsi1 == -crit:
            labell = 'Im($\epsilon_1$) = -Im($\epsilon_1$)$_c$'    
        elif im_epsi1 == 0:
            labell = 'Im($\epsilon_1$) = 0'  
        else:
            labell = 'Im($\epsilon_1$) = %.3f'%(im_epsi1)
            
        plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
    
    n = 10
    mini,maxi = np.min(list_Qabs_tot2),np.max(list_Qabs_tot2)
    eje_Lambda2 = np.linspace(mini,maxi,n)
    plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
    plt.yscale('log')
    plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))
    if save_graphs==1:
        os.chdir(path_g)
        plt.tight_layout(1)
        plt.savefig(name + '_grueso_modo%i_kz%.4f.png' %(modo,kz), format='png') 

#%%

from matplotlib.colors import SymLogNorm
import matplotlib.colors as colors

if graph_2D==1:
    print('Graficar '+ name + ' en 2D') 
    
    def Qscat2D(Omegac,Im_epsi1): 
        Epsi1 = re_epsi1 + 1j*Im_epsi1
        Qscatt = Cross_Section(kz,Omegac,Epsi1,nmax,R,hbaramu,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
        return (Qscatt)
       
    N = 500
    
    # if zoom==1:
    #     tol = 1e-4
    # else:
    #     tol1 = 5*1e-5
    #     tol2 = 3*1e-3
    
    # # omegac12,omegac22 = omegac0*(1-tol),omegac0*(1+tol)
    # # list_omegac = np.linspace(omegac12,omegac22,N)
    # # delta = (omegac22-omegac12)*0.5
    # # list_im_epsi1 = np.linspace(crit-delta,crit+delta,N)
    # crit1, crit2 = crit*(1+tol2),crit*(1- tol2)

    # list_im_epsi1 = np.linspace(crit - tol2,crit + tol2,N)
    # list_omegac = np.linspace(omegac0 - tol2,omegac0 + tol2,N)    


    if zoom==1:
        tol = 1e-3
    else:
        tol = 1e-2
        # tol = 5*1e-1 #para ver el otro polo
    
    omegac12,omegac22 = omegac0*(1-tol),omegac0*(1+tol)
    list_omegac = np.linspace(omegac12,omegac22,N)
    delta = (omegac22-omegac12)*0.5

    # list_im_epsi1 = np.linspace(im_epsi1c - delta,im_epsi1c + delta,N)
    # list_omegac = np.linspace(omegac0 - 3*tol, omegac0 + 3*tol,N)
    list_im_epsi1 = np.linspace(crit - 5*tol,crit + 5*tol,N)

    x = list_omegac
    y = list_im_epsi1
    X, Y = np.meshgrid(x, y, sparse=True)
    f = np.vectorize(Qscat2D)
    Z = f(X, Y)
        
    limits = [np.min(x) , np.max(x), np.min(y) , np.max(y)]
    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    vmin,vmax = np.min(Z), np.max(Z)
    maxlog=int(np.ceil( np.log10( np.abs(vmax) )))
    minlog=int(np.ceil( np.log10( np.abs(vmin) )))
    
    if vmin < 0 :
          tick_locations = ( [-(10.0**x) for x in np.linspace(minlog,-1,minlog+2)] 
                            + [0] 
                            + [(10.0**x) for x in np.linspace(-1,maxlog,maxlog+minlog+3)] )
    else:
          tick_locations = ( [(10.0**x) for x in np.linspace(minlog,maxlog,maxlog + np.abs(minlog) + 1) ])    
        
    pcm = plt.pcolormesh(X, Y, Z,
                          norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                              vmin=int(vmin), vmax=int(vmax)),cmap='RdBu_r')
    
    plt.plot(x,np.ones(N)*crit,'--',lw = lw,color = 'green')
    plt.plot(np.ones(N)*omegac0,y,'--', lw = lw,color = 'green')
    
#    im = plt.imshow(Z, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(pcm, extend='both')
    cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    # if paper == 0:
    #     cbar.set_label(name,fontsize=tamlegend)
    #plt.legend(loc='best',markerscale=2,fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_g)
        if paper == 1:
            plt.tight_layout()
        plt.savefig(name + '2D' + inf_fig, format='png') 
        if paper == 1:
            np.savetxt('info_' + name + '_modo%i_kz%.4f_zoom2D.txt' %(modo,kz), [inf_tot],fmt='%s')

#%%
