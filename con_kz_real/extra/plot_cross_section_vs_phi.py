#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

graficar Qscat Qabs o Qext (elegir en la variable cross_section)

ESTABAN MAL ESAS FORMULAS DEL CUADERNO (movidas a extra)
VER TESIS
(revisar cuentas)

variar el angulo phi de incidencia (ver overleaf "scattering con varios kz")
https://www.overleaf.com/read/tjgwynvjksvk

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

#%%

save_graphs = 1 #guardar los graficos 
zoom = 0        #graficos 2D con o sin zoom
paper = 0    #formato paper

list_cross_section = ['Qscat', 'Qabs', 'Qext'] 
cross_section = list_cross_section[0]

#%% 

#print('importar funciones para graficarlas')

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')

if save_graphs==1:
    path_g = path_basic + '/' + 'seccion_eficaz'
    if not os.path.exists(path_g):
        print('Creating folder to save graphs')
        os.mkdir(path_g)

#print('Importar modulos necesarios para este codigo')
err = 'cross_sections_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo'
err2 = 'path de la carpeta donde se encuentra cross_sections_conkz.py'
if cross_section == 'Qscat':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qscat2
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qscat2
        
    def Cross_Section(kz_var,rho,phi,z,omegac,epsi1,nmax,R,hbaramu,Ao,Bo):
        return Qscat2(kz_var,rho,phi,z,omegac,epsi1,nmax,R,hbaramu,Ao,Bo)

elif cross_section == 'Qabs':
    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qabs
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qabs
        
    def Cross_Section(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo):
        return Qabs(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo)

else: 

    try:
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qext
    except ModuleNotFoundError:
        print(err)
        path_basic = input(err2)
        sys.path.insert(1, path_basic)
        from cross_sections_conkz import Qext
        
    def Cross_Section(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo):
        return Qext(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo)

#%%   

print('Definir parametros del problema')

#valores de minimizo perdidas (ver header)
re_epsi1 = 3.9
R = 0.5 #micrones
hbaramu = 0.3        #eV mu_c
modo = 1
Ao,Bo = 1,1
rho = R/2
z = 0

nmax = 5
ind = 99
ind = 150

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

print('Importar datos de SPASER')

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/mu_%.1f' %(re_epsi1,hbaramu) 
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

info1 = 'kz = %.4f $\mu m^{-1}$, R = %.1f$\mu$m, nmax = %i, $\mu_c$ = %.1f eV' %(kz,R,nmax,hbaramu) +infAoBo 
info2 = '$\epsilon_1$ = %.1f - i%.5e y $\omega/c$ = %.5e 1/$\mu$m del modo = %i' %(re_epsi1,-crit,omegac0,modo)
info3 = r'$\rho$ = R/2, z = 0$\mu$m'

inf_tot = info1 + ', ' + info2  + ', ' + name_this_py
title = info1 +'\n' + info2 + '\n'  + info3 + ', ' + name_this_py
labelomegac = '$\omega/c$ = %.2f$\mu m^{-1}$' %(omegac0)
labelx = '$\omega/c$ [$\mu m^{-1}$]'
labely2 = 'Im($\epsilon_1$)'
labely1 = cross_section + '$_{ad}$'
name = cross_section
inf_fig = '_modo%i_kz%.4f.png' %(modo,kz)
if zoom == 1: 
    inf_fig = '_modo%i_kz%.4f_zoom.png' %(modo,kz)

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    if paper == 0:
        plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%
    
if re_epsi1 != 3.9:
    raise TypeError('Wrong value for Re(epsi1): Re(epsi1) = 3.9')
    
if R != 0.5:
    raise TypeError('Wrong value for R: R = 0.5 micrones')
    
if hbaramu != 0.3:
    raise TypeError('Wrong value for chemical potential: mu = 0.3')

#%%


print('')
print('Calcular '+ name + ' para diferentes angulos phi')

list_phi = [0,np.pi/4,np.pi/2,3*np.pi/4,np.pi,1.5*np.pi]
list_label = ['$\phi = 0$', '$\phi = \pi/4$', '$\phi = \pi/2$', '$\phi = 3\pi/4$','$\phi = \pi$' , '$\phi = 1.5\pi$']


N = int(1e3)               
omegac1,omegac2 = omegac0*0.98,omegac0*1.02
list_omegac = np.linspace(omegac1,omegac2,N)

list_Qabs_tot1 = []
for phi in list_phi:
    print(phi)
    list_Qabs1 = []
    for omeggac in list_omegac:
        epsi1 = re_epsi1 + 1j*crit
        Qabss = Cross_Section(kz,rho,phi,z,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo)
        list_Qabs1.append(Qabss)  
    list_Qabs_tot1.append(list_Qabs1)  

del list_Qabs1
            

colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']    
print('Graficar ' + name  + ' para diferentes Im(epsi1)')

graph(title,labelx,labely1,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)    
for j in range(len(list_Qabs_tot1)):
    list_Qabs2 = np.abs(list_Qabs_tot1[j])
    labell = list_label[j]
    plt.plot(list_omegac,list_Qabs2,'o',color = colores[j],ms = 4,alpha = 0.8,label = labell)
n = 10
mini,maxi = np.min(list_Qabs_tot1),np.max(list_Qabs_tot1)
eje_Lambda2 = np.linspace(mini,maxi,n)
plt.plot(omegac0*np.ones(n),eje_Lambda2,'-k',lw = 1,label = labelomegac)
plt.yscale('log')
plt.legend(loc='best',markerscale=2,fontsize=int(tamlegend*0.7))

if save_graphs==1:
    os.chdir(path_g)
    # plt.tight_layout(1)
    plt.savefig(name + '_modo%i_kz%.4f_vs_phi.png' %(modo,kz), format='png') 
    if paper == 0:
        np.savetxt('info_' + name + '_modo%i_kz%.4f_zoom1D.txt' %(modo,kz), [inf_tot],fmt='%s')

#%%
