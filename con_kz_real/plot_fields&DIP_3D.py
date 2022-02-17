#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila
mapas de color de los 
campos transversales + longitudinales:

    con un dipolo interior en el cilindro
para con y sin medio activo. El medio activo esta
modelado aca con una constante compleja
"""

import numpy as np
import sys
import os 
import matplotlib.pyplot as plt
from sympy import pi

#%% mirar celda 

save_graphs = True #guardar los graficos mapa de color del campo
non_active_medium = True #plotear campos con im(epsilon1) = 0

paper = True
normalizar = True #normalizar los campos

dipolo_alineado = False # si es True: el dipolo px,py,pz se alinea con el k incidente. 
                    # FALTABA NORMALIZAR EL P DEL DIPOLO. Sino el campo va a aumentar cada vez que aumente kz (porque si es proporcional a p, aumenta |p| y aumenta el |campo|)

list_field = ['Ez', 'Hz', 'Etot', 'Htot'] 
type_field = list_field[2]
if type_field == 'Ez' or type_field == 'Hz':
    modulo = False #if modulo == 1 graficar |Hz| o |Ez| y si vale 0 grafica la parte real

#%% mirar celda 

print('Definir parametros del problema y de los graficos')

med = 8  # poner el valor mas grande de grafico guardado
list_modos = [1,2,3,4] # barrido en modos 
list_ind = [0,50,500,1000,5000,-1]  # barrido en kz
list_ind = [1000]

Re_epsi1 = 3.9
Radio = 0.5              # micrones
hbaramu = 0.3        # eV mu_c

z_var = 0
nmaxx = 5 # sumatoria de campos
A0,B0 = 1,1

### datos del dipolo ####
thetaD = 0
zD = 0
rhoD = Radio/2

if dipolo_alineado == False: # elegir el momento dipolar 
    px_value = 1
    py_value = 1
    pz_value = 1

N = 150 # longitud de meshgrid
cota = 2*Radio # x e y van desde - cota hasta +cota
arrow = Radio*0.3 #longitud de la flecha del dipolo

#%% mirar celda 

if paper == 1: 
    tamfig = (4.5,3.5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 9
    labelpady = -1.5
    labelpadx = -0.5
    pad = 0.5
else:
    tamfig = (5,5)
    tamlegend = 10
    tamletra = 11
    tamtitle = 10
    tamnum = 7
    labelpady = -1.5
    labelpadx = -0.5
    pad = 0.5

#%%

if rhoD > Radio:
    raise TypeError('Invalid value for r_D: r_D must be < R (dipole inside cylinder)')    

if hbaramu!= 0.3:
    raise TypeError('Wrong value for chemical potential of graphene')

if Re_epsi1!= 3.9:
    raise TypeError('Wrong value for Re(epsi1)')
    
if Radio != 0.5 :
    raise TypeError('Wrong value for radius')

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save0 = path_basic + '/' + 'fields&DIP_3D'
path_save1 = path_save0 + '/' + type_field
path_graphene = path_basic.replace('/' + 'con_kz_real','') 

if save_graphs==True:    
    if not os.path.exists(path_save0):
        os.mkdir(path_save0)
    if not os.path.exists(path_save1):
        print('Creating folder to save graphs')
        os.mkdir(path_save1)    

#print('Importar modulos necesarios para este codigo')
try:
    sys.path.insert(1, path_basic)
    if type_field == 'Htot':
        from fields_conkz_DIP import Htot as fields
    elif type_field == 'Etot':
        from fields_conkz_DIP import Etot as fields
    elif type_field == 'Hz':
        from fields_conkz_DIP import Hz
        if modulo == True:
            def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
                rta = Hz(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
                return np.abs(rta)
        else:
            def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
                rta = Hz(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
                return rta.real
    elif type_field == 'Ez':
        from fields_conkz_DIP import Ez
        if modulo == True:
            def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
                rta = Ez(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
                return np.abs(rta)
        else:
            def fields(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo):
                rta = Ez(kz,omegac,epsi1,nmax,R,mu_c,rho,phi,z,p1,p2,pz,rho_D,theta_D,z_D,Ao,Bo)
                return rta.real                 
except ModuleNotFoundError:
    print('fields_conkz_DIP.py no se encuentra en ' + path_basic)
    print('reemplazar path_basic por el path donde se encuentra fields_conkz_DIP.py')

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en ' + path_graphene)
    print('reemplazar path_graphene por el path donde se encuentra constantes.py')

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()

#%%

info4 = r'$\rho_D$ = %.3f$\mu$m, $\theta_D$ = %s, $z_D$ = %i$\mu$m' %(rhoD,thetaD,zD)
thetaD = float(thetaD)
coordX_D, coordY_D = rhoD*np.cos(thetaD),rhoD*np.sin(thetaD)

def circle(radio):
    listx = []
    listy = []
    for theta in np.linspace(0,2*np.pi):
        x = radio*np.cos(theta)
        y = radio*np.sin(theta)
        listx.append(x)
        listy.append(y)
    return listx,listy

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,ptot,px,py):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    
    # plt.plot([coordX_D], [coordY_D] ,'x',color = 'blue')
    listx2,listy2 = circle(Radio)
    plt.plot(listx2,listy2,'-',color = 'blue')
    # p_theta = np.arctan(py/px)
    # px_theta = np.cos(p_theta)
    # py_theta = np.sin(p_theta)
    pvector = [px,py]
    pvector_normalize = pvector/np.linalg.norm(pvector)
    [px_val,py_val] = pvector_normalize
    
    if px + py !=0:
        plt.arrow(coordX_D, coordY_D , px_val*0.3, py_val*0.3, head_width = 0.05,
          width = 0.02,color = 'blue')
    else:
        plt.plot(coordX_D, coordY_D,'x',color = 'blue')
    
    if paper == False:
        plt.title(title,fontsize=int(tamtitle*0.9))
    if paper == True:
        plt.tight_layout()
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

    
if dipolo_alineado == True:
    # dipolo alineado al k incidente        
    def kt(kz_var,k0): #k transversal (unidad 1/micrometros)
        
        xz = kz_var/k0
        Rbarra = Radio*k0 #adimensional
        
        xt2rta = -xz**2 + mu2*epsi2
        
        inside = xt2rta + 0j
        xtmedio = (inside)**(1/2)   
        argument = xtmedio*Rbarra
    
        if argument.real>=0:
    	    xtmedio = xtmedio
        else:
    	    xtmedio = -xtmedio
    
        return xtmedio*k0
    
#%%
    
x = np.linspace(-cota,cota,N)
y = np.linspace(-cota,cota,N)
X, Y = np.meshgrid(x, y)
limits = [min(x) , max(x), min(y) , max(y)]

labelx,labely = 'x [$\mu$m]', 'y [$\mu$m]'
if type_field == 'Htot':
    labelz = '|Htot|$^2$'
    campo = 'Htot'
elif type_field == 'Etot':
    labelz = '|Etot|$^2$'
    campo = 'Etot'
else: 
    if modulo == True:
        labelz = '|' + type_field + '|'
        campo = labelz
    else:
        labelz = 'Re(' + type_field + ')'
        campo  =  'Re' + type_field 
#    name = type_field

#%%

path_load = path_basic  + '/' + 'real_freq' + '/' + 're_epsi1_%.2f_vs_kz/R_%.2f/mu_%.1f' %(Re_epsi1,Radio,hbaramu) 
    
if save_graphs == True:    
    if not os.path.exists(path_save0):
        print('Creating folder to save graphs')
        os.mkdir(path_save0)
    
def graficar_campos(modo,ind,med):
    
    os.chdir(path_load)
    print('Importar los valores de SPASER')
            
    name = 'opt_det_conkz_vs_kz_modo%i.txt' %(modo)        
    try:
        data = np.loadtxt(name,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(name) if l.startswith('#')):
            print('valores de ', name, ':', line)
    except OSError or IOError:
        print('El archivo ' + name + ' no se encuentra en ' + path_load)
        print('Cambiar path_load')
            
    data = np.transpose(data)
    [list_kz_opt,omegac_opt,epsi1_imag_opt,eq_det] = data

    ind = int(ind)
    kz = list_kz_opt[ind] #micrones
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = Re_epsi1 + 1j*im_epsi1
    
    if dipolo_alineado == True:
        a = kt(kz,omegac).real
        b = kt(kz,omegac).imag
        k_tot = kt(kz,omegac)**2 + kz**2
        p_cuadrado = (a**2 + b**2)/(1-(kz/np.abs(k_tot))**2)
        p_cuadrado = p_cuadrado.real
        
        px = a
        py = b
        pz = kz/np.abs(k_tot)
        pvector = [px,py,pz]
        pvector_normalize = pvector/np.linalg.norm(pvector)
        [px,py,pz] = pvector_normalize
    
    else:
        px = px_value
        py = py_value
        pz = pz_value
    
    ptot = px + py + pz 
    p1 = px + 1j*py # p+
    p2 = px - 1j*py # p-

    labelpx = 'px' + labelp(px)
    labelpy = 'py' + labelp(py)
    
    # labelp1 = 'p+' + labelp(p1)
    # labelp2 = 'p-' + labelp(p2)
    
    info1 = r'Ao = %i, Bo = %i, kz = %.4f 1/$\mu$m, z = %i$\mu$m' %(A0,B0,kz,z_var)
    info2 = r'R = %.2f$\mu$m, nmax = %i, $\mu_c$ = %.1f eV, $\nu$ = %i' %(Radio,nmaxx,hbaramu,modo)
    info3 = r'$\epsilon_1$ = %.1f - i%.4e, $\omega/c$ = %.4e 1/$\mu$m' %(Re_epsi1,-im_epsi1,omegac)
    infoDIP = labelpx + ', ' + labelpy + r', pz = %.2f' %(pz)
    title =  info1 + '\n' + info2 + '\n' + info3  + ', ' + info4 + '\n' + infoDIP
    info_tot = info1 + ', ' + info2 + ', ' + info3 + ', ' + info4 + ', ' + infoDIP + ', ' + name_this_py
    
    print('modo = ', modo)
    print(labelpx)
    print(labelpy)
    print('pz = ', pz)
    print('kz = ', kz)
                
    print('Graficar el campo ' + labelz + ' para el medio 1 y 2')
            
    if non_active_medium == True: #epsi1 = re_epsi1 ---> no hay medio activo
        info3_loss =  r'$\epsilon_1$ = %.1f, $\omega/c$ = %.4e 1/$\mu$m' %(Re_epsi1,omegac)
        title_loss =  info1 + '\n' + info2 + '\n' + info3_loss + ', ' + info4 + '\n' + infoDIP

        def Etot_2variable(x,y):   
            # phi = np.arctan2(y,x)
            phi = np.arctan2(y,x) # cuadrante [-pi,pi]
            if phi < 0:
                phi = phi + 2*np.pi # cuadrante [0,2pi]
            rho = (x**2+y**2)**(1/2)
            campo = fields(kz,omegac,Re_epsi1,nmaxx,Radio,hbaramu,rho,phi,z_var,p1,p2,pz,rhoD,thetaD,zD,A0,B0)
            return campo
        
        f2 = np.vectorize(Etot_2variable)
        Z2 = f2(X, Y)
        if normalizar == True:
            maxi = np.max(Z2)
            Z2 = Z2/maxi
        
        graph(title_loss,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,ptot,px,py)
        im = plt.imshow(Z2, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
        cbar = plt.colorbar(im)
        cbar.ax.tick_params(labelsize = tamnum)
        
        if paper == False:
            cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
        
        if save_graphs == True:
            os.chdir(path_save1)
            plt.savefig(campo + '_loss_med%i_modo%i.png' %(med,modo), format='png')
 
    np.savetxt('info_'+ campo + '_med%i.txt' %(med), [info_tot], fmt='%s') 
    def Etot_2variable(x,y):   
        phi = np.arctan2(y,x) # cuadrante [-pi,pi]
        if phi < 0:
            phi = phi + 2*np.pi # cuadrante [0,2pi]
        rho = (x**2+y**2)**(1/2)
        campo = fields(kz,omegac,epsi1,nmaxx,Radio,hbaramu,rho,phi,z_var,p1,p2,pz,rhoD,thetaD,zD,A0,B0)
        return campo
        
    f1 = np.vectorize(Etot_2variable)
    Z1 = f1(X, Y)
    if normalizar == True: #normarlizar TODOS los campos : con y sin medio activo
        maxi2 = np.max(Z1) 
        Z1 = Z1/maxi2
    
    graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad,ptot,px,py)
    im = plt.imshow(Z1, extent = limits, cmap=plt.cm.hot, interpolation='bilinear')
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize = tamnum)
    
    if paper == False :
        cbar.set_label(labelz,fontsize=tamlegend,labelpad = 1)
    
    if save_graphs == True:
        os.chdir(path_save1)
        plt.savefig(campo + '_med%i_modo%i.png' %(med,modo), format='png')
    
    del Z1,Z2
    print('')
    
#%%

for ind in list_ind: 
    for modo in list_modos:
        try:
            med = med + 1
        except:
            med = 1
        graficar_campos(modo,ind,med)

#%%
