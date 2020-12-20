#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 22:07:37 2020

@author: leila

Vamos a usar la formula que yo obtuve
del Qscat (que aparece en el overleaf
           en la seccion 4.8)

"""
import numpy as np
import sys
import os 
import matplotlib.pyplot as plt

save_graphs = 1 #guardar los graficos 2D del campo
graph_2D = 0    #graficos 2D

colores = ['coral','yellowgreen','midnightblue','green','darkred','aquamarine','hotpink','steelblue','purple']

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_save = path_basic + '/' + 'seccion_eficaz_bestpol'

if save_graphs==1:
    if not os.path.exists(path_save):
        print('Creating folder to save graphs')
        os.mkdir(path_save)

#print('Importar modulos necesarios para este codigo')

try:
    sys.path.insert(1, path_basic)
    from Qscat_conkz import Qscat
except ModuleNotFoundError:
    print('Qscat_conkz.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_basic = input('path de la carpeta donde se encuentra Qscat_conkz.py')
    sys.path.insert(1, path_basic)
    from Qscat_conkz import Qscat

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
R = 0.5             #micrones
hbaramu = 0.3        #eV mu_c

list_modos = [1,2,3,4]
list_ind = [1,5,15,25,40,50,100,-20,-1] #lista de indices
nmax = 10       #sumatoria desde -nmax hasta +nmax (se suman 2*nmax + 1 modos)

Ao,Bo = 0,1        #Ao ---> Hz, Bo ---> Ez

#%%

if re_epsi1 != 3.9:
    raise TypeError('Wrong value for Re(epsi1): Re(epsi1) = 3.9')
    
if R != 0.5:
    raise TypeError('Wrong value for R: R = 0.5 micrones')
    
if hbaramu != 0.3:
    raise TypeError('Wrong value for chemical potential: mu = 0.3')

#%%

if Ao*Bo != 0:
    path_save0 = path_save + '/' + '2pol'
    labelAoBo = 'Ao = 1/$\sqrt{2}$, Bo = 1/$\sqrt{2}$'
elif Bo == 0:
    path_save0 = path_save + '/' + 'polAo'
    labelAoBo = 'Ao = 1, Bo = 0'
elif Ao == 0:
    path_save0 = path_save + '/' + 'polBo'
    labelAoBo = 'Ao = 0, Bo = 1'

norm = np.linalg.norm([Ao,Bo])
Ao,Bo = Ao/norm, Bo/norm    

if save_graphs==1:    
    if not os.path.exists(path_save0):
        print('Creating folder to save graphs')
        os.mkdir(path_save0)

path_loadAoBO = path_basic  + '/' + 'best_pol2'

#%%

print('')
print('Importar los valores de SPASER y de su mejor polarizacion correspondiente (estan en un mismo .txt)')

def names(modo,ind):
    
    def labeltot(Aocte,Bocte):
        def label(var):
            if np.sign(var.imag) == 1 or var.imag ==0:
                label = '%.3f + i%.3f' %(var.real,var.imag)
            else:
                label = '%.3f - i%.3f' %(var.real,-var.imag)
            return label
        labelAoBo = 'Ao = ' + label(Aocte) + ', ' + 'Bo = ' + label(Bocte)
        return labelAoBo
    
    os.chdir(path_loadAoBO)

    nameAo1 = 'C1best_vs_kz_modo%i.txt' %(modo)
    nameBo1 = 'A1best_vs_kz_modo%i.txt' %(modo)
    nameAo2 = 'D2best_vs_kz_modo%i.txt' %(modo)
    nameBo2 = 'B2best_vs_kz_modo%i.txt' %(modo)        
        
    try:
        dataAo1 = np.loadtxt(nameAo1,delimiter = '\t', skiprows=1)
        for line in (l.strip() for l in open(nameAo1) if l.startswith('#')):
            print('values de ', nameAo1, ':', line)
    except OSError or IOError:
        print('El archivo ' + nameAo1 + ' no se encuentra en ' + path_loadAoBO)
    print('')
    
    dataAo1 = np.transpose(dataAo1)
    [list_re_Ao1,list_im_Ao1,list_kz_optAo1,omegac_opt,epsi1_imag_opt] = dataAo1
    
    dataBo1 = np.loadtxt(nameBo1,delimiter = '\t', skiprows=1)
    dataBo1 = np.transpose(dataBo1)
    [list_re_Bo1,list_im_Bo1,list_kz_optBo1,omegac_opt,epsi1_imag_opt] = dataBo1    

    dataAo2 = np.loadtxt(nameAo2,delimiter = '\t', skiprows=1)
    dataAo2 = np.transpose(dataAo2)
    [list_re_Ao2,list_im_Ao2,list_kz_optAo2,omegac_opt,epsi1_imag_opt] = dataAo2   

    dataBo2 = np.loadtxt(nameBo2,delimiter = '\t', skiprows=1)
    dataBo2 = np.transpose(dataBo2)
    [list_re_Bo2,list_im_Bo2,list_kz_optBo2,omegac_opt,epsi1_imag_opt] = dataBo2   
    
    if list_kz_optAo1[ind] != list_kz_optBo1[ind]:
        raise TypeError('Differents values of kz for medium 1')

    if list_kz_optAo2[ind] != list_kz_optBo2[ind]:
        raise TypeError('Differents values of kz for medium 2')
    
    kz = list_kz_optAo1[ind] #micrones
    im_epsi1 = epsi1_imag_opt[ind]
    omegac = omegac_opt[ind] 
    epsi1 = re_epsi1 + 1j*im_epsi1
    Ao1 = list_re_Ao1[ind] + 1j*list_im_Ao1[ind]
    Bo1 = list_re_Bo1[ind] + 1j*list_im_Bo1[ind]
    
    Ao2 = list_re_Ao2[ind] + 1j*list_im_Ao2[ind]
    Bo2 = list_re_Bo2[ind] + 1j*list_im_Bo2[ind]
    print('kz = ', kz)
    
    # Ao = coef_D2        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
    # Bo = coef_B2  
    
    # Ao = coef_C1        #Ao es el de Hz ---> coef_D2 de la formula (ver fieldsZ_conkz.py)
    # Bo = coef_A1  

    if kz < 0.13:
        path_save = path_save0 + '/' + 'kz_chico'
    else:
        path_save = path_save0 + '/' + 'kz_grande'
    
    if save_graphs==1:    
        if not os.path.exists(path_save):
            print('Creating folder to save graphs')
            os.mkdir(path_save)    

    norm1 = np.linalg.norm([Ao1,Bo1])
    Ao1,Bo1 = Ao1/norm1, Bo1/norm1
    
    norm2 = np.linalg.norm([Ao2,Bo2])
    Ao2,Bo2 = Ao2/norm2, Bo2/norm2    

    labelAo1Bo1 = labeltot(Ao1,Bo1)
    labelAo2Bo2 = labeltot(Ao2,Bo2)

    im_epsi1 = epsi1.imag

    info1 = 'modo = %i, kz = %.4f 1/$\mu$m, R = %.1f $\mu$m' %(modo,kz,R)
    info2 = 'nmax = %i, $\mu_c$ = %.1f eV, $\epsilon_1$ = %.1f - i%.2e' %(nmax,hbaramu,re_epsi1,-im_epsi1)
    info3 = ' y $\omega/c$ = %.2e 1/$\mu$m' %(omegac)
    title = info1 + '\n' + info2 + info3 + '\n' + name_this_py
    
    return kz,epsi1,omegac,Ao1,Bo1,Ao2,Bo2,labelAo1Bo1,labelAo2Bo2,path_save,title

#%%

print('')
print('Calcular Qscat para diferentes Im(epsi1)')

labelx = r'$\omega/c$ [1/$\mu$m]'
labely = 'Qscat/c'

tol = 1e-3
N = int(1e3) 
for ind in list_ind:
    for modo in list_modos:
        
        kz,epsi1,omegac0,Ao1,Bo1,Ao2,Bo2,labelAo1Bo1,labelAo2Bo2,path_save,title = names(modo,ind)
        
        # crit = epsi1.imag
        # list_im_epsi1_fino = [0,crit+3*tol,crit+2*tol,crit+tol,crit,crit-tol,crit-2*tol]
        
        if modo == 1 or modo == 2:            
            omegac1,omegac2 = omegac0*0.97,omegac0*1.03
        elif modo == 3: 
            omegac1,omegac2 = omegac0*0.9999,omegac0*1.0001
        else:
            omegac1,omegac2 = omegac0*0.99999,omegac0*1.00001
        
        list_omegac = np.linspace(omegac1,omegac2,N)
        
        list_Qscat0 = []
        list_Qscat1 = []
        list_Qscat2 = []
        for omeggac in list_omegac:
            # epsi1 = re_epsi1 + 1j*im_epsi1
            
            Qscatt0 = Qscat(kz,omeggac,epsi1,nmax,R,hbaramu,Ao,Bo)
            list_Qscat0.append(Qscatt0)  
            
            Qscatt1 = Qscat(kz,omeggac,epsi1,nmax,R,hbaramu,Ao1,Bo1)
            list_Qscat1.append(Qscatt1)  
    
        
            Qscatt2 = Qscat(kz,omeggac,epsi1,nmax,R,hbaramu,Ao2,Bo2)
            list_Qscat2.append(Qscatt2)  
            
        plt.figure(figsize=tamfig)
        plt.plot(list_omegac,list_Qscat0,'.r',ms=10,label = labelAoBo)
        plt.plot(list_omegac,list_Qscat1,'.b',ms=11,label = labelAo1Bo1)
        plt.plot(list_omegac,list_Qscat2,'.m',ms=10,label = labelAo2Bo2)
        plt.title(title,fontsize=tamtitle)
        plt.ylabel(labely,fontsize=tamletra)
        plt.xlabel(labelx,fontsize=tamletra)
        plt.tick_params(labelsize = tamnum)
        plt.yscale('log')
        plt.legend(loc='lower left',markerscale=2,fontsize = int(tamlegend*0.8))
        plt.grid(1)   
        
        if save_graphs==1:
            os.chdir(path_save)
            plt.savefig('Qscat_index%i_modo%i.png' %(ind,modo))
       
#%%