# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values.py: graficar optical gain y omega/c obtenidos
                            en formato paper (lindos) y para todos los modos
"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.interpolate import interp1d # para find_degenerations == 1
from scipy import optimize

#%%

save_graphs = 1
find_degenerations = 0 #encontrar degeneraciones en frecuencia: modos con el mismo omega/c

tamfig = (11,7)
tamlegend = 18
tamletra = 20
tamtitle = 20
tamnum = 16

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')
path_graphene = path_basic.replace('/sin_kz_inv/non-dispersive/real_freq','') 


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

re_epsi1 = 4.9
R = 0.5              #micrones
 
path_load = path_basic + '/re_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_nondispR_%.2f.txt' %(R), [inf_tot], fmt='%s')

list_color = ['darkred','yellowgreen','steelblue','coral']
labelx = '$\mu_c$ [eV]'

#%%

# if R != 0.5:
#     raise TypeError('Wrong value for radium')

#%%


# colors = sns.set_palette("RdBu_r")

# current_palette = sns.color_palette("RdBu_r", 4)
# list_color = current_palette

# sns.palplot(current_palette)
sns.set()

# hspace = 0.1
# wspace = 0.08

# fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
# plt.subplots_adjust(hspace =hspace,wspace = wspace)

# labely = 'Im($\epsilon_1$)'

# k = 0
# for modo in [1,2,3,4]:
#     os.chdir(path_load)
#     name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)
#     tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
#     tabla = np.transpose(tabla)
#     [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla

#     if modo == 1:
#         axs[0].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
#     else:
#         axs[1].plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
        
#     k = k + 1


# ticks1 = [-0.033,-0.035,-0.037,-0.039,-0.041]
# ticks2 = [-0.011,-0.017,-0.023,-0.029,-0.035]

# ticks = [ticks1,ticks2]

# for i in [0,1]:
#     axs[i].minorticks_on()
#     axs[i].tick_params(labelsize = tamnum)
#     # axs[i].set_yticks(ticks[i])
#     # axs[i].set_yticklabels(ticks[i])

# axs[0].set_ylabel(labely,fontsize = tamletra,labelpad=0) 
# axs[1].set_ylabel(labely,fontsize = tamletra,labelpad=6) 
# axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad=2)

# legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
# patchList2 = []
# for key in legend_mode2:
#         data_key = mpatches.Patch(color=legend_mode2[key], label=key)
#         patchList2.append(data_key)
        
# fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 

# s1,s2 = 'a)', 'b)'
# xs = 0.133
# fig.text(xs,0.83,s1,fontsize = tamletra) 
# fig.text(xs,0.422,s2,fontsize = tamletra) 


# if save_graphs==1:
#     os.chdir(path_load)
#     plt.savefig('Im_epsi1_vs_muR_%.2f.png' %(R),format = 'png')
    
#%% 

labely = 'Im($\epsilon_1$)'

plt.figure(figsize=tamfig)

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
    plt.plot(list_mu_opt,epsi1_imag_opt,'-',lw = 5,color = list_color[k])
        
    k = k + 1

plt.ylabel(labely,fontsize=tamletra,labelpad = 0)
plt.xlabel(labelx,fontsize=tamletra,labelpad=2)
plt.tick_params(labelsize = tamnum)

legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)
        
plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 

if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_muR_%.2f.png' %(R),format = 'png')
    
#%%

labely = '$\omega/c$'

plt.figure(figsize=tamfig)

k = 0
for modo in [1,2,3,4]:
    os.chdir(path_load)
    name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)
    tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
    tabla = np.transpose(tabla)
    [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
    plt.plot(list_mu_opt,omegac_opt,'-',lw = 5,color = list_color[k])
        
    k = k + 1

plt.ylabel(labely,fontsize=tamletra)
plt.xlabel(labelx,fontsize=tamletra,labelpad=2)
plt.tick_params(labelsize = tamnum)

legend_mode2 = {'mode 1' : list_color[0], 'mode 2' : list_color[1],'mode 3' : list_color[2], 'mode 4' : list_color[3] }
patchList2 = []
for key in legend_mode2:
        data_key = mpatches.Patch(color=legend_mode2[key], label=key)
        patchList2.append(data_key)
        
plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 


if save_graphs==1:
    os.chdir(path_load)
    plt.savefig('omegac_vs_muR_%.2f.png' %(R),format = 'png')
    
#%%
 
if find_degenerations == 1:
    
    for modo in [1,2,3]:

        os.chdir(path_load)
        name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)      
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla   
        
        name2 = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo + 1)      
        tabla2 = np.loadtxt(name2, delimiter='\t', skiprows=1)
        tabla2 = np.transpose(tabla2)
        [list_mu_opt2,omegac_opt2,epsi1_imag_opt2,eq_det2] = tabla2
        
        if list_mu_opt2[0] != list_mu_opt[0]:
            raise TypeError('Warning: mu del modo', modo, 'distinto del mu del modo', modo + 1)        
    
    def salto_modo1(modo):
        if type(modo) != int:
            raise TypeError('Modo debe ser entero')
        if modo < 4:
            rta = modo + 1
        elif modo == 4:
            rta = 1
        else:
            raise TypeError('Wrong value for mode: mode <= 4')
        return rta

    def salto_modo2(modo):
        if type(modo) != int:
            raise TypeError('Modo debe ser entero')
        if modo < 3:
            rta = modo + 2
        elif modo == 3:
            rta = 1
        elif modo == 4:
            rta = 2
        else:
            raise TypeError('Wrong value for mode: mode <= 4')
        return rta

    def omega_c(x,modo):
        """
        Parameters
        ----------
        x : mu_c obtenido para la condicion de spaser
        modo : mode
        
        Returns
        -------
        omega/c correspondiente al mu_c y al modo
        en condicion de spaser 
        """
        
        os.chdir(path_load)
        name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)      
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    
        
        f = interp1d(list_mu_opt,omegac_opt)
        return float(f(x))   

    def im_epsi1_crit(x,modo):
        """
        Parameters
        ----------
        x : mu_c obtenido para la condicion de spaser
        modo : mode
        
        Returns
        -------
        Im(epsi1) correspondiente al mu_c y al modo
        en condicion de spaser 
        """
        
        os.chdir(path_load)
        name = 'opt_det_sinkz_inv_vs_mu_modo%i.txt' %(modo)      
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla    
        f = interp1d(list_mu_opt,epsi1_imag_opt)
        return float(f(x)) 

    def resta_f(x,modo): 
        """
        Parameters
        ----------
        x : mu_c obtenido para la condicion de spaser
        modo : mode
        
        Returns
        -------
        minimizar la diferencia entre diferentes omega/c
        para hallar degeneraciones
        se usan saltos de a 1 modo

        """
        if modo > 4:
            raise TypeError('Wrong value for mode')
            
        f1 = omega_c(x,modo)
        
        modo2 = salto_modo1(modo)
        f2 = omega_c(x,modo2)
        
        return np.abs(f1 - f2)


    def resta_f2(x,modo): 
        """
        Parameters
        ----------
        x : mu_c obtenido para la condicion de spaser
        modo : mode
        
        Returns
        -------
        minimizar la diferencia entre diferentes omega/c
        para hallar degeneraciones
        se usan saltos de a 2 modos

        """

        if modo > 4:
            raise TypeError('Wrong value for mode')
            
        f1 = omega_c(x,modo)
        
        modo2 = salto_modo2(modo)
        f2 = omega_c(x,modo2)
        
        return np.abs(f1 - f2)
        
    list_modo = [1,2,3,4]
    cond_init = 0.31
    tol = 1e-9
    ite = 1150
    my_dict = {'modes': '','mu_c deg' : '','omega/c deg' : '', 'im_epsi1 del modo 1' : '', 'im_epsi1 del modo 2' : ''}
    
    for nu in list_modo :
        def F(x):
            return resta_f(x,nu)
        try:
            sol = optimize.root(F, cond_init, jac=None, method='hybr', tol=tol, 
                       options={'maxiter':ite})
        except ValueError:
            continue
        nu2 = salto_modo1(nu)
        sol1 = np.round(sol.x[0],9)
        print(nu,nu2, sol1)
        my_dict['modes'] = nu,nu2
        my_dict['mu_c deg'] = sol1
        my_dict['omega/c deg'] = omega_c(sol1,nu)
        
        key1_im_epsi1 = 'im_epsi1 del modo 1'
        my_dict[key1_im_epsi1] = im_epsi1_crit(sol1,nu)
        my_dict['im_epsi1 del modo %i' %(nu)] = my_dict.pop(key1_im_epsi1)
        
        key2_im_epsi1 = 'im_epsi1 del modo 2'
        my_dict[key2_im_epsi1] = im_epsi1_crit(sol1,nu2)
        my_dict['im_epsi1 del modo %i' %(nu2)] = my_dict.pop(key2_im_epsi1)

    print('')
    my_dict2 = {'modes': '','mu_c deg' : '','omega/c deg' : '', 'im_epsi1 del modo 1' : '', 'im_epsi1 del modo 2' : ''}
    for nu in list_modo :
        def F(x):
            return resta_f2(x,nu)
        try:
            sol = optimize.root(F, cond_init, jac=None, method='hybr', tol=tol, 
                       options={'maxiter':ite})
        except ValueError:
            continue
        nu2 = salto_modo2(nu)
        sol2 = np.round(sol.x[0],9)
        print(nu, nu2, sol2)    
        my_dict2['modes'] = nu,nu2
        my_dict2['mu_c deg'] = sol2
        my_dict2['omega/c deg'] = omega_c(sol2,nu)

        key1_im_epsi1 = 'im_epsi1 del modo 1'
        my_dict2[key1_im_epsi1] = im_epsi1_crit(sol2,nu)
        my_dict2['im_epsi1 del modo %i' %(nu)] = my_dict2.pop(key1_im_epsi1)
        
        key2_im_epsi1 = 'im_epsi1 del modo 2'
        my_dict2[key2_im_epsi1] = im_epsi1_crit(sol2,nu2)
        my_dict2['im_epsi1 del modo %i' %(nu2)] = my_dict2.pop(key2_im_epsi1)    
    
    if my_dict['omega/c deg'] != '' or my_dict2['omega/c deg'] != '' :
    
        deg_txt = 'info_degenerations_nondispR_%.2f.txt' %(R)
        np.savetxt(deg_txt, [], fmt='%s')    
        with open(deg_txt, 'w+') as file: 
            file.write('{}\n{}'.format(str(my_dict),str(my_dict2)))
            # file.writelines([inf_tot,str(my_dict),str(my_dict2)])
        file.close()
    
#%%  
