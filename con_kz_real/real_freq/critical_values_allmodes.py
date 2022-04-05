# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Sobre critical_values.py: graficar optical gain y omega/c obtenidos formato paper

"""

import numpy as np
import os 
import sys
import matplotlib.pyplot as plt
import seaborn as sns

save_graphs = 1
paper = 1
opcion1 = 1

#%%

if paper == 1: 

    if opcion1 == 1:
    
        tamfig = (3, 2)
        tamletra = 9
        tamnum = 8.5
        tamlegend = 8.5
        ticks_x = [0.4,0.6,0.8]
        ticks_y1_kz1 = [-0.036,-0.034,-0.032]
        ticks_y2_kz1 = [-0.012,-0.020,-0.028]
    
        ticks_y1_kz2 = [-0.761, -0.770, -0.779]
        ticks_y2_kz2 = [-0.012,-0.020,-0.028]
            
        columnspace = 0.5
        markerscale = 0.7
        loc1 = [0.275,0.9]  
        length_marker = 1
        labelpady = 1.1
        labelpadx = 0.7
        
    else:    
        tamfig = (3, 2)
        tamletra = 7
        tamnum = 6
        tamlegend = 6
        ticks_x = [0.3,0.4,0.5,0.6,0.7,0.8,0.9]
        ticks_y1_kz1 = [-0.036,-0.035,-0.034,-0.033,-0.032]
        ticks_y2_kz1 = [-0.010,-0.015,-0.020,-0.025,-0.030]
        
        ticks_y1_kz2 = [-0.76 , -0.765, -0.77 , -0.775, -0.78 ]
        ticks_y2_kz2 = [-0.010,-0.015,-0.020,-0.025,-0.030]

        columnspace = 3.5
        markerscale = 0.7
        loc1 = [0.052,1]
        length_marker = 3
        labelpady = 0.2
        labelpadx = 0.2

    lw = 1
    pad = -4        
    loc2 = [0.231,0.9]
    columnspace2 = 0.9

    dpi = 500
    
else:
    tamfig = (11,7)
    tamlegend = 18
    tamletra = 20
    tamtitle = 20
    tamnum = 16
    
    labelpady = 0
    labelpadx = 4
    lw = 3.5
    pad = 0
    loc1 = [0.15,0.88]
    loc2 = [0.025,1]

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_graphene = path_basic.replace('/con_kz_real/real_freq','') 


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

re_epsi1 = 3.9
R = 0.5              #micrones
list_kz = [0,0.0001,0.05,0.1,0.135,0.14,0.5] 
list_kz = [0.5] 
list_modos = [1,2,3,4]

path_load = path_basic + '/re_epsi1_%.2f_vs_mu/R_%.2f' %(re_epsi1,R)

info1 = 'R = %.1f $\mu$m, Re($\epsilon_1$) = %.2f' %(R,re_epsi1)
info2 = '$gamma_c$ = %.4f eV, $\mu_1$ = %i, $\mu_2$ = %i, $\epsilon_2$ = %i' %(hbargama,mu1,mu2,epsi2)
inf_tot = info1 + ',' + info2 + '. Ver ' + name_this_py

os.chdir(path_load)
np.savetxt('info_critical_values_conkz.txt', [inf_tot], fmt='%s')

colors = ['darkred','yellowgreen','coral','steelblue']
symbols = ['-','-','--','-.']
sns.set()

#%%

if R != 0.5:
    raise TypeError('Wrong value for radium')

#%%

hspace = -0.1
wspace = 0.03

for kz in list_kz: 
    
    if kz != 0.5:
        fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
        plt.subplots_adjust(hspace =hspace,wspace = wspace)
        
        labelx = '$\mu_c$ [eV]'
        labely = '[Im($\epsilon_1$)]$_c$'
        
        k = 0
        for modo in list_modos:
            os.chdir(path_load)
            name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
            tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
            tabla = np.transpose(tabla)
            [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        
            if modo == 1:
                line1, = axs[0].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
                line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
                if kz == 0.05:
                    axs[0].set_yticks(ticks_y1_kz1)
                elif kz == 0.14:
                    axs[0].set_yticks(ticks_y1_kz2)            
            else:
                axs[1].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
                if kz == 0.05:
                    axs[1].set_yticks(ticks_y2_kz1)
                elif kz == 0.14:
                    axs[1].set_yticks(ticks_y2_kz2)   
            k = k + 1
        
        
        for i in [0,1]:
            axs[i].minorticks_on()
            axs[i].tick_params(labelsize = tamnum,pad = pad)
    
        
        axs[0].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
        axs[1].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
        axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad =labelpadx) 
        
        #fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        fig.legend(loc = loc1,ncol = 4,markerscale=1,fontsize=tamlegend, 
                          columnspacing = columnspace,frameon=0,handletextpad=0.2, handlelength=length_marker)    
        # s1,s2 = 'a)', 'b)'
        # xs = 0.133
        # fig.text(xs,0.83,s1,fontsize = tamletra) 
        # fig.text(xs,0.422,s2,fontsize = tamletra) 
        
        
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_load)
            plt.savefig('Im_epsi1_vs_mu_kz%.4f.png'%(kz) ,format = 'png',dpi = dpi)
            
        labely = '$\omega/c$ [$\mu$m$^{-1}$]'
        
        fig,axs = plt.subplots(2,1, sharex=True, facecolor='w', figsize = tamfig)
        plt.subplots_adjust(hspace =hspace,wspace = wspace)
    
        k = 0
        for modo in list_modos:
            os.chdir(path_load)
            name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
            tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
            tabla = np.transpose(tabla)
            [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        
            if modo == 1:
                line1, = axs[0].plot(list_mu_opt,omegac_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
                line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
                # if kz == 0.05:
                #     axs[0].set_yticks(ticks_y1_kz1)
                # elif kz == 0.14:
                #     axs[0].set_yticks(ticks_y1_kz2)            
            else:
                axs[1].plot(list_mu_opt,omegac_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
                # if kz == 0.05:
                #     axs[1].set_yticks(ticks_y2_kz1)
                # elif kz == 0.14:
                #     axs[1].set_yticks(ticks_y2_kz2)   
            k = k + 1
        
        
        for i in [0,1]:
            axs[i].minorticks_on()
            axs[i].tick_params(labelsize = tamnum,pad = pad)
    
        
        axs[0].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
        axs[1].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
        axs[1].set_xlabel(labelx,fontsize = tamletra,labelpad =labelpadx) 
        
        #fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        fig.legend(loc = loc2,ncol = 4,markerscale=1,fontsize=tamlegend, 
                          columnspacing = columnspace2,frameon=0,handletextpad=0.2, handlelength=length_marker)    
        # s1,s2 = 'a)', 'b)'
        # xs = 0.133
        # fig.text(xs,0.83,s1,fontsize = tamletra) 
        # fig.text(xs,0.422,s2,fontsize = tamletra) 
        
        
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_load)
            plt.savefig('Omegac_vs_mu_kz%.4f.png'%(kz) ,format = 'png',dpi = dpi)
        
        plt.figure(figsize=tamfig)
        k = 0
        for modo in list_modos:
            os.chdir(path_load)
            name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
            tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
            tabla = np.transpose(tabla)
            [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
            line1, = plt.plot(list_mu_opt,omegac_opt,symbols[k],lw =lw,color = colors[k],label = 'modo %i'%(k+1))
            if k == 0:
                line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
                
            k = k + 1
        
        plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
        plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
        plt.tick_params(labelsize = tamnum,pad = pad)
        plt.xticks(ticks_x)  
                
        #plt.legend(handles=patchList2,loc=[0.02,0.99],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
    #    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
        fig.legend(loc = loc1,ncol = 4,markerscale=1,fontsize=tamlegend, 
                          columnspacing = columnspace,frameon=0,handletextpad=0.2, handlelength=length_marker)     
        if save_graphs==1:
            plt.tight_layout()
            os.chdir(path_load)
            plt.savefig('omegac_vs_mu_kz%.4f.png'%(kz) ,format = 'png',dpi = dpi)
            
#%%
loc4 = [0.012,0.98]
if kz == 0.5:
    
    labelx = '$\mu_c$ [eV]'
    labely = '[Im($\epsilon_1$)]$_c$'
         
    plt.figure(figsize=tamfig)

    k = 0
    for modo in list_modos:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        
        if modo == 1:
            line1, = plt.plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
            line1.set_dashes([2, 2, 10, 2])
        
        else:
            plt.plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
        
        k = k + 1

plt.legend(loc = loc4,ncol = 4,markerscale=1,fontsize=tamlegend, 
                      columnspacing = columnspace,frameon=0,handletextpad=0.2, handlelength=length_marker)    
            
plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
plt.tick_params(labelsize = tamnum,pad = pad)
if save_graphs==1:
    plt.tight_layout()    
    os.chdir(path_load)
    plt.savefig('Im_epsi1_vs_mu_kz%.4f.png'%(kz) ,format = 'png',dpi = dpi)


"""

if kz == 0.5:
    loc1 = [0.275,0.94]
    ticks_y1_kz3 = [-0.76,-0.77,-0.78]
    ticks_y2_kz3 = [-0.4224,-0.4225,-0.4226]
    ticks_y3_kz3 = [-0.889,-0.898,-0.906]
    ticks_y4_kz3 = [-0.275,-0.265,-0.255]
    
    tamfig2 = (3, 3.5)
    fig,axs = plt.subplots(4,1, sharex=True, facecolor='w', figsize = tamfig2)
    
    labelx = '$\mu_c$ [eV]'
    labely = '[Im($\epsilon_1$)]$_c$'
    tol = 2*1e-3
    k = 0
    for modo in list_modos:
        os.chdir(path_load)
        name = 'opt_det_conkz_vs_mu_kz%.4f_modo%i.txt' %(kz,modo)
        tabla = np.loadtxt(name, delimiter='\t', skiprows=1)
        tabla = np.transpose(tabla)
        [list_mu_opt,omegac_opt,epsi1_imag_opt,eq_det] = tabla
        maxi = np.max(epsi1_imag_opt)
        mini = np.min(epsi1_imag_opt)
    
        if modo == 1:
            line1, = axs[0].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
            line1.set_dashes([2, 2, 10, 2])  # 2pt line, 2pt break, 10pt line, 2pt break
            axs[0].set_ylim([mini-tol,maxi+tol])
            axs[0].set_yticks(ticks_y1_kz3)
            
        elif modo == 2:
            axs[1].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
            axs[1].set_ylim([mini-tol*1.5*1e-2,maxi+tol*1.5*1e-2])
            axs[1].set_yticks(ticks_y2_kz3)
            
        elif modo == 3:
            axs[2].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
            axs[2].set_ylim([mini-tol,maxi+tol])
            axs[2].set_yticks(ticks_y3_kz3)
        elif modo == 4:
            axs[3].plot(list_mu_opt,epsi1_imag_opt,symbols[k],lw =lw,color = colors[k],label = r'$n=$%i'%(k+1))
            axs[3].set_ylim([mini-tol,maxi+tol])
            axs[3].set_yticks(ticks_y4_kz3)
        k = k + 1
    
    
    for i in [0,1,2,3]:
        axs[i].minorticks_on()
        axs[i].tick_params(labelsize = tamnum,pad = pad)
        axs[i].set_ylabel(labely,fontsize = tamletra,labelpad =labelpady)
    
    axs[3].set_xlabel(labelx,fontsize = tamletra,labelpad =labelpadx) 
    
    #fig.legend(handles=patchList2,loc=[0.145,0.87],ncol=4,fontsize=tamlegend,frameon=0,handletextpad=0.5) 
#    fig.legend(loc = loc1, ncol = 4,markerscale=1,fontsize=tamlegend, columnspacing = 2,frameon=0,handletextpad=0.1)
    fig.legend(loc = loc1,ncol = 4,markerscale=1,fontsize=tamlegend, 
                      columnspacing = columnspace,frameon=0,handletextpad=0.2, handlelength=length_marker)    
    # s1,s2 = 'a)', 'b)'
    # xs = 0.133
    # fig.text(xs,0.83,s1,fontsize = tamletra) 
    # fig.text(xs,0.422,s2,fontsize = tamletra) 
    
#    fig.subplots_adjust(hspace=0.2)
    if save_graphs==1:
        plt.tight_layout(h_pad = 0.2)
        os.chdir(path_load)
        plt.savefig('Im_epsi1_vs_mu_kz%.4f.png'%(kz) ,format = 'png',dpi = dpi)

#%%    

"""