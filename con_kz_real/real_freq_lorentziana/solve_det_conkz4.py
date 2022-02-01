"""
Created on Wed May 20 08:58:35 2020
@author: leila
Diferencia con find_Lambda_conkz.py:
Cambiar las funciones de Bessel
(que sean J y Hankel, al igual que en gn)
determinante de 4x4
graficar el log|det(kz)| en mapa de color en funcion de 
omega en THz y del potencial quimico mu

resolver usando 3 variables: kz, omegaTHz y eta, para un mu0 fijo. Con el mu0 intentar llegar a soluciones cerca de 
freq1 y freq2 del medio activo.   
"""

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# from scipy.optimize import minimize  
from scipy.optimize import minimize_scalar

#%% 
save_data_opt = 1
save_graphs = 1

## minimizar con 3 variables: omegaTHz,kz y eta ############################

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq_lorentziana','')
path_graphene = path_basic2.replace('/' + 'con_kz_real','') 

try:
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante
except ModuleNotFoundError:
    print('det_conkz.py no se encuentra en ' + path_basic2)
    path_basic2 = input('path de la carpeta donde se encuentra det_conkz.py')
    sys.path.insert(1, path_basic2)
    from det_conkz import determinante

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
aux_cte = c*1e-12

#%%

print('Definir parametros del problema')

R = 0.5              #micrones
modo = 1
freq1 = 3.15 #en THz freq1 del medio activo
freq2 = 5.09 #en THz freq2 del medio activo
omegaTHz1 = freq1*2*np.pi
omegaTHz2 = freq2*2*np.pi
n = 251

#%%

########### parte 1 del barrido en mu ###################################################
list_mu0 = np.linspace(0.1,0.5,5) # para un mu0 y eta0 fijos
# list_mu0 = [0.53]
list_kz = np.linspace(0.1, 2.6, n)
list_omegaTHz = np.linspace(5,50,n)
list_cond_init = [20,2.15,0.8] # omegaTHz,kz y eta
########### parte 2 del barrido en mu ###################################################
# list_mu0 = np.linspace(0.6,0.9,4) # para un mu0 y eta0 fijos
# list_kz = np.linspace(0.2, 2.7, n)
# list_omegaTHz = np.linspace(5,60,n)
# list_cond_init = [30]
#######################################################################

#%%        

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16
labelpady = 0
labelpadx = 0
pad = 0
lw = 1
    
tol_NM = 1e-10
ite_NM = 4000

def graph(title,labelx,labely,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad):
    plt.figure(figsize=tamfig)
    plt.xlabel(labelx,fontsize=tamletra,labelpad =labelpadx)
    plt.ylabel(labely,fontsize=tamletra,labelpad =labelpady)
    plt.tick_params(labelsize = tamnum, pad = pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    return 

#%%

path_save = path_basic + '/' +  r'R_%.2f' %(R)
info = 'modo = %i, R = %.2f$\mu$m' %(modo,R)  

if save_data_opt==1:
    if not os.path.exists(path_save):
        print('Creating folder to save data')
        os.mkdir(path_save)
        
labelx = '$\omega$ [THz]'
labely1 = 'log|det|'

#%% 1 

labely2 = '$k_z$ [$\mu m^{-1}$]'
mu_opt = []
bnds = ((10,100),(0,4),(0.1,0.9))
for mu0 in list_mu0:
        
    inf = '$\mu_c$ = %.3f eV' %(mu0) + info + ', ' + name_this_py
    title = '$\mu_c$ = %.3f eV' %(mu0) + info + '\n' + name_this_py
    inf_fig = 'det_modo%i_mu%.4f' %(modo,mu0)
        
#
    print('Minimizar |det| con mu = %.1f eV' %(mu0))
    
    # list_cond_init  = [43, 59, 73]
    
    omegaTHz_opt = []
    kz_opt = []
    eta_opt = []
    eq_det = []

    # bnds = (15,100)
    
    def det_3variable(x):
        omegaTHz,kz,eta = x
        fTHz = omegaTHz/(2*np.pi)
        epsi1 = epsilon(fTHz, eta)
        omegac = omegaTHz/aux_cte
    
        rta = determinante(kz,omegac,epsi1,modo,R,mu0)
        return np.abs(rta)
    
#    res = minimize(det_3variable, list_cond_init, method='Nelder-Mead', tol=tol_NM, 
#                      options={'maxiter':ite_NM})
 
    
#    if res.message == 'Optimization terminated successfully.':

    res = minimize_scalar(det_3variable, method='bounded', bracket=None,  bounds= bnds, tol=tol_NM)

    if res.success == True:

        value = det_3variable(res.x)
        eq_det.append(value)
        mu_opt.append(mu0)
        omegaTHz_opt.append(res.x[0])
        kz_opt.append(res.x[1])
        eta_opt.append(res.x[2])
        
    cond_inicial = res.x
#

    print('Graficar |det| vs omega y kz con mu = %.1f eV' %(mu0))
    
    eta = res.x[2]
    def det_2variable(omegaTHz,kz):
        fTHz = omegaTHz/(2*np.pi)
        epsi1 = epsilon(fTHz, eta)
        omegac = omegaTHz/aux_cte
        rta = determinante(kz,omegac,epsi1,modo,R,mu0)
        return np.log10(np.abs(rta))
    
    #
    graph(title,labelx,labely2,tamfig,tamtitle,tamletra,tamnum,labelpadx,labelpady,pad)
    plt.title(title,fontsize=int(tamtitle*0.9))
    # im = plt.imshow(Z, extent = limits,  cmap='RdBu', interpolation='bilinear')
    
    X, Y = np.meshgrid(list_omegaTHz, list_kz, sparse=True)
    f = np.vectorize(det_2variable)
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
    plt.plot(omegaTHz_opt,kz_opt, 'g.',ms =11)  
    ejey1 = np.linspace(np.min(kz_opt),np.max(kz_opt),10)
    plt.plot(omegaTHz1*np.ones(10),ejey1,'-k',lw = 2.5)
    plt.plot(omegaTHz2*np.ones(10),ejey1,'-k',lw = 2.5)
    # cbar.set_ticks(tick_locations)
    cbar.ax.tick_params(labelsize=tamnum)
    cbar.set_label('log|det|',fontsize=tamlegend)
    if save_graphs==1:
        os.chdir(path_save)
        plt.savefig(inf_fig + '.png', format='png')


if save_data_opt == 1:
    os.chdir(path_save)
    print('Guardar data de minimizacion en .txt')

    tabla = np.array([mu_opt,omegaTHz_opt,kz_opt,eta_opt,eq_det])
    tabla = np.transpose(tabla)
    labeltxt = inf_fig + '.txt'
    header1 = inf
    header2 = 'mu [eV]    Omega [THz]    kz [1/micrones]    eta    |det|' + header1
    np.savetxt('opt_lorentz_conkz' + labeltxt, tabla, fmt='%1.11e', delimiter='\t', header = header2)

#%%

# del kz_opt, omegaTHz_opt, eta_opt, eq_det, list_kz, list_mu0

#%%

"""
print('')
print('Graficar |det| vs omega y mu para kz fijo para el modo = %i' %(modo))
       
def det_2variables(omegaTHz,mu):
    fTHz = omegaTHz/(2*np.pi)
    epsi1 = epsilon(fTHz, eta)
    omegac = omegaTHz/aux_cte
    
    rta = determinante(kz0,omegac,epsi1,modo,R,mu)
    return np.log10(np.abs(rta))
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
"""

#%%