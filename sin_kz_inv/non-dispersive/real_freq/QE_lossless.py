# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Las formulas de la seccion 7.1 (caso sin kz invisibilidad) del cuaderno corto
"Reescribir $\omega_n$ separando parte real e imaginaria"
son las ecuaciones de quasi approximation sin Im(omega/c) 
(sin perdidas)

+ formulas despreciando el gamma_c**2 

"""

import numpy as np
import os 
import sys

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')
path_graphene = path_basic.replace('/sin_kz/non-dispersive/real_freq','') 

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

# print('Aprox cuasi-estacionaria')    
# print('Definir la frecuencia del modo enesimo: Eq 16 (AA), medio 1 no dispersivo')

def im_epsi1_cuasi(omegac,modo,R,hbaramu):
 
    """
    Parameters
    viene de jn(x)--> n/x1**2
    y de hn(x)--> -n/x2**2
    ----------
    omegac : omega/c en 1/micrones
    modo: modo
    R: radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    im(epsilon1) para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 7.1 cuaderno corto
    plasmones sin kz - invisibilidad)

    """
    TK = 300               	# temperature in K
    hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K 
    
    omega = omegac*c
    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu))+np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02 = -4*pi*1j*intra*modo/R

    gamma_c = hbargama/hb
    num = omega02*omega*gamma_c
    den = (omega**2 + (gamma_c/2)**2)**2

    rta = -num/den
    if rta.imag == 0:
        rta = rta.real
    return rta 

def omegac_cuasi(modo,R,re_epsi1,hbaramu):
    """
    Parameters
    viene de jn(x)--> n/x1**2
    y de hn(x)--> -n/x2**2
    ----------
    modo : modo
    R : radio del cilindro en micrones
    re_epsi1 : re(epsilon1)
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    omega/c en 1/micrones para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 1.7 cuaderno corto)
    """
    TK = 300               	# temperature in K
    hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  

    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu))+np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02 = -4*pi*1j*intra*modo/R # no es compleja la frecuencia porque intra tiene otro 1j
    omega0 = omega02**(1/2)

#     cte = re_epsi1 - epsi2
#     gamma_c = hbargama/hb
#     gamma_c2 = (gamma_c)**2
    
#     term1 = -0.25*gamma_c2 
#     term2 = 0.5*omega0/cte
#     term3 = (omega02 + 2*gamma_c2*cte)**(1/2)
    
#     rta1 = term1 + term2*(-omega0 + term3) # -- > me quedo con la rama de +
# #    rta1 = term1 + term2*(-omega0 - term3) # -- > me quedo con la rama de -

#     omega_n = (rta1)**(1/2)
#     omegac_n = omega_n/c
#     if omegac_n.imag == 0:
#     	omegac_n = omegac_n.real
#     return omegac_n #omega/c

    cte = re_epsi1 - epsi2  #UNICO CAMBIO QUE HAY QUE HACER
    gamma_c = hbargama/hb
    gamma_c2 = (gamma_c)**2
    
    term1 = -0.25*gamma_c2 
    term2 = 0.5*omega0/cte
    term3 = (omega02 - 2*gamma_c2*cte)**(1/2)
    
    rta1 = term1 + term2*(omega0 + term3)
    # rta2 = term1 + term2*(omega0 - term3)

    omega_n = (rta1)**(1/2)
    omegac_n = omega_n/c
    if omegac_n.imag == 0:
    	omegac_n = omegac_n.real
    return omegac_n #omega/c
        
#%%

def im_epsi1_cuasi_aprox(omegac,modo,R,hbaramu):
 
    """
    Parameters
    viene de jn(x)--> n/x1**2
    y de hn(x)--> -n/x2**2
    ----------
    omegac : omega/c en 1/micrones
    modo: modo
    R: radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    im(epsilon1) para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 1.7 cuaderno corto)
    despreciando el gamma_c**2

    """
    TK = 300               	# temperature in K
    hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K 
    
    omega = omegac*c
    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu))+np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02 = -4*pi*1j*intra*modo/R

    gamma_c = hbargama/hb   #en unidades de frecuencia
    num = omega02*gamma_c
    den = omega**3

    rta = num/den
    if rta.imag == 0:
        rta = rta.real
    return rta  

def omegac_cuasi_aprox(modo,R,re_epsi1,hbaramu):
    """
    Parameters
    viene de jn(x)--> n/x1**2
    y de hn(x)--> -n/x2**2
    ----------
    modo : modo
    R : radio del cilindro en micrones
    re_epsi1 : re(epsilon1)
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    omega/c en 1/micrones para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 1.7 cuaderno corto)
    despreciando el gamma_c**2
    """
    TK = 300               	# temperature in K
#    hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  

    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu))+np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02 = -4*pi*1j*intra*modo/R
    
    num = -omega02
    den = re_epsi1 - epsi2
    # rta2 = term1 + term2*(omega0 - term3)

    omega_n = (num/den)**(1/2)  #no se puede despreciar el gamma**2 en el caso de invsibilidad
    omegac_n = omega_n/c
    if omegac_n.imag == 0:
    	omegac_n = omegac_n.real
    
    return omegac_n #omega/c

#%%
