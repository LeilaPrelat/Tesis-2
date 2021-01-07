# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 13:19:41 2020

@author: Usuario
"""

"""
Las formulas de la seccion 4.2 del cuaderno corto
"Reescribir $\omega_n$ separando parte real e imaginaria"
son las ecuaciones de quasi approximation sin Im(omega/c) 
(sin perdidas)

"""

import numpy as np
import os 
import sys

#%%

name_this_py = os.path.basename(__file__)
path = os.path.abspath(__file__) #path absoluto del .py actual
path_basic = path.replace('/' + name_this_py,'')
path_basic2 = path_basic.replace('/' + 'real_freq','')
path_graphene = path_basic.replace('/sin_kz_inv/dispersive/real_freq','') 

try:
    sys.path.insert(1, path_graphene)
    from constantes import constantes
except ModuleNotFoundError:
    print('constantes.py no se encuentra en el path_basic definido/carpeta de trabajo')
    path_graphene3 = input('path de la carpeta donde se encuentra constantes.py')
    sys.path.insert(1, path_graphene3)
    from constantes import constantes

pi,hb,c,alfac,hbargama,mu1,mu2,epsi2 = constantes()
 
#%%

# print('Aprox cuasi-estacionaria')    
# print('Definir la frecuencia del modo enesimo: Eq 16 (AA), medio 1 no dispersivo')

def im_epsi1_cuasi(omegac,Ep,epsiinf_DL,gamma_DL,modo,R,hbaramu):
 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrones
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    modo: modo
    R: radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    im(epsilon1) para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 7.1 cuaderno corto
    plasmones sin kz - invisibilidad---> RESULTA SER LO MISMO QUE LA SECCION 1 CASO SIN KZ)

    """
    TK = 300               	# temperature in K
    # hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K 

    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02n = -4*pi*1j*intra*modo/R #omega0n**2
    omega0n = omega02n**(1/2)
    E0n = hb*omega0n
    
    den = 2*(Ep**2 + E0n**2)
    num = gamma_DL*(Ep**2) + hbargama*(E0n**2)
    E2n = num/den

    cte = 2*(epsiinf_DL - epsi2)
    E = omegac*c*hb
    den = E2n**2 - E**2
    num = cte*E*E2n
    rta = num/den
    if rta.imag == 0:
    	rta = rta.real
    return rta

def omegac_cuasi(modo,Ep,epsiinf_DL,gamma_DL,R,hbaramu):
    """
    Parameters
    ----------
    Ep : hbar*omega_p (Drude-Lorentz) de permeabilidad del medio 1
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    modo : modo
    R : radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    omega/c en 1/micrones para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 7.1 cuaderno corto
    plasmones sin kz - invisibilidad---> RESULTA SER CASI LO MISMO: 
        cambiar epsiinf_DL + epsi2 por epsiinf_DL - epsi2)
    """
    TK = 300               	# temperature in K
    # hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  
    
    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02n = -4*pi*1j*intra*modo/R #omega0n**2
    omega0n = omega02n**(1/2)
    E0n = hb*omega0n
    # gamma_c = hbargama/hb
    den = 2*(Ep**2 + E0n**2)
    num = gamma_DL*(Ep**2) + hbargama*(E0n**2)
    E2n = num/den

    cte2 = 2*(epsiinf_DL - epsi2)
    
    term0 = (Ep**2 + E0n**2)
    term1 = -E2n**2 + term0/cte2
    term2 = (term0*(term0 - 4*cte2*E2n**2))**(1/2)
    rta_energy2 =  term1 + term2/cte2
    
    rta = rta_energy2**(1/2) #energy
    omegac_n = rta/(hb*c) #energy/hb da omega --> luego divido por c
    if omegac_n.imag == 0:
    	omegac_n = omegac_n.real
    return omegac_n #omega/c
        
#%%

def im_epsi1_cuasi_aprox(omegac,epsiinf_DL,gamma_DL,modo,R,hbaramu):
 
    """
    Parameters
    ----------
    omegac : omega/c en 1/micrones
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    modo: modo
    R: radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    im(epsilon1) para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 1.7 cuaderno corto)
    despreciando omega_p

    """
    # TK = 300               	# temperature in K
    # hbargama = 0.0001      	# collision frequency in eV
    # akb = 8.6173324E-5           ### Boltzman constant in eV/K 

    # TKmu = TK*akb/hbaramu
    # aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    # intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    # intra = intra*alfac*c
    # omega02n = -4*pi*1j*intra*modo/R #omega0n**2
    # omega0n = omega02n**(1/2)
    # E0n = hb*omega0n
    
    E2n = hbargama/2 #despreciando omega_p

    cte = 2*(epsiinf_DL - epsi2)
    E = omegac*c*hb
    den = E2n**2 - E**2
    num = cte*E*E2n
    rta = num/den
    if rta.imag == 0:
    	rta = rta.real
    	
    return rta

def omegac_cuasi_aprox(modo,epsiinf_DL,gamma_DL,R,hbaramu):
    """
    Parameters
    ----------
    epsiinf_DL: epsilon_infinito (Drude-Lorentz) de permeabilidad del medio 1
    gamma_DL : gamma (Drude-Lorentz) de permeabilidad del medio 1
            en unidades eV de energia
    modo : modo
    R : radio del cilindro en micrones
    hbaramu : potencial quimico del grafeno mu_c en eV

    Returns
    -------
    omega/c en 1/micrones para la aprox QE tal que no hay perdidas
    (se calculo el im(epsilon1) tal que im(omega/c) = 0 , ver seccion 4.2 cuaderno corto)
    """
    TK = 300               	# temperature in K
    # hbargama = 0.0001      	# collision frequency in eV
    akb = 8.6173324E-5           ### Boltzman constant in eV/K  
    
    TKmu = TK*akb/hbaramu
    aux = np.exp(1/(2*TKmu)) + np.exp(-1/(2*TKmu))
    intra = 2j*(1/pi)*akb*TK*np.log(aux)/(hb)
    intra = intra*alfac*c
    omega02n = -4*pi*1j*intra*modo/R #omega0n**2
    omega0n = omega02n**(1/2)
    E0n = hb*omega0n
    # gamma_c = hbargama/hb

    E2n = hbargama/2 #despreciando omega_p

    cte2 = 2*(epsiinf_DL - epsi2)
    
    term0 = (E0n**2)
    term1 = -E2n**2 + term0/cte2
    term2 = (term0*(term0 - 4*cte2*E2n**2))**(1/2)
    rta_energy2 =  term1 + term2/cte2
    
    rta = rta_energy2**(1/2) #energy
    omegac_n = rta/(hb*c) #energy/hb da omega --> luego divido por c
    if omegac_n.imag == 0:
    	omegac_n == omegac_n.real
    return omegac_n #omega/c
        
#%%
