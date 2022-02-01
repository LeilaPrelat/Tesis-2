"""
Created on Sun Jan 16 08:29:58 2022

calcula la permitividad de un material de tres niveles con el modelo lorentziano

parámetros y esquema obtenidos de la referencia:
Pavlov, S. G., Deßmann, N., Pohl, A., Zhukavin, R. K., Klaassen, T. O., Abrosimov, N. V., ... & Hübers, H. W. (2020). Terahertz transient stimulated emission from doped silicon. APL Photonics, 5(10), 106102.

usaremos los niveles 
1 ≡ 1s(A1), 2 ≡ 1s(E),3 ≡ 2p0
@author: npassarelli

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

usar_THz = 1 # siempre = 1 porque no funciona la otra opcion (habria que revisar cuentas)
graficar = 0

if usar_THz == 1:
    cte = 1e12
else:
    cte = 1

#%%

def epsilon(freq, eta):
    """
    Parameters
    ----------
    freq : freq en THz (entre 1 y 10). OJO es f = omega/(2*pi)
    eta : inversion de poblacion en THz (entre
                                         0 y 1)

    Returns
    -------
    epsilon : modelo lorenztiano para epsilon
    """

    # =============================================================================
    # parameters
    # =============================================================================
    f_12 = 3.15
    f_13 = 8.24 
    f_23 = f_13-f_12
    
    #tau_32=200e-12 # s
    #tau_21=10e-12 # s
    t_c = 10e-9 #s          tiempo de vida del foton en el resonador, hace falta Purcell correction (ver referencia 10)
    
    gamma_32=1 #THz # a ojo
    gamma_21=1 #THz
    gamma_13=1 #THz
    
    N = 3e15 # cm^-3   densidad de emisores
    
    sigma_13 = 4e-15 #cm^2
    sigma_12 = sigma_13
    sigma_23 = sigma_13
    
    eps_h = 16 # permitividad del medio (Si)
    
    h_bar = 1.0546e-27 # reduced plank constant in cm2 g s-1
    c = 2.99792458e10 # vel luz en cm/s
    
    # =============================================================================
    # precalculo
    # =============================================================================
    n_h = eps_h**.5
    
    lambda_32 = c/(f_23*cte) # en THz
    k_32=2*np.pi/lambda_32
    
    lambda_21 = c/(f_12*cte) # en THz
    k_21 = 2*np.pi/lambda_21
    
    lambda_13 = c/(f_13*cte) # en THz
    k_13 = 2*np.pi/lambda_13
    
    
    mu_32=np.sqrt(sigma_13*3*h_bar*eps_h*gamma_32*cte/(8*np.pi*n_h*k_32))
    mu_21=np.sqrt(sigma_12*3*h_bar*eps_h*gamma_21*cte/(8*np.pi*n_h*k_21))
    
    #mu_32=(c/(2*np.pi*f_23*1e12))**1.5 * (3*h_bar/(8*np.pi*eps_h**2*t_c))**.5
    #mu_21=(c/(2*np.pi*f_12*1e12))**1.5 * (3*h_bar/(8*np.pi*eps_h**2*t_c))**.5
    
    
    eps_L32=(8*np.pi*N*mu_32**2/(h_bar*gamma_32*cte))
    eps_L21=(8*np.pi*N*mu_21**2/(h_bar*gamma_21*cte))
    
    
    lambda_13 = c/(f_13*cte)
    
    refind_im = N*sigma_13*lambda_13/(4*np.pi) 
    eps_L13 = 2*refind_im*np.sqrt(eps_h)
    
    
    # =============================================================================
    # main
    # =============================================================================


    def chi_dot_epsi(freq, freq_res, gamma, epsi_l):
        gamma_div = gamma*0.5
        num = (freq_res - freq + 1j*gamma_div)*gamma_div
        den = (freq - freq_res)**2 + gamma_div**2
        rta = epsi_l*(num/den)
        return rta
    
    xi_abs = chi_dot_epsi(freq,f_13,gamma_13,eps_L13)
    xi_inv = chi_dot_epsi(freq,f_23,gamma_32,eps_L32)
    xi_srs = chi_dot_epsi(freq,f_12,gamma_21,eps_L21)
    
    epsilon = eps_h + xi_abs*(1 - eta) - xi_inv*eta - xi_srs*eta

    return epsilon

#%%

if graficar == 1:
    
    nfreq = 250
    neta = 250
    
    if usar_THz == 1:
        freq = np.linspace(1,10,nfreq)
    else:
        freq = np.linspace(1e12,1e13,nfreq)
    
    eta = np.linspace(0,1,neta)
    X, Y = np.meshgrid(freq, eta, sparse=True)
    
    f = np.vectorize(epsilon)
    Z = f(X, Y)
    
    
    fig, ax = plt.subplots(1,1)
    
    vmax = np.max(np.abs(np.imag(Z)))
    
    pcm = ax.pcolor(freq, eta, np.imag(Z),
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=-vmax, vmax=vmax),   
    				cmap='PiYG')
    
    cbar = fig.colorbar(pcm, ax=ax, extend='neither',fraction=.2)
    cbar.set_label('Im $\epsilon$')
    
    #ax.set_aspect('equal', 'box')
    
    #plt.xticks([])
    #plt.yticks([])
    ax.set_ylabel('population inversion')
    ax.set_xlabel('f=$\omega/2\pi$ [THz]')
    ax.set_ylim(0,1)
    # ax.set_xlim(1,10)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.90, top=0.95)
    plt.savefig('epsilon-imag.png')
    # plt.show()
    # plt.close()
    
    
    fig, ax = plt.subplots(1,1)
    
    vmin=np.min(np.real(Z))
    vmax=np.max(np.real(Z))
    
    pcm = ax.pcolor(freq, eta, np.real(Z),
    #                   norm=colors.LogNorm(vmin=np.quantile(absE3, 0), vmax=np.quantile(absE3, 1) ),
    				   norm=colors.Normalize(vmin=vmin, vmax=vmax),   
    				cmap='PiYG')
    
    
    cbar = fig.colorbar(pcm, ax=ax, extend='neither',fraction=.2)
    cbar.set_label('Re $\epsilon$')
    
    #ax.set_aspect('equal', 'box')
    
    #plt.xticks([])
    #plt.yticks([])
    ax.set_ylabel('population inversion')
    ax.set_xlabel('f=$\omega/2\pi$ [THz]')
    ax.set_ylim(0,1)
    # ax.set_xlim(1,10)
    plt.subplots_adjust(left=0.2, bottom=0.15, right=0.90, top=0.95)
    plt.savefig('epsilon-real.png')
    # plt.show()
    # plt.close()

#%%
