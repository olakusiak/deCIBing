import numpy as np
import matplotlib.pyplot as plt

"""
Constants
"""

TCMB = 2.726 #Kelvin
TCMB_uK = 2.726e6 #micro-Kelvin
hplanck=6.62606896e-34
kboltz=1.3806504e-23
clight=299792458.0


"""
conversion functions
"""

def tSZ_spectral_funct_at_nu_in_GHz(nu_in_GHz):
        T_cmb_uk = 2.726e-6
        frequency_in_Hz = nu_in_GHz*1e9
        x = hplanck*frequency_in_Hz/(kboltz*TCMB)
        gNU = x*(1./np.tanh(x/2.))-4.
        return T_cmb_uk * gNU * 1e6 # class_sz curves come in y*1.e-6 units, multiply by * 10e6

def dBnudT(nu_ghz):
    nu = 1.e9*np.asarray(nu_ghz)
    X = hplanck*nu/(kboltz*TCMB)
    return (2.*hplanck*nu**3.)/clight**2. * (np.exp(X))/(np.exp(X)-1.)**2. * X/TCMB_uK

# conversion from specific intensity to Delta T units (i.e., 1/dBdT|T_CMB)
#   i.e., from W/m^2/Hz/sr (1e-26 Jy/sr) --> uK_CMB
#   i.e., you would multiply a map in 1e-26 Jy/sr by this factor to get an output map in uK_CMB
def ItoDeltaT(nu_ghz):
    return 1./dBnudT(nu_ghz)

def Jysr_to_uKcmb(nu_ghz):
    return dBnudT(nu_ghz) / 1e-26
"""
Radio model
"""

#radio contribution parameters
#from Choi et al.
nu0_radio_ghz =150.0  #radio pivot frequency [GHz]
beta_radio =-0.5   #radio power-law index
ell0=3000
A_s=3.74

def sed_radio(nu_ghz, nu0_radio_ghz,beta_radio):
    nu = 1.e9*np.asarray(nu_ghz).astype(float)
    nu0_radio = nu0_radio_ghz*1.e9
    resp = (nu/nu0_radio)**beta_radio * (ItoDeltaT(np.asarray(nu_ghz).astype(float))/ItoDeltaT(nu0_radio_ghz))
    return resp

def model_radio(nu_ghz1,nu_ghz2, nu0_radio_ghz, beta_radio, A_s, ell0, ell_array):
    sed1 = sed_radio(nu_ghz1, nu0_radio_ghz,beta_radio)
    sed2 = sed_radio(nu_ghz2, nu0_radio_ghz,beta_radio)
    poisson_radio = A_s *(ell_array*(ell_array+1)/ell0/(ell0+1)) * sed1*sed2
    return poisson_radio
"""
CIB flux cut
"""

cib_flux = {}
# Planck flux cut, Table 1 in https://arxiv.org/pdf/1309.0382.pdf
cib_flux['100'] = 400
cib_flux['143'] = 350
cib_flux['217'] = 225
cib_flux['353'] = 315
cib_flux['545'] = 350
cib_flux['857'] = 710
cib_flux['3000'] = 1000
#SO
cib_flux['93'] = 7
cib_flux['145'] = 15
cib_flux['225'] = 20
cib_flux['280'] = 25

def make_flux_cut_list(cib_flux, nu_list):
    cib_flux_list = []
    keys = list(cib_flux.keys())
    for i,nu in enumerate(nu_list):
        if str(nu) in keys:
            cib_flux_list.append(cib_flux[str(nu)])
            print(cib_flux[str(nu)])
        else:
            cib_flux_list.append(0.3)
            print(0)
            print("Flux cuts: ", cib_flux_list)
    return cib_flux_list





def cib_spectral_response(freqs, Tdust_CIB,beta_CIB ): #input frequency in GHz
    # from pyilc
    # CIB = modified blackbody here
    # N.B. overall amplitude is not meaningful here; output ILC map (if you tried to preserve this component) would not be in sensible units

    TCMB = 2.726 #Kelvin
    TCMB_uK = 2.726e6 #micro-Kelvin
    hplanck=6.626068e-34 #MKS
    kboltz=1.3806503e-23 #MKS
    clight=299792458.0 #MKS

    # function needed for Planck bandpass integration/conversion following approach in Sec. 3.2 of https://arxiv.org/pdf/1303.5070.pdf
    # blackbody derivative
    # units are 1e-26 Jy/sr/uK_CMB
    def dBnudT(nu_ghz):
        nu = 1.e9*np.asarray(nu_ghz)
        X = hplanck*nu/(kboltz*TCMB)
        return (2.*hplanck*nu**3.)/clight**2. * (np.exp(X))/(np.exp(X)-1.)**2. * X/TCMB_uK

    # conversion from specific intensity to Delta T units (i.e., 1/dBdT|T_CMB)
    #   i.e., from W/m^2/Hz/sr (1e-26 Jy/sr) --> uK_CMB
    #   i.e., you would multiply a map in 1e-26 Jy/sr by this factor to get an output map in uK_CMB
    def ItoDeltaT(nu_ghz):
        return 1./dBnudT(nu_ghz)



#     Tdust_CIB = 24.0       #CIB effective dust temperature [K] (Table 9 of http://www.aanda.org/articles/aa/pdf/2014/11/aa22093-13.pdf)
#     beta_CIB = 1.2         #CIB modified blackbody spectral index (Table 9 of http://www.aanda.org/articles/aa/pdf/2014/11/aa22093-13.pdf ; Table 10 of that paper contains CIB monopoles)
    nu0_CIB_ghz = 353.0    #CIB pivot frequency [GHz]
    kT_e_keV = 5.0         #electron temperature for relativistic SZ evaluation [keV] (for reference, 5 keV is rough temperature for a 3x10^14 Msun/h cluster at z=0 (e.g., Arnaud+2005))
    nu0_radio_ghz = 150.0  #radio pivot frequency [GHz]
    beta_radio = -0.5      #radio power-law index

    nu_ghz = freqs
    nu = 1.e9*np.asarray(nu_ghz).astype(float)
    X_CIB = hplanck*nu/(kboltz*Tdust_CIB)
    nu0_CIB = nu0_CIB_ghz*1.e9
    X0_CIB = hplanck*nu0_CIB/(kboltz*Tdust_CIB)
    resp = (nu/nu0_CIB)**(3.0+(beta_CIB)) * ((np.exp(X0_CIB) - 1.0) / (np.exp(X_CIB) - 1.0)) * (ItoDeltaT(np.asarray(nu_ghz).astype(float))/ItoDeltaT(nu0_CIB_ghz))
    #resp[np.where(nu_ghz == None)] = 0. #this case is appropriate for HI or other maps that contain no CMB-relevant signals (and also no CIB); they're assumed to be denoted by None in nu_ghz
    return resp

def plotfunction(linear=False, **kwargs):
    plt.figure(figsize=(10,6))
    plt.title(r"", fontsize=30,  **kwargs )
    plt.xlabel(r"$\ell$", size=30)
    plt.yticks(fontsize=20)
    plt.xticks(fontsize=20)
    plt.grid(which='both',alpha=0.4)
    plt.yscale("log")
    plt.xscale("log")
    if linear == True:
        plt.yscale("linear")
        plt.xscale("linear")

"""
CIB shot noise
"""

cib_sn = {
'93':
    {'93x93':0.10,
     '93x100':0.12,
     '93x143':0.34,
     '93x145':0.36,
     '93x217':1.22,
     '93x225': 1.35,
     '93x280':2.46,
     '93x353':4.41,
     '93x545':9.72,
    },
'100':
    {'100x93':0.12,
     '100x100': 0.15,
     '100x143': 0.42,
     '100x145':0.44,
     '100x217': 1.50,
     '100x225':1.66,
     '100x280':3.02,
     '100x353': 5.40,
     '100x545':12,
    },

'143':
    {
     '143x93':0.34,
     '143x100':0.42,
     '143x143': 1.20,
     '143x145':1.25,
     '143x217':4.30,
     '143x225': 4.75,
     '143x280':8.54,
     '143x353':15.00,
     '143x545':35.00,
    },
'145':
    {
     '145x93': 0.36,
     '145x100':0.44,
     '145x143':1.25,
     '145x145':1.31,
     '145x217':4.49,
     '145x225': 4.96,
     '145x280': 8.92,
     '145x353':15.69,
        '145x545':36.58,
    },
'217':
    {
     '217x93':1.22,
     '217x100':1.50,
     '217x143': 4.30,
     '217x145':4.49,
     '217x217':16.00,
     '217x225':17.75,
     '217x280':32.59,
     '217x353':59.00,
    '217x545': 135.0,
    },
'225':
    {
     '225x93': 1.35,
     '225x100':1.66,
     '225x143': 4.75,
     '225x145':4.96,
     '225x217':17.75,
     '225x225': 19.69,
     '225x280':36.18,
     '225x353':65.57 ,
    '225x545':150.50,
    },
'280':
    {
     '280x93':2.46,
     '280x100':3.02,
     '280x143':8.54,
     '280x145':8.92,
     '280x217':32.59,
     '280x225':36.18,
     '280x280':66.85,
     '280x353':122.05,
    '280x545':286.25,
    },
'353':
    {
     '353x93':4.41,
     '353x100':5.40,
     '353x143': 15.00,
     '353x145':15.69,
     '353x217':59,
     '353x225': 65.57,
     '353x280':122.05,
     '353x353':225.00,
        '353x545':543.0,
    },

'545':
    {
     '545x93': 9.72,
     '545x100': 12.00,
     '545x143': 35.00,
     '545x145': 36.58,
     '545x217':135,
     '545x225': 150.50,
     '545x280':286.25,
     '545x353': 543,
     '545x545': 1454.00,
    },
}
