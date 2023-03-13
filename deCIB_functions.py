import numpy as np

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
