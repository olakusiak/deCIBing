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
