import numpy as np

#Planck
Nfreqs_Planck = 8
freqs_Planck = []
freqs_Planck.append('030')
freqs_Planck.append('044')
freqs_Planck.append('070')
freqs_Planck.append('100')
freqs_Planck.append('143')
freqs_Planck.append('217')
freqs_Planck.append('353')
freqs_Planck.append('545')
freqs_Planck_float = np.array([30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0, 545.0])
# Planck noise
noise_arr_Planck = np.zeros(Nfreqs_Planck)
noise_arr_Planck[0] = 195.079975053 #uK-arcmin, from Table 7 (first column) of https://arxiv.org/pdf/1502.01585.pdf -- converted via sqrt(3224.4*(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
noise_arr_Planck[1] = 226.090506617 # converted via sqrt(4331.0*(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
noise_arr_Planck[2] = 199.09525581 # converted via sqrt(3358.5*(4*Pi*(180/Pi)^2*60^2/(12*1024^2))) assuming Nside=1024
noise_arr_Planck[3] = 77.4 #uK-arcmin, from Table 6 of https://arxiv.org/pdf/1502.01587v2.pdf
noise_arr_Planck[4] = 33.0
noise_arr_Planck[5] = 46.8
noise_arr_Planck[6] = 153.6
noise_arr_Planck[7] = 0.78 * 1e-3 * 0.01723080316 * 1e6 * 60 #kJy/sr * deg --> converted to uK-arcmin
# Planck noise -- pol.
noise_arr_Planck_pol = np.zeros(Nfreqs_Planck)
noise_arr_Planck_pol[0] = 276.262025938 #uK-arcmin, from Table 7 (first column) of https://arxiv.org/pdf/1502.01585.pdf -- converted via sqrt( sqrt(6467.8*6465.1) *(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
noise_arr_Planck_pol[1] = 334.480764849 # converted via sqrt( sqrt(10088.0*8906.9) *(4*Pi*(180/Pi)^2*60^2/(12*1024^2)))
noise_arr_Planck_pol[2] = 282.141687884 # converted via sqrt( sqrt(6775.7*6713.7) *(4*Pi*(180/Pi)^2*60^2/(12*1024^2))) assuming Nside=1024
noise_arr_Planck_pol[3] = 117.6 #uK-arcmin, from Table 6 of https://arxiv.org/pdf/1502.01587v2.pdf
noise_arr_Planck_pol[4] = 70.2
noise_arr_Planck_pol[5] = 105.0
noise_arr_Planck_pol[6] = 438.6
# Planck beams
FWHM_arr_Planck = np.zeros(Nfreqs_Planck)
FWHM_arr_Planck[0] = 32.239
FWHM_arr_Planck[1] = 27.005
FWHM_arr_Planck[2] = 13.252
FWHM_arr_Planck[3] = 9.69 #arcmin, from Table 6 of https://arxiv.org/pdf/1502.01587v2.pdf
FWHM_arr_Planck[4] = 7.30
FWHM_arr_Planck[5] = 5.02
FWHM_arr_Planck[6] = 4.94
FWHM_arr_Planck[7] = 4.83 #arcmin

# convert to sigma in radians
sigma_arr_Planck = FWHM_arr_Planck / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
# Planck noise power spectra
MAX_NOISE=1.e9
ell_max=1e4
delta_ell=1
ell = np.arange(0,ell_max+1,delta_ell)
PS_noise_Planck = np.zeros((Nfreqs_Planck,np.int(ell_max)+1))
PS_noise_Planck_pol = np.zeros((Nfreqs_Planck,np.int(ell_max)+1))
for i in range(Nfreqs_Planck):
    PS_noise_Planck[i] = (noise_arr_Planck[i] * (1.0/60.0) * (np.pi/180.0))**2.0 * np.exp( ell*(ell+1)* sigma_arr_Planck[i]**2. ) #square to get the white-noise level -- see e.g. Eq. 2.32 (pg 15) of https://arxiv.org/pdf/1509.06770v4.pdf
    #PS_noise_Planck_pol[i] = (noise_arr_Planck_pol[i] * (1.0/60.0) * (np.pi/180.0))**2.0 * np.exp( ell*(ell+1)* sigma_arr_Planck[i]**2. ) #square to get the white-noise level -- see e.g. Eq. 2.32 (pg 15) of https://arxiv.org/pdf/1509.06770v4.pdf
    # handle overflow due to high noise at high ell
    PS_noise_Planck[i][(np.where(PS_noise_Planck[i] > MAX_NOISE))[0]] = MAX_NOISE
    #PS_noise_Planck_pol[i][(np.where(PS_noise_Planck_pol[i] > MAX_NOISE))[0]] = MAX_NOISE


#added below
for i in range(Nfreqs_Planck):
    with open(f'planck_noise/noise_{freqs_Planck[i]}GHz.txt', 'w') as f:
        for l in ell:
            f.write(f'{l} ')
        f.write('\n')
        for Cl in PS_noise_Planck[i]:
            f.write(f'{Cl} ')
