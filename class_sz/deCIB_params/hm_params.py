common_settings={}

# common_settings['redshift_epsabs']= 1.0e-100
# common_settings['redshift_epsrel']= 1e-6
# common_settings['mass_epsabs']= 1.0e-100
# common_settings['mass_epsrel']= 1e-6
# common_settings['ndim_masses']= 500
# common_settings['ndim_redshifts']= 500
# common_settings['epsabs_L_sat']=  1e-100
# common_settings['epsrel_L_sat']=  1e-6
# common_settings['n_z_L_sat'] = 350
# common_settings['n_m_L_sat'] = 350
# common_settings['n_nu_L_sat'] = 350
# common_settings['n_m_dndlnM']=700
# common_settings['n_z_dndlnM']=700
# common_settings['n_m_pressure_profile' ] = 500
# common_settings['n_z_pressure_profile' ] = 500

common_settings['z_min']= 0.005 #tbd
common_settings['z_max']= 12 #tbd
common_settings['ell_max']= 11e3
common_settings['ell_min']= 2.0
common_settings['M_min']= 7.0e8 #tbd
common_settings['M_max']= 3.5e15 #tbd

### Precision
common_settings['redshift_epsabs']= 1.0e-40
common_settings['redshift_epsrel']= 0.0005
common_settings['mass_epsabs']= 1.0e-40
common_settings['mass_epsrel']= 0.0005
common_settings['ndim_masses']= 150
common_settings['ndim_redshifts']= 150
common_settings['epsabs_L_sat']=  1e-40
common_settings['epsrel_L_sat']=  1e-6

common_settings['dell'] = 500 # to be changed, for testing now

common_settings['P_k_max_h/Mpc']= 50.0
common_settings['k_min_for_pk_class_sz']= 0.0001
common_settings['k_max_for_pk_class_sz']= 10.0
common_settings['k_per_decade_class_sz']= 20.0

common_settings['hm_consistency']= 1
#common_settings['class_sz_verbose']= 1 #it breaks when this is 1 or 0 sometimes

### This needs to confirmed  !!!!!!!!!!!!!!!!!!!!!!!!
common_settings['delta for galaxies'] = "200c"
common_settings['delta for matter density'] = "200c"
common_settings['delta for electron pressure'] = "200c"
common_settings['delta for cib'] = "200c"
common_settings['mass function'] = 'T08M200c'
common_settings['concentration parameter'] = 'B13'

common_settings['k_pivot'] = 0.05 # what about these
common_settings['N_ncdm'] = 1
common_settings['N_ur'] = 2.0328
common_settings['m_ncdm'] = 0.06


"""
Planck 2018 cosmology
"""
p18_pdict={}
p18_pdict['omega_b'] = 0.02242
p18_pdict['omega_cdm'] = 0.11933
p18_pdict['h'] = 0.6766
p18_pdict['tau_reio'] = 0.0561
p18_pdict['ln10^{10}A_s'] = 3.047
#p18_pdict['sigma8'] = 0.8102
p18_pdict['n_s'] = 0.9665


#### tSZ
Mmin_websky_msun   = 1.3e12 # approximately websky minimum M200m value
Mmax_websky_msun   = 1e16
websky_tsz_pdict = {}
websky_tsz_pdict['pressure profile'] = 'B12'  # check source/input.c for default parameter values of Battaglia et al profile (B12)
websky_tsz_pdict['units for tSZ spectrum'] = 'dimensionless'
websky_tsz_pdict['n_ell_pressure_profile' ] = 100
websky_tsz_pdict['n_m_pressure_profile' ] = 100
websky_tsz_pdict['n_z_pressure_profile' ] = 100
websky_tsz_pdict['x_outSZ'] = 4.
websky_tsz_pdict['truncate_wrt_rvir'] =0
websky_tsz_pdict['pressure_profile_epsrel'] =1e-3
websky_tsz_pdict['pressure_profile_epsabs'] =1e-40
websky_tsz_pdict['M_min_tSZ']= Mmin_websky_msun*p18_pdict['h']
websky_tsz_pdict['M_max_tSZ']= Mmax_websky_msun*p18_pdict['h']

#### CIB
websky_cib_pdict = {}
websky_cib_pdict['Redshift evolution of dust temperature' ] =  0.2
websky_cib_pdict['Dust temperature today in Kelvins' ] = 20.7
websky_cib_pdict['Emissivity index of sed' ] = 1.6
websky_cib_pdict['Power law index of SED at high frequency' ] = 1.7 # not given in WebSky paper actually not relevant since we dont use high freqs in websky.
websky_cib_pdict['Redshift evolution of L − M normalisation' ] = 1.28
websky_cib_pdict['Most efficient halo mass in Msun' ] = 10.**12.3
#websky_cib_pdict['Normalisation of L − M relation in [Jy MPc2/Msun]' ] =  4.461102571695613e-07 # L0_websky, not given in WebSky paper
websky_cib_pdict['Normalisation of L − M relation in [Jy MPc2/Msun]'] =  5.061102571695613e-07
websky_cib_pdict['Size of of halo masses sourcing CIB emission' ] = 0.3
websky_cib_pdict['z_plateau_cib' ] = 2.

# M_min_HOD_cib is the threshold above which nc = 1] =
websky_cib_pdict[ 'M_min_HOD_cib' ] = 10.**10
websky_cib_pdict['M_min_subhalo_in_Msun' ] = 1e11 # if M_min_subhalo_in_Msun is given, then M_min_HOD_cib is ignored
websky_cib_pdict['use_nc_1_for_all_halos_cib_HOD'] = 1
websky_cib_pdict['sub_halo_mass_function' ] = 'JvdB14'
#Mass bounds
websky_cib_pdict['use_redshift_dependent_M_min'] = 1
websky_cib_pdict['full_path_to_redshift_dependent_M_min'] ='/Users/aleksandra/software/class_sz/sz_auxiliary_files/websky_halo_mass_completion_z_Mmin_in_Msun_over_h.txt'
websky_cib_pdict['M_max_cib' ] = 1e16*p18_pdict['h']

#for the monopole computation, is it necessary ??
websky_cib_pdict['freq_min'] = 2e1
websky_cib_pdict['freq_max'] = 4e3
websky_cib_pdict['dlogfreq' ] = 0.05


fiona_alina = {}
fiona_alina['mass function'] = 'T10'
fiona_alina['concentration parameter'] = 'D08'
fiona_alina['delta for cib'] = '200m'
fiona_alina['hm_consistency'] = 1
fiona_alina['damping_1h_term'] = 0


# parameters for Cosmology
#mass bounds
fiona_alina['M_min_cib'] = 1e8*p18_pdict['h']
fiona_alina['M_max_cib'] = 1e16*p18_pdict['h']

# redshift bounds
fiona_alina['z_min'] = 0.07
fiona_alina['z_max'] = 6. # fiducial for MM20 : 6
fiona_alina['freq_min'] = 10.
fiona_alina['freq_max'] = 5e4

# HOD parameters for CIB
fiona_alina['M_min_HOD_cib'] = pow(10.,10)
#fiona_alina['M1_prime_HOD'] =pow(10.,125.1536196)*fiona_alina['h']

# CIB parametes see McCarthy & Madhavacheril 2020
fiona_alina['Redshift evolution of dust temperature'] =  0.36
fiona_alina['Dust temperature today in Kelvins'] = 24.4
fiona_alina['Emissivity index of sed'] = 1.75
fiona_alina['Power law index of SED at high frequency'] = 1.7
fiona_alina['Redshift evolution of L − M normalisation'] = 3.6
fiona_alina['Most efficient halo mass in Msun'] = pow(10.,12.6)
#fiona_alina['Normalisation of L − M relation in [Jy MPc2/Msun]'] = 6.4e-8 #original
fiona_alina['Normalisation of L − M relation in [Jy MPc2/Msun]'] = 7.0e-8 #this work, rescaled
fiona_alina['Size of of halo masses sourcing CIB emission'] = 0.5

# precision parameters
fiona_alina['pressure_profile_epsabs'] = 1.e-8
fiona_alina['pressure_profile_epsrel'] = 1.e-3
# precision for redshift integal
fiona_alina['redshift_epsabs'] = 1e-40#1.e-40
fiona_alina['redshift_epsrel'] = 1e-4#1.e-10 # fiducial value 1e-8
# precision for mass integal
fiona_alina['mass_epsabs'] = 1e-40 #1.e-40
fiona_alina['mass_epsrel'] = 1e-4#1e-10
# precision for Luminosity integral (sub-halo mass function)
fiona_alina['L_sat_epsabs'] = 1e-40 #1.e-40
fiona_alina['L_sat_epsrel'] = 1e-3#1e-10

fiona_alina['z_max_pk'] = fiona_alina['z_max']
