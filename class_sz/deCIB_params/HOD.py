s_blue =  0.455
s_green =  0.648
s_red = 0.842

#### HOD PARAMS
###  mean
blue_hod_pdict = {}

# [0 : minimize] Parameter values at minimum:
#    weight  minuslogpost  alpha_s_HOD  sigma_log10M_HOD  logM1_prime  logM_min_HOD  x_out_truncated_nfw_profile_satellite_galaxies  logA_shot_noise  M1_prime_HOD     M_min_HOD  A_shot_noise  minuslogprior  minuslogprior__0         chi2  chi2__soliket.GXG_GXM_MXM_LOGSN_PS_Likelihood
# 0     1.0    694.825789      1.79262          0.958794     13.38958     12.137475                                        1.684816         -6.96772  2.452337e+13  1.372384e+12  1.077159e-07       2.980262          2.980262  1383.691053                                    1383.691053
blue_hod_pdict['galaxy_sample'] ='unwise'
blue_hod_pdict['galaxy_sample_id'] ='blue'
blue_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
blue_hod_pdict['M0_HOD']= 0  # Msun/h
blue_hod_pdict['M0 equal M_min (HOD)'] = 'no'
blue_hod_pdict['x_out_truncated_nfw_profile']= 1.0
blue_hod_pdict['M_min_gal']= 7.0e8
blue_hod_pdict['M_max_gal']= 3.5e15

blue_hod_pdict['alpha_s_HOD']=  1.061623
blue_hod_pdict['sigma_log10M_HOD']=   0.019955
blue_hod_pdict['M1_prime_HOD']= 10** 12.608837
blue_hod_pdict['M_min_HOD']= 10** 11.692076
blue_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=     1.8
A_shot_noise_blue = 10**(-7.062235)


#### HOD PARAMS
green_hod_pdict = {}
green_hod_pdict['galaxy_sample'] ='unwise'
green_hod_pdict['galaxy_sample_id'] ='green'
green_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
green_hod_pdict['M0_HOD']= 0
green_hod_pdict['M0 equal M_min (HOD)']= 'no'
green_hod_pdict['x_out_truncated_nfw_profile']= 1.0
green_hod_pdict['M_min_gal']= 7.0e8
green_hod_pdict['M_max_gal']= 3.5e15

green_hod_pdict['alpha_s_HOD']=    1.138029
green_hod_pdict['sigma_log10M_HOD']=   0.032215
green_hod_pdict['M1_prime_HOD']= 10**  13.171083
green_hod_pdict['M_min_HOD']= 10** 12.226724
green_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=      0.944354
A_shot_noise_green = 10**(  -6.815523)



#### HOD PARAMS
red_hod_pdict={}
red_hod_pdict['galaxy_sample'] ='unwise'
red_hod_pdict['galaxy_sample_id'] ='red'
red_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
red_hod_pdict['M0_HOD']= 0
red_hod_pdict['M0 equal M_min (HOD)'] : 'no'
red_hod_pdict['x_out_truncated_nfw_profile']= 1.0
red_hod_pdict['M_min_gal']= 7.0e8
red_hod_pdict['M_max_gal']= 3.5e15

red_hod_pdict['alpha_s_HOD']=    1.758605
red_hod_pdict['sigma_log10M_HOD']=        0.071076
red_hod_pdict['M1_prime_HOD']= 10** 13.567938
red_hod_pdict['M_min_HOD']= 10** 12.556522
red_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=  1.018972
A_shot_noise_red = 10**(-5.539891)


# #### HOD PARAMS
# ###  50 data points
# blue_hod_pdict = {}
#
# # [0 : minimize] Parameter values at minimum:
# #    weight  minuslogpost  alpha_s_HOD  sigma_log10M_HOD  logM1_prime  logM_min_HOD  x_out_truncated_nfw_profile_satellite_galaxies  logA_shot_noise  M1_prime_HOD     M_min_HOD  A_shot_noise  minuslogprior  minuslogprior__0         chi2  chi2__soliket.GXG_GXM_MXM_LOGSN_PS_Likelihood
# # 0     1.0    694.825789      1.79262          0.958794     13.38958     12.137475                                        1.684816         -6.96772  2.452337e+13  1.372384e+12  1.077159e-07       2.980262          2.980262  1383.691053                                    1383.691053
# blue_hod_pdict['galaxy_sample'] ='unwise'
# blue_hod_pdict['galaxy_sample_id'] ='blue'
# blue_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
# blue_hod_pdict['M0_HOD']= 0  # Msun/h
# blue_hod_pdict['M0 equal M_min (HOD)'] = 'no'
# blue_hod_pdict['x_out_truncated_nfw_profile']= 1.0
#
# blue_hod_pdict['alpha_s_HOD']=   1.0641872
# blue_hod_pdict['sigma_log10M_HOD']=  0.019417893
# blue_hod_pdict['M1_prime_HOD']= 10** 12.6061
# blue_hod_pdict['M_min_HOD']= 10**11.691596
# blue_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=  1.7994517
# A_shot_noise_blue = 10**(-7.0661348 )
#
#
# # blue_hod_pdict['alpha_s_HOD']= 1.79262
# # blue_hod_pdict['sigma_log10M_HOD']=  0.958794
# # blue_hod_pdict['M1_prime_HOD']= 10**13.38958
# # blue_hod_pdict['M_min_HOD']= 10**12.137475
# # blue_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']= 1.684816
# # A_shot_noise_blue = 1.077159 #e-07
#
# # ## ellmax=3000, 0covamt
# # blue_hod_pdict['alpha_s_HOD']=  1.333957
# # blue_hod_pdict['sigma_log10M_HOD']=  0.791019
# # blue_hod_pdict['M1_prime_HOD']= 10**13.144242
# # blue_hod_pdict['M_min_HOD']= 10**12.206059
# # blue_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']= 0.959523
# # A_shot_noise_blue = 0.6540711597277759 #e-07
#
# #### HOD PARAMS
# green_hod_pdict = {}
# green_hod_pdict['galaxy_sample'] ='unwise'
# green_hod_pdict['galaxy_sample_id'] ='green'
# green_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
# green_hod_pdict['M0_HOD']= 0
# green_hod_pdict['M0 equal M_min (HOD)']= 'no'
# green_hod_pdict['x_out_truncated_nfw_profile']= 1.0
#
# green_hod_pdict['alpha_s_HOD']=       1.131729
# green_hod_pdict['sigma_log10M_HOD']=     0.035472
# green_hod_pdict['M1_prime_HOD']= 10**  13.171437
# green_hod_pdict['M_min_HOD']= 10**  12.226522
# green_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=  0.935355
# A_shot_noise_green = 10**( -6.813347 )
#
# # green_hod_pdict['sigma_log10M_HOD'] = 0.60500869
# # green_hod_pdict['alpha_s_HOD'] =  1.2348499
# # green_hod_pdict['M1_prime_HOD']=  10**12.867959
# # green_hod_pdict['M_min_HOD']= 10**12.387252
# # green_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=  2.4961516
# # A_shot_noise_green = 1.3465079
#
#
#
# #### HOD PARAMS
# red_hod_pdict={}
# red_hod_pdict['galaxy_sample'] ='unwise'
# red_hod_pdict['galaxy_sample_id'] ='red'
# red_hod_pdict['UNWISE_dndz_file'] =  "/Users/aleksandra/software/class_sz/sz_auxiliary_files/UNWISE_galaxy_distributions/normalised_dndz_cosmos.txt"
# red_hod_pdict['M0_HOD']= 0
# red_hod_pdict['M0 equal M_min (HOD)'] : 'no'
# red_hod_pdict['x_out_truncated_nfw_profile']= 1.0
#
# red_hod_pdict['alpha_s_HOD']=      1.666696
# red_hod_pdict['sigma_log10M_HOD']=      0.041204
# red_hod_pdict['M1_prime_HOD']= 10**  13.173948
# red_hod_pdict['M_min_HOD']= 10**    12.372546
# red_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=      1.560796
# A_shot_noise_red= 10**(   -5.54813)
# # red_hod_pdict['alpha_s_HOD'] = 1.1849386
# # red_hod_pdict['sigma_log10M_HOD'] =  0.75082121
# # red_hod_pdict['M1_prime_HOD']=  10**13.201484
# # red_hod_pdict['M_min_HOD']= 10**13.234943
# # A_shot_noise_red= 27.954418
# # red_hod_pdict['x_out_truncated_nfw_profile_satellite_galaxies']=  1.2956410
