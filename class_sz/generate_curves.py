import numpy as np
from classy_sz import Class
import time
import sys
sys.path.append('/moto/home/akk2175/ILC_params/')
from HOD import *
from ilc_params import *
from ilc_functions import *

"""
Choose unWISE color, frequencies, and CIB dictionary
"""

color = "blue"
nu_list = 93,100,143,145,217,220,225,280,353,545
nu_list_str = '93,100,143,145,217,220,225,280,353,545'
CIB_pdict = websky_cib_pdict #fiona_alina  #

print(CIB_pdict)
#save?
path_save = "/moto/home/akk2175/class_sz_curves/results/ilc-moto-3-3_websky/"
save_to_file = "yes"

"""
Dictionaries
"""

if color=="blue":
    hod_pdict = blue_pdict
    s = s_blue
if color=="green":
    hod_pdict = green_pdict
    s = s_green
if color=="red":
    hod_pdict = red_pdict
    s = s_red

print(common_settings['dell'])
print(hod_pdict['galaxy_sample_id'])
print(hod_pdict)

#get the CIB flux cut values
cib_flux_list = make_flux_cut_list(cib_flux, nu_list)

"""
Primary lensed CMB
"""

M = Class()
M.set(common_settings_cmb)
M.set(p18_pdict)
M.compute()
cl_tot = M.raw_cl(10000)
cl_lensed = M.lensed_cl(10000)
ell = cl_tot['ell']
M.struct_cleanup()  # clean output
M.empty()           # clean input

factor = TCMB_uK**2 * ell*(ell+1.)/2./np.pi #class output is unitless, so multpiply by TCMB_uk ^2


"""
Compute components with class_sz
"""

N = Class()
N.set(hod_pdict)
N.set(p18_pdict)
N.set(websky_tsz_pdict)
N.set(CIB_pdict)
N.set(common_settings)

N.set({
        'output':'cib_cib_1h,cib_cib_2h,tSZ_cib_1h,tSZ_cib_2h,tSZ_1h,tSZ_2h,gal_cib_1h,gal_cib_2h,tSZ_gal_1h, tSZ_gal_2h,gal_gal_1h, gal_gal_2h, gal_lensmag_1h,gal_lensmag_2h, lensmag_lensmag_1h, lensmag_lensmag_2h, cib_lensmag_1h,cib_lensmag_2h, tSZ_lensmag_1h,tSZ_lensmag_2h',
         'class_sz_verbose':0,
        })
N.set({
        'cib_frequency_list_num' : len(nu_list),
        'cib_frequency_list_in_GHz' : nu_list_str,

        'cib_Snu_cutoff_list [mJy]': str(list(cib_flux_list))[1:-1],
        'has_cib_flux_cut': 1
      })
N.compute()

cl_cib_cib = N.cl_cib_cib()
cl_tsz_cib = N.cl_tSZ_cib()
cl_sz = N.cl_sz()
cl_cib_g = N.cl_gal_cib()
cl_tsz_g = N.cl_yg()
cl_gg = N.cl_gg()
cl_cib_mu = N.cl_cib_m()
cl_tsz_mu = N.cl_ym()
cl_gm = N.cl_gm()
cl_mm = N.cl_mm()

print(cl_cib_mu['93'])

cl_gg_ell = np.asarray(cl_gg['ell'])
cl_gg_1h = np.asarray(cl_gg['1h'])
cl_gg_2h = np.asarray(cl_gg['2h'])

cl_to_dl_gg = cl_gg_ell*(cl_gg_ell+1)/2/np.pi

yy = np.asarray(cl_sz['1h'])+np.asarray(cl_sz['2h'])


"""
Save
"""
#CMB
if save_to_file == "yes":
    np.savetxt(path_save+"ell_dl_CMB_lensed.txt", (ell, factor*cl_lensed['tt']))
    np.savetxt(path_save+"ell_dl_CMB.txt", (ell, factor*cl_tot['tt']))
    #np.savetxt(path_save+"ell_dl_gg_"+color+".txt", (cl_gg_ell, cl_gg_1h+cl_gg_2h))

ell_cib = np.asarray(cl_cib_cib[str(nu_list[0])+'x'+str(nu_list[0])]['ell'])
cls_to_dls = ell_cib*(ell_cib+1.)/2./np.pi
for (i,nu) in enumerate(nu_list):
    ## tSZ
    #plt.plot(cls_tSZ[i]['ell'],cls_tSZ[i]['1h']+cls_tSZ[i]['2h'], color="blue", label=r'tSZ')
    tsz_tsz = yy*abs(tSZ_spectral_funct_at_nu_in_GHz(nu)**2)

    ## CIB
    cls_cib_1h = np.asarray(cl_cib_cib[str(nu_list[i])+'x'+str(nu_list[i])]['1h'])
    cls_cib_2h = np.asarray(cl_cib_cib[str(nu_list[i])+'x'+str(nu_list[i])]['2h'])
    CIB_uK = (cls_cib_1h+cls_cib_2h) / Jysr_to_uKcmb(nu) /Jysr_to_uKcmb(nu)
    ##Save the curves
    if save_to_file == "yes":
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_CIBxCIB.txt", (cl_cib_cib['217x217']['ell'], CIB_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_tSZxtSZ.txt", (cl_sz['ell'],tsz_tsz))


ell_cib = np.asarray(cl_cib_cib[str(nu_list[0])+'x'+str(nu_list[0])]['ell'])
cls_to_dls = ell_cib*(ell_cib+1.)/2./np.pi
for (i,nu) in enumerate(nu_list):
    # tSZ x g
    yg = np.asarray(cl_tsz_g['1h'])+np.asarray(cl_tsz_g['2h'])
    ym = np.asarray(cl_tsz_mu['1h'])+np.asarray(cl_tsz_mu['2h'])
    tSZg_uK = yg*tSZ_spectral_funct_at_nu_in_GHz(nu)
    tSZm_uK = ym*tSZ_spectral_funct_at_nu_in_GHz(nu)

    tSZg_uK = yg*tSZ_spectral_funct_at_nu_in_GHz(nu)

    #CIB x g
    cl_cib_g_1h = np.asarray(cl_cib_g[str(nu_list[i])]['1h'])
    cl_cib_g_2h = np.asarray(cl_cib_g[str(nu_list[i])]['2h'])
    cl_cib_mu_1h = np.asarray(cl_cib_mu[str(nu_list[i])]['1h'])
    cl_cib_mu_2h = np.asarray(cl_cib_mu[str(nu_list[i])]['2h'])
    CIBg_uK = (cl_cib_g_1h + cl_cib_g_2h)/Jysr_to_uKcmb(nu)
    CIBm_uK = (cl_cib_mu_1h + cl_cib_mu_2h)/Jysr_to_uKcmb(nu)



    #tSZ x CIB
    ell_tsz_cib = np.asarray(cl_tsz_cib[str(nu_list[0])]['ell'])
    cl_tsz_cib_1h = np.asarray(cl_tsz_cib[str(nu)]['1h'])
    cl_tsz_cib_2h = np.asarray(cl_tsz_cib[str(nu)]['2h'])
    CIBtSZ_uK = (cl_tsz_cib_1h+cl_tsz_cib_2h)*tSZ_spectral_funct_at_nu_in_GHz(nu)/ Jysr_to_uKcmb(nu)
    ##Save the curves
    if save_to_file == "yes":
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_tSZxCIB.txt", (ell_tsz_cib, CIBtSZ_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_tSZxg_"+color+".txt", (cl_tsz_g['ell'], tSZg_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_CIBxg_"+color+".txt", (cl_cib_g[str(nu_list[1])]['ell'], CIBg_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_tSZxg_wLensmag_"+color+".txt", (cl_tsz_g['ell'], tSZg_uK+(5*s-2)*tSZm_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_CIBxg_wLensmag_"+color+".txt", (cl_cib_g[str(nu_list[1])]['ell'], CIBg_uK +(5*s-2)*CIBm_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_tSZxm_"+color+".txt", (cl_tsz_g['ell'], tSZm_uK))
        np.savetxt(path_save+'ell_dl_'+str(nu)+'x'+str(nu)+"_GHz_CIBxm_"+color+".txt", (cl_cib_g[str(nu_list[1])]['ell'],CIBm_uK))

#radio x radio

ell_radio = np.asarray(cl_sz['ell'])
for (i,nu1) in enumerate(nu_list):
    for (j,nu2) in enumerate(nu_list):
        cl_radio = model_radio(nu1,nu2, nu0_radio_ghz, beta_radio, A_s, ell0, ell_radio) #for the same frequency nu1 x nu2
        if save_to_file == "yes":
            np.savetxt(path_save+'ell_dl_'+str(nu1)+"x"+str(nu2)+"_GHz_radioxradio.txt", ( ell_radio, cl_radio ) )



#CIB x CIB
for (i,nu1) in enumerate(nu_list):
    for (j,nu2) in enumerate(nu_list):
        if nu1!=nu2:
            nu1_nu2 = str(nu1)+"x"+str(nu2)
            nu2_nu1 = str(nu2)+"x"+str(nu1)
            if nu1_nu2 in cl_cib_cib:
                print(nu1_nu2)
                cl_cib_1h = np.asarray(cl_cib_cib[nu1_nu2]['1h'])
                cl_cib_2h = np.asarray(cl_cib_cib[nu1_nu2]['2h'])
                cl_cib_uK = (cl_cib_1h+cl_cib_2h) / Jysr_to_uKcmb(nu1) /Jysr_to_uKcmb(nu2)
                if save_to_file == "yes":
                    np.savetxt(path_save+'ell_dl_'+str(nu1)+"x"+str(nu2)+"_GHz_CIBxCIB.txt", ( cl_cib_cib['217x217']['ell'], cl_cib_uK ) )
                    np.savetxt(path_save+'ell_dl_'+str(nu2)+"x"+str(nu1)+"_GHz_CIBxCIB.txt", ( cl_cib_cib['217x217']['ell'], cl_cib_uK ) )


#tsz x tsz
for (i,nu1) in enumerate(nu_list):
    for (j,nu2) in enumerate(nu_list):
        if nu1!=nu2:
            cl_tsz_cross = ( np.asarray(cl_sz['1h'])+np.asarray(cl_sz['2h']) )*tSZ_spectral_funct_at_nu_in_GHz(nu1)*tSZ_spectral_funct_at_nu_in_GHz(nu2)
            if save_to_file == "yes":
                np.savetxt(path_save+'ell_dl_'+str(nu1)+"x"+str(nu2)+"_GHz_tSZxtSZ.txt", ( cl_sz['ell'], cl_tsz_cross ) )
                np.savetxt(path_save+'ell_dl_'+str(nu2)+"x"+str(nu1)+"_GHz_tSZxtSZ.txt", ( cl_sz['ell'], cl_tsz_cross ) )


# tsz x CIB

for (i,nu_cib) in enumerate(nu_list):
    ell_tsz_cib = np.asarray(cl_tsz_cib[str(nu_cib)]['ell'])
    cl_tsz_cib_1h = np.asarray(cl_tsz_cib[str(nu_cib)]['1h'])
    cl_tsz_cib_2h = np.asarray(cl_tsz_cib[str(nu_cib)]['2h'])
    CIBy_uK = (cl_tsz_cib_1h+cl_tsz_cib_2h)/ Jysr_to_uKcmb(nu_cib)

    for (j,nu_sz) in enumerate(nu_list):
        if nu_sz!=nu_cib:
            CIBtSZ_uK = CIBy_uK*tSZ_spectral_funct_at_nu_in_GHz(nu_sz)
            if save_to_file == "yes":
                np.savetxt(path_save+"ell_dl_"+str(nu_sz)+"x"+str(nu_cib)+"_GHz_tSZxCIB.txt", ( ell_tsz_cib, (CIBtSZ_uK)) )
                np.savetxt(path_save+"ell_dl_"+str(nu_cib)+"x"+str(nu_sz)+"_GHz_CIBxtSZ.txt", ( ell_tsz_cib, (CIBtSZ_uK)) )

for (i,nu_cib) in enumerate(nu_list):
    ell_tsz_cib = np.asarray(cl_tsz_cib[str(nu_cib)]['ell'])
    cl_tsz_cib_1h = np.asarray(cl_tsz_cib[str(nu_cib)]['1h'])
    cl_tsz_cib_2h = np.asarray(cl_tsz_cib[str(nu_cib)]['2h'])
    CIBy_uK = (cl_tsz_cib_1h+cl_tsz_cib_2h)/ Jysr_to_uKcmb(nu_cib)

    for (j,nu_sz) in enumerate(nu_list):
        if nu_sz!=nu_cib:
            CIBtSZ_uK = CIBy_uK*tSZ_spectral_funct_at_nu_in_GHz(nu_sz)
            if save_to_file == "yes":
                np.savetxt(path_save+"ell_dl_"+str(nu_sz)+"x"+str(nu_cib)+"_GHz_tSZxCIB.txt", ( ell_tsz_cib, (CIBtSZ_uK)) )
print("Done")
print(path_save)
