from classy_sz import Class
import numpy as np
from scipy.interpolate import interp1d
import sys
sys.path.append('/Users/aleksandra/Desktop/Research/ILC_params/')

from HOD import *
from ilc_params import *


def chi2(model,data,cov):
    diff = data-model
    inv = np.linalg.inv(cov)
    chi2 = np.dot(diff, np.dot(inv, diff))
    return chi2

def full_chi2(model,data,cov):
    log_det = np.linalg.slogdet(cov)[1]
    const = np.log(2 * np.pi) * (-len(data) / 2) + log_det * (-1 / 2)
    inv = np.linalg.inv(cov)
    diff = data-model
    chi2 = np.dot(diff, inv.dot(diff))
    print("chi2 w sigmas", chi2)
    return -0.5 * chi2 + const

def binning(ell_class, dl_class, ell_data, bpwf, Nellbins=55, conv2cl=True):
    """
    Interpolate class dl's, convert to cl's, and bin, according to alex binning scheme
    in log
    bwf- bandpower window function
    """
    #interpolate and to cl's (Alex data is in cl's)
    #print(ell_class)
    dl_class = np.log(dl_class)
    f_kg = interp1d(ell_class, dl_class)
    new_ell = np.arange(2, ell_data[Nellbins], 1) # up to 1051.5
    inter_dls = np.asarray(f_kg(new_ell))
    inter_dls = np.exp(inter_dls)
    if conv2cl==True:
        inter_cls = inter_dls*(2.0*np.pi)/(new_ell)/(new_ell+1.0)
    #binning / bandpower WF from Alex
    clbinned = np.zeros(Nellbins)
    for i in range (Nellbins):
        wi = bpwf[i]
        #print("wi shape:", wi.shape)
        #print("interp dls :", inter_dls.shape)
        # wi starts from ell=2 according to Alex, email 1-9-22; could add ell=0,1, but would contribute nothing to the sum
        ci_binned = np.sum(wi[2:len(inter_cls)+2]*inter_cls)
        clbinned[i]=ci_binned
    #print("clbinned:", clbinned)
    return ell_data, clbinned, inter_cls



def fiducial_gg_full(params_dict,s, shot_noise, transfer_funct, ell_data, cl_data, cov_data, bpwf_gg, pixwind_bin,ng, Npoints):
    """
    params_dict,
    s, shot_noise, transfer_funct,
    ell_data,
    bpwf_gg - bandpower window function ot bin the theory cl's
    pixwind_bin -- binned pixel window function ("pixel_window_bin_nside2048.txt")
    ng - in steredians
    Npoints
    """
    M = Class()
    M.set(params_dict)
    M.set(p18_pdict)
    M.set(common_settings)
    M.set(websky_tsz_pdict)
    M.set({
        "ell_max":6400,
        'dell':10,
        'output':'gal_gal_1h, gal_gal_2h, gal_lensmag_1h,gal_lensmag_2h, lensmag_lensmag_1h, lensmag_lensmag_2h',
        })
    M.compute()

    # Cl_gxg
    theory = M.cl_gg()
    cl_ell_theory = theory['ell']
    dl_1h_theory = theory['1h']
    dl_2h_theory = theory['2h']
    dl_gg_theory_1h = np.asarray(list(dl_1h_theory))
    dl_gg_theory_2h =  np.asarray(list(dl_2h_theory))
    #print("gg 1h", dl_gg_theory_1h)
    #print("gg theory", dl_gg_theory_2h+dl_gg_theory_1h)
    #print(dl_gg_theory_1h)
    #print(cl_ell_theory)
    ell_theory = np.asarray(list(cl_ell_theory))
    ell_gg_binned, cl_gg_binned_1h, k = binning(ell_theory, dl_gg_theory_1h, ell_data, bpwf_gg, Nellbins=Npoints)
    ell_gg_binned, cl_gg_binned_2h, k = binning(ell_theory, dl_gg_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints)
    ell_gg_binned, cl_gg_binned_tot, inter_cls_gg = binning(ell_theory, dl_gg_theory_1h+dl_gg_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints)
    #print("gg tot", cl_gg_binned_tot)

    # Cl_kxmu
    theory_km =M.cl_gm()# self.theory.get_Cl_kxmu()
    cl_ell_theory_gm = theory_km['ell']
    dl_1h_theory_gm = theory_km['1h']
    dl_2h_theory_gm = theory_km['2h']
    dl_gm_theory_1h = np.asarray(list(dl_1h_theory_gm))
    dl_gm_theory_2h = np.asarray(list(dl_2h_theory_gm))
    ell_theory = np.asarray(list(cl_ell_theory_gm))
    #print("gMu 1h", dl_gm_theory_1h)
    #print("gMu 2h", dl_gm_theory_2h)
    ell_gm_binned, cl_gm_binned_1h,k = binning(ell_theory, dl_gm_theory_1h, ell_data, bpwf_gg, Nellbins=Npoints)
    ell_gm_binned, cl_gm_binned_2h, k = binning(ell_theory, dl_gm_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints,)
    ell_gm_binned, cl_gm_binned_tot, inter_cls_gm = binning(ell_theory, dl_gm_theory_1h+dl_gm_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints,)

    # Cl_muxmu
    theory_mm =M.cl_mm()# self.theory.get_Cl_kxmu()
    cl_ell_theory_mm = theory_mm['ell']
    dl_1h_theory_mm = theory_mm['1h']
    dl_2h_theory_mm = theory_mm['2h']
    dl_mm_theory_1h = np.asarray(list(dl_1h_theory_mm))
    dl_mm_theory_2h = np.asarray(list(dl_2h_theory_mm))
    ell_theory = np.asarray(list(cl_ell_theory_mm))
    #print("Mu 1h", dl_mm_theory_1h)
    #print("Mu 2h", dl_mm_theory_2h)
    ell_mm_binned, cl_mm_binned_1h, k = binning(ell_theory, dl_mm_theory_1h, ell_data, bpwf_gg, Nellbins=Npoints, )
    ell_mm_binned, cl_mm_binned_2h, k = binning(ell_theory, dl_mm_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints, )
    ell_mm_binned, cl_mm_binned_tot, inter_cls_mm = binning(ell_theory, dl_mm_theory_1h+dl_mm_theory_2h, ell_data, bpwf_gg, Nellbins=Npoints, )

    #print("ell_gg: ", ell_gg_binned)
    #print("cl_gg: ", cl_gg_binned_tot)
    #print("cl_gm: ", cl_gm_binned_tot)
    #print("cl_mm: ", cl_mm_binned_tot)
    # km multiplied by 5s-2
    print(" A_shot_noise:",  shot_noise)

    #cl_gg_theory +  A_SN* p_ell^-2 + 1/ng * (1- p_ell^-2)
    p_ell = pixwind_bin[:Npoints]
    #print(p_ell)
    cl_bin_tot  =  cl_gg_binned_tot + 2*(5*s-2)*(cl_gm_binned_tot)  + (5*s-2)*(5*s-2)*(cl_mm_binned_tot) + shot_noise*p_ell**(-2) + 1/ng * (1-p_ell**(-2))
    cl_theory_tot = (dl_gg_theory_1h+dl_gg_theory_2h) + 2*(5*s-2)*(dl_gm_theory_1h+dl_gm_theory_2h) + (5*s-2)*(5*s-2)*(dl_mm_theory_1h+dl_mm_theory_2h)
    cl_inter_tot  =  inter_cls_gg + 2*(5*s-2)*(inter_cls_gm)  + (5*s-2)*(5*s-2)*(inter_cls_mm) + shot_noise*1.e-7

    #apply the transfer function
    #trans function
    trans = np.append(transfer_funct, np.ones(Npoints-len(transfer_funct)))
    cl_bin_tot  =  cl_bin_tot * trans
    #print("cl_bin_tot", cl_bin_tot)

    #cut the first data point in gg
    ell_binned, cl_bin_tot = ell_gm_binned[1:Npoints], cl_bin_tot[1:Npoints]
    chi_sq = chi2(cl_bin_tot, cl_data[1:Npoints], cov_data[1:Npoints,1:Npoints])
    chi_sq_full = full_chi2(cl_bin_tot, cl_data[1:Npoints], cov_data[1:Npoints,1:Npoints])
    M.struct_cleanup()
    M.empty()
    return ell_binned, cl_bin_tot, cl_ell_theory, cl_theory_tot, chi_sq, chi_sq_full, cl_gg_binned_1h[1:], cl_gg_binned_2h[1:], cl_gm_binned_1h[1:],cl_gm_binned_2h[1:],  cl_mm_binned_1h[1:],  cl_mm_binned_2h[1:]


def computed_gg_full(ell_class, cl_class, shot_noise, transfer_funct, ell_data, cl_data, cov_data, bpwf_gg, pixwind_bin,ng, Npoints):
    """
    params_dict,
    s, shot_noise, transfer_funct,
    ell_data,
    bpwf_gg - bandpower window function ot bin the theory cl's
    pixwind_bin -- binned pixel window function ("pixel_window_bin_nside2048.txt")
    ng - in steredians
    Npoints
    """

    print(" A_shot_noise:",  shot_noise)

    #cl_gg_theory +  A_SN* p_ell^-2 + 1/ng * (1- p_ell^-2)
    p_ell = pixwind_bin[:Npoints]
    #print(p_ell)

    ell_class_bin, cl_class_bin, k = binning(ell_class, cl_class, ell_data, bpwf_gg, Nellbins=Npoints)

    cl_bin_tot  =  cl_class_bin + shot_noise*p_ell**(-2) + 1/ng * (1-p_ell**(-2))
    # cl_theory_tot = (dl_gg_theory_1h+dl_gg_theory_2h) + 2*(5*s-2)*(dl_gm_theory_1h+dl_gm_theory_2h) + (5*s-2)*(5*s-2)*(dl_mm_theory_1h+dl_mm_theory_2h)
    # cl_inter_tot  =  inter_cls_gg + 2*(5*s-2)*(inter_cls_gm)  + (5*s-2)*(5*s-2)*(inter_cls_mm) + shot_noise*1.e-7

    #apply the transfer function
    #trans function
    trans = np.append(transfer_funct, np.ones(Npoints-len(transfer_funct)))
    cl_bin_tot  =  cl_bin_tot * trans
    #print("cl_bin_tot", cl_bin_tot)

    #cut the first data point in gg
    print("ell_class_bin: ", ell_class_bin)
    print(cl_bin_tot)
    ell_binned, cl_bin_tot = ell_class_bin[1:Npoints], cl_bin_tot[1:Npoints]
    chi_sq = chi2(cl_bin_tot, cl_data[1:Npoints], cov_data[1:Npoints,1:Npoints])
    chi_sq_full = full_chi2(cl_bin_tot, cl_data[1:Npoints], cov_data[1:Npoints,1:Npoints])

    return ell_binned, cl_bin_tot, chi_sq, chi_sq_full
