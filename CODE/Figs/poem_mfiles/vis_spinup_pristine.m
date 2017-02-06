% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% dpath = [datap 'NoPDc_NoAct_TrefO_flev1e4/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_flev1e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev4e4/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_flev4e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev8e4/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_flev8e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_NoWgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_1e4_NoWgt/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4/'];
% fpath = [figp 'NoPDc_NoMetab_TrefO_1e4/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4_HalfC/'];
% fpath = [figp 'NoPDc_NoMetab_TrefO_1e4_HalfC/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_1e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_4e4_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_4e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_8e4_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_8e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e5_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_1e5_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_1e6_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_Cmax_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoMet_TrefO_Cmax_C&Mwgt/'];
% fpath = [figp 'NoPDc_NoMet_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_noC&Mwgt/'];
% fpath = [figp 'NoPDc_NoAct_TrefO_1e6_noC&Mwgt/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab/'];
% fpath = [figp 'NoPDc_TrefO_KHparams_cmax-metab/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS/'];
% fpath = [figp 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% fpath = [figp 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFprefZ/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFprefZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit15/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit15/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc025/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ00/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ00/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ00/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ00/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ01/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ01/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ05/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ05/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_05/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_05/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_075/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_075/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_125/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_125/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_15/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_15/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Tmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Lmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_LTmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_LTmort/'];
% dpath = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% fpath = [figp 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ025_NOnmort/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05_NOnmort/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ075_NOnmort/';
% npath5 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ1_NOnmort/';
% npath6 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ0_NOnmort/';
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MFMZ05_MPMZ01_NOnmort/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MFMZ1_MPMZ05_NOnmort/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MFMZ1_MPMZ01_NOnmort/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_fcrit10_MZ01_MF2xMP_NOnmort/';
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_nmort/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_Tmort/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_Lmort/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_LTmort/';
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro10/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro20/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro30/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro40/';
% npath5 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro50/';
% npath6 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro60/';
% npath7 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro70/';
% npath8 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro80/';
% npath9 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro90/';
% npath1 = 'NoPDc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/';
% npath2 = 'PDc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/';
% npath1 = 'NoPDc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/';
% npath2 = 'PDc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/';
% npath4 = 'NoPDc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/';
% npath5 = 'PDc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/';
% npath6 = 'Dc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/';
% npath7 = 'NoPDc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/';
% npath8 = 'PDc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/';
% npath9 = 'Dc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath6 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath7 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath8 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath9 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath10 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath11 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath12 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath13 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath14 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath15 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath1 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath2 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath3 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath4 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath1 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath2 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath3 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath4 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath5 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';
% npath2 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE15_BP75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE25_BP75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE15_BP75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP50/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_BP75/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP25/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP50/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25_BP75/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP25/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP50/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP75/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP25/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP50/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_BP75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP50/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP75/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05_BP100/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE10_BP100/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE15_BP100/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE20_BP100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE075/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE075_BP100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt125/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt130/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt135/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt140/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-7/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-3/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-1/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e1/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e3/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_3xZ/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_BE05/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ00_MFMZ01_NOnmort_BE05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ01_MFMZ05_NOnmort_BE05/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ01_MFMZ10_NOnmort_BE05/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ00_NOnmort_BE05/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ05_NOnmort_BE05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ10_NOnmort_BE05/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_S00_NOnmort_BE05/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_S05_NOnmort_BE05/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE001/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE09/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE09/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath3 = 'Dc_TrefO_Hartvig_metab_MFeqMP_MZ01_NOnmort_BE05/';
% npath4 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_NOnmort_BE05/';
% npath5 = 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath6 = 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath7 = 'Dc_TrefO_mizer_metab_MFeqMP_MZ01_NOnmort_BE05/';
% npath8 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_NOnmort_BE05/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/';
% npath3 = 'Dc_TrefO_Hartvig_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/';
% npath4 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_NOnmort_BE05_RE01/';
% npath5 = 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/';
% npath6 = 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/';
% npath7 = 'Dc_TrefO_mizer_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/';
% npath8 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_NOnmort_BE05_RE01/';
% npath9 = 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath10 = 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath11 = 'Dc_TrefO_JC_metab_MFeqMP_MZ01_NOnmort_BE05/';
% npath12 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_NOnmort_BE05/';
% npath13 = 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/';
% npath14 = 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/';
% npath15 = 'Dc_TrefO_JC_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/';
% npath16 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_NOnmort_BE05_RE01/';
npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_Hartvignmort_BE05_RE01/';
npath2 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort_BE05_RE01/';
npath3 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_mizernmort_BE05_RE01/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_mizernmort_BE05_RE01/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_Hartvignmort_BE05_RE01/';
npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_Hartvignmort_BE05_RE01/';
npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_JCnmort_BE05_RE01/';
npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_JCnmort_BE05_RE01/';
npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_JCnmort_BE05_RE01/';
npath12 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01/';

npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500/';

dp = {npath0};
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% dp = {npath6;npath7;npath8;npath9;npath10};
% dp = {npath14;npath15;npath16};
% dp = {npath1};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
%     npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20;npath21};

%spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
spots = {'SPac'};

cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

fplot=1;

%%
if (fplot==1)
    for i=1:length(dp)
        %%
        dpath = [datap char(dp(i))];
        fpath = [figp char(dp(i))];
        
        mclev=NaN*ones(length(spots),8);
        Zcon=NaN*ones(length(spots),3);
        
        %%
        for s=1:length(spots)
            %%
            close all
            loc = spots{s};
            sname = 'Spinup_';
            lname = [loc '_'];
            %lname = ['phen_' loc '_'];
            SP = csvread([dpath sname lname 'Sml_p.csv']);
            SF = csvread([dpath sname lname 'Sml_f.csv']);
            SD = csvread([dpath sname lname 'Sml_d.csv']);
            MP = csvread([dpath sname lname 'Med_p.csv']);
            MF = csvread([dpath sname lname 'Med_f.csv']);
            MD = csvread([dpath sname lname 'Med_d.csv']);
            LP = csvread([dpath sname lname 'Lrg_p.csv']);
            LD = csvread([dpath sname lname 'Lrg_d.csv']);
            C = csvread([dpath sname lname 'Cobalt.csv']);
            z(:,1) = C(:,2);
            z(:,2) = C(:,3);
            z(:,3) = C(:,4);
            z=floor(z);
            
            %% Plots over time
            x=1:length(SP);
            y=x/365;
            lstd=length(SP);
            
            %% Mean consumption level
            c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
            mclev(s,:) = nanmean(c);
            
            % Zoop overconsumption
            
            Zcon(s,:) = nansum(z)/lstd;
            
            %% PLOTS
            
            %% Piscivore
            figure(1)
            subplot(4,1,1)
            plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
            plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
            plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title(['Spinup Pelagic Piscivores ' loc])
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            legend('Larvae','Juveniles','Adults')
            stamp(cfile)
            
            subplot(4,1,2)
            plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Larvae')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,3)
            plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Juveniles')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,4)
            plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Adults')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            print('-dpng',[fpath sname lname 'oneloc_pisc_time.png'])
            
            %% Planktivore
            figure(2)
            subplot(3,1,1)
            plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
            plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title(['Spinup Forage Fishes ' loc])
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            legend('Immature','Adults')
            stamp(cfile)
            
            subplot(3,1,2)
            plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Immature')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(3,1,3)
            plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Adults')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            print('-dpng',[fpath sname lname 'oneloc_plan_time.png'])
            
            %% Detritivore
            figure(3)
            subplot(4,1,1)
            plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
            plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
            plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title(['Spinup Demersal Piscivores ' loc])
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            legend('Larvae','Juveniles','Adults')
            stamp(cfile)
            
            subplot(4,1,2)
            plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Larvae')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,3)
            plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Juveniles')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,4)
            plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Adults')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            print('-dpng',[fpath sname lname 'oneloc_detr_time.png'])
            
            %% All biomass in subplots
            %SP
            figure(4)
            subplot(3,3,2)
            plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title({loc; 'SP'})
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            stamp(cfile)
            
            subplot(3,3,5)
            plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MP')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(3,3,8)
            plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('LP')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            %FF
            subplot(3,3,1)
            plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SF')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(3,3,4)
            plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MF')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            %Detritivore
            subplot(3,3,3)
            plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SD')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            stamp(cfile)
            
            subplot(3,3,6)
            plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MD')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(3,3,9)
            plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('LD')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            print('-dpng',[fpath sname lname 'oneloc_all_sizes_sub.png'])
            
            %% All size classes of all
            
            figure(5)
            plot(y,log10(SP(:,1)),'Linewidth',1); hold on;
            plot(y,log10(MP(:,1)),'Linewidth',1); hold on;
            plot(y,log10(LP(:,1)),'Linewidth',1); hold on;
            plot(y,log10(SF(:,1)),'Linewidth',1); hold on;
            plot(y,log10(MF(:,1)),'Linewidth',1); hold on;
            plot(y,log10(SD(:,1)),'Linewidth',1); hold on;
            plot(y,log10(MD(:,1)),'Linewidth',1); hold on;
            plot(y,log10(LD(:,1)),'Linewidth',1); hold on;
            legend('SP','MP','LP','SF','MF','SD','MD','LD')
            legend('location','eastoutside')
            xlim([y(1) y(end)])
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            title(['Spinup ' loc])
            stamp(cfile)
            print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
            
            %% Final mean biomass size spectrum
            t=1:length(SP);
            lyr=t((end-365+1):end);
            
            SP_sum=sum(SP(lyr,1));
            SF_sum=sum(SF(lyr,1));
            SD_sum=sum(SD(lyr,1));
            MP_sum=sum(MP(lyr,1));
            MF_sum=sum(MF(lyr,1));
            MD_sum=sum(MD(lyr,1));
            LP_sum=sum(LP(lyr,1));
            LD_sum=sum(LD(lyr,1));
            
            SP_mean=mean(SP(lyr,1));
            SF_mean=mean(SF(lyr,1));
            SD_mean=mean(SD(lyr,1));
            MP_mean=mean(MP(lyr,1));
            MF_mean=mean(MF(lyr,1));
            MD_mean=mean(MD(lyr,1));
            LP_mean=mean(LP(lyr,1));
            LD_mean=mean(LD(lyr,1));
            
            P_sum=[SP_sum;MP_sum;LP_sum];
            F_sum=[SF_sum;MF_sum];
            D_sum=[SD_sum;MD_sum;LD_sum];
            P_mean=[SP_mean;MP_mean;LP_mean];
            F_mean=[SF_mean;MF_mean];
            D_mean=[SD_mean;MD_mean;LD_mean];
            
            Pwgt = [0.0025; 2.5298; 2.5298e3];
            Fwgt = [0.0025; 2.5298];
            Dwgt = [0.0025; 2.5298; 2.5298e3];
            
            figure(6)
            subplot(2,3,1)
            bar(log10(P_sum),'k')
            xlim([0 4])
            title('Pel Pisc')
            ylabel('log10 Total Biomass (g m^-^2)')
            subplot(2,3,4)
            bar(log10(P_mean),'k')
            xlim([0 4])
            ylabel('log10 Mean Biomass (g m^-^2)')
            stamp(cfile)
            
            subplot(2,3,2)
            bar(log10(F_sum),'b')
            xlim([0 3])
            title({loc; 'Forage Fishes'})
            xlabel('Stage')
            subplot(2,3,5)
            bar(log10(F_mean),'b')
            xlim([0 3])
            xlabel('Stage')
            
            subplot(2,3,3)
            bar(log10(D_sum),'r')
            xlim([0 4])
            title('Dem Pisc')
            subplot(2,3,6)
            bar(log10(D_mean),'r')
            xlim([0 4])
            print('-dpng',[fpath sname lname 'oneloc_all_biomass_spec.png'])
            
            %% Reproduction
            rep(:,1)=MF(:,1).*MF(:,18);
            rep(:,2)=LD(:,1).*LD(:,18);
            rep(:,3)=LP(:,1).*LP(:,18);
            
            figure(7)
            subplot(4,1,1)
            plot(y,log10(rep),'Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title(['Spinup Reproduction ' loc])
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            legend('F','D','P')
            stamp(cfile)
            
            subplot(4,1,2)
            plot(y,log10(rep(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Forage Fishes')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,3)
            plot(y,log10(rep(:,2)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Demersal Piscivores')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            
            subplot(4,1,4)
            plot(y,log10(rep(:,3)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('Pelagic Piscivores')
            xlabel('Time (y)')
            ylabel('log10 Biomass (g m^-^2)')
            print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
            
            %% Maturation
            m(:,1)=MF(:,19);
            m(:,2)=MD(:,19);
            m(:,3)=MP(:,19);
            m(:,4)=LD(:,19);
            m(:,5)=LP(:,19);
            
            figure(8)
            subplot(3,2,1)
            plot(y,log10(m(:,1)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title([loc ' log10 Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
            ylabel('Forage Fishes')
            stamp(cfile)
            
            subplot(3,2,3)
            plot(y,log10(m(:,2)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('M')
            ylabel('Demersal Piscivores')
            
            subplot(3,2,5)
            plot(y,log10(m(:,3)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('M')
            ylabel('Pelagic Piscivores')
            xlabel('Time (y)')
            
            subplot(3,2,4)
            plot(y,log10(m(:,4)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('L')
            ylabel('Demersal Piscivores')
            
            subplot(3,2,6)
            plot(y,log10(m(:,5)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('L')
            xlabel('Time (y)')
            ylabel('Pelagic Piscivores')
            print('-dpng',[fpath sname lname 'oneloc_matur_time.png'])
            
            %% Predation mortality
            %SP
            figure(9)
            subplot(2,3,2)
            plot(y,log10(SP(:,17)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title({loc; 'SP'})
            stamp(cfile)
            
            subplot(2,3,5)
            plot(y,log10(MP(:,17)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MP')
            xlabel('Time (y)')
            
            %FF
            subplot(2,3,1)
            plot(y,log10(SF(:,17)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SF')
            ylabel('log10 Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
            
            subplot(2,3,4)
            plot(y,log10(MF(:,17)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MF')
            xlabel('Time (y)')
            %ylabel('Biomass eaten by predators (g m^-^2)')
            
            %Detritivore
            subplot(2,3,3)
            plot(y,log10(SD(:,17)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SD')
            
            subplot(2,3,6)
            plot(y,log10(MD(:,17)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MD')
            xlabel('Time (y)')
            
            print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub.png'])
            
            %% All consumption in subplots
            %SP
            figure(10)
            subplot(3,3,2)
            plot(y,(SP(:,14)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title({loc; 'SP'})
            stamp(cfile)
            
            subplot(3,3,5)
            plot(y,(MP(:,14)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MP')
            
            subplot(3,3,8)
            plot(y,(LP(:,14)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('LP')
            xlabel('Time (y)')
            
            %FF
            subplot(3,3,1)
            plot(y,(SF(:,14)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SF')
            
            subplot(3,3,4)
            plot(y,(MF(:,14)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MF')
            xlabel('Time (y)')
            ylabel('Biomass Consumed (g g^-^1 m^-^2)')
            
            %Detritivore
            subplot(3,3,3)
            plot(y,(SD(:,14)),'b','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('SD')
            
            subplot(3,3,6)
            plot(y,(MD(:,14)),'r','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('MD')
            
            subplot(3,3,9)
            plot(y,(LD(:,14)),'k','Linewidth',1); hold on;
            xlim([y(1) y(end)])
            title('LD')
            xlabel('Time (y)')
            print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub.png'])
            
            
            %% Recruitment
            yr=x;
            SFL=m(yr,3);
            SDL=m(yr,4);
            SPL=m(yr,5);
            FA=MF(yr,1);
            DA=LD(yr,1);
            PA=LP(yr,1);
            
            st=1:365:length(yr);
            en=365:365:length(yr);
            SPy = NaN*ones(100,1);
            SFy = SPy;
            SDy = SPy;
            PAy = SPy;
            FAy = SPy;
            DAy = SPy;
            for n=1:100
                SPy(n) = nansum(SPL(st(n):en(n)));
                SFy(n) = nansum(SFL(st(n):en(n)));
                SDy(n) = nansum(SDL(st(n):en(n)));
                PAy(n) = nansum(PA(st(n):en(n)));
                FAy(n) = nansum(FA(st(n):en(n)));
                DAy(n) = nansum(DA(st(n):en(n)));
            end
            
            %
            figure(11)
            subplot(3,1,3)
            plot(1:100,log10(SPy),'k','Linewidth',2); hold on;
            xlim([1 100])
            ylabel('log10 Recruits (g m^-^2)')
            title('Pelagic piscivores')
            stamp(cfile)
            
            subplot(3,1,1)
            plot(1:100,log10(SFy),'b','Linewidth',2); hold on;
            xlim([1 100])
            ylabel('log10 Recruits (g m^-^2)')
            title({loc; 'Forage fishes'})
            
            subplot(3,1,2)
            plot(1:100,log10(SDy),'r','Linewidth',2); hold on;
            xlim([1 100])
            ylabel('log10 Recruits (g m^-^2)')
            title('Demersal piscivores')
            print('-dpng',[fpath sname lname 'oneloc_recruitment.png'])
            
        end
        
        if (length(lname) > 5)
            save([dpath sname lname(1:5) 'consump.mat'],'mclev','Zcon');
            csvwrite([dpath sname lname(1:5) 'clevel.csv'],mclev);
            csvwrite([dpath sname lname(1:5) 'Zconsump.csv'],Zcon);
        else
            save([dpath sname 'consump.mat'],'mclev','Zcon');
            csvwrite([dpath sname 'clevel.csv'],mclev);
            csvwrite([dpath sname 'Zconsump.csv'],Zcon);
        end
        
    end
else
    for i=1:length(dp)
        
        dpath = [datap char(dp(i))];
        fpath = [figp char(dp(i))];
        
        mclev=NaN*ones(length(spots),8);
        Zcon=NaN*ones(length(spots),3);
        
        %%
        for s=1:length(spots)
            %%
            close all
            loc = spots{s};
            sname = 'Spinup_';
            lname = [loc '_'];
            %lname = ['phen_' loc '_'];
            SP = csvread([dpath sname lname 'Sml_p.csv']);
            SF = csvread([dpath sname lname 'Sml_f.csv']);
            SD = csvread([dpath sname lname 'Sml_d.csv']);
            MP = csvread([dpath sname lname 'Med_p.csv']);
            MF = csvread([dpath sname lname 'Med_f.csv']);
            MD = csvread([dpath sname lname 'Med_d.csv']);
            LP = csvread([dpath sname lname 'Lrg_p.csv']);
            LD = csvread([dpath sname lname 'Lrg_d.csv']);
            C = csvread([dpath sname lname 'Cobalt.csv']);
            z(:,1) = C(:,2);
            z(:,2) = C(:,3);
            z(:,3) = C(:,4);
            z=floor(z);
            
            %% Mean consumption level
            c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
            mclev(s,:) = nanmean(c);
            
            % Zoop overconsumption
            lstd=length(C);
            Zcon(s,:) = nansum(z)/lstd;
            
        end
        
        if (length(lname) > 5)
            save([dpath sname lname(1:5) 'consump.mat'],'mclev','Zcon');
            csvwrite([dpath sname lname(1:5) 'clevel.csv'],mclev);
            csvwrite([dpath sname lname(1:5) 'Zconsump.csv'],Zcon);
        else
            save([dpath sname 'consump.mat'],'mclev','Zcon');
            csvwrite([dpath sname 'clevel.csv'],mclev);
            csvwrite([dpath sname 'Zconsump.csv'],Zcon);
        end
        
    end
end

%% DEBUG ---------------------------------------------------
% Find NaNs
mos=8894:8897;
tend=length(mos);
bio=NaN*ones(9,tend);
con=NaN*ones(8,tend);
die=NaN*ones(8,tend);
gamma=NaN*ones(8,tend);
rec=NaN*ones(8,tend);
nu=NaN*ones(8,tend);

for k=1:tend
    n=mos(k);
    %bio
    id = isnan(SP(n,1));
    bio(1,k) = sum(id);
    id = isnan(SF(n,1));
    bio(2,k) = sum(id);
    id = isnan(SD(n,1));
    bio(3,k) = sum(id);
    id = isnan(MP(n,1));
    bio(4,k) = sum(id);
    id = isnan(MF(n,1));
    bio(5,k) = sum(id);
    id = isnan(MD(n,1));
    bio(6,k) = sum(id);
    id = isnan(LP(n,1));
    bio(7,k) = sum(id);
    id = isnan(LD(n,1));
    bio(8,k) = sum(id);
    id = isnan(C(n,1));
    bio(9,k) = sum(id);
    %con
    id = isnan(SP(n,14));
    con(1,k) = sum(id);
    id = isnan(SF(n,14));
    con(2,k) = sum(id);
    id = isnan(SD(n,14));
    con(3,k) = sum(id);
    id = isnan(MP(n,14));
    con(4,k) = sum(id);
    id = isnan(MF(n,14));
    con(5,k) = sum(id);
    id = isnan(MD(n,14));
    con(6,k) = sum(id);
    id = isnan(LP(n,14));
    con(7,k) = sum(id);
    id = isnan(LD(n,14));
    con(8,k) = sum(id);
    %die
    id = isnan(SP(n,17));
    die(1,k) = sum(id);
    id = isnan(SF(n,17));
    die(2,k) = sum(id);
    id = isnan(SD(n,17));
    die(3,k) = sum(id);
    id = isnan(MP(n,17));
    die(4,k) = sum(id);
    id = isnan(MF(n,17));
    die(5,k) = sum(id);
    id = isnan(MD(n,17));
    die(6,k) = sum(id);
    id = isnan(LP(n,17));
    die(7,k) = sum(id);
    id = isnan(LD(n,17));
    die(8,k) = sum(id);
    %gamma/rep
    id = isnan(SP(n,16));
    gamma(1,k) = sum(id);
    id = isnan(SF(n,16));
    gamma(2,k) = sum(id);
    id = isnan(SD(n,16));
    gamma(3,k) = sum(id);
    id = isnan(MP(n,16));
    gamma(4,k) = sum(id);
    id = isnan(MF(n,18));
    gamma(5,k) = sum(id);
    id = isnan(MD(n,16));
    gamma(6,k) = sum(id);
    id = isnan(LP(n,18));
    gamma(7,k) = sum(id);
    id = isnan(LD(n,18));
    gamma(8,k) = sum(id);
    %rec
    id = isnan(SP(n,19));
    rec(1,k) = sum(id);
    id = isnan(SF(n,19));
    rec(2,k) = sum(id);
    id = isnan(SD(n,19));
    rec(3,k) = sum(id);
    id = isnan(MP(n,19));
    rec(4,k) = sum(id);
    id = isnan(MF(n,19));
    rec(5,k) = sum(id);
    id = isnan(MD(n,19));
    rec(6,k) = sum(id);
    id = isnan(LP(n,19));
    rec(7,k) = sum(id);
    id = isnan(LD(n,19));
    rec(8,k) = sum(id);
    %nu
    id = isnan(SP(n,15));
    nu(1,k) = sum(id);
    id = isnan(SF(n,15));
    nu(2,k) = sum(id);
    id = isnan(SD(n,15));
    nu(3,k) = sum(id);
    id = isnan(MP(n,15));
    nu(4,k) = sum(id);
    id = isnan(MF(n,15));
    nu(5,k) = sum(id);
    id = isnan(MD(n,15));
    nu(6,k) = sum(id);
    id = isnan(LP(n,15));
    nu(7,k) = sum(id);
    id = isnan(LD(n,15));
    nu(8,k) = sum(id);
    
end
%%
figure
bar(bio')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD','B')
legend('location','northwest')
print('-dpng',[fpath 'nan_biotest.png'])

figure
bar(con')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[fpath 'nan_contest.png'])

figure
bar(die')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[fpath 'nan_dietest.png'])

figure
bar(gamma')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[fpath 'nan_gammatest.png'])

figure
bar(rec')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[fpath 'nan_rectest.png'])

figure
bar(nu')
set(gca,'XTickLabel',mos)
legend('SP','SF','SD','MP','MF','MD','LP','LD')
legend('location','northwest')
print('-dpng',[fpath 'nan_nutest.png'])





