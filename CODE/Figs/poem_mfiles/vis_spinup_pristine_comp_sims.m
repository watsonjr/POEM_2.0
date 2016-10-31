%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% dpath = [datap 'NoPDc_NoAct_TrefO_flev1e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev4e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_flev8e4/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_NoWgt/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4/'];
% dpath = [datap 'NoPDc_NoMetab_TrefO_1e4_HalfC/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_4e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_8e4_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e5_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoMet_TrefO_Cmax_C&Mwgt/'];
% dpath = [datap 'NoPDc_NoAct_TrefO_1e6_noC&Mwgt/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS/'];
% dpath = [datap 'NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MprefLZoverMZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFprefZ/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/'];
% dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/'];
% dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit15/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/'];
% dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann_Tmort/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann_Tmort/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
% dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
% dpath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
% dpath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
% dpath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
% dpath12 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
% dpath13 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
% dpath14 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
% dpath15 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ00/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ00/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ01/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MPMZ05/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_05/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_075/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_125/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_McannInc_15/'];
% dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_NOnmort/'];
% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Tmort/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Lmort/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_LTmort/'];
% npath0 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish025/'];
% npath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05/'];
% npath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10/'];
% npath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20/'];
% npath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30/'];
% npath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40/'];
% npath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish45/'];
% npath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/'];
% npath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish55/'];
% npath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60/'];
% npath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70/'];
% npath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05_halfM/'];
% npath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/'];
% npath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/'];
% npath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/'];
% npath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/'];
% npath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/'];
% npath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/'];
% npath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/'];
% npath0 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ0_NOnmort/'];
% npath10 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath2 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ025_NOnmort/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ05_NOnmort/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ075_NOnmort/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ1_NOnmort/'];
% npath0 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish025/'];
% npath1 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05/'];
% npath2 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish075/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30/'];
% npath6 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40/'];
% npath7 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/'];
% npath8 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60/'];
% npath9 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70/'];
% npath2 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_nmort/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_Tmort/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_Lmort/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_LTmort/'];
% npath2 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/'];
% npath6 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/'];
% npath7 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/'];
% npath8 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/'];
% npath9 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/'];
% npath10 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish80_halfM/'];
% npath11 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish90_halfM/'];
% npath1 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro10/'];
% npath2 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro20/'];
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro30/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro40/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro50/'];
% npath6 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro60/'];
% npath7 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro70/'];
% npath8 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro80/'];
% npath9 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_repro90/'];
% npath1 = [datap 'NoPDc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/'];
% npath2 = [datap 'Dc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/'];
% npath3 = [datap 'PDc_TrefO_KHparams_all_resp025_MFeqMP_MZ01_NOnmort/'];
% npath4 = [datap 'NoPDc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/'];
% npath5 = [datap 'Dc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/'];
% npath6 = [datap 'PDc_TrefO_KHparams_all_resp05_MFeqMP_MZ01_NOnmort/'];
% npath7 = [datap 'NoPDc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/'];
% npath8 = [datap 'Dc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/'];
% npath9 = [datap 'PDc_TrefO_KHparams_all_resp075_MFeqMP_MZ01_NOnmort/'];
% npath10 = [datap 'NoPDc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/'];
% npath11 = [datap 'Dc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/'];
% npath12 = [datap 'PDc_TrefO_KHparams_all_MFeqMP_MZ01_NOnmort/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath8 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath9 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath10 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath11 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath12 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath13 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath14 = [datap 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath15 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/'];
% npath16 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/'];
% npath17 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath18 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/'];
% npath19 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath20 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/'];
% npath21 = [datap 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish05/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish10/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish20/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish30/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish40/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish05/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish10/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish20/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish30/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish40/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish05/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish10/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish20/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish30/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish40/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish05/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish10/'];
% npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish20/'];
% npath19 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish30/'];
% npath20 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish40/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ00_NOnmort_BE05/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_BE05/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ05_NOnmort_BE05/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ10_NOnmort_BE05/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_S00_NOnmort_BE05/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_S05_NOnmort_BE05/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ00_MFMZ01_NOnmort_BE05/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ01_MFMZ05_NOnmort_BE05/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_fcrit10_MPMZ01_MFMZ10_NOnmort_BE05/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE001/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE02/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE03/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE04/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE05/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE06/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE07/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE08/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE09/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE02/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE03/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE04/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE05/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE06/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE07/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE08/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE09/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath3 = [datap 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath6 = [datap 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_metab_MFeqMP_MZ01_NOnmort_BE05/'];
% npath8 = [datap 'Dc_TrefO_mizer_metab_MFeqMP_MZ01_NOnmort_BE05/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_NOnmort_BE05/'];
% npath10 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_NOnmort_BE05/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath3 = [datap 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath6 = [datap 'Dc_TrefO_mizer_cmax-metab40_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath8 = [datap 'Dc_TrefO_mizer_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath10 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath11 = [datap 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
% npath12 = [datap 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
% npath13 = [datap 'Dc_TrefO_JC_metab_MFeqMP_MZ01_NOnmort_BE05/'];
% npath14 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_NOnmort_BE05/'];
% npath15 = [datap 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath16 = [datap 'Dc_TrefO_JC_cmax-metab25_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath17 = [datap 'Dc_TrefO_JC_metab_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath18 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_NOnmort_BE05_RE01/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish005/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish01/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish02/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish03/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish04/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish05/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish06/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_JCnmort_BE05_RE01/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_Hartvignmort_BE05_RE01/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_mizernmort_BE05_RE01/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_JCnmort_BE05_RE01/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_Hartvignmort_BE05_RE01/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_mizernmort_BE05_RE01/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_BE05/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_JCnmort_BE05_RE01/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_Hartvignmort_BE05_RE01/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01/'];
% npath12 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort_BE05_RE01/'];
% npath14 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish01/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish02/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish03/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish04/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish05/'];
% npath0 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish01/'];
% npath2 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish02/'];
% npath3 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish03/'];
% npath4 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish04/'];
% npath5 = [datap 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish05/'];
% npath0 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish01/'];
% npath2 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish02/'];
% npath3 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish03/'];
% npath4 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish04/'];
% npath5 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish05/'];
npath0 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001/'];
npath1 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish01/'];
npath2 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish02/'];
npath3 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish03/'];
npath4 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish04/'];
npath5 = [datap 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish05/'];

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% dp = {dpath1; dpath2; dpath3; dpath4; dpath5; dpath6};
% sims = {'f05','f10','f15','f20','f30','f40'};
% cfile = 'PDc_MFeqMP_fcrit_comp';
% dp = {dpath3;dpath4;dpath2;dpath5;dpath6;dpath7;dpath8};
% sims = {'MP0.50','MP0.75','Meq','MF1.5','MF2.0','MF2.5','MF3.0'};
% cfile = 'PDc_MFbetter_enc_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5;dpath6;dpath7};
% sims = {'con','T','L','T&L','simpQ','compQ1','compQ1'};
% cfile = 'PDc_MFeqMP_mort_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5};
% sims = {'full','redSM','T-redSM','redS','T-redS'};
% cfile = 'PDc_MFeqMP_cannibal_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5;dpath6;dpath7;dpath8};
% sims = {'0.0+N','0.01+N','0.025+N','0.05+N','0.10+N','0.20+N','0.30+N','0.40+N'};
% cfile = 'PDc_MFeqMP_fishing_nmort_comp';
% dp = {dpath1;dpath9;dpath10;dpath11;dpath12;dpath13;dpath14;dpath15};
% sims = {'0.0','0.01','0.025','0.05','0.10','0.20','0.30','0.40'};
% cfile = 'PDc_MFeqMP_fishing_comp';
% dp = {dpath2;dpath3;dpath4;dpath1;dpath5;dpath6;dpath7};
% sims = {'0.0','0.1','0.5','1.0','MP0.0','MP0.1','MP0.5'};
% cfile = 'PDc_MFeqMP_MZ_comp';
% dp = {dpath2;dpath3;dpath1;dpath4;dpath5};
% sims = {'0.5','0.75','1.0','1.25','1.5'};
% cfile = 'PDc_MFeqMP_Mcann_comp';
% dp = {dpath1;dpath2;dpath3;dpath4;dpath5};
% sims = {'none','const','T','L','T&L'};
% cfile = 'PDc_MFeqMP_mort2_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;...
%     npath10;npath11};
% sims = {'0.0MZ01','0.025MZ01','0.05MZ01','0.10MZ01','0.20MZ01','0.30MZ01',...
%     '0.40MZ01','0.45MZ01','0.50MZ01','0.55MZ01','0.60MZ01','0.70MZ01'};
% cfile = 'PDc_MFeqMP_MZ01_fishing_comp';
% dp = {npath0;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% sims = {'0.0','0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70'};
% cfile = 'PDc_MFeqMP_MZ01_halfM_fishing_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5};
% sims = {'0.0','0.1','0.25','0.5','0.75','1.0'};
% cfile = 'Dc_MFeqMP_NOnmort_MZpref_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;...
%     npath9};
% sims = {'0.025','0.05','0.075','0.10','0.20','0.30','0.40','0.50',...
%     '0.60','0.70'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_comp';
% dp = {npath1;npath2;npath3;npath4;npath5};
% sims = {'none','const','T','L','T&L'};
% cfile = 'Dc_MFeqMP_MZ01_mort_comp';
% dp = {npath1;npath2;npath3;npath4;npath5};
% sims = {'0.1','0.2','0.3','0.4','0.5'};
% cfile = 'Dc_MFeqMP_MZ01_fcrit_comp';
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11};
% sims = {'0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_half_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'};
% cfile = 'Dc_MFeqMP_MZ01_reproeff_comp';
% dp = {npath1;npath2;npath3};
% sims = {'noPDc','PDc','Dc'};
% cfile = 'K&Hparams_all_MFeqMP_MZ01_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12};
% sims = {'no25','D25','PD25','no50','D50','PD50','no75','D75','PD75','no100','D100','PD100'};
% cfile = 'K&Hparams_all_MFeqMP_MZ01_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath8;npath9;npath10;npath11;npath12;npath13;npath14};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'PDc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath15;npath16;npath17;npath18;npath19;npath20;npath21};
% sims = {'0.05','0.075','0.1','0.2','0.3','0.4','0.5'};
% cfile = 'NoPDc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit_comp';
% dp = {npath17;npath3;npath10};
% sims = {'noPDc','Dc','PDc'};
% cfile = 'Hartvig_cmax-metab_MFeqMP_MZ01_fcrit10_PDc_comp';
% dp = {npath21;npath7;npath14};
% sims = {'noPDc','Dc','PDc'};
% cfile = 'Hartvig_cmax-metab_MFeqMP_MZ01_fcrit50_PDc_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6};
% sims = {'0.025','0.05','0.075','0.1','0.15','0.2'};
% cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_BentEff_comp';
% dp = {npath1;npath2;npath3;npath4;npath5};
% sims = {'.05','.1','.2','.3','.4'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit30_BE20_fishing_comp';
% dp = {npath6;npath7;npath8;npath9;npath10};
% sims = {'.05','.1','.2','.3','.4'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit30_BE30_fishing_comp';
% dp = {npath11;npath12;npath13;npath14;npath15};
% sims = {'.05','.1','.2','.3','.4'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit40_BE20_fishing_comp';
% dp = {npath16;npath17;npath18;npath19;npath20};
% sims = {'.05','.1','.2','.3','.4'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit40_BE30_fishing_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8};
% sims = {'MZ0','MZ01','MZ05','MZ1','S0','S05','P0F01','P01F05','P01F1'};
% cfile = 'Dc_Hartvig_fcrit10_BE05_Mpref_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'.01','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit30_BE05_RepEff_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit40_BE05_RepEff_comp';
% dp = {npath1;npath10};
% sims = {'85-10C','60-15C'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_fcrit30_BE05_cmaxH_comp';
% dp = {npath1;npath2;npath4;npath5;npath7;npath9};
% sims = {'60-15C-30','85-10C-30','60-15C-40','85-10C-40','metab','all'};
% cfile = 'Dc_Hartvig_MFeqMP_MZ01_BE05_eqs_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'60-15C-30','85-10C-30','40-10C-30','60-15C-40','85-10C-40','40-10C-40',...
%     'Hmetab','Mmetab','Hall','Mall'};
% cfile = 'Dc_Hartvig_mizer_MFeqMP_MZ01_BE05_eqs_comp';
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'60-15C-30','85-10C-30','40-10C-30','60-15C-40','85-10C-40','40-10C-40',...
%     'Hmetab','Mmetab','Hall','Mall'};
% cfile = 'Dc_Hartvig_mizer_MFeqMP_MZ01_BE05_RE01_eqs_comp';
% dp = {npath1;npath2;npath3;npath11;npath4;npath5;npath6;npath12;npath7;npath8;npath13;npath9;npath10;npath14};
% sims = {'60-15C-30','85-10C-30','40-10C-30','25-10C-30','60-15C-40','85-10C-40','40-10C-40',...
% '25-10C-40','Hmetab','Mmetab','Jmetab','Hall','Mall','Jall'};
% cfile = 'Dc_Hartvig_mizer_JC_MFeqMP_MZ01_BE05_eqs_comp';
% dp = {npath1;npath2;npath3;npath15;npath4;npath5;npath6;npath16;npath7;npath8;npath17;npath9;npath10;npath18};
% sims = {'60-15C-30','85-10C-30','40-10C-30','25-10C-30','60-15C-40','85-10C-40','40-10C-40',...
% '25-10C-40','Hmetab','Mmetab','Jmetab','Hall','Mall','Jall'};
% cfile = 'Dc_Hartvig_mizer_JC_MFeqMP_MZ01_BE05_RE01_eqs_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7};
% sims = {'0','.05','.1','.2','.3','.4','.5','.6'};
% cfile = 'Dc_MFeqMP_fcrit40_MZ01_BE05_RE01_fishing_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;...
%     npath10;npath11;npath12;npath13;npath14};
% sims = {'40N','40J','40H','40M','30N','30J','30H','30M','10N','10J','10H','10M',...
%     'Jall','Hall','Mall'};
% cfile = 'Dc_Hartvig_mizer_JC_MFeqMP_MZ01_BE05_RE01_mort_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5};
% sims = {'0','.1','.2','.3','.4','.5'};
% cfile = 'Dc_MFeqMP_fcrit10_MZ01_BE05_RE01_mizernmort_fishing_comp';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;};
% sims = {'0','.1','.2','.3','.4','.5'};
% cfile = 'Dc_mizer_all_MFeqMP_MZ01_BE05_RE01_fishing_comp';
dp = {npath0;npath1;npath2;npath3;npath4;npath5;};
sims = {'0','.1','.2','.3','.4','.5'};
cfile = 'Dc_JC_all_MFeqMP_MZ01_BE05_RE001_fishing_comp';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

%%

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met'};
cols=cols';

ndp = length(dp);

load('cmap_ppt_angles.mat')

%%
for i=1:length(spots)
    close all
    
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:ndp
        
        dpath = char(dp(s));
        load([dpath sname sname2 'consump.mat'],'mclev','Zcon');
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %% Logmean biomass
        f1 = figure(1);
        subplot(3,1,1)
        plot(s-0.25,log10(squeeze(all_mean(1,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(1,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(1,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,log10(squeeze(all_mean(2,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(2,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(2,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title('M')
        
        subplot(3,1,3)
        plot(s,log10(squeeze(all_mean(3,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(3,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        title('L')
        xlabel('Sim')
        
        %% Logmean biomass
        f21 = figure(21);
        subplot(3,1,1)
        plot(s-0.25,log10(squeeze(all_mean(1,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(1,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(1,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,log10(squeeze(all_mean(2,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(2,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(2,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title('M')
        
        subplot(3,1,3)
        plot(s,log10(squeeze(all_mean(3,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(3,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('L')
        xlabel('Sim')
        
        
        %% Feeding level
        f2=figure(2);
        subplot(3,3,1)
        plot(s,mclev(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title(stage(1))
        
        subplot(3,3,2)
        plot(s,mclev(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title([loc ' ' stage(2)])
        
        subplot(3,3,3)
        plot(s,mclev(i,3),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title(stage(3))
        
        subplot(3,3,4)
        plot(s,mclev(i,4),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            ylabel('Feeding level')
        end
        title(stage(4))
        
        subplot(3,3,5)
        plot(s,mclev(i,5),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title(stage(5))
        
        subplot(3,3,6)
        plot(s,mclev(i,6),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title(stage(6))
        
        subplot(3,3,7)
        plot(s,mclev(i,7),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title(stage(7))
        
        subplot(3,3,8)
        plot(s,mclev(i,8),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title(stage(8))
        
        %% Growth rate (nu - energy for biomass production)
        f3 = figure(3);
        subplot(3,1,1)
        plot(s-0.25,Fmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dmgr(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmgr(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        ylabel('Mean growth/repro rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pmgr(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmgr(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        title('L')
        
        %% Consump per biomass (I)
        f4 = figure(4);
        subplot(3,1,1)
        plot(s-0.25,Fcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dcon(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dcon(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        ylabel('Mean consumption rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pcon(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dcon(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        title('L')
        
        %% Fraction zoop losses consumed
        f5 = figure(5);
        subplot(3,1,1)
        plot(s,z(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' Med zoo'])
        
        subplot(3,1,2)
        plot(s,z(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('Large zoo')
        ylabel('Fraction of flux consumed')
        
        subplot(3,1,3)
        plot(s,z(i,3),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('Detritus')
        
        
        %% Zoop con
        f12 = figure(12);
        subplot(2,1,1)
        plot(s,Zcon(i,1),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' Med zoo'])
        ylabel('Fraction of times overconsumed')
        
        subplot(2,1,2)
        plot(s,Zcon(i,2),'.k','MarkerSize',25); hold on;
        ylim([0 1])
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title([loc ' Large zoo'])
        ylabel('Fraction of times overconsumed')
        
        
        %% Size spectrum (sum stages)
        spec = nansum(all_mean(:,:,i),2);
        
        f7 = figure(7);
        stamp(cfile)
        plot(1:2:6,log10(spec),'LineWidth',2); hold on;
        xlim([0 6])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==ndp)
            legend(sims)
            legend('location','northwest')
            stamp(cfile)
        end
        ylabel('log Mean Biom (g m^-^2) in final year')
        xlabel('Size class')
        title(loc)
        
        %% Production (= nu * biom)
        f8 = figure(8);
        subplot(3,1,1)
        plot(s-0.25,Fprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dprod(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dprod(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pprod(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dprod(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('L')
        
        %% Reproduction
        f9 = figure(9);
        subplot(2,1,1)
        plot(s-0.25,Frep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Prep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Drep(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        ylabel('Mean repro rate (g g^-^1 d^-^1)')
        title(loc)
        
        subplot(2,1,2)
        plot(s-0.25,(Frep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Prep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Drep(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        ylabel('Mean biom reproduced (g d^-^1)')
        
        %% Metabolism
        f10 = figure(10);
        subplot(3,1,1)
        plot(s-0.25,Fmet(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pmet(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dmet(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,(Fmet(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pmet(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmet(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        ylabel('Mean metabolic rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(3,1,3)
        plot(s,(Pmet(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmet(3,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('L')
        
        %% Predation
        f11 = figure(11);
        subplot(2,1,1)
        plot(s-0.25,Fpred(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Ppred(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dpred(1,i),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(2,1,2)
        plot(s-0.25,(Fpred(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Ppred(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dpred(2,i)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,0,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        ylabel('Mean predation rate (g g^-^1 d^-^1) in final year',...
            'HorizontalAlignment','left')
        title('M')
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f16=figure(16);
        plot(s-0.1,log10(fishsp(1,i)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,log10(fishsp(2,i)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.1,log10(fishsp(3,i)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 2])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All stages'])
        
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        f15=figure(15);
        plot(s,log10(sumspec(i)),'k.','MarkerSize',25); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All fishes and stages'])
        
    end
    
    print(f1,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'Logmean_biomass.png'])
    print(f21,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'Logmean_biomass_axes.png'])
    print(f2,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'con_level.png'])
    print(f3,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'nu.png'])
    print(f4,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'consump.png'])
    print(f5,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'frac_zoop_loss.png'])
    print(f7,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'size_spec.png'])
    print(f8,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'prod.png'])
    print(f9,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'rep.png'])
    print(f10,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'met.png'])
    print(f11,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'pred.png'])
    print(f12,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'zoo_con.png'])
    print(f15,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'tot_mean_biomass_spec.png'])
    print(f16,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'tot_mean_biomass_type.png'])
    
end

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:ndp
        
        dpath = char(dp(s));
        load([dpath sname sname2 'consump.mat'],'mclev','Zcon');
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f17=figure(17);
        subplot(3,3,i)
        plot(s-0.1,log10(fishsp(1,i)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,log10(fishsp(2,i)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.1,log10(fishsp(3,i)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 2])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile)
        end
        if (i==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title([loc ' All stages'])
        
        
    end 
    
end
print(f17,'-dpng',[fpath sname sname2 cfile '_tot_mean_biomass_type_all_locs.png'])
    

