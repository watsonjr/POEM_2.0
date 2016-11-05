%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
datap = '/Volumes/GFDL/CSV/';

% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
% dpath = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];
% fpath = [figp 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];
% npath0 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath1 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish025/';
% npath2 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05/';
% npath3 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10/';
% npath4 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20/';
% npath5 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30/';
% npath6 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40/';
% npath7 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/';
% npath8 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60/';
% npath9 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70/';
% npath10 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish45/';
% npath11 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish55/';
% npath2 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05_halfM/';
% npath3 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/';
% npath4 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/';
% npath5 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/';
% npath6 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/';
% npath7 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/';
% npath8 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/';
% npath9 = 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/';
% npath0 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish025/';
% npath1 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05/';
% npath2 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish075/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20/';
% npath5 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30/';
% npath6 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40/';
% npath7 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/';
% npath8 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60/';
% npath9 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70/';
% npath10 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish80/';
% npath3 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/';
% npath4 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/';
% npath5 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/';
% npath6 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/';
% npath7 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/';
% npath8 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/';
% npath9 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/';
% npath10 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish80_halfM/';
% npath11 = 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish90_halfM/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish10/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish20/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish30/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20_fish40/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish05/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish10/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish20/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish30/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30_fish40/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish05/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish10/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish20/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish30/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20_fish40/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish05/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish10/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish20/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish30/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30_fish40/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fishM43L12/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish160/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish170/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish180/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish190/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05_halfM/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10_halfM/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15_halfM/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20_halfM/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30_halfM/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40_halfM/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50_halfM/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60_halfM/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70_halfM/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80_halfM/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90_halfM/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100_halfM/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110_halfM/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120_halfM/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130_halfM/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140_halfM/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150_halfM/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05_halfL/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10_halfL/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15_halfL/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20_halfL/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30_halfL/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40_halfL/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50_halfL/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60_halfL/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70_halfL/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80_halfL/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90_halfL/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100_halfL/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110_halfL/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120_halfL/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130_halfL/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140_halfL/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150_halfL/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish005/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish01/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish02/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish03/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish04/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish05/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish06/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01_fish05/';
% npath1 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish01/';
% npath2 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish02/';
% npath3 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish03/';
% npath4 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish04/';
% npath5 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort_BE05_RE01_fish05/';
% npath1 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish01/';
% npath2 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish02/';
% npath3 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish03/';
% npath4 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish04/';
% npath5 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE01_fish05/';
% npath0 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001/';
% npath1 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish01/';
% npath2 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish02/';
% npath3 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish03/';
% npath4 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish04/';
% npath5 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort_BE05_RE001_fish05/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish01/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish02/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish03/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish04/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish05/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish06/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish07/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MP_fish08/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish01/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish02/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish03/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish04/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish05/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish06/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish07/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_MD_fish08/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish01/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish02/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish03/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish04/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish05/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish06/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish07/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LP_fish08/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish01/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish02/';
% npath35 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish03/';
% npath36 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish04/';
% npath37 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish05/';
% npath38 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish06/';
% npath39 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish07/';
% npath40 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_LD_fish08/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_M_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish01/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish02/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish03/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish04/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish05/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish06/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish07/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_L_fish08/';
npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish01/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish02/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish03/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish04/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish05/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MF_fish06/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish01/';
npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish02/';
npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish03/';
npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish04/';
npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish05/';
npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MP_fish06/';
npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish01/';
npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish02/';
npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish03/';
npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish04/';
npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish05/';
npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_MD_fish06/';
npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish01/';
npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish02/';
npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish03/';
npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish04/';
npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish05/';
npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LP_fish06/';
npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish01/';
npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish02/';
npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish03/';
npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish04/';
npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish05/';
npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_LD_fish06/';
npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish01/';
npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish02/';
npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish03/';
npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish04/';
npath35 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish05/';
npath36 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_M_fish06/';
npath37 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish01/';
npath38 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish02/';
npath39 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish03/';
npath40 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish04/';
npath41 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish05/';
npath42 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01_L_fish06/';

% dp = {npath3};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% dp = {npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;...
%     npath9};
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% dp = {npath10};
dp = {npath37;npath38;npath39;npath40;npath41;npath42};
% dp = {npath9;npath10;npath11;npath12;npath13;npath14;npath15;npath16};
% dp = {npath20};

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

%%
for i=1:length(dp)
    close all
    
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
    load([dpath sname sname2 'consump.mat'],'mclev','Zcon');

    %% Zoop con
    
    figure(20)
    subplot(2,1,1)
    bar(Zcon(:,1),'k')
    ylim([0 1])
    set(gca,'XTickLabel',spots);
    title('Med zoo')
    ylabel('Fraction of times overconsumed')
    
    subplot(2,1,2)
    bar(Zcon(:,2),'k')
    ylim([0 1])
    set(gca,'XTickLabel',spots);
    title('Large zoo')
    ylabel('Fraction of times overconsumed')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_zoo_con.png'])
    
    %%
    Psum = NaN*ones(3,length(spots));
    Fsum = NaN*ones(2,length(spots));
    Dsum = Psum;
    Pmean = Psum;
    Fmean = Fsum;
    Dmean = Psum;
    Pmgr = Psum;
    Fmgr = Fsum;
    Dmgr = Psum;
    Pgge = Psum;
    Fgge = Fsum;
    Dgge = Psum;
    Pprod = Psum;
    Fprod = Fsum;
    Dprod = Psum;
    Pcon = Psum;
    Fcon = Fsum;
    Dcon = Psum;
    Prep = Fsum;
    Frep = Fsum;
    Drep = Fsum;
    Pmet = Psum;
    Fmet = Fsum;
    Dmet = Psum;
    Ppred = Psum;
    Fpred = Fsum;
    Dpred = Psum;
    Pnat = Psum;
    Fnat = Fsum;
    Dnat = Psum;
    Pfish = Psum;
    Ffish = Fsum;
    Dfish = Psum;
    Ptotcatch = Psum;
    Ftotcatch = Fsum;
    Dtotcatch = Psum;
    all_mean=NaN*ones(3,3,length(spots));
    z = NaN*ones(length(spots),3);
    
    %%
    for s=1:length(spots)
        %%
        loc = spots{s};
        lname = [sname2 loc '_'];
        SP = csvread([dpath sname lname 'Sml_p.csv']);
        SF = csvread([dpath sname lname 'Sml_f.csv']);
        SD = csvread([dpath sname lname 'Sml_d.csv']);
        MP = csvread([dpath sname lname 'Med_p.csv']);
        MF = csvread([dpath sname lname 'Med_f.csv']);
        MD = csvread([dpath sname lname 'Med_d.csv']);
        LP = csvread([dpath sname lname 'Lrg_p.csv']);
        LD = csvread([dpath sname lname 'Lrg_d.csv']);
        C = csvread([dpath sname lname 'Cobalt.csv']);
        
        t=1:length(SP);
        lyr=t((end-365+1):end);
        
        %% Final mean biomass in each size
        
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
        
        Psum(:,s) = P_sum;
        Fsum(:,s) = F_sum;
        Dsum(:,s) = D_sum;
        Pmean(:,s) = P_mean;
        Fmean(:,s) = F_mean;
        Dmean(:,s) = D_mean;
        
        
        all_mean(1:2,1,s) = F_mean;
        all_mean(:,2,s) = P_mean;
        all_mean(:,3,s) = D_mean;
        
        f1 = figure(1);
        subplot(3,3,s)
        plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 6])
        %ylim([-20 10])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==4)
            %legend('F','P','D')
            %legend('location','southeast')
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title(loc)
        xlabel('Stage')
        if (s==3)
            stamp(cfile)
        end
        
        f21 = figure(21);
        subplot(3,3,s)
        plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 6])
        ylim([-5 2])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title(loc)
        xlabel('Stage')
        if (s==3)
            stamp(cfile)
        end
        
        %% Feeding level
        f2=figure(2);
        subplot(3,3,s)
        bar(mclev(s,:),'k')
        ylim([0 1])
        set(gca,'XTickLabel',[]);
        for n=1:8
            text(n-0.5,-0.2,stage{n},'Rotation',45)
        end
        title(spots{s})
        if (s==4)
            ylabel('Feeding level')
        end
        stamp(cfile)
        
        %% Growth rate (nu - energy for biomass production)
        SP_mgr=nanmean(SP(lyr,15));
        SF_mgr=nanmean(SF(lyr,15));
        SD_mgr=nanmean(SD(lyr,15));
        MP_mgr=nanmean(MP(lyr,15));
        MF_mgr=nanmean(MF(lyr,15));
        MD_mgr=nanmean(MD(lyr,15));
        LP_mgr=nanmean(LP(lyr,15));
        LD_mgr=nanmean(LD(lyr,15));
        
        P_mgr=[SP_mgr;MP_mgr;LP_mgr];
        F_mgr=[SF_mgr;MF_mgr];
        D_mgr=[SD_mgr;MD_mgr;LD_mgr];
        
        Pmgr(:,s) = P_mgr;
        Fmgr(:,s) = F_mgr;
        Dmgr(:,s) = D_mgr;
        
        f3 = figure(3);
        subplot(1,3,1)
        plot(s-0.25,F_mgr(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_mgr(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_mgr(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean growth rate (g g^-^1 d^-^1) in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(F_mgr(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_mgr(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_mgr(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean growth/repro rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(1,3,3)
        plot(s,(P_mgr(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_mgr(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha3=gca;
            for n=1:9
                text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
        title('L')
        if (s==3)
            stamp(cfile)
        end
        
        %% Consump per biomass (I)
        SP_con=nanmean(SP(lyr,14));
        SF_con=nanmean(SF(lyr,14));
        SD_con=nanmean(SD(lyr,14));
        MP_con=nanmean(MP(lyr,14));
        MF_con=nanmean(MF(lyr,14));
        MD_con=nanmean(MD(lyr,14));
        LP_con=nanmean(LP(lyr,14));
        LD_con=nanmean(LD(lyr,14));
        
        P_con=[SP_con;MP_con;LP_con];
        F_con=[SF_con;MF_con];
        D_con=[SD_con;MD_con;LD_con];
        
        Pcon(:,s) = P_con;
        Fcon(:,s) = F_con;
        Dcon(:,s) = D_con;
        
        f4 = figure(4);
        subplot(3,3,s)
        plot(0.5:2:3.5,(F_con),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(1:2:6,(P_con),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(1.5:2:6.5,(D_con),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 6])
        %ylim([-25 15])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==4)
            %legend('F','P','D')
            %legend('location','southeast')
            ylabel('Mean consumption (g g^-^1 d^-^1) in final year')
        end
        title(loc)
        xlabel('Stage')
        if (s==3)
            stamp(cfile)
        end
        
        %% Fraction zoop losses consumed
        z(s,1) = nanmean(C(lyr,2));
        z(s,2) = nanmean(C(lyr,3));
        z(s,3) = nanmean(C(lyr,4));
        
        f5 = figure(5);
        subplot(3,3,s)
        bar(z(s,:)); hold on;
        xlim([0 4])
        ylim([0 1])
        set(gca,'XTick',1:3,'XTickLabel',{'MZ','LZ','Det'})
        if (s==4)
            ylabel('Fraction flux consumed')
        end
        title(loc)
        if (s==3)
            stamp(cfile)
        end
        
        %% Size spectrum (sum stages)
        spec = nansum(all_mean(:,:,s),2);
        
        f6 = figure(6);
        subplot(3,3,s)
        plot(1:2:6,log10(spec),'sk',...
            'MarkerFaceColor','k',...
            'MarkerSize',15); hold on;
        xlim([0 6])
        %ylim([-4 4])clo
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        title(loc)
        if (s==4)
            ylabel('log Mean Biom (g m^-^2) in final year')
        end
        xlabel('Size')
        if (s==3)
            stamp(cfile)
        end
        
        f7 = figure(7);
        stamp(cfile)
        %     plot(1:2:6,log10(spec),'color',cmap_ppt(s,:),...
        %         'LineWidth',2); hold on;
        plot(1:2:6,log10(spec),'LineWidth',2); hold on;
        xlim([0 6])
        %ylim([-25 15])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==9)
            legend(spots)
            legend('location','northwest')
        end
        ylabel('log Mean Biom (g m^-^2) in final year')
        xlabel('Size class')
        if (s==1)
            stamp(cfile)
        end
        
        %% Production (= nu * biom)
        SP_prod=nanmean(SP(lyr,24));
        SF_prod=nanmean(SF(lyr,24));
        SD_prod=nanmean(SD(lyr,24));
        MP_prod=nanmean(MP(lyr,24));
        MF_prod=nanmean(MF(lyr,24));
        MD_prod=nanmean(MD(lyr,24));
        LP_prod=nanmean(LP(lyr,24));
        LD_prod=nanmean(LD(lyr,24));
        
        P_prod=[SP_prod;MP_prod;LP_prod];
        F_prod=[SF_prod;MF_prod];
        D_prod=[SD_prod;MD_prod;LD_prod];
        
        Pprod(:,s) = P_prod;
        Fprod(:,s) = F_prod;
        Dprod(:,s) = D_prod;
        
        f8 = figure(8);
        subplot(1,3,1)
        plot(s-0.25,F_prod(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_prod(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_prod(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(F_prod(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_prod(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_prod(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(1,3,3)
        plot(s,(P_prod(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_prod(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha3=gca;
            for n=1:9
                text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
        title('L')
        if (s==1)
            stamp(cfile)
        end
        
        %% Reproduction
        F_rep(1,1)=nanmean(MF(lyr,18));
        D_rep(1,1)=nanmean(LD(lyr,18));
        P_rep(1,1)=nanmean(LP(lyr,18));
        F_rep(2,1)=nanmean(MF(lyr,1).*MF(lyr,18));
        D_rep(2,1)=nanmean(LD(lyr,1).*LD(lyr,18));
        P_rep(2,1)=nanmean(LP(lyr,1).*LP(lyr,18));
        
        Prep(:,s) = P_rep;
        Frep(:,s) = F_rep;
        Drep(:,s) = D_rep;
        
        f9 = figure(9);
        subplot(1,2,1)
        plot(s-0.25,F_rep(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_rep(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_rep(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean repro rate (g g^-^1 d^-^1) in final year')
        
        subplot(1,2,2)
        plot(s-0.25,(F_rep(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_rep(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_rep(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean biom reproduced (g d^-^1) in final year')
        if (s==1)
            stamp(cfile)
        end
        
        %% Metabolism
        SP_met=nanmean(SP(lyr,27));
        SF_met=nanmean(SF(lyr,27));
        SD_met=nanmean(SD(lyr,27));
        MP_met=nanmean(MP(lyr,27));
        MF_met=nanmean(MF(lyr,27));
        MD_met=nanmean(MD(lyr,27));
        LP_met=nanmean(LP(lyr,27));
        LD_met=nanmean(LD(lyr,27));
        
        P_met=[SP_met;MP_met;LP_met];
        F_met=[SF_met;MF_met];
        D_met=[SD_met;MD_met;LD_met];
        
        Pmet(:,s) = P_met;
        Fmet(:,s) = F_met;
        Dmet(:,s) = D_met;
        
        f10 = figure(10);
        subplot(1,3,1)
        plot(s-0.25,F_met(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_met(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_met(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(F_met(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_met(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_met(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
        title('M')
        
        subplot(1,3,3)
        plot(s,(P_met(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_met(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha3=gca;
            for n=1:9
                text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean metabolism (g g^-^1 d^-^1) in final year')
        title('L')
        if (s==3)
            stamp(cfile)
        end
        
        %% Predation
        SP_pred=nanmean(SP(lyr,25));
        SF_pred=nanmean(SF(lyr,25));
        SD_pred=nanmean(SD(lyr,25));
        MP_pred=nanmean(MP(lyr,25));
        MF_pred=nanmean(MF(lyr,25));
        MD_pred=nanmean(MD(lyr,25));
        LP_pred=nanmean(LP(lyr,25));
        LD_pred=nanmean(LD(lyr,25));
        
        P_pred=[SP_pred;MP_pred;LP_pred];
        F_pred=[SF_pred;MF_pred];
        D_pred=[SD_pred;MD_pred;LD_pred];
        
        Ppred(:,s) = P_pred;
        Fpred(:,s) = F_pred;
        Dpred(:,s) = D_pred;
        
        f11 = figure(11);
        subplot(1,2,1)
        plot(s-0.25,F_pred(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_pred(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_pred(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean predation rate (d^-^1) in final year')
        title('S')
        
        subplot(1,2,2)
        plot(s-0.25,(F_pred(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_pred(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_pred(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean predation rate (d^-^1) in final year')
        title('M')
        if (s==3)
            stamp(cfile)
        end
        
        %% Natural mortality
        Pnat(1,s)=nanmean(SP(lyr,26));
        Fnat(1,s)=nanmean(SF(lyr,26));
        Dnat(1,s)=nanmean(SD(lyr,26));
        Pnat(2,s)=nanmean(MP(lyr,26));
        Fnat(2,s)=nanmean(MF(lyr,26));
        Dnat(2,s)=nanmean(MD(lyr,26));
        Pnat(3,s)=nanmean(LP(lyr,26));
        Dnat(3,s)=nanmean(LD(lyr,26));
        
        %% Fishing
        MP_fish=nanmean(MP(lyr,28));
        MF_fish=nanmean(MF(lyr,28));
        MD_fish=nanmean(MD(lyr,28));
        LP_fish=nanmean(LP(lyr,28));
        LD_fish=nanmean(LD(lyr,28));
        
        P_fish=[0;MP_fish;LP_fish];
        F_fish=[0;MF_fish];
        D_fish=[0;MD_fish;LD_fish];
        
        Pfish(:,s) = P_fish;
        Ffish(:,s) = F_fish;
        Dfish(:,s) = D_fish;
        
        MP_totcatch=nansum(MP(lyr,28));
        MF_totcatch=nansum(MF(lyr,28));
        MD_totcatch=nansum(MD(lyr,28));
        LP_totcatch=nansum(LP(lyr,28));
        LD_totcatch=nansum(LD(lyr,28));
        
        P_totcatch=[0;MP_totcatch;LP_totcatch];
        F_totcatch=[0;MF_totcatch];
        D_totcatch=[0;MD_totcatch;LD_totcatch];
        
        Ptotcatch(:,s) = P_totcatch;
        Ftotcatch(:,s) = F_totcatch;
        Dtotcatch(:,s) = D_totcatch;
        
        f12 = figure(12);
        subplot(1,2,1)
        plot(s-0.25,F_totcatch(2),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_totcatch(2),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_totcatch(2),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Total catch (g) in final year')
        
        subplot(1,2,2)
        plot(s,P_totcatch(3),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_totcatch(3),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Total catch (g) in final year')
        if (s==1)
            stamp(cfile)
        end
        
        %% Total mortality w/o fishing
        %model lengths
        L(1) = 10^((log10(2)+log10(20))/2);
        L(2) = 10^((log10(20)+log10(200))/2);
        L(3) = 10^((log10(200)+log10(2000))/2);
        %model mass in grams
        M = 0.01 .* (0.1.*L).^3;
        %Andersen & Beyer mortality rate per year (natural + predation)
        %physiol mort * growth constant * M^-0.25
        AB = (0.35 .* 4.5 .* M.^(-0.25)) ./365;
        
        Fmort = Fpred + Fnat;
        Pmort = Ppred + Pnat;
        Dmort = Dpred + Dnat;
        
        f13=figure(13);
        subplot(1,3,1)
        plot(s-0.25,Fmort(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pmort(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dmort(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(1)*ones(11,1),'--k'); hold on;
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean mortality rate w/o fishing (d^-^1) in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(Fmort(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pmort(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmort(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(2)*ones(11,1),'--k'); hold on;
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        title('M')
        
        subplot(1,3,3)
        plot(s,(Pmort(3,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmort(3,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(3)*ones(11,1),'--k'); hold on;
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        title('L')
        
        %% Total mortality w/ fishing
        Fmortf = Fpred + Fnat + Ffish;
        Pmortf = Ppred + Pnat + Pfish;
        Dmortf = Dpred + Dnat + Dfish;
        
        f14=figure(14);
        subplot(1,3,1)
        plot(s-0.25,Fmortf(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,Pmortf(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,Dmortf(1,s),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(1)*ones(11,1),'--k'); hold on;
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean mortality rate w/fishing (d^-^1) in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(Fmortf(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(Pmortf(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmortf(2,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(2)*ones(11,1),'--k'); hold on;
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        title('M')
        
        subplot(1,3,3)
        plot(s,(Pmortf(3,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(Dmortf(3,s)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            plot(0:10,AB(3)*ones(11,1),'--k'); hold on;
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        title('L')
        
        %% Gross growth efficiency (= nu/consump) 
        SP_gge=nanmean(SP(lyr,15)./SP(lyr,14));
        SF_gge=nanmean(SF(lyr,15)./SF(lyr,14));
        SD_gge=nanmean(SD(lyr,15)./SD(lyr,14));
        MP_gge=nanmean(MP(lyr,15)./MP(lyr,14));
        MF_gge=nanmean(MF(lyr,15)./MF(lyr,14));
        MD_gge=nanmean(MD(lyr,15)./MD(lyr,14));
        LP_gge=nanmean(LP(lyr,15)./LP(lyr,14));
        LD_gge=nanmean(LD(lyr,15)./LD(lyr,14));
        
        P_gge=[SP_gge;MP_gge;LP_gge];
        F_gge=[SF_gge;MF_gge];
        D_gge=[SD_gge;MD_gge;LD_gge];
        
        Pgge(:,s) = P_gge;
        Fgge(:,s) = F_gge;
        Dgge(:,s) = D_gge;
        
        f15 = figure(15);
        subplot(1,3,1)
        plot(s-0.25,F_gge(1),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,P_gge(1),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,D_gge(1),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha1=gca;
            for n=1:9
                text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean gross growth efficiency in final year')
        title('S')
        
        subplot(1,3,2)
        plot(s-0.25,(F_gge(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,(P_gge(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_gge(2)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha2=gca;
            for n=1:9
                text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean gross growth efficiency in final year')
        title('M')
        
        subplot(1,3,3)
        plot(s,(P_gge(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.25,(D_gge(3)),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 10])
        ylim([-0.5 1])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if(s==9)
            ha3=gca;
            for n=1:9
                text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
            end
        end
        ylabel('Mean gross growth efficiency in final year')
        title('L')
        if (s==3)
            stamp(cfile)
        end
        
    end
    print(f1,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass.png'])
    print(f21,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass_axes.png'])
    print(f2,'-dpng',[fpath sname sname2 'All_oneloc_con_level.png'])
    print(f3,'-dpng',[fpath sname sname2 'All_oneloc_nu.png'])
    print(f4,'-dpng',[fpath sname sname2 'All_oneloc_consump.png'])
    print(f5,'-dpng',[fpath sname sname2 'All_oneloc_frac_zoop_loss.png'])
    print(f6,'-dpng',[fpath sname sname2 'All_oneloc_size_spec_sub.png'])
    print(f7,'-dpng',[fpath sname sname2 'All_oneloc_size_spec.png'])
    print(f8,'-dpng',[fpath sname sname2 'All_oneloc_prod.png'])
    print(f9,'-dpng',[fpath sname sname2 'All_oneloc_rep.png'])
    print(f10,'-dpng',[fpath sname sname2 'All_oneloc_met.png'])
    print(f11,'-dpng',[fpath sname sname2 'All_oneloc_pred.png'])
    print(f12,'-dpng',[fpath sname sname2 'All_oneloc_catch.png'])
    print(f13,'-dpng',[fpath sname sname2 'All_oneloc_mort_nof.png'])
    print(f14,'-dpng',[fpath sname sname2 'All_oneloc_mort_f.png'])
    print(f15,'-dpng',[fpath sname sname2 'All_oneloc_gge.png'])
    
    save([dpath sname sname2 'lastyr_sum_mean_biom'],'Psum','Fsum',...
        'Dsum','Pmean','Fmean','Dmean','all_mean',...
        'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
        'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
        'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
        'Dtotcatch','Pgge','Fgge','Dgge');
    
    %% Sum mean biom over stages
    fishsp = squeeze(nansum(all_mean));
    
    figure(16);
    plot((1-0.1):9,log10(fishsp(1,:)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:9,log10(fishsp(2,:)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot((1+0.1):10,log10(fishsp(3,:)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    ylim([-2 2])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-2.2,spots{n},'HorizontalAlignment','center')
    end
    ylabel('log10 Mean Biom (g m^-^2) in final year')
    title('All stages')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_type.png'])
    
    sumspec = squeeze(nansum(nansum(all_mean)));
    
    figure(17);
    subplot(2,1,1)
    plot(1:9,log10(sumspec),'k.','MarkerSize',25); hold on;
    xlim([0 10])
    %ylim([-2 1])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-0.6,spots{n},'HorizontalAlignment','center')
    end
    ylabel('log10 Mean Biom (g m^-^2) in final year')
    title('All fishes and stages')
    
    subplot(2,1,2)
    plot(1:9,(sumspec),'k.','MarkerSize',25); hold on;
    xlim([0 10])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-1,spots{n},'HorizontalAlignment','center')
    end
    ylabel('Mean Biom (g m^-^2) in final year')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec.png'])
    
    %% All on one
    figure(18)
    subplot(2,2,1)
    bar(log(Fmean))
    xlim([0 4])
    ylim([-10 3])
    xlabel('Stage')
    title('Forage')
    ylabel('log Mean Biomass (g m^-^2) in final year',...
        'HorizontalAlignment','right')
    
    subplot(2,2,2)
    bar(log(Pmean))
    xlim([0 4])
    ylim([-10 3])
    xlabel('Stage')
    title('Pel Pisc')
    
    subplot(2,2,3)
    bar(log(Dmean))
    xlim([0 4])
    ylim([-10 3])
    xlabel('Stage')
    title('Dem Pisc')
    
    subplot(2,2,4)
    bar(log(Fsum))
    xlim([0 0.5])
    ylim([-0.5 0.5])
    legend(spots)
    legend('location','west')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_biomass_spec.png'])
    
end




