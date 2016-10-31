%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together
%Fishing rate vs. catch

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
% dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
% dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
% dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
% dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
% dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
% dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
% dpath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
% dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_NOnmort/'];
% dpath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
% dpath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
% dpath12 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
% dpath13 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
% dpath14 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
% dpath15 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
% dpath16 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];
% npath0 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/'];
% npath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish025/'];
% npath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05/'];
% npath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10/'];
% npath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20/'];
% npath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30/'];
% npath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40/'];
% npath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50/'];
% npath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60/'];
% npath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70/'];
% npath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish45/'];
% npath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish55/'];
% npath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish05_halfM/'];
% npath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/'];
% npath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/'];
% npath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/'];
% npath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/'];
% npath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/'];
% npath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/'];
% npath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/'];
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
% npath3 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish10_halfM/'];
% npath4 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish20_halfM/'];
% npath5 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish30_halfM/'];
% npath6 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish40_halfM/'];
% npath7 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish50_halfM/'];
% npath8 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish60_halfM/'];
% npath9 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish70_halfM/'];
% npath10 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish80_halfM/'];
% npath11 = [datap 'Dc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort_fish90_halfM/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fishM43L12/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140/'];
% npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150/'];
% npath19 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish160/'];
% npath20 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish170/'];
% npath21 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish180/'];
% npath22 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish190/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05_halfM/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10_halfM/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15_halfM/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20_halfM/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30_halfM/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40_halfM/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50_halfM/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60_halfM/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70_halfM/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80_halfM/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90_halfM/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100_halfM/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110_halfM/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120_halfM/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130_halfM/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140_halfM/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150_halfM/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish05_halfL/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish10_halfL/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish15_halfL/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish20_halfL/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish30_halfL/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish40_halfL/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish50_halfL/'];
% npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish60_halfL/'];
% npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish70_halfL/'];
% npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish80_halfL/'];
% npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish90_halfL/'];
% npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish100_halfL/'];
% npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish110_halfL/'];
% npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish120_halfL/'];
% npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish130_halfL/'];
% npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish140_halfL/'];
% npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_fish150_halfL/'];
% npath0 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01/'];
% npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish005/'];
% npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish01/'];
% npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish02/'];
% npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish03/'];
% npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish04/'];
% npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish05/'];
% npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_RE01_fish06/'];
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

% dp = {dpath1;dpath10;dpath11;dpath12;dpath13;dpath14;dpath15;dpath16};
% sims = {'0.0','0.01','0.025','0.05','0.10','0.20','0.30','0.40'};
% cfile = 'MFeqMP_fishing_catch';
% dp = {dpath2;dpath3;dpath4;dpath5;dpath6;dpath7;dpath8;dpath9};
% sims = {'0.0+N','0.01+N','0.025+N','0.05+N','0.10+N','0.20+N','0.30+N','0.40+N'};
% cfile = 'MFeqMP_fishing_nmort_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath10;npath7;...
%     npath11;npath8;npath9};
% sims = {'0.0','0.025','0.05','0.10','0.20','0.30',...
%     '0.40','0.45','0.50','0.55','0.60','0.70'};
% cfile = 'MFeqMP_MZ01_fishing_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath10;npath7;...
%     npath11;npath8;npath9};
% sims = {'0.0','0.025','0.05','0.10','0.20','0.30',...
%     '0.40','0.45','0.50','0.55','0.60','0.70'};
% cfile = 'MFeqMP_MZ01_fishing_catch';
% dp = {npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9};
% sims = {'0.05','0.10','0.20','0.30','0.40','0.50','0.60','0.70'};
% cfile = 'MFeqMP_MZ01_halfM_fishing_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;...
%     npath9};
% sims = {'0.025','0.05','0.075','0.10','0.20','0.30','0.40','0.50',...
%     '0.60','0.70'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_catch';
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11};
% sims = {'0.10','0.20','0.30','0.40','0.50','0.60','0.70','0.80','0.90'};
% cfile = 'Dc_MFeqMP_MZ01_fishing_half_catch';
% dp = {npath2;npath4;npath5;npath6;npath7;npath9;npath10;npath11;...
%     npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;npath20;...
%     npath21;npath22};
% sims = [0.1:0.1:1.9];
% cfile = 'Dc_MFeqMP_MZ01_BE05_fishing_catch';
% dp = {npath2;npath4;npath5;npath6;npath7;npath9;npath10;npath11;...
%     npath12;npath13;npath14;npath15;npath16;npath17};
% sims = [0.1:0.1:1.5];
% cfile = 'Dc_MFeqMP_MZ01_BE05_fishing_halfM_catch';
% dp = {npath2;npath4;npath5;npath6;npath7;npath9;npath10;npath11;...
%     npath12;npath13;npath14;npath15;npath16;npath17};
% sims = [0.1:0.1:1.5];
% cfile = 'Dc_MFeqMP_MZ01_BE05_fishing_halfL_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7};
% sims = {'0','.05','.1','.2','.3','.4','.5','.6'};
% cfile = 'Dc_MFeqMP_fcrit40_MZ01_BE05_RE01_fishing_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;};
% sims = {'0','.1','.2','.3','.4','.5'};
% cfile = 'Dc_MFeqMP_fcrit10_MZ01_BE05_RE01_mizernmort_fishing_catch';
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;};
% sims = {'0','.1','.2','.3','.4','.5'};
% cfile = 'Dc_mizer_all_MFeqMP_MZ01_BE05_RE01_fishing_catch';
dp = {npath0;npath1;npath2;npath3;npath4;npath5;};
sims = {'0','.1','.2','.3','.4','.5'};
cfile = 'Dc_JC_all_MFeqMP_MZ01_BE05_RE001_fishing_catch';

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

load('cmap_ppt_angles.mat')

%%
for d=1:length(dp)
    
    dpath = char(dp(d));
    load([dpath sname sname2 'lastyr_sum_mean_biom']);
    
    Ftot = sum(Ftotcatch(:));
    Ptot = sum(Ptotcatch(:));
    Dtot = sum(Dtotcatch(:));
    Tot = Ftot+Ptot+Dtot;
    
    f1 = figure(1);
    plot(d,Ftot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total F catch (g) in final year')
    
    f2 = figure(2);
    plot(d,Ptot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total P catch (g) in final year')
    
    f3 = figure(3);
    plot(d,Dtot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total D catch (g) in final year')
    
    f4 = figure(4);
    plot(d,Tot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total catch (g) in final year')
    
    
end
print(f1,'-dpng',[fpath sname sname2 cfile '_allF.png'])
print(f2,'-dpng',[fpath sname sname2 cfile '_allP.png'])
print(f3,'-dpng',[fpath sname sname2 cfile '_allD.png'])
print(f4,'-dpng',[fpath sname sname2 cfile '_allTotal.png'])


%%
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    close all
    %%
    for d=1:length(dp)
        
        dpath = char(dp(d));
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        
        %% Fishing
        %         MP_fish=nanmean(MP(lyr,28));
        %         MP_totcatch=nansum(MP(lyr,28));
        Ftot = sum(Ftotcatch(:,s));
        Ptot = sum(Ptotcatch(:,s));
        Dtot = sum(Dtotcatch(:,s));
        Tot = Ftot+Ptot+Dtot;
        
        f5 = figure(5);
        plot(d,Ftotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title([loc ' Total catch (g) in final year'])
            ylabel('MF')
            stamp(cfile)
        end
        
        f6 = figure(6);
        subplot(3,1,1)
        plot(d,Ptotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title([loc ' Total catch (g) in final year'])
            ylabel('MP')
            stamp(cfile)
        end
        subplot(3,1,2)
        plot(d,Ptotcatch(3,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('LP')
        end
        subplot(3,1,3)
        plot(d,Ptot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('All P')
        end
        
        f7 = figure(7);
        subplot(3,1,1)
        plot(d,Dtotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title([loc ' Total catch (g) in final year'])
            ylabel('MD')
            stamp(cfile)
        end
        subplot(3,1,2)
        plot(d,Dtotcatch(3,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('LD')
        end
        subplot(3,1,3)
        plot(d,Dtot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('All D')
        end
        
        f8 = figure(8);
        plot(d,Tot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total catch (g) in final year')
            stamp(cfile)
        end
        
    end
    
    print(f5,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'F.png'])
    print(f6,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'P.png'])
    print(f7,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'D.png'])
    print(f8,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'All.png'])
    
end


%%
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    
    %%
    for d=1:length(dp)
        
        dpath = char(dp(d));
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        
        %% Fishing
        Ftot = sum(Ftotcatch(:,s));
        Ptot = sum(Ptotcatch(:,s));
        Dtot = sum(Dtotcatch(:,s));
        Tot = Ftot+Ptot+Dtot;
        
        f9 = figure(9);
        subplot(3,3,s)
        plot(d,Tot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total catch (g) in final year')
            stamp(cfile)
        end
        
        f10 = figure(10);
        subplot(3,3,s)
        plot(d,Ftot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total F catch (g) in final year')
            stamp(cfile)
        end
        
        f11 = figure(11);
        subplot(3,3,s)
        plot(d,Ptot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total P catch (g) in final year')
            stamp(cfile)
        end
        
        f12 = figure(12);
        subplot(3,3,s)
        plot(d,Dtot,'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total D catch (g) in final year')
            stamp(cfile)
        end
        
    end
    
end
print(f9,'-dpng',[fpath sname sname2 cfile '_all_locs.png'])
print(f10,'-dpng',[fpath sname sname2 cfile '_allF_locs.png'])
print(f11,'-dpng',[fpath sname sname2 cfile '_allP_locs.png'])
print(f12,'-dpng',[fpath sname sname2 cfile '_allD_locs.png'])


