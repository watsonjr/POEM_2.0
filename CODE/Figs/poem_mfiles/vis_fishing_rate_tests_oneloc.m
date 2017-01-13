% Visualize output of POEM
% Fishing spinup at one location
% 50 years

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish09/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0001_LD_fish14/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish01/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish02/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish03/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish04/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish05/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish06/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LP_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0005_LD_fish14/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish01/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish02/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish03/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish04/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish05/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish06/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LP_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish01/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish02/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish03/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish04/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish05/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish06/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LP_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0025_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish01/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish02/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish03/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish04/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish05/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish06/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LP_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE005_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish11/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish12/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish13/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish14/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish15/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish16/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish17/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish18/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish19/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_MF_fish20/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00005_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish11/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish12/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_MF_fish13/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00003_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish11/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish12/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish13/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish14/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish15/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish16/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish17/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish18/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish19/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_MF_fish20/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE00001_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE1000_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0500_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0100_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0050_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0025_LD_fish07/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE0010_LD_fish14/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00075_LD_fish14/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00050_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00025_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_MF_fish07/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00010_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00005_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish10/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish11/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish12/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish13/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish14/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish15/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish16/';
% npath31 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish17/';
% npath32 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish18/';
% npath33 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish19/';
% npath34 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00003_LD_fish20/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortH2_BE05_RE00001_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0050_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0025_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortM2_BE05_RE0010_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0050_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0025_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE0010_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00050_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00025_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00010_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00005_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00001_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_MF_fish10/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish01/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish02/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish03/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish04/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish05/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish06/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish07/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish08/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish09/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortJC2_BE05_RE00000_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0050/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE1000/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE1000/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE0100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab85_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE01/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE0010/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE0010/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE001/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmortH2_BE05_RE1000/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmortH2_BE05_RE1000/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE1000/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmortH2_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmortH2_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmortH2_BE05_RE0100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE0100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmortH2_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmortH2_BE05_RE0010/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmortH2_BE05_RE0010/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE1000/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE1000/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0010/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0010/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE1000/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE1000/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE1000/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0010/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0010/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_LD_fish10/';
 
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060_LD_fish10/';  

% npath0 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish10/';
  
% npath0 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010/';
% npath1 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish01/';
% npath2 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish02/';
% npath3 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish03/';
% npath4 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish04/';
% npath5 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish05/';
% npath6 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish06/';
% npath7 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish07/';
% npath8 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish08/';
% npath9 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish09/';
% npath10 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_MF_fish10/';
% npath11 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish01/';
% npath12 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish02/';
% npath13 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish03/';
% npath14 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish04/';
% npath15 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish05/';
% npath16 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish06/';
% npath17 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish07/';
% npath18 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish08/';
% npath19 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish09/';
% npath20 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LP_fish10/';
% npath21 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish01/';
% npath22 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish02/';
% npath23 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish03/';
% npath24 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish04/';
% npath25 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish05/';
% npath26 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish06/';
% npath27 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish07/';
% npath28 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish08/';
% npath29 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish09/';
% npath30 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010_LD_fish10/';

% npath0 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish01/';
% npath2 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish02/';
% npath3 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish03/';
% npath4 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish04/';
% npath5 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish05/';
% npath6 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish06/';
% npath7 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish07/';
% npath8 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish08/';
% npath9 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish09/';
% npath10 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fish10/';
% npath11 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish01/';
% npath12 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish02/';
% npath13 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish03/';
% npath14 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish04/';
% npath15 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish05/';
% npath16 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish06/';
% npath17 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish07/';
% npath18 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish08/';
% npath19 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish09/';
% npath20 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fish10/';
% npath21 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish01/';
% npath22 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish02/';
% npath23 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish03/';
% npath24 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish04/';
% npath25 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish05/';
% npath26 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish06/';
% npath27 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish07/';
% npath28 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish08/';
% npath29 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish09/';
% npath30 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_MF_fish10/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish01/';
% npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish02/';
% npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish03/';
% npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish04/';
% npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish05/';
% npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish06/';
% npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish07/';
% npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish08/';
% npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish09/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LP_fish10/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim_LD_fish10/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish01/';
% npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish02/';
% npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish03/';
% npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish04/';
% npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish05/';
% npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish06/';
% npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish07/';
% npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish08/';
% npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish09/';
% npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim_LD_fish10/';
npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim/';
npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish01/';
npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish02/';
npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish03/';
npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish04/';
npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish05/';
npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish06/';
npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish07/';
npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish08/';
npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish09/';
npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fish10/';

% dp = {npath0;npath1;npath2};

% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;...
%     npath11;npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;...
%     npath20;npath21};
% dp = {npath0};

%FORAGE
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fishing_catch';

%PELAGICS
% dp = {npath0;npath11;npath12;npath13;npath14;npath15;npath16;npath17;...
%     npath18;npath19;npath20};
% sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fishing_catch';

%DEMERSALS
dp = {npath0;npath21;npath22;npath23;npath24;npath25;npath26;npath27;npath28;npath29;npath30}; 
sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim_LD_fishing_catch';

%BASELINE
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;...
%      npath11;npath12;npath13;npath14;npath15;npath16;npath17};
% sims = {'.25-1','.5-1','.75-1','1-1','.25-.1','.5-.1','.75-.1','1-.1',...
%     '.25-.01','.5-.01','.75-.01','1-.01'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort2_BE05_LDpref_REtests';


sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

%%
for i=1:length(dp)
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
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
        %lyr=t((end-365+1):end);
        lyr=1:365;
        
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
        
        %% Fraction zoop losses consumed
        z(s,1) = nanmean(C(lyr,2));
        z(s,2) = nanmean(C(lyr,3));
        z(s,3) = nanmean(C(lyr,4));
        
        %% Size spectrum (sum stages)
        spec = nansum(all_mean(:,:,s),2);
        
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
        
        %% Total mortality w/ fishing
        Fmortf = Fpred + Fnat + Ffish;
        Pmortf = Ppred + Pnat + Pfish;
        Dmortf = Dpred + Dnat + Dfish;
        
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
        
    end
    
    save([dpath sname sname2 'lastyr_sum_mean_biom'],'Psum','Fsum',...
        'Dsum','Pmean','Fmean','Dmean','all_mean',...
        'Pmgr','Fmgr','Dmgr','Pcon','Fcon','Dcon','z','Pprod','Fprod','Dprod',...
        'Prep','Frep','Drep','Pmet','Fmet','Dmet','Ppred','Fpred','Dpred',...
        'Pnat','Fnat','Dnat','Pfish','Ffish','Dfish','Ptotcatch','Ftotcatch',...
        'Dtotcatch','Pgge','Fgge','Dgge');
    
end

%%
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    close all
    %%
    for d=1:length(dp)
        
        dpath = char(dp(d));
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
        
        %% Fishing
        Ftot = sum(Ftotcatch(:,s));
        Ptot = sum(Ptotcatch(:,s));
        Dtot = sum(Dtotcatch(:,s));
        Tot = Ftot+Ptot+Dtot;
        
        f5 = figure(5);
        plot(d,Ftotcatch(2,s),'.k','MarkerSize',20); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title([loc ' Total catch (g) in final year'])
            ylabel('MF')
            stamp(cfile2)
        end
        
        f6 = figure(6);
%         subplot(3,1,1)
%         plot(d,Ptotcatch(2,s),'.k','MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (d==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             title([loc ' Total catch (g) in final year'])
%             ylabel('MP')
%             stamp(cfile2)
%         end
%         subplot(3,1,2)
%         plot(d,Ptotcatch(3,s),'.k','MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (d==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             ylabel('LP')
%         end
%         subplot(3,1,3)
        plot(d,Ptot,'.k','MarkerSize',20); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('All P')
            stamp(cfile2)
        end
        
        f7 = figure(7);
%         subplot(3,1,1)
%         plot(d,Dtotcatch(2,s),'.k','MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (d==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             title([loc ' Total catch (g) in final year'])
%             ylabel('MD')
%             stamp(cfile2)
%         end
%         subplot(3,1,2)
%         plot(d,Dtotcatch(3,s),'.k','MarkerSize',15); hold on;
%         xlim([0 length(dp)+1])
%         if (d==length(dp))
%             set(gca,'XTick',1:length(dp),'XTickLabel',sims);
%             ylabel('LD')
%         end
%         subplot(3,1,3)
        plot(d,Dtot,'.k','MarkerSize',20); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('All D')
            stamp(cfile2)
        end
        
        f8 = figure(8);
        plot(d,Tot,'.k','MarkerSize',20); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total catch (g) in final year')
            stamp(cfile2)
        end
        
    end
    
    print(f5,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'F.png'])
    print(f6,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'P.png'])
    print(f7,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'D.png'])
    print(f8,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'All.png'])
    
end

%%
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    
    %%
    for d=1:length(dp)
        
        dpath = char(dp(d));
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
        
        %% Fishing
        Ftot = sum(Ftotcatch(:,s));
        Ptot = sum(Ptotcatch(:,s));
        Dtot = sum(Dtotcatch(:,s));
        Tot = Ftot+Ptot+Dtot;
        
        f9 = figure(9);
        subplot(4,3,s)
        plot(d,Tot,'.k','MarkerSize',25); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            if (s==4)
                ylabel('Total catch (g) in final year')
            end
            stamp(cfile2)
        end
        
        f10 = figure(10);
        subplot(4,3,s)
        plot(d,Ftot,'.k','MarkerSize',25); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            if (s==4)
                ylabel('Total F catch (g) in final year')
            end
            stamp(cfile2)
        end
        
        f11 = figure(11);
        subplot(4,3,s)
        plot(d,Ptot,'.k','MarkerSize',25); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            if (s==4)
                ylabel('Total P catch (g) in final year')
            end
            stamp(cfile2)
        end
        
        f12 = figure(12);
        subplot(4,3,s)
        plot(d,Dtot,'.k','MarkerSize',25); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            if (s==4)
                ylabel('Total D catch (g) in final year')
            end
            stamp(cfile2)
        end
        
    end
    
end
print(f9,'-dpng',[figp sname sname2 cfile2 '_all_locs.png'])
print(f10,'-dpng',[figp sname sname2 cfile2 '_allF_locs.png'])
print(f11,'-dpng',[figp sname sname2 cfile2 '_allP_locs.png'])
print(f12,'-dpng',[figp sname sname2 cfile2 '_allD_locs.png'])
