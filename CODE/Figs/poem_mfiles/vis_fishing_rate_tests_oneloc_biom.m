%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

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

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0050/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortI2_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100_LD_fish10/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish01/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish02/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish03/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish04/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish05/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish06/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish07/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish08/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish09/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100_LD_fish10/';

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

npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim/';
npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish01/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish02/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish03/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish04/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish05/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish06/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish07/';
npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish08/';
npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish09/';
npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_MF_fish10/';
npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish01/';
npath12 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish02/';
npath13 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish03/';
npath14 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish04/';
npath15 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish05/';
npath16 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish06/';
npath17 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish07/';
npath18 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish08/';
npath19 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish09/';
npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LP_fish10/';
npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish01/';
npath22 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish02/';
npath23 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish03/';
npath24 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish04/';
npath25 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish05/';
npath26 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish06/';
npath27 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish07/';
npath28 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish08/';
npath29 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish09/';
npath30 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fish10/';


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
cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim_LD_fishing_catch';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
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
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
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
            stamp(cfile2)
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
            stamp(cfile2)
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
            stamp(cfile2)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All stages'])
        
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        f15=figure(15);
        plot(s,log10(sumspec(i)),'k.','MarkerSize',25); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile2)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All fishes and stages'])
        
    end
    
    print(f1,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'Logmean_biomass.png'])
    print(f21,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'Logmean_biomass_axes.png'])
    print(f15,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'tot_mean_biomass_spec.png'])
    print(f16,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'tot_mean_biomass_type.png'])
    
end

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:ndp
        
        dpath = char(dp(s));
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
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
        ylim([-2 1])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
        end
        if (i==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title([loc ' All stages'])
        
        
    end 
    
end
print(f17,'-dpng',[figp sname sname2 cfile2 '_tot_mean_biomass_type_all_locs.png'])
    

