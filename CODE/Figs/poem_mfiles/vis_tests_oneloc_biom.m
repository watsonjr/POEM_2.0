%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE1000/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE1000/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE1000/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE0100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD025_nmort0_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD050_nmort0_BE05_RE0010/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD075_nmort0_BE05_RE0010/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0010/';
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
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE1000/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmort0_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmort0_BE05_RE0010/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmort0_BE05_RE0010/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmort0_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE1000/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE1000/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE1000/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE0100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D025_nmortH2_BE05_RE0010/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_RE0010/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D075_nmortH2_BE05_RE0010/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_LD100_nmortH2_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortM2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortJC2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D050_nmort0_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortM2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortJC2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D030_nmort0_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortM2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortJC2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmortH2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D010_nmort0_BE05_RE0010/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_A100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_A100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_A100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_A100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_A100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_A100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_A100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_A100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_A100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_A100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_A100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_A100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_A100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_A100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_A100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_A100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_A100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_A100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_A100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_A100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_A100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_A100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_A100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_A100/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_HA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_HA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_HA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_HA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_HA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_HA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_HA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_HA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_HA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_HA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_HA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_HA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_HA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_HA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_HA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_HA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_HA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_HA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_HA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_HA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_HA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_HA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_HA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_HA100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_MA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_MA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_MA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_MA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_MA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_MA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_MA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_MA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_MA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_MA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_MA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_MA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_MA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_MA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_MA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_MA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_MA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_MA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_MA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_MA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_MA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_MA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_MA100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_JCA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_JCA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_JCA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_JCA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_JCA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_JCA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_JCA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_JCA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_JCA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_JCA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_JCA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_JCA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_JCA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_JCA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_JCA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_JCA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_JCA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_JCA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_JCA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_JCA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_JCA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_JCA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_JCA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_JCA100/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_HA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_HA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_HA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_HA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_HA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_HA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_HA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_HA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_HA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_HA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_HA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_HA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_HA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_HA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_HA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_HA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_HA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_HA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_HA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_HA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_HA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_HA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_HA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_HA100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_MA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_MA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_MA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_MA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_MA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_MA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_MA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_MA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_MA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_MA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_MA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_MA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_MA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_MA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_MA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_MA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_MA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_MA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_MA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_MA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_MA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_MA100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_JCA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_JCA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_JCA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_JCA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_JCA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_JCA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_JCA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_JCA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_JCA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_JCA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_JCA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_JCA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_JCA100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_JCA100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_JCA100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_JCA100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_JCA100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_JCA100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_JCA100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_JCA100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_JCA100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_JCA100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_JCA100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_JCA100/';

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

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE1000_BA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0500_BA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0100_BA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0090_BA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0080_BA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0070_BA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0060_BA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0050_BA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0040_BA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0030_BA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0020_BA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_nmort0_BE05_RE0010_BA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_BA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_BA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_BA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_BA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_BA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_BA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_BA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_BA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_BA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_BA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_BA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE1000_BA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0100_BA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0090_BA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0080_BA/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0070_BA/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0060_BA/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0050_BA/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0040_BA/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0030_BA/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0020_BA/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_BA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_nmort0_BE05_RE1000_BA/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_nmort0_BE05_RE0500_BA/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_nmort0_BE05_RE0100_BA/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_nmort0_BE05_RE0090_BA/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_nmort0_BE05_RE0080_BA/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_Bassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_Bassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_Bassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_Bassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_Bassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_Bassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_Bassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_Bassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_Bassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_Bassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_Bassim/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_BAassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmort0_BE05_RE0010_BAassim/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_BAassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_BAassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_BAassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_BAassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_BAassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAassim/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_BAboth/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_BAboth/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_BAboth/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_BAboth/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_BAboth/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_BAboth/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_BAboth/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_BAboth/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_BAboth/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_BAboth/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_BAboth/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_BAboth/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_BAboth/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_BAboth/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_BAboth/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_BAboth/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_BAboth/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_BAboth/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_BAboth/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_BAboth/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_BAboth/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_BAboth/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_BAboth/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_BAboth/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_BAboth/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_BAboth/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_BAboth/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_BAboth/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_BAboth/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_BAboth/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_BAboth/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_BAboth/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_BAboth/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_BAboth/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAboth/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_BAboth/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_BAboth/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_BAboth/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_BAboth/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_BAboth/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_BAboth/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_BAboth/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_BAboth/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_BAboth/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_BAboth/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_BAboth/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_BAboth/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAboth/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_BAbent/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_BAbent/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_BAbent/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_BAbent/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_BAbent/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_BAbent/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_BAbent/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_BAbent/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_BAbent/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_BAbent/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_BAbent/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_BAbent/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_BAbent/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_BAbent/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_BAbent/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_BAbent/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_BAbent/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_BAbent/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_BAbent/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_BAbent/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_BAbent/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_BAbent/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_BAbent/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_BAbent/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_BAbent/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_BAbent/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_BAbent/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_BAbent/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_BAbent/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_BAbent/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_BAbent/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_BAbent/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_BAbent/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_BAbent/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_BAbent/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_BAbent/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_Bassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_Bassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_Bassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_Bassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_Bassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_Bassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_Bassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_Bassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_Bassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_Bassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_Bassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_Bassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_Bassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_Bassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_Bassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_Bassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_Bassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_Bassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_Bassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_Bassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_Bassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_Bassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_Bassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_Bassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_Bassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_Bassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_Bassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_Bassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_Bassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_Bassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_Bassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_Bassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_Bassim/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_Bassim/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_Bassim/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_Bassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_Bassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_Bassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_Bassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_Bassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_Bassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_Bassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_Bassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_Bassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_Bassim/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE1000/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0500/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0090/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0080/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0070/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0060/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0050/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0040/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0030/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0020/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_nmortPW2_BE05_RE0010/';

% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00010_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00010_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00090_BAassim/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00080_BAassim/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00070_BAassim/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00060_BAassim/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00050_BAassim/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_BAassim/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00030_BAassim/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00020_BAassim/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_BAassim/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_K100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_K100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_K100/';   
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_K100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_K100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_K100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_K100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_K100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_K100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_K100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_K100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_K100/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_K100/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_K100/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_K100/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_K75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_K75/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_K75/';           
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_K75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_K75/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_K75/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_K75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_K75/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_K75/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_K75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_K75/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_K75/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_K75/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_K75/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_K75/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_K50/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_K50/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_K50/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_K50/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_K50/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_K50/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_K50/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_K50/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_K50/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_K50/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_K50/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_K50/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_K50/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_K50/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_K50/';

% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE1000_K25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0500_K25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0100_K25/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0090_K25/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0080_K25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0070_K25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0060_K25/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0050_K25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0040_K25/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0030_K25/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0020_K25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE0010_K25/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE1000_K25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0500_K25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0100_K25/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0090_K25/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0080_K25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0070_K25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0060_K25/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0050_K25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0040_K25/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0030_K25/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0020_K25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE0010_K25/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE1000_K25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0500_K25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0100_K25/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0090_K25/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0080_K25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0070_K25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0060_K25/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0050_K25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0040_K25/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0030_K25/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0020_K25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE0010_K25/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE1000_K25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0500_K25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0100_K25/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0090_K25/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0080_K25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0070_K25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0060_K25/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0050_K25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0040_K25/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0030_K25/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0020_K25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE0010_K25/';
% npath0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE1000_K25/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0500_K25/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0100_K25/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0090_K25/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0080_K25/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0070_K25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0060_K25/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0050_K25/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0040_K25/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0030_K25/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0020_K25/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE0010_K25/';

% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_K100/';   
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00010_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00010_K100/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00090_K100/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00080_K100/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00070_K100/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00060_K100/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00050_K100/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_K100/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00030_K100/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00020_K100/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_K100/';

% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_K75/';   
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00010_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00010_K75/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00090_K75/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00080_K75/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00070_K75/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00060_K75/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00050_K75/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_K75/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00030_K75/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00020_K75/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_K75/';

% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_RE00010_K50/';   
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortH2_BE05_RE00010_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortJC2_BE05_RE00010_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortM2_BE05_RE00010_K50/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00090_K50/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00080_K50/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00070_K50/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00060_K50/';
% npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00050_K50/';
% npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00040_K50/';
% npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00030_K50/';
% npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00020_K50/';
% npath11 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_RE00010_K50/';

npath1='Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D100_nmort0_BE05_CC275_RE0500/';
npath2='Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE05_CC275_RE0500/';
dp = {npath1;npath2};
sims = {'.30','.40'};
cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_D100_nmort0_BE05_CC275_RE0500_fcrit_tests';

% npath1 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100/';
% npath2 = 'Dc_TrefO_JC_all_MFeqMP_MZ01_nmort2_BE05_RE0010/';
% npath3 = 'Dc_TrefO_mizer_all_MFeqMP_MZ01_nmort2_BE05_RE0100/';
% dp = {npath1;npath2;npath3};
% sims = {'Hartvig','J&C','mizer'};
% cfile2 = 'Dc_TrefO_MFeqMP_MZ01_nmort2_BE05_all_param_tests';

%BASELINE D pref
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;...
%      npath11};
% sims = {'.25-1','.5-1','.75-1','1-1','.25-.1','.5-.1','.75-.1','1-.1',...
%     '.25-.01','.5-.01','.75-.01','1-.01'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_D050_nmortH2_BE05_LDpref_REtests';

%BASELINE RE & mort
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11};
% sims = {'1','.5','.1','.09','.08','.07','.06','.05','.04','.03','.02','.01'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmortPW2_BE05_K25_REtests';

%lower RE & mort
% dp = {npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11};
% sims = {'.009','.008','.007','.006','.005','.004','.003','.002','.001'};
% cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_nmort0_BE05_K50_lowREtests';


%FORAGE
% dp = {npath0;npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10};
% sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile2 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_MF_fishing_comp';

%PELAGICS
% dp = {npath0;npath11;npath12;npath13;npath14;npath15;npath16;npath17;...
%     npath18;npath19;npath20};
% sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile2 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LP_fishing_comp';

%DEMERSALS
% dp = {npath0;npath21;npath22;npath23;npath24;npath25;npath26;npath27;npath28;npath29;npath30}; %;...
% %npath31;npath32;npath33;npath34};
% sims = {'0','.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
% cfile2 = 'Dc_TrefO_Hartvig_all_MFeqMP_MZ01_nmort2_BE05_RE0100_LD_fishing_comp';


sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

ndp = length(dp);

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
        subplot(4,3,i)
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
    

