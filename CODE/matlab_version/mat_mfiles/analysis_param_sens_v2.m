% POEM output at all locations
% Only last 12 months of 150 years saved 

clear all
close all

%GFDL/NC/Matlab_new_size/Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100/param_sens/
cfile = 'Dc_enc50-b210_m4-b210-k060_c50-b210_D075_J075_A075_Sm025_nmort1_BE08_noCC_RE00100';
nfile = ['/Volumes/GFDL/NC/Matlab_new_size/',cfile,'/param_sens/'];

load([nfile 'Climatol_All_fish03_means_param_sens_v2.mat'])
ptext{29}='kap75';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);

%% Outputs
F = SF + MF;
P = SP + MP + LP;
D = SD + MD + LD;
All = F + P + D;
% Relative biomass of each group
v1 = log10(nansum(F)) - log10(nansum(F(:,1)));
v2 = log10(nansum(P)) - log10(nansum(P(:,1)));
v3 = log10(nansum(D)) - log10(nansum(D(:,1)));
% Total biomass
v4 = log10(nansum(All)) - log10(nansum(All(:,1)));
% Weighted latitudinal center of biomass
vlat = lat(ID);
nhem = find(vlat>0);
shem = find(vlat<0);
totb = repmat(nansum(All),length(ID),1);
nh = All(nhem,:) .* vlat(nhem) ./ totb(nhem,:);
NH = nansum(nh);
sh = All(shem,:) .* vlat(shem) ./ totb(shem,:);
SH = nansum(sh);
latc = (NH-SH)/2;
v5 = log10(latc) - log10(latc(:,1));

% Combined
v = sqrt(v1.^2 + v2.^2 + v3.^2 + v4.^2 + v5.^2);




