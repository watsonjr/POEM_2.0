% Catch at different F rates

clear all
close all

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')

%% 10km grid
hpath = '/Volumes/GFDL/NEMURO/10km/';
load([hpath 'gridspec_10km.mat'],'LON','LAT');
load([hpath 'Data_grid_10km_hist.mat']);

GRD10 = GRD;
tot_area10_km2 = nansum(GRD10.AREA) * 1e-6;

%% FEISTY Output
cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
pp = ['/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/'];
fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/CalCurr/'];
ppath = [pp cfile '/CC/'];

F = 1:15;
mcatch10=NaN*ones(length(F),4);
for f=1:length(F)
    if f<10
        harv = ['All_fish0' num2str(F(f))];
    else
        harv = ['All_fish' num2str(F(f))];
    end
    
    load([fpath 'Means_Historic10km_' harv '_' cfile '.mat'],...
        'mf_my','mp_my','md_my',...
        'lp_my','ld_my');
    mf_my10 = mf_my;
    mp_my10 = mp_my;
    md_my10 = md_my;
    lp_my10 = lp_my;
    ld_my10 = ld_my;
    
    %% FEISTY LME biomass in MT
    % g/m2/d --> total g
    Amf_catch10 = mf_my10 .* GRD10.AREA * 365; %mean fish catch per yr
    Amp_catch10 = mp_my10 .* GRD10.AREA * 365;
    Amd_catch10 = md_my10 .* GRD10.AREA * 365;
    Alp_catch10 = lp_my10 .* GRD10.AREA * 365;
    Ald_catch10 = ld_my10 .* GRD10.AREA * 365;
    
    %%
    AF_catch10 = Amf_catch10;
    AP_catch10 = Amp_catch10 + Alp_catch10;
    AD_catch10 = Amd_catch10 + Ald_catch10;
    All_catch10 = AF_catch10 + AP_catch10 + AD_catch10;
    
    %% g --> MT
    mcatch10(f,1) = nansum(AF_catch10) * 1e-6;
    mcatch10(f,2) = nansum(AP_catch10) * 1e-6;
    mcatch10(f,3) = nansum(AD_catch10) * 1e-6;
    mcatch10(f,4) = nansum(All_catch10) * 1e-6;
    
end

% MT/km2
mcatch10_km2 = mcatch10 ./ tot_area10_km2;



%% Plots
frate = 0.1:0.1:1.5;

figure(1)
plot(frate,mcatch10(:,1),'k','LineWidth',2)
ylabel('Mean catch (MT)')
xlabel('Fishing rate (yr^-^1)')
title('F')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist10km_yield_curve_F.png'])

figure(2)
plot(frate,mcatch10(:,2),'k','LineWidth',2)
ylabel('Mean catch (MT)')
xlabel('Fishing rate (yr^-^1)')
title('P')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist10km_yield_curve_P.png'])

figure(3)
plot(frate,mcatch10(:,3),'k','LineWidth',2)
ylabel('Mean catch (MT)')
xlabel('Fishing rate (yr^-^1)')
title('D')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist10km_yield_curve_D.png'])

figure(4)
plot(frate,mcatch10(:,4),'k','LineWidth',2)
ylabel('Mean catch (MT)')
xlabel('Fishing rate (yr^-^1)')
title('All')
stamp([harv '_' cfile])
print('-dpng',[ppath 'Hist10km_yield_curve_All.png'])

%% RAM stocks in California or P Coast (30 total)
%sum(CCram$mLandings,na.rm=T)
% 56,754.98 MT
%sum(CCram$mCatch,na.rm=T)
% 447,356.9 MT
%sum(CCram$mBiomass,na.rm=T)
% 2,596,265
%mean(CCram$mnF,na.rm=T)
% 0.06590057
%mean(CCram$mdF,na.rm=T)
% 0.060325
%mean(CCram$mnER,na.rm=T)
% 0.08126313
%mean(CCram$mdER,na.rm=T)
% 0.07224141

% F=0.1 MAY BE A MORE REALISTIC FISHING RATE


