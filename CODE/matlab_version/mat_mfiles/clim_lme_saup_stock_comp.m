% Compare POEM vs. SAUP
% to Stock et al PNAS model vs. SAUP

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

cfile = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_noCC_RE00100';
harv = 'All_fish03_Juve09';

ppath = [pp cfile '/'];
dpath = [dp cfile '/'];

%Reconciling Fisheries Catch and Ocean Productivity
%***TEMPLATE FOR FEEDBACK, PENDING FINAL CHECKS***
%model: 4
%alpha: 0.14
%TE_p: 0.13
%TE_b: 0.40
%f_T: 0.74
%T_100,warm: 19.99
%All fluxes in g C m-2 yr-1, Temperature in degrees celsius
%cols = LME  NPP   MESOZP  FDET   TLeq     T  modcatch SAUcatch
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'])

% POEM
%tab  = MT/km2
%tab2 = gC/m2
%cols = 'LME','ZP','Det','Bent','T','Pcatch','Scatch'
load([dpath 'LME_prod_catch_SAUP_' cfile '.mat'],'tab','tab2')

% Temps w/color
load([cpath 'LME_clim_temp_zoop_det.mat']);
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

%% Temp
figure(1)
plot(-2:30,-2:30,'--k'); hold on;
plot(StockPNAS(:,6),tab(:,5),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME temp')

% Zoop
figure(2)
plot(0:40,0:40,'--k'); hold on;
plot(StockPNAS(:,3),tab(:,2),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME mesozoop')

% Det
figure(3)
plot(0:50,0:50,'--k'); hold on;
plot(StockPNAS(:,4),tab(:,3),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME bottom detritus')

% Modeled catch
figure(4)
plot(0:0.1:1.2,0:0.1:1.2,'--k'); hold on;
plot(StockPNAS(:,7),tab2(:,6),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME model catch')

% SAUP catch
figure(5)
plot(0:0.1:1.2,0:0.1:1.2,'--k'); hold on;
plot(StockPNAS(:,8),tab2(:,7),'.k','MarkerSize',20)
xlabel('Stock PNAS')
ylabel('POEM')
title('LME SAU catch')
