% Visualize Matlab diffusion test cases

clear all
close all

dpath = '/Volumes/GFDL/CSV/advect_tests/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/advect_tests/';

load('/Volumes/GFDL/CSV/advect_tests/POEM_diff_test_gradT_Atl_dt30min');
cname = 'Atl_dt30min_k600_b100';

grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');


%% plot info
% Land
surf_tmask = tmask(:,:,1);
lmask = surf_tmask;
lmask(lmask==0) = 999;
lmask(lmask==1) = NaN;

%axes labels
xt = -250:50:50;
xl = xt;
xl(xl<-180) = xl(xl<-180) + 350;


%% Global flat
figure
    surf(geolon_t,geolat_t,TF2);
    view(2);
    shading flat;
    colorbar;
    caxis([0 1e2]);
    colormap('jet')
    title(['Day 365 Year 1'])
    %print('-dpng',[fpath 'diff_test_' cname '_' num2str(t(n)) '.png'])


%% Arctic projection
nd = size(gT,3);
%t = 1:35:nd;
%t = 1:3:15;
%t = round(t);
t = 1:2189:nd;
d = t/48;
d = round(d);

for n=1:length(t)
    B1 = gT(:,:,t(n));
    
    figure
    m_proj('stereographic','lat',90,'long',30,'radius',30);
    m_pcolor(geolon_t,geolat_t,B1);
    shading flat
    colorbar
    colormap('jet')
    caxis([-1e-3 1e-3])
    m_grid('xtick',12,'tickdir','out','ytick',[70 80],'linest','-');
    m_coast('patch',[.7 .7 .7],'edgecolor','k');
    title(['Day ' num2str(d(n)) ' Year 1'])
    print('-dpng',[fpath 'Mat_diff_test_' cname '_arcticproj_' num2str(d(n)) '.png'])
end

