% Visualize advection test cases

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

bio = csvread([dpath 'bio_2Dadvect_test_Pac.csv']);
grid = csvread('grid_csv.csv');
load('gridspec_forecast.mat');

totb = sum(bio,2);

%%
B1=NaN*ones(size(geolat_t));
B2=B1;
B3=B1;
B4=B1;
B5=B1;

%for n=1:size(bio,1);
    B1(grid(:,1))=bio(1,:);
    B2(grid(:,1))=bio(73,:);
    B3(grid(:,1))=bio(146,:);
    B4(grid(:,1))=bio(219,:);
    B5(grid(:,1))=bio(365,:);
%end

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


%% Start
figure(1)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=1')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B1); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([0 1e5])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_Pac1.png'])

figure(2)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=73')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B2); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([0 1e5])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_Pac2.png'])

figure(3)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=146')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B3); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([0 1e5])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_Pac3.png'])

figure(4)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=219')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B4); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([0 1e5])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_Pac4.png'])

%
figure(5)
% Create two axes
ax1 = axes;
pcolor(ax1,geolon_t,geolat_t,lmask); hold on;
shading flat
title('Day=365')
view(2)
ax2 = axes;
pcolor(ax2,geolon_t,geolat_t,B5); hold on;
shading flat
% Link them together
linkaxes([ax1,ax2])
% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
% Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')
caxis([0 1e5])
% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.1 .11 .75 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0675 .815]);
set(ax1,'XTick',xt,'XTickLabel',xl)
print('-dpng',[fpath 'advec_test_Pac5.png'])