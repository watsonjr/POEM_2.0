% P:D ratio by LME 
% Climatology
% 150 years
% Saved as mat files
% Compare to SAU & Reg

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
cdir='/Volumes/GFDL/GCM_DATA/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp_zoop_det.mat']);

tlme = lme_mask_onedeg;
AREA_OCN = max(area,1);

% "lfile" never changes, has lme areas
lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = ['/Volumes/GFDL/NC/Matlab_new_size/' lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%% SAUP
%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Dlme_wcatch','Plme_wcatch');

Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);
lme_catch_all = nansum(lme_catch,3);

%1950-2006 SAUP average
id = find(yr>1950 & yr<=2006);

slme_mcatch = nanmean(lme_catch_all(id,:));
slme_mcatch = slme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));

%Top 10 yrs by LME SAUP 
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch_all(:,i),'descend');
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
sAlme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
sPlme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
sDlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

sFracPD = sPlme_mcatch10 ./ (sPlme_mcatch10 + sDlme_mcatch10);

%on grid
tlme = lme_mask_onedeg;
sFracPD_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    sFracPD_grid(lid) = sFracPD(L);
end

clear yr totcatch lme_catch Dlme_wcatch Plme_wcatch

%% Reg
%use weighted catches
load([spath 'Reg_LME_Catch_annual.mat']);

% sum of functional types in each
Flme_land_all = nansum(Flme_wland,3);
Plme_land_all = nansum(Plme_wland,3);
Dlme_land_all = nansum(Dlme_wland,3);
Flme_catch_all = nansum(Flme_wall,3);
Plme_catch_all = nansum(Plme_wall,3);
Dlme_catch_all = nansum(Dlme_wall,3);

%1950-2006 Reg average
id = find(yr>1950 & yr<=2006);

slme_mland = nanmean(lme_land(id,:));
slme_mland = slme_mland';
Fslme_mland = nanmean(Flme_land_all(id,:));
Fslme_mland = Fslme_mland';
Pslme_mland = nanmean(Plme_land_all(id,:));
Pslme_mland = Pslme_mland';
Dslme_mland = nanmean(Dlme_land_all(id,:));
Dslme_mland = Dslme_mland';

slme_mcatch = nanmean(lme_all(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

%Top 10 yrs by total
slme_mland10 = NaN*ones(size(slme_mland));
Flme_mland10 = NaN*ones(size(slme_mland));
Plme_mland10 = NaN*ones(size(slme_mland));
Dlme_mland10 = NaN*ones(size(slme_mland));
slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));
for i=1:66
    [sort_lme_land,ix] = sort(lme_land(:,i),'descend');
    sort_Flme_land = Flme_land_all(ix,i);
    sort_Plme_land = Plme_land_all(ix,i);
    sort_Dlme_land = Dlme_land_all(ix,i);
    slme_mland10(i) = nanmean(sort_lme_land(1:10));
    Flme_mland10(i) = nanmean(sort_Flme_land(1:10));
    Plme_mland10(i) = nanmean(sort_Plme_land(1:10));
    Dlme_mland10(i) = nanmean(sort_Dlme_land(1:10));
    
    [sort_lme_all,ix] = sort(lme_all(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_all(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
Alme_mland10 = slme_mland10 ./ lme_area_km2;
Flme_mland10 = Flme_mland10 ./ lme_area_km2;
Plme_mland10 = Plme_mland10 ./ lme_area_km2;
Dlme_mland10 = Dlme_mland10 ./ lme_area_km2;

Alme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

lFracPD = Plme_mland10 ./ (Plme_mland10 + Dlme_mland10);
cFracPD = Plme_mcatch10 ./ (Plme_mcatch10 + Dlme_mcatch10);

%on grid
tlme = lme_mask_onedeg;
lFracPD_grid = NaN*ones(180,360);
cFracPD_grid = NaN*ones(180,360);
for L=1:66
    lid = find(tlme==L);
    lFracPD_grid(lid) = lFracPD(L);
    cFracPD_grid(lid) = cFracPD(L);
end

%% Comparison stats
keep=[1:61,63];
load(['/Users/cpetrik/Dropbox/Princeton/POEM_other/poem_ms/',...
    'Stock_PNAS_catch_oceanprod_output.mat'],'notLELC')
notLELC = notLELC(notLELC<=63);

diffL = sFracPD_grid - lFracPD_grid;
diffC = sFracPD_grid - cFracPD_grid;

sAlme_mcatch10 = log10(sAlme_mcatch10);
sPlme_mcatch10 = log10(sPlme_mcatch10);
sDlme_mcatch10 = log10(sDlme_mcatch10);
Alme_mcatch10 = log10(Alme_mcatch10);
Plme_mcatch10 = log10(Plme_mcatch10);
Dlme_mcatch10 = log10(Dlme_mcatch10);

%r
rLfrac=corr(sFracPD(keep),lFracPD(keep));
rCfrac=corr(sFracPD(keep),cFracPD(keep));
rCA=corr(sAlme_mcatch10(keep),Alme_mcatch10(keep));
rCP=corr(sPlme_mcatch10(keep),Plme_mcatch10(keep));
rCD=corr(sDlme_mcatch10(keep),Dlme_mcatch10(keep));

%root mean square error
o=sFracPD(keep);
p=lFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseLfrac = sqrt(num/n);

o=sFracPD(keep);
p=cFracPD(keep);
n = length(o);
num=nansum((p-o).^2);
rmseCfrac = sqrt(num/n);

o=sAlme_mcatch10(keep);
p=Alme_mcatch10(keep);
n = length(o);
num=nansum((p-o).^2);
rmseCA = sqrt(num/n);

o=sPlme_mcatch10(keep);
p=Plme_mcatch10(keep);
n = length(o);
num=nansum((p-o).^2);
rmseCP = sqrt(num/n);

o=sDlme_mcatch10(keep);
p=Dlme_mcatch10(keep);
n = length(o);
num=nansum((p-o).^2);
rmseCD = sqrt(num/n);

%Fmed
FLfrac=10^(median(sFracPD(keep)-lFracPD(keep)));
FCfrac=10^(median(sFracPD(keep)-cFracPD(keep)));
FCA=10^(median(sAlme_mcatch10(keep)-Alme_mcatch10(keep)));
FCP=10^(median(sPlme_mcatch10(keep)-Plme_mcatch10(keep)));
FCD=10^(median(sDlme_mcatch10(keep)-Dlme_mcatch10(keep)));

% Table
fish_stat(1,1) = rLfrac;
fish_stat(2,1) = rmseLfrac;
fish_stat(3,1) = FLfrac;
fish_stat(1,2) = rCfrac;
fish_stat(2,2) = rmseCfrac;
fish_stat(3,2) = FCfrac;
fish_stat(1,3) = rCA;
fish_stat(2,3) = rmseCA;
fish_stat(3,3) = FCA;
fish_stat(1,4) = rCP;
fish_stat(2,4) = rmseCP;
fish_stat(3,4) = FCP;
fish_stat(1,5) = rCD;
fish_stat(2,5) = rmseCD;
fish_stat(3,5) = FCD;

Fstat = array2table(fish_stat,'RowNames',{'r','RMSE','Fmed'},...
    'VariableNames',{'PvsDLandings','PvsDCatch','AllCatch','PCatch','DCatch'});
writetable(Fstat,[spath 'LME_SAU_vs_Reg_stats.csv'],'Delimiter',',','WriteRowNames',true)
save([spath 'LME_SAU_vs_Reg_stats.mat'],'fish_stat')

%% Plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac
% ENTER -100 TO MAP ORIGIN LONG

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

revamp=colormap(cmocean('amp'));
revamp=flipud(revamp);
% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);
close all

x=-7:0.1:2;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Correlation
% P vs D landings
figure(1)
plot(x,x,'--k');hold on;
plot(x,x+0.1,':b'); hold on;
plot(x,x-0.1,':b'); hold on;
plot(x,x+0.2,':r'); hold on;
plot(x,x-0.2,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),sFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(lFracPD(lme),sFracPD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rLfrac)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseLfrac)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FLfrac)])
axis([-0.05 1.05 -0.05 1.05])
xlabel('Reg landings')
ylabel('SAU')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_SAU_vs_Reg_landings_comp_fracP.png'])

% P vs D catch
figure(2)
plot(x,x,'--k');hold on;
plot(x,x+0.1,':b'); hold on;
plot(x,x-0.1,':b'); hold on;
plot(x,x+0.2,':r'); hold on;
plot(x,x-0.2,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),sFracPD(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(cFracPD(lme),sFracPD(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rCfrac)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseCfrac)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FCfrac)])
axis([-0.05 1.05 -0.05 1.05])
xlabel('Reg catch')
ylabel('SAU')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_SAU_vs_Reg_catch_comp_fracP.png'])

%% All catch
figure(3)
plot(x,x,'--k');hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(Alme_mcatch10(lme),sAlme_mcatch10(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(Alme_mcatch10(lme),sAlme_mcatch10(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-3.5,1.75,['r = ' sprintf('%2.2f',rCA)])
text(-3.5,1.5,['RMSE = ' sprintf('%2.2f',rmseCA)])
text(-3.5,1.25,['Fmed = ' sprintf('%2.2f',FCA)])
axis([-4 2 -4 2])
xlabel('Reg catch')
ylabel('SAU')
title('log10 All catch (MT km^-^2)')
stamp('')
print('-dpng',[spath 'LME_SAU_vs_Reg_catch_comp_fracP.png'])

% P catch
figure(4)
plot(x,x,'--k');hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(Plme_mcatch10(lme),sPlme_mcatch10(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(Plme_mcatch10(lme),sPlme_mcatch10(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-6.5,0.75,['r = ' sprintf('%2.2f',rCP)])
text(-6.5,0.5,['RMSE = ' sprintf('%2.2f',rmseCP)])
text(-6.5,0.25,['Fmed = ' sprintf('%2.2f',FCP)])
axis([-7 1 -7 1])
xlabel('Reg catch')
ylabel('SAU')
title('log10 P catch (MT km^-^2)')
stamp('')
print('-dpng',[spath 'LME_SAU_vs_Reg_catch_comp_fracP.png'])

%% D catch
figure(5)
plot(x,x,'--k');hold on;
plot(x,x2h,':b'); hold on;
plot(x,x2l,':b'); hold on;
plot(x,x5h,':r'); hold on;
plot(x,x5l,':r'); hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(Dlme_mcatch10(lme),sDlme_mcatch10(lme),'o','MarkerSize',15,'color',tmap(tid(lme,2),:)); hold on;
    text(Dlme_mcatch10(lme),sDlme_mcatch10(lme),num2str(lme),...
        'Color','k','HorizontalAlignment','center'); hold on;
end
text(-3.5,1.75,['r = ' sprintf('%2.2f',rCD)])
text(-3.5,1.5,['RMSE = ' sprintf('%2.2f',rmseCD)])
text(-3.5,1.25,['Fmed = ' sprintf('%2.2f',FCD)])
axis([-4 2 -4 2])
xlabel('Reg catch')
ylabel('SAU')
title('log10 D catch (MT km^-^2)')
stamp('')
print('-dpng',[spath 'LME_SAU_vs_Reg_catch_comp_fracP.png'])


%% Subplot with maps and corr
figure(6)
% SAU
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('SAU Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffL)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('SAU - Reg difference')

%Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,lFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Reg Fraction LP landings')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(lFracPD(lme),sFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rLfrac)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseLfrac)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FLfrac)])
axis([0 1.05 0 1.05])
xlabel('Reg landings')
ylabel('SAU')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_fracPD_SAU_vs_Reg_landings_comp_subplot.png'])


figure(7)
% SAU
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,sFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
colorbar('Position',[0.01 0.5 0.5 0.05],'orientation','horizontal')                   %Set color axis if needed
set(gcf,'renderer','painters')
title('SAU Fraction LP catch')

% Diff
subplot('Position',[0.0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,diffC)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([-1 1]);
set(gcf,'renderer','painters')
title('SAU - Reg difference')

%Reg
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(geolat_t,geolon_t,cFracPD_grid)
cmocean('balance')
load coast;                     %decent looking coastlines
h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([0 1]);
set(gcf,'renderer','painters')
title('Reg Fraction LP catch')

%Corr
subplot('Position',[0.55 0.075 0.4 0.4])
plot(x,x,'--k');hold on;
for i=1:length(keep)
    lme=keep(i);
    plot(cFracPD(lme),sFracPD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
end
text(0.05,0.95,['r = ' sprintf('%2.2f',rLfrac)])
text(0.05,0.9,['RMSE = ' sprintf('%2.2f',rmseLfrac)])
text(0.05,0.85,['Fmed = ' sprintf('%2.2f',FLfrac)])
axis([0 1.05 0 1.05])
xlabel('Reg catch')
ylabel('SAU')
title('Fraction Large Pelagics')
stamp('')
print('-dpng',[spath 'LME_fracPD_SAU_vs_Reg_catch_comp_subplot.png'])



