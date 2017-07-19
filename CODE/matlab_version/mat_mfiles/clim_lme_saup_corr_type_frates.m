%POEM catch vs. SAUP catch by LME
%Find fishing selectivity that gives highest corr for each type

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/';

Pdir = '/Volumes/GFDL/POEM_JLD/esm26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

%use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

%Colormap
load('MyColormaps.mat')
load('cmap_ppt_angles.mat')
cmap1(1,:)=[1 1 1];
cmap1(2,:)=cmap_ppt(1,:);
cmap1(3,:)=cmap_ppt(3,:);
cmap1(4,:)=cmap_ppt(5,:);

cmap2(1,:)=cmap_ppt(1,:);
cmap2(2,:)=cmap_ppt(3,:);
cmap2(3,:)=cmap_ppt(5,:);

cmap3(1,:)=[0 0 0];
cmap3(2,:)=cmap_ppt(1,:);
cmap3(3,:)=cmap_ppt(3,:);
cmap3(4,:)=cmap_ppt(5,:);

cmap4(1,:)=cmap_ppt(3,:);
cmap4(2,:)=cmap_ppt(1,:);
cmap4(3,:)=cmap_ppt(5,:);

AREA_OCN = max(area,1);

lfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100';
lpath = [dp lfile '/'];
load([lpath 'LME_clim_fished03_' lfile '.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%%
Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));
%Top 10 yrs SAUP
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

%% Assign a color to each LME based on temp
tmap=colormap(jet(66));
lme_ptemp(:,2)=1:length(lme_ptemp);
[B,I] = sort(lme_ptemp(:,1));
I(:,2)=1:length(lme_ptemp);
[B2,I2] = sort(I(:,1));
tid = I(I2,:);

keep = [1:9,11:34,36:38,46:54,59:60,62,65];

x=-6:0.5:8;
x2h = x+log10(2);
x2l = x-log10(2);
x5h = x+log10(5);
x5l = x-log10(5);

%% Loop over POEM params
fqs = 0.5:0.5:3.5;
pqs = 0.25:0.25:1.5;
dqs = 1:6;
frate = 0.1;

rall = NaN*ones(length(fqs),length(pqs),length(dqs));
rF = rall;
rP = rall;
rD = rall;
rmse = rall;
rmseF = rall;
rmseP = rall;
rmseD = rall;

sname = 'Clim_';
dp = 'Dc_enc70-b200_cm20_m-b175-k09_fcrit20_c-b250_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
datap = '/Volumes/GFDL/CSV/Matlab_new_size/';
dpath = [datap char(dp) '/'];
ppath = [pp char(dp) '/'];

for fq=1%:length(fqs)
    for pq=1:3;%length(pqs)
        for dq=1:length(dqs)
            MFsel = fqs(fq);
            LPsel = pqs(pq);
            LDsel = dqs(dq);
            tF = num2str(1000+int64(100*frate*MFsel));
            tP = num2str(1000+int64(100*frate*LPsel));
            tD = num2str(1000+int64(100*frate*LDsel));
            
            charv = ['fish_F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
            harv = ['F',tF(2:end),'_P',tP(2:end),'_D',tD(2:end)];
            cfile = [harv,'_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];

            load([dpath 'LME_clim_fished_',harv,'_loop_' dp '.mat'],'lme_mcatch');
            
            close all 
            
            %% POEM LME biomass in MT
            plme_mcatch = nansum(lme_mcatch,2) * 1e-6;
            plme_Fmcatch = (lme_mcatch(:,1)) * 1e-6;
            plme_Pmcatch = (lme_mcatch(:,2)+lme_mcatch(:,4)) * 1e-6;
            plme_Dmcatch = (lme_mcatch(:,3)+lme_mcatch(:,5)) * 1e-6;
            % MT/km2
            plme_mcatch = plme_mcatch ./ lme_area_km2;
            plme_Fmcatch = plme_Fmcatch ./ lme_area_km2;
            plme_Pmcatch = plme_Pmcatch ./ lme_area_km2;
            plme_Dmcatch = plme_Dmcatch ./ lme_area_km2;
            
            %log10 Difference
            diff_catch = (log10(plme_mcatch) - log10(slme_mcatch10)) ./ log10((slme_mcatch10+eps));
            Fdiff_catch = (log10(plme_Fmcatch) - log10(Flme_mcatch10)) ./ log10((Flme_mcatch10+eps));
            Pdiff_catch = (log10(plme_Pmcatch) - log10(Plme_mcatch10)) ./ log10((Plme_mcatch10+eps));
            Ddiff_catch = (log10(plme_Dmcatch) - log10(Dlme_mcatch10)) ./ log10((Dlme_mcatch10+eps));
            
            code = [1:66]';
            T = table(code,slme_mcatch10,plme_mcatch,Flme_mcatch10,plme_Fmcatch,...
                Plme_mcatch10,plme_Pmcatch,Dlme_mcatch10,plme_Dmcatch,'VariableNames',...
                {'lme','saupAll','poemAll','saupF','poemF','saupP','poemP','saupD','poemD'});
            writetable(T,[dpath 'LME_saup_catch_clim_fished',harv,'.csv']);
            
            %% Plots
            l10s=log10(slme_mcatch10+eps);
            l10p=log10(plme_mcatch);
            l10sF=log10(Flme_mcatch10+eps);
            l10pF=log10(plme_Fmcatch);
            l10sP=log10(Plme_mcatch10+eps);
            l10pP=log10(plme_Pmcatch);
            l10sD=log10(Dlme_mcatch10+eps);
            l10pD=log10(plme_Dmcatch);
            
%             figure(1)
%             subplot(2,2,1)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:66
%                 plot(l10s(i),l10p(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
%             end
%             axis([-6 2 -6 2])
%             xlabel('SAUP mean of top 10 years')
%             ylabel('POEM mean of Climatology')
%             title('Mean catch')
%             
%             subplot(2,2,2)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:66
%                 plot(l10sF(i),l10pF(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
%             end
%             axis([-6 2 -6 2])
%             xlabel('SAUP F catch (log10 MT km^-^2)')
%             ylabel('POEM F catch (log10 MT km^-^2)')
%             title('Mean F catch')
%             
%             subplot(2,2,3)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:66
%                 plot(l10sP(i),l10pP(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
%             end
%             axis([-6 2 -6 2])
%             xlabel('SAUP P catch (log10 MT km^-^2)')
%             ylabel('POEM P catch (log10 MT km^-^2)')
%             title('Mean P catch')
%             
%             subplot(2,2,4)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:66
%                 plot(l10sD(i),l10pD(i),'.','MarkerSize',25,'color',tmap(tid(i,2),:)); hold on;
%             end
%             axis([-6 2 -6 2])
%             xlabel('SAUP D catch (log10 MT km^-^2)')
%             ylabel('POEM D catch (log10 MT km^-^2)')
%             title('Mean D catch')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished',harv,'_SAUP_log10catch_comp_types.png'])
%             
%             
%             
            %% Drop Arctic, Antarctic, Hawaii, Australia -------------------------
             
            % Stats
            %r
            rall(fq,pq,dq)=corr(l10s(keep),l10p(keep));
            rF(fq,pq,dq)=corr(l10sF(keep),l10pF(keep));
            rP(fq,pq,dq)=corr(l10sP(keep),l10pP(keep));
            rD(fq,pq,dq)=corr(l10sD(keep),l10pD(keep));
            
            %root mean square error
            o=l10s(keep);
            p=l10p(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmse(fq,pq,dq) = sqrt(num/n);
            
            o=l10sF(keep);
            p=l10pF(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseF(fq,pq,dq) = sqrt(num/n);
            
            o=l10sP(keep);
            p=l10pP(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseP(fq,pq,dq) = sqrt(num/n);
            
            o=l10sD(keep);
            p=l10pD(keep);
            n = length(o);
            num=nansum((p-o).^2);
            rmseD(fq,pq,dq) = sqrt(num/n);
%             
%             figure(11)
%             subplot(2,2,1)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:length(keep)
%                 lme=keep(i);
%                 plot(l10s(lme),l10p(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
%             end
%             text(-3.5,1.5,['r = ' num2str(rall(fq,pq,dq))])
%             text(-3.5,1.0,['RMSE = ' num2str(rmse(fq,pq,dq))])
%             axis([-4 2 -4 2])
%             xlabel('SAUP mean of top 10 years')
%             ylabel('POEM total mean of Climatology')
%             title('Mean catch without Polar and Australia')
%             
%             subplot(2,2,2)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:length(keep)
%                 lme=keep(i);
%                 plot(l10sF(lme),l10pF(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
%             end
%             text(-5.5,1.5,['r = ' num2str(rF(fq,pq,dq))])
%             text(-5.5,1.0,['RMSE = ' num2str(rmseF(fq,pq,dq))])
%             axis([-6 2 -6 2])
%             xlabel('SAUP F catch (log10 MT km^-^2)')
%             ylabel('POEM F catch (log10 MT km^-^2)')
%             title('Mean F catch')
%             
%             subplot(2,2,3)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:length(keep)
%                 lme=keep(i);
%                 plot(l10sP(lme),l10pP(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
%             end
%             text(-5.5,1.5,['r = ' num2str(rP(fq,pq,dq))])
%             text(-5.5,1.0,['RMSE = ' num2str(rmseP(fq,pq,dq))])
%             axis([-6 2 -6 2])
%             xlabel('SAUP P catch (log10 MT km^-^2)')
%             ylabel('POEM P catch (log10 MT km^-^2)')
%             title('Mean P catch')
%             
%             subplot(2,2,4)
%             plot(x,x,'--k'); hold on;
%             plot(x,x2h,'--b'); hold on;
%             plot(x,x2l,'--b'); hold on;
%             plot(x,x5h,'--r'); hold on;
%             plot(x,x5l,'--r'); hold on;
%             for i=1:length(keep)
%                 lme=keep(i);
%                 plot(l10sD(lme),l10pD(lme),'.k','MarkerSize',25,'color',tmap(tid(lme,2),:)); hold on;
%             end
%             text(-5.5,1.5,['r = ' num2str(rD(fq,pq,dq))])
%             text(-5.5,1.0,['RMSE = ' num2str(rmseD(fq,pq,dq))])
%             axis([-6 2 -6 2])
%             xlabel('SAUP D catch (log10 MT km^-^2)')
%             ylabel('POEM D catch (log10 MT km^-^2)')
%             title('Mean D catch')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished',harv,'_SAUP_log10catch_comp_types_LELC.png'])
%             
%             
%             %% MAPS
%             
%             plme = NaN*ones(ni,nj);
%             plmeF = NaN*ones(ni,nj);
%             plmeP = NaN*ones(ni,nj);
%             plmeD = NaN*ones(ni,nj);
%             
%             slme = NaN*ones(ni,nj);
%             slmeF = NaN*ones(ni,nj);
%             slmeP = NaN*ones(ni,nj);
%             slmeD = NaN*ones(ni,nj);
%             
%             tlme = lme_mask_onedeg;
%             
%             for L=1:66
%                 lid = find(tlme==L);
%                 
%                 plme(lid) = plme_mcatch(L);
%                 plmeF(lid) = plme_Fmcatch(L);
%                 plmeP(lid) = plme_Pmcatch(L);
%                 plmeD(lid) = plme_Dmcatch(L);
%                 
%                 slme(lid) = slme_mcatch10(L);
%                 slmeF(lid) = Flme_mcatch10(L);
%                 slmeP(lid) = Plme_mcatch10(L);
%                 slmeD(lid) = Dlme_mcatch10(L);
%             end
%             
%             dlme = ((log10(plme)-log10(slme))./log10(slme)) .* 100;
%             dlmeF = ((log10(plmeF)-log10(slmeF))./log10(slmeF)) .* 100;
%             dlmeP = ((log10(plmeP)-log10(slmeP))./log10(slmeP)) .* 100;
%             dlmeD = ((log10(plmeD)-log10(slmeD))./log10(slmeD)) .* 100;
%             
%             %% Catch
%             
%             % all
%             figure(21)
%             axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%                 'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%             surfm(geolat_t,geolon_t,dlme)
%             cmocean('balance')
%             load coast;                     %decent looking coastlines
%             h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%             caxis([-50 50]);
%             hcb = colorbar('h');
%             ylim(hcb,[-50 50])                   %Set color axis if needed
%             set(gcf,'renderer','painters')
%             title('Climatology difference from SAUP total annual catch (MT)')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished_LME_SAUP_catch_diff.png'])
%             
%             %% F
%             figure(22)
%             axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%                 'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%             surfm(geolat_t,geolon_t,dlmeF)
%             cmocean('balance')
%             load coast;                     %decent looking coastlines
%             h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%             caxis([-50 50]);
%             hcb = colorbar('h');
%             ylim(hcb,[-50 50])                   %Set color axis if needed
%             set(gcf,'renderer','painters')
%             title('Climatology difference from SAUP total annual F catch (MT)')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished_LME_SAUP_catch_diffF.png'])
%             
%             % all
%             figure(23)
%             axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%                 'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%             surfm(geolat_t,geolon_t,dlmeP)
%             cmocean('balance')
%             load coast;                     %decent looking coastlines
%             h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%             caxis([-50 50]);
%             hcb = colorbar('h');
%             ylim(hcb,[-50 50])                   %Set color axis if needed
%             set(gcf,'renderer','painters')
%             title('Climatology difference from SAUP total annual P catch (MT)')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished_LME_SAUP_catch_diffP.png'])
%             
%             % all
%             figure(24)
%             axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%                 'Grid','off','FLineWidth',1,'origin',[0 -100 0])
%             surfm(geolat_t,geolon_t,dlmeD)
%             cmocean('balance')
%             load coast;                     %decent looking coastlines
%             h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);
%             caxis([-50 50]);
%             hcb = colorbar('h');
%             ylim(hcb,[-50 50])                   %Set color axis if needed
%             set(gcf,'renderer','painters')
%             title('Climatology difference from SAUP total annual D catch (MT)')
%             stamp(cfile)
%             print('-dpng',[ppath 'Clim_fished_LME_SAUP_catch_diffD.png'])
%             
        end
    end
end

save([dp 'Clim_fished',harv,'_LME_SAUP_catch_comp_qs.mat'],'rall','rF','rP',...
    'rD','rmse','rmseF','rmseP','rmseD')

%% Plots of r and RMSE
cfile2 = ['_cm20_m-b175-k09_D075_J100_A050_Sm025_nmort1_BE05_CC100_RE00100'];

fqs2 = 0.5:0.5:4;
pqs2 = 0.25:0.25:1.75;
dqs2 = 1:7;
[pgrid1,fgrid1]=meshgrid(pqs2,fqs2);
[dgrid1,fgrid2]=meshgrid(dqs2,fqs2);
[dgrid2,pgrid2]=meshgrid(dqs2,pqs2);
nf = length(fqs);
np = length(pqs);
nd = length(dqs);

Rall= NaN*ones(nf+1,np+1,nd+1);
RF  = NaN*ones(nf+1,np+1,nd+1);
RP  = NaN*ones(nf+1,np+1,nd+1);
RD  = NaN*ones(nf+1,np+1,nd+1);
Rall(1:nf,1:np,1:nd)=rall;
RF(1:nf,1:np,1:nd)=rF;
RP(1:nf,1:np,1:nd)=rP;
RD(1:nf,1:np,1:nd)=rD;

for j=1:nd
    figure(2)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,Rall(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('Corr all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_FP'])
    
    figure(3)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RF(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('Corr F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_FP'])
    
    figure(4)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RP(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('Corr P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_FP'])
    
    figure(5)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RD(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('Corr D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_FP'])
end
%%
for j=1:np
    figure(6)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(Rall(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_FD'])
    
    figure(7)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RF(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_FD'])
    
    figure(8)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RP(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_FD'])
    
    figure(9)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RD(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('Corr D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_FD'])
end
for j=1:nf
    figure(10)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(Rall(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrAll_frates_PD'])
    
    figure(11)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RF(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrF_frates_PD'])
    
    figure(12)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RP(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrP_frates_PD'])
    
    figure(13)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RD(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('Corr D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_corrD_frates_PD'])
    
end

%%
Rmseall  = NaN*ones(nf+1,np+1,nd+1);
RmseF  = NaN*ones(nf+1,np+1,nd+1);
RmseP  = NaN*ones(nf+1,np+1,nd+1);
RmseD  = NaN*ones(nf+1,np+1,nd+1);
Rmseall(1:nf,1:np,1:nd)=rmse;
RmseF(1:nf,1:np,1:nd)=rmseF;
RmseP(1:nf,1:np,1:nd)=rmseP;
RmseD(1:nf,1:np,1:nd)=rmseD;

for j=1:nd
    figure(14)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,Rmseall(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('RMSE all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_FP'])
    
    figure(15)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RmseF(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('RMSE F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_FP'])
    
    figure(16)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RmseP(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('RMSE P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_FP'])
    
    figure(17)
    subplot(2,3,j)
    pcolor(fgrid1,pgrid1,RmseD(:,:,j))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('P sel')
    title('RMSE D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_FP'])
end
for j=1:np
    figure(18)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(Rmseall(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_FD'])
    
    figure(19)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseF(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_FD'])
    
    figure(20)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseP(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_FD'])
    
    figure(21)
    subplot(2,3,j)
    pcolor(fgrid2,dgrid1,squeeze(RmseD(:,j,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('F sel')
    ylabel('D sel')
    title('RMSE D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_FD'])
end
for j=1:nf
    figure(22)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(Rmseall(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE all fish')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseAll_frates_PD'])
    
    figure(23)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseF(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE F')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseF_frates_PD'])
    
    figure(24)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseP(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE P')
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseP_frates_PD'])
    
    figure(25)
    subplot(3,3,j)
    pcolor(pgrid2,dgrid2,squeeze(RmseD(j,:,:)))
    colorbar
    caxis([0 1])
    %set(gca,'XTickLabel',bees(1:2:end),'YTickLabel',kays(2:2:end))
    xlabel('P sel')
    ylabel('D sel')
    title('RMSE D')
    stamp(cfile2)
    %print('-dpng',[pp 'Clim_',harv,'_LME_SAUP_catch_rmseD_frates_PD'])
    
end

%%

corr_all = normalize(rall)+normalize(rF)+normalize(rP)+normalize(rD);
max(corr_all(:)) %k=0.11;b=0.25

rmse_all = normalize(rmse)+normalize(rmseF)+normalize(rmseP)+normalize(rmseD);
min(rmse_all(:)) %k=0.11;b=0.225

corr_PD = normalize(rP)+normalize(rD);
max(corr_PD(:)) %k=0.11;b=0.25

rmse_PD = normalize(rmseP)+normalize(rmseD);
min(rmse_PD(:)) %k=0.11;b=0.225
