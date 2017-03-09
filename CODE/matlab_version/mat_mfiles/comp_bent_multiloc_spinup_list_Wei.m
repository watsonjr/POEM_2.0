% Calculate different skill metrics for each oneloc bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Bent_CC/';
datap = '/Volumes/GFDL/CSV/Matlab_big_size/';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};

sname = 'Spinup_';

np1 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE1000';
np2 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE0500';
np3 = 'Dc_TrefO_Hold_cmax-metab_MFeqMP_fcrit40_D100_nmort2_BE05_CC050_RE0100';
np4 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0500';
np5 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0100';
np6 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0050';
np7 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC050_RE0010';
np8 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE10_CC050_RE0100';
np9 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE10_CC050_RE0050';
np10 = 'Dc_TrefO_cmax-metab4_enc4_MFeqMP_fcrit40_D100_nmort3_BE05_CC100_RE0100';
dp = {np1;np2;np3;np4;np5;np6;np7;np8;np9;np10};

sims = {'212-05-05-10';'212-05-05-05';'212-05-05-01';'443-05-05-05';'443-05-05-01';'443-05-05-005';...
    '443-05-05-001';'443-10-05-01';'443-10-05-005';'443-05-10-01'};

%% Wei benthic biomass
seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';

% all mean biomasses in log10 mg C/m2
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

%% put on same grid as POEM output
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/hindcast_gridspec.mat',...
    'geolon_t','geolat_t','tmask');
grid = csvread([cpath 'grid_csv.csv']);

mask = tmask(:,:,1);

test=seafl(:,2);
id=find(test>80);
test(id)=test(id)-360;

Zi = griddata(test,seafl(:,1),invert,geolon_t,geolat_t);
Zf = griddata(test,seafl(:,1),fish,geolon_t,geolat_t);
Zi = Zi.*mask;
Zf = Zf.*mask;

[ni,nj]=size(geolon_t);

%% ocean cells
% ocean=zeros(ni,nj);
% ocean(grid(:,1))=ones(size(grid(:,1)));
% figure(55)
% surf(geolon_t,geolat_t,ocean); view(2); hold on;

%% Maps
% Inv
figure(1)
surf(geolon_t,geolat_t,log10(Zi)); view(2); hold on;
shading flat
title('log10 mean benthic invert biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_inverts.png')


% Fish
figure(2)
surf(geolon_t,geolat_t,log10(Zf)); view(2); hold on;
shading flat
title('log10 mean benthic fish biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-2.5 0.5])
print('-dpng','/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Wei_fish.png')

%% colors
cm={[1 0.5 0],...   %orange
    [0.5 0.5 0],... %tan/army
    [0 0.7 0],...   %g
    [0 1 1],...     %c
    [0 0 0.75],...  %b
    [0.5 0 1],...   %purple
    [1 0 1],...     %m
    [1 0 0],...     %r
    [0.5 0 0],...   %maroon
    [0.75 0.75 0.75],... %lt grey
    [0.5 0.5 0.5],...    %med grey
    [49/255 79/255 79/255],... %dk grey
    [0 0 0],...      %black
    [1 1 0],...      %yellow
    [127/255 255/255 0],... %lime green
    [0 0.5 0],...    %dk green
    [0/255 206/255 209/255],... %turq
    [0 0.5 0.75],...   %med blue
    [188/255 143/255 143/255],... %rosy brown
    [255/255 192/255 203/255],... %pink
    [255/255 160/255 122/255]}; %peach

%% Skill metrics

% Stowe et al 2009
% Univariate (model vs observ)
% 1. correlation coefficient [-1 1]; want close to 1
% 2. root mean square error; want close to 0; considering the magnitude rather than the direction
% 3. average error (bias); want close to 0; neg & pos discrepancies can cancel each other
% 4. average absolute error; want close to 0; considering the magnitude rather than the direction
% 5. modeling efficiency; want close to 1; <0 means avg obs better than model; how well a model predicts relative to avg obs

metrics={'r','RMSE','AE','AAE','MEF'};
metrics=metrics';

%Sometimes it is appropriate to log-transform the observations and predictions before
%calculating goodness-of-fit statistics so that differences between predicted and observed
%values will not be highly skewed and dominated by a small proportion of high values.

r = NaN*ones(length(dp),1);
RMSE = r;
AE = r;
AAE = r;
MEF = r;
NSTD = r;
RMSD = r;

%%
iskill=NaN*ones(7,length(dp));
    fskill=NaN*ones(7,length(dp));
    ibest=cell(5,1);
    fbest=cell(5,1);
    
for i=1:length(dp)
    close all
    
    %%
    cfile = dp{i};
    fpath=['/Volumes/GFDL/NC/Matlab_big_size/' cfile '/'];
    ppath = [pp cfile];
    load([fpath 'Means_spinup_' cfile '.mat'],'md_mean','ld_mean','b_mean');
    
    %% Model bent & dem
    Zmd=NaN*ones(ni,nj);
    Zld=NaN*ones(ni,nj);
    Zb=NaN*ones(ni,nj);
    
    Zmd(grid(:,1))=md_mean;
    Zld(grid(:,1))=ld_mean;
    Zb(grid(:,1))=b_mean;
    
    Zd = Zmd+Zld;
    
    % Bent
    figure(3)
    surf(geolon_t,geolat_t,log10(Zb)); view(2); hold on;
    shading flat
    title('log10 mean benthic biomass (g m^-^2)')
    colormap('jet')
    colorbar('h')
    caxis([-2.5 0.5])
    stamp(cfile)
    print('-dpng',[ppath '_Spinup_global_Bent.png'])
    
    % D
    figure(4)
    surf(geolon_t,geolat_t,log10(Zd)); view(2); hold on;
    shading flat
    title('log10 mean model M&L D biomass (g m^-^2)')
    colormap('jet')
    colorbar('h')
    caxis([-2.5 0.5])
    stamp(cfile)
    print('-dpng',[ppath '_Spinup_global_Demersal.png'])
    
    
    %% Inverts
    obs = real(log(Zi(grid(:,1))));
    model = real(log(b_mean));
    
    n = length(obs);
    
    o=obs;
    p=model;
    %             o(isnan(o))=0;
    %             p(isnan(p))=0;
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    iskill(1,i) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    iskill(2,i) = sqrt(num/n);
    
    % average error
    iskill(3,i) = nansum(p-o) / n;
    
    % average absolute error
    iskill(4,i) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    iskill(5,i) = (num1-num2)/num1;
    
    % Taylor normalized std
    iskill(6,i) = psig/osig;
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    iskill(7,i) = sqrt(q1/n);
    
    
    %% Fish
    obs = real(log(Zf(grid(:,1))));
    model = real(log(md_mean+ld_mean));
    
    n = length(obs);
    
    o=obs;
    p=model;
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    fskill(1,i) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    fskill(2,i) = sqrt(num/n);
    
    % average error
    fskill(3,i) = nansum(p-o) / n;
    
    % average absolute error
    fskill(4,i) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    fskill(5,i) = (num1-num2)/num1;
    
    % Taylor normalized std
    fskill(6,i) = psig/osig;
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    fskill(7,i) = sqrt(q1/n);
    
end %sims

%% Plot results
cfile2 = 'Dc_TrefO_cmax-metab_enc_MFeqMP_fcrit40_D100_BentCCtests';

% Bar graphs
figure(5)
subplot(2,2,1)
bar(iskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([-0.1 1])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.15,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
stamp(cfile2)
title('Inverts skill','HorizontalAlignment','left')

subplot(2,2,2)
bar(iskill(2,:),'k')
ylabel('Root mean square error')
ylim([0 2])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,3)
bar(iskill(3,:),'k')
ylabel('Average error')
ylim([-0.2 1.2])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.25,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,4)
bar(iskill(5,:),'k')
ylabel('Modeling efficiency')
ylim([-1.5 0.5])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-1.55,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
print('-dpng',[pp cfile2 '_Spinup_skill_Wei_inverts'])

%%
figure(6)
subplot(2,2,1)
bar(fskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([-0.1 1])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.15,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
stamp(cfile2)
title('Fish skill','HorizontalAlignment','left')

subplot(2,2,2)
bar(fskill(2,:),'k')
ylabel('Root mean square error')
ylim([0 3.5])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,3)
bar(fskill(3,:),'k')
ylabel('Average error')
ylim([0 3])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-0.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end

subplot(2,2,4)
bar(fskill(5,:),'k')
ylabel('Modeling efficiency')
ylim([-5 0.5])
xlim([0 length(dp)+1])
set(gca,'XTick',1:length(dp),'XTickLabel',[]);
for t=1:length(dp)
    text(t,-5.05,num2str(sims{t}),'Rotation',45,'HorizontalAlignment','right')
end
print('-dpng',[pp cfile2 '_Spinup_skill_Wei_fish'])

%% Taylor diagrams
% INVERTS
[rmsd,it]=sort(iskill(7,:),'descend');
theta=acos(iskill(1,it));    %corr coeff
rho=iskill(6,it);            %stdev
simtex=sims(it);
simtex{length(sims)+1}='obs';

% % Just those <2
% id = find(rho<2);
% theta = theta(id);
% rho = rho(id);
% simtex=sims(id);
% simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(7)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:length(dp)
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Inverts Taylor diagram')
legend([' '; simtex])
legend('location','northeast')
print('-dpng',[pp cfile2 '_Spinup_global_invert_Taylor_Wei'])

%% FISH

[rmsd,it]=sort(fskill(7,:),'descend');
theta=acos(fskill(1,it));    %corr coeff
rho=fskill(6,it);            %stdev
simtex=sims(it);
simtex{length(sims)+1}='obs';

% Just those <1.5
% id = find(rho<1.5);
% theta = theta(id);
% rho = rho(id);
% simtex=sims(id);
% simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(8)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:length(dp)
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Fish Taylor diagram')
legend([' '; simtex])
legend('location','northeast')
print('-dpng',[pp cfile2 '_Spinup_global_fish_Taylor_Wei'])

%% Best
one=[1;5];
zer=[2;3;4];
for z=1:5
    Z1=sum(z==one);
    Z2=sum(z==zer);
    if (Z1>0)
        ibest{z,1}=sims(find(min(abs(1-(iskill(z,:))))==abs(1-(iskill(z,:)))));
        fbest{z,1}=sims(find(min(abs(1-(fskill(z,:))))==abs(1-(fskill(z,:)))));
    elseif (Z2>0)
        ibest{z,1}=sims(find(min(abs(0-(iskill(z,:))))==abs(0-(iskill(z,:)))));
        fbest{z,1}=sims(find(min(abs(0-(fskill(z,:))))==abs(0-(fskill(z,:)))));
    end
end

Ti=table(metrics,ibest);
Tf=table(metrics,fbest);

%%
save([datap 'Bent_CC_tests/' cfile2 '_global_skill_Wei.mat'],...
    'invert','fish','metrics','iskill','fskill','Ti','Tf');



