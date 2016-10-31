% Calculate different skill metrics for each bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_mizernmort_BE05_RE01';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'Data_spinup_pristine_' cfile '.mat'],'MD','LD','BENT');
seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei2010_Global_seafloor_biomass.csv',1,0);

%% Wei benthic biomass
Wcol = {'latitude','longitude','depth','bact.biom.mean','meio.biom.mean',...
    'macro.biom.mean','mega.biom.mean','inv.biom.mean','fis.biom.mean'};
Wcol = Wcol';
% all mean biomasses in log10 mg C/m2
invert = seafl(:,8);
fish = seafl(:,9);
% convert to g WW/m2
invert = 10.^(invert) * 1e-3 * 9.0;
fish = 10.^(fish) * 1e-3 * 9.0;

%%
[loc,days]=size(LD.bio);
x=1:days;

lyr=x((end-365+1):end);

md_mean=mean(MD.bio(:,lyr),2);
ld_mean=mean(LD.bio(:,lyr),2);
b_mean=mean(BENT.bio(:,lyr),2);

% Plots in space
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);

Zmd=griddata(grid(:,2),grid(:,3),md_mean,X,Y);
Zld=griddata(grid(:,2),grid(:,3),ld_mean,X,Y);
Zb=griddata(grid(:,2),grid(:,3),b_mean,X,Y);
Zi=griddata(seafl(:,2),seafl(:,1),invert,X,Y);
Zf=griddata(seafl(:,2),seafl(:,1),fish,X,Y);

Zd = Zmd+Zld;

%% maps 
% Bent
figure(1)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zb))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean model benthic biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_BENT.png'])

% D
figure(2)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zd))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean model M&L D biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Spinup_global_Demersal.png'])

% Inv
figure(3)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zi))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Wei benthic invert biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Wei_inverts.png'])

% D
figure(4)
m_proj('miller','lat',82);
m_pcolor(X,Y,real(log10(Zf))); hold on;
shading flat
m_coast('patch',[.5 .5 .5],'edgecolor','none');
m_grid;
title('log10 mean Wei benthic fish biomass (g m^-^2)')
colormap('jet')
colorbar('h')
caxis([-3 1])
stamp(cfile)
print('-dpng',[ppath 'Wei_fish.png'])

%% Skill metrics 

metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

% REDO Use log biomass

obs(:,1) = real(log(Zi(:)));
obs(:,2) = real(log(Zf(:)));
model(:,1) = real(log(Zb(:)));
model(:,2) = real(log(Zd(:)));

n = length(Zi(:));

skill=NaN*ones(14,2);
for j=1:2
    o=obs(:,j);
    p=model(:,j);
    o(o==0)=1e-6;
    p(p==0)=1e-6;
    
    omean=repmat(nanmean(o),n,1);
    pmean=repmat(nanmean(p),n,1);
    osig=nanstd(o);
    psig=nanstd(p);
    
    % corr coeff
    num=nansum((o-omean).*(p-pmean));
    d1=nansum((o-omean).^2);
    d2=nansum((p-pmean).^2);
    den=sqrt(d1*d2);
    skill(1,j) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    skill(2,j) = sqrt(num/n);
    
    % reliability index  %%What about predicted zero values?
    q1=nansum((log(o./p)).^2);
    skill(3,j) = exp(sqrt(q1/n));
    
    % average error
    skill(4,j) = nansum(p-o) / n;
    
    % average absolute error
    skill(5,j) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    skill(6,j) = (num1-num2)/num1;
    
    % Taylor R
    num=nansum((o-omean).*(p-pmean));
    skill(7,j) = num/(n*osig*psig);
    
    % Taylor normalized std
    skill(8,j) = psig/osig;
    
    % unbiased root mean square difference
    % sign tells difference between model and obs std
    q1=nansum(((p-pmean)-(o-omean)).^2);
    if (psig>=osig)
        skill(9,j) = sqrt(q1/n);
    else
        skill(9,j) = -1*sqrt(q1/n);
    end
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    skill(10,j) = sqrt(q1/n);
    
    % normalized bias
    skill(11,j) = (pmean(1,1)-omean(1,1));%./osig;
    
    % Joliff et al 2009 S1
    R=skill(7,j);
    s=skill(8,j);
    num=2*(1+R);
    den=(s + (1/s)).^2;
    skill(12,j) = 1 - (num/den);
    
    % Joliff et al 2009 S2
    num=(1+R).^4;
    den=4*(s + (1/s)).^2;
    skill(13,j) = 1 - (num/den);
    
    % Joliff et al 2009 S3
    q1=exp(-((s-1).^2) / 0.18);
    q2=(1+R)/2;
    skill(14,j) = 1 - (q1*q2);
    
end

%% Plot results

simtex = {'inverts','fish'};

% Bar graphs
figure(5)
subplot(2,2,1)
bar(skill(1,:),'k')
title('Correlation coefficient')
set(gca,'XTickLabel',simtex)
stamp(cfile)
%print('-depsc2',[ppath 'corr_coeff_Wei'])

subplot(2,2,2)
bar(skill(2,:),'k')
title('Root mean square error')
set(gca,'XTickLabel',simtex)
%print('-depsc2',[ppath 'RMSE_Wei'])

subplot(2,2,3)
bar(skill(6,:),'k')
title('Modeling efficiency')
set(gca,'XTickLabel',simtex)
%print('-depsc2',[ppath 'MEF_Wei'])
print('-dpng',[ppath 'corr_rmse_mef_Wei'])


%% Taylor diagram using Joliff corr coeff
% cm={[1 0.5 0],...   %orange
%     [0.5 0.5 0],... %tan/army
%     [0 0.7 0],...   %g
%     [0 1 1],...     %c
%     [0 0 0.75],...  %b
%     [0.5 0 1],...   %purple
%     [1 0 1],...     %m
%     [1 0 0],...     %r
%     [0.5 0 0],...   %maroon
%     [0.35 0.35 0.35]}; %grey

cm={'r','b'}; 

[rmsd,it]=sort(skill(10,:),'descend');
theta=acos(skill(7,:));
rho=skill(8,:);
simtex{3}='obs';

%%
tr=0;
rr=1;
figure(6)
h0=polar(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:2
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 2 0 2.1])
title('Joliff Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[ppath 'Taylor_Joliff_Wei'])

%% Taylor diagram using corr coeff calculated
theta2=acos(skill(1,:));

figure(7)
h0=polar(0,2,'.'); hold on;        %0.05
set(h0,'color','w');
for s=1:2
    h=polar(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
end
%set(h,'clim',[0 1.5]);
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 2 0 2.1])
title('Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[ppath 'Taylor_Wei'])


%% Best
one=[1;3;6;7;8];
zer=[2;4;5;9;10;11];
high=[12;13;14];
for z=1:14
    Z1=sum(z==one);
    Z2=sum(z==zer);
    Z3=sum(z==high);
    if (Z1>0)
        best{z,1}=simtex{find(min(abs(1-(skill(z,:))))==abs(1-(skill(z,:))))};
    elseif (Z2>0)
        best{z,1}=simtex{find(min(abs(0-(skill(z,:))))==abs(0-(skill(z,:))))};
    elseif (Z3>0)
        best{z,1}=simtex{find(max(skill(z,:))==skill(z,:))};
    end
end

T=table(metrics,best);

save([dpath 'skill_Wei.mat'],'model','obs','metrics','skill');

%% Save Wei biomasses at locations
% glon = seafl(:,2);
% glat = seafl(:,1);
% 
% % Georges Bank (NEP)
% lon=find(glon<=-66 & glon>=-67);
% lat=find(glat<=42 & glat>=41);
% gid=intersect(lon,lat);
% 
% % Eastern Bering Sea (Pribs)
% lon=find(glon<=-169 & glon>=-170);
% lat=find(glat<=57 & glat>=56);
% eid=intersect(lon,lat);
% 
% % Subarctic Pacific Gyre (Ocean Station Papa)
% lon=find(glon<=-145 & glon>=-146);
% lat=find(glat<=51 & glat>=50);
% pid=intersect(lon,lat);
% 
% % HOT
% lon=find(glon<=-157 & glon>=-158);
% lat=find(glat<=23 & glat>=22);
% hid=intersect(lon,lat);
% 
% % BATS
% lon=find(glon<=-64 & glon>=-65);
% lat=find(glat<=32 & glat>=31);
% bid=intersect(lon,lat);
% 
% % North Sea
% lon=find(glon<=4 & glon>=3);
% lat=find(glat<=57 & glat>=56);
% nid=intersect(lon,lat);
% 
% % E Eq Pac
% lon=find(glon<=-110 & glon>=-111);
% lat=find(glat<=5.5 & glat>=5);
% qid=intersect(lon,lat);
% 
% % Subpolar W Pac station K2: 47oN, 160oE
% lon=find(glon<=161 & glon>=160);
% lat=find(glat<=48 & glat>=47);
% kid=intersect(lon,lat);
% 
% % Subtropical W Pac station S1: 30oN, 145oE
% lon=find(glon<=146 & glon>=145);
% lat=find(glat<=31 & glat>=30);
% sid=intersect(lon,lat);
% 
% names={'Georges Bank','Eastern Bering Sea','Ocean Station Papa',...
%     'HOT','BATS','North Sea','Eastern Equatorial Pacific','K2','S1'};
% ids=[gid;eid;pid;hid;bid;nid;qid;kid;sid];
% B = invert(ids);
% F = fish(ids);
% 
% %
% T=table(names',B,F,...
%     'VariableNames',{'Location','Inverts','Fish'});
% writetable(T,'Wei_inverts_fish_gWWm2_locs.csv','Delimiter',',');
% save('Wei_inverts_fish_gWWm2_locs.mat','T');

