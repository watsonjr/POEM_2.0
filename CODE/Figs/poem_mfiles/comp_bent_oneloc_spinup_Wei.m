% Calculate different skill metrics for each oneloc bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
datap = '/Volumes/GFDL/CSV/';

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

npath1 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/'];
npath2 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP25/'];
npath3 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP50/'];
npath4 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_BP75/'];
npath5 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/'];
npath6 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/'];
npath7 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP25/'];
npath8 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP50/'];
npath9 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10_BP75/'];
npath10 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/'];
npath11 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP25/'];
npath12 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP50/'];
npath13 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_BP75/'];
npath14 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/'];
npath15 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP25/'];
npath16 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP50/'];
npath17 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10_BP75/'];
npath18 = [datap 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/'];

dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;npath11;npath12;...
    npath13;npath14;npath15;npath16;npath17;npath18};
sims = {'30-025-1','30-05-1/4','30-05-1/2','30-05-3/4','30-05-1','30-075-1','30-10-1/4',...
    '30-10-1/2','30-10-3/4','30-10-1','40-05-1/4','40-05-1/2','40-05-3/4','40-05-1','40-10-1/4',...
    '40-10-1/2','40-10-3/4','40-10-1'};
%cfile = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_BentPref_comp';

seafl = csvread('/Users/cpetrik/Dropbox/Princeton/POEM_other/Wei_inverts_fish_gWWm2_locs.csv',1,1);

sname = 'Spinup_';
sname2 = '';

%% Wei benthic biomass
spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
% in g WW/m2
invert = seafl(:,1);
fish = seafl(:,2);

%% Model bent & dem
ndp = length(dp);
modi = NaN*ones(length(spots),ndp);
modf = NaN*ones(length(spots),ndp);
for s=1:ndp
    dpath = char(dp(s));
    load([dpath sname sname2 'lastyr_bent_biom.mat'],'Bmean');
    load([dpath sname sname2 'lastyr_sum_mean_biom.mat'],'all_mean');
    
    %Demersal fishes (M&L)
    dem = nansum(squeeze(all_mean(2:3,3,:)));
    
    modi(:,s) = Bmean;
    modf(:,s) = dem';
    
end


%% Skill metrics

% Stowe et al 2009
% Univariate (model vs observ)
% 1. correlation coefficient [-1 1]; want close to 1
% 2. root mean square error; want close to 0; considering the magnitude rather than the direction
% 3. reliability index; want close to 1; avg factor model differ from obs; e.g. RI=2.0 is model predicts obs within a multiplicative factor of 2
% 4. average error (bias); want close to 0; neg & pos discrepancies can cancel each other
% 5. average absolute error; want close to 0; considering the magnitude rather than the direction
% 6. modeling efficiency; want close to 1; <0 means avg obs better than model; how well a model predicts relative to avg obs

metrics={'r','RMSE','RI','AE','AAE','MEF','R','norm std','unb RMSD',...
    'tot RMSD','bias','S1','S2','S3'};
metrics=metrics';

%Sometimes it is appropriate to log-transform the observations and predictions before
%calculating goodness-of-fit statistics so that differences between predicted and observed
%values will not be highly skewed and dominated by a small proportion of high values.

%% Inverts
obs = real(log(invert));
model = real(log(modi));

n = length(obs);

iskill=NaN*ones(14,2);
for j=1:ndp
    o=obs;
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
    iskill(1,j) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    iskill(2,j) = sqrt(num/n);
    
    % reliability index  %%What about predicted zero values?
    q1=nansum((log(o./p)).^2);
    iskill(3,j) = exp(sqrt(q1/n));
    
    % average error
    iskill(4,j) = nansum(p-o) / n;
    
    % average absolute error
    iskill(5,j) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    iskill(6,j) = (num1-num2)/num1;
    
    % Taylor R
    num=nansum((o-omean).*(p-pmean));
    iskill(7,j) = num/(n*osig*psig);
    
    % Taylor normalized std
    iskill(8,j) = psig/osig;
    
    % unbiased root mean square difference
    % sign tells difference between model and obs std
    q1=nansum(((p-pmean)-(o-omean)).^2);
    if (psig>=osig)
        iskill(9,j) = sqrt(q1/n);
    else
        iskill(9,j) = -1*sqrt(q1/n);
    end
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    iskill(10,j) = sqrt(q1/n);
    
    % normalized bias
    iskill(11,j) = (pmean(1,1)-omean(1,1));%./osig;
    
    % Joliff et al 2009 S1
    R=iskill(7,j);
    s=iskill(8,j);
    num=2*(1+R);
    den=(s + (1/s)).^2;
    iskill(12,j) = 1 - (num/den);
    
    % Joliff et al 2009 S2
    num=(1+R).^4;
    den=4*(s + (1/s)).^2;
    iskill(13,j) = 1 - (num/den);
    
    % Joliff et al 2009 S3
    q1=exp(-((s-1).^2) / 0.18);
    q2=(1+R)/2;
    iskill(14,j) = 1 - (q1*q2);
    
end

%% Fish
obs = real(log(fish));
model = real(log(modf));

n = length(obs);

fskill=NaN*ones(14,2);
for j=1:ndp
    o=obs;
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
    fskill(1,j) = num/den;
    
    % root mean square error
    num=nansum((p-o).^2);
    fskill(2,j) = sqrt(num/n);
    
    % reliability index  %%What about predicted zero values?
    q1=nansum((log(o./p)).^2);
    fskill(3,j) = exp(sqrt(q1/n));
    
    % average error
    fskill(4,j) = nansum(p-o) / n;
    
    % average absolute error
    fskill(5,j) = nansum(abs(p-o)) / n;
    
    % modeling efficiency
    num1=nansum((o-omean).^2);
    num2=nansum((p-o).^2);
    fskill(6,j) = (num1-num2)/num1;
    
    % Taylor R
    num=nansum((o-omean).*(p-pmean));
    fskill(7,j) = num/(n*osig*psig);
    
    % Taylor normalized std
    fskill(8,j) = psig/osig;
    
    % unbiased root mean square difference
    % sign tells difference between model and obs std
    q1=nansum(((p-pmean)-(o-omean)).^2);
    if (psig>=osig)
        fskill(9,j) = sqrt(q1/n);
    else
        fskill(9,j) = -1*sqrt(q1/n);
    end
    
    % total root mean square difference
    q1=nansum((o-p).^2);
    fskill(10,j) = sqrt(q1/n);
    
    % normalized bias
    fskill(11,j) = (pmean(1,1)-omean(1,1));%./osig;
    
    % Joliff et al 2009 S1
    R=fskill(7,j);
    s=fskill(8,j);
    num=2*(1+R);
    den=(s + (1/s)).^2;
    fskill(12,j) = 1 - (num/den);
    
    % Joliff et al 2009 S2
    num=(1+R).^4;
    den=4*(s + (1/s)).^2;
    fskill(13,j) = 1 - (num/den);
    
    % Joliff et al 2009 S3
    q1=exp(-((s-1).^2) / 0.18);
    q2=(1+R)/2;
    fskill(14,j) = 1 - (q1*q2);
    
end

%% Plot results

% Bar graphs
figure(4)
subplot(3,1,1)
bar(iskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([-0.25 0.75])
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.3,sims{t},'HorizontalAlignment','right','Rotation',45)
end
title('Inverts fcrit-BEff-Bpref')

subplot(3,1,2)
bar(iskill(2,:),'k')
ylabel('Root mean square error')
set(gca,'XTick',1:ndp,'XTickLabel',[]);
ylim([0 2.5])
for t=1:ndp
    text(t,-0.1,sims{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,1,3)
bar(iskill(6,:),'k')
ylabel('Modeling efficiency')
ylim([-1.5 0.5])
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-1.6,sims{t},'HorizontalAlignment','right','Rotation',45)
end
print('-dpng',[fpath sname 'locs_corr_rmse_mef_Wei_inverts'])

figure(5)
subplot(3,1,1)
bar(fskill(1,:),'k')
ylabel('Correlation coefficient')
ylim([0 1])
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.1,sims{t},'HorizontalAlignment','right','Rotation',45)
end
title('Fish fcrit-BEff-Bpref')

subplot(3,1,2)
bar(fskill(2,:),'k')
ylabel('Root mean square error')
set(gca,'XTick',1:ndp,'XTickLabel',[]);
ylim([0 2.5])
for t=1:ndp
    text(t,-0.1,sims{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,1,3)
bar(fskill(6,:),'k')
ylabel('Modeling efficiency')
ylim([-1 1])
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-1.1,sims{t},'HorizontalAlignment','right','Rotation',45)
end
print('-dpng',[fpath sname 'locs_corr_rmse_mef_Wei_fish'])


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


%% INVERTS
% Taylor diagram using Joliff corr coeff
[rmsd,it]=sort(iskill(10,:),'descend');
theta=acos(iskill(7,:));
rho=iskill(8,:);

% Just those <2

id = find(rho<2);
theta = theta(id);
rho = rho(id);
simtex=sims(id);
simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(6)
h0=polar(0,2,'.'); hold on;
set(h0,'color','w');
for s=1:length(id)
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 2 0 2])
title('Inverts Jolliff Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[fpath sname 'locs_invert_Taylor_Joliff_Wei2'])

% Taylor diagram using corr coeff calculated
theta2=acos(iskill(1,:));

% Just those <2
id = find(rho<2);
theta2 = theta2(id);
rho = rho(id);
simtex=sims(id);
simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(7)
h0=polar(0,2,'.'); hold on;
set(h0,'color','w');
for s=1:length(id)
    h=polar(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 2 0 2])
title('Inverts Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[fpath sname 'locs_invert_Taylor_Wei2'])


%% FISH
% Taylor diagram using Joliff corr coeff

[rmsd,it]=sort(fskill(10,:),'descend');
theta=acos(fskill(7,:));
rho=fskill(8,:);

%% Just those <1.5

id = find(rho<1.5);
theta = theta(id);
rho = rho(id);
simtex=sims(id);
simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(8)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:length(id)
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Fish Jolliff Taylor diagram')
legend([' ' simtex])
legend('location','westoutside')
print('-dpng',[fpath sname 'locs_fish_Taylor_Joliff_Wei2'])

%% Taylor diagram using corr coeff calculated
theta2=acos(fskill(1,:));

% Just those <1.5
id = find(rho<1.5);
theta2 = theta2(id);
rho = rho(id);
simtex=sims(id);
simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(9)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:length(id)
    h=polar(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Fish Taylor diagram')
legend([' ' simtex])
legend('location','westoutside')
print('-dpng',[fpath sname 'locs_fish_Taylor_Wei2'])


%% Best
one=[1;3;6;7;8];
zer=[2;4;5;9;10;11];
high=[12;13;14];
for z=1:14
    Z1=sum(z==one);
    Z2=sum(z==zer);
    Z3=sum(z==high);
    if (Z1>0)
        ibest{z,1}=sims{find(min(abs(1-(iskill(z,:))))==abs(1-(iskill(z,:))))};
        fbest{z,1}=sims{find(min(abs(1-(fskill(z,:))))==abs(1-(fskill(z,:))))};
    elseif (Z2>0)
        ibest{z,1}=sims{find(min(abs(0-(iskill(z,:))))==abs(0-(iskill(z,:))))};
        fbest{z,1}=sims{find(min(abs(0-(fskill(z,:))))==abs(0-(fskill(z,:))))};
    elseif (Z3>0)
        ibest{z,1}=sims{find(max(iskill(z,:))==iskill(z,:))};
        fbest{z,1}=sims{find(max(fskill(z,:))==fskill(z,:))};
    end
end

Ti=table(metrics,ibest);
Tf=table(metrics,fbest);

save([fpath sname 'locs_skill_Wei.mat'],'modi','modf','invert','fish','metrics','iskill',...
    'fskill','Ti','Tf');

%%
oi = abs(1-iskill(1,:));
zi = abs(0-iskill(2,:));
[oi,id1] = sort(oi);
ir = sims(id1);
[zi,id2] = sort(zi);
irmse = sims(id2);
[xi,id3] = sort(iskill(6,:),'descend');
imef = sims(id3);

ir(1:3)'
irmse(1:3)'
imef(1:3)'

of = abs(1-fskill(1,:));
zf = abs(0-fskill(2,:));
[of,fd1] = sort(of);
fr = sims(fd1);
[zf,fd2] = sort(zf);
frmse = sims(fd2);
[xf,fd3] = sort(fskill(6,:),'descend');
fmef = sims(fd3);

fr(1:3)'
frmse(1:3)'
fmef(1:3)'

