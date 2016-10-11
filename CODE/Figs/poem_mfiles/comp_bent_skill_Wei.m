% Compare different skill metrics of each bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

dp = '/Volumes/GFDL/NC/';
sp = '/Volumes/GFDL/CSV/comps/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

cfile0 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';
cfile1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE075';
cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10';
cfile3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort';
cfile4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25';
cfile5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit35_MZ01_NOnmort_BE05';
cfile6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';
cfile7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10';
cfile8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MFMZ1_NOnmort_BE05';

simtex={'BE05_40','BE075_40','BE10_40','BE15_40','BE25_40','BE05_35',...
    'BE05_30','BE10_30','BE05_30_MFMZ1'};
simtex2={'BE05.40','BE075.40','BE10.40','BE15.40','BE25.40','BE05.35',...
    'BE05.30','BE10.30','BE05.30.MFMZ1'};

dpath0 = [dp cfile0 '/'];
dpath1 = [dp cfile1 '/'];
dpath2 = [dp cfile2 '/'];
dpath3 = [dp cfile3 '/'];
dpath4 = [dp cfile4 '/'];
dpath5 = [dp cfile5 '/'];
dpath6 = [dp cfile6 '/'];
dpath7 = [dp cfile7 '/'];
dpath8 = [dp cfile8 '/'];
%ppath = [pp cfile '/'];

load([dpath0 'skill_Wei.mat'],'skill');
skill0 = skill;
clear skill
load([dpath1 'skill_Wei.mat'],'skill','metrics');
skill1 = skill;
clear skill
load([dpath2 'skill_Wei.mat'],'skill');
skill2 = skill;
clear skill
load([dpath3 'skill_Wei.mat'],'skill');
skill3 = skill;
clear skill
load([dpath4 'skill_Wei.mat'],'skill');
skill4 = skill;
clear skill
load([dpath5 'skill_Wei.mat'],'skill');
skill5 = skill;
clear skill
load([dpath6 'skill_Wei.mat'],'skill');
skill6 = skill;
clear skill
load([dpath7 'skill_Wei.mat'],'skill');
skill7 = skill;
clear skill
load([dpath8 'skill_Wei.mat'],'skill');
skill8 = skill;
clear skill

%%
iskill(:,1) = skill0(:,1);
iskill(:,2) = skill1(:,1);
iskill(:,3) = skill2(:,1);
iskill(:,4) = skill3(:,1);
iskill(:,5) = skill4(:,1);
iskill(:,6) = skill5(:,1);
iskill(:,7) = skill6(:,1);
iskill(:,8) = skill7(:,1);
iskill(:,9) = skill8(:,1);

fskill(:,1) = skill0(:,2);
fskill(:,2) = skill1(:,2);
fskill(:,3) = skill2(:,2);
fskill(:,4) = skill3(:,2);
fskill(:,5) = skill4(:,2);
fskill(:,6) = skill5(:,2);
fskill(:,7) = skill6(:,2);
fskill(:,8) = skill7(:,2);
fskill(:,9) = skill8(:,2);

%% Best

one=[1;3;6;7;8];
zer=[2;4;5;9;10;11];
high=[12;13;14];
for z=1:14
    Z1=sum(z==one);
    Z2=sum(z==zer);
    Z3=sum(z==high);
    if (Z1>0)
        ibest{z,1}=simtex{find(min(abs(1-(iskill(z,:))))==abs(1-(iskill(z,:))))};
        fbest{z,1}=simtex{find(min(abs(1-(fskill(z,:))))==abs(1-(fskill(z,:))))};
    elseif (Z2>0)
        ibest{z,1}=simtex{find(min(abs(0-(iskill(z,:))))==abs(0-(iskill(z,:))))};
        fbest{z,1}=simtex{find(min(abs(0-(fskill(z,:))))==abs(0-(fskill(z,:))))};
    elseif (Z3>0)
        ibest{z,1}=simtex{find(max(iskill(z,:))==iskill(z,:))};
        fbest{z,1}=simtex{find(max(fskill(z,:))==fskill(z,:))};
    end
end


Ti=table(metrics,ibest);
Tf=table(metrics,fbest);

save([sp 'skill_Wei.mat'],'Ti','Tf','metrics','iskill','fskill','simtex');
%%
oi = abs(1-iskill(1,:));
zi = abs(0-iskill(2,:));
[oi,id1] = sort(oi);
ir = simtex(id1);
[zi,id2] = sort(zi);
irmse = simtex(id2);
[xi,id3] = sort(iskill(6,:),'descend');
imef = simtex(id3);

ir(1:3)'
irmse(1:3)'
imef(1:3)'
%%
of = abs(1-fskill(1,:));
zf = abs(0-fskill(2,:));
[of,fd1] = sort(of);
fr = simtex(fd1);
[zf,fd2] = sort(zf);
frmse = simtex(fd2);
[xf,fd3] = sort(fskill(6,:),'descend');
fmef = simtex(fd3);

fr(1:3)'
frmse(1:3)'
fmef(1:3)'

%% Bar graphs
ndp=length(simtex2);

figure(5)
orient portrait
subplot(3,2,1)
bar(iskill(1,:))
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end
title('Inverts')
ylabel('Correlation coefficient')

subplot(3,2,3)
bar(iskill(2,:))
ylabel('Root mean square error')
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,2,5)
bar(iskill(6,:))
ylabel('Modeling efficiency')
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,2,2)
bar(fskill(1,:))
title('Fish')
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,2,4)
bar(fskill(2,:))
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-0.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end

subplot(3,2,6)
bar(fskill(6,:))
set(gca,'XTick',1:ndp,'XTickLabel',[]);
for t=1:ndp
    text(t,-1.05,simtex2{t},'HorizontalAlignment','right','Rotation',45)
end
print('-dpng',[pp 'Spinup_global_corr_rmse_mef_Wei'])

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


% INVERTS
% Taylor diagram using Joliff corr coeff
[rmsd,it]=sort(iskill(10,:),'descend');
theta=acos(iskill(7,:));    %Joliff corr coeff
rho=iskill(8,:);            %stdev
simtex=simtex2;
simtex{ndp+1}='obs';

% % Just those <2
% id = find(rho<2);
% theta = theta(id);
% rho = rho(id);
% simtex=sims(id);
% simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(6)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:ndp
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Inverts Jolliff Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[pp 'Spinup_global_invert_Taylor_Joliff_Wei2'])

% Taylor diagram using corr coeff calculated
theta2=acos(iskill(1,:)); %corr coeff

% Just those <2
% id = find(rho<2);
% theta2 = theta2(id);
% rho = rho(id);
% simtex=sims(id);
% simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(7)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:ndp
    h=polar(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Inverts Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[pp 'Spinup_global_invert_Taylor_Wei2'])


%% FISH
% Taylor diagram using Joliff corr coeff

[rmsd,it]=sort(fskill(10,:),'descend');
theta=acos(fskill(7,:));
rho=fskill(8,:);

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
for s=1:ndp
    h=polar(theta(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Fish Jolliff Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[pp 'Spinup_global_fish_Taylor_Joliff_Wei2'])

% Taylor diagram using corr coeff calculated
theta2=acos(fskill(1,:));

% Just those <1.5
% id = find(rho<1.5);
% theta2 = theta2(id);
% rho = rho(id);
% simtex=sims(id);
% simtex{length(id)+1}='obs';

tr=0;
rr=1;
figure(9)
h0=polar(0,1.5,'.'); hold on;
set(h0,'color','w');
for s=1:ndp
    h=polar(theta2(s),rho(s),'.'); hold on;
    set(h,'color',cm{s},'MarkerSize',25);
    %     set(h,'MarkerSize',25);
end
h2=polar(tr,rr,'k*');
set(h2,'MarkerSize',10);
axis([0 1.5 0 1.5])
title('Fish Taylor diagram')
legend([' ' simtex])
legend('location','northeast')
print('-dpng',[pp 'Spinup_global_fish_Taylor_Wei2'])





