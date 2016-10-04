% Compare different skill metrics of each bent eff simulation

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';

dp = '/Volumes/GFDL/NC/';
sp = '/Volumes/GFDL/CSV/comps/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

cfile1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05';
cfile2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10';
cfile3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort';
cfile4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25';
cfile5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10';

simtex={'BE05_40','BE10_40','BE15_40','BE25_40','BE10_30'};

dpath1 = [dp cfile1 '/'];
dpath2 = [dp cfile2 '/'];
dpath3 = [dp cfile3 '/'];
dpath4 = [dp cfile4 '/'];
dpath5 = [dp cfile5 '/'];
%ppath = [pp cfile '/'];

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

%%
iskill(:,1) = skill1(:,1);
iskill(:,2) = skill2(:,1);
iskill(:,3) = skill3(:,1);
iskill(:,4) = skill4(:,1);
iskill(:,5) = skill5(:,1);

fskill(:,1) = skill1(:,2);
fskill(:,2) = skill2(:,2);
fskill(:,3) = skill3(:,2);
fskill(:,4) = skill4(:,2);
fskill(:,5) = skill5(:,2);

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

%save([dpath 'skill_Wei.mat'],'model','obs','metrics','skill');

%% Bar graphs
simtex2={'BE05.4','BE10.4','BE15.4','BE25.4','BE10.3'};

figure(5)
orient portrait
subplot(3,2,1)
bar(iskill(1,:))
title('Inverts Correlation coefficient')
set(gca,'XTickLabel',simtex2)
%print('-depsc2',[ppath 'corr_coeff_Wei'])

subplot(3,2,3)
bar(iskill(2,:))
title('Root mean square error')
set(gca,'XTickLabel',simtex2)
%print('-depsc2',[ppath 'RMSE_Wei'])

subplot(3,2,5)
bar(iskill(6,:))
title('Modeling efficiency')
set(gca,'XTickLabel',simtex2)

subplot(3,2,2)
bar(fskill(1,:))
title('Fish Correlation coefficient')
set(gca,'XTickLabel',simtex2)
%print('-depsc2',[ppath 'corr_coeff_Wei'])

subplot(3,2,4)
bar(fskill(2,:))
title('Root mean square error')
set(gca,'XTickLabel',simtex2)
%print('-depsc2',[ppath 'RMSE_Wei'])

subplot(3,2,6)
bar(fskill(6,:))
title('Modeling efficiency')
set(gca,'XTickLabel',simtex2)
print('-dpng',[pp 'corr_rmse_mef_Wei'])