%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together
%Fishing rate vs. catch

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

dpath2 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/'];
dpath3 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010/'];
dpath4 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025/'];
dpath5 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05/'];
dpath6 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10/'];
dpath7 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20/'];
dpath8 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30/'];
dpath9 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40/'];
dpath1 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_NOnmort/'];
dpath10 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish010_NOnmort/'];
dpath11 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish025_NOnmort/'];
dpath12 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish05_NOnmort/'];
dpath13 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish10_NOnmort/'];
dpath14 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish20_NOnmort/'];
dpath15 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish30_NOnmort/'];
dpath16 = [datap 'PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_fish40_NOnmort/'];

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

% dp = {dpath1;dpath10;dpath11;dpath12;dpath13;dpath14;dpath15;dpath16};
% sims = {'0.0','0.01','0.025','0.05','0.10','0.20','0.30','0.40'};
% cfile = 'MFeqMP_fishing_catch';

dp = {dpath2;dpath3;dpath4;dpath5;dpath6;dpath7;dpath8;dpath9};
sims = {'0.0+N','0.01+N','0.025+N','0.05+N','0.10+N','0.20+N','0.30+N','0.40+N'};
cfile = 'MFeqMP_fishing_nmort_catch';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

%%

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met'};
cols=cols';

load('cmap_ppt_angles.mat')

%%
for d=1:length(dp)
    
    dpath = char(dp(d));
    load([dpath sname sname2 'lastyr_sum_mean_biom']);
    
    Ftot = sum(Ftotcatch(:));
    Ptot = sum(Ptotcatch(:));
    Dtot = sum(Dtotcatch(:));
    Tot = Ftot+Ptot+Dtot;
    
    f1 = figure(1);
    plot(d,Ftot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total F catch (g) in final year')
    
    f2 = figure(2);
    plot(d,Ptot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total P catch (g) in final year')
    
    f3 = figure(3);
    plot(d,Dtot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total D catch (g) in final year')
    
    f4 = figure(4);
    plot(d,Tot,'.k','MarkerSize',15); hold on;
    xlim([0 length(dp)+1])
    if (d==length(dp))
        set(gca,'XTick',1:length(dp),'XTickLabel',sims);
        stamp(cfile)
    end
    ylabel('Total catch (g) in final year')
    
    
end
print(f1,'-dpng',[fpath sname sname2 cfile '_allF.png'])
print(f2,'-dpng',[fpath sname sname2 cfile '_allP.png'])
print(f3,'-dpng',[fpath sname sname2 cfile '_allD.png'])
print(f4,'-dpng',[fpath sname sname2 cfile '_allTotal.png'])


%%
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    close all
    %%
    for d=1:length(dp)
        
        dpath = char(dp(d));
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        
        %% Fishing
        %         MP_fish=nanmean(MP(lyr,28));
        %         MP_totcatch=nansum(MP(lyr,28));
        
        f5 = figure(5);
        plot(d,Ftotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total MF catch (g) in final year')
            stamp(cfile)
        end
        
        f6 = figure(6);
        subplot(2,1,1)
        plot(d,Ptotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total MP catch (g) in final year')
            stamp(cfile)
        end
        subplot(2,1,2)
        plot(d,Ptotcatch(3,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('Total LP catch (g) in final year')
        end
        
        f7 = figure(7);
        subplot(2,1,1)
        plot(d,Dtotcatch(2,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            title(loc)
            ylabel('Total MD catch (g) in final year')
            stamp(cfile)
        end
        subplot(2,1,2)
        plot(d,Dtotcatch(3,s),'.k','MarkerSize',15); hold on;
        xlim([0 length(dp)+1])
        if (d==length(dp))
            set(gca,'XTick',1:length(dp),'XTickLabel',sims);
            ylabel('Total LD catch (g) in final year')
        end
        
    end
    
    print(f5,'-dpng',[fpath sname sname2 cfile '_' lname 'F.png'])
    print(f6,'-dpng',[fpath sname sname2 cfile '_' lname 'P.png'])
    print(f7,'-dpng',[fpath sname sname2 cfile '_' lname 'D.png'])
    
end



