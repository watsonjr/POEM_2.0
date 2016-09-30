%Visualize benthic biomass output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath6 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath7 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath8 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath9 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath10 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath11 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit10_MZ01_NOnmort/';
% npath12 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit20_MZ01_NOnmort/';
% npath13 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
% npath14 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath15 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit50_MZ01_NOnmort/';
% npath16 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath17 = 'NoPDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath18 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath19 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
% npath20 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit05_MZ01_NOnmort/';
% npath21 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit075_MZ01_NOnmort/';
npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE025/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE075/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE10/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE20/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE30/';
% npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';
% npath1 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
% npath2 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE10/';
% npath3 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort/';
% npath4 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE20/';
% npath5 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE25/';
% npath6 = 'PDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE30/';

dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7};%;npath8;npath9};
% dp = {npath6;npath7;npath8;npath9;npath10};
% dp = {npath11;npath12;npath13;npath14;npath15};
% dp = {npath1;npath2;npath3;npath4;npath5;npath6;npath7;npath8;npath9;npath10;...
%     npath11;npath12;npath13;npath14;npath15;npath16;npath17;npath18;npath19;...
%     npath20;npath21};

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';
sims = {'0.025','0.05','0.075','0.10','0.15','0.20','0.30'};
cfile2 = 'Dc_Hartvig_cmax-metab_MFeqMP_MZ01_fcrit30_BentEff_comp';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/grid_360x200_id_locs_area_dep.mat','depth')
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/cobalt_zoop_1980.mat')

%% Predicted biomass
% Using linear equations based on depth from Wei et al. 2010
% log10 biomass in mg C /m2
%non-log10
meio = 10.^(2.18-(2.39e-4)*depth);
macro = 10.^(3.05-(5.15e-4)*depth);
mega = 10.^(1.81-(3.07e-4)*depth);
%convert from mg C to g wet weight
meio = meio * 9e-3;
macro = macro * 9e-3;
mega = mega * 9e-3;

allb = meio+macro+mega;

%%
for i=1:length(dp)
    %%
    dpath = [datap char(dp(i))];
    fpath = [figp char(dp(i))];
    cfile = char(dp(i));
    
    load([dpath sname sname2 'lastyr_sum_mean_biom.mat']);
    
    Bsum = NaN*ones(length(spots),1);
    Bmean = Bsum;
    zfspec = NaN*ones(4,length(spots));
    
    %%
    for s=1:length(spots)
        %%
        loc = spots{s};
        lname = [sname2 loc '_'];
        
        % Benthic biomass
        C = csvread([dpath sname lname 'Cobalt.csv']);
        
        t=1:length(C);
        lyr=t((end-365+1):end);
        
        Bsum(s,1) = nansum(C(lyr,1));
        Bmean(s,1) = nanmean(C(lyr,1));
        
        %% Size spectrum (sum stages) w/ zoop
        spec = nansum(all_mean(:,:,s),2);
        zfspec(1,s) = Zmean(1,s);
        zfspec(2:4,s) = spec;
        zfspec(2,s) = zfspec(2,s) + Zmean(2,s);
        
        %subplot all locations
        f6 = figure(6);
        subplot(3,3,s)
        plot(1:4,log10(zfspec(:,s)),'sk','MarkerFaceColor','k','MarkerSize',15); hold on;
        xlim([0 5])
        %ylim([-4 4])
        set(gca,'XTick',1:4,'XTickLabel',{'ZM','ZL+S','M','L'})
        title(loc)
        if (s==4)
            ylabel('log Mean Biom (g m^-^2) in final year')
        end
        xlabel('Size')
        if (s==3)
            stamp(cfile)
        end
        
        %all locations on one
        f7 = figure(7);
        plot(1:4,log10(zfspec(:,s)),'LineWidth',2); hold on;
        xlim([0 5])
        %ylim([-2 2])
        set(gca,'XTick',1:4,'XTickLabel',{'ZM','ZL+S','M','L'})
        if (s==9)
            legend(spots)
            legend('location','southwest')
        end
        ylabel('log Mean Biom (g m^-^2) in final year')
        xlabel('Size class')
        if (s==1)
            stamp(cfile)
        end
        
        %all locations on one from all sims
        f1 = figure(1);
        subplot(3,3,i)
        plot(1:4,log10(zfspec(:,s)),'LineWidth',2); hold on;
        xlim([0 5])
        ylim([-3 1])
        set(gca,'XTick',1:4,'XTickLabel',{'ZM','ZL+S','M','L'})
        title(sims{i})
        if (i==4)
            ylabel('log Mean Biom (g m^-^2) in final year')
        end
        if (s==9)
            stamp(cfile)
        end
        
        
    end
    print(f6,'-dpng',[fpath sname sname2 'All_oneloc_size_spec_zoop_sub.png'])
    print(f7,'-dpng',[fpath sname sname2 'All_oneloc_size_spec_zoop.png'])
    close(f6)
    close(f7)
    
    save([dpath sname sname2 'lastyr_bent_biom.mat'],'Bsum','Bmean','zfspec');
    
    %% Bent biom at all locations
    
    %Demersal fishes (M&L)
    dem = nansum(squeeze(all_mean(2:3,3,:)));
    BF = Bmean + dem';
    
    f15=figure(15);
    plot(1:9,log10(allb),'.r','MarkerSize',25); hold on;
    plot(1:9,log10(Bmean),'.k','MarkerSize',25); hold on;
    plot(1:9,log10(BF),'.b','MarkerSize',25); hold on;
    xlim([0 10])
    ylim([-2 2])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-2.2,spots{n},'HorizontalAlignment','center')
    end
    ylabel('log10 Mean Benthic Biom (g m^-^2) in final year')
    title('Benthic biomass')
    legend('predict','modelB','modelBF')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_bent_biomass.png'])
    
    
    f16=figure(16);
    plot(1:9,log10(allb),'.r','MarkerSize',25); hold on;
    plot(1:9,log10(Bsum),'.k','MarkerSize',25); hold on;
    xlim([0 10])
    ylim([-2 2])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-2.2,spots{n},'HorizontalAlignment','center')
    end
    ylabel('log10 Sum Benthic Biom (g m^-^2) in final year')
    title('Benthic biomass')
    legend('predict','model')
    stamp(cfile)
    print('-dpng',[fpath sname sname2 'All_oneloc_tot_sum_bent_biomass.png'])
    
    %all locations on one from all sims
    f2 = figure(2);
    subplot(3,3,i)
    plot(1:9,log10(allb),'.r','MarkerSize',25); hold on;
    plot(1:9,log10(Bmean),'.k','MarkerSize',25); hold on;
    plot(1:9,log10(BF),'.b','MarkerSize',25); hold on;
    xlim([0 10])
    ylim([-2 2])
    set(gca,'XTick',1:9,'XTickLabel',[])
    for n=1:9
        text(n,-2.2,spots{n},'HorizontalAlignment','center')
    end
    title(sims{i})
    if (i==4)
        ylabel('log10 Mean Benthic Biom (g m^-^2) in final year')
    end
    if (i==length(dp))
        stamp(cfile)
        legend('p','B','BF')
    end
    
    close(f15)
    close(f16)
    
end
print(f1,'-dpng',[cpath sname sname2 cfile2 '_size_spec_zoop_sims.png'])
print(f2,'-dpng',[cpath sname sname2 cfile2 '_mean_bent_biom_sims.png'])


