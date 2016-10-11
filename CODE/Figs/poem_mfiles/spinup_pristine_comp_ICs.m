%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

%datap = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

npath1 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt140/';
npath2 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt135/';
npath3 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt130/';
npath4 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_cobalt125/';
npath5 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05/';
npath6 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-7/';
npath7 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-3/';
npath8 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e-1/';
npath9 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e1/';
npath10 = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_MZ01_NOnmort_BE05_init1e3/';

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

dp = {npath1;npath2;npath3;npath4;npath5};
sims = {'2000','1995','1990','1985','1980'};
cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_init_year_comp';

% dp = {npath6;npath5;npath7;npath8;npath9;npath10};
% sims = {'1e-7','1e-5','1e-3','1e-1','1e1','1e3'};
% cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_init_num_comp';

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

ndp = length(dp);

load('cmap_ppt_angles.mat')

%%
for i=1:length(spots)
    close all
    
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:ndp
        
        dpath = [datap char(dp(s))];
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %% Logmean biomass
        f1 = figure(1);
        subplot(3,1,1)
        plot(s-0.25,log10(squeeze(all_mean(1,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(1,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(1,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,log10(squeeze(all_mean(2,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(2,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(2,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title('M')
        
        subplot(3,1,3)
        plot(s,log10(squeeze(all_mean(3,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(3,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        title('L')
        xlabel('Sim')
        
        %% Logmean biomass
        f21 = figure(21);
        subplot(3,1,1)
        plot(s-0.25,log10(squeeze(all_mean(1,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(1,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(1,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        title([loc ' S'])
        
        subplot(3,1,2)
        plot(s-0.25,log10(squeeze(all_mean(2,1,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(3,:),'MarkerSize',15); hold on;
        plot(s,log10(squeeze(all_mean(2,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(2,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title('M')
        
        subplot(3,1,3)
        plot(s,log10(squeeze(all_mean(3,2,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(1,:),'MarkerSize',15); hold on;
        plot(s+0.25,log10(squeeze(all_mean(3,3,i))),'sk',...
            'MarkerFaceColor',cmap_ppt(2,:),'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 2])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
        end
        title('L')
        xlabel('Sim')
        
        
        %% Size spectrum (sum stages)
        spec = nansum(all_mean(:,:,i),2);
        
        f7 = figure(7);
        stamp(cfile)
        plot(1:2:6,log10(spec),'LineWidth',2); hold on;
        xlim([0 6])
        set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
        if (s==ndp)
            legend(sims)
            legend('location','northwest')
            stamp(cfile)
        end
        ylabel('log Mean Biom (g m^-^2) in final year')
        xlabel('Size class')
        title(loc)
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f16=figure(16);
        plot(s-0.1,log10(fishsp(1,i)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,log10(fishsp(2,i)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.1,log10(fishsp(3,i)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 2])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All stages'])
        
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        f15=figure(15);
        plot(s,log10(sumspec(i)),'k.','MarkerSize',25); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All fishes and stages'])
        
    end
    
    print(f1,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'Logmean_biomass.png'])
    print(f21,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'Logmean_biomass_axes.png'])
    print(f7,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'size_spec.png'])
    print(f15,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'tot_mean_biomass_spec.png'])
    print(f16,'-dpng',[fpath loc '/' sname sname2 cfile '_' lname 'tot_mean_biomass_type.png'])
    
end

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:ndp
        
        dpath = [datap char(dp(s))];
        load([dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f17=figure(17);
        subplot(3,3,i)
        plot(s-0.1,log10(fishsp(1,i)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,log10(fishsp(2,i)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.1,log10(fishsp(3,i)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 2])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile)
        end
        if (i==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title([loc ' All stages'])
        
        
    end 
    
end
print(f17,'-dpng',[fpath sname sname2 cfile '_tot_mean_biomass_type_all_locs.png'])
    

