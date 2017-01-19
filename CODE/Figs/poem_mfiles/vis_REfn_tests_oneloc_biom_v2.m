% Visualize output of POEM
% Fishing spinup at one location
% 50 years

clear all
close all

datap = '/Volumes/GFDL/CSV/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

sims = {'poly3','piece','exp','2poly3','2piece'};
rfrac = {'poly3','piece','exp','2poly3','2piece'};
fcrit = 40;
nmort = 'M2';
kad = 100;
pref = 'D050';

sname = 'Spinup_';
sname2 = '';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','catch'};
cols=cols';

load('cmap_ppt_angles.mat')

ndp=length(rfrac);

%%
for i=1:length(spots)
    close all
    loc = spots{i};
    lname = [sname2 loc '_'];
    for s=1:length(rfrac)
        RE = rfrac{s};
%         dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_' pref '_nmort'  nmort '_BE05_RE' RE '_LD_fish000'];
        dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_MZ01_nmort'  nmort '_BE05_RE' RE '_BAassim_LD_fish000'];
        cfile = char(dp);
%         cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_' pref '_nmort'  nmort '_BE05_REfn_tests'];
        cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_MZ01_nmort'  nmort '_BE05_BAassim_REfn_tests'];
        dpath = [char(dp) '/'];
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
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
            stamp(cfile2)
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
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
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
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
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
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
        end
        title('L')
        xlabel('Sim')
        
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
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All stages'])
        
        sumspec = squeeze(nansum(nansum(all_mean)));
        
        f15=figure(15);
        plot(s,log10(sumspec(i)),'k.','MarkerSize',25); hold on;
        xlim([0 ndp+1])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',sims);
            stamp(cfile2)
        end
        ylabel('log10 Mean Biom (g m^-^2) in final year')
        title([loc ' All fishes and stages'])
        
    end
    
    print(f1,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'Logmean_biomass.png'])
    print(f21,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'Logmean_biomass_axes.png'])
    print(f15,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'tot_mean_biomass_spec.png'])
    print(f16,'-dpng',[figp loc '/' sname sname2 cfile2 '_' lname 'tot_mean_biomass_type.png'])
    
end

%%
for i=1:length(spots)
    loc = spots{i};
    lname = [sname2 loc '_'];
    
    %%
    for s=1:length(rfrac)
        RE = rfrac{s};
%         dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_' pref '_nmort'  nmort '_BE05_RE' RE '_LD_fish000'];
        dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_MZ01_nmort'  nmort '_BE05_RE' RE '_BAassim_LD_fish000'];
        cfile = char(dp);
%         cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%             '_' pref '_nmort'  nmort '_BE05_REfn_tests'];
        cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
            '_MZ01_nmort'  nmort '_BE05_BAassim_REfn_tests'];
        dpath = [char(dp) '/'];
        load([datap dpath sname sname2 'lastyr_sum_mean_biom']);
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f17=figure(17);
        subplot(4,3,i)
        plot(s-0.1,log10(fishsp(1,i)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(s,log10(fishsp(2,i)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(s+0.1,log10(fishsp(3,i)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-2 1])
        set(gca,'XTick',1:9,'XTickLabel',[])
        if (s==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
        end
        if (i==4)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title([loc ' All stages'])
        
        
    end
    
end
print(f17,'-dpng',[figp sname sname2 cfile2 '_tot_mean_biomass_type_all_locs.png'])


