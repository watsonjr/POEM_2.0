%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev1e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev1e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev4e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev4e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_flev8e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_flev8e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e4_NoWgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e4_NoWgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMetab_TrefO_1e4/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMetab_TrefO_1e4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMetab_TrefO_1e4_HalfC/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMetab_TrefO_1e4_HalfC/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_4e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_4e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_8e4_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_8e4_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e5_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e5_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e6_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e6_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_Cmax_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_Cmax_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMet_TrefO_Cmax_C&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMet_TrefO_Cmax_C&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoAct_TrefO_1e6_noC&Mwgt/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoAct_TrefO_1e6_noC&Mwgt/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab_MFeatS/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab_MFeatS/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ/';

cfile = 'PDc_TrefO_KHparams_cmax-metab_MFeatS_MeatMZ';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

load([dpath sname sname2 'consump.mat'],'mclev','Zcon');

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S'};
cols=cols';

load('cmap_ppt_angles.mat')
load([dpath sname sname2 'lastyr_sum_mean_biom']);

%%
for s=1:length(spots) 
    loc = spots{s};
    %% Final mean biomass in each size 
    f1 = figure(1);
    subplot(3,3,s)
    plot(0.5:2:5.5,log10(squeeze(all_mean(:,1,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log10(squeeze(all_mean(:,2,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,log10(squeeze(all_mean(:,3,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-20 10])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        %legend('F','P','D')
        %legend('location','southeast')
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    title(loc)
    xlabel('Stage')
    if (s==7)
        stamp(cfile)
    end
    
    %% Feeding level
    f2=figure(2);
    subplot(3,3,s)
    bar(mclev(s,:),'k')
    ylim([0 1])
    set(gca,'XTickLabel',[]);
    for n=1:8
        text(n-0.5,-0.2,stage{n},'Rotation',45)
    end
    title(spots{s})
    if (s==4)
        ylabel('Feeding level')
    end
    if (s==7)
        stamp(cfile)
    end
    
    %% Growth rate (nu - energy for biomass production)
    f3 = figure(3);
    subplot(1,3,1)
    plot(s-0.25,Fmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmgr(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-0.1 0.1])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha1=gca;
        for n=1:7
            text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('S')
    
    subplot(1,3,2)
    plot(s-0.25,(Fmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmgr(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 5])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha2=gca;
        for n=1:7
            text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('M')
    
    subplot(1,3,3)
    plot(s,(Pmgr(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmgr(3,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 8])
    %ylim([-2 7])
    set(gca,'XTick',1:7,'XTickLabel',[])
    if(s==7)
        ha3=gca;
        for n=1:7
            text(n-0.5,ha3.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean biom prod rate (g g^-^1 d^-^1) in final year')
    title('L')
    if (s==7)
        stamp(cfile)
    end
    
    %% Consump per biomass (I)
    f4 = figure(4);
    subplot(3,3,s)
    plot(0.5:2:3.5,(Fcon(:,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,(Pcon(:,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.5:2:6.5,(Dcon(:,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        %legend('F','P','D')
        %legend('location','southeast')
        ylabel('Mean consumption (g g^-^1 d^-^1) in final year')
    end
    title(loc)
    xlabel('Stage')
    if (s==7)
        stamp(cfile)
    end
    
    %% Fraction zoop losses consumed
    f5 = figure(5);
    subplot(3,3,s)
    bar(z(s,:)); hold on;
    xlim([0 4])
    ylim([0 1])
    set(gca,'XTick',1:3,'XTickLabel',{'MZ','LZ','Det'})
    if (s==4)
        ylabel('Fraction flux consumed')
    end
    title(loc)
    if (s==7)
        stamp(cfile)
    end
    
    %% Size spectrum (sum stages)
    spec = nansum(all_mean(:,:,s),2);
    
    f6 = figure(6);
    subplot(3,3,s)
    plot(1:2:6,log10(spec),'sk',...
        'MarkerFaceColor','k',...
        'MarkerSize',15); hold on;
    xlim([0 6])
    %ylim([-4 4])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    title(loc)
    if (s==4)
        ylabel('log Mean Biom (g m^-^2) in final year')
    end
    xlabel('Size')
    if (s==7)
        stamp(cfile)
    end
    
    f7 = figure(7);
    %     plot(1:2:6,log10(spec),'color',cmap_ppt(s,:),...
    %         'LineWidth',2); hold on;
    plot(1:2:6,log10(spec),'LineWidth',2); hold on;
    xlim([0 6])
    %ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==7)
        legend(spots)
        legend('location','northwest')
    end
    ylabel('log Mean Biom (g m^-^2) in final year')
    xlabel('Size class')
    if (s==7)
        stamp(cfile)
    end
    
end
print(f1,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass.png'])
print(f2,'-dpng',[fpath sname sname2 'All_oneloc_con_level.png'])
print(f3,'-dpng',[fpath sname sname2 'All_oneloc_nu.png'])
print(f4,'-dpng',[fpath sname sname2 'All_oneloc_consump.png'])
print(f5,'-dpng',[fpath sname sname2 'All_oneloc_frac_zoop_loss.png'])
print(f6,'-dpng',[fpath sname sname2 'All_oneloc_size_spec_sub.png'])
print(f7,'-dpng',[fpath sname sname2 'All_oneloc_size_spec.png'])


%% Sum mean biom over stages
sumspec = squeeze(nansum(nansum(all_mean)));

figure(8);
subplot(2,1,1)
plot(1:7,log10(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 8])
%ylim([-2 1])
set(gca,'XTick',1:7,'XTickLabel',[])
% for n=1:7
%     text(n,-2.1,spots{n},'HorizontalAlignment','center')
% end
ylabel('log10 Mean Biom (g m^-^2) in final year')
title('All fishes and stages')

subplot(2,1,2)
plot(1:7,(sumspec),'k.','MarkerSize',25); hold on;
xlim([0 8])
set(gca,'XTick',1:7,'XTickLabel',[])
for n=1:7
    text(n,-0.05,spots{n},'HorizontalAlignment','center')
end
ylabel('Mean Biom (g m^-^2) in final year')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_tot_mean_biomass_spec.png'])

%% All on one
figure(9)
subplot(2,2,1)
bar(log(Fmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Forage')
ylabel('log Mean Biomass (g m^-^2) in final year',...
    'HorizontalAlignment','right')

subplot(2,2,2)
bar(log(Pmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Pel Pisc')

subplot(2,2,3)
bar(log(Dmean))
xlim([0 4])
ylim([-10 3])
xlabel('Stage')
title('Dem Pisc')

subplot(2,2,4)
bar(log(Fsum))
xlim([0 0.5])
ylim([-0.5 0.5])
legend(spots)
legend('location','west')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_biomass_spec.png'])

%% Zoop con
figure(10)
subplot(2,1,1)
bar(Zcon(:,1),'k')
ylim([0 1])
set(gca,'XTickLabel',spots);
title('Med zoo')
ylabel('Fraction of times overconsumed')

subplot(2,1,2)
bar(Zcon(:,2),'k')
ylim([0 1])
set(gca,'XTickLabel',spots);
title('Large zoo')
ylabel('Fraction of times overconsumed')
stamp(cfile)
print('-dpng',[fpath sname sname2 'All_oneloc_zoo_con.png'])




