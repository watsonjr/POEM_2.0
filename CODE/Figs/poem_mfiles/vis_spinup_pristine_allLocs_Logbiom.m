%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit01/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFmet/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_NoMFpred/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFbetterMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LTmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_simpQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_compQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_bioQmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF50/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP4_fcrit10_LencF75/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_sameA/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA1_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_FdiffA2_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_Fenc2x_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA1_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFdiffA2_Lmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit20/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit30/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit40/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit05_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc15/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc20/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc25/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MFenc30/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc075/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc050/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_fcrit10_MPenc025/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_cann_Tmort/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/PDc_TrefO_KHparams_cmax-metab_MFeqMP_fcrit10_Scann_Tmort/';

fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

cfile = 'MFeqMP_fcrit10_Scann_Tmort';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

load([dpath sname sname2 'consump.mat'],'mclev','Zcon');

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod'};
cols=cols';

load('cmap_ppt_angles.mat')
load([dpath sname sname2 'lastyr_sum_mean_biom']);

%%
f21 = figure(21);
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
     
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
    if (s==1)
        ylim([-30 2])
    elseif (s==2)
        ylim([-20 2])
    elseif (s==3)
        ylim([-10 2])
    elseif (s==4)
        ylim([-100 2])
    elseif (s==5)
        ylim([-20 2])
    elseif (s==6)
        ylim([-20 2])
    elseif (s==7)
        ylim([-150 2])
    else
        ylim([-20 2])
    end
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
    if (s==4)
        ylabel('log10 Mean Biom (g m^-^2) in final year')
    end
    title(loc)
    xlabel('Stage')
end
time=clock;
chour=num2str(time(4));
cmin=num2str(time(5));
cstring=[strrep(cfile,'_','\_'),' ', date, ' ', chour,':',cmin];
subplot(3,3,6)
text(7,5,cstring,'Rotation',270)
print(f21,'-dpng',[fpath sname sname2 'All_oneloc_Logmean_biomass_axes' cfile '.png'])

%% Predation
%model lengths
L(1) = 10^((log10(2)+log10(20))/2);
L(2) = 10^((log10(20)+log10(200))/2);
L(3) = 10^((log10(200)+log10(2000))/2);
%model mass in grams
M = 0.01 .* (0.1.*L).^3;
%Andersen & Beyer mortality rate per year (natural + predation)
%physiol mort * growth constant * M^-0.25
AB = (0.35 .* 4.5 .* M.^(-0.25)) ./365;

Fmort = Fpred + Fnat;
Pmort = Ppred + Pnat;
Dmort = Dpred + Dnat;

f11 = figure(11);
for s=1:length(spots)
    
    loc = spots{s};
    lname = [sname2 loc '_'];
    
    subplot(1,2,1)
    plot(s-0.25,Fmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,Pmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,Dmort(1,s),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    set(gca,'XTick',1:9,'XTickLabel',[])
    if(s==9)
        plot(0:10,AB(1)*ones(11,1),'--k'); hold on;
        ha1=gca;
        for n=1:9
            text(n-0.5,ha1.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean mortality rate (d^-^1) in final year')
    title('S')
    
    subplot(1,2,2)
    plot(s-0.25,(Fmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(s,(Pmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(s+0.25,(Dmort(2,s)),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 10])
    set(gca,'XTick',1:9,'XTickLabel',[])
    if(s==9)
        plot(0:10,AB(2)*ones(11,1),'--k'); hold on;
        ha2=gca;
        for n=1:9
            text(n-0.5,ha2.YLim(1),spots{n},'Rotation',45)
        end
    end
    ylabel('Mean mortality rate (d^-^1) in final year')
    title('M')

end
time=clock;
chour=num2str(time(4));
cmin=num2str(time(5));
cstring=[strrep(cfile,'_','\_'),' ', date, ' ', chour,':',cmin];
subplot(1,2,2)
text(11,0.01,cstring,'Rotation',270,'HorizontalAlignment','center')
print(f11,'-dpng',[fpath sname sname2 'All_oneloc_mort_comp' cfile '.png'])



