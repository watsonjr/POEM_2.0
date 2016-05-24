%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together
%Plots of constant and phenology together
%ppt coloring
%Recruitment variability

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefPD/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefPD/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

sname = 'Oneloc_hist_';

load([dpath 'onelocs_hist_recruitment_PV.mat']);

load('cmap_ppt_angles.mat')

spots = {'GB','EBS','OSP','HOT','BATS','NS'};


%% Each on same axes
for s=1:length(spots)
    loc = spots{s};
    f9 = figure(9);
    subplot(2,3,s)
    plot(0.9,F(2*s -1),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1.9,P(2*s -1),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(2.9,D(2*s -1),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    plot(1.1,F(2*s),'^k',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(2.1,P(2*s),'^k',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(3.1,D(2*s),'^k',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 4])
    ylim([0 1])
    set(gca,'XTick',1:3,'XTickLabel',{'F','P','D'})
    title(loc)
    ylabel('Recruitment variability (PV)')
    
    
end

print(f9,'-dpng',[fpath sname 'All_oneloc_rec_PV.png'])
