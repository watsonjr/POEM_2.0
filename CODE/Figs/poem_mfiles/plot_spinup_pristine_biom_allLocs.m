%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together
%Plots of constant and phenology together
%ppt coloring

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefPD/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefPD/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

sname = 'Spinup_';
cname = '';
pname = 'phen_';

load([dpath sname cname 'lastyr_sum_mean_biom'],'all_sum');
Call_sum=all_sum;
clear all_sum

load([dpath sname pname 'lastyr_sum_mean_biom'],'all_sum');
Pall_sum=all_sum;
clear all_sum

load('cmap_ppt_angles.mat')

spots = {'GB','EBS','OSP','HOT','BATS','NS'};

%% test color order
figure
plot(1:5,ones(5,1),'color',cmap_ppt(1,:),'LineWidth',2); hold on;
plot(1:5,2*ones(5,1),'color',cmap_ppt(2,:),'LineWidth',2); hold on;
plot(1:5,3*ones(5,1),'color',cmap_ppt(3,:),'LineWidth',2); hold on;
plot(1:5,4*ones(5,1),'color',cmap_ppt(4,:),'LineWidth',2); hold on;
plot(1:5,5*ones(5,1),'color',cmap_ppt(5,:),'LineWidth',2); hold on;
plot(1:5,6*ones(5,1),'color',cmap_ppt(6,:),'LineWidth',2); hold on;
plot(1:5,7*ones(5,1),'color',cmap_ppt(7,:),'LineWidth',2); hold on;
ylim([0 8])


%% Each on same axes
for s=1:length(spots)
    loc = spots{s};
    f9 = figure(9);
    subplot(2,3,s)
    plot(0.75:2:5.75,log(squeeze(Call_sum(:,1,s))),'s','color',...
        cmap_ppt(3,:),'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log(squeeze(Call_sum(:,2,s))),'s','color',...
        cmap_ppt(1,:),'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:2:6.25,log(squeeze(Call_sum(:,3,s))),'s','color',...
        cmap_ppt(2,:),'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    plot(0.75:2:5.75,log(squeeze(Pall_sum(:,1,s))),'^','color',...
        cmap_ppt(3,:),'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log(squeeze(Pall_sum(:,2,s))),'^','color',...
        cmap_ppt(1,:),'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:2:6.25,log(squeeze(Pall_sum(:,3,s))),'^','color',...
        cmap_ppt(2,:),'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
%     if (s==4)
%         legend('F con','P con','D con','F seas','P seas','D seas')
%         legend('location','southeast')
%     end
    title(loc)
    ylabel('log Tot Biom (g m^-^2) in final year')
    xlabel('Stage')
    
    f10 = figure(10);
    subplot(2,3,s)
    plot(0.75:2:5.75,log(squeeze(Call_sum(:,1,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log(squeeze(Call_sum(:,2,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:2:6.25,log(squeeze(Call_sum(:,3,s))),'sk',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    plot(0.75:2:5.75,log(squeeze(Pall_sum(:,1,s))),'^k',...
        'MarkerFaceColor',cmap_ppt(3,:),...
        'MarkerSize',15); hold on;
    plot(1:2:6,log(squeeze(Pall_sum(:,2,s))),'^k',...
        'MarkerFaceColor',cmap_ppt(1,:),...
        'MarkerSize',15); hold on;
    plot(1.25:2:6.25,log(squeeze(Pall_sum(:,3,s))),'^k',...
        'MarkerFaceColor',cmap_ppt(2,:),...
        'MarkerSize',15); hold on;
    xlim([0 6])
    ylim([-25 15])
    set(gca,'XTick',1:2:5,'XTickLabel',{'S','M','L'})
%     if (s==4)
%         legend('F con','P con','D con','F seas','P seas','D seas')
%         legend('location','southeast')
%     end
    title(loc)
    ylabel('log Tot Biom (g m^-^2) in final year')
    xlabel('Stage')
    
%     f11 = figure(11);
%     %subplot(2,3,s)
%     plot(1:3,log(squeeze(Call_sum(:,1,s))),'s','color',...
%         cmap_ppt(1,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     plot(1:3,log(squeeze(Call_sum(:,2,s))),'s','color',...
%         cmap_ppt(2,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     plot(1:3,log(squeeze(Call_sum(:,3,s))),'s','color',...
%         cmap_ppt(3,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     plot(1:3,log(squeeze(Pall_sum(:,1,s))),'^','color',...
%         cmap_ppt(1,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     plot(1:3,log(squeeze(Pall_sum(:,2,s))),'^','color',...
%         cmap_ppt(2,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     plot(1:3,log(squeeze(Pall_sum(:,3,s))),'^','color',...
%         cmap_ppt(3,:),'LineWidth',2,...
%         'MarkerSize',15); hold on;
%     xlim([0 4])
%     ylim([-25 15])
%     if (s==4)
%         legend('F','P','D')
%         legend('location','southeast')
%     end
%     title(loc)
%     ylabel('log Tot Biom (g m^-^2) in final year')
%     xlabel('Stage')
end
%%
print(f9,'-dpng',[fpath sname 'All_oneloc_Logtot_biomass_con_phen_line1.png'])
print(f10,'-dpng',[fpath sname 'All_oneloc_Logtot_biomass_con_phen_line2.png'])
