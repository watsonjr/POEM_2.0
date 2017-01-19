% Visualize output of POEM
% Spinup at one location
% 100 years, but only last year saved

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ_TrefOrig/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ_TrefOrig/';

spots = {'GB','EBS','OSP','HOT','BATS','NS'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Oneloc_hist_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),2);
R=NaN*ones(30,3,length(spots));
for s=1:length(spots)
    %%
    close all
    loc = spots{s};
    %lname = [loc '_'];
    lname = ['phen_' loc '_'];
    
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    C = csvread([dpath sname lname 'Cobalt.csv']);
    z(:,1) = C(:,2);
    z(:,2) = C(:,3);
    
    %% Plots over time
    x=1:length(SP);
    y=x/365;
    lstd=length(SP);
    id1 = 0:365:(lstd-1);
    id2 = 365:365:(lstd);
    ID  = [id1 id2];
    
    %% Mean consumption level
    c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
    mclev(s,:) = nanmean(c);
    
    % Zoop overconsumption
    Zcon(s,:) = nansum(z)/lstd;
    
    
    %% Large Pelagic
    figure(1)
    subplot(4,1,1)
    plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    %print('-dpng',[fpath sname lname 'oneloc_P_time.png'])
    
    %% Forage Fish
    figure(2)
    subplot(3,1,1)
    plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
    plot(y,MF(:,1),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Forage Fishes ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Immature','Adults')
    
    subplot(3,1,2)
    plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Immature')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,1,3)
    plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    %print('-dpng',[fpath sname lname 'oneloc_F_time.png'])
    
    %% Demersal
    figure(3)
    subplot(4,1,1)
    plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
    plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
    plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Demersal Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    %print('-dpng',[fpath sname lname 'oneloc_D_time.png'])
    
    %% All biomass in subplots
    %SP
    figure(4)
    subplot(3,3,2)
    plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,5)
    plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,8)
    plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    %FF
    subplot(3,3,1)
    plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,4)
    plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,6)
    plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(3,3,9)
    plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    %print('-dpng',[fpath sname lname 'oneloc_all_sizes_sub.png'])
    
    %% All size classes of all
    
    figure(5)
    plot(y,log10(SP(:,1)),'Linewidth',2); hold on;
    plot(y,log10(MP(:,1)),'Linewidth',2); hold on;
    plot(y,log10(LP(:,1)),'Linewidth',2); hold on;
    plot(y,log10(SF(:,1)),'Linewidth',2); hold on;
    plot(y,log10(MF(:,1)),'Linewidth',2); hold on;
    plot(y,log10(SD(:,1)),'Linewidth',2); hold on;
    plot(y,log10(MD(:,1)),'Linewidth',2); hold on;
    plot(y,log10(LD(:,1)),'Linewidth',2); hold on;
    legend('SP','MP','LP','SF','MF','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    title(['Spinup ' loc])
    %print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
    
    %% Reproduction
    rep(:,1)=MF(:,1).*MF(:,18);
    rep(:,2)=LD(:,1).*LD(:,18);
    rep(:,3)=LP(:,1).*LP(:,18);
    
    figure(6)
    subplot(4,1,1)
    plot(y,rep(:,:),'Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Reproduction ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    legend('F','D','SP')
    
    subplot(4,1,2)
    plot(y,rep(:,1),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,3)
    plot(y,rep(:,2),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    
    subplot(4,1,4)
    plot(y,rep(:,3),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('Pelagic Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g m^-^2)')
    %print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
    
    %% Maturation
    m(:,1)=MF(:,19);
    m(:,2)=MD(:,19);
    m(:,3)=MP(:,19);
    m(:,4)=LD(:,19);
    m(:,5)=LP(:,19);
    
    figure(7)
    subplot(3,2,1)
    plot(y,m(:,1),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title([loc ' Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
    ylabel('Forage Fishes')
    
    subplot(3,2,3)
    plot(y,m(:,2),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,5)
    plot(y,m(:,3),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('M')
    ylabel('Pelagic Piscivores')
    xlabel('Time (y)')
    
    subplot(3,2,4)
    plot(y,m(:,4),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    ylabel('Demersal Piscivores')
    
    subplot(3,2,6)
    plot(y,m(:,5),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('L')
    xlabel('Time (y)')
    ylabel('Pelagic Piscivores')
    %print('-dpng',[fpath sname lname 'oneloc_matur_time.png'])
    
    %% Predation mortality
    %SP
    figure(8)
    subplot(2,3,2)
    plot(y,SP(:,17),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    
    subplot(2,3,5)
    plot(y,MP(:,17),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    xlabel('Time (y)')
    
    %FF
    subplot(2,3,1)
    plot(y,SF(:,17),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    ylabel('Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
    
    subplot(2,3,4)
    plot(y,MF(:,17),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    %ylabel('Biomass eaten by predators (g m^-^2)')
    
    %Detritivore
    subplot(2,3,3)
    plot(y,SD(:,17),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(2,3,6)
    plot(y,MD(:,17),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    xlabel('Time (y)')
    
    %print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub.png'])
    
    %% All consumption in subplots
    %SP
    figure(9)
    subplot(3,3,2)
    plot(y,SP(:,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title({loc; 'SP'})
    
    subplot(3,3,5)
    plot(y,MP(:,14),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MP')
    
    subplot(3,3,8)
    plot(y,LP(:,14),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LP')
    xlabel('Time (y)')
    
    %FF
    subplot(3,3,1)
    plot(y,SF(:,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SF')
    
    subplot(3,3,4)
    plot(y,MF(:,14),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MF')
    xlabel('Time (y)')
    ylabel('Biomass Consumed (g m^-^2)')
    
    %Detritivore
    subplot(3,3,3)
    plot(y,SD(:,14),'b','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('SD')
    
    subplot(3,3,6)
    plot(y,MD(:,14),'r','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('MD')
    
    subplot(3,3,9)
    plot(y,LD(:,14),'k','Linewidth',1); hold on;
    xlim([y(1) y(end)])
    title('LD')
    xlabel('Time (y)')
    %print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub.png'])
    
end

% if (length(lname) > 5)
%     save([dpath sname lname(1:5) 'consump.mat'],'mclev','Zcon');
%     csvwrite([dpath sname lname(1:5) 'clevel.csv'],mclev);
%     csvwrite([dpath sname lname(1:5) 'Zconsump.csv'],Zcon);
%     csvwrite([dpath sname lname(1:5) 'Recruits.csv'],R);
%     save([dpath sname lname(1:5) 'Recruits.mat'],'R');
% else
%     save([dpath sname 'consump.mat'],'mclev','Zcon');
%     csvwrite([dpath sname 'clevel.csv'],mclev);
%     csvwrite([dpath sname 'Zconsump.csv'],Zcon);
%     csvwrite([dpath sname 'Recruits.csv'],R);
%     save([dpath sname 'Recruits.mat'],'R');
% end
% 

