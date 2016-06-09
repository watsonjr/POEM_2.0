% Visualize output of POEM
% Pristine historical at one location
% 145 years

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
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/NoPDc_NoMetab_TrefO_1e4/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/NoPDc_NoMetab_TrefO_1e4/';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP'};

cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S'};
cols=cols';

fplot=1;

mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);
for s=1:length(spots)
    %%
    close all
    loc = spots{s};
    sname = 'Spinup_';
    lname = [loc '_'];
    %lname = ['phen_' loc '_'];
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
    z(:,3) = C(:,4);
    z=floor(z);
    
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
    
    %% PLOTS
    if (fplot==1)
        %% Piscivore
        figure(1)
        subplot(4,1,1)
        plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
        plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
        plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title(['Spinup Pelagic Piscivores ' loc])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        legend('Larvae','Juveniles','Adults')
        
        subplot(4,1,2)
        plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Larvae')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,3)
        plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Juveniles')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,4)
        plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Adults')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        print('-dpng',[fpath sname lname 'oneloc_pisc_time.png'])
        
        %% Planktivore
        figure(2)
        subplot(3,1,1)
        plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
        plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title(['Spinup Forage Fishes ' loc])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        legend('Immature','Adults')
        
        subplot(3,1,2)
        plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Immature')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,1,3)
        plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Adults')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        print('-dpng',[fpath sname lname 'oneloc_plan_time.png'])
        
        %% Detritivore
        figure(3)
        subplot(4,1,1)
        plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
        plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
        plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title(['Spinup Demersal Piscivores ' loc])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        legend('Larvae','Juveniles','Adults')
        
        subplot(4,1,2)
        plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Larvae')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,3)
        plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Juveniles')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,4)
        plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Adults')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        print('-dpng',[fpath sname lname 'oneloc_detr_time.png'])
        
        %% All biomass in subplots
        %SP
        figure(4)
        subplot(3,3,2)
        plot(y,log10(SP(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title({loc; 'SP'})
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,3,5)
        plot(y,log10(MP(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MP')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,3,8)
        plot(y,log10(LP(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('LP')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        %FF
        subplot(3,3,1)
        plot(y,log10(SF(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SF')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,3,4)
        plot(y,log10(MF(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MF')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        %Detritivore
        subplot(3,3,3)
        plot(y,log10(SD(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SD')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,3,6)
        plot(y,log10(MD(:,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MD')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(3,3,9)
        plot(y,log10(LD(:,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('LD')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        print('-dpng',[fpath sname lname 'oneloc_all_sizes_sub.png'])
        
        %% All size classes of all
        
        figure(5)
        plot(y,log10(SP(:,1)),'Linewidth',1); hold on;
        plot(y,log10(MP(:,1)),'Linewidth',1); hold on;
        plot(y,log10(LP(:,1)),'Linewidth',1); hold on;
        plot(y,log10(SF(:,1)),'Linewidth',1); hold on;
        plot(y,log10(MF(:,1)),'Linewidth',1); hold on;
        plot(y,log10(SD(:,1)),'Linewidth',1); hold on;
        plot(y,log10(MD(:,1)),'Linewidth',1); hold on;
        plot(y,log10(LD(:,1)),'Linewidth',1); hold on;
        legend('SP','MP','LP','SF','MF','SD','MD','LD')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        title(['Spinup ' loc])
        print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
        
        %% Final mean biomass size spectrum
        t=1:length(SP);
        lyr=t((end-365+1):end);
        
        SP_sum=sum(SP(lyr,1));
        SF_sum=sum(SF(lyr,1));
        SD_sum=sum(SD(lyr,1));
        MP_sum=sum(MP(lyr,1));
        MF_sum=sum(MF(lyr,1));
        MD_sum=sum(MD(lyr,1));
        LP_sum=sum(LP(lyr,1));
        LD_sum=sum(LD(lyr,1));
        
        SP_mean=mean(SP(lyr,1));
        SF_mean=mean(SF(lyr,1));
        SD_mean=mean(SD(lyr,1));
        MP_mean=mean(MP(lyr,1));
        MF_mean=mean(MF(lyr,1));
        MD_mean=mean(MD(lyr,1));
        LP_mean=mean(LP(lyr,1));
        LD_mean=mean(LD(lyr,1));
        
        P_sum=[SP_sum;MP_sum;LP_sum];
        F_sum=[SF_sum;MF_sum];
        D_sum=[SD_sum;MD_sum;LD_sum];
        P_mean=[SP_mean;MP_mean;LP_mean];
        F_mean=[SF_mean;MF_mean];
        D_mean=[SD_mean;MD_mean;LD_mean];
        
        Pwgt = [0.0025; 2.5298; 2.5298e3];
        Fwgt = [0.0025; 2.5298];
        Dwgt = [0.0025; 2.5298; 2.5298e3];
        
        figure(6)
        subplot(2,3,1)
        bar(P_sum,'k')
        xlim([0 4])
        title('Pel Pisc')
        ylabel('Total Biomass (g m^-^2)')
        subplot(2,3,4)
        bar(P_mean,'k')
        xlim([0 4])
        ylabel('Mean Biomass (g m^-^2)')
        
        subplot(2,3,2)
        bar(F_sum,'b')
        xlim([0 3])
        title({loc; 'Forage Fishes'})
        xlabel('Stage')
        subplot(2,3,5)
        bar(F_mean,'b')
        xlim([0 3])
        xlabel('Stage')
        
        subplot(2,3,3)
        bar(D_sum,'r')
        xlim([0 4])
        title('Dem Pisc')
        subplot(2,3,6)
        bar(D_mean,'r')
        xlim([0 4])
        print('-dpng',[fpath sname lname 'oneloc_all_biomass_spec.png'])
        
        %% Reproduction
        rep(:,1)=MF(:,1).*MF(:,18);
        rep(:,2)=LD(:,1).*LD(:,18);
        rep(:,3)=LP(:,1).*LP(:,18);
        
        figure(7)
        subplot(4,1,1)
        plot(y,log10(rep),'Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title(['Spinup Reproduction ' loc])
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        legend('F','D','SP')
        
        subplot(4,1,2)
        plot(y,log10(rep(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Forage Fishes')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,3)
        plot(y,log10(rep(:,2)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Demersal Piscivores')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,4)
        plot(y,log10(rep(:,3)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Pelagic Piscivores')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
        
        %% Maturation
        m(:,1)=MF(:,19);
        m(:,2)=MD(:,19);
        m(:,3)=MP(:,19);
        m(:,4)=LD(:,19);
        m(:,5)=LP(:,19);
        
        figure(8)
        subplot(3,2,1)
        plot(y,log10(m(:,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title([loc ' log10 Maturation Biomass (g m^-^2)'],'HorizontalAlignment','left')
        ylabel('Forage Fishes')
        
        subplot(3,2,3)
        plot(y,log10(m(:,2)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('M')
        ylabel('Demersal Piscivores')
        
        subplot(3,2,5)
        plot(y,log10(m(:,3)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('M')
        ylabel('Pelagic Piscivores')
        xlabel('Time (y)')
        
        subplot(3,2,4)
        plot(y,log10(m(:,4)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('L')
        ylabel('Demersal Piscivores')
        
        subplot(3,2,6)
        plot(y,log10(m(:,5)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('L')
        xlabel('Time (y)')
        ylabel('Pelagic Piscivores')
        print('-dpng',[fpath sname lname 'oneloc_matur_time.png'])
        
        %% Predation mortality
        %SP
        figure(9)
        subplot(2,3,2)
        plot(y,log10(SP(:,17)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title({loc; 'SP'})
        
        subplot(2,3,5)
        plot(y,log10(MP(:,17)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MP')
        xlabel('Time (y)')
        
        %FF
        subplot(2,3,1)
        plot(y,log10(SF(:,17)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SF')
        ylabel('log10 Biomass eaten by predators (g m^-^2)','HorizontalAlignment','right')
        
        subplot(2,3,4)
        plot(y,log10(MF(:,17)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MF')
        xlabel('Time (y)')
        %ylabel('Biomass eaten by predators (g m^-^2)')
        
        %Detritivore
        subplot(2,3,3)
        plot(y,log10(SD(:,17)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SD')
        
        subplot(2,3,6)
        plot(y,log10(MD(:,17)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MD')
        xlabel('Time (y)')
        
        print('-dpng',[fpath sname lname 'oneloc_all_sizes_pred_sub.png'])
        
        %% All consumption in subplots
        %SP
        figure(10)
        subplot(3,3,2)
        plot(y,log10(SP(:,14)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title({loc; 'SP'})
        
        subplot(3,3,5)
        plot(y,log10(MP(:,14)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MP')
        
        subplot(3,3,8)
        plot(y,log10(LP(:,14)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('LP')
        xlabel('Time (y)')
        
        %FF
        subplot(3,3,1)
        plot(y,log10(SF(:,14)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SF')
        
        subplot(3,3,4)
        plot(y,log10(MF(:,14)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MF')
        xlabel('Time (y)')
        ylabel('log10 Biomass Consumed (g m^-^2)')
        
        %Detritivore
        subplot(3,3,3)
        plot(y,log10(SD(:,14)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('SD')
        
        subplot(3,3,6)
        plot(y,log10(MD(:,14)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('MD')
        
        subplot(3,3,9)
        plot(y,log10(LD(:,14)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('LD')
        xlabel('Time (y)')
        print('-dpng',[fpath sname lname 'oneloc_all_sizes_consump_sub.png'])
        
        
        %% Recruitment
        yr=x;
        SFL=m(yr,3);
        SDL=m(yr,4);
        SPL=m(yr,5);
        FA=MF(yr,1);
        DA=LD(yr,1);
        PA=LP(yr,1);
        
        st=1:365:length(yr);
        en=365:365:length(yr);
        SPy = NaN*ones(100,1);
        SFy = SPy;
        SDy = SPy;
        PAy = SPy;
        FAy = SPy;
        DAy = SPy;
        for n=1:100
            SPy(n) = nansum(SPL(st(n):en(n)));
            SFy(n) = nansum(SFL(st(n):en(n)));
            SDy(n) = nansum(SDL(st(n):en(n)));
            PAy(n) = nansum(PA(st(n):en(n)));
            FAy(n) = nansum(FA(st(n):en(n)));
            DAy(n) = nansum(DA(st(n):en(n)));
        end
        
        %
        figure(11)
        subplot(3,1,3)
        plot(1:100,log10(SPy),'k','Linewidth',2); hold on;
        xlim([1 100])
        ylabel('log10 Recruits (g m^-^2)')
        title('Pelagic piscivores')
        
        subplot(3,1,1)
        plot(1:100,log10(SFy),'b','Linewidth',2); hold on;
        xlim([1 100])
        ylabel('log10 Recruits (g m^-^2)')
        title({loc; 'Forage fishes'})
        
        subplot(3,1,2)
        plot(1:100,log10(SDy),'r','Linewidth',2); hold on;
        xlim([1 100])
        ylabel('log10 Recruits (g m^-^2)')
        title('Demersal piscivores')
        print('-dpng',[fpath sname lname 'oneloc_recruitment.png'])
    end
end

if (length(lname) > 5)
    save([dpath sname lname(1:5) 'consump.mat'],'mclev','Zcon');
    csvwrite([dpath sname lname(1:5) 'clevel.csv'],mclev);
    csvwrite([dpath sname lname(1:5) 'Zconsump.csv'],Zcon);
else
    save([dpath sname 'consump.mat'],'mclev','Zcon');
    csvwrite([dpath sname 'clevel.csv'],mclev);
    csvwrite([dpath sname 'Zconsump.csv'],Zcon);
end



