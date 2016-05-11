% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

spots = {'GB','EBS','OSP','HOT','BATS','NS'};

mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),2);
for s=1:length(spots)
    %%
    close all
    loc = spots{s};
    sname = 'Spinup_';
    lname = [loc '_'];
    %lname = ['phen_' loc '_'];
    P(:,1) = csvread([dpath sname lname 'Sml_p.csv']);
    P(:,2) = csvread([dpath sname lname 'Med_p.csv']);
    P(:,3) = csvread([dpath sname lname 'Lrg_p.csv']);
    F(:,1) = csvread([dpath sname lname 'Sml_f.csv']);
    F(:,2) = csvread([dpath sname lname 'Med_f.csv']);
    D(:,1) = csvread([dpath sname lname 'Sml_d.csv']);
    D(:,2) = csvread([dpath sname lname 'Med_d.csv']);
    D(:,3) = csvread([dpath sname lname 'Lrg_d.csv']);
    r(:,1) = csvread([dpath sname 'Rep_' lname 'Med_f.csv']);
    r(:,2) = csvread([dpath sname 'Rep_' lname 'Lrg_d.csv']);
    r(:,3) = csvread([dpath sname 'Rep_' lname 'Lrg_p.csv']);
    m(:,1) = csvread([dpath sname 'Rec_' lname 'Med_f.csv']);
    m(:,2) = csvread([dpath sname 'Rec_' lname 'Lrg_d.csv']);
    m(:,3) = csvread([dpath sname 'Rec_' lname 'Lrg_p.csv']);
    c(:,1) = csvread([dpath sname 'Clev_' lname 'Sml_p.csv']);
    c(:,2) = csvread([dpath sname 'Clev_' lname 'Med_p.csv']);
    c(:,3) = csvread([dpath sname 'Clev_' lname 'Lrg_p.csv']);
    c(:,4) = csvread([dpath sname 'Clev_' lname 'Sml_f.csv']);
    c(:,5) = csvread([dpath sname 'Clev_' lname 'Med_f.csv']);
    c(:,6) = csvread([dpath sname 'Clev_' lname 'Sml_d.csv']);
    c(:,7) = csvread([dpath sname 'Clev_' lname 'Med_d.csv']);
    c(:,8) = csvread([dpath sname 'Clev_' lname 'Lrg_d.csv']);
    z(:,1) = csvread([dpath sname lname 'ZMcon.csv']);
    z(:,2) = csvread([dpath sname lname 'ZLcon.csv']);
    
    %% Plots over time
    x=1:length(P);
    y=x/365;
    lstd=length(P);
    id1 = 0:365:(lstd-1);
    id2 = 365:365:(lstd);
    ID  = [id1 id2];
    
    %% Mean consumption level
    mclev(s,:) = nanmean(c);
    
    %% Zoop overconsumption
    Zcon(s,:) = nansum(z)/lstd;
    
    %% Piscivore
    figure(1)
    subplot(4,1,1)
    plot(y,P,'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Pelagic Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y,P(:,1),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,3)
    plot(y,P(:,2),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,4)
    plot(y,P(:,3),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_pisc_time.png'])
    
    %% Planktivore
    figure(2)
    subplot(3,1,1)
    plot(y,F,'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Forage Fishes ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    legend('Immature','Adults')
    
    subplot(3,1,2)
    plot(y,F(:,1),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Immature')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(3,1,3)
    plot(y,F(:,2),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    print('-dpng',[fpath sname lname 'oneloc_plan_time.png'])
    
    %% Detritivore
    figure(3)
    subplot(4,1,1)
    plot(y,D,'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Demersal Piscivores ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    legend('Larvae','Juveniles','Adults')
    
    subplot(4,1,2)
    plot(y,D(:,1),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Larvae')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,3)
    plot(y,D(:,2),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Juveniles')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,4)
    plot(y,D(:,3),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Adults')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_detr_time.png'])
    
    %% All size classes of all
    
    figure(4)
    plot(y,P,'Linewidth',2); hold on;
    plot(y,F,'Linewidth',2); hold on;
    plot(y,D,'Linewidth',2); hold on;
    legend('SP','MP','LP','SF','MF','SD','MD','LD')
    legend('location','eastoutside')
    xlim([y(1) y(end)])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    title(['Spinup ' loc])
    print('-dpng',[fpath sname lname 'oneloc_all_sizes.png'])
    
    %% Final mean biomass size spectrum
    t=1:length(P);
    lyr=t((end-(30*365)+1):end);
    %lyr=1:365;
    SP_sum=sum(P(lyr,1));
    SF_sum=sum(F(lyr,1));
    SD_sum=sum(D(lyr,1));
    MP_sum=sum(P(lyr,2));
    MF_sum=sum(F(lyr,2));
    MD_sum=sum(D(lyr,2));
    LP_sum=sum(P(lyr,3));
    LD_sum=sum(D(lyr,3));
    
    SP_mean=mean(P(lyr,1));
    SF_mean=mean(F(lyr,1));
    SD_mean=mean(D(lyr,1));
    MP_mean=mean(P(lyr,2));
    MF_mean=mean(F(lyr,2));
    MD_mean=mean(D(lyr,2));
    LP_mean=mean(P(lyr,3));
    LD_mean=mean(D(lyr,3));
    
    P_sum=[SP_sum;MP_sum;LP_sum];
    F_sum=[SF_sum;MF_sum];
    D_sum=[SD_sum;MD_sum;LD_sum];
    P_mean=[SP_mean;MP_mean;LP_mean];
    F_mean=[SF_mean;MF_mean];
    D_mean=[SD_mean;MD_mean;LD_mean];
    
    Pwgt = [0.0025; 2.5298; 2.5298e3];
    Fwgt = [0.0025; 2.5298];
    Dwgt = [0.0025; 2.5298; 2.5298e3];
    
    figure(5)
    subplot(2,3,1)
    bar(P_sum,'k')
    xlim([0 4])
    title('Pelagic Piscivores')
    ylabel('Total Biomass (g km^-^2)')
    subplot(2,3,4)
    bar(P_mean,'k')
    xlim([0 4])
    ylabel('Mean Biomass (g km^-^2)')
    
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
    title('Demersal Fishes')
    subplot(2,3,6)
    bar(D_mean,'r')
    xlim([0 3])
    print('-dpng',[fpath sname lname 'oneloc_all_biomass_spec.png'])
    
    %% Reproduction
    figure(10)
    subplot(4,1,1)
    plot(y,r,'Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title(['Spinup Reproduction ' loc])
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    legend('F','D','P')
    
    subplot(4,1,2)
    plot(y,r(:,1),'b','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Forage Fishes')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,3)
    plot(y,r(:,2),'r','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Demersal Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    
    subplot(4,1,4)
    plot(y,r(:,3),'k','Linewidth',2); hold on;
    xlim([y(1) y(end)])
    title('Pelagic Piscivores')
    xlabel('Time (y)')
    ylabel('Biomass (g km^-^2)')
    print('-dpng',[fpath sname lname 'oneloc_rep_time.png'])
    
    
    %% Recruitment
    yr=x;
    SFL=m(yr,1);
    SDL=m(yr,2);
    SPL=m(yr,3);
    FA=F(yr,2);
    DA=D(yr,3);
    PA=P(yr,3);
    
    st=1:365:length(yr);
    en=365:365:length(yr);
    SPy = NaN*ones(50,1);
    SFy = SPy;
    SDy = SPy;
    PAy = SPy;
    FAy = SPy;
    DAy = SPy;
    for n=1:50
        SPy(n) = nansum(SPL(st(n):en(n)));
        SFy(n) = nansum(SFL(st(n):en(n)));
        SDy(n) = nansum(SDL(st(n):en(n)));
        PAy(n) = nansum(PA(st(n):en(n)));
        FAy(n) = nansum(FA(st(n):en(n)));
        DAy(n) = nansum(DA(st(n):en(n)));
    end
    
    %
    figure(7)
    subplot(3,1,1)
    plot(1:50,SPy,'Linewidth',2); hold on;
    xlim([1 50])
    ylabel( 'Recruits (g km^-^2)')
    title({loc;'Pelagic piscivores'})
    
    subplot(3,1,2)
    plot(1:50,SFy,'Linewidth',2); hold on;
    xlim([1 50])
    ylabel('Recruits (g km^-^2)')
    title('Forage fishes')
    
    subplot(3,1,3)
    plot(1:50,SDy,'Linewidth',2); hold on;
    xlim([1 50])
    ylabel('Recruits (g km^-^2)')
    title('Demersal piscivores')
    print('-dpng',[fpath sname lname 'oneloc_recruitment.png'])
end

