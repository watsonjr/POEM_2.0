% Visualize output of POEM
% Spinup at one location
% 100 years, but only last year saved

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_test_runs/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Comparisons/';

sims = {'.05','.1','.15','.2','.25','.3'};
RE = [1.0,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
CarCap = [0.25:0.25:5.0];
benteff = [0.05:0.05:0.3];
fcrit = 40;
nmort = '0';
kad = 100;
pref = 'D100';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')

%%
for i=1:length(benteff)
    BE = benteff{i};
%     dp = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
%         '_' pref '_nmort'  nmort '_BE' BE '_RE0010'];
    dp = 'NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit40_D100_nmort0_BE30_CC25_RE0010';
    dpath = [datap char(dp) '/'];
    fpath = [figp char(dp) '/'];
    cfile = char(dp);
    cfile2 = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
        '_' pref '_nmort'  nmort '_BE' BE '_RE0010_BentMortTests'];
    
    %%
    Bmean = NaN*ones(1,length(spots));
    Psum = NaN*ones(3,length(spots));
    Fsum = NaN*ones(2,length(spots));
    Pmean = Psum;
    Fmean = Fsum;
    Dmean = Psum;
    Pcon = Psum;
    Fcon = Fsum;
    Dcon = Psum;
    all_mean=NaN*ones(3,4,length(spots));
    Zcon = NaN*ones(length(spots),3);
    mclev = NaN*ones(length(spots),8);
    
    %%
    for s=1:length(spots)
        %%
        close all
        loc = spots{s};
        lname = [loc '_'];
        %lname = ['phen_' loc '_'];
        
        %     SF = Spinup_Sml_f;
        %     SP = Spinup_Sml_p;
        %     SD = Spinup_Sml_d;
        %     MF = Spinup_Med_f;
        %     MP = Spinup_Med_p;
        %     MD = Spinup_Med_d;
        %     LP = Spinup_Lrg_p;
        %     LD = Spinup_Lrg_d;
        %     C  = Spinup_Cobalt;
        
        SP = csvread([dpath sname lname 'Sml_p.csv']);
        SF = csvread([dpath sname lname 'Sml_f.csv']);
        SD = csvread([dpath sname lname 'Sml_d.csv']);
        MP = csvread([dpath sname lname 'Med_p.csv']);
        MF = csvread([dpath sname lname 'Med_f.csv']);
        MD = csvread([dpath sname lname 'Med_d.csv']);
        LP = csvread([dpath sname lname 'Lrg_p.csv']);
        LD = csvread([dpath sname lname 'Lrg_d.csv']);
        C = csvread([dpath sname lname 'Cobalt.csv']);
        z(:,1) = C(:,3);
        z(:,2) = C(:,4);
        z(:,3) = C(:,5);
        
        
        %% Plots over time
        x=182:365:length(SP);
        %x=1:90:length(SP);
        y=x/365;
        lstd=length(SP);
        t=1:lstd;
        lyr=t((end-365+1):end);
        
        %% Mean biomasses
        
        SP_mean=mean(SP(lyr,1));
        SF_mean=mean(SF(lyr,1));
        SD_mean=mean(SD(lyr,1));
        MP_mean=mean(MP(lyr,1));
        MF_mean=mean(MF(lyr,1));
        MD_mean=mean(MD(lyr,1));
        LP_mean=mean(LP(lyr,1));
        LD_mean=mean(LD(lyr,1));
        B_mean =mean(C(lyr,1));
        
        P_mean=[SP_mean;MP_mean;LP_mean];
        F_mean=[SF_mean;MF_mean];
        D_mean=[SD_mean;MD_mean;LD_mean];
        
        Pmean(:,s) = P_mean;
        Fmean(:,s) = F_mean;
        Dmean(:,s) = D_mean;
        Bmean(:,s) = B_mean;
        
        all_mean(1:2,1,s) = F_mean;
        all_mean(:,2,s) = P_mean;
        all_mean(:,3,s) = D_mean;
        all_mean(1,4,s) = B_mean;
        
        %% Mean consumption level
        c=[SF(:,21) SP(:,21) SD(:,21) MF(:,21) MP(:,21) MD(:,21) LP(:,21) LD(:,21)];
        mclev(s,:) = nanmean(c);
        
        % Zoop overconsumption
        Zcon(s,:) = nansum(z)/lstd;
        
        %% Benthic invert biomass
        figure(11)
        subplot(4,1,1)
        plot(y,log10(C(x,1)),'b','Linewidth',1); hold on;
        plot(y,log10(C(x,2)),'k','Linewidth',1); hold on;
        plot(y,log10(C(x,5)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Benthic inverts')
        xlabel('Time (y)')
        ylabel('log10')
        legend('biomass','eaten','frac')
        
        subplot(4,1,2)
        plot(y,log10(C(x,1)),'b','Linewidth',1); hold on;
        %plot(y,C(x,1),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Benthic inverts biomass')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,3)
        plot(y,log10(C(x,2)),'k','Linewidth',1); hold on;
        %plot(y,C(x,2),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Benthic inverts consumed')
        xlabel('Time (y)')
        ylabel('log10 Biomass (g m^-^2)')
        
        subplot(4,1,4)
        %plot(y,log10(C(x,5)),'r','Linewidth',1); hold on;
        plot(y,C(x,5),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Fraction of total benthic invertsconsumed')
        xlabel('Time (y)')
        ylabel('Fraction')
        print('-dpng',[dpath sname lname 'oneloc_B_time.png'])
        
        %% Demersal
        figure(3)
        subplot(4,1,1)
        plot(y,log10(SD(x,1)),'b','Linewidth',1); hold on;
        plot(y,log10(MD(x,1)),'r','Linewidth',1); hold on;
        plot(y,log10(LD(x,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title(['Spinup Demersal Piscivores ' loc])
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        legend('Larvae','Juveniles','Adults')
        
        subplot(4,1,2)
        plot(y,log10(SD(x,1)),'b','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Larvae')
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        
        subplot(4,1,3)
        plot(y,log10(MD(x,1)),'r','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Juveniles')
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        
        subplot(4,1,4)
        plot(y,log10(LD(x,1)),'k','Linewidth',1); hold on;
        xlim([y(1) y(end)])
        title('Adults')
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        print('-dpng',[dpath sname lname 'oneloc_D_time.png'])
        
        %% All size classes of all
        
        figure(5)
        plot(y,log10(SP(x,1)),'Linewidth',2); hold on;
        plot(y,log10(MP(x,1)),'Linewidth',2); hold on;
        plot(y,log10(LP(x,1)),'Linewidth',2); hold on;
        plot(y,log10(SF(x,1)),'Linewidth',2); hold on;
        plot(y,log10(MF(x,1)),'Linewidth',2); hold on;
        plot(y,log10(SD(x,1)),'Linewidth',2); hold on;
        plot(y,log10(MD(x,1)),'Linewidth',2); hold on;
        plot(y,log10(LD(x,1)),'Linewidth',2); hold on;
        legend('SP','MP','LP','SF','MF','SD','MD','LD')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-10 2])
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        title(['Spinup ' loc])
        print('-dpng',[dpath sname lname 'oneloc_all_sizes.png'])
        
        %% All fish summed over sizes
        figure(15)
        BE = SF(x,1)+MF(x,1);
        P = SP(x,1)+MP(x,1)+LP(x,1);
        D = SD(x,1)+MD(x,1)+LD(x,1);
        
        plot(y,log10(P),'k','Linewidth',2); hold on;
        plot(y,log10(BE),'b','Linewidth',2); hold on;
        plot(y,log10(D),'r','Linewidth',2); hold on;
        legend('P','F','D')
        legend('location','eastoutside')
        xlim([y(1) y(end)])
        ylim([-10 2])
        xlabel('Time (y)')
        ylabel('Biomass (g m^-^2)')
        title(['Spinup ' loc])
        print('-dpng',[dpath sname lname 'oneloc_all_type.png'])
        
        
    end
    save([dpath sname 'lastyr_sum_mean_biom'],'Pmean','Fmean','Dmean','Bmean',...
        'all_mean','z','mclev','Zcon');
    
end

%%
ndp = length(benteff);
for s=1:length(spots)
    
    loc = spots{s};
    lname = [loc '_'];
    
    %%
    for i=1:length(benteff)
        BE = benteff{i};
        dp = ['NoDc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
        '_' pref '_nmort'  nmort '_BE' BE '_RE0010'];
        dpath = [datap char(dp) '/'];
        load([dpath sname 'lastyr_sum_mean_biom']);
        
        %% Sum mean biom over stages
        fishsp = squeeze(nansum(all_mean));
        
        f17=figure(17);
        subplot(3,1,s)
        plot(i-0.1,log10(fishsp(1,s)),'sk','MarkerFaceColor',cmap_ppt(3,:),...
            'MarkerSize',15); hold on;
        plot(i,log10(fishsp(2,s)),'sk','MarkerFaceColor',cmap_ppt(1,:),...
            'MarkerSize',15); hold on;
        plot(i+0.1,log10(fishsp(3,s)),'sk','MarkerFaceColor',cmap_ppt(2,:),...
            'MarkerSize',15); hold on;
        xlim([0 ndp+1])
        ylim([-5 0])
        if (i==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-5.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
        end
        if (s==2)
            ylabel('log10 Mean Biom (g m^-^2) in final year')
        end
        title([loc ' All stages'])
        
        
        %% Mean bent biom
        f19=figure(19);
        subplot(3,1,s)
        plot(i,log10(all_mean(1,4,s)),'k.','MarkerSize',25); hold on;
        xlim([0 ndp+1])
        ylim([-2 0])
        set(gca,'XTick',1:ndp,'XTickLabel',[])
        if (i==ndp)
            set(gca,'XTick',1:ndp,'XTickLabel',[]);
            for t=1:ndp
                text(t,-2.1,sims{t},'Rotation',45,'HorizontalAlignment','right')
            end
            stamp(cfile2)
        end
        if (s==2)
            ylabel('log10 Mean Bent Biom (g m^-^2) in final year')
        end
        title([loc ' Bent'])
        
    end
    
end
print(f17,'-dpng',[figp 'Matlab_' sname cfile2 '_tot_mean_biomass_type_all_locs.png'])
print(f19,'-dpng',[figp 'Matlab_' sname cfile2 '_mean_bent_biomass_all_locs.png'])

