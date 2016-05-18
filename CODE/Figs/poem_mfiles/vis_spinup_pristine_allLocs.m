%Visualize output of POEM
%Spinup at one location
%100 years
%Plots of all locations together

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/Megrey_swim_encounter_beta1/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Megrey_swim_encounter_beta1/';
% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/beta_flev5000/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/beta_flev5000/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ/';

sname = 'Spinup_';
sname2 = '';
%sname2 = 'phen_';

load([fpath sname sname2 'consump.mat'],'mclev','Zcon');

spots = {'GB','EBS','OSP','HOT','BATS','NS'};

%%
Psum = NaN*ones(3,length(spots));
Fsum = NaN*ones(2,length(spots));
Dsum = Psum;
Pmean = Psum;
Fmean = Fsum;
Dmean = Psum;
for s=1:length(spots)
    %%
    loc = spots{s};
    lname = [sname2 loc '_'];
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    
    %% Plots over time
    x=1:length(SP);
    y=x/365;
    lstd=length(SP);
    id1 = 0:365:(lstd-1);
    id2 = 365:365:(lstd);
    ID  = [id1 id2];
    
    
    %% Final mean biomass in each size
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
    
    Psum(:,s) = P_sum;
    Fsum(:,s) = F_sum;
    Dsum(:,s) = D_sum;
    Pmean(:,s) = P_mean;
    Fmean(:,s) = F_mean;
    Dmean(:,s) = D_mean;
    
    %% Total
    %Pisc
    f1=figure(1);
    subplot(2,3,s)
    bar(P_sum,'k')
    xlim([0 4])
    if (s==2)
        title({'Pel Pisc'; loc})
    else
        title(loc)
    end
    ylabel('Total Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %Forage
    f2=figure(2);
    subplot(2,3,s)
    bar(F_sum,'b')
    xlim([0 3])
    if (s==2)
        title({'Forage'; loc})
    else
        title(loc)
    end
    ylabel('Total Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %Dem
    f3=figure(3);
    subplot(2,3,s)
    bar(D_sum,'r')
    xlim([0 4])
    if (s==2)
        title({'Dem Pisc'; loc})
    else
        title(loc)
    end
    ylabel('Total Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %% Mean
    %Pisc
    f4=figure(4);
    subplot(2,3,s)
    bar(P_mean,'k')
    xlim([0 4])
    if (s==2)
        title({'Pel Pisc'; loc})
    else
        title(loc)
    end
    ylabel('Mean Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %Forage
    f5=figure(5);
    subplot(2,3,s)
    bar(F_mean,'b')
    xlim([0 3])
    if (s==2)
        title({'Forage'; loc})
    else
        title(loc)
    end
    ylabel('Mean Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %Dem
    f6=figure(6);
    subplot(2,3,s)
    bar(D_mean,'r')
    xlim([0 4])
    if (s==2)
        title({'Dem Pisc'; loc})
    else
        title(loc)
    end
    ylabel('Mean Biomass (g m^-^2) in final year')
    xlabel('Stage')
    
    %% Each on same axes
    all_sum=NaN*ones(3,3);
    all_sum(1:2,1) = F_sum;
    all_sum(:,2) = P_sum;
    all_sum(:,3) = D_sum;
    
    f10 = figure(10);
    subplot(2,3,s)
    bar(log(all_sum))
    xlim([0 4])
    ylim([-20 10])
    if (s==4)
        legend('F','P','D')
        legend('location','southeast')
    end
    title(loc)
    ylabel('log Tot Biom (g m^-^2) in final year')
    xlabel('Stage')
    
    f11 = figure(11);
    subplot(2,3,s)
    bar(all_sum)
    xlim([0 4])
    if (s==4)
        legend('F','P','D')
        legend('location','northeast')
    end
    title(loc)
    ylabel('Tot Biom (g m^-^2) in final year')
    xlabel('Stage')
    
end

print(f1,'-dpng',[fpath sname sname2 'All_oneloc_tot_biomass_spec_Pisc.png'])
print(f2,'-dpng',[fpath sname sname2 'All_oneloc_tot_biomass_spec_Forage.png'])
print(f3,'-dpng',[fpath sname sname2 'All_oneloc_tot_biomass_spec_Dem.png'])
print(f4,'-dpng',[fpath sname sname2 'All_oneloc_mean_biomass_spec_Pisc.png'])
print(f5,'-dpng',[fpath sname sname2 'All_oneloc_mean_biomass_spec_Forage.png'])
print(f6,'-dpng',[fpath sname sname2 'All_oneloc_mean_biomass_spec_Dem.png'])
print(f10,'-dpng',[fpath sname sname2 'All_oneloc_Logtot_biomass_spec_FPD.png'])
print(f11,'-dpng',[fpath sname sname2 'All_oneloc_tot_biomass_spec_FPD.png'])

%% All on one
figure(7)
subplot(2,3,1)
bar(log(Fsum))
xlim([0 4])
ylabel('log Total Biomass (g m^-^2) in final year')
legend(spots)
title('Forage')

subplot(2,3,2)
bar(log(Psum))
title('Pel Pisc')

subplot(2,3,3)
bar(log(Dsum))
title('Dem Pisc')

subplot(2,3,4)
bar(log(Fmean))
xlim([0 4])
xlabel('Stage')
ylabel('log Mean Biomass (g m^-^2) in final year')

subplot(2,3,5)
bar(log(Pmean))
xlabel('Stage')

subplot(2,3,6)
bar(log(Dmean))
xlabel('Stage')
print('-dpng',[fpath sname sname2 'All_oneloc_biomass_spec.png'])

%% Feeding level

stage={'SF','SP','SD','MF','MP','MD','LP','LD'};
for s = 1:length(spots)
    figure(8)
    subplot(2,3,s)
    bar(mclev(s,:),'k')
    set(gca,'XTickLabel',[]);
    for n=1:8
        text(n-0.5,-0.1,stage{n},'Rotation',45)
    end
    title(spots{s})
    if (s==1)
        ylabel('Feeding level')
    end
    if (s==4)
        ylabel('Feeding level')
    end
end
print('-dpng',[fpath sname sname2 'All_oneloc_con_level.png'])

%% Zoop con

figure(9)
subplot(2,1,1)
bar(Zcon(:,1),'k')
set(gca,'XTickLabel',spots);
title('Med zoo')
ylabel('Fraction of times overconsumed')

subplot(2,1,2)
bar(Zcon(:,2),'k')
set(gca,'XTickLabel',spots);
title('Large zoo')
ylabel('Fraction of times overconsumed')
print('-dpng',[fpath sname sname2 'All_oneloc_zoo_con.png'])




