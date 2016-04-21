% Visualize output of POEM
% Spinup at one location
% 100 years

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Spinup_BATS_';
loc = 'BATS';
SP = csvread([dpath sname 'Sml_p.csv']);
SF = csvread([dpath sname 'Sml_f.csv']);
SD = csvread([dpath sname 'Sml_d.csv']);
MP = csvread([dpath sname 'Med_p.csv']);
MF = csvread([dpath sname 'Med_f.csv']);
MD = csvread([dpath sname 'Med_d.csv']);
LP = csvread([dpath sname 'Lrg_p.csv']);

%% Plots over time
x=1:length(SP);
x=x/365;

%% Piscivore
figure(1)
subplot(1,2,1)
plot(x,SP,'b'); hold on;
plot(x,MP,'r'); hold on;
plot(x,LP,'k'); hold on;
xlim([x(1) x(end)])
title('Spinup')
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('Larvae','Juveniles','Adults')
subplot(2,2,2)
plot(x(1:730),SP(1:730),'b','Linewidth',2); hold on;
plot(x(1:730),MP(1:730),'r','Linewidth',2); hold on;
plot(x(1:730),LP(1:730),'k','Linewidth',2); hold on;
xlim([x(1) x(730)])
title(['Pelagic Piscivore ' loc])
xlabel('Time (y)')
subplot(2,2,4)
plot(x((end-731):end),SP((end-731):end),'b','Linewidth',2); hold on;
plot(x((end-731):end),MP((end-731):end),'r','Linewidth',2); hold on;
plot(x((end-731):end),LP((end-731):end),'k','Linewidth',2); hold on;
xlim([x(end-731) x(end)])
xlabel('Time (y)')
print('-dpng',[fpath sname 'oneloc_pisc_time.png'])

% Planktivore
figure(2)
subplot(1,2,1)
plot(x,SF,'b'); hold on;
plot(x,MF,'r'); hold on;
xlim([x(1) x(end)])
title('Spinup')
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('Immature','Adults')
subplot(2,2,2)
plot(x(1:730),SF(1:730),'b','Linewidth',2); hold on;
plot(x(1:730),MF(1:730),'r','Linewidth',2); hold on;
xlim([x(1) x(730)])
title(['Forage fishes ' loc])
subplot(2,2,4)
plot(x((end-731):end),SF((end-731):end),'b','Linewidth',2); hold on;
plot(x((end-731):end),MF((end-731):end),'r','Linewidth',2); hold on;
xlim([x(end-731) x(end)])
print('-dpng',[fpath sname 'oneloc_plan_time.png'])

% Detritivore
figure(3)
subplot(1,2,1)
plot(x,SD,'b'); hold on;
plot(x,MD,'r'); hold on;
xlim([x(1) x(end)])
title('Spinup')
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
legend('Immature','Adults')
subplot(2,2,2)
plot(x(1:730),SD(1:730),'b','Linewidth',2); hold on;
plot(x(1:730),MD(1:730),'r','Linewidth',2); hold on;
xlim([x(1) x(730)])
title(['Demersal ' loc])
subplot(2,2,4)
plot(x((end-731):end),SD((end-731):end),'b','Linewidth',2); hold on;
plot(x((end-731):end),MD((end-731):end),'r','Linewidth',2); hold on;
xlim([x(end-731) x(end)])
print('-dpng',[fpath sname 'oneloc_detr_time.png'])

%% All size classes of all

figure(4)
plot(x((end-731):end),SP((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),MP((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),LP((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),SF((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),MF((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),SD((end-731):end,:),'Linewidth',2); hold on;
plot(x((end-731):end),MD((end-731):end,:),'Linewidth',2); hold on;
legend('SP','MP','LP','SF','MF','SD','MD')
legend('location','eastoutside')
xlim([x(end-731) x(end)])
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
title(['Spinup ' loc])
print('-dpng',[fpath sname 'oneloc_all_sizes.png'])

%% Final mean biomass size spectrum
t=1:length(SP);
lyr=t((end-365):end);
SP_sum=sum(SP(lyr));
SF_sum=sum(SF(lyr));
SD_sum=sum(SD(lyr));
MP_sum=sum(MP(lyr));
MF_sum=sum(MF(lyr));
MD_sum=sum(MD(lyr));
LP_sum=sum(LP(lyr));

SP_mean=mean(SP(lyr));
SF_mean=mean(SF(lyr));
SD_mean=mean(SD(lyr));
MP_mean=mean(MP(lyr));
MF_mean=mean(MF(lyr));
MD_mean=mean(MD(lyr));
LP_mean=mean(LP(lyr));

P_sum=[SP_sum;MP_sum;LP_sum];
F_sum=[SF_sum;MF_sum];
D_sum=[SD_sum;MD_sum];
P_mean=[SP_mean;MP_mean;LP_mean];
F_mean=[SF_mean;MF_mean];
D_mean=[SD_mean;MD_mean];

Pwgt = [0.0025; 2.5298; 2.5298e3];
Fwgt = [0.0025; 2.5298];
Dwgt = [0.0025; 2.5298];

figure(5)
subplot(2,3,1)
bar(P_sum,'k')
xlim([0 4])
title('Pelagic Piscivores')
ylabel('Total Annual Biomass (g km^-^2)')
subplot(2,3,4)
bar(P_mean,'k')
xlim([0 4])
ylabel('Mean Annual Biomass (g km^-^2)')

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
xlim([0 3])
title('Demersal Fishes')
subplot(2,3,6)
bar(D_mean,'r')
xlim([0 3])
print('-dpng',[fpath sname 'oneloc_all_biomass_spec.png'])

%% log scale with weight

figure(6)
subplot(2,1,1)
plot(log(Pwgt),log(P_sum),'k','Linewidth',2); hold on;
plot(log(Fwgt),log(F_sum),'b','Linewidth',2); hold on;
plot(log(Dwgt),log(D_sum),'r','Linewidth',2); hold on;
xlabel('log Weight of size class (g)')
ylabel('log Total Annual Biomass (g km^-^2)')

subplot(2,1,2)
plot(log(Pwgt),log(P_mean),'k','Linewidth',2); hold on;
plot(log(Fwgt),log(F_mean),'b','Linewidth',2); hold on;
plot(log(Dwgt),log(D_mean),'r','Linewidth',2); hold on;
xlabel('log Weight of size class (g)')
ylabel('log Mean Annual Biomass (g km^-^2)')
%print('-dpng',[fpath sname 'oneloc_all_logbiomass_spec.png'])

%%
figure(7)
subplot(3,1,1)
plot(log(Pwgt),log(P_mean),'k','Linewidth',2); hold on;
xlim([-6 8])
%ylim([-40 10])
title({loc; 'Pelagic Piscivores'})

subplot(3,1,2)
plot(log(Fwgt),log(F_mean),'b','Linewidth',2); hold on;
xlim([-6 8])
title('Forage Fishes')
ylabel('log Mean Annual Biomass (g km^-^2)')

subplot(3,1,3)
plot(log(Dwgt),log(D_mean),'r','Linewidth',2); hold on;
xlim([-6 8])
title('Demersal Fishes')
xlabel('log Weight of size class (g)')
print('-dpng',[fpath sname 'oneloc_each_logbiomass_spec.png'])






