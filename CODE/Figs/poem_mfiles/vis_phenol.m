% Visualize phenology output of POEM

clear all
close all

dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

sname = 'Oneloc_pris_';
lname = 'NS_phenol_';
loc = 'NS';

SP = csvread([dpath sname lname 'Sml_p.csv']);
SF = csvread([dpath sname lname 'Sml_f.csv']);
SD = csvread([dpath sname lname 'Sml_d.csv']);
MP = csvread([dpath sname lname 'Med_p.csv']);
MF = csvread([dpath sname lname 'Med_f.csv']);
MD = csvread([dpath sname lname 'Med_d.csv']);
LP = csvread([dpath sname lname 'Lrg_p.csv']);
ddMF = csvread([dpath sname 'DD_' lname 'Med_f.csv']);
ddMD = csvread([dpath sname 'DD_' lname 'Med_d.csv']);
ddLP = csvread([dpath sname 'DD_' lname 'Lrg_p.csv']);
kMF = csvread([dpath sname 'K_' lname 'Med_f.csv']);
kMD = csvread([dpath sname 'K_' lname 'Med_d.csv']);
kLP = csvread([dpath sname 'K_' lname 'Lrg_p.csv']);
rMF = csvread([dpath sname 'Rep_' lname 'Med_f.csv']);
rMD = csvread([dpath sname 'Rep_' lname 'Med_d.csv']);
rLP = csvread([dpath sname 'Rep_' lname 'Lrg_p.csv']);

%% Plots over time
x=1:length(SP);
yfrac=x/365;
y=1861:(1/365):(1861+yfrac(end));
%y=1861:(1/365):2005;
y=y(1:(end-1));
lstd=length(SP);
id1 = 0:365:(lstd-1);
id2 = 365:365:(lstd);
ID  = [id1 id2];

figure(1)
plot(y((end-(30*365)):end),SP((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),MP((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),LP((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),SF((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),MF((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),SD((end-(30*365)):end,:),'Linewidth',2); hold on;
plot(y((end-(30*365)):end),MD((end-(30*365)):end,:),'Linewidth',2); hold on;
legend('SP','MP','LP','SF','MF','SD','MD')
legend('location','eastoutside')
xlim([y(end-(30*365)) y(end)])
xlabel('Time (y)')
ylabel('Biomass (g km^-^2)')
title(['Historic ' loc])

%%
figure(2)
subplot(3,1,1)
plot(y((end-(30*365)):end),rLP((end-(30*365)):end,:),'Linewidth',2); hold on;
xlim([y(end-(30*365)) y(end)])
ylabel('LP Repro (g km^-^2)')
title(['Historic ' loc])
subplot(3,1,2)
plot(y((end-(30*365)):end),rMF((end-(30*365)):end,:),'Linewidth',2); hold on;
xlim([y(end-(30*365)) y(end)])
ylabel('MF Repro (g km^-^2)')
subplot(3,1,3)
plot(y((end-(30*365)):end),rMD((end-(30*365)):end,:),'Linewidth',2); hold on;
xlim([y(end-(30*365)) y(end)])
xlabel('Time (y)')
ylabel('MD Repro (g km^-^2)')



