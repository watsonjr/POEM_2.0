% Visualize output of POEM
% Spinup at 100 locations
% 30 years
% Saved as mat files

clear all
close all

dpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/NPZ/';

load([dpath 'Spinup_pristine_PISC.mat']);
load([dpath 'Spinup_pristine_PLAN.mat']);
load([dpath 'Spinup_pristine_DETR.mat']);

% CONSTANT OVER TIME
% OUTPUT PROBLEM?

%%
years=2006:2006+29;
days=1:365;
[dy,yr]=meshgrid(days,years);
%%
pi=NaN*ones(length(years)*365,100,10);
pl=pi;
de=pi;
n=0;
for y=1:length(years)
    for d=1:365
        n=n+1;
        for i=1:10
            pi(n,:,i)=A_PISC{y,d}(:,i);
            pl(n,:,i)=A_PLAN{y,d}(:,i);
            de(n,:,i)=A_DETR{y,d}(:,i);
        end
    end
end

%% Plots over time
x=1:(length(years)*365);

%Piscivore
figure(1)
%subplot(1,2,1)
plot(x,squeeze(pi(:,1,:)));
xlim([x(1) x(end)])
title('Piscivore')
xlabel('Time (d)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
%%
subplot(2,2,2)
plot(x(1:730),pi(1:730,:),'Linewidth',2)
xlim([1 730])
subplot(2,2,4)
plot(x((end-731):end),pi((end-731):end,:),'Linewidth',2)
xlim([x(end-731) x(end)])

%% Planktivore
figure(2)
subplot(1,2,1)
plot(x,squeeze(pl(:,1,:)))
xlim([x(1) x(end)])
title('Planktivore')
xlabel('Time (d)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
%%
subplot(2,2,2)
plot(x(1:730),pl(1:730,:),'Linewidth',2)
xlim([1 730])
subplot(2,2,4)
plot(x((end-366):end),pl((end-366):end,:),'Linewidth',2)
xlim([x(end-366) x(end)])

%% Detritivore
figure(3)
subplot(1,2,1)
plot(x,squeeze(de(:,1,:)))
xlim([x(1) x(end)])
title('Detritivore')
xlabel('Time (d)')
ylabel('Biomass (g km^-^2)')
legend('1','2','3','4','5','6','7','8','9','10')
%%
subplot(2,2,2)
plot(x(1:730),de(1:730,:),'Linewidth',2)
xlim([1 730])
subplot(2,2,4)
plot(x((end-366):end),de((end-366):end,:),'Linewidth',2)
xlim([x(end-366) x(end)])



%% Plots in space

lyr=x((end-365):end);
pi_sum=squeeze(sum(pi(lyr,:,:)));
pi_mean=squeeze(mean(pi(lyr,:,:)));

% pl_sum=sum(pl(lyr,:));
% pl_mean=mean(pl(lyr,:));
% de_sum=sum(de(lyr,:));
% de_mean=mean(de(lyr,:));

cpath = '/Users/Colleen/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
grid = csvread([cpath 'grid_csv.csv']);
%fix lon shift
id=find(grid(:,2)<-180);
grid(id,2)=grid(id,2)+360;

x=-180:180;
y=-90:90;
[X,Y]=meshgrid(x,y);
%%
for n=1%:10
    b=NaN*ones(length(grid),1);
    b(grid(13700:(13700+100-1),1))=pi_mean(:,n);
    Z=griddata(grid(:,2),grid(:,3),b,X,Y);
    figure
    m_proj('miller','lat',82);
    m_pcolor(X,Y,Z); hold on;
    m_coast('patch',[.5 .5 .5],'edgecolor','none');
    m_grid;
end





