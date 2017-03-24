function plotPoem_1D_1B(param, result)

warning('off')

R = result.R;
R(R<0) = eps;
B = result.B;
t = result.t;
y = nanmean(result.y((t(end)-364):t(end),:));
y(y<0) = eps;
f = nanmean(result.f((t(end)-364):t(end),:));
mortpred = nanmean(result.mortpred((t(end)-364):t(end),:));

xlimit = [min(param.wc) max(param.wu)];
wc = param.wc;

%Feeding level
fl = NaN*ones(3,3);
fl(1:2,1) = f(4:5);
fl(1:3,2) = f(6:8);
fl(1:3,3) = f(9:11);

%Total mortality
morttot = mortpred(param.ixFish) + param.mort0 + param.F;

gre = [0 0.75 0.5];
purp = [0.75 0 0.75];

close all
subplot(4,1,1)
semilogy(t, R(:,1:2),'-','color',gre); hold on;
semilogy(t, R(:,3),'-','color',purp)
hold on
semilogy(t, B(:,1),'r-','linewidth',0.5)
semilogy(t, B(:,2),'r-','linewidth',1)
semilogy(t, B(:,3),'b-','linewidth',0.5)
semilogy(t, B(:,4),'b-','linewidth',1)
semilogy(t, B(:,5),'b-','linewidth',1.5)
semilogy(t, B(:,6),'k-','linewidth',0.5)
semilogy(t, B(:,7),'k-','linewidth',1)
semilogy(t, B(:,8),'k-','linewidth',1.5)
xlabel('Time (days)')
ylabel('Biomass')
xlim([t(1) t(end)])
ylim([1e-6 1]*max(B(:)));
hold off
title(['Repro effic = ' num2str(param.RE)])
    

subplot(4,1,2)
loglog(param.wc(1:2), y(:,1:2),'-','color',gre, 'linewidth',2)
hold on
loglog(param.wc(3), y(:,3),'.','color',purp, 'MarkerSize',25)
loglog(param.wc(4:5), y(:,4:5),'r-', 'linewidth',2)
loglog(param.wc(6:8), y(:,6:8),'b-', 'linewidth',2)
loglog(param.wc(9:11), y(:,9:11),'k-', 'linewidth',2)
xlim(xlimit)
%ylim([1e-6 1]*max(y(:)));
ylim([1e-3 1e3]);
%set(gca,'XTickLabel','')
ylabel('Biomass')
xlabel('Weight (g)')


subplot(4,1,3)
loglog(param.wc(1:2), mortpred(1:2),'-','color',gre, 'linewidth',2)
hold on
loglog(param.wc(3), mortpred(3),'.','color',purp, 'MarkerSize',25)
loglog(param.wc(4:5), mortpred(4:5),'r-', 'linewidth',2)
loglog(param.wc(6:8), mortpred(6:8),'b--', 'linewidth',2)
loglog(param.wc(9:11), mortpred(9:11),'k-', 'linewidth',2)
hold on
%plot(wc(param.ixFish), mortpred(param.ixFish)'+param.mort0+param.F)
loglog(param.wc(4:5), morttot(1:2),'r-', 'linewidth',1)
loglog(param.wc(6:8), morttot(3:5),'b--', 'linewidth',1)
loglog(param.wc(9:11), morttot(6:8),'k-', 'linewidth',1)
xlim(xlimit)
ylim([1e-4 100])
ylabel('Mortality (d^{-1})')
xlabel('Weight (g)')

%Feeding level
subplot(4,1,4)
bar(fl); hold on;
colormap([1 0 0; 0 0 1; 0 0 0]);
plot(0:4, param.fc*ones(5,1), 'k--')
ylim([0 1])
xlim([0 4])
ylabel('Feeding level')
set(gca,'XTickLabel',{'S','M','L'})

% semilogx(param.wc(6:7), f(6:7), 'r-', 'linewidth',2)
% hold on
% semilogx(param.wc(8:10), f(8:10), 'b-', 'linewidth',2)
% semilogx(param.wc(11:13), f(11:13), 'k-', 'linewidth',2)
% semilogx(wc([1 end]), param.fc*[1 1], 'c--')
% ylim([0 1])
% xlim(xlimit)
% ylabel('Feeding lvl.')
% set(gca,'XTickLabel','')


%defaultAxesVertical('double')