function plotPoem(param, result)

y = result.y;
R = result.R;
B = result.B;
t = result.t;

xlimit = [min(param.wc) max(param.wu)];
[f, mortpred] = calcEncounter(y(end,:)', param);
wc = param.wc;

%Feeding level
fl = NaN*ones(3,3);
fl(1:2,1) = f(6:7);
fl(1:3,2) = f(8:10);
fl(1:3,3) = f(11:13);

%Total mortality
morttot = mortpred(param.ixFish)'+param.mort0+param.F;

gre = [0 0.75 0.5];
purp = [0.75 0 0.75];

clf
subplot(4,1,1)
semilogy(t, R(:,1:2),'-','color',gre)
semilogy(t, R(:,3:5),'-','color',purp)
hold on
semilogy(t, B(:,1),'r-','linewidth',0.5)
semilogy(t, B(:,2),'r-','linewidth',1)
semilogy(t, B(:,3),'b-','linewidth',0.5)
semilogy(t, B(:,4),'b-','linewidth',1)
semilogy(t, B(:,5),'b-','linewidth',1.5)
semilogy(t, B(:,6),'k-','linewidth',0.5)
semilogy(t, B(:,7),'k-','linewidth',1)
semilogy(t, B(:,8),'k-','linewidth',1.5)
xlabel('Time (yrs)')
ylabel('Biomass')
ylim([1e-6 1]*max(y(:)));
hold off
title(['Repro effic = ' num2str(param.RE)])
    

subplot(4,1,2)
loglog(param.wc(1:2), y(end,1:2),'-','color',gre, 'linewidth',2)
hold on
loglog(param.wc(3:5), y(end,3:5),'-','color',purp, 'linewidth',2)
loglog(param.wc(6:7), y(end,6:7),'r-', 'linewidth',2)
loglog(param.wc(8:10), y(end,8:10),'b-', 'linewidth',2)
loglog(param.wc(11:13), y(end,11:13),'k-', 'linewidth',2)
xlim(xlimit)
%ylim([1e-6 1]*max(y(:)));
ylim([1e-3 1e3]);
%set(gca,'XTickLabel','')
ylabel('Biomass')
xlabel('Weight (g)')


subplot(4,1,3)
loglog(param.wc(1:2), mortpred(1:2),'-','color',gre, 'linewidth',2)
hold on
loglog(param.wc(3:5), mortpred(3:5),'-','color',purp, 'linewidth',2)
loglog(param.wc(6:7), mortpred(6:7),'r-', 'linewidth',2)
loglog(param.wc(8:10), mortpred(8:10),'b--', 'linewidth',2)
loglog(param.wc(11:13), mortpred(11:13),'k-', 'linewidth',2)
hold on
%plot(wc(param.ixFish), mortpred(param.ixFish)'+param.mort0+param.F)
loglog(param.wc(6:7), morttot(1:2),'r-', 'linewidth',1)
loglog(param.wc(8:10), morttot(3:5),'b--', 'linewidth',1)
loglog(param.wc(11:13), morttot(6:8),'k-', 'linewidth',1)
xlim(xlimit)
ylim([1e-2 100])
ylabel('Mortality (yr^{-1})')
xlabel('Weight (g)')

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