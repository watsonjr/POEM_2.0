function plotPoem(param, result)

y = result.y;
R = result.R;
B = result.B;
t = result.t;

xlimit = [min(param.wc) max(param.wu)];
[f, mortpred] = calcEncounter(y(end,:)', param);
wc = param.wc;

clf
subplot(4,1,1)
semilogy(t, R, 'r-')
hold on
semilogy(t, B(:,1),'b-','linewidth',0.5)
semilogy(t, B(:,2),'b-','linewidth',1)
semilogy(t, B(:,3),'m-','linewidth',0.5)
semilogy(t, B(:,4),'m-','linewidth',1)
semilogy(t, B(:,5),'m-','linewidth',1.5)
semilogy(t, B(:,6),'g-','linewidth',0.5)
semilogy(t, B(:,7),'g-','linewidth',1)
semilogy(t, B(:,8),'g-','linewidth',1.5)
xlabel('Time (yrs)')
ylabel('Biomass')
ylim([1e-6 1]*max(y(:)));
hold off


subplot(4,1,2)
loglog(param.wc(1:2), y(end,1:2),'b-', 'linewidth',2)
hold on
loglog(param.wc(3:5), y(end,3:5),'g-', 'linewidth',2)
loglog(param.wc(6:7), y(end,6:7),'b-', 'linewidth',2)
loglog(param.wc(8:10), y(end,8:10),'m-', 'linewidth',2)
loglog(param.wc(11:13), y(end,11:13),'g-', 'linewidth',2)
xlim(xlimit)
ylim([1e-6 1]*max(y(:)));
set(gca,'XTickLabel','')
ylabel('Biomass')


subplot(4,1,3)
semilogx(param.wc(6:7), f(6:7), 'b-', 'linewidth',2)
hold on
semilogx(param.wc(8:10), f(8:10), 'm-', 'linewidth',2)
semilogx(param.wc(11:13), f(11:13), 'g-', 'linewidth',2)
semilogx(wc([1 end]), param.fc*[1 1], 'k--')
ylim([0 1])
xlim(xlimit)
ylabel('Feeding lvl.')
set(gca,'XTickLabel','')


subplot(4,1,4)
loglog(param.wc(1:2), mortpred(1:2),'b-', 'linewidth',2)
hold on
loglog(param.wc(3:5), mortpred(3:5),'g-', 'linewidth',2)
loglog(param.wc(6:7), mortpred(6:7),'b-', 'linewidth',2)
loglog(param.wc(8:10), mortpred(8:10),'m-', 'linewidth',2)
loglog(param.wc(11:13), mortpred(11:13),'g-', 'linewidth',2)
hold on
%plot(wc(param.ixFish), mortpred(param.ixFish)'+param.mort0+param.F)
xlim(xlimit)
ylim([0.1 100])
ylabel('Mortality (yr^{-1})')
xlabel('Weight (g)')


defaultAxesVertical('double')