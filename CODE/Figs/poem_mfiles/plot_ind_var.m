figure
plot(SD(:,5),'k'); hold on;
plot(SD(:,11),'r'); hold on;
plot(SD(:,14),'--b'); hold on;
plot(C(:,2),'m'); hold on;
xlim([0 365])

%% consump
figure
%plot(SP(:,5),'k'); hold on;
plot(SP(:,11),'r'); hold on;
plot(SP(:,14),'--b'); hold on;
plot(C(:,2),'m'); hold on;
xlim([365 700])

%% nu
figure
plot(SP(:,15),'k'); hold on;
xlim([365 700])

%% zoop flux consumed
figure
plot(C(:,2),'m'); hold on;
%xlim([365 700])

%% clev
figure
plot(SP(:,21),'m'); hold on;
%xlim([0 365])
xlim([365 700])