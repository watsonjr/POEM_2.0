%POEM catch vs. SAUP catch by LME

clear all
close all

spath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/SAUP/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
dp = '/Volumes/GFDL/NC/';
pp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';

cfile = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05';

dpath = [dp cfile '/'];
ppath = [pp cfile '/'];

load([dpath 'LME_hist_fished_' cfile '.mat']);
load([spath 'LME_Catch_annual.mat']);
load([cpath 'hindcast_gridspec.mat'],'dat','geolat_t','geolon_t');
load([cpath 'lme_mask_esm2m.mat']);
grid = csvread([cpath 'grid_csv.csv']);

%% 
%1950-2000 SAUP average
id = find(yr>=1950 & yr<2000);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';

%Top 10 yrs SAUP
tot_catch = sum(lme_catch);
[sort_lme_catch,id10] = sort(tot_catch,'descend');
yrs10 = yr(id10(1:10));
slme_mcatch10 = nanmean(lme_catch(id10(1:10),:));
slme_mcatch10 = slme_mcatch10';

%Top 20 yrs SAUP
yrs20 = yr(id10(1:20));
slme_mcatch20 = nanmean(lme_catch(id10(1:20),:));
slme_mcatch20 = slme_mcatch20';

%POEM LME biomass in MT
plme_mcatch = nansum(lme_mcatch00,2) * 1e-6;

%Difference
diff_catch = plme_mcatch - slme_mcatch;

code = [1:66]';
T = table(code,slme_mcatch,plme_mcatch,'VariableNames',{'lme','saup','poem'});
writetable(T,[dpath 'LME_saup_catch_hist_fished_' cfile '.csv']);

%% Fit line
fraw = fit(slme_mcatch,plme_mcatch,'poly1');
mraw = fraw.p1;
braw = fraw.p2;
yraw = mraw * slme_mcatch + braw;

l10s=log10(slme_mcatch);
l10p=log10(plme_mcatch);
ii = find(l10s~=-Inf);
l10s = l10s(ii);
l10p = l10p(ii);

flog = fit(l10s,l10p,'poly1');
mlog = flog.p1;
blog = flog.p2;
ylog = mlog * l10s + blog;

%%
x=1:0.5:8;

figure(1)
plot(10.^x,10.^x,'--k'); hold on;
plot(slme_mcatch,plme_mcatch,'.k','MarkerSize',25); hold on;
plot(slme_mcatch,yraw,'r'); hold on;
%plot(fraw,slme_mcatch,plme_mcatch); hold on;
axis([0 9e6 0 9e6])
xlabel('SAUP catch (MT)')
ylabel('POEM catch (MT)')
title('1950-2000 mean catch')
print('-dpng',[ppath 'hist_fished_SAUP_catch_comp.png'])

figure(2)
plot(x,x,'--k'); hold on;
plot(l10s,l10p,'.k','MarkerSize',25); hold on;
plot(l10s,ylog,'r'); hold on;
%plot(flog,l10s,l10p); hold on;
axis([1 7 1 7])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
title('1950-2000 mean catch')
print('-dpng',[ppath 'hist_fished_SAUP_log10catch_comp.png'])

figure(3)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10),log10(plme_mcatch),'.k','MarkerSize',25); hold on;
axis([1 7 1 7])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('1950-2000 mean catch')
print('-dpng',[ppath 'hist_fished_SAUP10_log10catch_comp.png'])

figure(4)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch20),log10(plme_mcatch),'.k','MarkerSize',25); hold on;
axis([1 7 1 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('1950-2000 mean catch')
print('-dpng',[ppath 'hist_fished_SAUP20_log10catch_comp.png'])


%% Plot by region

nwa = [6:9,18];
nea = [19:26,59:60];
nep = [1:4,10:11,54:55,64:65];
nwp = [47:53,56];
sam = 13:17;
afr = [27:31,33];
aus = 39:46;
ind = [32,34:38];
crb = [5,12];

figure(5)
subplot(2,2,1)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nwa)),log10(plme_mcatch(nwa)),'.c','MarkerSize',25); hold on;
axis([4 7 4 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Atl 1950-2000 mean catch')
for i=1:length(nwa)
    text(log10(slme_mcatch10(nwa(i))),log10(plme_mcatch(nwa(i))),num2str(nwa(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,2)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nea)),log10(plme_mcatch(nea)),'.c','MarkerSize',25); hold on;
axis([4 7 4 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Atl')
for i=1:length(nea)
    text(log10(slme_mcatch10(nea(i))),log10(plme_mcatch(nea(i))),num2str(nea(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,3)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nwp)),log10(plme_mcatch(nwp)),'.c','MarkerSize',25); hold on;
axis([3 7 3 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('NW Pac')
for i=1:length(nwp)
    text(log10(slme_mcatch10(nwp(i))),log10(plme_mcatch(nwp(i))),num2str(nwp(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,2,4)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(nep)),log10(plme_mcatch(nep)),'.c','MarkerSize',25); hold on;
axis([2 7 2 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('NE Pac')
for i=1:length(nep)
    text(log10(slme_mcatch10(nep(i))),log10(plme_mcatch(nep(i))),num2str(nep(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'hist_fished_SAUP10_log10catch_comp_AtlPac.png'])

figure(6)
subplot(2,3,1)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(afr)),log10(plme_mcatch(afr)),'.c','MarkerSize',25); hold on;
axis([4 7 4 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('Africa 1950-2000 mean catch')
for i=1:length(afr)
    text(log10(slme_mcatch10(afr(i))),log10(plme_mcatch(afr(i))),num2str(afr(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,2)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(aus)),log10(plme_mcatch(aus)),'.c','MarkerSize',25); hold on;
axis([4 7 4 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('Australia')
for i=1:length(aus)
    text(log10(slme_mcatch10(aus(i))),log10(plme_mcatch(aus(i))),num2str(aus(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,3)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(ind)),log10(plme_mcatch(ind)),'.c','MarkerSize',25); hold on;
axis([4 7 4 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('India-Indonesia')
for i=1:length(ind)
    text(log10(slme_mcatch10(ind(i))),log10(plme_mcatch(ind(i))),num2str(ind(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,4)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(crb)),log10(plme_mcatch(crb)),'.c','MarkerSize',25); hold on;
axis([5 7 5 7])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('Caribbean')
for i=1:length(crb)
    text(log10(slme_mcatch10(crb(i))),log10(plme_mcatch(crb(i))),num2str(crb(i)),...
        'Color','k','HorizontalAlignment','center')
end

subplot(2,3,5)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(sam)),log10(plme_mcatch(sam)),'.c','MarkerSize',25); hold on;
axis([5 8 5 8])
xlabel('log10 SAUP catch (MT) top 20 yrs')
ylabel('log10 POEM catch (MT)')
title('S Amer')
for i=1:length(sam)
    text(log10(slme_mcatch10(sam(i))),log10(plme_mcatch(sam(i))),num2str(sam(i)),...
        'Color','k','HorizontalAlignment','center')
end
print('-dpng',[ppath 'hist_fished_SAUP10_log10catch_comp_Other.png'])

%% Drop Arctic, Antarctic, Hawaii, Australia

keep = [1:9,11:38,46:54,59:60,62,65];

figure(7)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)),'.k','MarkerSize',25); hold on;
axis([4.5 7.5 4.5 7.5])
xlabel('log10 SAUP catch (MT) top 10 yrs')
ylabel('log10 POEM catch (MT)')
title('1950-2000 mean catch without Polar and Australia')
% for i=1:length(keep)
%     text(log10(slme_mcatch10(keep(i))),log10(plme_mcatch(keep(i))),num2str(keep(i)),...
%         'Color','k','HorizontalAlignment','center')
% end
print('-dpng',[ppath 'hist_fished_SAUP10_log10catch_comp_out.png'])

%% exclude zero catch because log=-Inf
nn=[1:63,65:66];
[rraw,praw] = corrcoef(slme_mcatch,plme_mcatch)
[rlog,plog] = corrcoef(l10s,l10p)
[rlog10,plog10] = corrcoef(log10(slme_mcatch10(nn)),log10(plme_mcatch(nn)))
[rlog20,plog20] = corrcoef(log10(slme_mcatch20(nn)),log10(plme_mcatch(nn)))
[rlog10O,plog10O] = corrcoef(log10(slme_mcatch10(keep)),log10(plme_mcatch(keep)))
[rlogO,plogO] = corrcoef(log10(slme_mcatch(keep)),log10(plme_mcatch(keep)))


figure(8)
plot(x,x,'--k'); hold on;
plot(l10s,l10p,'.k','MarkerSize',25); hold on;
axis([1 7 1 7])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
title('1950-2000 mean catch')
text(2,6,['r = ' num2str(rlog(2,1))])
print('-dpng',[ppath 'hist_fished_SAUP_log10catch_comp_rp.png'])

figure(9)
plot(x,x,'--k'); hold on;
plot(log10(slme_mcatch(keep)),log10(plme_mcatch(keep)),'.k','MarkerSize',25); hold on;
axis([1 7 1 7])
xlabel('log10 SAUP catch (MT)')
ylabel('log10 POEM catch (MT)')
text(2,6,['r = ' num2str(rlogO(2,1))])
title('1950-2000 mean catch without Polar and Australia')
print('-dpng',[ppath 'hist_fished_SAUP_log10catch_comp_out_rp.png'])

