%Consolodate spinup sim output into tables

clear all
close all

datap = '/Volumes/GFDL/CSV/';
d = dir(datap);
d2 = {d.name};
d2 = d2';
d3 = d2(4:161);

sname = 'Spinup_';
locs = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1'};
stgs={'SF','SP','SD','MF','MP','MD','LP','LD'};

Smean = NaN*ones(length(d3),length(locs),length(stgs));
Lmean = NaN*ones(length(d3),length(stgs),length(locs));
%%
for i=1:length(d3)
    
    dpath = [datap char(d3(i)) '/'];
    load([dpath sname 'lastyr_sum_mean_biom.mat'],'all_mean');
    
    for L = 1:size(all_mean,3)
        
        Smean(i,L,1) = log10(all_mean(1,1,L));
        Smean(i,L,2) = log10(all_mean(1,2,L));
        Smean(i,L,3) = log10(all_mean(1,3,L));
        Smean(i,L,4) = log10(all_mean(2,1,L));
        Smean(i,L,5) = log10(all_mean(2,2,L));
        Smean(i,L,6) = log10(all_mean(2,3,L));
        Smean(i,L,7) = log10(all_mean(3,2,L));
        Smean(i,L,8) = log10(all_mean(3,3,L));
        
        Lmean(i,1,L) = log10(all_mean(1,1,L));
        Lmean(i,2,L) = log10(all_mean(1,2,L));
        Lmean(i,3,L) = log10(all_mean(1,3,L));
        Lmean(i,4,L) = log10(all_mean(2,1,L));
        Lmean(i,5,L) = log10(all_mean(2,2,L));
        Lmean(i,6,L) = log10(all_mean(2,3,L));
        Lmean(i,7,L) = log10(all_mean(3,2,L));
        Lmean(i,8,L) = log10(all_mean(3,3,L));
        
    end
    
end

%%
spath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';

sf = table(d3,Smean(:,1,1),Smean(:,2,1),Smean(:,3,1),Smean(:,4,1),...
    Smean(:,5,1),Smean(:,6,1),Smean(:,7,1),Smean(:,8,1),Smean(:,9,1),...
    'VariableNames',['Sim' locs]);

sp = table(d3,Smean(:,1,2),Smean(:,2,2),Smean(:,3,2),Smean(:,4,2),...
    Smean(:,5,2),Smean(:,6,2),Smean(:,7,2),Smean(:,8,2),Smean(:,9,2),...
    'VariableNames',['Sim' locs]);

sd = table(d3,Smean(:,1,3),Smean(:,2,3),Smean(:,3,3),Smean(:,4,3),...
    Smean(:,5,3),Smean(:,6,3),Smean(:,7,3),Smean(:,8,3),Smean(:,9,3),...
    'VariableNames',['Sim' locs]);

mf = table(d3,Smean(:,1,4),Smean(:,2,4),Smean(:,3,4),Smean(:,4,4),...
    Smean(:,5,4),Smean(:,6,4),Smean(:,7,4),Smean(:,8,4),Smean(:,9,4),...
    'VariableNames',['Sim' locs]);

mp = table(d3,Smean(:,1,5),Smean(:,2,5),Smean(:,3,5),Smean(:,4,5),...
    Smean(:,5,5),Smean(:,6,5),Smean(:,7,5),Smean(:,8,5),Smean(:,9,5),...
    'VariableNames',['Sim' locs]);

md = table(d3,Smean(:,1,6),Smean(:,2,6),Smean(:,3,6),Smean(:,4,6),...
    Smean(:,5,6),Smean(:,6,6),Smean(:,7,6),Smean(:,8,6),Smean(:,9,6),...
    'VariableNames',['Sim' locs]);

lp = table(d3,Smean(:,1,7),Smean(:,2,7),Smean(:,3,7),Smean(:,4,7),...
    Smean(:,5,7),Smean(:,6,7),Smean(:,7,7),Smean(:,8,7),Smean(:,9,7),...
    'VariableNames',['Sim' locs]);



ld = table(d3,Smean(:,1,8),Smean(:,2,8),Smean(:,3,8),Smean(:,4,8),...
    Smean(:,5,8),Smean(:,6,8),Smean(:,7,8),Smean(:,8,8),Smean(:,9,8),...
    'VariableNames',['Sim' locs]);

%
gb = table(d3,Lmean(:,1,1),Lmean(:,2,1),Lmean(:,3,1),Lmean(:,4,1),...
    Lmean(:,5,1),Lmean(:,6,1),Lmean(:,7,1),Lmean(:,8,1),...
    'VariableNames',['Sim' stgs]);

ebs = table(d3,Lmean(:,1,2),Lmean(:,2,2),Lmean(:,3,2),Lmean(:,4,2),...
    Lmean(:,5,2),Lmean(:,6,2),Lmean(:,7,2),Lmean(:,8,2),...
    'VariableNames',['Sim' stgs]);

osp = table(d3,Lmean(:,1,3),Lmean(:,2,3),Lmean(:,3,3),Lmean(:,4,3),...
    Lmean(:,5,3),Lmean(:,6,3),Lmean(:,7,3),Lmean(:,8,3),...
    'VariableNames',['Sim' stgs]);

hot = table(d3,Lmean(:,1,4),Lmean(:,2,4),Lmean(:,3,4),Lmean(:,4,4),...
    Lmean(:,5,4),Lmean(:,6,4),Lmean(:,7,4),Lmean(:,8,4),...
    'VariableNames',['Sim' stgs]);

bats = table(d3,Lmean(:,1,5),Lmean(:,2,5),Lmean(:,3,5),Lmean(:,4,5),...
    Lmean(:,5,5),Lmean(:,6,5),Lmean(:,7,5),Lmean(:,8,5),...
    'VariableNames',['Sim' stgs]);

ns = table(d3,Lmean(:,1,6),Lmean(:,2,6),Lmean(:,3,6),Lmean(:,4,6),...
    Lmean(:,5,6),Lmean(:,6,6),Lmean(:,7,6),Lmean(:,8,6),...
    'VariableNames',['Sim' stgs]);

eep = table(d3,Lmean(:,1,7),Lmean(:,2,7),Lmean(:,3,7),Lmean(:,4,7),...
    Lmean(:,5,7),Lmean(:,6,7),Lmean(:,7,7),Lmean(:,8,7),...
    'VariableNames',['Sim' stgs]);

k2 = table(d3,Lmean(:,1,8),Lmean(:,2,8),Lmean(:,3,8),Lmean(:,4,8),...
    Lmean(:,5,8),Lmean(:,6,8),Lmean(:,7,8),Lmean(:,8,8),...
    'VariableNames',['Sim' stgs]);

s1 = table(d3,Lmean(:,1,9),Lmean(:,2,9),Lmean(:,3,9),Lmean(:,4,9),...
    Lmean(:,5,9),Lmean(:,6,9),Lmean(:,7,9),Lmean(:,8,9),...
    'VariableNames',['Sim' stgs]);

%%
writetable(sf,[spath 'logmean_bio_SF.csv']);
writetable(sp,[spath 'logmean_bio_SP.csv']);
writetable(sd,[spath 'logmean_bio_SD.csv']);
writetable(mf,[spath 'logmean_bio_MF.csv']);
writetable(mp,[spath 'logmean_bio_MP.csv']);
writetable(md,[spath 'logmean_bio_MD.csv']);
writetable(lp,[spath 'logmean_bio_LP.csv']);
writetable(ld,[spath 'logmean_bio_LD.csv']);
writetable(k2,[spath 'logmean_bio_K2.csv']);
writetable(ebs,[spath 'logmean_bio_EBS.csv']);
writetable(gb,[spath 'logmean_bio_GB.csv']);
writetable(ns,[spath 'logmean_bio_NS.csv']);
writetable(osp,[spath 'logmean_bio_OSP.csv']);
writetable(s1,[spath 'logmean_bio_S1.csv']);
writetable(hot,[spath 'logmean_bio_HOT.csv']);
writetable(bats,[spath 'logmean_bio_BATS.csv']);
writetable(eep,[spath 'logmean_bio_EEP.csv']);

