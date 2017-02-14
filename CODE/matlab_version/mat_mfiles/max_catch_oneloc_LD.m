% Find RE of max LD catch for each loc and mort

clear all
close all

datap = '/Volumes/GFDL/CSV/Matlab_big_size/';
figp = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_Big_sizes/Fishing/';

fsim = {'.1','.2','.3','.4','.5','.6','.7','.8','.9','1'};
frate = {'01','02','03','04','05','06','07','08','09','10'};
Frate = 1:10;
RE = {'1000','0500','0100','0050','0010'};
reff = [1.0,0.5,0.1,0.05,0.01];
sreff = {'1.0','0.5','0.1','0.05','0.01'};
nmrt = [0,2:5];
Mort = {'None','Hartvig','mizer','J&C','P&W'};
fcrit = 40;
kad = 100;
pref = 'D100';
BE = '05';
CC = '050';

spots = {'GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1','Aus','PUp'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';

sname = 'Spinup_';
mclev=NaN*ones(length(spots),8);
Zcon=NaN*ones(length(spots),3);

load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/poem_mfiles/cmap_ppt_angles.mat')
cmap3=cmap_ppt([5,1,3],:);

%%
maxC = NaN*ones(length(frate),length(RE),length(spots));
for n = 1:length(nmrt)
    nmort = num2str(nmrt(n));
    for r = 1:length(RE)
        
        rfrac = RE{r};
        
        %%
        for i=1:length(frate)
            F = frate{i};
            dp = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_RE' rfrac '_LD_fish' F];
            dpath = [char(dp) '/'];
            load([datap dpath sname 'lastyr_sum_mean_biom']);
            cfile2 = ['Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit' num2str(fcrit) ...
                '_' pref '_nmort' nmort '_BE' BE '_CC' CC '_LD_fishing_catch'];
            
            %% Fishing
            Dtot = sum(Dtotcatch);
            maxC(i,r,:) = Dtot;
            
        end %frate
    end %RE
    if (n==1)
        m0_maxC = maxC;
    elseif (n==2)
        m2_maxC = maxC;
    elseif (n==3)
        m3_maxC = maxC;
    elseif (n==4)
        m4_maxC = maxC;
    elseif (n==5)
        m5_maxC = maxC;
    end
end %nmort

%% Plot
leg=1:5;
leg=repmat(leg,10,1);

for s = 1:length(spots)
    f1=figure(1);
    subplot(4,3,s)
    plot(1:10,m0_maxC(:,:,s))
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{1}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg)
xlim([0 length(frate)+1])
legend(sreff)
stamp(cfile2)
print(f1,'-dpng',[figp cfile2 '_allRE_nmort0.png'])

for s = 1:length(spots)
    f2=figure(2);
    subplot(4,3,s)
    plot(1:10,m2_maxC(:,:,s))
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{2}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg)
xlim([0 length(frate)+1])
legend(sreff)
stamp(cfile2)
print(f2,'-dpng',[figp cfile2 '_allRE_nmort2.png'])

for s = 1:length(spots)
    f3=figure(3);
    subplot(4,3,s)
    plot(1:10,m3_maxC(:,:,s))
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{3}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg)
xlim([0 length(frate)+1])
legend(sreff)
stamp(cfile2)
print(f3,'-dpng',[figp cfile2 '_allRE_nmort3.png'])

for s = 1:length(spots)
    f4=figure(4);
    subplot(4,3,s)
    plot(1:10,m4_maxC(:,:,s))
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{4}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg)
xlim([0 length(frate)+1])
legend(sreff)
stamp(cfile2)
print(f4,'-dpng',[figp cfile2 '_allRE_nmort4.png'])

for s = 1:length(spots)
    f5=figure(5);
    subplot(4,3,s)
    plot(1:10,m5_maxC(:,:,s))
    xlim([0 length(frate)+1])
    if (s==2)
        str = {['M=' Mort{5}], spots{s}};
        title(str)
    else
        title(spots{s})
    end
    if (s==4)
        ylabel('Total D catch (g) in final year')
    end
end
subplot(4,3,12)
plot(1:10,leg)
xlim([0 length(frate)+1])
legend(sreff)
stamp(cfile2)
print(f5,'-dpng',[figp cfile2 '_allRE_nmort5.png'])


%% Find best
[Rmesh,Fmesh] = meshgrid(reff,Frate);
Cmax = nan*ones(length(spots),2,length(nmrt));
for s=1:11
    maxD0 = m0_maxC(:,:,s);
    m0 = find(maxD0==max(maxD0(:)));
    Cmax(s,1,1) = Rmesh(m0);
    Cmax(s,2,1) = Fmesh(m0);
    
    maxD2 = m2_maxC(:,:,s);
    m2 = find(maxD2==max(maxD2(:)));
    Cmax(s,1,2) = Rmesh(m2);
    Cmax(s,2,2) = Fmesh(m2);
    
    maxD3 = m3_maxC(:,:,s);
    m3 = find(maxD3==max(maxD3(:)));
    Cmax(s,1,3) = Rmesh(m3);
    Cmax(s,2,3) = Fmesh(m3);
    
    maxD4 = m4_maxC(:,:,s);
    m4 = find(maxD4==max(maxD4(:)));
    Cmax(s,1,4) = Rmesh(m4);
    Cmax(s,2,4) = Fmesh(m4);
    
    maxD5 = m5_maxC(:,:,s);
    m5 = find(maxD5==max(maxD5(:)));
    Cmax(s,1,5) = Rmesh(m5);
    Cmax(s,2,5) = Fmesh(m5);
end

%%
save([figp cfile2 '_allRE_allMort.mat'],'Cmax','m0_maxC','m2_maxC',...
    'm3_maxC','m4_maxC','m5_maxC');

