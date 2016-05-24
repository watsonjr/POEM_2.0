% Visualize output of POEM
% Pristine historical at one location
% 145 years

clear all
close all

% dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/';
% fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/';
dpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/CSV/No_PD_coupling_no_activ/';
fpath = '/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Figs/PNG/No_PD_coupling_no_activ/';

spots = {'GB','EBS','OSP','HOT','BATS','NS'};
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','egg','clev','DD','S'};
cols=cols';

sname = 'Oneloc_hist_';

%% Continuous
sp=NaN*ones(145*365,length(spots));
sf = sp;
sd = sp;
mp = sp;
mf = sp;
md = sp;
lp = sp;
ld = sp;
Rmf = sp;
Rld = sp;
Rlp = sp;
Mmf = sp;
Mld = sp;
Mlp = sp;

for s=1:length(spots)
    loc = spots{s};
    lname = [loc '_'];
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    
    sp(:,s) = SP(:,1);
    sf(:,s) = SF(:,1);
    sd(:,s) = SD(:,1);
    mp(:,s) = MP(:,1);
    mf(:,s) = MF(:,1);
    md(:,s) = MD(:,1);
    ld(:,s) = LD(:,1);
    lp(:,s) = LP(:,1);
    
    Rmf(:,s) = MF(:,1).*MF(:,18);
    Rld(:,s) = LD(:,1).*LD(:,18);
    Rlp(:,s) = LP(:,1).*LP(:,18);
    
    Mmf(:,s) = MF(:,19);
    Mld(:,s) = LD(:,19);
    Mlp(:,s) = LP(:,19);
    
end

save([dpath 'Oneloc_hist_all.mat'],'sp','sf','sd','mp',...
    'mf','md','ld','lp','Rmf','Rld','Rlp','Mmf','Mld','Mlp',...
    'spots')

%% Phenology
psp=NaN*ones(145*365,length(spots));
psf = psp;
psd = psp;
pmf = psp;
pmd = psp;
pmp = psp;
pld = psp;
plp = psp;
pDDmf = psp;
pDDld = psp;
pDDlp = psp;
pKmf = psp;
pKld = psp;
pKlp = psp;
pRmf = psp;
pRld = psp;
pRlp = psp;
pMmf = psp;
pMld = psp;
pMlp = psp;
for s=1:length(spots)
    close all
    loc = spots{s};
    lname = ['phen_' loc '_'];
    
    SP = csvread([dpath sname lname 'Sml_p.csv']);
    SF = csvread([dpath sname lname 'Sml_f.csv']);
    SD = csvread([dpath sname lname 'Sml_d.csv']);
    MP = csvread([dpath sname lname 'Med_p.csv']);
    MF = csvread([dpath sname lname 'Med_f.csv']);
    MD = csvread([dpath sname lname 'Med_d.csv']);
    LP = csvread([dpath sname lname 'Lrg_p.csv']);
    LD = csvread([dpath sname lname 'Lrg_d.csv']);
    
    psp(:,s) = SP(:,1);
    psf(:,s) = SF(:,1);
    psd(:,s) = SD(:,1);
    pmp(:,s) = MP(:,1);
    pmf(:,s) = MF(:,1);
    pmd(:,s) = MD(:,1);
    pld(:,s) = LD(:,1);
    plp(:,s) = LP(:,1);
    
    pRmf(:,s) = MF(:,1).*MF(:,18);
    pRld(:,s) = LD(:,1).*LD(:,18);
    pRlp(:,s) = LP(:,1).*LP(:,18);
    
    pMmf(:,s) = MF(:,19);
    pMld(:,s) = LD(:,19);
    pMlp(:,s) = LP(:,19);
    
    pDDmf(:,s) = MF(:,22);
    pDDld(:,s) = LD(:,22);
    pDDlp(:,s) = LP(:,22);
    
    pKmf(:,s) = MF(:,23);
    pKld(:,s) = LD(:,23);
    pKlp(:,s) = LP(:,23);
    
end

save([dpath 'Oneloc_hist_phenol_all.mat'],'psp','psf','psd',...
    'pmp','pmf','pmd','pld','plp','pDDmf','pDDld','pDDlp',...
    'pKmf','pKld','pKlp','pRmf','pRld','pRlp','pMmf','pMld',...
    'pMlp','spots')
