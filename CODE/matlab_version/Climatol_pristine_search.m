%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_pristine_search()

global DAYS GRD NX ID
global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl
global Z_s Z_m Z_l Lambda K_l K_j K_a fcrit h gam kt bpow
global bent_eff rfrac CC D J Sm A benc bcmx
global Tu_s Tu_m Tu_l Nat_mrt MORT
global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE
global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE
global MFsel MPsel MDsel LPsel LDsel efn cfn
global tstep K CGRD ni nj

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
Pdrpbx = '/Users/cpetrik/Dropbox/';
Fdrpbx = '/Users/Colleen/Dropbox/';
load('/Volumes/GFDL/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_grid.mat');
NX = length(GRD.Z);
ID = 1:NX;
%2D Grid for advect-diff
load([Pdrpbx 'Princeton/POEM_2.0/CODE/Data/Hindcast_cgrid_cp2D.mat']);
[ni,nj] = size(CGRD.mask);


%%%%%%%%%%%%%%% Initialize Model Variables
%! Feeding preferences
Sm = 0.25;  %Feeding 2 sizes down
J = 1.0;    %Juvenile feeding reduction
D = 0.75;   %Demersal feeding in pelagic reduction
A = 0.50;    %Adult predation reduction
% kays = 0.0405:0.01:0.125;
% bees = 0.125:0.025:0.25; %bees = 0.1:0.05:0.35;
% bees = 0.175:0.005:0.195;
% Fish = 0.1:0.1:0.6;
fqs = 0.5:0.5:3.5;
pqs = 0.25:0.25:1.5;
dqs = 1:7;

for fq=1:length(fqs)
    for pq=4;%1:length(pqs)
        for dq=7;%1:length(dqs)
            MFsel = fqs(fq);
            LPsel = pqs(pq);
            LDsel = dqs(dq);
            
            %for F=1%:length(Fish)
            %! Set fishing rate
            frate = 0.1; %Fish(F);
            dfrate = frate/365.0;
            %0=no fishing; 1=fishing
            if (frate>0)
                harv = 1;
            else
                harv = 0;
            end
            
            
            %! Make core parameters/constants (global)
            make_parameters()
            
            
            %! Create a directory for output
            tfcrit = num2str(int64(100*fcrit));
            td = num2str(1000+int64(100*LD_phi_MP));
            tj = num2str(1000+int64(100*MP_phi_S));
            tsm = num2str(1000+int64(100*MF_phi_MZ));
            ta = num2str(1000+int64(100*LP_phi_MF));
            tbe = num2str(100+int64(100*bent_eff));
            tmort = num2str(MORT);
            tcc = num2str(1000+int64(100*CC));
            tre = num2str(100000+int64(round(10000*rfrac)));
            tre2 = num2str(100000+int64(round(10000*rfrac*1)));
            if (frate >= 0.1)
                tfish = num2str(100+int64(10*frate));
                tF = num2str(1000+int64(100*frate*MFsel));
                tP = num2str(1000+int64(100*frate*LPsel));
                tD = num2str(1000+int64(100*frate*LDsel));
            else
                tfish = num2str(1000+int64(100*frate));
                tF = num2str(1000+int64(100*frate*MFsel));
                tP = num2str(1000+int64(100*frate*LPsel));
                tD = num2str(1000+int64(100*frate*LDsel));
            end
            if (MFsel > 0)
                if (LPsel > 0 && LDsel > 0)
                    sel='All';
                else
                    sel='F';
                end
            else
                if (LPsel > 0 && LDsel > 0)
                    sel = 'L';
                elseif (LPsel > 0)
                    sel = 'P';
                elseif (LDsel > 0)
                    sel = 'D';
                end
            end
            if (pdc == 0)
                coup = 'NoDc';
            elseif (pdc == 1)
                coup = 'Dc';
            elseif (pdc == 2)
                coup = 'PDc';
            end
            tcfn = num2str(h);
            tefn = num2str(round(gam));
            tkfn = num2str(100+int64(100*kt));
            tbfn = num2str(1000+int64(1000*bpow));
            tbenc = num2str(1000+int64(1000*benc));
            tbcmx = num2str(1000+int64(1000*bcmx));
            %simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
            %simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
            %simname = ['Diff_',coup,'_enc',tefn,'_cmax-metab',tcfn,'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_RE',tre(2:end)];
            %simname = [coup,'_enc',tefn,'_cmax-metab',tcfn,'_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
            %simname = [coup,'_enc',tefn,'_cm',tcfn,'_m-b200_c-b',tbfn(2:end),'_m-k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
            %simname = [coup,'_enc',tefn,'-b',tbfn(2:end),'_cm',tcfn,'_m-b175-k',tkfn(2:end),'_fcrit',tfcrit,'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
            simname = [coup,'_enc',tefn,'-b',tbenc(2:end),'_cm',tcfn,'_m-b',tbfn(2:end),'-k',tkfn(2:end),'_fcrit',tfcrit,'_c-b',tbcmx(2:end),'_D',td(2:end),'_J',tj(2:end),'_A',ta(2:end),'_Sm',tsm(2:end),'_nmort',tmort,'_BE',tbe(2:end),'_CC',tcc(2:end),'_lgRE',tre(2:end),'_mdRE',tre2(2:end)];
            if (~isdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname]))
                mkdir(['/Volumes/GFDL/NC/Matlab_new_size/',simname])
            end
            % if (~isdir([Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname]))
            %     mkdir([Pdrpbx 'Princeton/POEM_2.0/CODE/Figs/PNG/Matlab_New_sizes/',simname])
            % end
            
            %! Storage variables
            % Dims
            nt = 12;
            
            S_Bent_bio = zeros(NX,DAYS);
            
            S_Sml_f = zeros(NX,DAYS);
            S_Sml_p = zeros(NX,DAYS);
            S_Sml_d = zeros(NX,DAYS);
            S_Med_f = zeros(NX,DAYS);
            S_Med_p = zeros(NX,DAYS);
            S_Med_d = zeros(NX,DAYS);
            S_Lrg_p = zeros(NX,DAYS);
            S_Lrg_d = zeros(NX,DAYS);
            
            S_Med_f_fish = zeros(NX,DAYS);
            S_Med_p_fish = zeros(NX,DAYS);
            S_Med_d_fish = zeros(NX,DAYS);
            S_Lrg_p_fish = zeros(NX,DAYS);
            S_Lrg_d_fish = zeros(NX,DAYS);
            
            Spinup_Sml_f.bio = NaN*ones(NX,nt);
            Spinup_Sml_p.bio = NaN*ones(NX,nt);
            Spinup_Sml_d.bio = NaN*ones(NX,nt);
            Spinup_Med_f.bio = NaN*ones(NX,nt);
            Spinup_Med_p.bio = NaN*ones(NX,nt);
            Spinup_Med_d.bio = NaN*ones(NX,nt);
            Spinup_Lrg_p.bio = NaN*ones(NX,nt);
            Spinup_Lrg_d.bio = NaN*ones(NX,nt);
            Spinup_Bent.bio = NaN*ones(NX,nt);
            
            Spinup_Med_f.yield = NaN*ones(NX,nt);
            Spinup_Med_p.yield = NaN*ones(NX,nt);
            Spinup_Med_d.yield = NaN*ones(NX,nt);
            Spinup_Lrg_p.yield = NaN*ones(NX,nt);
            Spinup_Lrg_d.yield = NaN*ones(NX,nt);
            
            
            %! Initialize
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
            Med_d.td(1:NX) = 0.0;
            Lrg_d.td(1:NX) = 0.0;
            ENVR = sub_init_env(ID);
            
            
            %% %%%%%%%%%%%%%%%%%%%% Run the Model
            %! Run model with no fishing
            MNT=0;
            for YR = 1:YEARS % years
                
                for DAY = 1:DT:DAYS % days
                    
                    %%%! Future time step
                    DY = int64(ceil(DAY));
                    [num2str(YR),' , ', num2str(mod(DY,365))]
                    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                        sub_futbio(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,...
                        Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,dfrate,CC);
                    
                    if (YR==YEARS)
                        S_Bent_bio(:,DY) = BENT.mass;
                        
                        S_Sml_f(:,DY) = Sml_f.bio;
                        S_Sml_p(:,DY) = Sml_p.bio;
                        S_Sml_d(:,DY) = Sml_d.bio;
                        S_Med_f(:,DY) = Med_f.bio;
                        S_Med_p(:,DY) = Med_p.bio;
                        S_Med_d(:,DY) = Med_d.bio;
                        S_Lrg_p(:,DY) = Lrg_p.bio;
                        S_Lrg_d(:,DY) = Lrg_d.bio;
                        
                        S_Med_f_fish(:,DY) = Med_f.caught;
                        S_Med_p_fish(:,DY) = Med_p.caught;
                        S_Med_d_fish(:,DY) = Med_d.caught;
                        S_Lrg_p_fish(:,DY) = Lrg_p.caught;
                        S_Lrg_d_fish(:,DY) = Lrg_d.caught;
                    end
                    
                end %Days
                
            end %Years
            
            %! Calculate monthly means and save
            aa = (cumsum(MNTH)+1);
            a = [1,aa(1:end-1)]; % start of the month
            b = cumsum(MNTH); % end of the month
            for i = 1:12
                %! Put vars of netcdf file
                Spinup_Bent.bio(:,i) = mean(S_Bent_bio(:,a(i):b(i)),2);
                Spinup_Sml_f.bio(:,i) = mean(S_Sml_f(:,a(i):b(i)),2);
                Spinup_Sml_p.bio(:,i) = mean(S_Sml_p(:,a(i):b(i)),2);
                Spinup_Sml_d.bio(:,i) = mean(S_Sml_d(:,a(i):b(i)),2);
                Spinup_Med_f.bio(:,i) = mean(S_Med_f(:,a(i):b(i)),2);
                Spinup_Med_p.bio(:,i) = mean(S_Med_p(:,a(i):b(i)),2);
                Spinup_Med_d.bio(:,i) = mean(S_Med_d(:,a(i):b(i)),2);
                Spinup_Lrg_p.bio(:,i) = mean(S_Lrg_p(:,a(i):b(i)),2);
                Spinup_Lrg_d.bio(:,i) = mean(S_Lrg_d(:,a(i):b(i)),2);
                
                Spinup_Med_f.yield(:,i) = mean(S_Med_f_fish(:,a(i):b(i)),2);
                Spinup_Med_p.yield(:,i) = mean(S_Med_p_fish(:,a(i):b(i)),2);
                Spinup_Med_d.yield(:,i) = mean(S_Med_d_fish(:,a(i):b(i)),2);
                Spinup_Lrg_p.yield(:,i) = mean(S_Lrg_p_fish(:,a(i):b(i)),2);
                Spinup_Lrg_d.yield(:,i) = mean(S_Lrg_d_fish(:,a(i):b(i)),2);
                
            end %Monthly mean
            %%% Save
            if harv==1
                save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,...
                    '/Clim_means_fish_lF',tF(2:end),'_qP',tP(2:end),'_qD',tD(2:end),'.mat'],...
                    'Spinup_Sml_f','Spinup_Sml_p','Spinup_Sml_d','Spinup_Med_f',...
                    'Spinup_Med_p','Spinup_Med_d','Spinup_Lrg_p','Spinup_Lrg_d',...
                    'Spinup_Bent')
            else
                save(['/Volumes/GFDL/CSV/Matlab_new_size/',simname,'/Clim_means.mat'],...
                    'Spinup_Sml_f','Spinup_Sml_p','Spinup_Sml_d','Spinup_Med_f',...
                    'Spinup_Med_p','Spinup_Med_d','Spinup_Lrg_p','Spinup_Lrg_d',...
                    'Spinup_Bent')
            end
            
        end
    end
end

end
