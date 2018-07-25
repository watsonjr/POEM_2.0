% POEM output at all locations

clear all
close all

% kays = 0.0405:0.01:0.0916;
% bees = 0.125:0.025:0.25; %bees = 0.1:0.05:0.35;
% for n = 1:length(kays)
%     kt = kays(n);
%     
%     for g = 1:length(bees)
%         bpow = bees(g);
%         tkfn = num2str(100+int64(100*kt));
%         tbfn = num2str(100+int64(100*bpow));
%         
%         cfile = ['Dc_enc70_cmax-metab20_b',tbfn(2:end),'_k',tkfn(2:end),'_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC050_lgRE00100_mdRE00100'];
%         
        cfile = 'Dc_enc70_cmax-metab20_b18_k09_fcrit20_D075_J100_A050_Sm025_nmort1_BE05_CC100_lgRE00100_mdRE00100';
        fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];
        
        harv = '03';
        
        %% SP
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_sml_p.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        SP.bio = biomass;
        clear biomass
        
        % SF
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_sml_f.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        SF.bio = biomass;
        clear biomass
        
        % SD
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_sml_d.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        SD.bio = biomass;
        clear biomass
        
        % MP
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_med_p.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        MP.bio = biomass;
        MP.yield = yield;
        clear biomass yield
        
        %% MF
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_med_f.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        MF.bio = biomass;
        MF.yield = yield;
        clear biomass yield
        
        %% MD
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_med_d.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        MD.bio = biomass;
        MD.yield = yield;
        clear biomass yield
        
        % LP
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_lrg_p.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        LP.bio = biomass;
        LP.yield = yield;
        clear biomass yield
        
        % LD
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_lrg_d.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        LD.bio = biomass;
        LD.yield = yield;
        clear biomass yield
        
        % Benthic material
        ncid = netcdf.open([fpath 'Hindcast_fished' harv '_bent.nc'],'NC_NOWRITE');
        [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
        for i = 1:nvars
            varname = netcdf.inqVar(ncid, i-1);
            eval([ varname ' = netcdf.getVar(ncid,i-1);']);
            eval([ varname '(' varname ' == 99999) = NaN;']);
        end
        netcdf.close(ncid);
        
        BENT.bio = biomass;
        clear biomass
        
        %% Take means
        nt = length(time);
        
        %Time
        sp_tmean=mean(SP.bio,1);
        sf_tmean=mean(SF.bio(:,1:nt),1);
        sd_tmean=mean(SD.bio,1);
        mp_tmean=mean(MP.bio,1);
        mf_tmean=mean(MF.bio,1);
        md_tmean=mean(MD.bio,1);
        lp_tmean=mean(LP.bio,1);
        ld_tmean=mean(LD.bio,1);
        b_tmean=mean(BENT.bio,1);
        
        mf_tmy=mean(MF.yield,1);
        mp_tmy=mean(MP.yield,1);
        md_tmy=mean(MD.yield,1);
        lp_tmy=mean(LP.yield,1);
        ld_tmy=mean(LD.yield,1);
        
        % Last 10 years
        lyr=time((end-(10*12)+1):end);
        sp_mean=mean(SP.bio(:,lyr),2);
        sf_mean=mean(SF.bio(:,lyr),2);
        sd_mean=mean(SD.bio(:,lyr),2);
        mp_mean=mean(MP.bio(:,lyr),2);
        mf_mean=mean(MF.bio(:,lyr),2);
        md_mean=mean(MD.bio(:,lyr),2);
        lp_mean=mean(LP.bio(:,lyr),2);
        ld_mean=mean(LD.bio(:,lyr),2);
        b_mean=mean(BENT.bio(:,lyr),2);
        
        mf_my=mean(MF.yield(:,lyr),2);
        mp_my=mean(MP.yield(:,lyr),2);
        md_my=mean(MD.yield(:,lyr),2);
        lp_my=mean(LP.yield(:,lyr),2);
        ld_my=mean(LD.yield(:,lyr),2);
        
        %
        save([fpath 'Means_Hindcast_fished_' harv '_' cfile '.mat'],...
            'sf_mean','sp_mean','sd_mean','mf_mean','mp_mean','md_mean','b_mean',...
            'lp_mean','ld_mean','sf_tmean','sp_tmean','sd_tmean','mf_tmean','mp_tmean',...
            'md_tmean','b_tmean','lp_tmean','ld_tmean','time','lyr',...
            'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
            'mf_my','mp_my','md_my','lp_my','ld_my');
        
%     end
% end





