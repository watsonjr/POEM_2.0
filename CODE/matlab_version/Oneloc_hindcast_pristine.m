%%%%!! RUN HINDCAST FOR ONE LOCATION
function Oneloc_hindcast_pristine()
	%! Make parameters
	make_parameters(0) % make core parameters/constants
	%! Setup
	%! Load COBALT and grid data
	Tref = readdlm('./Data/grid_phenol_T0raw_NOflip.csv',','); %min temp for each yr at each location
	%global TrefP = readdlm('./Data/grid_phenol_T0p_clim_min_NOflip.csv',','); %1901-1950 climatological min temp at each location for upper 100m
	%global TrefB = readdlm('./Data/grid_phenol_T0b_clim_min_NOflip.csv',','); %1901-1950 climatological min temp at each location for bottom
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm('./Data/grid_phenol_DTraw_NOflip.csv',',');
	global Sp = readdlm('./Data/Gaussian_spawn_2mo.csv',',');
	global GRD = load('./Data/Data_grid_hindcast_NOTflipped.jld')

	XY = zeros(Int,360,200);
  XY(GRD('ID')) = collect(1:GRD('N'))
	%ID = 40319 %30181 % Georges Bank
  %ID = 42639 %15105 % Eastern Bering Sea
  %ID = 41782 %19526 % Ocean Station Papa
  %ID = 36334 %17377 % Hawaii
	%ID = 38309 %30335 % Bermuda
  %ID = 42744 %40403 % North Sea
	ids = (40319,42639,41782,36334,38309,42744,30051)
	names = ('GB','EBS','OSP','HOT','BATS','NS','EEP')

	for L = 1:7
		ID = ids(L)
		loc = names(L)

		const global NX = length(ID)
		const YEARS = 145; % integration period in years
		const global DAYS = 365; % number of days
		const global MNTH = collect((31,28,31,30,31,30,31,31,30,31,30,31)) % days in month
		%! Initialize
		phen=1;
		Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
		Med_d.td(1) = 0.0;
		Lrg_d.td(1) = 0.0;
		ENVR = sub_init_env(ID);
		%READ INITIAL VALUES FROM PRE-INDUSTRIAL RUN

		%! Storage
		if (phen==1)
			Oneloc_hist_Sml_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Sml_f.csv'),'w')
			Oneloc_hist_Sml_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Sml_p.csv'),'w')
			Oneloc_hist_Sml_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Sml_d.csv'),'w')
			Oneloc_hist_Med_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Med_f.csv'),'w')
			Oneloc_hist_Med_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Med_p.csv'),'w')
			Oneloc_hist_Med_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Med_d.csv'),'w')
			Oneloc_hist_Lrg_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Lrg_p.csv'),'w')
			Oneloc_hist_Lrg_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Lrg_d.csv'),'w')
			Oneloc_hist_Cobalt = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_phen_',loc,'_Cobalt.csv'),'w')
		else
			Oneloc_hist_Sml_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Sml_f.csv'),'w')
			Oneloc_hist_Sml_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Sml_p.csv'),'w')
			Oneloc_hist_Sml_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Sml_d.csv'),'w')
			Oneloc_hist_Med_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Med_f.csv'),'w')
			Oneloc_hist_Med_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Med_p.csv'),'w')
			Oneloc_hist_Med_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Med_d.csv'),'w')
			Oneloc_hist_Lrg_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Lrg_p.csv'),'w')
			Oneloc_hist_Lrg_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Lrg_d.csv'),'w')
			Oneloc_hist_Cobalt = open(string('/Volumes/GFDL/CSV/',simname, '/Oneloc_hist_',loc,'_Cobalt.csv'),'w')
		end

		%%%%%%%%%%%%%%%%%% RUN MODEL
		%! Iterate Model forward in time
		for YR = 1:YEARS % years
			%! Load a year's COBALT data
			ti = string(1860+YR)
			COBALT = load(string('/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_',ti(1:end),'.jld'));
			%reset spawning flag
			if (phen == 1)
				Med_f.S = zeros(Float64,NX,DAYS)
				Lrg_d.S = zeros(Float64,NX,DAYS)
				Lrg_p.S = zeros(Float64,NX,DAYS)
			end
			for DAY = 1:DT:DAYS % days
				%%%! ticker
				DY  = Int(ceil(DAY))
				println(YR,' , ', mod(DY,365))
				%%%! Future time step
				sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
				DY+=1

				%%%! Daily storage
				%! Save
				writecsv(Oneloc_hist_Sml_f,(Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S(DY-1) Sml_f.prod Sml_f.pred Sml_f.nmort Sml_f.met Sml_f.caught))
				writecsv(Oneloc_hist_Sml_p,(Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S(DY-1) Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught))
				writecsv(Oneloc_hist_Sml_d,(Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S(DY-1) Sml_d.prod Sml_d.pred Sml_d.nmort Sml_d.met Sml_d.caught))
				writecsv(Oneloc_hist_Med_f,(Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S(DY-1) Med_f.prod Med_f.pred Med_f.nmort Med_f.met Med_f.caught))
				writecsv(Oneloc_hist_Med_p,(Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S(DY-1) Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught))
				writecsv(Oneloc_hist_Med_d,(Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S(DY-1) Med_d.prod Med_d.pred Med_d.nmort Med_d.met Med_d.caught))
				writecsv(Oneloc_hist_Lrg_p,(Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S(DY-1) Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught))
				writecsv(Oneloc_hist_Lrg_d,(Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S(DY-1) Lrg_d.prod Lrg_d.pred Lrg_d.nmort Lrg_d.met Lrg_d.caught))
				writecsv(Oneloc_hist_Cobalt,(BENT.mass ENVR.fZm ENVR.fZl ENVR.fB))

			end %Days
		end %Years
		%%% close save
	  close(Oneloc_hist_Sml_f)
	  close(Oneloc_hist_Sml_p)
	  close(Oneloc_hist_Sml_d)
	  close(Oneloc_hist_Med_f)
	  close(Oneloc_hist_Med_p)
	  close(Oneloc_hist_Med_d)
	  close(Oneloc_hist_Lrg_p)
		close(Oneloc_hist_Lrg_d)
		close(Oneloc_hist_Cobalt)
	end %Locations
end
