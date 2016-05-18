
####!! RUN SPINUP FOR ONE LOCATION
function Testoneloc()

	#! Make parameters
	make_parameters(0) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
  COBALT = load("./Data/JLD/Data_hindcast_molCm2_000120.jld"); # 1980
	#! Add phenology params from csv file with ID as row
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	YEARS = 100
  DAYS = 365

	#! choose where to run the model
	global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
	#XY = zeros(Int,200,360); # choose a particulat place or everywhere
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii OT
	#ID = 38309 #30335 # Bermuda ATS
  #ID = 42744 #40403 # North Sea
	const global NX = length(ID)
	phen=1;
	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	Med_d.td[1] = 0.0;
	Lrg_d.td[1] = 0.0;
	ENVR = sub_init_env(ID);

	#! Storage
	if (phen==1)
		Spinup_Sml_f  = open("./Data/CSV/Spinup_phen_BATS_Sml_f.csv","w")
		Spinup_Sml_p  = open("./Data/CSV/Spinup_phen_BATS_Sml_p.csv","w")
		Spinup_Sml_d  = open("./Data/CSV/Spinup_phen_BATS_Sml_d.csv","w")
		Spinup_Med_f  = open("./Data/CSV/Spinup_phen_BATS_Med_f.csv","w")
		Spinup_Med_p  = open("./Data/CSV/Spinup_phen_BATS_Med_p.csv","w")
		Spinup_Med_d  = open("./Data/CSV/Spinup_phen_BATS_Med_d.csv","w")
		Spinup_Lrg_p  = open("./Data/CSV/Spinup_phen_BATS_Lrg_p.csv","w")
		Spinup_Lrg_d  = open("./Data/CSV/Spinup_phen_BATS_Lrg_d.csv","w")
		Spinup_Cobalt = open("./Data/CSV/Spinup_phen_BATS_Cobalt.csv","w")
	else
		Spinup_Sml_f  = open("./Data/CSV/Spinup_BATS_Sml_f.csv","w")
		Spinup_Sml_p  = open("./Data/CSV/Spinup_BATS_Sml_p.csv","w")
		Spinup_Sml_d  = open("./Data/CSV/Spinup_BATS_Sml_d.csv","w")
		Spinup_Med_f  = open("./Data/CSV/Spinup_BATS_Med_f.csv","w")
		Spinup_Med_p  = open("./Data/CSV/Spinup_BATS_Med_p.csv","w")
		Spinup_Med_d  = open("./Data/CSV/Spinup_BATS_Med_d.csv","w")
		Spinup_Lrg_p  = open("./Data/CSV/Spinup_BATS_Lrg_p.csv","w")
		Spinup_Lrg_d  = open("./Data/CSV/Spinup_BATS_Lrg_d.csv","w")
		Spinup_Cobalt = open("./Data/CSV/Spinup_BATS_Cobalt.csv","w")
	end

	#! Iterate forward in time with NO fishing
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days

			###! ticker
			DY  = Int(ceil(DAY))
			println(YR," , ", mod(DY,365))

			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Tref,Dthresh,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
			DY+=1

			#! Save
			writecsv(Spinup_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1]])
			writecsv(Spinup_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1]])
			writecsv(Spinup_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1]])
			writecsv(Spinup_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1]])
			writecsv(Spinup_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1]])
			writecsv(Spinup_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1]])
			writecsv(Spinup_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1]])
			writecsv(Spinup_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1]])
			writecsv(Spinup_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl])
		end
	end

	### close save
  close(Spinup_Sml_f)
  close(Spinup_Sml_p)
  close(Spinup_Sml_d)
  close(Spinup_Med_f)
  close(Spinup_Med_p)
  close(Spinup_Med_d)
  close(Spinup_Lrg_p)
	close(Spinup_Lrg_d)
	close(Spinup_Cobalt)

end


####!! RUN HINDCAST FOR ONE LOCATION
function Oneloc_hindcast_pristine()
	#! Make parameters
	make_parameters(0) # make core parameters/constants
	#! Setup
	#! Load COBALT and grid data
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii
	ID = 38309 #30335 # Bermuda
  #ID = 42744 #40403 # North Sea
	const global NX = length(ID)
	const global YEARS = 145; # integration period in years
	const global DAYS = 365; # number of days
	const global MNTH = collect([31,28,31,30,31,30,31,31,30,31,30,31]) # days in month
	#! Initialize
	phen=1;
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	Med_d.td[1] = 0.0;
	Lrg_d.td[1] = 0.0;
	ENVR = sub_init_env(ID);
	#! Storage
	if (phen==1)
		Oneloc_hist_Sml_f  = open("./Data/CSV/Oneloc_hist_phen_BATS_Sml_f.csv","w")
		Oneloc_hist_Sml_p  = open("./Data/CSV/Oneloc_hist_phen_BATS_Sml_p.csv","w")
		Oneloc_hist_Sml_d  = open("./Data/CSV/Oneloc_hist_phen_BATS_Sml_d.csv","w")
		Oneloc_hist_Med_f  = open("./Data/CSV/Oneloc_hist_phen_BATS_Med_f.csv","w")
		Oneloc_hist_Med_p  = open("./Data/CSV/Oneloc_hist_phen_BATS_Med_p.csv","w")
		Oneloc_hist_Med_d  = open("./Data/CSV/Oneloc_hist_phen_BATS_Med_d.csv","w")
		Oneloc_hist_Lrg_p  = open("./Data/CSV/Oneloc_hist_phen_BATS_Lrg_p.csv","w")
		Oneloc_hist_Lrg_d  = open("./Data/CSV/Oneloc_hist_phen_BATS_Lrg_d.csv","w")
		Oneloc_hist_Cobalt = open("./Data/CSV/Oneloc_hist_phen_BATS_Cobalt.csv","w")
	else
		Oneloc_hist_Sml_f  = open("./Data/CSV/Oneloc_hist_BATS_Sml_f.csv","w")
		Oneloc_hist_Sml_p  = open("./Data/CSV/Oneloc_hist_BATS_Sml_p.csv","w")
		Oneloc_hist_Sml_d  = open("./Data/CSV/Oneloc_hist_BATS_Sml_d.csv","w")
		Oneloc_hist_Med_f  = open("./Data/CSV/Oneloc_hist_BATS_Med_f.csv","w")
		Oneloc_hist_Med_p  = open("./Data/CSV/Oneloc_hist_BATS_Med_p.csv","w")
		Oneloc_hist_Med_d  = open("./Data/CSV/Oneloc_hist_BATS_Med_d.csv","w")
		Oneloc_hist_Lrg_p  = open("./Data/CSV/Oneloc_hist_BATS_Lrg_p.csv","w")
		Oneloc_hist_Lrg_d  = open("./Data/CSV/Oneloc_hist_BATS_Lrg_d.csv","w")
		Oneloc_hist_Cobalt = open("./Data/CSV/Oneloc_hist_BATS_Cobalt.csv","w")
	end

	################## RUN MODEL
	#! Iterate Model forward in time
	for YR = 1:YEARS # years
		#! Load a year's COBALT data
		ti = string(YR+1000000)
		COBALT = load(string("./Data/JLD/Data_hindcast_molCm2_",ti[2:end],".jld"));
		#reset spawning flag
		if (phen == 1)
			Med_f.S = zeros(Float64,NX,DAYS)
			Med_d.S = zeros(Float64,NX,DAYS)
			Lrg_p.S = zeros(Float64,NX,DAYS)
		end
		for DAY = 1:DT:DAYS # days
			###! ticker
			DY  = Int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Tref,Dthresh,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
			DY+=1

			###! Daily storage
			#! Save
			writecsv(Oneloc_hist_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1]])
			writecsv(Oneloc_hist_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1]])
			writecsv(Oneloc_hist_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1]])
			writecsv(Oneloc_hist_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1]])
			writecsv(Oneloc_hist_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1]])
			writecsv(Oneloc_hist_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1]])
			writecsv(Oneloc_hist_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1]])
			writecsv(Oneloc_hist_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1]])
			writecsv(Oneloc_hist_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl])

		end
	end
	### close save
  close(Oneloc_hist_Sml_f)
  close(Oneloc_hist_Sml_p)
  close(Oneloc_hist_Sml_d)
  close(Oneloc_hist_Med_f)
  close(Oneloc_hist_Med_p)
  close(Oneloc_hist_Med_d)
  close(Oneloc_hist_Lrg_p)
	close(Oneloc_hist_Lrg_d)
	close(Oneloc_hist_Cobalt)
end



####!! RUN FORECAST FOR ONE LOCATION
function Oneloc_forecast_pristine()
	#! Make parameters
	make_parameters(0) # make core parameters/constants
	#! Setup
	#! Load COBALT and grid data
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global GRD = load("./Data/Data_grid_forecast_NOTflipped.jld")
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii OS
	#ID = 38309 #30335 # Bermuda ATS
  #ID = 42744 #40403 # North Sea
	const global NX = length(ID)
	const global YEARS = 95; # integration period in years
	const global DAYS = 365; # number of days
	const global MNTH = collect([31,28,31,30,31,30,31,31,30,31,30,31]) # days in month
	#! Initialize
	phen=0;
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	Med_d.td[1] = 0.0;
	Lrg_d.td[1] = 0.0;
	ENVR = sub_init_env(ID);
	#! Storage
	if (phen==1)
		Oneloc_fore_Sml_f  = open("./Data/CSV/Oneloc_fore_phen_BATS_Sml_f.csv","w")
		Oneloc_fore_Sml_p  = open("./Data/CSV/Oneloc_fore_phen_BATS_Sml_p.csv","w")
		Oneloc_fore_Sml_d  = open("./Data/CSV/Oneloc_fore_phen_BATS_Sml_d.csv","w")
		Oneloc_fore_Med_f  = open("./Data/CSV/Oneloc_fore_phen_BATS_Med_f.csv","w")
		Oneloc_fore_Med_p  = open("./Data/CSV/Oneloc_fore_phen_BATS_Med_p.csv","w")
		Oneloc_fore_Med_d  = open("./Data/CSV/Oneloc_fore_phen_BATS_Med_d.csv","w")
		Oneloc_fore_Lrg_p  = open("./Data/CSV/Oneloc_fore_phen_BATS_Lrg_p.csv","w")
		Oneloc_fore_Lrg_d  = open("./Data/CSV/Oneloc_fore_phen_BATS_Lrg_d.csv","w")
		Oneloc_fore_Cobalt = open("./Data/CSV/Oneloc_fore_phen_BATS_Cobalt.csv","w")
	else
		Oneloc_fore_Sml_f  = open("./Data/CSV/Oneloc_fore_BATS_Sml_f.csv","w")
		Oneloc_fore_Sml_p  = open("./Data/CSV/Oneloc_fore_BATS_Sml_p.csv","w")
		Oneloc_fore_Sml_d  = open("./Data/CSV/Oneloc_fore_BATS_Sml_d.csv","w")
		Oneloc_fore_Med_f  = open("./Data/CSV/Oneloc_fore_BATS_Med_f.csv","w")
		Oneloc_fore_Med_p  = open("./Data/CSV/Oneloc_fore_BATS_Med_p.csv","w")
		Oneloc_fore_Med_d  = open("./Data/CSV/Oneloc_fore_BATS_Med_d.csv","w")
		Oneloc_fore_Lrg_p  = open("./Data/CSV/Oneloc_fore_BATS_Lrg_p.csv","w")
		Oneloc_fore_Lrg_d  = open("./Data/CSV/Oneloc_fore_BATS_Lrg_d.csv","w")
		Oneloc_fore_Cobalt = open("./Data/CSV/Oneloc_fore_BATS_Cobalt.csv","w")
	end

	################## RUN MODEL
	#! Iterate Model forward in time
	for YR = 1:YEARS # years
		#! Load a year's COBALT data
		ti = string(YR+1000000)
		COBALT = load(string("./Data/JLD/Data_forecast_molCm2_",ti[2:end],".jld"));
		#reset spawning flag
		if (phen == 1)
			Med_f.S = ones(Float64,NX,DAYS)
			Med_d.S = ones(Float64,NX,DAYS)
			Lrg_p.S = ones(Float64,NX,DAYS)
		end
		for DAY = 1:DT:DAYS # days
			###! ticker
			DY  = Int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Tref,Dthresh,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
			DY+=1
			###! Daily storage
			#! Save
			writecsv(Oneloc_fore_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1]])
			writecsv(Oneloc_fore_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1]])
			writecsv(Oneloc_fore_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1]])
			writecsv(Oneloc_fore_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1]])
			writecsv(Oneloc_fore_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1]])
			writecsv(Oneloc_fore_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1]])
			writecsv(Oneloc_fore_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1]])
			writecsv(Oneloc_fore_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1]])
			writecsv(Oneloc_fore_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl])

		end
	end
	### close save
  close(Oneloc_fore_Sml_f)
  close(Oneloc_fore_Sml_p)
  close(Oneloc_fore_Sml_d)
  close(Oneloc_fore_Med_f)
  close(Oneloc_fore_Med_p)
  close(Oneloc_fore_Med_d)
  close(Oneloc_fore_Lrg_p)
	close(Oneloc_fore_Lrg_d)
	close(Oneloc_fore_Cobalt)
end



####!! RUN SPINUP FOR ALL LOCATION
function Spinup_pristine()

	############### Initialize Model Variables
	#! Make parameters
	make_parameters(0) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # if on laptop
	#COBALT = load("./Data/JLD/Data_000001.jld"); # if at school

	#! choose where and when to run the model
	const global YEARS = 20; # integration period in years
	const global NX = 48111
	const global ID = collect(1:NX)

	#! Storage variables
	S_Sml_f = zeros(NX,1)
	S_Sml_p = zeros(NX,1)
	S_Sml_d = zeros(NX,1)
	S_Med_f = zeros(NX,1)
	S_Med_p = zeros(NX,1)
	S_Med_d = zeros(NX,1)
	S_Lrg_p = zeros(NX,1)

	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);

	############### Setup NetCDF save
	#! Init netcdf file for storage
	#varatts = {"longname" => "Biomass","units" => "kg/m^2"}
	#X_atts = {"longname" => "Space", "units" => "grid cell"}
	#timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}
	#Use "Dict{Any,Any}(a=>b, ...)" instead.
	varatts = Dict("longname" => "Biomass",
           "units"    => "kg/m^2")
	X_atts = Dict("longname" => "Space",
			"units"    => "grid cell")
	timatts = Dict("longname" => "Time",
			"units"    => "hours since 01-01-2000 00:00:00")

	#! Init dims of netcdf file
	X   = collect(1:NX)
	tim = collect(1)

	#! setup netcdf path to store to
	file_sml_f = "./Data/NC/Data_spinup_pristine_sml_f.nc"
	file_sml_p = "./Data/NC/Data_spinup_pristine_sml_p.nc"
	file_sml_d = "./Data/NC/Data_spinup_pristine_sml_d.nc"
	file_med_f = "./Data/NC/Data_spinup_pristine_med_f.nc"
	file_med_p = "./Data/NC/Data_spinup_pristine_med_p.nc"
	file_med_d = "./Data/NC/Data_spinup_pristine_med_d.nc"
	file_lrg_p = "./Data/NC/Data_spinup_pristine_lrg_p.nc"

	#! remove if already in existence
	isfile(file_sml_f) ? rm(file_sml_f) : nothing
	isfile(file_sml_p) ? rm(file_sml_p) : nothing
	isfile(file_sml_d) ? rm(file_sml_d) : nothing
	isfile(file_med_f) ? rm(file_med_f) : nothing
	isfile(file_med_p) ? rm(file_med_p) : nothing
	isfile(file_med_d) ? rm(file_med_d) : nothing
	isfile(file_lrg_p) ? rm(file_lrg_p) : nothing

	#! create netcdf files
	nccreate(file_sml_f,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_sml_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_sml_d,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_med_f,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_med_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_med_d,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_lrg_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=varatts)

	#! Initializing netcdf files
	println("Initializing file system (takes about 2 minutes)")
	ncwrite(zeros(NX,1),file_sml_f,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"biomass",[1,1])

	###################### Run the Model
	#! Run model with no fishing
	for YR = 1:YEARS # years

		for DAY = 1:DT:365 # days

			###! Future time step
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

		end

	end

	##################### Clean up
	#! Store
	for i = 1:NX
		S_Sml_f[i,1] = Sml_f.bio[i]
		S_Sml_p[i,1] = Sml_p.bio[i]
		S_Sml_d[i,1] = Sml_d.bio[i]
		S_Med_f[i,1] = Med_f.bio[i]
		S_Med_p[i,1] = Med_p.bio[i]
		S_Med_d[i,1] = Med_d.bio[i]
		S_Lrg_p[i,1] = Lrg_p.bio[i]
	end

	#! Save
	ncwrite(S_Sml_f,file_sml_f,"biomass",[1,1])
	ncwrite(S_Sml_p,file_sml_p,"biomass",[1,1])
	ncwrite(S_Sml_d,file_sml_d,"biomass",[1,1])
	ncwrite(S_Med_f,file_med_f,"biomass",[1,1])
	ncwrite(S_Med_p,file_med_p,"biomass",[1,1])
	ncwrite(S_Med_d,file_med_d,"biomass",[1,1])
	ncwrite(S_Lrg_p,file_lrg_p,"biomass",[1,1])

	#! Close save
  ncclose(file_sml_f)
  ncclose(file_sml_p)
  ncclose(file_sml_d)
  ncclose(file_med_f)
  ncclose(file_med_p)
  ncclose(file_med_d)
  ncclose(file_lrg_p)
end







####!! RUN SPINUP FOR ALL LOCATION
function Spinup_fished()

	############### Initialize Model Variables
	#! Make parameters
	make_parameters(1) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # if on laptop
	#COBALT = load("./Data/JLD/Data_000001.jld"); # if at school

	#! choose where and when to run the model
	const global YEARS = 1; # integration period in years
	const global NX = 48111
	const global ID = collect(1:NX)

	#! Storage variables
	S_PISC = zeros(48111,PI_N,1)
	S_PLAN = zeros(48111,PL_N,1)
	S_DETR = zeros(48111,DE_N,1)
	S_BENT = zeros(48111,BE_N,1)

	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	ENVR = sub_init_env(ID);


	############### Setup NetCDF save
	#! Init netcdf file for storage
	#varatts = {"longname" => "Biomass","units" => "kg/m^2"}
	#X_atts = {"longname" => "Space", "units" => "grid cell"}
	#S_atts = {"longname" => "Size classes", "units"  => "g"}
	#timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}
	#Use "Dict{Any,Any}(a=>b, ...)" instead.
	varatts = Dict("longname" => "Biomass",
           "units"    => "kg/m^2")
	X_atts = Dict("longname" => "Space",
			"units"    => "grid cell")
	S_atts = Dict("longname" => "Size classes",
			"units"    => "g")
	timatts = Dict("longname" => "Time",
			"units"    => "hours since 01-01-2000 00:00:00")

	#! Init dims of netcdf file
	S_pi=collect(1:PI_N)
	S_pl=collect(1:PL_N)
	S_de=collect(1:DE_N)
	S_be=collect(1:BE_N)
	X=collect(1:NX)
	tim=collect(1)

	#! setup netcdf path to store to
	file_pisc = "./Data/NC/Data_spinup_fished_pisc.nc"
	file_plan = "./Data/NC/Data_spinup_fished_plan.nc"
	file_detr = "./Data/NC/Data_spinup_fished_detr.nc"
	file_bent = "./Data/NC/Data_spinup_fished_bent.nc"

	#! remove if already in existence
	isfile(file_pisc) ? rm(file_pisc) : nothing
	isfile(file_plan) ? rm(file_plan) : nothing
	isfile(file_detr) ? rm(file_detr) : nothing
	isfile(file_bent) ? rm(file_bent) : nothing

	#! create netcdf files
	nccreate(file_pisc,"biomass","X",X,X_atts,"S",S_pi,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_plan,"biomass","X",X,X_atts,"S",S_pl,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_detr,"biomass","X",X,X_atts,"S",S_de,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_bent,"biomass","X",X,X_atts,"S",S_be,S_atts,"time",tim,timatts,atts=varatts)

	#! Initializing netcdf files
	println("Initializing file system (takes about 2 minutes)")
	ncwrite(zeros(NX,PI_N,1),file_pisc,"biomass",[1,1,1])
	ncwrite(zeros(NX,PL_N,1),file_plan,"biomass",[1,1,1])
	ncwrite(zeros(NX,DE_N,1),file_detr,"biomass",[1,1,1])
	ncwrite(zeros(NX,BE_N,1),file_bent,"biomass",[1,1,1])



	###################### Run the Model
	#! Run model with no fishing
	for YR = 1:YEARS # years

		for DAY = 1:DT:365 # days

			###! Future time step
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

		end

	end


	##################### Clean up
	#! Store
	for i = 1:NX
		S_PISC[i,:,1] = PISC.bio[i]'
		S_PLAN[i,:,1] = PLAN.bio[i]'
		S_DETR[i,:,1] = DETR.bio[i]'
		S_BENT[i,:,1] = BENT.bio[i]'
	end

	#! Save
	ncwrite(S_PISC,file_pisc,"biomass",[1,1,1])
	ncwrite(S_PLAN,file_plan,"biomass",[1,1,1])
	ncwrite(S_DETR,file_detr,"biomass",[1,1,1])
	ncwrite(S_BENT,file_bent,"biomass",[1,1,1])

	#! Close save
    ncclose(file_pisc)
    ncclose(file_plan)
    ncclose(file_detr)
    ncclose(file_bent)

end





####!! RUN FORECAST FOR ALL LOCATIONS
function Forecast_pristine()

	################ SETUP
	#! Load COBALT and grid data
	GRD = load("./Data/Data_grid.jld"); # spatial information

	#! Make parameters
	make_parameters(0) # make core parameters/constants
	const global YEARS = 2; # integration period in years
	const global DAYS = 365; # number of days
	const global NX = 48111
	const global ID = collect(1:NX)
	const global MNTH = collect([31,28,31,30,31,30,31,31,30,31,30,31]) # days in month

	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	ENVR = sub_init_env(ID);
	pisc = ncread("./Data/NC/Data_spinup_pristine_pisc.nc","biomass")
	plan = ncread("./Data/NC/Data_spinup_pristine_plan.nc","biomass")
	detr = ncread("./Data/NC/Data_spinup_pristine_detr.nc","biomass")
	bent = ncread("./Data/NC/Data_spinup_pristine_bent.nc","biomass")
	for i = 1:NX
		PISC.bio[i] = squeeze(pisc[i,:],1)
		PLAN.bio[i] = squeeze(plan[i,:],1)
		DETR.bio[i] = squeeze(detr[i,:],1)
		BENT.bio[i] = squeeze(bent[i,:],1)
	end


	################ STORAGE
	#! Storage arrays (daily)
	DAY_PISC = zeros(NX,PI_N,DAYS);
	DAY_PLAN = zeros(NX,PL_N,DAYS);
	DAY_DETR = zeros(NX,DE_N,DAYS);
	DAY_BENT = zeros(NX,BE_N,DAYS);

	#! Init netcdf file for storage
	#Use "Dict{Any,Any}(a=>b, ...)" instead.
	varatts = Dict("longname" => "Biomass",
           "units"    => "kg/m^2")
	X_atts = Dict("longname" => "Space",
			"units"    => "grid cell")
	S_atts = Dict("longname" => "Size classes",
			"units"    => "g")
	timatts = Dict("longname" => "Time",
			"units"    => "hours since 01-01-2000 00:00:00")

	#! Init dims of netcdf file
	S_pi=collect([1:PI_N]);	S_pl=collect([1:PL_N]) ; S_de=collect([1:DE_N]);
	S_be=collect([1:BE_N]); X=collect([1:NX]); tim=collect([1:12*YEARS])

	#! setup netcdf path to store to
	file_pisc = "./Data/NC/Data_forecast_pristine_pisc_adv.nc"
	file_plan = "./Data/NC/Data_forecast_pristine_plan_adv.nc"
	file_detr = "./Data/NC/Data_forecast_pristine_detr_adv.nc"
	file_bent = "./Data/NC/Data_forecast_pristine_bent_adv.nc"

	#! remove if already in existence
	isfile(file_pisc) ? rm(file_pisc) : nothing
	isfile(file_plan) ? rm(file_plan) : nothing
	isfile(file_detr) ? rm(file_detr) : nothing
	isfile(file_bent) ? rm(file_bent) : nothing

	#! create netcdf files
	nccreate(file_pisc,"biomass","X",X,X_atts,"S",S_pi,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_plan,"biomass","X",X,X_atts,"S",S_pl,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_detr,"biomass","X",X,X_atts,"S",S_de,S_atts,"time",tim,timatts,atts=varatts)
	nccreate(file_bent,"biomass","X",X,X_atts,"S",S_be,S_atts,"time",tim,timatts,atts=varatts)

	#! Initializing netcdf files
	println("Initializing file system (takes about 2 minutes)")
	ncwrite(zeros(NX,PI_N,1),file_pisc,"biomass",[1,1,1])
	ncwrite(zeros(NX,PL_N,1),file_plan,"biomass",[1,1,1])
	ncwrite(zeros(NX,DE_N,1),file_detr,"biomass",[1,1,1])
	ncwrite(zeros(NX,1,1),file_bent,"biomass",[1,1,1])


	################## RUN MODEL
	#! Iterate Model forward in time
	MNT = 0; # monthly ticker
	for YR = 1:YEARS # years

		#! Load a year's COBALT data
		ti = string(YR+1000000)
		COBALT = load("./Data/Data_000001.jld"); # testing one year
		#COBALT = load(string("./Data/Data_",ti[2:end],".jld")); #  laptop
		#COBALT = load(string("./Data/JLD/Data_",ti[2:end],".jld")); # src comp

		for DAY = 1:DT:DAYS # days

			###! ticker
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))

			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

			###! Daily storage
			for i = 1:NX
				DAY_PISC[i,:,DAY] = PISC.bio[i]
				DAY_PLAN[i,:,DAY] = PLAN.bio[i]
				DAY_DETR[i,:,DAY] = DETR.bio[i]
				DAY_BENT[i,:,DAY] = BENT.bio[i]
			end

		end

		#! Calculate monthly means and save
		a = [1;(cumsum(MNTH)+1)[1:end-1]] # start of the month
		b = cumsum(MNTH) # end of the month
		for i = 1:12
			MNT += 1 # Update monthly ticker
			ncwrite(mean(DAY_PISC[:,:,a[i]:b[i]],3),file_pisc,"biomass",[1,1,MNT])
			ncwrite(mean(DAY_PLAN[:,:,a[i]:b[i]],3),file_plan,"biomass",[1,1,MNT])
			ncwrite(mean(DAY_DETR[:,:,a[i]:b[i]],3),file_detr,"biomass",[1,1,MNT])
			ncwrite(mean(DAY_BENT[:,:,a[i]:b[i]],3),file_bent,"biomass",[1,1,MNT])
		end

	end

	#! Close save
	ncclose(file_pisc)
	ncclose(file_plan)
	ncclose(file_detr)
	ncclose(file_bent)

end
