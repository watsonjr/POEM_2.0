
####!! RUN SPINUP FOR ONE LOCATION
function Testoneloc()

	#! Make parameters
	make_parameters(0) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
  COBALT = load("./Data/JLD/Data_hindcast_PC_000120.jld"); # 1980
	#! Add phenology params from csv file with ID as row
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	#global TrefP = readdlm("./Data/grid_phenol_T0p_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for upper 100m
	#global TrefB = readdlm("./Data/grid_phenol_T0b_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for bottom
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	YEARS = 100
  DAYS = 365

	#! choose where to run the model
	global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
	#XY = zeros(Int,200,360); # choose a particulat place or everywhere
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii OT
	#ID = 38309 #30335 # Bermuda ATS
  #ID = 42744 #40403 # North Sea
	ids = [40319,42639,41782,36334,38309,42744,30051]
	names = ["GB","EBS","OSP","HOT","BATS","NS","EEP"]

	for L = 1:7
		ID = ids[L]
		loc = names[L]

		const global NX = length(ID)
		phen=0;
		#! Initialize
		Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
		Med_d.td[1] = 0.0;
		Lrg_d.td[1] = 0.0;
		ENVR = sub_init_env(ID);

		#! Storage
		if (phen==1)
			Spinup_Sml_f  = open(string("./Data/CSV/Spinup_phen_",loc,"_Sml_f.csv"),"w")
			Spinup_Sml_p  = open(string("./Data/CSV/Spinup_phen_",loc,"_Sml_p.csv"),"w")
			Spinup_Sml_d  = open(string("./Data/CSV/Spinup_phen_",loc,"_Sml_d.csv"),"w")
			Spinup_Med_f  = open(string("./Data/CSV/Spinup_phen_",loc,"_Med_f.csv"),"w")
			Spinup_Med_p  = open(string("./Data/CSV/Spinup_phen_",loc,"_Med_p.csv"),"w")
			Spinup_Med_d  = open(string("./Data/CSV/Spinup_phen_",loc,"_Med_d.csv"),"w")
			Spinup_Lrg_p  = open(string("./Data/CSV/Spinup_phen_",loc,"_Lrg_p.csv"),"w")
			Spinup_Lrg_d  = open(string("./Data/CSV/Spinup_phen_",loc,"_Lrg_d.csv"),"w")
			Spinup_Cobalt = open(string("./Data/CSV/Spinup_phen_",loc,"_Cobalt.csv"),"w")
		else
			Spinup_Sml_f  = open(string("./Data/CSV/Spinup_",loc,"_Sml_f.csv"),"w")
			Spinup_Sml_p  = open(string("./Data/CSV/Spinup_",loc,"_Sml_p.csv"),"w")
			Spinup_Sml_d  = open(string("./Data/CSV/Spinup_",loc,"_Sml_d.csv"),"w")
			Spinup_Med_f  = open(string("./Data/CSV/Spinup_",loc,"_Med_f.csv"),"w")
			Spinup_Med_p  = open(string("./Data/CSV/Spinup_",loc,"_Med_p.csv"),"w")
			Spinup_Med_d  = open(string("./Data/CSV/Spinup_",loc,"_Med_d.csv"),"w")
			Spinup_Lrg_p  = open(string("./Data/CSV/Spinup_",loc,"_Lrg_p.csv"),"w")
			Spinup_Lrg_d  = open(string("./Data/CSV/Spinup_",loc,"_Lrg_d.csv"),"w")
			Spinup_Cobalt = open(string("./Data/CSV/Spinup_",loc,"_Cobalt.csv"),"w")
		end

		#! Iterate forward in time with NO fishing
		for YR = 1:YEARS # years

			#reset spawning flag
			if (phen == 1)
				Med_f.S = zeros(Float64,NX,DAYS)
				Lrg_d.S = zeros(Float64,NX,DAYS)
				Lrg_p.S = zeros(Float64,NX,DAYS)
			end

			for DAY = 1:DT:DAYS # days

				###! ticker
				DY  = Int(ceil(DAY))
				println(YR," , ", mod(DY,365))

				###! Future time step
				sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
				DY+=1

				#! Save
				writecsv(Spinup_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1] Sml_f.prod])
				writecsv(Spinup_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1] Sml_p.prod])
				writecsv(Spinup_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1] Sml_d.prod])
				writecsv(Spinup_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1] Med_f.prod])
				writecsv(Spinup_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1] Med_p.prod])
				writecsv(Spinup_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1] Med_d.prod])
				writecsv(Spinup_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1] Lrg_p.prod])
				writecsv(Spinup_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1] Lrg_d.prod])
				writecsv(Spinup_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl ENVR.fB])
			end #Days
		end #Years

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

	end #Locations

end


####!! RUN HINDCAST FOR ONE LOCATION
function Oneloc_hindcast_pristine()
	#! Make parameters
	make_parameters(0) # make core parameters/constants
	#! Setup
	#! Load COBALT and grid data
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	#global TrefP = readdlm("./Data/grid_phenol_T0p_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for upper 100m
	#global TrefB = readdlm("./Data/grid_phenol_T0b_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for bottom
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii
	#ID = 38309 #30335 # Bermuda
  #ID = 42744 #40403 # North Sea
	ids = [40319,42639,41782,36334,38309,42744,30051]
	names = ["GB","EBS","OSP","HOT","BATS","NS","EEP"]

	for L = 1:7
		ID = ids[L]
		loc = names[L]

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
			Oneloc_hist_Sml_f  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Sml_f.csv"),"w")
			Oneloc_hist_Sml_p  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Sml_p.csv"),"w")
			Oneloc_hist_Sml_d  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Sml_d.csv"),"w")
			Oneloc_hist_Med_f  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Med_f.csv"),"w")
			Oneloc_hist_Med_p  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Med_p.csv"),"w")
			Oneloc_hist_Med_d  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Med_d.csv"),"w")
			Oneloc_hist_Lrg_p  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Lrg_p.csv"),"w")
			Oneloc_hist_Lrg_d  = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Lrg_d.csv"),"w")
			Oneloc_hist_Cobalt = open(string("./Data/CSV/Oneloc_hist_phen_",loc,"_Cobalt.csv"),"w")
		else
			Oneloc_hist_Sml_f  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Sml_f.csv"),"w")
			Oneloc_hist_Sml_p  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Sml_p.csv"),"w")
			Oneloc_hist_Sml_d  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Sml_d.csv"),"w")
			Oneloc_hist_Med_f  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Med_f.csv"),"w")
			Oneloc_hist_Med_p  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Med_p.csv"),"w")
			Oneloc_hist_Med_d  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Med_d.csv"),"w")
			Oneloc_hist_Lrg_p  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Lrg_p.csv"),"w")
			Oneloc_hist_Lrg_d  = open(string("./Data/CSV/Oneloc_hist_",loc,"_Lrg_d.csv"),"w")
			Oneloc_hist_Cobalt = open(string("./Data/CSV/Oneloc_hist_",loc,"_Cobalt.csv"),"w")
		end

		################## RUN MODEL
		#! Iterate Model forward in time
		for YR = 1:YEARS # years
			#! Load a year's COBALT data
			ti = string(YR+1000000)
			COBALT = load(string("./Data/JLD/Data_hindcast_PC_",ti[2:end],".jld"));
			#reset spawning flag
			if (phen == 1)
				Med_f.S = zeros(Float64,NX,DAYS)
				Lrg_d.S = zeros(Float64,NX,DAYS)
				Lrg_p.S = zeros(Float64,NX,DAYS)
			end
			for DAY = 1:DT:DAYS # days
				###! ticker
				DY  = Int(ceil(DAY))
				println(YR," , ", mod(DY,365))
				###! Future time step
				sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
				DY+=1

				###! Daily storage
				#! Save
				writecsv(Oneloc_hist_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1] Sml_f.prod])
				writecsv(Oneloc_hist_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1] Sml_p.prod])
				writecsv(Oneloc_hist_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1] Sml_d.prod])
				writecsv(Oneloc_hist_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1] Med_f.prod])
				writecsv(Oneloc_hist_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1] Med_p.prod])
				writecsv(Oneloc_hist_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1] Med_d.prod])
				writecsv(Oneloc_hist_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1] Lrg_p.prod])
				writecsv(Oneloc_hist_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1] Lrg_d.prod])
				writecsv(Oneloc_hist_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl ENVR.fB])

			end #Days
		end #Years
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
	end #Locations
end



####!! RUN FORECAST FOR ONE LOCATION
function Oneloc_forecast_pristine()
	#! Make parameters
	make_parameters(0) # make core parameters/constants
	#! Setup
	#! Load COBALT and grid data
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	#global TrefP = readdlm("./Data/grid_phenol_T0p_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for upper 100m
	#global TrefB = readdlm("./Data/grid_phenol_T0b_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for bottom
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	global GRD = load("./Data/Data_grid_forecast_NOTflipped.jld")
	XY = zeros(Int,360,200);
  XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = 40319 #30181 # Georges Bank
  #ID = 42639 #15105 # Eastern Bering Sea
  #ID = 41782 #19526 # Ocean Station Papa
  #ID = 36334 #17377 # Hawaii OS
	#ID = 38309 #30335 # Bermuda ATS
  #ID = 42744 #40403 # North Sea
	ids = [40319,42639,41782,36334,38309,42744,30051]
	names = ["GB","EBS","OSP","HOT","BATS","NS","EEP"]

	for L = 1:7
		ID = ids[L]
		loc = names[L]

		const global NX = length(ID)
		const global YEARS = 95; # integration period in years
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
			Oneloc_fore_Sml_f  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Sml_f.csv"),"w")
			Oneloc_fore_Sml_p  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Sml_p.csv"),"w")
			Oneloc_fore_Sml_d  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Sml_d.csv"),"w")
			Oneloc_fore_Med_f  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Med_f.csv"),"w")
			Oneloc_fore_Med_p  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Med_p.csv"),"w")
			Oneloc_fore_Med_d  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Med_d.csv"),"w")
			Oneloc_fore_Lrg_p  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Lrg_p.csv"),"w")
			Oneloc_fore_Lrg_d  = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Lrg_d.csv"),"w")
			Oneloc_fore_Cobalt = open(string("./Data/CSV/Oneloc_fore_phen_",loc,"_Cobalt.csv"),"w")
		else
			Oneloc_fore_Sml_f  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Sml_f.csv"),"w")
			Oneloc_fore_Sml_p  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Sml_p.csv"),"w")
			Oneloc_fore_Sml_d  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Sml_d.csv"),"w")
			Oneloc_fore_Med_f  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Med_f.csv"),"w")
			Oneloc_fore_Med_p  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Med_p.csv"),"w")
			Oneloc_fore_Med_d  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Med_d.csv"),"w")
			Oneloc_fore_Lrg_p  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Lrg_p.csv"),"w")
			Oneloc_fore_Lrg_d  = open(string("./Data/CSV/Oneloc_fore_",loc,"_Lrg_d.csv"),"w")
			Oneloc_fore_Cobalt = open(string("./Data/CSV/Oneloc_fore_",loc,"_Cobalt.csv"),"w")
		end

		################## RUN MODEL
		#! Iterate Model forward in time
		for YR = 1:YEARS # years
			#! Load a year's COBALT data
			ti = string(YR+1000000)
			COBALT = load(string("./Data/JLD/Data_forecast_PC_",ti[2:end],".jld"));
			#reset spawning flag
			if (phen == 1)
				Med_f.S = zeros(Float64,NX,DAYS)
				Lrg_d.S = zeros(Float64,NX,DAYS)
				Lrg_p.S = zeros(Float64,NX,DAYS)
			end
			for DAY = 1:DT:DAYS # days
				###! ticker
				DY  = Int(ceil(DAY))
				println(YR," , ", mod(DY,365))
				###! Future time step
				sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
				DY+=1
				###! Daily storage
				#! Save
				writecsv(Oneloc_fore_Sml_f,[Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep Sml_f.rec Sml_f.egg Sml_f.clev Sml_f.DD Sml_f.S[DY-1] Sml_f.prod])
				writecsv(Oneloc_fore_Sml_p,[Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.egg Sml_p.clev Sml_p.DD Sml_p.S[DY-1] Sml_p.prod])
				writecsv(Oneloc_fore_Sml_d,[Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.egg Sml_d.clev Sml_d.DD Sml_d.S[DY-1] Sml_d.prod])
				writecsv(Oneloc_fore_Med_f,[Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.egg Med_f.clev Med_f.DD Med_f.S[DY-1] Med_f.prod])
				writecsv(Oneloc_fore_Med_p,[Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.egg Med_p.clev Med_p.DD Med_p.S[DY-1] Med_p.prod])
				writecsv(Oneloc_fore_Med_d,[Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.egg Med_d.clev Med_d.DD Med_d.S[DY-1] Med_d.prod])
				writecsv(Oneloc_fore_Lrg_p,[Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.egg Lrg_p.clev Lrg_p.DD Lrg_p.S[DY-1] Lrg_p.prod])
				writecsv(Oneloc_fore_Lrg_d,[Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.egg Lrg_d.clev Lrg_d.DD Lrg_d.S[DY-1] Lrg_d.prod])
				writecsv(Oneloc_fore_Cobalt,[BENT.mass ENVR.fZm ENVR.fZl ENVR.fB])

			end #Days
		end #Years
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
	end #Locations
end



####!! RUN SPINUP FOR ALL LOCATIONS
function Spinup_pristine()

	############### Initialize Model Variables
	#! Make parameters
	make_parameters(0) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/JLD/Data_hindcast_PC_000120.jld"); # 1980
	#! Add phenology params from csv file with ID as row
	Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
	#global TrefP = readdlm("./Data/grid_phenol_T0p_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for upper 100m
	#global TrefB = readdlm("./Data/grid_phenol_T0b_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for bottom
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
	global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
	global GRD = load("./Data/Data_grid_forecast_NOTflipped.jld")
	YEARS = 2 #50
  DAYS = 365

	#! choose where and when to run the model
	const global NX = 48111
	const global ID = collect(1:NX);

	#! Storage variables
	S_Sml_f = zeros(NX,DAYS);
	S_Sml_p = zeros(NX,DAYS);
	S_Sml_d = zeros(NX,DAYS);
	S_Med_f = zeros(NX,DAYS);
	S_Med_p = zeros(NX,DAYS);
	S_Med_d = zeros(NX,DAYS);
	S_Lrg_p = zeros(NX,DAYS);
	S_Lrg_d = zeros(NX,DAYS);

	S_Sml_f_rec = zeros(NX,DAYS);
	S_Sml_p_rec = zeros(NX,DAYS);
	S_Sml_d_rec = zeros(NX,DAYS);
	S_Med_f_rec = zeros(NX,DAYS);
	S_Med_p_rec = zeros(NX,DAYS);
	S_Med_d_rec = zeros(NX,DAYS);
	S_Lrg_p_rec = zeros(NX,DAYS);
	S_Lrg_d_rec = zeros(NX,DAYS);

	S_Sml_f_con = zeros(NX,DAYS);
	S_Sml_p_con = zeros(NX,DAYS);
	S_Sml_d_con = zeros(NX,DAYS);
	S_Med_f_con = zeros(NX,DAYS);
	S_Med_p_con = zeros(NX,DAYS);
	S_Med_d_con = zeros(NX,DAYS);
	S_Lrg_p_con = zeros(NX,DAYS);
	S_Lrg_d_con = zeros(NX,DAYS);

	S_Sml_f_nu = zeros(NX,DAYS);
	S_Sml_p_nu = zeros(NX,DAYS);
	S_Sml_d_nu = zeros(NX,DAYS);
	S_Med_f_nu = zeros(NX,DAYS);
	S_Med_p_nu = zeros(NX,DAYS);
	S_Med_d_nu = zeros(NX,DAYS);
	S_Lrg_p_nu = zeros(NX,DAYS);
	S_Lrg_d_nu = zeros(NX,DAYS);

	S_Sml_f_prod = zeros(NX,DAYS);
	S_Sml_p_prod = zeros(NX,DAYS);
	S_Sml_d_prod = zeros(NX,DAYS);
	S_Med_f_prod = zeros(NX,DAYS);
	S_Med_p_prod = zeros(NX,DAYS);
	S_Med_d_prod = zeros(NX,DAYS);
	S_Lrg_p_prod = zeros(NX,DAYS);
	S_Lrg_d_prod = zeros(NX,DAYS);

	S_Sml_f_gamma = zeros(NX,DAYS);
	S_Sml_p_gamma = zeros(NX,DAYS);
	S_Sml_d_gamma = zeros(NX,DAYS);
	S_Med_f_gamma = zeros(NX,DAYS);
	S_Med_p_gamma = zeros(NX,DAYS);
	S_Med_d_gamma = zeros(NX,DAYS);
	S_Lrg_p_gamma = zeros(NX,DAYS);
	S_Lrg_d_gamma = zeros(NX,DAYS);

	S_Sml_f_rep = zeros(NX,DAYS);
	S_Sml_p_rep = zeros(NX,DAYS);
	S_Sml_d_rep = zeros(NX,DAYS);
	S_Med_f_rep = zeros(NX,DAYS);
	S_Med_p_rep = zeros(NX,DAYS);
	S_Med_d_rep = zeros(NX,DAYS);
	S_Lrg_p_rep = zeros(NX,DAYS);
	S_Lrg_d_rep = zeros(NX,DAYS);

	S_Sml_f_egg = zeros(NX,DAYS);
	S_Sml_p_egg = zeros(NX,DAYS);
	S_Sml_d_egg = zeros(NX,DAYS);
	S_Med_f_egg = zeros(NX,DAYS);
	S_Med_p_egg = zeros(NX,DAYS);
	S_Med_d_egg = zeros(NX,DAYS);
	S_Lrg_p_egg = zeros(NX,DAYS);
	S_Lrg_d_egg = zeros(NX,DAYS);

	S_Sml_f_die = zeros(NX,DAYS);
	S_Sml_p_die = zeros(NX,DAYS);
	S_Sml_d_die = zeros(NX,DAYS);
	S_Med_f_die = zeros(NX,DAYS);
	S_Med_p_die = zeros(NX,DAYS);
	S_Med_d_die = zeros(NX,DAYS);
	S_Lrg_p_die = zeros(NX,DAYS);
	S_Lrg_d_die = zeros(NX,DAYS);

	S_Sml_f_clev = zeros(NX,DAYS);
	S_Sml_p_clev = zeros(NX,DAYS);
	S_Sml_d_clev = zeros(NX,DAYS);
	S_Med_f_clev = zeros(NX,DAYS);
	S_Med_p_clev = zeros(NX,DAYS);
	S_Med_d_clev = zeros(NX,DAYS);
	S_Lrg_p_clev = zeros(NX,DAYS);
	S_Lrg_d_clev = zeros(NX,DAYS);

	S_Sml_f_S = zeros(NX,DAYS);
	S_Sml_p_S = zeros(NX,DAYS);
	S_Sml_d_S = zeros(NX,DAYS);
	S_Med_f_S = zeros(NX,DAYS);
	S_Med_p_S = zeros(NX,DAYS);
	S_Med_d_S = zeros(NX,DAYS);
	S_Lrg_p_S = zeros(NX,DAYS);
	S_Lrg_d_S = zeros(NX,DAYS);

	S_Sml_f_DD = zeros(NX,DAYS);
	S_Sml_p_DD = zeros(NX,DAYS);
	S_Sml_d_DD = zeros(NX,DAYS);
	S_Med_f_DD = zeros(NX,DAYS);
	S_Med_p_DD = zeros(NX,DAYS);
	S_Med_d_DD = zeros(NX,DAYS);
	S_Lrg_p_DD = zeros(NX,DAYS);
	S_Lrg_d_DD = zeros(NX,DAYS);

	#! Initialize
	phen=0;
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
	Med_d.td[1:NX] = 0.0;
	Lrg_d.td[1:NX] = 0.0;
	ENVR = sub_init_env(ID);

	############### Setup NetCDF save
	# #! Init netcdf file for storage
	# #biomatts = {"longname" => "Biomass","units" => "kg/m^2"}
	# #X_atts = {"longname" => "Space", "units" => "grid cell"}
	# #timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}
	# #Use "Dict{Any,Any}(a=>b, ...)" instead.
	biomatts = Dict("longname" => "Biomass",
	         "units"    => "g/m^2")
	X_atts = Dict("longname" => "Space",
			"units"    => "grid cell")
	timatts = Dict("longname" => "Time",
			"units"    => "days since 01-01-1980 00:00:00")
	specatts = Dict("longname" => "Biomass rate",
			         "units"    => "g/g/day")
	fracatts = Dict("longname" => "Fraction",
			"units"    => "unitless")
	DDatts = Dict("longname" => "Cumulative degree days",
			"units"    => "degrees Celsius")

	# #! Init dims of netcdf file
	X   = collect(1:NX);
	tim = collect(1:DAYS);

	# #! setup netcdf path to store to
	file_sml_f = "/Volumes/GFDL/NC/Data_spinup_pristine_sml_f.nc"
	file_sml_p = "/Volumes/GFDL/NC/Data_spinup_pristine_sml_p.nc"
	file_sml_d = "/Volumes/GFDL/NC/Data_spinup_pristine_sml_d.nc"
	file_med_f = "/Volumes/GFDL/NC/Data_spinup_pristine_med_f.nc"
	file_med_p = "/Volumes/GFDL/NC/Data_spinup_pristine_med_p.nc"
	file_med_d = "/Volumes/GFDL/NC/Data_spinup_pristine_med_d.nc"
	file_lrg_p = "/Volumes/GFDL/NC/Data_spinup_pristine_lrg_p.nc"
	file_lrg_d = "/Volumes/GFDL/NC/Data_spinup_pristine_lrg_d.nc"

	# #! remove if already in existence
	isfile(file_sml_f) ? rm(file_sml_f) : nothing
	isfile(file_sml_p) ? rm(file_sml_p) : nothing
	isfile(file_sml_d) ? rm(file_sml_d) : nothing
	isfile(file_med_f) ? rm(file_med_f) : nothing
	isfile(file_med_p) ? rm(file_med_p) : nothing
	isfile(file_med_d) ? rm(file_med_d) : nothing
	isfile(file_lrg_p) ? rm(file_lrg_p) : nothing
	isfile(file_lrg_d) ? rm(file_lrg_d) : nothing

	nccreate(file_sml_f,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_d,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_f,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_d,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_p,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_d,"biomass","X",X,X_atts,"time",tim,timatts,atts=biomatts);

	nccreate(file_sml_f,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_p,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_d,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_f,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_p,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_d,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_p,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_d,"prod","X",X,X_atts,"time",tim,timatts,atts=biomatts);

	nccreate(file_sml_f,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_p,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_sml_d,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_f,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_p,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_med_d,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_p,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);
	nccreate(file_lrg_d,"rec","X",X,X_atts,"time",tim,timatts,atts=biomatts);

	nccreate(file_sml_f,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"con","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"nu","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"gamma","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"rep","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"egg","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_p,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_sml_d,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_f,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_p,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_med_d,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_p,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);
	nccreate(file_lrg_d,"die","X",X,X_atts,"time",tim,timatts,atts=specatts);

	nccreate(file_sml_f,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_sml_p,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_sml_d,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_med_f,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_med_p,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_med_d,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_lrg_p,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);
	nccreate(file_lrg_d,"clev","X",X,X_atts,"time",tim,timatts,atts=fracatts);

	nccreate(file_sml_f,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_sml_p,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_sml_d,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_med_f,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_med_p,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_med_d,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_lrg_p,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)
	nccreate(file_lrg_d,"S","X",X,X_atts,"time",tim,timatts,atts=fracatts)

	nccreate(file_sml_f,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_sml_p,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_sml_d,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_med_f,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_med_p,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_med_d,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_lrg_p,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)
	nccreate(file_lrg_d,"DD","X",X,X_atts,"time",tim,timatts,atts=DDatts)

	# #! Initializing netcdf files
	println("Initializing file system (takes about 5 minutes)")
	ncwrite(zeros(NX,1),file_sml_f,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"biomass",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"biomass",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"prod",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"prod",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"prod",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"prod",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"prod",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"prod",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"prod",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"prod",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"con",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"con",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"con",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"con",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"con",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"con",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"con",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"con",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"rec",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"rec",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"rec",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"rec",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"rec",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"rec",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"rec",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"rec",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"nu",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"nu",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"nu",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"nu",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"nu",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"nu",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"nu",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"nu",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"gamma",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"gamma",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"rep",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"rep",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"rep",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"rep",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"rep",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"rep",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"rep",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"rep",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"egg",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"egg",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"egg",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"egg",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"egg",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"egg",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"egg",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"egg",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"die",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"die",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"die",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"die",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"die",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"die",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"die",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"die",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"clev",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"clev",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"clev",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"clev",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"clev",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"clev",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"clev",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"clev",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"S",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"S",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"S",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"S",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"S",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"S",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"S",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"S",[1,1])

	ncwrite(zeros(NX,1),file_sml_f,"DD",[1,1])
	ncwrite(zeros(NX,1),file_sml_p,"DD",[1,1])
	ncwrite(zeros(NX,1),file_sml_d,"DD",[1,1])
	ncwrite(zeros(NX,1),file_med_f,"DD",[1,1])
	ncwrite(zeros(NX,1),file_med_p,"DD",[1,1])
	ncwrite(zeros(NX,1),file_med_d,"DD",[1,1])
	ncwrite(zeros(NX,1),file_lrg_p,"DD",[1,1])
	ncwrite(zeros(NX,1),file_lrg_d,"DD",[1,1])

	###################### Run the Model
	#! Run model with no fishing
	for YR = 1:YEARS # years

		#reset spawning flag
		if (phen == 1)
			Med_f.S = zeros(Float64,NX,DAYS)
			Lrg_d.S = zeros(Float64,NX,DAYS)
			Lrg_p.S = zeros(Float64,NX,DAYS)
		end

		for DAY = 1:DT:DAYS # days

			###! Future time step
			DY  = Int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

			if (YR == YEARS)
				##################### Clean up
				#! Store last year of spinup
				for i = 1:NX
					S_Sml_f[i,DY] = Sml_f.bio[i]
					S_Sml_p[i,DY] = Sml_p.bio[i]
					S_Sml_d[i,DY] = Sml_d.bio[i]
					S_Med_f[i,DY] = Med_f.bio[i]
					S_Med_p[i,DY] = Med_p.bio[i]
					S_Med_d[i,DY] = Med_d.bio[i]
					S_Lrg_p[i,DY] = Lrg_p.bio[i]
					S_Lrg_d[i,DY] = Lrg_d.bio[i]

					S_Sml_f_rec[i,DY] = Sml_f.rec[i]
					S_Sml_p_rec[i,DY] = Sml_p.rec[i]
					S_Sml_d_rec[i,DY] = Sml_d.rec[i]
					S_Med_f_rec[i,DY] = Med_f.rec[i]
					S_Med_p_rec[i,DY] = Med_p.rec[i]
					S_Med_d_rec[i,DY] = Med_d.rec[i]
					S_Lrg_p_rec[i,DY] = Lrg_p.rec[i]
					S_Lrg_d_rec[i,DY] = Lrg_d.rec[i]

					S_Sml_f_con[i,DY] = Sml_f.I[i]
					S_Sml_p_con[i,DY] = Sml_p.I[i]
					S_Sml_d_con[i,DY] = Sml_d.I[i]
					S_Med_f_con[i,DY] = Med_f.I[i]
					S_Med_p_con[i,DY] = Med_p.I[i]
					S_Med_d_con[i,DY] = Med_d.I[i]
					S_Lrg_p_con[i,DY] = Lrg_p.I[i]
					S_Lrg_d_con[i,DY] = Lrg_d.I[i]

					S_Sml_f_nu[i,DY] = Sml_f.nu[i]
					S_Sml_p_nu[i,DY] = Sml_p.nu[i]
					S_Sml_d_nu[i,DY] = Sml_d.nu[i]
					S_Med_f_nu[i,DY] = Med_f.nu[i]
					S_Med_p_nu[i,DY] = Med_p.nu[i]
					S_Med_d_nu[i,DY] = Med_d.nu[i]
					S_Lrg_p_nu[i,DY] = Lrg_p.nu[i]
					S_Lrg_d_nu[i,DY] = Lrg_d.nu[i]

					S_Sml_f_prod[i,DY] = Sml_f.prod[i]
					S_Sml_p_prod[i,DY] = Sml_p.prod[i]
					S_Sml_d_prod[i,DY] = Sml_d.prod[i]
					S_Med_f_prod[i,DY] = Med_f.prod[i]
					S_Med_p_prod[i,DY] = Med_p.prod[i]
					S_Med_d_prod[i,DY] = Med_d.prod[i]
					S_Lrg_p_prod[i,DY] = Lrg_p.prod[i]
					S_Lrg_d_prod[i,DY] = Lrg_d.prod[i]

					S_Sml_f_gamma[i,DY] = Sml_f.gamma[i]
					S_Sml_p_gamma[i,DY] = Sml_p.gamma[i]
					S_Sml_d_gamma[i,DY] = Sml_d.gamma[i]
					S_Med_f_gamma[i,DY] = Med_f.gamma[i]
					S_Med_p_gamma[i,DY] = Med_p.gamma[i]
					S_Med_d_gamma[i,DY] = Med_d.gamma[i]
					S_Lrg_p_gamma[i,DY] = Lrg_p.gamma[i]
					S_Lrg_d_gamma[i,DY] = Lrg_d.gamma[i]

					S_Sml_f_rep[i,DY] = Sml_f.rep[i]
					S_Sml_p_rep[i,DY] = Sml_p.rep[i]
					S_Sml_d_rep[i,DY] = Sml_d.rep[i]
					S_Med_f_rep[i,DY] = Med_f.rep[i]
					S_Med_p_rep[i,DY] = Med_p.rep[i]
					S_Med_d_rep[i,DY] = Med_d.rep[i]
					S_Lrg_p_rep[i,DY] = Lrg_p.rep[i]
					S_Lrg_d_rep[i,DY] = Lrg_d.rep[i]

					S_Sml_f_egg[i,DY] = Sml_f.egg[i]
					S_Sml_p_egg[i,DY] = Sml_p.egg[i]
					S_Sml_d_egg[i,DY] = Sml_d.egg[i]
					S_Med_f_egg[i,DY] = Med_f.egg[i]
					S_Med_p_egg[i,DY] = Med_p.egg[i]
					S_Med_d_egg[i,DY] = Med_d.egg[i]
					S_Lrg_p_egg[i,DY] = Lrg_p.egg[i]
					S_Lrg_d_egg[i,DY] = Lrg_d.egg[i]

					S_Sml_f_die[i,DY] = Sml_f.die[i]
					S_Sml_p_die[i,DY] = Sml_p.die[i]
					S_Sml_d_die[i,DY] = Sml_d.die[i]
					S_Med_f_die[i,DY] = Med_f.die[i]
					S_Med_p_die[i,DY] = Med_p.die[i]
					S_Med_d_die[i,DY] = Med_d.die[i]
					S_Lrg_p_die[i,DY] = Lrg_p.die[i]
					S_Lrg_d_die[i,DY] = Lrg_d.die[i]

					S_Sml_f_clev[i,DY] = Sml_f.clev[i]
					S_Sml_p_clev[i,DY] = Sml_p.clev[i]
					S_Sml_d_clev[i,DY] = Sml_d.clev[i]
					S_Med_f_clev[i,DY] = Med_f.clev[i]
					S_Med_p_clev[i,DY] = Med_p.clev[i]
					S_Med_d_clev[i,DY] = Med_d.clev[i]
					S_Lrg_p_clev[i,DY] = Lrg_p.clev[i]
					S_Lrg_d_clev[i,DY] = Lrg_d.clev[i]

					S_Sml_f_S[i,DY] = Sml_f.S[i]
					S_Sml_p_S[i,DY] = Sml_p.S[i]
					S_Sml_d_S[i,DY] = Sml_d.S[i]
					S_Med_f_S[i,DY] = Med_f.S[i]
					S_Med_p_S[i,DY] = Med_p.S[i]
					S_Med_d_S[i,DY] = Med_d.S[i]
					S_Lrg_p_S[i,DY] = Lrg_p.S[i]
					S_Lrg_d_S[i,DY] = Lrg_d.S[i]

					S_Sml_f_DD[i,DY] = Sml_f.DD[i]
					S_Sml_p_DD[i,DY] = Sml_p.DD[i]
					S_Sml_d_DD[i,DY] = Sml_d.DD[i]
					S_Med_f_DD[i,DY] = Med_f.DD[i]
					S_Med_p_DD[i,DY] = Med_p.DD[i]
					S_Med_d_DD[i,DY] = Med_d.DD[i]
					S_Lrg_p_DD[i,DY] = Lrg_p.DD[i]
					S_Lrg_d_DD[i,DY] = Lrg_d.DD[i]

				end #Grid cells
			end #If last year

		end #Days

	end #Years

	#! Save
	ncwrite(S_Sml_f,file_sml_f,"biomass",[1,1])
	ncwrite(S_Sml_p,file_sml_p,"biomass",[1,1])
	ncwrite(S_Sml_d,file_sml_d,"biomass",[1,1])
	ncwrite(S_Med_f,file_med_f,"biomass",[1,1])
	ncwrite(S_Med_p,file_med_p,"biomass",[1,1])
	ncwrite(S_Med_d,file_med_d,"biomass",[1,1])
	ncwrite(S_Lrg_p,file_lrg_p,"biomass",[1,1])
	ncwrite(S_Lrg_d,file_lrg_d,"biomass",[1,1])

	ncwrite(S_Sml_f_rec,file_sml_f,"rec",[1,1])
	ncwrite(S_Sml_p_rec,file_sml_p,"rec",[1,1])
	ncwrite(S_Sml_d_rec,file_sml_d,"rec",[1,1])
	ncwrite(S_Med_f_rec,file_med_f,"rec",[1,1])
	ncwrite(S_Med_p_rec,file_med_p,"rec",[1,1])
	ncwrite(S_Med_d_rec,file_med_d,"rec",[1,1])
	ncwrite(S_Lrg_p_rec,file_lrg_p,"rec",[1,1])
	ncwrite(S_Lrg_d_rec,file_lrg_d,"rec",[1,1])

	ncwrite(S_Sml_f_con,file_sml_f,"con",[1,1])
	ncwrite(S_Sml_p_con,file_sml_p,"con",[1,1])
	ncwrite(S_Sml_d_con,file_sml_d,"con",[1,1])
	ncwrite(S_Med_f_con,file_med_f,"con",[1,1])
	ncwrite(S_Med_p_con,file_med_p,"con",[1,1])
	ncwrite(S_Med_d_con,file_med_d,"con",[1,1])
	ncwrite(S_Lrg_p_con,file_lrg_p,"con",[1,1])
	ncwrite(S_Lrg_d_con,file_lrg_d,"con",[1,1])

	ncwrite(S_Sml_f_nu,file_sml_f,"nu",[1,1])
	ncwrite(S_Sml_p_nu,file_sml_p,"nu",[1,1])
	ncwrite(S_Sml_d_nu,file_sml_d,"nu",[1,1])
	ncwrite(S_Med_f_nu,file_med_f,"nu",[1,1])
	ncwrite(S_Med_p_nu,file_med_p,"nu",[1,1])
	ncwrite(S_Med_d_nu,file_med_d,"nu",[1,1])
	ncwrite(S_Lrg_p_nu,file_lrg_p,"nu",[1,1])
	ncwrite(S_Lrg_d_nu,file_lrg_d,"nu",[1,1])

	ncwrite(S_Sml_f_prod,file_sml_f,"prod",[1,1])
	ncwrite(S_Sml_p_prod,file_sml_p,"prod",[1,1])
	ncwrite(S_Sml_d_prod,file_sml_d,"prod",[1,1])
	ncwrite(S_Med_f_prod,file_med_f,"prod",[1,1])
	ncwrite(S_Med_p_prod,file_med_p,"prod",[1,1])
	ncwrite(S_Med_d_prod,file_med_d,"prod",[1,1])
	ncwrite(S_Lrg_p_prod,file_lrg_p,"prod",[1,1])
	ncwrite(S_Lrg_d_prod,file_lrg_d,"prod",[1,1])

	ncwrite(S_Sml_f_gamma,file_sml_f,"gamma",[1,1])
	ncwrite(S_Sml_p_gamma,file_sml_p,"gamma",[1,1])
	ncwrite(S_Sml_d_gamma,file_sml_d,"gamma",[1,1])
	ncwrite(S_Med_f_gamma,file_med_f,"gamma",[1,1])
	ncwrite(S_Med_p_gamma,file_med_p,"gamma",[1,1])
	ncwrite(S_Med_d_gamma,file_med_d,"gamma",[1,1])
	ncwrite(S_Lrg_p_gamma,file_lrg_p,"gamma",[1,1])
	ncwrite(S_Lrg_d_gamma,file_lrg_d,"gamma",[1,1])

	ncwrite(S_Sml_f_rep,file_sml_f,"rep",[1,1])
	ncwrite(S_Sml_p_rep,file_sml_p,"rep",[1,1])
	ncwrite(S_Sml_d_rep,file_sml_d,"rep",[1,1])
	ncwrite(S_Med_f_rep,file_med_f,"rep",[1,1])
	ncwrite(S_Med_p_rep,file_med_p,"rep",[1,1])
	ncwrite(S_Med_d_rep,file_med_d,"rep",[1,1])
	ncwrite(S_Lrg_p_rep,file_lrg_p,"rep",[1,1])
	ncwrite(S_Lrg_d_rep,file_lrg_d,"rep",[1,1])

	ncwrite(S_Sml_f_egg,file_sml_f,"egg",[1,1])
	ncwrite(S_Sml_p_egg,file_sml_p,"egg",[1,1])
	ncwrite(S_Sml_d_egg,file_sml_d,"egg",[1,1])
	ncwrite(S_Med_f_egg,file_med_f,"egg",[1,1])
	ncwrite(S_Med_p_egg,file_med_p,"egg",[1,1])
	ncwrite(S_Med_d_egg,file_med_d,"egg",[1,1])
	ncwrite(S_Lrg_p_egg,file_lrg_p,"egg",[1,1])
	ncwrite(S_Lrg_d_egg,file_lrg_d,"egg",[1,1])

	ncwrite(S_Sml_f_die,file_sml_f,"die",[1,1])
	ncwrite(S_Sml_p_die,file_sml_p,"die",[1,1])
	ncwrite(S_Sml_d_die,file_sml_d,"die",[1,1])
	ncwrite(S_Med_f_die,file_med_f,"die",[1,1])
	ncwrite(S_Med_p_die,file_med_p,"die",[1,1])
	ncwrite(S_Med_d_die,file_med_d,"die",[1,1])
	ncwrite(S_Lrg_p_die,file_lrg_p,"die",[1,1])
	ncwrite(S_Lrg_d_die,file_lrg_d,"die",[1,1])

	ncwrite(S_Sml_f_clev,file_sml_f,"clev",[1,1])
	ncwrite(S_Sml_p_clev,file_sml_p,"clev",[1,1])
	ncwrite(S_Sml_d_clev,file_sml_d,"clev",[1,1])
	ncwrite(S_Med_f_clev,file_med_f,"clev",[1,1])
	ncwrite(S_Med_p_clev,file_med_p,"clev",[1,1])
	ncwrite(S_Med_d_clev,file_med_d,"clev",[1,1])
	ncwrite(S_Lrg_p_clev,file_lrg_p,"clev",[1,1])
	ncwrite(S_Lrg_d_clev,file_lrg_d,"clev",[1,1])

	ncwrite(S_Sml_f_S,file_sml_f,"S",[1,1])
	ncwrite(S_Sml_p_S,file_sml_p,"S",[1,1])
	ncwrite(S_Sml_d_S,file_sml_d,"S",[1,1])
	ncwrite(S_Med_f_S,file_med_f,"S",[1,1])
	ncwrite(S_Med_p_S,file_med_p,"S",[1,1])
	ncwrite(S_Med_d_S,file_med_d,"S",[1,1])
	ncwrite(S_Lrg_p_S,file_lrg_p,"S",[1,1])
	ncwrite(S_Lrg_d_S,file_lrg_d,"S",[1,1])

	ncwrite(S_Sml_f_DD,file_sml_f,"DD",[1,1])
	ncwrite(S_Sml_p_DD,file_sml_p,"DD",[1,1])
	ncwrite(S_Sml_d_DD,file_sml_d,"DD",[1,1])
	ncwrite(S_Med_f_DD,file_med_f,"DD",[1,1])
	ncwrite(S_Med_p_DD,file_med_p,"DD",[1,1])
	ncwrite(S_Med_d_DD,file_med_d,"DD",[1,1])
	ncwrite(S_Lrg_p_DD,file_lrg_p,"DD",[1,1])
	ncwrite(S_Lrg_d_DD,file_lrg_d,"DD",[1,1])

	#! Close save
  ncclose(file_sml_f)
  ncclose(file_sml_p)
  ncclose(file_sml_d)
  ncclose(file_med_f)
  ncclose(file_med_p)
  ncclose(file_med_d)
  ncclose(file_lrg_p)
	ncclose(file_lrg_d)

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
	#biomatts = {"longname" => "Biomass","units" => "kg/m^2"}
	#X_atts = {"longname" => "Space", "units" => "grid cell"}
	#S_atts = {"longname" => "Size classes", "units"  => "g"}
	#timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}
	#Use "Dict{Any,Any}(a=>b, ...)" instead.
	biomatts = Dict("longname" => "Biomass",
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
	nccreate(file_pisc,"biomass","X",X,X_atts,"S",S_pi,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_plan,"biomass","X",X,X_atts,"S",S_pl,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_detr,"biomass","X",X,X_atts,"S",S_de,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_bent,"biomass","X",X,X_atts,"S",S_be,S_atts,"time",tim,timatts,atts=biomatts)

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
	biomatts = Dict("longname" => "Biomass",
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
	nccreate(file_pisc,"biomass","X",X,X_atts,"S",S_pi,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_plan,"biomass","X",X,X_atts,"S",S_pl,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_detr,"biomass","X",X,X_atts,"S",S_de,S_atts,"time",tim,timatts,atts=biomatts)
	nccreate(file_bent,"biomass","X",X,X_atts,"S",S_be,S_atts,"time",tim,timatts,atts=biomatts)

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
