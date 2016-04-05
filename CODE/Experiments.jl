
####!! RUN SPINUP FOR ONE LOCATION
function Testoneloc()

	#! Make parameters
	make_parameters(0) # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # first year's data 
	YEARS = 10
	DAYS = 365

	#! choose where to run the model
	GRD = load("./Data/Data_grid.jld")
	XY = zeros(Int,360,200); # choose a particulat place or everywhere
	XY[GRD["ID"]] = collect(1:GRD["N"])
	#ID = XY[195,102] # Humboldt
	#ID = XY[270,156] # Iberian location
	#ID = XY[265,156] # Iberian location off shore
	#ID = XY[260,156] # Iberian location further off shore
	#ID = XY[250,156] # Iberian location # way off shore
	ID = 1 # 70 m depth
	const global NX = length(ID)

	#! Initialize
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);

	#! Storage
	Spinup_Sml_f = open("./Data/CSV/Spinup_Sml_f.csv","w")
    Spinup_Sml_p = open("./Data/CSV/Spinup_Sml_p.csv","w")
    Spinup_Sml_d = open("./Data/CSV/Spinup_Sml_d.csv","w")
    Spinup_Med_f = open("./Data/CSV/Spinup_Med_f.csv","w")
    Spinup_Med_p = open("./Data/CSV/Spinup_Med_p.csv","w")
    Spinup_Med_d = open("./Data/CSV/Spinup_Med_d.csv","w")
    Spinup_Lrg_p = open("./Data/CSV/Spinup_Lrg_p.csv","w")

	#! Iterate forward in time with NO fishing
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days 

			###! ticker
			DY  = Int(ceil(DAY))
			println(YR," , ", mod(DY,365))

			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,BENT);
			DY+=1

			#! Save
			writecsv(Spinup_Sml_f,Sml_f.bio)
			writecsv(Spinup_Sml_p,Sml_p.bio)
			writecsv(Spinup_Sml_d,Sml_d.bio)
			writecsv(Spinup_Med_f,Med_f.bio)
			writecsv(Spinup_Med_p,Med_p.bio)
			writecsv(Spinup_Med_d,Med_d.bio)
			writecsv(Spinup_Lrg_p,Lrg_p.bio)

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
	Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);

	############### Setup NetCDF save
	#! Init netcdf file for storage
	varatts = {"longname" => "Biomass","units" => "kg/m^2"}
	X_atts = {"longname" => "Space", "units" => "grid cell"}
	timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}

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
			sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,BENT);

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
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);


	############### Setup NetCDF save
	#! Init netcdf file for storage
	varatts = {"longname" => "Biomass","units" => "kg/m^2"}
	X_atts = {"longname" => "Space", "units" => "grid cell"}
	S_atts = {"longname" => "Size classes", "units"  => "g"}
	timatts = {"longname" => "Time", "units" => "hours since 01-01-2000 00:00:00"}

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
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

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
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
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
	varatts = {"longname" => "Biomass",
           "units"    => "kg/m^2"}
	X_atts = {"longname" => "Space",
			"units"    => "grid cell"}
	S_atts = {"longname" => "Size classes",
			"units"    => "g"}
	timatts = {"longname" => "Time",
			"units"    => "hours since 01-01-2000 00:00:00"}

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
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

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


