
####!! RUN SPINUP FOR ONE LOCATION
function run_testoneloc()

	#! Make parameters
	make_parameters() # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # first year's data 

	#! choose where to run the model
	XY = zeros(360,200); # choose a particulat place or everywhere
	XY[GRD["ID"]] =[1:GRD["N"]]
	#ID = XY[195,102] # Humboldt
	ID = XY[270,156] # Iberian location
	#ID = XY[265,156] # Iberian location off shore
	#ID = XY[260,156] # Iberian location further off shore
	#ID = XY[250,156] # Iberian location # way off shore
	const global NX = length(ID)

	#! Initialize
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);

	#! Storage
	Spinup_PISC = open("./Data/CSV/Spinup_PISC.csv","w")
    Spinup_PLAN = open("./Data/CSV/Spinup_PLAN.csv","w")
    Spinup_DETR = open("./Data/CSV/Spinup_DETR.csv","w")
    Spinup_BENT = open("./Data/CSV/Spinup_BENT.csv","w")

	#! Iterate forward in time with NO fishing
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days 

			###! ticker
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))

			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

			#! Save
			writecsv(Spinup_PISC,PISC.bio[1]')
			writecsv(Spinup_PLAN,PLAN.bio[1]')
			writecsv(Spinup_DETR,DETR.bio[1]')
			writecsv(Spinup_BENT,BENT.bio[1]')

		end
	end
	### close save
    close(Spinup_PISC)
    close(Spinup_PLAN)
    close(Spinup_DETR)
    close(Spinup_BENT)

end



####!! RUN SPINUP FOR ALL LOCATION
function run_spinup_raw()

	#! Make parameters
	make_parameters() # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # first year's data 

	#! choose where to run the model
	const global NX = 48111
	const global ID = [1:NX]

	#! Storage
	S_PISC = zeros(48111,PI_N)
	S_PLAN = zeros(48111,PL_N)
	S_DETR = zeros(48111,DE_N)
	S_BENT = zeros(48111,1)

	#! Initialize
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID);

	#! Run model with no fishing
	const global FISHING = 00000000000 / (365/DT) # 0MT per year
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days 

			###! Future time step
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

		end

	end

	#! Store
	for i = 1:NX
		S_PISC[i,:] = PISC.bio[i]'
		S_PLAN[i,:] = PLAN.bio[i]'
		S_DETR[i,:] = DETR.bio[i]'
		S_BENT[i,:] = BENT.bio[i]'
	end

	#! Save
	npzwrite("./Data/NPZ/Spinup_raw_PISC.npy",S_PISC)
	npzwrite("./Data/NPZ/Spinup_raw_PLAN.npy",S_PLAN)
	npzwrite("./Data/NPZ/Spinup_raw_DETR.npy",S_DETR)
	npzwrite("./Data/NPZ/Spinup_raw_BENT.npy",S_BENT)

end

####!! RUN SPINUP FOR ALL LOCATION WITH FISHING
function run_spinup_fishing()

	#! Make parameters
	make_parameters() # make core parameters/constants

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # first year's data 

	#! choose where to run the model
	const global NX = 48111
	const global ID = [1:NX]

	#! Storage
	S_PISC = zeros(48111,PI_N)
	S_PLAN = zeros(48111,PL_N)
	S_DETR = zeros(48111,DE_N)
	S_BENT = zeros(48111,1)

	#! Initialize
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
	PISC = npzread("./Data/NPZ/Spinup_raw_PISC.npy")
	PLAN = npzread("./Data/NPZ/Spinup_raw_PLAN.npy")
	DETR = npzread("./Data/NPZ/Spinup_raw_DETR.npy")
	BENT = npzread("./Data/NPZ/Spinup_raw_BENT.npy")
	ENVR = sub_init_env(ID);

	#! Run model with no fishing
	const global FISHING = 80000000000 / (365/DT) # 80MT per year
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days 

			###! Future time step
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

		end

	end

	#! Store
	for i = 1:NX
		S_PISC[i,:] = PISC.bio[i]'
		S_PLAN[i,:] = PLAN.bio[i]'
		S_DETR[i,:] = DETR.bio[i]'
		S_BENT[i,:] = BENT.bio[i]'
	end

	#! Save
	npzwrite("./Data/NPZ/Spinup_fishing_PISC.npy",S_PISC)
	npzwrite("./Data/NPZ/Spinup_fishing_PLAN.npy",S_PLAN)
	npzwrite("./Data/NPZ/Spinup_fishing_DETR.npy",S_DETR)
	npzwrite("./Data/NPZ/Spinup_fishing_BENT.npy",S_BENT)

end

####!! RUN FORECAST FOR ALL LOCATIONS
function run_forecast()

	#! Load COBALT and grid data
	COBALT = load("./Data/Data_000001.jld"); # first year's data 
	GRD = load("./Data/Data_grid.jld"); # spatial information

	#! Make parameters
	make_parameters() # make core parameters/constants
	const global DT = 1.; # time step
	const global YEARS = 1; # integration period in years
	const global DAYS = 365; # number of days 
	const global ID = [1:size(GRD["ID"])[1]]
	const global NX = length(ID)
	const global MNTH = [31,28,31,30,31,30,31,31,30,31,30,31]

	#! Initialize
	PISC,PLAN,DETR,BENT = sub_init_fish(ID);
	ENVR = sub_init_env(ID); 

	#! Storage arrays (daily and monthly)
	DAY_PISC = Array(Float64,NX,PI_N,DAYS);
	DAY_PLAN = Array(Float64,NX,PL_N,DAYS);
	DAY_DETR = Array(Float64,NX,DE_N,DAYS);
	DAY_BENT = Array(Float64,NX,1,DAYS);
	MNT_PISC = Array(Float64,NX,PI_N,12);
	MNT_PLAN = Array(Float64,NX,PL_N,12);
	MNT_DETR = Array(Float64,NX,DE_N,12);
	MNT_BENT = Array(Float64,NX,1,12);

	#! Init netcdf file for storage
	var_atts = @Compat.Dict("longname" => "Fish biomass",
			"units"    => "kg/m^2")
	space_atts = @Compat.Dict("longname" => "Space",
			"units"    => "grid cell")
	time_atts = @Compat.Dict("longname" => "Time",
			"units"    => "months since 01-01-2006 00:00:00")
	pisc_atts = @Compat.Dict("longname" => "Sizeclass",
			"units"    => "months since 01-01-2006 00:00:00")

	latdim = NcDim("lat",lat,latatts)
	londim = NcDim("lon",lon,lonatts)
	timdim = NcDim("time",tim,timatts)


	#! Iterate forward in time
	for YR = 1:YEARS # years

		for DAY = 1:DT:DAYS # days 

			###! ticker
			DY  = int(ceil(DAY))
			println(YR," , ", mod(DY,365))

			###! Future time step
			sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT);

		end
	
		#! Calculate monthly means
		a = [1,(cumsum(MNTH)+1)[1:end-1]] # start of the month
		b = cumsum(MNTH) # end of the month
		for i = 1:12
			MNT_PISC[i] = mean(DAY_PISC[:,:,a[i]:b[i]])
			MNT_PLAN[i] = mean(DAY_PLAN[:,:,a[i]:b[i]])
			MNT_DETR[i] = mean(DAY_DETR[:,:,a[i]:b[i]])
			MNT_BENT[i] = mean(DAY_BENT[:,:,a[i]:b[i]])
		end

		#! Save monthly means
		npzwrite("./Data/NPZ/Spinup_PISC.npz",["x"=>1,"bio"=>MNT_PISC])
		npzwrite("./Data/NPZ/Spinup_PLAN.npz",["x"=>1,"bio"=>MNT_PLAN])
		npzwrite("./Data/NPZ/Spinup_DETR.npz",["x"=>1,"bio"=>MNT_DETR])
		npzwrite("./Data/NPZ/Spinup_BENT.npz",["x"=>1,"bio"=>MNT_BENT])
	end

end


