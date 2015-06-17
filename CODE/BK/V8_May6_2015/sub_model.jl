
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ
include("Parameters.jl")
include("sub_init.jl")
include("sub_functions.jl")
include("sub_routines.jl")

####!! RUN SPINUP
function run_testoneloc()

	#! Make parameters
	make_parameters() # make core parameters/constants
	const global DT = 1.0; # time step
	const global YEARS = 10; # integration period in years
	const global DAYS = 365; # number of days 

	#! setup spinup (loop first year of COBALT)
	COBALT = load("./Data/Data_000001.jld"); # first year's data 

	#! choose where to run the model
	GRD = load("./Data/Data_grid.jld"); # spatial information
	XY = zeros(360,200); # choose a particulat place or everywhere
	XY[GRD["ID"]] =[1:GRD["N"]]
	ID = XY[272,156] # Iberian location
	NX = length(ID)

	#! Initialize
	PISC,PLAN,DETR,W = sub_init_fish(ID);
	ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,ENV_dZm,ENV_dZl,ENV_det = sub_init_env(ID); 

	#! Storage
	Spinup_PISC = open("./Data/CSV/Spinup_PISC.csv","w")
    Spinup_PLAN = open("./Data/CSV/Spinup_PLAN.csv","w")
    Spinup_DETR = open("./Data/CSV/Spinup_DETR.csv","w")
    Spinup_W    = open("./Data/CSV/Spinup_W.csv","w")

	#! Iterate forward in time
	for YR = 1:YEARS # years

		for DAY = 1:DAYS/DT # days 

			###! ticker
			DY  = mod(int(ceil(DAY)),365)+1
			println(YR," , ", mod(DY,365))

			###! COBALT information
			get_COBALT!(COBALT,ID,DY,ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,
						ENV_dZm,ENV_dZl,ENV_det);

			###! DEMOGRAPHIC CALCULATIONS
			#! temperature multiplier
			PISC.tmet = sub_tmet(PISC.tmet,ENV_Tp);
			PLAN.tmet = sub_tmet(PLAN.tmet,ENV_Tp);
			DETR.tmet = sub_tmet(DETR.tmet,ENV_Tb);

			#! metabolism
			map(sub_metabolism_pi!,PISC.met,PISC.tmet);
			map(sub_metabolism_pl!,PLAN.met,PLAN.tmet);
			map(sub_metabolism_de!,DETR.met,DETR.tmet);

			#! fraction of time piscivore spends in pelagic
			PISC.tdif = map(sub_tdif,PLAN.bio,DETR.bio);

			#! handling times
			map(sub_tau_pi!,PISC.tau,PISC.met);
			map(sub_tau_pl!,PLAN.tau,PLAN.met);
			map(sub_tau_de!,DETR.tau,DETR.met);

			#! encounter rates
			map(sub_enc_pipi!,PISC.enc_pi,PISC.bio,PISC.tdif);
			map(sub_enc_pipl!,PISC.enc_pl,PLAN.bio,PISC.tdif);
			map(sub_enc_pide!,PISC.enc_de,DETR.bio,PISC.tdif);
			map(sub_enc_piz!,PISC.enc_z,ENV_Zm,ENV_Zl,PISC.tdif);
			map(sub_enc_plz!,PLAN.enc_z,ENV_Zm,ENV_Zl);
			map(sub_enc_dede!,DETR.enc_de,DETR.bio);
			map(sub_enc_dew!,DETR.ENC_w,W.bio);
			
			#! total biomass encountered of each group
			PISC.ENC_pi = map((x) -> sum(x,1),PISC.enc_pi);
			PISC.ENC_pl = map((x) -> sum(x,1),PISC.enc_pl);
			PISC.ENC_de = map((x) -> sum(x,1),PISC.enc_de);
			PISC.ENC_z  = map((x) -> sum(x,1),PISC.enc_z);
			PLAN.ENC    = map((x) -> sum(x,1),PLAN.enc_z);
			DETR.ENC_de = map((x) -> sum(x,1),DETR.enc_de);

			#! total biomass encountered
			PISC.ENC = PISC.ENC_pi + PISC.ENC_pl + PISC.ENC_de + PISC.ENC_z;
			DETR.ENC = DETR.ENC_de + DETR.ENC_w;

			#! reset consumption and mortality
			PISC.I,PLAN.I,DETR.I = sub_reset_con_fish(PISC.I,PLAN.I,DETR.I);
			PISC.I_z,PLAN.I_z    = sub_reset_con_zoo(PISC.I_z,PLAN.I_z);
			PISC.d,PLAN.d,DETR.d,W.d = sub_reset_mort(PISC.d,PLAN.d,DETR.d,W.d)

			#! total biomass consumed and lost to predation SPEED UP
			map(sub_consume_pipi!,PISC.I,PISC.d,PISC.bio,PISC.enc_pi,PISC.ENC,PISC.tau);
			map(sub_consume_pipl!,PISC.I,PLAN.d,PISC.bio,PISC.enc_pl,PISC.ENC,PISC.tau);
			map(sub_consume_pide!,PISC.I,DETR.d,PISC.bio,PISC.enc_de,PISC.ENC,PISC.tau);
			map(sub_consume_dede!,DETR.I,DETR.d,DETR.bio,DETR.enc_de,DETR.ENC,DETR.tau);
			map(sub_consume_dew!,DETR.I,W.d,DETR.bio,DETR.ENC_w,DETR.ENC,DETR.tau);

			#! zooplankton consumption
			map(sub_consume_piz!,PISC.I_z,PISC.bio,PISC.enc_z,PISC.ENC,PISC.tau);
			map(sub_consume_plz!,PLAN.I_z,PLAN.bio,PLAN.enc_z,PLAN.ENC,PLAN.tau);
			
			#! OFFLINE Coupling
			map(sub_offline!,PISC.I_z,PLAN.I_z,ENV_dZm,ENV_dZl)

			#! total energy available for growth nu
			map(sub_nu_pi!,PISC.nu,PISC.bio,PISC.I,PISC.met);
			map(sub_nu_pl!,PLAN.nu,PLAN.bio,PLAN.I,PLAN.met);
			map(sub_nu_de!,DETR.nu,DETR.bio,DETR.I,DETR.met);

			#! somatic growth
			map(sub_gamma_pi!,PISC.gamma,PISC.nu,PISC.bio,PISC.d);
			map(sub_gamma_pl!,PLAN.gamma,PLAN.nu,PLAN.bio,PLAN.d);
			map(sub_gamma_de!,DETR.gamma,DETR.nu,DETR.bio,DETR.d);

			#! egg production
			map(sub_rep_pi!,PISC.REP,PISC.nu,PISC.bio);
			map(sub_rep_pl!,PLAN.REP,PLAN.nu,PLAN.bio);
			map(sub_rep_de!,DETR.REP,DETR.nu,DETR.bio);

			#! total biomass somatic growth
			map(sub_grw_pi!,PISC.GRW,PISC.nu,PISC.bio);
			map(sub_grw_pl!,PLAN.GRW,PLAN.nu,PLAN.bio);
			map(sub_grw_de!,DETR.GRW,DETR.nu,DETR.bio);

			#! total biomass maturing
			map(sub_mat_pi!,PISC.MAT,PISC.gamma,PISC.bio);
			map(sub_mat_pl!,PLAN.MAT,PLAN.gamma,PLAN.bio);
			map(sub_mat_de!,DETR.MAT,DETR.gamma,DETR.bio);

			#! total biomass lost to natural mortality
			map(sub_mrt_pi!,PISC.MRT,PISC.bio);
			map(sub_mrt_pl!,PLAN.MRT,PLAN.bio);
			map(sub_mrt_de!,DETR.MRT,DETR.bio);

			#! Mass balance SPEED UP
			map(sub_update_pi!,PISC.bio,PISC.REP,PISC.GRW,PISC.MAT,PISC.d,PISC.MRT);
			map(sub_update_pl!,PLAN.bio,PLAN.REP,PLAN.GRW,PLAN.MAT,PLAN.d,PLAN.MRT);
			map(sub_update_de!,DETR.bio,DETR.REP,DETR.GRW,DETR.MAT,DETR.d,DETR.MRT);
			map(sub_update_w!,W.bio,W.d,ENV_det);

			#! Forward Euler checks
			map(sub_check_pi!,PISC.bio);
			map(sub_check_pl!,PLAN.bio);
			map(sub_check_de!,DETR.bio);
			map(sub_check_w!,W.bio);

			#! Save
			writecsv(Spinup_PISC,PISC.bio[1]')
			writecsv(Spinup_PLAN,PLAN.bio[1]')
			writecsv(Spinup_DETR,DETR.bio[1]')
			writecsv(Spinup_W,W.bio[1]')

		end
	end
	### close save
    close(Spinup_PISC)
    close(Spinup_PLAN)
    close(Spinup_DETR)
    close(Spinup_W)

end

##### RUN MODEL
run_testoneloc()


