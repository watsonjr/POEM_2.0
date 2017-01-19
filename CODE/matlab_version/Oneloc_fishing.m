%%%%!! RUN SPINUP WITH DIFF FISHING RATES FOR ONE LOCATION
function Oneloc_fishing()

	%! setup spinup (loop first year of COBALT)
  COBALT = load('/Volumes/GFDL/POEM_JLD/Data_hindcast_PC_000120.jld'); % 120=1980
	%! Add phenology params from csv file with ID as row
	Tref = readdlm('./Data/grid_phenol_T0raw_NOflip.csv',','); %min temp for each yr at each location
	global TrefP = Tref
	global TrefB = Tref
	global Dthresh = readdlm('./Data/grid_phenol_DTraw_NOflip.csv',',');
	global Sp = readdlm('./Data/Gaussian_spawn_2mo.csv',',');
	YEARS = 100
  global DAYS = 365

	%! choose where to run the model
	global GRD = load('./Data/Data_grid_hindcast_NOTflipped.jld')
	XY = zeros(Int,360,200);
  XY(GRD('ID')) = collect(1:GRD('N'))
	ids = (40319,42639,41782,36334,38309,42744,30051,41284,38003)
	names = ('GB','EBS','OSP','HOT','BATS','NS','EEP','K2','S1')

	simname = 'Dc_TrefO_Hartvig_cmax-metab_MFeqMP_fcrit30_MZ01_NOnmort_BE05_RE0075_MF';

	%! Storage variables
	S_Sml_f = zeros(DAYS,7);
	S_Sml_p = zeros(DAYS,7);
	S_Sml_d = zeros(DAYS,7);
	S_Med_f = zeros(DAYS,7);
	S_Med_p = zeros(DAYS,7);
	S_Med_d = zeros(DAYS,7);
	S_Lrg_p = zeros(DAYS,7);
	S_Lrg_d = zeros(DAYS,7);

	S_Sml_f_catch = zeros(DAYS,7);
	S_Sml_p_catch = zeros(DAYS,7);
	S_Sml_d_catch = zeros(DAYS,7);
	S_Med_f_catch = zeros(DAYS,7);
	S_Med_p_catch = zeros(DAYS,7);
	S_Med_d_catch = zeros(DAYS,7);
	S_Lrg_p_catch = zeros(DAYS,7);
	S_Lrg_d_catch = zeros(DAYS,7);

	%! Make parameters
	harv = 1 %0=no fishing; 1=fishing
	fishing = collect(0.1:0.1:0.7)
	phen=0;

	for L = 6 %1:9
		ID = ids(L)
		loc = names(L)
		const global NX = length(ID)


		%! Storage
		Spinup_Sml_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_f.csv'),'w')
		Spinup_Sml_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_p.csv'),'w')
		Spinup_Sml_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_d.csv'),'w')
		Spinup_Med_f  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_f.csv'),'w')
		Spinup_Med_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_p.csv'),'w')
		Spinup_Med_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_d.csv'),'w')
		Spinup_Lrg_p  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Lrg_p.csv'),'w')
		Spinup_Lrg_d  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Lrg_d.csv'),'w')

		Spinup_Sml_f_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_f_catch.csv'),'w')
		Spinup_Sml_p_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_p_catch.csv'),'w')
		Spinup_Sml_d_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Sml_d_catch.csv'),'w')
		Spinup_Med_f_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_f_catch.csv'),'w')
		Spinup_Med_p_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_p_catch.csv'),'w')
		Spinup_Med_d_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Med_d_catch.csv'),'w')
		Spinup_Lrg_p_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Lrg_p_catch.csv'),'w')
		Spinup_Lrg_d_catch  = open(string('/Volumes/GFDL/CSV/',simname, '/Spinup_',loc,'_Lrg_d_catch.csv'),'w')

		for r = 1:7
			frate = fishing(r)
			make_parameters(harv, frate) % make core parameters/constants
			%! Initialize
			Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
			Med_d.td(1) = 0.0;
			Lrg_d.td(1) = 0.0;
			ENVR = sub_init_env(ID);
			%! Iterate forward in time with NO fishing
			for YR = 1:YEARS % years
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
					%! Store
					S_Sml_f(DY,r) = Sml_f.bio(1)
					S_Sml_p(DY,r) = Sml_p.bio(1)
					S_Sml_d(DY,r) = Sml_d.bio(1)
					S_Med_f(DY,r) = Med_f.bio(1)
					S_Med_p(DY,r) = Med_p.bio(1)
					S_Med_d(DY,r) = Med_d.bio(1)
					S_Lrg_p(DY,r) = Lrg_p.bio(1)
					S_Lrg_d(DY,r) = Lrg_d.bio(1)
					S_Sml_f_catch(DY,r) = Sml_f.caught(1)
					S_Sml_p_catch(DY,r) = Sml_p.caught(1)
					S_Sml_d_catch(DY,r) = Sml_d.caught(1)
					S_Med_f_catch(DY,r) = Med_f.caught(1)
					S_Med_p_catch(DY,r) = Med_p.caught(1)
					S_Med_d_catch(DY,r) = Med_d.caught(1)
					S_Lrg_p_catch(DY,r) = Lrg_p.caught(1)
					S_Lrg_d_catch(DY,r) = Lrg_d.caught(1)
					DY+=1
				end %Days
				%! Save
				if (YR == YEARS)
					writecsv(Spinup_Sml_f,S_Sml_f)
					writecsv(Spinup_Sml_p,S_Sml_p)
					writecsv(Spinup_Sml_d,S_Sml_d)
					writecsv(Spinup_Med_f,S_Med_f)
					writecsv(Spinup_Med_p,S_Med_p)
					writecsv(Spinup_Med_d,S_Med_d)
					writecsv(Spinup_Lrg_p,S_Lrg_p)
					writecsv(Spinup_Lrg_d,S_Lrg_d)
					writecsv(Spinup_Sml_f_catch,S_Sml_f_catch)
					writecsv(Spinup_Sml_p_catch,S_Sml_p_catch)
					writecsv(Spinup_Sml_d_catch,S_Sml_d_catch)
					writecsv(Spinup_Med_f_catch,S_Med_f_catch)
					writecsv(Spinup_Med_p_catch,S_Med_p_catch)
					writecsv(Spinup_Med_d_catch,S_Med_d_catch)
					writecsv(Spinup_Lrg_p_catch,S_Lrg_p_catch)
					writecsv(Spinup_Lrg_d_catch,S_Lrg_d_catch)
				end
			end %Years
		end %fishing rates
		%%% close save
		close(Spinup_Sml_f)
		close(Spinup_Sml_p)
		close(Spinup_Sml_d)
		close(Spinup_Med_f)
		close(Spinup_Med_p)
		close(Spinup_Med_d)
		close(Spinup_Lrg_p)
		close(Spinup_Lrg_d)
		close(Spinup_Sml_f_catch)
		close(Spinup_Sml_p_catch)
		close(Spinup_Sml_d_catch)
		close(Spinup_Med_f_catch)
		close(Spinup_Med_p_catch)
		close(Spinup_Med_d_catch)
		close(Spinup_Lrg_p_catch)
		close(Spinup_Lrg_d_catch)
	end %Locations
end
