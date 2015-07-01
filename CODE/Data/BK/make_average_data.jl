

###### OFFLINE DATA REQUIRED BY POEM 2.0
###### DAILY AVERAGES FOR SPINUP
#! Zm: medium zooplankton biomass (g m-2)
#! Zl: large zooplankton biomass (g m-2)
#! dZm: medium zooplankton mortality rate (g m-2 day-1)
#! dZl: large zooplankton mortality rate (g m-2 day-1)
#! dDet: detrital flux to the benthos (g m-2 day-1)
#! Tp: pelagic temperature averaged over the top 200m (deg C) 
#! Tb: bottom temperature (deg C)
#! Time and Lon and Lat
using HDF5, JLD

###### LOAD COBALT DATA
#! Prepare daily average files
Tp_mu  = zeros(48111,365)
Tb_mu  = zeros(48111,365)
Zm_mu  = zeros(48111,365)
Zl_mu  = zeros(48111,365)
dZm_mu = zeros(48111,365)
dZl_mu = zeros(48111,365)
det_mu = zeros(48111,365)

##! Load in daily data
for i = 1:10; # number of years to average over
	##! ticker
	id = string(1000000+i);

	##! Load
	Zm  = load("./JLD/Data_"id[2:end]".jld","D_Zm");
	Zl  = load("./JLD/Data_"id[2:end]".jld","D_Zl");
	dZm = load("./JLD/Data_"id[2:end]".jld","D_dZm");
	dZl = load("./JLD/Data_"id[2:end]".jld","D_dZl");
	Tp  = load("./JLD/Data_"id[2:end]".jld","D_Tp");
	Tb  = load("./JLD/Data_"id[2:end]".jld","D_Tb");
	Det = load("./JLD/Data_"id[2:end]".jld","D_det");
	
	##! Add
	Tp_mu  = Tp_mu + Tp;
	Tb_mu  = Tb_mu + Tb;
	Zm_mu  = Zm_mu + Zm;
	Zl_mu  = Zl_mu + Zl;
	dZm_mu = dZm_mu + dZm;
	dZl_mu = dZl_mu + dZl;
	det_mu = det_mu + Det;

	##! ticker
	println(i)

end

##! Calculate average
Tp_mu  = Tp_mu / 94;
Tb_mu  = Tb_mu / 94;
Zm_mu  = Zm_mu / 94;
Zl_mu  = Zl_mu / 94;
dZm_mu = dZm_mu / 94;
dZl_mu = dZl_mu / 94;
det_mu = det_mu / 94;

###! Pull out time series for one place(Iberian sea)
GRD = load("./Data_grid.jld")
XY = zeros(360,200);
XY[GRD["ID"]] =[1:GRD["N"]]
id = XY[270,156] # Iberian location
id = XY[265,156] # further off shore ...
id = XY[260,156]
id = XY[255,156]

##! get data
#Tp_x = squeeze(Tp_mu[id,:],1);
#Tb_x = squeeze(Tb_mu[id,:],1);
#Zm_x = squeeze(Zm_mu[id,:],1);
#Zl_x = squeeze(Zl_mu[id,:],1);
#dZm_x = squeeze(dZm_mu[id,:],1);
#dZl_x = squeeze(dZl_mu[id,:],1);

##! Save
save("./JLD/Data_daily_averages.jld","Tp",Tp_mu,"Tb",Tb_mu,
		"Zm",Zm_mu,"Zl",Zl_mu,"dZm",dZm_mu,"dZl",dZl_mu,"det",det_mu);
#save("./JLD/Data_loc_average.jld","Tp",Tp_x,"Tb",Tb_x,
#        "Zm",Zm_x,"Zl",Zl_x,"dZm",dZm_x,"dZl",dZl_x);

##! Look
#using PyPlot
#GRD = load("./JLD/Data_grid.jld")
#A = zeros(360,200);
#A[GRD["GRD_ID"]] = Tp_mu[:,1]-273;
#A[272,156] = 0.; # Iberian location
#pcolormesh(A)
#plt.clim([0 30]);


