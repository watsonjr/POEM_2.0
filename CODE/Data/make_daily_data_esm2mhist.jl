

###### OFFLINE DATA REQUIRED BY POEM 2.0
#! Zm: medium zooplankton biomass (g m-2)
#! Zl: large zooplankton biomass (g m-2)
#! dZm: medium zooplankton mortality rate (g m-2 day-1)
#! dZl: large zooplankton mortality rate (g m-2 day-1)
#! dDet: detrital flux to the benthos (g m-2 day-1)
#! Tp: pelagic temperature averaged over the top 200m (deg C)
#! Tb: bottom temperature (deg C)
#! Time and Lon and Lat

#! 1) data to be interpolated to daily resolution,
#! and saved in monthly chunks
#! 2) daily average values will be saved too as a separate data file


###### LOAD COBALT DATA
#! ncinfo to look into netCDF files, ncread to load
#! e.g. ncinfo("./GCM/ocean.200601-210012.temp_100_avg.nc");
using NetCDF, HDF5, JLD, Grid, MAT

#! Time
TIME = ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc",
    "average_T1"); # time

#! Physical Scalers: (Pelagic Temp, Bottom Temp (potential temp)), deg C
Tb = ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_btm.186101-200512.btm_temp.nc","btm_temp");
Tp = ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean.186101-200512.temp100.nc","TEMP100");

#! Zooplankton abundances: medium and large (mol C m-2)
Zm=ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_biomass_100.186101-200512.nmdz_100.nc",
	"nmdz_100");
Zl=ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_biomass_100.186101-200512.nlgz_100.nc",
	"nlgz_100");

#! Zooplankton mortality rates: medium and large size: (mol C m-2 s-1)
dZm=ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_miscflux_100.186101-200512.jhploss_nmdz_100.nc",
	"jhploss_nmdz_100");
dZl=ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_miscflux_100.186101-200512.jhploss_nlgz_100.nc",
	"jhploss_nlgz_100");

#! Detrital flux at the sea floor (mol C m-2 s-1)
dDet=ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean_cobalt_btm.186101-200512.fndet_btm.nc","fndet_btm");

#! Ocean currents
#U = ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean.186101-200512.u_100_avg.nc","U_100");
#V = ncread("/Volumes/GFDL/GCM_DATA/ESM2M_hist/ocean.186101-200512.v_100_avg.nc","V_100");

###### INTERPOLATE DATA TO SIZE-BASED MODEL TIME SCALES
#! Save in annual chunks (365 days)
#! load grid data (for pressure to calc bottom temp)
Pr = load("/Volumes/GFDL/POEM_JLD/Data_grid_hindcast.jld","Pr")
R_Cp = 0.11 # R/Cp: gas constant over specific heat capacity
Po = 1000 # millibars (air pressure at sea surface)

#! index of water cells
WID = find(Zm[:,:,1] .!= -1.0e10); # spatial index of water cells
NID = length(WID); # number of water cells

#! months in a year
lstd = Int(TIME[length(TIME)])+31 - Int(TIME[1])
id1 = collect(TIME[1]:365:(lstd-1))
id2 = collect(TIME[1]+365:365:(lstd))
ID  = [id1 id2];
nyr = Int(lstd/365)
#! pull out annual information
#! transform to size-based model units (g, day-1, m-2)
#! x 12.01  mol C --> grams C
#! / 0.32 grams C --> dry weight.
#! *60 *60 *24 --> per day (if flux)
for i = 1:nyr
	id = map(Float64,ID[i,:])
	I = find(id[1].<=TIME.<=id[2])

	#! pull raw data out
	time = TIME[I];
	tp  = Tp[:,:,I];
	tb  = Tb[:,:,I];
	zm  = Zm[:,:,I];
	zl  = Zl[:,:,I];
	dzm = dZm[:,:,I];
	dzl = dZl[:,:,I];
	det = dDet[:,:,I];
  # u = U[:,:,I];
  # v = V[:,:,I];

	#! setup POEM data files
	D_Tp  = zeros(NID,365);
	D_Tb  = zeros(NID,365);
	D_Zm  = zeros(NID,365);
 	D_Zl  = zeros(NID,365);
	D_dZm = zeros(NID,365);
 	D_dZl = zeros(NID,365);
 	D_det = zeros(NID,365);
  # D_u = zeros(NID,365);
  # D_v = zeros(NID,365);

  #! NaN velocities
  # u[find(u.==minimum(u))] = 0.0
  # v[find(v.==minimum(v))] = 0.0

	#! interpolate to daily resolution
	for j = 1:NID
		#! indexes
		m,n = ind2sub((360,200),WID[j]); # spatial index of water cell

    #! v currents from m/s to m/d
		# Y = zeros(size(time))
		# Y[:] = v[m,n,:];
		# yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		# Xi = collect(time[1]:1:time[end]);
		# Yi = yi[Xi[1:end-1]];
		# D_v[j,1:length(Yi)] = Yi * 60 *60 * 24; # m d-1
    #
		# #! u currents from m/s to m/d
		# Y = zeros(size(time))
		# Y[:] = u[m,n,:];
		# yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		# Xi = collect(time[1]:1:time[end]);
		# Yi = yi[Xi[1:end-1]];
		# D_u[j,1:length(Yi)] = Yi * 60 *60 * 24; # m d-1

		#! pelagic temperature (from Kelvin to Celcius)
		Y = zeros(size(time))
		Y[:] = tp[m,n,:] - 273;
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
		Yi = yi[Xi[1:end-1]];
		D_Tp[j,1:length(Yi)] = Yi;

		#! bottom temperature (in Celcius)
		Y = zeros(size(time))
		Y[:] = tb[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
		Yi = yi[Xi[1:end-1]];
		D_Tb[j,1:length(Yi)] = Yi; ## FIX THIS LATER FOR POT TEMP
		#D_Tb[j,:] = ((Yi+273) / ((Po/Pr[j])^R_Cp) ) - 273

		#! medium zoo: from mol N m-2 to g(WW) m-2
    # 106/16 mol C in 1 mol N
    # 12.01 g C in 1 mol C
    # 1 g dry W in 9 g wet W (Pauly & Christiansen)
		Y = zeros(size(time))
		Y[:] = zm[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
    Yi = yi[Xi[1:end-1]];
		D_Zm[j,1:length(Yi)] = Yi * (106.0/16.0) * 12.01 * 9.0;

		#! large zoo: from mol N m-2 to g(WW) m-2
    # 106/16 mol C in 1 mol N
    # 12.01 g C in 1 mol C
    # 1 g dry W in 9 g wet W (Pauly & Christiansen)
		Y = zeros(size(time))
		Y[:] = zl[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
    Yi = yi[Xi[1:end-1]];
		D_Zl[j,1:length(Yi)] = Yi * (106.0/16.0) * 12.01 * 9.0;

		#! medium zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
    # 106/16 mol C in 1 mol N
    # 12.01 g C in 1 mol C
    # 1 g dry W in 9 g wet W (Pauly & Christiansen)
		Y = zeros(size(time))
		Y[:] = dzm[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
		Yi = yi[Xi[1:end-1]];
		D_dZm[j,1:length(Yi)] = Yi * (106.0/16.0) * 12.01 * 9.0 * 60 * 60 *24 ;

		#! large zoo mortality: from mol N m-2 s-1 to g(WW) m-2 d-1
    # 106/16 mol C in 1 mol N
    # 12.01 g C in 1 mol C
    # 1 g dry W in 9 g wet W (Pauly & Christiansen)
		Y = zeros(size(time))
		Y[:] = dzl[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
    Yi = yi[Xi[1:end-1]];
		D_dZl[j,1:length(Yi)] = Yi * (106.0/16.0) * 12.01 * 9.0 *60 *60 *24;

		#! detrital flux to benthos: from mol C m-2 s-1 to g(WW) m-2 d-1
    # 106/16 mol C in 1 mol N
    # 12.01 g C in 1 mol C
    # 1 g dry W in 9 g wet W (Pauly & Christiansen)
		Y = zeros(size(time))
		Y[:] = det[m,n,:];
		yi = InterpIrregular(time, Y, BCnil, InterpLinear);
		Xi = collect(time[1]:1:time[end]);
    Yi = yi[Xi[1:end-1]];
		D_det[j,1:length(Yi)] = Yi * (106.0/16.0) * 12.01 * 9.0 *60 *60 *24;

	end

	#! convert to single precision to save space
	D_Tp  = map(Float64,D_Tp);
	D_Tb  = map(Float64,D_Tb);
	D_Zm  = map(Float64,D_Zm);
	D_Zl  = map(Float64,D_Zl);
	D_dZm = map(Float64,D_dZm);
	D_dZl = map(Float64,D_dZl);
	D_det = map(Float64,D_det);

	#! save
	println(i)
	#ti = string(1000000+i);
  ti = string(1860+i);
  di = "/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_";
	# save(string(di,ti[2:end],".jld"), "Zm",D_Zm,"Zl",D_Zl,"dZm",D_dZm,"dZl",D_dZl,
	# 								  "Tp",D_Tp,"Tb",D_Tb,"det",D_det,"U",D_u,"V",D_v);
  save(string(di,ti[1:end],".jld"), "Zm",D_Zm,"Zl",D_Zl,"dZm",D_dZm,"dZl",D_dZl,
                  									  "Tp",D_Tp,"Tb",D_Tb,"det",D_det);
end
