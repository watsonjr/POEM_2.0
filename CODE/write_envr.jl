using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT #,MATLAB
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")


make_parameters(0,0) # make core parameters/constants

COBALT = load("/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_2005.jld"); # 120=1980

#! Add phenology params from csv file with ID as row
Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
global TrefP = Tref
global TrefB = Tref
global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
global DAYS = 365
const global MNTH = collect([31,28,31,30,31,30,31,31,30,31,30,31]) # days in month

#! choose where to run the model
global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld")
const global NX = 48111
const global ID = collect(1:NX);

simname = "ESM2M_2005_Julia_interp";
Det_2005_mo  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_det_mo.csv"),"w")
Zm_2005_mo  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_Zm_mo.csv"),"w")
Zl_2005_mo  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_Zl_mo.csv"),"w")
dZm_2005_mo  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_dZm_mo.csv"),"w")
dZl_2005_mo  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_dZl_mo.csv"),"w")
Det_2005_day  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_det_day.csv"),"w")
Zm_2005_day  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_Zm_day.csv"),"w")
Zl_2005_day  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_Zl_day.csv"),"w")
dZm_2005_day  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_dZm_day.csv"),"w")
dZl_2005_day  = open(string("/Volumes/GFDL/POEM_JLD/esm2m_hist/",simname, "_dZl_day.csv"),"w")

writecsv(Det_2005_day,COBALT["det"])
writecsv(Zm_2005_day,COBALT["Zm"])
writecsv(Zl_2005_day,COBALT["Zl"])
writecsv(dZm_2005_day,COBALT["dZm"])
writecsv(dZl_2005_day,COBALT["dZl"])

#! Run model with no fishing
S_Det = zeros(NX,DAYS);
S_Zm = zeros(NX,DAYS);
S_Zl = zeros(NX,DAYS);
S_dZm = zeros(NX,DAYS);
S_dZl = zeros(NX,DAYS);

ENVR = sub_init_env(ID);

MNT = 0
YR=1;
for DAY = 1:DT:DAYS # days

  ###! Future time step
  DY  = Int(ceil(DAY))
  println(YR," , ", mod(DY,365))
  get_COBALT!(COBALT,ID,DY,ENVR)
	sub_neg!(ENVR.det);
	sub_neg!(ENVR.Zm);
	sub_neg!(ENVR.Zl);
	sub_neg!(ENVR.dZm);
	sub_neg!(ENVR.dZl);

  #! Store
  for i = 1:NX
    S_Det[i,DY] = ENVR.det[i]
    S_Zm[i,DY] = ENVR.Zm[i]
    S_Zl[i,DY] = ENVR.Zl[i]
    S_dZm[i,DY] = ENVR.dZm[i]
    S_dZl[i,DY] = ENVR.dZl[i]
  end #Grid cells

end #Days

#! Calculate monthly means and save
a = [1;(cumsum(MNTH)+1)[1:end-1]] # start of the month
b = cumsum(MNTH) # end of the month
for i = 1:12
  MNT += 1 # Update monthly ticker
  writecsv(Det_2005_mo,[mean(S_Det[:,a[i]:b[i]],2)])
  writecsv(Zm_2005_mo,[mean(S_Zm[:,a[i]:b[i]],2)])
  writecsv(Zl_2005_mo,[mean(S_Zl[:,a[i]:b[i]],2)])
  writecsv(dZm_2005_mo,[mean(S_dZm[:,a[i]:b[i]],2)])
  writecsv(dZl_2005_mo,[mean(S_dZl[:,a[i]:b[i]],2)])

end #Monthly mean

### close save
close(Det_2005_day)
close(Zm_2005_day)
close(Zl_2005_day)
close(dZm_2005_day)
close(dZl_2005_day)
close(Det_2005_mo)
close(Zm_2005_mo)
close(Zl_2005_mo)
close(dZm_2005_mo)
close(dZl_2005_mo)
