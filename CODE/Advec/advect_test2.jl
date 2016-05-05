# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID")
GRD = load("./Data/Data_grid_cp_2D.jld")
COBALT = load("./Data/JLD/Data_hindcast_000120.jld"); # 1980

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
bio[ID] = 1.0e3*ones(Float64,size(ID));
U[ID] = COBALT["U"][:,1];
V[ID] = COBALT["V"][:,1];

ni, nj = size(U);

biot = sub_advection(GRD,bio,U,V,ni,nj)

### ADVECTION ###
function sub_advection(GRD,Bio_in,U,V,ni,nj)
	# time steps within one day
	ntime = 24.0
	nt = Int(ntime)
	dtime = 1.0/ntime
	# biol concentration
	Tfield = zeros(Float64,360,200,nt);
	Tfield[:,:,1] = Bio_in;
	Ttendency = zeros(Float64,360,200,nt);
  # grid size
	isd = 1
	jsd = 1
	ied = ni
	jed = nj

	# time loop
	for time=1:nt-1
		t = time
		wrk1 = zeros(Float64,ni,nj)
	  wrk1 = -horz_advect_tracer_upwind(U,V,Tfield[:,:,t],ni,nj)
		for j=jsd:jed
			for i=isd:ied
					Ttendency[i,j,t] = Ttendency[i,j,t] + wrk1[i,j]
			end
		end

		for j=jsd:jed
	    for i=isd:ied
	      Tfield[i,j,time+1]=(Tfield[i,j,time] + dtime.*Ttendency[i,j,time])
	    end
	  end
	end

	# return
	return Tfield[:,:,nt]
end


function horz_advect_tracer_upwind(uvel,vvel,Tracer_field,ni,nj)
	isd = 1
	jsd = 1
	ied = ni
	jed = nj

	fe = zeros(Float64,ni,nj)
	fn = zeros(Float64,ni,nj)
	upwind = zeros(Float64,ni,nj)

	# i-flux
	for j=jsd+1:jed
		for i=isd:ied-1
			velocity = 0.5*uvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[i+1,j]) .*GRD["lmask"][i,j].*GRD["lmask"][i+1,j]
		end
	end
	# j-flux
	for j=jsd:jed-1
		for i=isd+1:ied
			velocity = 0.5*vvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer_field[i,j] + uneg.*Tracer_field[i,j+1]) .*GRD["lmask"][i,j].*GRD["lmask"][i,j+1]
		end
	end
	# combined
	for j=jsd+1:jed
		for i=isd+1:ied
			upwind[i,j] = GRD["lmask"][i,j].*(fe[i,j]-fe[i-1,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
		end
	end

	# return
	return upwind
end
