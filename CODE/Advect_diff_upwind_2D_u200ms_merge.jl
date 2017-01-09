###! DIFFUSION & ADVECTION TIME LOOP###-------------------------------
function sub_advec_diff(GRD,Bio_in,K,U,V,ni,nj)
	# K = diffusivity in m2/s
	# U & V = velocities in m/s
	# ntime = time steps in a day
	# dtime = # seconds in ntime
	dtime = 60.0*60.0*12.0
	ntime = (60.0*60.0*24.0) / dtime
	nt = Int(ntime)
	# biol concentration
	Tfield = zeros(Float64,360,200,nt);
	Tfield[:,:,1] = Bio_in;
	# grid size
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	# time loop
	for time=1:nt-1
		t = time
		wrk1 = zeros(Float64,ni,nj);
	  wrk1 = -horz_advect_tracer_upwind(U,V,Tfield[:,:,t],ni,nj)
		wrk2 = zeros(Float64,ni,nj);
	  wrk2 = horz_diff_upwind(K,Tfield[:,:,t],ni,nj)
		Tfield[:,:,time+1] = (Tfield[:,:,time] + dtime.*wrk1 + dtime.*wrk2)
	end

	# return
	return Tfield[:,:,nt]

end





function horz_diff_upwind(K,Tracer,ni,nj)
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	fe = zeros(Float64,ni,nj);
	fn = zeros(Float64,ni,nj);
	upwind = zeros(Float64,ni,nj);

	# Gradient i
	gradTi = zeros(Float64,ni,nj);
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied
			if (i == ied)
				gradTi[i,j] = (Tracer[isd,j] - Tracer[i,j]) ./ GRD["dyte"][i,j]
			else
				gradTi[i,j] = (Tracer[i+1,j] - Tracer[i,j]) ./ GRD["dyte"][i,j]
			end
		end
	end
	# Gradient j
	gradTj = zeros(Float64,ni,nj);
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied
			if (j < jed)
				gradTj[i,j] = (Tracer[i,j+1] - Tracer[i,j]) ./ GRD["dxtn"][i,j]
			else
				gradTj[i,j] = (Tracer[ni-i+1,j] - Tracer[i,j]) ./ GRD["dxtn"][i,j]
			end
		end
	end
	gradT = gradTi + gradTj

	diffusiv = 0.5*K
	upos     = diffusiv + abs(diffusiv)
	uneg     = diffusiv - abs(diffusiv)
	# i-flux
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied
			if (i == ied)
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*gradT[i,j] + uneg.*gradT[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
			else
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*gradT[i,j] + uneg.*gradT[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
			end
		end
	end

	# j-flux
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied
			if (j < jed)
				fn[i,j]  = GRD["dxtn"][i,j].*(upos.*gradT[i,j] + uneg.*gradT[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
			else
				fn[i,j]  = GRD["dxtn"][i,j].*(upos.*gradT[i,j] + uneg.*gradT[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
			end
		end
	end

	# combined
	for j=jsd:jed
		for i=isd:ied
			if (j > 1)
				if (i > 1)
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[i-1,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				else
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[ied,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				end
			# no need for else at South Pole
			end
		end
	end

	# return
	return upwind
end





function horz_advect_tracer_upwind(uvel,vvel,Tracer,ni,nj)
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	fe = zeros(Float64,ni,nj);
	fn = zeros(Float64,ni,nj);
	upwind = zeros(Float64,ni,nj);

	# i-flux
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied

			velocity = 0.5*uvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (i == ied)
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
			else
				fe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
			end
		end
	end

	# j-flux
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied

			velocity = 0.5*vvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (j < jed)
					fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
			else
					fn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
			end
		end
	end

	# combined
	for j=jsd:jed
		for i=isd:ied
			if (j > 1)
				if (i > 1)
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[i-1,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				else
					upwind[i,j] = GRD["lmask"][i,j,1].*(fe[i,j]-fe[ied,j]+fn[i,j]-fn[i,j-1]).*GRD["datr"][i,j]
				end
			# no need for else at South Pole
			end
		end
	end

	# return
	return upwind
end
