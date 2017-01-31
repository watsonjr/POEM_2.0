###! DIFFUSION & ADVECTION TIME LOOP###-------------------------------
function sub_advec_diff(GRD,Bio_in,K,U,V,ni,nj)
	# K = diffusivity in m2/s
	# U & V = velocities in m/s
	# ntime = time steps in a day
	# dtime = # seconds in ntime
	dtime = 60.0*60.0*0.25
	ntime = (60.0*60.0*24.0) / dtime
	nt = Int(ntime)
	# biol concentration
	Tfield = zeros(Float64,360,200,nt+1);
	Tfield[:,:,1] = Bio_in;
	# grid size
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	# time loop
	for time=1:nt
		t = time
		wrk1 = zeros(Float64,ni,nj);
		wrk2 = zeros(Float64,ni,nj);
	  wrk1, wrk2 = horz_advect_diff_upwind(K,U,V,Tfield[:,:,t],ni,nj)
		Tfield[:,:,time+1] = (Tfield[:,:,time] - dtime.*wrk1 + dtime.*wrk2)
	end

	# return
	return Tfield[:,:,nt+1]

end





function horz_advect_diff_upwind(K,uvel,vvel,Tracer,ni,nj)
	isd = 1
	jsd = 2 #ignore j=1 b/c land (Antarctica)
	ied = ni
	jed = nj

	dfe = zeros(Float64,ni,nj);
	dfn = zeros(Float64,ni,nj);
	dupwind = zeros(Float64,ni,nj);
	afe = zeros(Float64,ni,nj);
	afn = zeros(Float64,ni,nj);
	aupwind = zeros(Float64,ni,nj);

	# Gradient i
	gradTi = zeros(Float64,ni,nj);
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied
			if (i == ied)
				gradTi[i,j] = (Tracer[isd,j] - Tracer[i,j]) ./ GRD["dyte"][i,j] .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
			else
				gradTi[i,j] = (Tracer[i+1,j] - Tracer[i,j]) ./ GRD["dyte"][i,j] .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
			end
		end
	end
	# Gradient j
	gradTj = zeros(Float64,ni,nj);
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied
			if (j < jed)
				gradTj[i,j] = (Tracer[i,j+1] - Tracer[i,j]) ./ GRD["dxtn"][i,j] .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
			else
				gradTj[i,j] = (Tracer[ni-i+1,j] - Tracer[i,j]) ./ GRD["dxtn"][i,j] .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
			end
		end
	end
	gradT = gradTi + gradTj

	diffusiv = 0.5*K
	kpos     = diffusiv + abs(diffusiv)
	kneg     = diffusiv - abs(diffusiv)
	# i-flux
	for j=jsd:jed
		for i=isd:ied #i=isd-1:ied
			velocity = 0.5*uvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (GRD["lmask"][i,j,1] > 0)
				if (i == ied)
					dfe[i,j]  = GRD["dyte"][i,j].*(kpos.*gradT[i,j] + kneg.*gradT[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
					afe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[isd,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][isd,j,1]
				else
					dfe[i,j]  = GRD["dyte"][i,j].*(kpos.*gradT[i,j] + kneg.*gradT[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
					afe[i,j]  = GRD["dyte"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i+1,j,1]
				end
			end
		end
	end

	# j-flux
	for j=jsd:jed #j=jsd-1:jed
		for i=isd:ied
			velocity = 0.5*vvel[i,j]
			upos     = velocity + abs(velocity)
			uneg     = velocity - abs(velocity)
			if (GRD["lmask"][i,j,1] > 0)
				if (j < jed)
					dfn[i,j]  = GRD["dxtn"][i,j].*(kpos.*gradT[i,j] + kneg.*gradT[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
					afn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[i,j+1]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][i,j+1,1]
				else
					dfn[i,j]  = GRD["dxtn"][i,j].*(kpos.*gradT[i,j] + kneg.*gradT[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
					afn[i,j]  = GRD["dxtn"][i,j].*(upos.*Tracer[i,j] + uneg.*Tracer[ni-i+1,j]) .*GRD["lmask"][i,j,1] .*GRD["lmask"][ni-i+1,j,1]
				end
			end
		end
	end

	# combined
	for j=jsd:jed
		for i=isd:ied
			if (j > 1)
				if (i > 1)
					dupwind[i,j] = GRD["lmask"][i,j,1].*(dfe[i,j]-dfe[i-1,j]+dfn[i,j]-dfn[i,j-1]).*GRD["datr"][i,j]
					aupwind[i,j] = GRD["lmask"][i,j,1].*(afe[i,j]-afe[i-1,j]+afn[i,j]-afn[i,j-1]).*GRD["datr"][i,j]
				else
					dupwind[i,j] = GRD["lmask"][i,j,1].*(dfe[i,j]-dfe[ied,j]+dfn[i,j]-dfn[i,j-1]).*GRD["datr"][i,j]
					aupwind[i,j] = GRD["lmask"][i,j,1].*(afe[i,j]-afe[ied,j]+afn[i,j]-afn[i,j-1]).*GRD["datr"][i,j]
				end
			# no need for else at South Pole
			end
		end
	end

	# return
	return aupwind, dupwind
end
