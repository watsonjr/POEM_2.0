###! Diffusion
function sub_diffuse(Bio_in,Ns,A,DX,DY)
	Bio_out = similar(Bio_in)
	for I = 1:NX
		# Biomasses in Neighboring cells
		bio = Bio_in[GRD["Neigh"][I]]
		bio_out = zeros(Ns)
		#!
		bm = bio[1]
		bu = bio[2]
		bd = bio[3]
		bl = bio[4]
		br = bio[5]

		# Prefactors
		Alpha = A * dtime / (DX[I]^2)
		Beta  = A * dtime / (DY[I]^2)

		# calculation (2D Forward in Time Centered Scheme)
		DIFF =  (Alpha*(bu+bd)) + (Beta*(bl+br)) + (bm*(1-(2*Alpha)-(2*Beta)))

		# Update biomass
		if DIFF < 0
			bio_out = 0
		else
			bio_out = DIFF
		end
		Bio_out[I] = bio_out
	end
	return Bio_out
end


### ADVECTION ###
function sub_advec2D(Bio_in,Nu,U,V,DX,DY,Tp,Tb,tdif,wgt)
	#Bio_in = .bio; Nu = .nu; DX = GRD["dxtn"]; DY = GRD["dyte"]; U = ENVR.U; V = ENVR.V;
	#Tp = ENVR.Tp; Tb = ENVR.Tb; tdif = Med_f.td; wgt = M_m;

	#! 2D bio
	bio2 = zeros(360,200);
	bio2[GRD["ID"]] = Bio_in;
	nu = zeros(360,200);
	nu[GRD["ID"]] = Nu;
	dx = zeros(360,200);
	dx[GRD["ID"]] = DX;
	dy = zeros(360,200);
	dy[GRD["ID"]] = DY;
	U2 = zeros(360,200);
	U2[GRD["ID"]] = U;
	V2 = zeros(360,200);
	V2[GRD["ID"]] = V;

	#! Calc swimming speed (m/d)
	#T = (Tp.*tdif) + (Tb.*(1.0-tdif))
	#w = ((3.9*wgt.^0.13 * exp(0.149*T)) /100*60*60*24)
	Q = zeros(360,200);
	#Q[GRD["ID"]] = w;

	#! Find direction to swim in
	KK = zeros(Int64,360,200);
	LID = zeros(360,200);
	LID[GRD["ID"]] = GRD["lmask"];
	ni, nj = size(KK);
	for i = 1:ni
		for j = 1:nj
			if (j==1)
				if (i==1)
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[ni-i+1,j],nu[360,j],nu[i+1,j]])[2]
				elseif (i==ni)
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[ni-i+1,j],nu[i-1,j],nu[1,j]])[2]
				else
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[ni-i+1,j],nu[i-1,j],nu[i+1,j]])[2]
				end
			elseif (j==nj)
				if (i==1)
					KK[i,j] = findmax([nu[i,j],nu[ni-i+1,j],nu[i,j-1],nu[360,j],nu[i+1,j]])[2]
				elseif (i==ni)
					KK[i,j] = findmax([nu[i,j],nu[ni-i+1,j],nu[i,j-1],nu[i-1,j],nu[1,j]])[2]
				else
					KK[i,j] = findmax([nu[i,j],nu[ni-i+1,j],nu[i,j-1],nu[i-1,j],nu[i+1,j]])[2]
				end
			else
				if (i==1)
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[i,j-1],nu[360,j],nu[i+1,j]])[2]
				elseif (i==ni)
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[i,j-1],nu[i-1,j],nu[1,j]])[2]
				else
					KK[i,j] = findmax([nu[i,j],nu[i,j+1],nu[i,j-1],nu[i-1,j],nu[i+1,j]])[2]
				end
			end
		end
	end

	#! Calculate effective speed
	u = U2
	v = V2
	I1 = find(KK .== 1)
	I2 = find(KK .== 2)
	I3 = find(KK .== 3)
	I4 = find(KK .== 4)
	I5 = find(KK .== 5)

	u[I1] = U2[I1]
	v[I1] = V2[I1]

	u[I2] = U2[I2]
	v[I2] = V2[I2] + Q[I2]

	u[I3] = U2[I3]
	v[I3] = V2[I3] - Q[I3]

	u[I4] = U2[I4] - Q[I4]
	v[I4] = V2[I4]

	u[I5] = U2[I5] + Q[I5]
	v[I5] = V2[I5]

	#! Land
	# I1 = find(v.>0. + LID[:,2].==0.) 	#up
	# I2 = find(v.<0. + LID[:,3].==0.)	#down
	# I3 = find(u.<0. + LID[:,4].==0.)	#left
	# I4 = find(u.>0. + LID[:,5].==0.)	#right
	# v[I1] *= -1
	# v[I2] *= -1
	# u[I3] *= -1
	# u[I4] *= -1
	u = u.*LID;
	v = v.*LID;

	#! U advection
	Bio_minus = similar(bio2)
	Bio_minus = zeros(360,200);
	for i = 1:ni
		for j = 1:nj
			UU = u[i,j] * dtime / dx[i,j] #Produces NaN when dx=0 at land cells
			if (i==1)
				BB = bio2[i,j] - ((UU/2)*(bio2[i+1,j]-bio2[360,j])) -
								(((UU^2)/2)*(bio2[i+1,j]-(2*bio2[i,j])+bio2[360,j]))
			elseif (i==ni)
				BB = bio2[i,j] - ((UU/2)*(bio2[1,j]-bio2[i-1,j])) -
								(((UU^2)/2)*(bio2[1,j]-(2*bio2[i,j])+bio2[i-1,j]))
			else
				BB = bio2[i,j] - ((UU/2)*(bio2[i+1,j]-bio2[i-1,j])) -
								(((UU^2)/2)*(bio2[i+1,j]-(2*bio2[i,j])+bio2[i-1,j]))
			end
			if BB < 0
				Bio_minus[i,j] = 0
			elseif isnan(BB)
				Bio_minus[i,j] = 0
			else
				Bio_minus[i,j] = BB
			end
		end
	end

	#! V advection
	Bio_plus = similar(bio2)
	Bio_plus = zeros(360,200);
	bio = Bio_minus
	for i = 1:ni
		for j = 1:nj
			VV = v[i,j] * dtime / dy[i,j]
			if (j==1)
				BB = bio[i,j] - ((VV/2)*(bio[i,j+1]-bio[ni-i+1,j])) -
							(((VV^2)/2)*(bio[i,j+1]-(2*bio[i,j])+bio[ni-i+1,j]))
			elseif (j==nj)
				BB = bio[i,j] - ((VV/2)*(bio[ni-i+1,j]-bio[i,j-1])) -
							(((VV^2)/2)*(bio[ni-i+1,j]-(2*bio[i,j])+bio[i,j-1]))
			else
				BB = bio[i,j] - ((VV/2)*(bio[i,j+1]-bio[i,j-1])) -
							(((VV^2)/2)*(bio[i,j+1]-(2*bio[i,j])+bio[i,j-1]))
			end
			if BB < 0
				Bio_plus[i,j] = 0
			elseif isnan(BB)
				Bio_plus[i,j] = 0
			else
				Bio_plus[i,j] = BB
			end
		end
	end

	#! Vectorize
	biov = zeros(NX);
	biov = Bio_plus[GRD["ID"]];

	# return
	return biov
end
