### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 069375b9-3662-4281-9180-253392d5370a
begin
	LOAD_PATH = ["@", "@stdlib"]
	using Pkg
	Pkg.activate(".")
	import DelimitedFiles, Mmap, FFTW, HDF5
end

# ╔═╡ b0209ff4-73e1-4cda-b672-26d201054178
md"""# Processing DNS Snapshots"""

# ╔═╡ aec4a0c4-22ca-11ed-1949-a38e3057bf79
const cfg = (
	# set here or in environment to overwrite HDF5 files
	writeh5 = lowercase(get(ENV, "WRITE", "false")) in ("true", "yes", "1"),
	jlpath = "../../data/dns-julia",
	lmpath = "../../data/dns-lee-moser-2015",
	outfile = "../../data/dns-profiles.h5",
 	subset = 41:1:200, # takes about 5 minutes
	parameters = (η = 0.97, ν = 3.5e-4, uτ = 0.0637309),
);

# ╔═╡ 605883de-441e-41dd-b63a-ec67367e3c41
md"## Write profiles to HDF5 file"

# ╔═╡ 25b96178-052d-4b11-a5ba-2b83ce52c753
h5write(fn, lmdata, jldata) = HDF5.h5open(fn, "w") do h5
	h5lm = HDF5.create_group(h5, "lee-moser-2015")
	h5jl = HDF5.create_group(h5, "abl-jl")
	for k in keys(lmdata)
		h5lm[string(k)] = collect(lmdata[k])
	end
	for k in keys(jldata)
		h5jl[string(k)] = collect(jldata[k])
	end
	"wrote `$fn`"
end

# ╔═╡ 554b22fc-4e6c-4716-b778-826e26be3ae4
md"## Read profiles of Lee & Moser (2015) data"

# ╔═╡ 7ac73a2f-215d-4b4a-bb43-4ba881d26f9f
function load_lee_moser_data(dir)
	
	readdat(fn) =  DelimitedFiles.readdlm("$dir/$fn", comments=true, comment_char='%')
	
	cfg = readlines("$dir/LM_Channel_0180_mean_prof.dat")[33:43]
	@assert contains(cfg[7], "Kinematic Viscosity")
	ν = parse(Float64, split(cfg[7], '=')[end])
	@assert contains(cfg[10], "Friction vel")
	uτ = parse(Float64, split(cfg[10], '=')[end])
	
	# y/delta, y^+, U, dU/dy, W, P
	meanprof = readdat("LM_Channel_0180_mean_prof.dat")

	# y/delta, y^+, u'u', v'v', w'w', u'v', u'w', v'w', k
	flucprof = readdat("LM_Channel_0180_vel_fluc_prof.dat")
	
	# y/delta, y^+, Production, Turbulent_Transport, Viscous_Transport,
	# Pressure_Strain, Pressure_Transport, Viscous_Dissipation, Balance
	tkeprof = readdat("LM_Channel_0180_RSTE_k_prof.dat")
	
	E1, E2, k1, k2 = HDF5.h5open("$dir/LM_Channel_0180_1d_energy_spectra.h5") do f
		k = "kx"
		# E1 seems to be defined too large, adds up to 4× TKE → divide by 2
		E1 = (read(f, "Euu_$k") .+ read(f, "Evv_$k") .+ read(f, "Eww_$k")) / 2
		k1 = read(f, k)
		k = "kz"
		E2 = read(f, "Euu_$k") .+ read(f, "Evv_$k") .+ read(f, "Eww_$k")
		k2 = read(f, k)
		E1, E2, k1, k2
	end
	
	(	x3p = meanprof[:,2],
		u1 = meanprof[:,3],
		τdiff = meanprof[:,4], # in plus-units → contains ν
		τturb = -flucprof[:,6],
		tke = flucprof[:,3] .+ flucprof[:,4] .+ flucprof[:,5],
		prod = tkeprof[:,3],
		diss = tkeprof[:,8],
		E1=E1/uτ^2, E2=E2/uτ^2, k1p=k1*ν/uτ, k2p=k2*ν/uτ,
		ν=ν, uτ=uτ,
	)
end

# ╔═╡ eedf8ccb-156c-4544-b5a4-b429aab9b768
md"## Build profiles of Julia data"

# ╔═╡ 67f06956-7ef2-4a89-8348-a64bfabdc0a0
function spectra(vel, include_mke = false)
	# compute spectra of ui·ui along κ1 and κ2

	N = size(vel)
	κ1max = div(N[1]-1,2)
	κ2max = div(N[2]-1,2)
	uhat = FFTW.rfft(vel, (1,2)) / prod(N[1:2])
	Ek1 = zeros(div(N[1]+1,2), N[3])
	Ek2 = zeros(div(N[2]+1,2), N[3])
	for i3 = 1:N[3]
		
		# handle MKE: κ1=κ2=0
		if include_mke
			mke = abs2(uhat[1,1,i3])
			Ek1[1,i3] += mke
			Ek2[1,i3] += mke
		end
		
		# E(κ1): handle κ1=0, contributions from  κ2>0 (skipping MKE)
		for κ2=1:κ2max
			Ek1[1,i3] += abs2(uhat[1,1+κ2,i3]) # κ2 positive
			Ek1[1,i3] += abs2(uhat[1,end+1-κ2,i3]) # κ2 negative
		end
		
		# E(κ1): handle κ1≠0
		for i1=2:size(Ek1,1)
			Ek1[i1,i3] = 2*abs2(uhat[i1,1,i3]) # ×2 for ±κ1, κ2=0
			for κ2=1:κ2max
				Ek1[i1,i3] += 2*abs2(uhat[i1,1+κ2,i3]) # ×2 for ±κ1
				Ek1[i1,i3] += 2*abs2(uhat[i1,end+1-κ2,i3]) # ×2 for ±κ1
			end
		end

		# E(κ2): handle κ2=0, contributions from  κ1>0 (skipping MKE)
		for κ1=1:κ1max
			Ek2[1,i3] += 2*abs2(uhat[1+κ1,1,i3]) # ×2 for ±κ1
		end

		# E(κ2): handle κ2≠0
		for κ2=1:κ2max
			Ek2[1+κ2,i3] += abs2(uhat[1,1+κ2,i3]) # κ1=0
			Ek2[1+κ2,i3] += abs2(uhat[1,end+1-κ2,i3]) # κ1=0
			for κ1=1:κ1max
				Ek2[1+κ2,i3] += 2*abs2(uhat[1+κ1,1+κ2,i3]) # ×2 for ±κ1
				Ek2[1+κ2,i3] += 2*abs2(uhat[1+κ1,end+1-κ2,i3]) # ×2 for ±κ1
			end
		end
			
	end
	Ek1, Ek2
end

# ╔═╡ 3e38ac64-709a-4e68-9aa1-de47c06b11e1
function dh(u; l1=nothing, l2=nothing)
	uhat = FFTW.rfft(u, (1,2)) / prod(size(u)[1:2])
	kxmax, kymax = div.(size(u)[1:2] .- 1, 2)
	kx, ky = [0:kxmax; 0], [0:kymax; 0; -kymax:-1] # nyquist set to zero
	ux = FFTW.brfft(uhat .* (2*π*1im/l1) .* reshape(kx, :, 1, 1), size(u,1), (1,2))
	uy = FFTW.brfft(uhat .* (2*π*1im/l2) .* reshape(ky, 1, :, 1), size(u,1), (1,2))
	ux, uy
end

# ╔═╡ da1e9ab7-40ba-4b1d-835f-74d0a9f72b99
function mmap_cbd(fn; verbose = false)
	verbose && println("mapping file '$fn'...")
	isfile(fn) || error("Cannot open file $fn")
	open(fn, "r") do f
		id = read(f, UInt64)
		T = id == 288230376151834571 ? Float32 :
			id == 576460752303546315 ? Float64 :
			error("File ID does not match number type")
		N = Tuple(read(f, UInt64) for i=1:3)
		xmin = Tuple(read(f, Float64) for i=1:3)
		xmax = Tuple(read(f, Float64) for i=1:3)
		x1, x2, x3 = Tuple(read!(f, zeros(Float64, n)) for n=N)
		x1, x2, x3, Mmap.mmap(f, Array{T,3}, N)
	end
end

# ╔═╡ 40cea051-e8f7-40a6-8b14-29373d237dee
function load_julia_dns_data(path, snaps; η = nothing, ν = nothing, uτ = nothing)
	path = joinpath(path, "snapshots")
	ids = sort(map(f -> split(f, '-')[end], readdir(path)))[snaps]

	x1, x2, x3c = mmap_cbd("$path/state-$(ids[1])/vel1.cbd")[1:3]
	x3i = mmap_cbd("$path/state-$(ids[1])/vel3.cbd")[3]
	N = length.((x1, x2, x3c))

	lτ = ν / uτ
	
	# vertical derivative
	x1length, x2length, x3length = (4π, 2π, 2)
	Dx3 = ζ -> x3length * η*π/2 / sin(η*π/2) * cos(η*(2*ζ-1)*π/2)
	ζc = LinRange(0, 1, 2*N[3]+1)[2:2:end-1]
	ζi = LinRange(0, 1, 2*N[3]+1)[1:2:end]
	Δζ = 1/N[3]

	profiles = (
		u1 = zeros(N[3]),
		u13i = zeros(N[3]+1),
		u1u3 = zeros(N[3]),
		u1u3i = zeros(N[3]+1),
		u1u1 = zeros(N[3]),
		u2u2 = zeros(N[3]),
		u3u3i = zeros(N[3]+1),
		dissc = zeros(N[3]),
		dissi = zeros(N[3]+1),
		E11k1 = zeros(div(N[1]+1,2), N[3]),
		E22k1 = zeros(div(N[1]+1,2), N[3]),
		E33k1 = zeros(div(N[1]+1,2), N[3]),
		E11k2 = zeros(div(N[2]+1,2), N[3]),
		E22k2 = zeros(div(N[2]+1,2), N[3]),
		E33k2 = zeros(div(N[2]+1,2), N[3]),
	)
	
	uij = zeros(N[1], N[2], N[3])
	dx3(u) = begin # c→c, assumes homog. Dirichlet BCs
		@. @views uij[:,:,2:end-1] = (u[:,:,3:end] - u[:,:,1:end-2]) / (2*Δζ)
		# u'(Δζ/2) = (-4 u(0) + 3 u(Δζ/2) + u(3Δζ/2)) / (3Δζ)
		@. @views uij[:,:,1] = (3*u[:,:,1] + u[:,:,2]) / (3*Δζ)
		@. @views uij[:,:,end] = -(3*u[:,:,end] + u[:,:,end-1]) / (3*Δζ)
		uij ./= reshape(Dx3.(ζc), 1, 1, :)
	end
	
	uiji = zeros(N[1], N[2], N[3]+1)
	dx3i(u) = begin # c→i, assumes homog. Dirichlet BCs
		@. @views uiji[:,:,2:end-1] = (u[:,:,2:end] - u[:,:,1:end-1]) / Δζ
		# u'(0) = (-8 u(0) + 9 u(Δζ/2) - u(3Δζ/2)) / (3Δζ)
		@. @views uiji[:,:,1] = (9*u[:,:,1] - u[:,:,2]) / (3*Δζ)
		@. @views uiji[:,:,end] = -(9*u[:,:,end] - u[:,:,end-1]) / (3*Δζ)
		uiji ./= reshape(Dx3.(ζi), 1, 1, :)
	end

	uijc = zeros(N[1], N[2], N[3])
	dx3c(u) = begin # i→c assumes homog. Dirichlet BCs
		@. uijc[:,:,end] = 0
		@. uijc[:,:,1:end-1] = u / Δζ # above
		@. uijc[:,:,2:end] -= u / Δζ # below
		uijc ./= reshape(Dx3.(ζc), 1, 1, :)
	end
	
	havg(x) = sum(x, dims=(1,2))[:] / (N[1] * N[2])
	u3c = zeros(N[1], N[2], N[3])
	u3i = zeros(N[1], N[2], N[3]+1)
	u1i = zeros(N[1], N[2], N[3]+1)

	for id in ids
		dir = "$path/state-$id"
		u1, u2, u3 = [mmap_cbd("$dir/vel$i.cbd")[end] for i=1:3]

		# mean velocity
		profiles.u1 .+= havg(u1)

		@. u3c[:,:,1] = 0 # initialize lowest level
		@. u3c[:,:,2:end] = u3 / 2 # u3 above
		@. u3c[:,:,1:end-1] += u3 / 2 # u3 below
		@. u3i[:,:,2:end-1] = u3
		@. @views u1i[:,:,2:end-1] = u1[:,:,1:end-1]/2 # below
		@. @views u1i[:,:,2:end-1] += u1[:,:,2:end]/2 # above

		# kinetic energy
		profiles.u1u1 .+= havg(u1 .* u1)
		profiles.u2u2 .+= havg(u2 .* u2)
		profiles.u3u3i .+= havg(u3i .* u3i)
		
		# turbulent transport
		profiles.u1u3 .+= havg(u1 .* u3c)
		profiles.u1u3i .+= havg(u1i .* u3i)

		# dissipation i-nodes
		u13i = dx3i(u1) # uses uiji
		profiles.u13i .+= havg(u13i)
		profiles.dissi .+= havg(u13i.^2)
		u23i = dx3i(u2) # uses uiji
		profiles.dissi .+= havg(u23i.^2)
		u31i, u32i = dh(u3i, l1=x1length, l2=x2length)
		profiles.dissi .+= havg(u31i.^2)
		profiles.dissi .+= havg(u32i.^2)
		
		# dissipation on c-nodes
		u11, u12 = dh(u1, l1=x1length, l2=x2length)
		profiles.dissc .+= havg(u11.^2)
		profiles.dissc .+= havg(u12.^2)
		u21, u22 = dh(u2, l1=x1length, l2=x2length)
		profiles.dissc .+= havg(u21.^2)
		profiles.dissc .+= havg(u22.^2)
		u33 = dx3c(u3) # uses uijc
		profiles.dissc .+= havg(u33.^2)

		# spectra
		s = spectra(u1)
		profiles.E11k1 .+= s[1]
		profiles.E11k2 .+= s[2]
		s = spectra(u2)
		profiles.E22k1 .+= s[1]
		profiles.E22k2 .+= s[2]
		s = spectra(u3)
		profiles.E33k1[:,1:end-1] .+= s[1] ./ 2
		profiles.E33k1[:,2:end] .+= s[1] ./ 2
		profiles.E33k2[:,1:end-1] .+= s[2] ./ 2
		profiles.E33k2[:,2:end] .+= s[2] ./ 2
	end

	wrap(a, mirror = false) = let N = length(a), n = div(N, 2), α = mirror ? -1 : 1
		iseven(N) ? (a[1:n] .+ α * a[end:-1:n+1]) / 2 :
			       	[(a[1:n] .+ α * a[end:-1:n+2]) / 2; a[n+1]]
	end

	# normalize by number of samples
	for k in keys(profiles)
		profiles[k] ./= length(snaps)
	end

	# compute/interpolate fields
	i2c(x) = @. x[1:end-1]/2 + x[2:end]/2
	tke = (profiles.u1u1 .+ profiles.u2u2 .+ i2c(profiles.u3u3i) .- profiles.u1.^2)/2
	prod = - i2c(profiles.u1u3i .* profiles.u13i)
	diss = ν * (profiles.dissc .+ i2c(profiles.dissi .- profiles.u13i.^2))
	E1 = (profiles.E11k1 .+ profiles.E22k1 .+ profiles.E33k1)
	E2 = (profiles.E11k2 .+ profiles.E22k2 .+ profiles.E33k2)

	# evaluate prod/diss budget
	Δx3c = (1/N[3]) * Dx3.(ζc)
	Δx3i = (1/N[3]) * Dx3.(ζi)
	println("(prod, diss)=", (sum(prod .* Δx3c), sum(diss .* Δx3c)))

	(   x3p = x3c[1:div(N[3], 2)] / lτ,
		k1p = (0:div(N[1]-1, 2)) * (2π/x1length*lτ),
		k2p = (0:div(N[2]-1, 2)) * (2π/x2length*lτ),
		u1 = wrap(profiles.u1) / uτ,
		τturb = wrap(-i2c(profiles.u1u3i), true) / uτ^2,
		τdiff = wrap(ν * i2c(profiles.u13i), true) / uτ^2,
		tke = wrap(tke) / uτ^2,
		prod = wrap(prod) / (uτ^3/lτ),
		diss = wrap(diss) / (uτ^3/lτ),
		E1 = (E1[:,1:div(N[3],2)] .+ E1[:,end:-1:div(N[3],2)+1])' / (2*uτ^2),
		E2 = (E2[:,1:div(N[3],2)] .+ E2[:,end:-1:div(N[3],2)+1])' / (2*uτ^2),
	)
end

# ╔═╡ 22891b8a-381e-4495-a6b4-c6b9163feb18
cfg.writeh5 ? let
	lmdata = load_lee_moser_data(cfg.lmpath)
	jldata = load_julia_dns_data(cfg.jlpath, cfg.subset; cfg.parameters...)
	h5write(cfg.outfile, lmdata, jldata)
end : "skipped writing HDF5 files"

# ╔═╡ Cell order:
# ╟─b0209ff4-73e1-4cda-b672-26d201054178
# ╠═aec4a0c4-22ca-11ed-1949-a38e3057bf79
# ╟─069375b9-3662-4281-9180-253392d5370a
# ╟─605883de-441e-41dd-b63a-ec67367e3c41
# ╠═25b96178-052d-4b11-a5ba-2b83ce52c753
# ╠═22891b8a-381e-4495-a6b4-c6b9163feb18
# ╟─554b22fc-4e6c-4716-b778-826e26be3ae4
# ╟─7ac73a2f-215d-4b4a-bb43-4ba881d26f9f
# ╟─eedf8ccb-156c-4544-b5a4-b429aab9b768
# ╟─40cea051-e8f7-40a6-8b14-29373d237dee
# ╟─67f06956-7ef2-4a89-8348-a64bfabdc0a0
# ╟─3e38ac64-709a-4e68-9aa1-de47c06b11e1
# ╟─da1e9ab7-40ba-4b1d-835f-74d0a9f72b99
