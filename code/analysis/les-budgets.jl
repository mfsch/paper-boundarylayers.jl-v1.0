### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 60a0bfb9-4aab-4831-a07e-a98a4f5ac20c
begin
	LOAD_PATH = ["@", "@stdlib"]
	using Pkg
	Pkg.activate(".")
	import Mmap, PyPlot, FFTW, JSON3, HDF5, BoundaryLayerDynamics
	const BLD = BoundaryLayerDynamics
end;

# ╔═╡ 668e802b-eb57-474c-98cb-86fa661668d9
md"""
# Processing LES Snapshots
"""

# ╔═╡ 67a565f2-45ad-466a-ba1e-e5fb3e53560e
const cfg = (
	# set here or in environment to overwrite HDF5 files
	writeh5 = lowercase(get(ENV, "WRITE", "false")) in ("true", "yes", "1"),
	jlpath = "../../data/les-julia",
	ftpath = "../../data/les-fortran",
	outfile = "../../data/les-profiles.h5",
	dims = (64, 64, 64),
 	subset = 41:200, # takes about 20s
);

# ╔═╡ 76de4f32-bbfe-449c-9a23-b007fcb005e9
md"""## Examine last run"""

# ╔═╡ 1410f160-e7fc-4961-927c-d2de6c2ee835
md"""
## Write profiles to HDF5 files
"""

# ╔═╡ 78e8a13a-4a6f-4153-939f-26be6a98ec37
h5write(fn, ftdata, jldata) = HDF5.h5open(fn, "w") do h5
	h5ft = HDF5.create_group(h5, "albertson-parlange")
	h5jl = HDF5.create_group(h5, "abl-jl")
	for k in keys(ftdata)
		h5ft[string(k)] = collect(ftdata[k])
	end
	for k in keys(jldata)
		h5jl[string(k)] = collect(jldata[k])
	end
	"wrote `$fn`"
end

# ╔═╡ 6b769ef4-8d5e-4638-8ebf-e018d6b3ce4f
md"""
## Build mean profiles from Fortran data
"""

# ╔═╡ 025e3866-3f96-41fc-822e-0fc78bc9eafb
md"""
#### Build mean profiles from Julia data
"""

# ╔═╡ 12068898-667d-4947-a48f-b489317a5cab
function add_vel1(umean, u)
	umean .+= sum(u, dims=(1,2))[:] / prod(size(u)[1:2])
end

# ╔═╡ 53e4fe37-3f17-4244-bd7e-56637ec4dffa
function add_vel1vel3(uwmean, u, w)
	nh = prod(size(u)[1:2])
	uwmean[2:end]   .+= sum(u[:,:,2:end]   .* w, dims=(1,2))[:] / (2*nh) # u*wbelow
	uwmean[1:end-1] .+= sum(u[:,:,1:end-1] .* w, dims=(1,2))[:] / (2*nh) # u*wabove
end

# ╔═╡ e67c7981-3bc3-4bc8-9400-97f99aa1eef9
function add_ke(kemean, u, v, w)
	nh = prod(size(u)[1:2])
	wh = zeros(size(u))
	wh[:,:,1:end-1] .+= w/2 # w above
	wh[:,:,2:end] .+= w/2 # w below
	kemean .+= sum(u .^ 2, dims=(1,2))[:] / (2*nh)
	kemean .+= sum(v .^ 2, dims=(1,2))[:] / (2*nh)
	kemean .+= sum(wh .^ 2, dims=(1,2))[:] / (2*nh)
end

# ╔═╡ 02ea41ac-a1d2-4d98-8aa5-9a69a3f400ef
struct TermBuffer
	terms
	buffer
	TermBuffer(terms) = new(terms, Dict{Symbol,Array}())
end

# ╔═╡ 82a546e8-75d8-4c7c-abfa-4d4d2412ae00
 function BLD.Logging.log_sample!(log::TermBuffer, (field, values)::Pair, t; kwargs...)
	 field in log.terms || return
	 n3 = size(values, 3) + (haskey(kwargs, :bcs) ? 2 : 0)
	 if !haskey(log.buffer, field)
		 log.buffer[field] = zeros(size(values)[1:2]..., n3)
	 end
	 vals(bc::BLD.BoundaryConditions.ConstantValue) = bc.value
	 vals(bc::BLD.BoundaryConditions.DynamicValues) = bc.values
	 if haskey(kwargs, :bcs)
		 log.buffer[field][:,:,1] .= vals(kwargs[:bcs][1].type)
		 log.buffer[field][:,:,2:end-1] .= values
		 log.buffer[field][:,:,end] .= vals(kwargs[:bcs][2].type)
	 else
		 log.buffer[field] .= values # only retain first log
	 end
 end

# ╔═╡ 8df38d6b-43b8-43f6-a507-f766d5fef53d
function init_sgs(; dims=(64,64,64), l1 = 2π, l2 = 4/3*π, z0 = 1e-4)
	dom = BLD.Domain((l1, l2, 1), BLD.RoughWall(z0), BLD.FreeSlipBoundary())
	proc = [BLD.StaticSmagorinskyModel()]
	abl = BLD.Model(dims, dom, proc)
	rhs = deepcopy(abl.state)
	output = TermBuffer((:sgs11, :sgs12, :sgs13, :sgs22, :sgs23, :sgs33))

	interp(a) = if size(a, 3) == dims[3]
		a
	elseif size(a, 3) == dims[3] + 1
		(a[:,:,2:end] .+ a[:,:,1:end-1]) / 2
	else
		ai = zeros(dims)
		ai[:,:,1:end-1] .+= a/2 # above
		ai[:,:,2:end] .+= a/2 # below
		ai
	end
	hmean(a) = sum(a, dims=(1,2))[:] / prod(size(a)[1:2])
	hprod(a, b) = sum(interp(a) .* interp(b), dims=(1,2))[:] / prod(size(a)[1:2])

	function sgs(sgs13, diss, u13, u, v, w)
		# computes all instantaneous contributions from sgs stresses
		BLD.initialize!(abl, vel1 = u, vel2 = v, vel3 = w)
		BLD.evolve!(abl, 0.0, dt = 0.0, output = output, method = BLD.Euler())
		sgs = output.buffer
		terms = abl.physical_spaces[size(u)[1:2]].terms
		sgs13 .+= hmean(sgs[:sgs13])
		u13 .+= hmean(terms.vel1_3.values)
		diss .-= hprod(sgs[:sgs11], terms.vel1_1.values) .+
				 hprod(sgs[:sgs22], terms.vel2_2.values) .+
				 hprod(sgs[:sgs33], terms.vel3_3.values) .+
				 hprod(sgs[:sgs12], terms.vel1_2.values) .+
				 hprod(sgs[:sgs12], terms.vel2_1.values) .+
				 hprod(sgs[:sgs13], terms.vel3_1.values) .+
				 hprod(sgs[:sgs23], terms.vel3_2.values)

		# compute terms with dui/dx3 gradients separately
		dissi3 = hprod(sgs[:sgs13], terms.vel1_3.values) .+
				 hprod(sgs[:sgs23], terms.vel2_3.values)
		# overwrite wall value (at z=dz/2)
		@assert size(sgs[:sgs13], 3) == dims[3] + 1
		sgs13w = (sgs[:sgs13][:,:,1] .+ sgs[:sgs13][:,:,2])/2
		sgs23w = (sgs[:sgs23][:,:,1] .+ sgs[:sgs23][:,:,2])/2
		z1 = 1/(2*dims[3])
		dudzw = u[:,:,1] / (z1 * log(z1/z0))
		dvdzw = v[:,:,1] / (z1 * log(z1/z0))
		dissi3[1] = hmean(sgs13w .* dudzw)[] .+ hmean(sgs23w .* dvdzw)[]
		diss .-= dissi3 # add to mean
	end
end

# ╔═╡ a11e0947-5451-4e25-8894-0e0be324ef8d
md"""
## Examine mean profiles
"""

# ╔═╡ c3a61806-dd18-4ea6-b267-a5bf2ee07ad0
md"""
## Shared helper functions
"""

# ╔═╡ c4c06d18-ce1d-45f2-9376-a6309e3e4f5a
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
		data = Mmap.mmap(f, Array{T,3}, N)
	end
end

# ╔═╡ 814af8ef-64ba-4a46-ab47-9ef4d17909c3
function dh(u; l1=nothing, l2=nothing)
	uhat = FFTW.rfft(u, (1,2)) / prod(size(u)[1:2])
	kxmax, kymax = div.(size(u)[1:2] .- 1, 2)
	kx, ky = [0:kxmax; 0], [0:kymax; 0; -kymax:-1] # nyquist set to zero
	ux = FFTW.brfft(uhat .* (2*π*1im/l1) .* reshape(kx, :, 1, 1), size(u,1), (1,2))
	uy = FFTW.brfft(uhat .* (2*π*1im/l2) .* reshape(ky, 1, :, 1), size(u,1), (1,2))
	ux, uy
end

# ╔═╡ 6633423f-8961-43fe-8329-33192b1b142d
function spectra(vel)

	N = size(vel)
	κ1max = div(N[1]-1,2)
	κ2max = div(N[2]-1,2)
	uhat = FFTW.rfft(vel, (1,2)) / prod(N[1:2])
	Ek1 = zeros(div(N[1]+1,2), N[3])
	Ek2 = zeros(div(N[2]+1,2), N[3])
	for i3 = 1:N[3]
		
		# handle κ1 = 0 → remove MKE
		#E11[1,i3] = abs2(uhat[1,1,i3]) # contribution from κ2=0
		# contribution from  κ2>0
		for κ2=1:κ2max
			Ek1[1,i3] += abs2(uhat[1,1+κ2,i3]) # κ1=0
			Ek1[1,i3] += abs2(uhat[1,end+1-κ2,i3]) # κ1=0
		end
		
		for i1=2:size(Ek1,1)
			κ1 = i1-1
			Ek1[i1,i3] = 2*abs2(uhat[i1,1,i3]) # ×2 for ±κ1, κ2=0
			for κ2=1:κ2max
				Ek1[i1,i3] += 2*abs2(uhat[i1,1+κ2,i3]) # ×2 for ±κ1
				Ek1[i1,i3] += 2*abs2(uhat[i1,end+1-κ2,i3]) # ×2 for ±κ1
			end
		end

		# handle κ2 = 0 → remove MKE
		for κ1=1:κ1max
			Ek2[1,i3] += 2*abs2(uhat[1+κ1,1,i3]) # ×2 for ±κ1
		end

		for κ2=1:κ2max
			# κ1=0, κ2≠0
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

# ╔═╡ 3359d7ed-c99e-4a1e-9575-67cd327ac763
function load_profiles_ft(path, subset; dims = nothing, z0 = 1e-4, l1=2π, l2=4π/3, verbose=true)
	fpath = "$path/instantaneous-fields"
	snapshots = map(x -> x[3:5], filter(startswith("u-"), readdir(fpath)))[subset]

	profiles = (
		u = zeros(dims[3]),
		uw = zeros(dims[3]),
		sgs13 = zeros(dims[3]+1),
		tke = zeros(dims[3]),
		prod = zeros(dims[3]),
		diss = zeros(dims[3]),
		Ek1 = zeros(div(dims[1]-1, 2), dims[3]),
		Ek2 = zeros(div(dims[2]-1, 2), dims[3]),
	)

	havg(x) = sum(x, dims=(1,2))[:] / prod(size(x)[1:2])

	for s in snapshots
		u, v, w, t11, t12, t13, t22, t23, t33 = (mmap_cbd("$fpath/$f-$s.cbd")
			for f in ("u", "v", "w", "txx", "txy", "txz", "tyy", "tyz", "tzz"))
		wext = zeros(dims[1:2]..., dims[3]+1)
		wext[:,:,1:end-1] .= w

		# simple profiles
		profiles.u .+= havg(u)
		wc = (wext[:,:,1:end-1] .+ wext[:,:,2:end]) / 2
		profiles.uw .+= havg(u .* wc)
		# tij has negative sign in LES output
		profiles.sgs13[1:end-1] .-= havg(t13)
		profiles.tke .+= (havg(u.^2) .+ havg(v.^2) .+ havg(wc.^2))/2

		# spectra
		E11k1, E11k2 = spectra(u)
		E22k1, E22k2 = spectra(v)
		E33k1, E33k2 = spectra(wc)
		profiles.Ek1 .+= E11k1[2:end,:] / 2
		profiles.Ek1 .+= E22k1[2:end,:] / 2
		profiles.Ek1 .+= E33k1[2:end,:] / 2
		profiles.Ek2 .+= E11k2[2:end,:] / 2
		profiles.Ek2 .+= E22k2[2:end,:] / 2
		profiles.Ek2 .+= E33k2[2:end,:] / 2

		# horizontal derivatives
		u11, u12 = dh(u, l1=l1, l2=l2)
		u21, u22 = dh(v, l1=l1, l2=l2)
		u31, u32 = dh(wext, l1=l1, l2=l2)
		u31c = (u31[:,:,1:end-1] .+ u31[:,:,2:end]) / 2
		u32c = (u32[:,:,1:end-1] .+ u32[:,:,2:end]) / 2

		# vertical derivatives
		dz = 1/dims[3]
		u13 = (u[:,:,2:end] .- u[:,:,1:end-1]) / dz
		u13c = zeros(size(u))
		u13c[:,:,1:end-1] .+= u13/2 # above
		u13c[:,:,2:end]   .+= u13/2 # below
		u13c[:,:,1] .= u[:,:,1] / (0.5*dz * log(0.5*dz/z0))
		u23 = (v[:,:,2:end] .- v[:,:,1:end-1]) / dz
		u23c = zeros(size(v))
		u23c[:,:,1:end-1] .+= u23/2 # above
		u23c[:,:,2:end]   .+= u23/2 # below
		u23c[:,:,1] .= v[:,:,1] / (0.5*dz * log(0.5*dz/z0))
		u33c = (wext[:,:,2:end] .- wext[:,:,1:end-1]) / dz

		# interpolated stresses
		t13c = t13/2 # t13 below
		t13c[:,:,1:end-1] .+= t13[:,:,2:end]/2 # t13 above
		t23c = t23/2 # t23 below
		t23c[:,:,1:end-1] .+= t23[:,:,2:end]/2 # t23 above

		# dissipation term (tij has negative sign in LES output)
		profiles.diss .+= havg(u11  .* t11 ) .+
						  havg(u12  .* t12 ) .+
						  havg(u13c .* t13c) .+
						  havg(u21  .* t12 ) .+
						  havg(u22  .* t22 ) .+
						  havg(u23c .* t23c) .+
						  havg(u31c .* t13c) .+
						  havg(u32c .* t23c) .+
						  havg(u33c .* t33 )
	end

	for k in keys(profiles)
		k in (:zc, :zi) && continue
		profiles[k] ./= length(subset)
	end

	dz = 1/dims[3]
	dudz = [0.0; diff(profiles.u) / dz; 0.0]
	dudz = (dudz[1:end-1] .+ dudz[2:end])/2
	dudz[1] = profiles.u[1] / (0.5*dz * log(0.5*dz/z0))
	profiles.prod .= - profiles.uw .* dudz
	profiles.prod .*= z0
	profiles.tke .-= profiles.u.^2/2
	t13c = -(profiles.sgs13[1:end-1] .+ profiles.sgs13[2:end])/2
	profiles.diss .-= t13c .* dudz
	profiles.diss .*= z0
	
	if verbose
		println("FT dudz wall: ", (1/0.4 * 2/dz, dudz[1]))
		println("FT prod, diss: ", (sum(profiles.prod), sum(profiles.diss)))
	end

	profiles
end

# ╔═╡ 294cbe47-4d5c-48e2-8599-77ab84931824
function load_profiles_jl(path, subset; dims = nothing, z0 = 1e-4, verbose=true)
	snapshots = readdir("$path/snapshots")[subset]

	profiles = (
		zc = LinRange(0, 1, 1+2*dims[3])[2:2:end] / z0,
		zi = LinRange(0, 1, 1+2*dims[3])[1:2:end] / z0,
		k1 = 1:div(dims[1]-1, 2),
		k2 = 1:div(dims[2]-1, 2),
		u = zeros(dims[3]),
		uw = zeros(dims[3]),
		sgs13 = zeros(dims[3]+1),
		tke = zeros(dims[3]),
		prod = zeros(dims[3]),
		diss = zeros(dims[3]),
		u13 = zeros(dims[3]-1),
		Ek1 = zeros(div(dims[1]-1, 2), dims[3]),
		Ek2 = zeros(div(dims[2]-1, 2), dims[3]),
	)

	add_sgs = init_sgs()
	
	for s in snapshots
		spath = joinpath(path, "snapshots", s)
		u, v, w = (mmap_cbd(joinpath(spath, vel))
			for vel in readdir(spath))
		add_vel1(profiles.u, u)
		add_vel1vel3(profiles.uw, u, w)
		add_ke(profiles.tke, u, v, w)
		add_sgs(profiles.sgs13, profiles.diss, profiles.u13, u, v, w)
		wc = zeros(size(u))
		wc[:,:,1:end-1] .+= w/2 # above
		wc[:,:,2:end] .+= w/2 # below
		E11k1, E11k2 = spectra(u)
		E22k1, E22k2 = spectra(v)
		E33k1, E33k2 = spectra(wc)
		profiles.Ek1 .+= E11k1[2:end,:] / 2
		profiles.Ek1 .+= E22k1[2:end,:] / 2
		profiles.Ek1 .+= E33k1[2:end,:] / 2
		profiles.Ek2 .+= E11k2[2:end,:] / 2
		profiles.Ek2 .+= E22k2[2:end,:] / 2
		profiles.Ek2 .+= E33k2[2:end,:] / 2
	end
	
	for k in keys(profiles)
		k in (:zc, :zi, :k1, :k2) && continue
		profiles[k] ./= length(subset)
	end

	# compute/subtract mean products
	interp(p) = (p[2:end] .+ p[1:end-1]) / 2
	zc1 = 1 / (2*dims[3])
	dudz_wall = profiles.u[1] / (zc1 * log(zc1/z0))
	dudzc = [dudz_wall; interp([profiles.u13; 0])]
	profiles.prod .= - profiles.uw .* dudzc
	profiles.diss .-= - dudzc .* interp(profiles.sgs13)
	profiles.tke .-= profiles.u .^ 2 / 2

	# normalization, assuming ut=1
	profiles.prod .*= z0
	profiles.diss .*= z0
	
	if verbose
		println("JL dudz wall: ", (1/0.4 * 1/zc1, dudz_wall))
		println("JL prod, diss: ", (sum(profiles.prod), sum(profiles.diss)))
	end

	profiles
end

# ╔═╡ 4c41c8ae-e41e-44a0-a54f-df0363d47e70
cfg.writeh5 ? let
	ftdata = load_profiles_ft(cfg.ftpath, cfg.subset, dims = cfg.dims)
	jldata = load_profiles_jl(cfg.jlpath, cfg.subset, dims = cfg.dims)
	h5write(cfg.outfile, ftdata, jldata)
end : "skipped writing HDF5 files"

# ╔═╡ 9ccdea97-69ec-4992-9f1c-613cddbcefe2
function plt(f::Function)
    PyPlot.close("all") # clean up if the last call has crashed
	
	# font setup
	rcp = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcp["font.family"] = "monospace"
	rcp["font.monospace"] = ["Source Code Pro"]

	f(PyPlot)
    figs = [PyPlot.figure(num) for num in PyPlot.get_fignums()] # force order
    PyPlot.close("all")
    length(figs) == 1 ? figs[1] : figs
end

# ╔═╡ 926c81c7-5aef-449b-ad9f-4bf4b3e00064
function examine_ft_log(path)
	l = readlines(path)
	vals(marker, field) = map(filter(l -> contains(l, marker), l)) do l
		parse(Float64, split(l)[field])
	end
	t = vals("dt,tt", 4)
	ustar = vals("u_star_havg", 3)
	cfl = vals("cfl,cfl_old", 2)
	ke = vals("divu,ke", 6)

	plt() do P
		P.figure(figsize=(8,3))
		P.plot(t, ustar / maximum(ustar))
		P.plot(t, ke / maximum(ke))
		P.plot(t, cfl * 10)
		P.xlim(0, t[end])
		P.ylim(0,1.1)
	end
end

# ╔═╡ fcea0662-c01d-4a6c-a7cd-7ba3dcc02d7e
let fn = cfg.ftpath * "/.cache/output.log"
	isfile(fn) && examine_ft_log(fn)
end

# ╔═╡ 6fc8046f-fb19-46a0-be95-9e70fba8f11b
function examine_jl_log(path)
	l = readlines(path)
	vals(marker, field) = map(filter(l -> contains(l, marker), l)) do l
		parse(Float64, split(split(l)[field], ',')[1])
	end
	t = vals("Simulation Time", 4)
	ustar = sqrt.(-vals("τ₁₃", 5))
	c = vals("Courant", 5)
	tke = vals("Turbulent KE", 4)
	mke = vals("Mean KE", 4)
	ke = mke .+ tke

	plt() do P
		P.figure(figsize=(8,3))
		P.plot(t, ustar / maximum(ustar))
		P.plot(t, ke / maximum(ke))
		P.plot(t, c * 10)
		P.plot(t, tke / maximum(tke))
		P.xlim(0, t[end])
		P.ylim(0,1.1)
	end
end

# ╔═╡ c12e6133-e6df-4371-89b4-3104dedaa114
let fn = cfg.jlpath * "/.cache/output.log"
	isfile(fn) && examine_jl_log(fn)
end

# ╔═╡ 3780b8cf-1b3d-479b-bfd6-b9a1f516b66a
plt() do P

	subset = 41:1:79
	pf = load_profiles_ft(cfg.ftpath, subset, dims = cfg.dims)
	pj = load_profiles_jl(cfg.jlpath, subset, dims = cfg.dims)

	P.plot(pf.diss, pj.zc, ".-k", lw=1)
	P.plot(pj.diss, pj.zc, ".-C1", lw=1)
	P.plot(pf.prod, pj.zc, ".-k", lw=1)
	P.plot(pj.prod, pj.zc, ".-C0", lw=1)
	#P.plot(pf.prod ./ pf.diss, pj.zc, ".-k", lw=1)
	#P.plot(pj.prod ./ pj.diss, pj.zc, ".-C0", lw=1)
	#P.gca().set_xscale(:log)
	P.gca().set_yscale(:log)
end

# ╔═╡ 7050f5b7-d22e-4f38-9a08-e1b9bd50b6e2
plt() do P

	plot = :spectra
	subset_jl = 30:48
	subset_ft = subset_jl

	p = load_profiles_jl(cfg.jlpath, subset_jl, dims = (64, 64, 64))
	p2 = load_profiles_ft(cfg.ftpath, subset_ft, dims = (64, 64, 64))
	
	if plot == :spectra
		for i=1:2:size(p.Ek1, 2)
			P.plot(p.k1 .* p.zc[i], p.Ek1[:,i] ./ p.zc[i], "-k", lw=0.2)
			#P.plot(p.k1, p.Ek1[:,i] .* p.k1 ./ p.zc[i], "-k", lw=0.2)
		end
		P.plot(10.0.^(2:5), 1e0*(10.0.^(2:5)).^(-1), "--k")
		P.plot(10.0.^(4:6), 1e3*(10.0.^(4:6)).^(-5/3), "--k")
	else
		P.contour(p.k2, p.zc, p.Ek2' .* reshape(p.k2, 1, :), 0:0.6:2, colors="C1")
		P.contour(p.k2, p.zc, p2.Ek2' .* reshape(p.k2, 1, :), 0:0.6:2, colors="k")
		#P.contour(p.k1, p.zc, p.Ek1' .* reshape(p.k1, 1, :), 0:0.2:2, colors="C1")
		#P.contour(p.k1, p.zc, p2.Ek1' .* reshape(p.k1, 1, :), 0:0.2:2, colors="k")
	end
	P.gca().set_xscale(:log)
	P.gca().set_yscale(:log)
end

# ╔═╡ Cell order:
# ╟─668e802b-eb57-474c-98cb-86fa661668d9
# ╠═67a565f2-45ad-466a-ba1e-e5fb3e53560e
# ╟─60a0bfb9-4aab-4831-a07e-a98a4f5ac20c
# ╟─76de4f32-bbfe-449c-9a23-b007fcb005e9
# ╟─fcea0662-c01d-4a6c-a7cd-7ba3dcc02d7e
# ╟─c12e6133-e6df-4371-89b4-3104dedaa114
# ╟─926c81c7-5aef-449b-ad9f-4bf4b3e00064
# ╟─6fc8046f-fb19-46a0-be95-9e70fba8f11b
# ╟─1410f160-e7fc-4961-927c-d2de6c2ee835
# ╟─78e8a13a-4a6f-4153-939f-26be6a98ec37
# ╟─4c41c8ae-e41e-44a0-a54f-df0363d47e70
# ╟─6b769ef4-8d5e-4638-8ebf-e018d6b3ce4f
# ╟─3359d7ed-c99e-4a1e-9575-67cd327ac763
# ╟─025e3866-3f96-41fc-822e-0fc78bc9eafb
# ╟─294cbe47-4d5c-48e2-8599-77ab84931824
# ╟─12068898-667d-4947-a48f-b489317a5cab
# ╟─53e4fe37-3f17-4244-bd7e-56637ec4dffa
# ╟─e67c7981-3bc3-4bc8-9400-97f99aa1eef9
# ╟─02ea41ac-a1d2-4d98-8aa5-9a69a3f400ef
# ╟─82a546e8-75d8-4c7c-abfa-4d4d2412ae00
# ╟─8df38d6b-43b8-43f6-a507-f766d5fef53d
# ╟─a11e0947-5451-4e25-8894-0e0be324ef8d
# ╟─3780b8cf-1b3d-479b-bfd6-b9a1f516b66a
# ╟─7050f5b7-d22e-4f38-9a08-e1b9bd50b6e2
# ╟─c3a61806-dd18-4ea6-b267-a5bf2ee07ad0
# ╠═c4c06d18-ce1d-45f2-9376-a6309e3e4f5a
# ╠═814af8ef-64ba-4a46-ab47-9ef4d17909c3
# ╠═6633423f-8961-43fe-8329-33192b1b142d
# ╟─9ccdea97-69ec-4992-9f1c-613cddbcefe2
