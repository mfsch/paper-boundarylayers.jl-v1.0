### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 770498c7-7321-4b46-9bce-9f76ba71387f
begin
	LOAD_PATH = ["@", "@stdlib"]
	using Pkg
	Pkg.activate(".")
	import PyPlot, HDF5, DelimitedFiles, JSON3, BenchmarkTools, FFTW, Statistics
end

# ╔═╡ 5f67622d-28f6-487d-ab3c-dbe5687f08bc
md"# Figures for Julia Code Paper"

# ╔═╡ 8a86c512-6a37-4849-bd5a-2f5741dbeca9
const cfg = (
	# set here or in environment to overwrite HDF5 files
	writefigs = lowercase(get(ENV, "WRITE", "false")) in ("true", "yes", "1"),
	writeh5 = lowercase(get(ENV, "WRITEH5", "false")) in ("true", "yes", "1"),
	codedir = abspath("."),
	datadir = abspath("../../data"),
	figsdir = abspath("../../figures"),
);

# ╔═╡ d0bfa853-cf87-4f1b-b49f-b14b3f1b9ed3
md"## Validation Laminar"

# ╔═╡ 9ddd165f-6bf6-4aa3-9fc2-6c6adce3a0ea
md"## Validation DNS"

# ╔═╡ 5dd0205a-5903-4141-b6cf-99649d6f5379
md"## Validation LES"

# ╔═╡ f0359c18-7687-4aba-9634-57115dcaffed
md"## Performance & Scaling"

# ╔═╡ ccba1953-37db-4cfa-a072-99403cae3cc3
md"## Write performance data to HDF5 file"

# ╔═╡ 9bc89a83-8001-4af3-8d35-97543bcb703e
function load_benchmark_suites(dir; arch = "icx", nh = 256)
	dir = joinpath(dir, "benchmarks")
	bm = BenchmarkTools.BenchmarkGroup()
	
	function get!(bg, k)
		if !haskey(bg, k)
			bg[k] = BenchmarkTools.BenchmarkGroup()
		end
		bg[k]
	end
	
	for fn in filter(startswith("suite-$arch"), readdir(dir))
		np, nh, nv = parse.(Int, split(fn[1:end-5], '-')[3:end])
		nvpp = Int(nv/np)
		bg = get!(get!(bm, nh), nvpp)
		bg[np] = BenchmarkTools.load("$dir/$fn")[]
	end
	bm[nh]
end

# ╔═╡ dceb48af-f0fc-4f0a-ac90-e661c0a67192
function load_integration_benchmarks(dir, arch = "icx")

	jldir = joinpath(dir, "benchmarks")
	jlfns = filter(startswith("integration-$arch"), readdir(jldir))
	bm_jl = map(filter(endswith("json"), jlfns)) do fn
		np, nh, nv = parse.(Int, split(fn, '-')[3:5])
		sim = fn[end-7:end-5]
		data = JSON3.read(read(joinpath(jldir, fn), String))
		("julia-$sim", np) => Statistics.median(diff(data["wallTime"]) ./ data["frequency"])
	end
	
	ftdir = joinpath(dir, "logs")
	ftfns = filter(startswith("integration-$arch"), readdir(ftdir))
	bm_ft = map(filter(endswith("log"), ftfns)) do fn
		np = parse(Int, split(fn, '-')[3])
		code = split(fn, '-')[4]
		sim = fn[end-6:end-4]
		lines = readlines(joinpath(ftdir, fn))
		steps = map(filter(contains("SIM"), lines)) do line
			parse(Int, split(line, ' ')[end])
		end
		times = map(filter(contains("time since beginning"), lines)) do line
			parse(Float64, split(strip(line), ' ')[end])
		end
		("$code-$sim", np) => Statistics.median(diff(times) ./ diff(steps))
	end

	np = sort(unique(last.(vcat(first.(bm_ft), first.(bm_jl)))))
	times = Dict(bm_jl..., bm_ft...)
	(processes = np,
		fortran_dns = [times[("fortran-dns", np)] for np=np],
		fortran_les = [times[("fortran-les", np)] for np=np],
		julia_dns = [times[("julia-dns", np)] for np=np],
		julia_les = [times[("julia-les", np)] for np=np],
	)
end

# ╔═╡ e713ebe2-fe50-4a73-98cd-dc20d6baea86
benchmarks = cfg.writeh5 ? HDF5.h5open(cfg.datadir * "/performance.h5", "w") do h5
	bms = load_benchmark_suites(cfg.datadir * "/performance") |>
		BenchmarkTools.median |> BenchmarkTools.time
	bmi = load_integration_benchmarks(cfg.datadir * "/performance")
	h5s = HDF5.create_group(h5, "suites")
	h5i = HDF5.create_group(h5, "integration")

	for term in ("continuity", "diffusion", "sgs_stress", "advection")
		g1 = HDF5.create_group(h5s, term)
		for nvpp in keys(bms)
			g2 = HDF5.create_group(g1, string(nvpp))
			for np in keys(bms[nvpp])
				g2[string(np)] = bms[nvpp][np][term]
			end
		end
	end
	for k in keys(bmi)
		h5i[string(k)] = collect(bmi[k])
	end
	"wrote performance data"
end : "skipped writing HDF5 file"

# ╔═╡ 85bfa2a1-f6f1-4a85-bd60-c2c88cd6af8e
md"## Helper Functions"

# ╔═╡ 6bec8977-07bf-47e6-86af-c7e49309a1bd
h5read(fn, group = nothing) = HDF5.h5open(fn) do h5
	h5 = isnothing(group) ? h5 : h5[group]
	NamedTuple(Symbol(k) => read(h5, k) for k in keys(h5))
end

# ╔═╡ f1d6e064-b036-11eb-1fad-730b31b597a7
function plt(f::Function)
    PyPlot.close("all") # clean up if the last call has crashed
	
	# font setup
	rcp = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
	rcp["font.family"] = "monospace"
	rcp["font.monospace"] = ["Source Code Pro"]
	rcp["text.usetex"] = false

	f(PyPlot)
    figs = [PyPlot.figure(num) for num in PyPlot.get_fignums()] # force order
    PyPlot.close("all")
    length(figs) == 1 ? figs[1] : figs
end

# ╔═╡ 340ed995-8a10-48a5-8562-a41814a0db01
plt() do p
	HDF5.h5open(joinpath(cfg.datadir, "errors-laminar.h5")) do h5

		imin = 2

		nx = read(h5, "nx")[imin:end]
		εp = read(h5, "err-poiseuille")[imin:end]
		εc = read(h5, "err-couette")[imin:end]
		εtg = read(h5, "err-taylorgreen")[imin:end]

		nt = read(h5, "nt")[imin:end]
		εeul = read(h5, "err-tg-euler")[imin:end]
		εab2 = read(h5, "err-tg-ab2")[imin:end]
		εrk2 = read(h5, "err-tg-ssprk22")[imin:end]
		εrk3 = read(h5, "err-tg-ssprk33")[imin:end]

		fig, ax = p.subplots(1, 4, figsize=(12,4), sharey=true)
		xr = 2.0 .^ [3, 7]
		yr = 2.0 .^ [-29, -4]

		# plot grid
		gdpad = 1.2
		gdpows = -24:2:24
		for iax=1:4
			for i=gdpows
				xs = [2^0, 2^8]
				gdcol = "#eeeeee"
				ord = [2, 1, 2, 3][iax]
				ax[iax].plot(xs, [2.0^i, 2.0^(i-ord*diff(log.(2,xs))[])],
					color=gdcol, lw=1)
			end
		end

		ax[1].plot(nx, εtg, ".:k", ms=10, label="Taylor-Green vortex")
		ax[1].plot(nx, εp, ".--k", ms=10, label="Poiseuille flow")
		ax[1].plot(nx, εc, ".-k", ms=10, label="Couette flow")

		ax[2].plot(nt, εeul, ".-k", ms=10, label="forward Euler")
		ax[3].plot(nt, εab2, ".--k", ms=10, label="Adams-Bashforth 2")
		ax[3].plot(nt, εrk2, ".-k", ms=10, label="SSPRK (2,2)")
		ax[4].plot(nt, εrk3, ".-k", ms=10, label="SSPRK (3,3)")

		# configure plot details
		[ax.legend(frameon=false, loc=3) for ax=ax]
		ax[1].set_xlabel(PyPlot.L"N_3")
		ax[2].set_xlabel(PyPlot.L"N_t")
		ax[3].set_xlabel(PyPlot.L"N_t")
		ax[4].set_xlabel(PyPlot.L"N_t")
		ax[1].set_ylabel("error")

		for ax=ax
			ax.set_xscale("log", base=2)
			ax.set_yscale("log", base=2)
			ax.set_xlim(xr[1]/gdpad, xr[2]*gdpad)
		end
		ax[1].set_ylim(yr[1]/gdpad, yr[2])
		[ax.set_xticks(2.0.^(3:1:7)) for ax=ax]
		ax[1].set_yticks(2.0.^(-28:4:-4))

		p.subplots_adjust(wspace=0)
		cfg.writefigs && p.savefig("$(cfg.figsdir)/validation-laminar.pdf", bbox_inches="tight", metadata=Dict("CreationDate" => nothing))
	end
end

# ╔═╡ bfbd6dcb-6a87-4b6b-b67b-c2847c84ad45
plt() do p
	fig, ax = p.subplots(1, 4, figsize=(12,4), sharey=true)

	ms = 10
	leg = ax[1].legend([p.plt.Line2D([0], [0],
				color=("k", "C1")[i],
				marker=("", ".")[i],
				ls=("-", "")[i], ms=ms) for i=1:2],
		("Lee & Moser (2015)", "new code"), frameon=false, loc=2, numpoints=2)

	# set up y-axis
	ax[1].set_yscale("log", base=10)
	ax[1].set_ylim(5e-1, 185)
	ax[1].set_ylabel(PyPlot.L"x_3^+")
	ax[1].set_yticks([1, 10, 100, 182])
	ax[1].yaxis.set_major_formatter(p.matplotlib.ticker.ScalarFormatter())
	
	#lm, jl = lmdns, jldns
	lm = h5read(cfg.datadir * "/dns-profiles.h5", "lee-moser-2015")
	jl = h5read(cfg.datadir * "/dns-profiles.h5", "abl-jl")
	
	# plot velocities
	ax[1].set_xlim(0, 19)
	ax[1].plot(lm.u1, lm.x3p, "-k")
	ax[1].plot(jl.u1, jl.x3p, ".C1", ms=ms)
	ax[1].legend([p.plt.Line2D([0], [0], color="k"), ],
		("streamwise velocity", ), frameon=false, loc=4)
	ax[1].set_xlabel(PyPlot.L"u_1^+")
	
	# plot momentum balance
	ax[2].set_xlim(-0.05, 1.05)
	ax[2].set_xticks((0, 0.5, 1))
	ax[2].set_xticklabels(("0", "0.5", "1"))
	ax[2].plot(lm.τdiff, lm.x3p, ":k")
	ax[2].plot(jl.τdiff, jl.x3p, ".C3", ms=ms)
	ax[2].plot(lm.τturb, lm.x3p, "-k")
	ax[2].plot(jl.τturb, jl.x3p, ".C1", ms=ms)
	ax[2].legend([p.plt.Line2D([0], [0], color="k",
		ls=("-","--",":")[i]) for i=1:2:3],
		("advective transp.", "total transp.", "diffusive transp.")[1:2:3],
		frameon=false, loc=8)
	ax[2].set_xlabel("\$\\tau_{13}^+\$")
	
	# plot energy budget terms
	ax[3].set_xticks((0, 0.1, 0.2))
	ax[3].set_xticklabels(("0", "0.1", "0.2"))
	ax[3].plot(lm.diss, lm.x3p, ":k")
	ax[3].plot(jl.diss, jl.x3p, ".C3", ms=ms)
	ax[3].plot(lm.prod, lm.x3p, "-k")
	ax[3].plot(jl.prod, jl.x3p, ".C1", ms=ms)
	ax[3].legend([p.plt.Line2D([0], [0], color="k", ls=("-",":")[i]) for i=1:2],
		("production", "dissipation"), frameon=false, loc=1)
	ax[3].set_xlabel("\$\\mathcal{P}^+\\!,\\; \\varepsilon^+\$")
	
	# plot energy spectra
	contours = 1.2 * (0.75:1:10) ./ 182
	contours = (1:1.25:100) ./ 182
	cratio = 1/5 # premultiplied spectra have different magnitudes
	cms = ms/1.75 # marker for contour plot lines
	clines = [(0,(0.01,1.8))]
	ax[4].set_xscale("log", base=10)
	ax[4].set_xlim((0.5, 20) ./ 182) # kmin = 0.5 = 1 * 2π / (4π)
	ax[4].contour(lm.k2p, lm.x3p, lm.E2 .* reshape(lm.k2p, 1, :), 
		contours, colors="k", linestyles=":")
	ax[4].contour(lm.k1p, lm.x3p, lm.E1 .* reshape(lm.k1p, 1, :),
		cratio*contours, colors="k")
	cs = ax[4].contour(jl.k2p, jl.x3p, jl.E2 .* reshape(jl.k2p, 1, :),
		contours, colors="C3", linestyles=clines, linewidths=cms)
	[ls.set_capstyle("round") for ls in cs.collections]
	cs = ax[4].contour(jl.k1p, jl.x3p, jl.E1 .* reshape(jl.k1p, 1, :),
		cratio*contours, colors="C1", linestyles=clines, linewidths=cms)
	[ls.set_capstyle("round") for ls in cs.collections]
	ax[4].legend([p.plt.Line2D([0], [0], color="k", ls=("-",":")[i]) for i=1:2],
		(PyPlot.L"k_1^+ E^+_{ii}\!(\,k_1^+)", PyPlot.L"k_2^+ E^+_{ii}\!(\,k_2^+)"),
		frameon=false, loc=4, ncol=1)
	ax[4].set_xlabel(PyPlot.L"k_1^+\!,\; k_2^+")
	
	ax[1].add_artist(leg)
	ax[4].yaxis.set_tick_params(which="both", right=true)
	p.subplots_adjust(wspace=0)
	cfg.writefigs && p.savefig("$(cfg.figsdir)/validation-dns.pdf",
		bbox_inches="tight", metadata=Dict("CreationDate" => nothing))
end

# ╔═╡ 05c20f4b-50ad-412c-89ad-bfc132fa99f0
plt() do p
	fig, ax = p.subplots(1, 4, figsize=(12,4), sharey=true)

	κ = 0.40
	z0 = 1e-4
	ms = 10
	l1, l2 = (2π, 4π/3)
	leg = ax[1].legend([p.plt.Line2D([0], [0], 
				color=("k", "C1")[i],
				marker=("", ".")[i],
				ls=("-", "")[i], ms=ms) for i=1:2],
		("Fortran code", "new code"), frameon=false, loc=2, numpoints=2)
	yrange = (50, 1e4) # first dot at ½Δz/z0=78.125

	# set up y-axis
	ax[1].set_yscale("log", base=10)
	ax[1].set_ylim(26, 1e4)
	ax[1].set_ylabel(PyPlot.L"x_3 z_0^{-1}")

	ft = h5read("$(cfg.datadir)/les-profiles.h5", "albertson-parlange")
	jl = h5read("$(cfg.datadir)/les-profiles.h5", "abl-jl")
	
	# plot velocities
	ax[1].set_xlim(10, 27.5)
	ax[1].plot(log.(yrange)./κ, yrange, "--k")
	ax[1].plot(ft.u, jl.zc, "+-k", lw=1, ms=ms)
	ax[1].plot(jl.u, jl.zc, ".C1", ms=ms)
	ax[1].legend([p.plt.Line2D([0], [0], color="k", ls=("","--")[i], marker=("+","")[i], ms=ms) for i=1:2],
		("streamwise velocity","log-law with \$\\kappa = 0.4\$"), frameon=false, loc=4, numpoints=2)
	ax[1].set_xlabel(PyPlot.L"u_1 u_\tau^{-1}")
	
	
	# plot momentum balance
	ax[2].set_xlim(-0.05, 1.05)
	ax[2].set_xticks((0, 0.5, 1))
	ax[2].set_xticklabels(("0", "0.5", "1"))
	interp(u) = (u[1:end-1] .+ u[2:end]) / 2
	ax[2].plot(-interp(ft.sgs13), jl.zc, "x-k", lw=1, ms=ms/sqrt(2))
	ax[2].plot(-interp(jl.sgs13), jl.zc, ".C3", ms=ms)
	ax[2].plot(-ft.uw, jl.zc, "+-k", lw=1, ms=ms)
	ax[2].plot(-jl.uw, jl.zc, ".C1", ms=ms)
	ax[2].legend([p.plt.Line2D([0], [0], color="k", ls="", marker=("+","x")[i],
		ms=(ms, ms/sqrt(2))[i]) for i=1:2],
		("resolved transp.", "subgrid-scale transp."),
		frameon=false, loc=3, numpoints=2)
	ax[2].set_xlabel(PyPlot.L"\tau_{13} u_\tau^{-2}")
	
	# plot energy balance
	ax[3].set_xticks((0, 0.005))
	ax[3].set_xticklabels(("0", "0.005"))
	ax[3].plot(ft.diss, jl.zc, "x-k", lw=1, ms=ms/sqrt(2))
	ax[3].plot(jl.diss, jl.zc, ".C3", ms=ms)
	ax[3].plot(ft.prod, jl.zc, "+-k", lw=1, ms=ms)
	ax[3].plot(jl.prod, jl.zc, ".C1", ms=ms)
	ax[3].legend([p.plt.Line2D([0], [0], color="k", ls="", marker=("+","x")[i], ms=(ms, ms/sqrt(2))[i]) for i=1:2],
		("production", "dissipation"), frameon=false, loc=1, numpoints=2)
	ax[3].set_xlabel(PyPlot.L"\mathcal{P}\, z_0 u_\tau^{-3} \!,\; \varepsilon\,z_0 u_\tau^{-3}")
	
	# plot energy spectra
	contours = (0:0.03:10, 0:0.125:10)
	ax[4].set_xscale("log", base=10)
	#ax[4].set_xticks((0.5, 1, 10))
	#ax[4].set_xticklabels((" ½ ", "1", "10")) # ½ doesn’t render without spaces
	#ax[4].xaxis.set_major_formatter(p.matplotlib.ticker.ScalarFormatter())
	ax[4].set_xlim(2*π/l1*z0, 2*π*31/l1*z0)
	ax[4].contour(2*π*jl.k1/l1*z0, jl.zc, ft.Ek1' .* reshape(jl.k1/l1, 1, :),
		contours[1], colors="k")
	cs = ax[4].contour(2*π*jl.k1/l1*z0, jl.zc, jl.Ek1' .* reshape(jl.k2/l1, 1, :),
		contours[1], colors="C1", linestyles=[(0,(0.01,1.5))], linewidths=ms/2)
	[ls.set_capstyle("round") for ls in cs.collections]
	if false # include E(k2)
		ax[4].contour(2*π*jl.k2/l2*z0, jl.zc, ft.Ek2' .* reshape(jl.k2/l2, 1, :),
			contours[2], colors="k", linestyles=":")
		cs = ax[4].contour(2*π*jl.k2/l2*z0, jl.zc, jl.Ek2' .* reshape(jl.k2/l2, 1, :),
			contours[2], colors="C3", linestyles=[(0,(0.01,1.5))], linewidths=ms/2)
		[ls.set_capstyle("round") for ls in cs.collections]
	end
	ax[4].legend([p.plt.Line2D([0], [0], color="k", ls=("-",":")[i]) for i=1:1],
		(PyPlot.L"k_1 E_{ii} z_0 u_\tau^{-2}", PyPlot.L"k_2 E_{ii} z_0 u_\tau^{-2}"), frameon=false, loc=8, ncol=2)
	ax[4].set_xlabel(PyPlot.L"k_1 z_0")
	
	ax[1].add_artist(leg)
	#ax[4].yaxis.set_tick_params(which="both", right=true)
	p.subplots_adjust(wspace=0)
	cfg.writefigs && p.savefig("$(cfg.figsdir)/validation-les.pdf", bbox_inches="tight", metadata=Dict("CreationDate" => nothing))
end

# ╔═╡ a1604e8a-d5b1-41b8-9e69-047f74645afc
plt() do p

	nvmain = 1280
	tref = 1e9 # 1 s = 1e9 ns
	gdpows = -16:1:4
	ms = 10
	benchmarks # to enforce correct order in Pluto notebook
	
	bms = h5read(cfg.datadir * "/performance.h5", "suites")
	nps = sort(parse.(Int, unique(vcat(collect.(keys.(values(first(bms))))...))))
	
	f, ax = p.subplots(1, 5, sharey=true, figsize=(12,4))
	[ax[i].set_title(("pressure", "advection (resolved)", "advection (SGS)",
		"diffusion", "full time step")[i], size=11, weight="semibold",
		#ha="left", x=0.01) for i=1:5]
		ha="right", va="top", x=0.97, y=0.95) for i=1:5]
	
	# plot grid with perfect scaling in background
	for ax=ax
		xs = extrema(nps)
		xs = [80, 1280]
		xs = 2 .^ ([-0.5, 0.5] .+ log.(2, xs)) # power of 2 padding
		ax.set_xlim(xs...)
		for i=gdpows
			gdcol = "#eeeeee"
			ax.plot(xs, [2.0, 2.0] .^ i, color=gdcol, lw=1)
			ax.plot(xs, [2.0^i, 2.0^(i-diff(log.(2,xs))[])], color=gdcol, lw=1)
		end
	end
	
	# weak scaling (lines only)
	nv = Set()
	for nvpp in reverse(sort(parse.(Int, collect(keys(first(bms))))))
		col = "k"
		np = sort(parse.(Int, collect(keys(first(bms)[string(nvpp)]))))
		push!(nv, (nvpp * np)...) # collect nv for strong scaling below
		t(bm) = [bm[string(nvpp)][string(np)]/tref for np=np]
		ax[1].plot(np, t(bms.continuity), ":k", c=col)
		ax[2].plot(np, t(bms.advection), ":k", c=col)
		ax[3].plot(np, t(bms.sgs_stress), ":k", c=col)
		ax[4].plot(np, t(bms.diffusion), ":k", c=col)
	end

	# strong scaling (markers and lines)
	for nv in sort(collect(nv))
		col = nv == nvmain ? "C1" : "k"
		nvpp = parse.(Int, keys(first(bms)))
		np = map(nvpp) do nvpp
			np = Int(nv / nvpp)
			haskey(first(bms)[string(nvpp)], string(np)) ? np : nothing
		end |> nps -> filter(np -> !isnothing(np), nps) |> sort
		t(bm) = [bm[string(div(nv,np))][string(np)]/tref for np=np]
		ax[1].plot(np, t(bms.continuity), ".--k", ms=ms, c=col)
		ax[2].plot(np, t(bms.advection), ".--k", ms=ms, c=col)
		ax[3].plot(np, t(bms.sgs_stress), ".--k", ms=ms, c=col)
		ax[4].plot(np, t(bms.diffusion), ".--k", ms=ms, c=col)
	end
	
	# overall performance
	mew = 2.5
	bmi = h5read(cfg.datadir * "/performance.h5", "integration")
	np = bmi.processes
	#ax[5].plot(np, bmi.oldfortran_dns, "+--", color="C2", ms=ms, mew=mew)
	#ax[5].plot(np, bmi.oldfortran_les, "x--", color="C2", ms=0.8*ms, mew=mew)
	ax[5].plot(np, bmi.fortran_dns, "+--", color="C4", ms=ms, mew=mew)
	ax[5].plot(np, bmi.fortran_les, "x--", color="C4", ms=0.8*ms, mew=mew)
	ax[5].plot(np, bmi.julia_dns, "+--", color="C1", ms=ms, mew=mew)
	ax[5].plot(np, bmi.julia_les, "x--", color="C1", ms=0.8*ms, mew=mew)
	
	# format x-axis
	f.text(.5, .05, "number of MPI processes", ha="center")
	for ax=ax
		ax.set_xscale("log", base=2)
		ax.set_xticks(nps)
		ax.xaxis.set_major_formatter(p.matplotlib.ticker.ScalarFormatter())
	end

	# format y-axis
	[ax.yaxis.set_tick_params(left=false) for ax=ax[2:end]]
	ax[end].yaxis.set_tick_params(right=true)
	ax[1].set_ylabel("time per evaluation (s)")
	ax[1].set_yscale("log", base=2)
	ax[1].set_yticks(2.0.^gdpows)
	ax[1].set_ylim(2.0^(-9), 2^0)
	
	p.tight_layout(rect=(0, 0.05, 1, 1))
	p.subplots_adjust(wspace=0)
	cfg.writefigs && p.savefig(cfg.figsdir * "/performance.pdf", bbox_inches="tight", metadata=Dict("CreationDate" => nothing))
end

# ╔═╡ Cell order:
# ╟─5f67622d-28f6-487d-ab3c-dbe5687f08bc
# ╠═8a86c512-6a37-4849-bd5a-2f5741dbeca9
# ╟─770498c7-7321-4b46-9bce-9f76ba71387f
# ╟─d0bfa853-cf87-4f1b-b49f-b14b3f1b9ed3
# ╟─340ed995-8a10-48a5-8562-a41814a0db01
# ╟─9ddd165f-6bf6-4aa3-9fc2-6c6adce3a0ea
# ╟─bfbd6dcb-6a87-4b6b-b67b-c2847c84ad45
# ╟─5dd0205a-5903-4141-b6cf-99649d6f5379
# ╟─05c20f4b-50ad-412c-89ad-bfc132fa99f0
# ╟─f0359c18-7687-4aba-9634-57115dcaffed
# ╟─a1604e8a-d5b1-41b8-9e69-047f74645afc
# ╟─ccba1953-37db-4cfa-a072-99403cae3cc3
# ╟─e713ebe2-fe50-4a73-98cd-dc20d6baea86
# ╟─9bc89a83-8001-4af3-8d35-97543bcb703e
# ╟─dceb48af-f0fc-4f0a-ac90-e661c0a67192
# ╟─85bfa2a1-f6f1-4a85-bd60-c2c88cd6af8e
# ╟─6bec8977-07bf-47e6-86af-c7e49309a1bd
# ╟─f1d6e064-b036-11eb-1fad-730b31b597a7
