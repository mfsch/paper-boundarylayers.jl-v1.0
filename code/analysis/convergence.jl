### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 68359980-d40a-4e7a-b38a-acc31984385a
begin
	LOAD_PATH = ["@", "@stdlib"]
	using Pkg
	Pkg.activate(".")
	import PyPlot, HDF5, Random
	using BoundaryLayerDynamics
	include("../BoundaryLayerDynamics.jl/test/laminar_flow_problems.jl")
	global_maximum = maximum # required in laminar_flow_problems.jl
end;

# ╔═╡ 411f9a8c-bdf6-406c-bfba-57b869160007
const cfg = (
	# set here or in environment to overwrite HDF5 files
	writeh5 = lowercase(get(ENV, "WRITE", "false")) in ("true", "yes", "1"),
	pow2max = 5, # 7 for figures, takes about 6 minutes
	seed = 12298,
);

# ╔═╡ 0287875c-5b58-4928-b1f7-19cd165e6394
function runtests(Nx, Nts; T=Float64, seed=rand(UInt))

	Random.seed!(seed)

	# parameters for Poiseuille & Couette flow
	ν  = one(T) / 4 * 3 + rand(T) / 2
	ex = (θ = rand(T) * 2 * π; (cos(θ), sin(θ)))
	δ  = one(T) / 4 * 3 + rand(T) / 2
	uτ = one(T) / 4 * 3 + rand(T) / 2
 	t  = (δ^2 / ν) / 6
	Nt = Nt_viscous(T, Nx[end], t=t, δ=δ, ν=ν, Cmax=1/4)

	εp = poiseuille_error.(T, 3, Nx, Nt, t=t, ν=ν, δ=δ, uτ=uτ, dir=ex, η=nothing, method=SSPRK33())
	εc = couette_error.(T, 3, Nx, Nt, t=t, ν=ν, δ=δ, uτ=uτ, dir=ex, η=nothing, method=SSPRK33())

	# parameters for Taylor-Green vortex
    α  = one(T) / 4 * 3 + rand(T) / 2
    β  = one(T) / 4 * 3 + rand(T) / 2
    γ  = one(T) / 4 * 3 + rand(T) / 2
    U  = one(T) / 4 * 3 + rand(T) / 2
    t  = 1 / ((α^2 + β^2 + γ^2) * ν)
    Nt = Nt_viscous(T, Nx[end], δ=T(π)/(2*β), ν=ν, t=t, Cmax=1/4)

    εtgv = taylor_green_vortex_error.(T, 3, Nx, Nt, t=t, ν=ν, α=α, β=β, γ=γ, U=U, method=SSPRK33())
    εtgh = [taylor_green_vortex_error.(T, 3, 3, Nts,  t=t, ν=ν, α=α, β=β, γ=γ, U=U, W=zero(T), method=m) for m=(Euler(), AB2(), SSPRK22(), SSPRK33())]

	εp, εc, εtgv, εtgh
end

# ╔═╡ bac705f6-e4f2-4d5f-a8e7-2b3429652e31
begin
	Nx = [2^i for i=2:cfg.pow2max]
	Nt = [2^i for i=2:cfg.pow2max]
	εp, εc, εtgv, εtgh = runtests(Nx, Nt, seed=cfg.seed)
end

# ╔═╡ 4aca4118-14c4-49a0-a5ce-d748c7298c25
if cfg.writeh5
	HDF5.h5open("../../data/errors-laminar.h5", "w") do h5
		h5["nx"] = Nx
		h5["nt"] = Nt
		h5["err-poiseuille"] = εp
		h5["err-couette"] = εc
		h5["err-taylorgreen"] = εtgv
		h5["err-tg-euler"] = εtgh[1]
		h5["err-tg-ab2"] = εtgh[2]
		h5["err-tg-ssprk22"] = εtgh[3]
		h5["err-tg-ssprk33"] = εtgh[4]
	end
end

# ╔═╡ fcc49851-122e-480d-bdd6-1a4eff1d2843
function plt(f::Function)
    PyPlot.close("all") # clean up if the last call has crashed
    f(PyPlot)
    figs = [PyPlot.figure(num) for num in PyPlot.get_fignums()] # force order
    PyPlot.close("all")
    length(figs) == 1 ? figs[1] : figs
end

# ╔═╡ b0e90eda-d293-4cfb-a10f-1e942074f2b8
plt() do p
	fig, ax = p.subplots(1, 2, figsize=(12, 6), sharey=true)
	ax[1].loglog(Nx, εp, ".-", basex=2, basey=2)
	ax[1].loglog(Nx, εc, ".-", basex=2, basey=2)
	ax[1].loglog(Nx, εtgv, ".-", basex=2, basey=2)
	ax[2].loglog(Nt, εtgh[1], ".-", basex=2, basey=2)
	ax[2].loglog(Nt, εtgh[2], ".-", basex=2, basey=2)
	ax[2].loglog(Nt, εtgh[3], ".-", basex=2, basey=2)
	ax[2].loglog(Nt, εtgh[4], ".-", basex=2, basey=2)
	ax[1].grid(true)
	ax[2].grid(true)
	ax[2].set_yticks([2.0^i for i=-(5+3*cfg.pow2max):0])
end

# ╔═╡ Cell order:
# ╠═411f9a8c-bdf6-406c-bfba-57b869160007
# ╠═68359980-d40a-4e7a-b38a-acc31984385a
# ╠═0287875c-5b58-4928-b1f7-19cd165e6394
# ╠═bac705f6-e4f2-4d5f-a8e7-2b3429652e31
# ╠═4aca4118-14c4-49a0-a5ce-d748c7298c25
# ╠═b0e90eda-d293-4cfb-a10f-1e942074f2b8
# ╠═fcc49851-122e-480d-bdd6-1a4eff1d2843
