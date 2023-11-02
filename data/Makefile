analysis := ../code/analysis

.DEFAULT_GOAL := profiles
.PHONY: profiles

profiles: errors-laminar.h5 dns-profiles.h5 les-profiles.h5 performance.h5

errors-laminar.h5: $(analysis)/convergence.jl
	cd "$(analysis)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(analysis)" && env WRITE=true julia --project convergence.jl

dns-profiles.h5: $(analysis)/dns-budgets.jl dns-julia/snapshots dns-lee-moser-2015
	cd "$(analysis)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(analysis)" && env WRITE=true julia --project dns-budgets.jl

les-profiles.h5: $(analysis)/les-budgets.jl les-julia/snapshots
	cd "$(analysis)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(analysis)" && env WRITE=true julia --project les-budgets.jl

performance.h5: $(analysis)/figures.jl performance/benchmarks performance/logs
	cd "$(analysis)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(analysis)" && env WRITEH5=true julia --project figures.jl
