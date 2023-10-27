analysis := ../code/analysis

.DEFAULT_GOAL := profiles
.PHONY: profiles

profiles: errors-laminar.h5 dns-profiles.h5 les-profiles.h5 performance.h5

les-profiles.h5: $(analysis)/les-budgets.jl les-julia/snapshots
	cd "$(analysis)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(analysis)" && env WRITE=true julia --project les-budgets.jl
