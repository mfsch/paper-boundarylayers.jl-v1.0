simdir := ../../code/validation
tmpdir := .cache
juliacmd := julia --project="$(simdir)"
mpiprocs := 64

.DEFAULT_GOAL := snapshots

snapshots: $(simdir)/les-julia.jl $(simdir)/Project.toml $(simdir)/Manifest.toml
	mkdir -p $(tmpdir)
	$(juliacmd) -e 'using BoundaryLayerDynamics' # make sure precompilation has run
	mpiexec -np $(mpiprocs) $(juliacmd) $< 2>&1 | tee $(tmpdir)/output.log
