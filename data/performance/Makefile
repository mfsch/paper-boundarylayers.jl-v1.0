juliaversion = 1.9.0
tooldir := .cache
juliabin = $(tooldir)/julia-$(juliaversion)/bin/julia
juliadepot = $(tooldir)/julia-depot
juliacmd = env JULIA_DEPOT_PATH=$(tooldir)/julia-depot $(juliabin) --project=../../code/performance
modulecmd = module load impi
fortrandir = $(tooldir)/fortran-les
account = TG-ATM180022
skxsizes = 48 96 192 384 768
icxsizes = 80 160 320 640 1280

.DEFAULT: all
.PHONY: all julia update-patch clean launch-skx launch-icx \
    $(addprefix launch-skx-,$(skxsizes)) \
    $(addprefix launch-icx-,$(icxsizes))
.PRECIOUS: $(fortrandir)/bin/les-mpi-%

all: launch-skx launch-icx
launch-skx: $(addprefix launch-skx-,$(skxsizes))
launch-icx: $(addprefix launch-icx-,$(icxsizes))
clean:
	rm -rf $(tooldir)

# launch a single run on either set of nodes
$(addprefix launch-skx-,$(skxsizes)): launch-skx-%: $(fortrandir)/bin/les-mpi-% $(juliabin) $(juliadepot)
	mkdir -p $(fortrandir)/{dns,les}/{input,output}
	ln -rsf ../../code/performance/config-dns.yaml $(fortrandir)/dns/input/config.yaml
	ln -rsf ../../code/performance/config-les.yaml $(fortrandir)/les/input/config.yaml
	$(modulecmd) && env JULIACMD="$(juliacmd)" sbatch --partition=skx-normal --nodes=$$(($*/48)) --ntasks-per-node=48 --account=$(account) ../../code/performance/launch.sh
$(addprefix launch-icx-,$(icxsizes)): launch-icx-%: $(fortrandir)/bin/les-mpi-% $(juliabin) $(juliadepot)
	mkdir -p $(fortrandir)/{dns,les}/{input,output}
	ln -rsf ../../code/performance/config-dns.yaml $(fortrandir)/dns/input/config.yaml
	ln -rsf ../../code/performance/config-les.yaml $(fortrandir)/les/input/config.yaml
	$(modulecmd) && env JULIACMD="$(juliacmd)" sbatch --partition=icx-normal --nodes=$$(($*/80)) --ntasks-per-node=80 --account=$(account) ../../code/performance/launch.sh

# build a new patch based on the current code in the cache directory
update-patch:
	diff -Naur ../../code/albertson-parlange-les/src $(fortrandir)/src | tee ../../code/performance/setup.patch

# build the fortran code with the specified nproc (separate recipes for skx &
# icx nodes, since they require different vertical sizes)
$(addprefix $(fortrandir)/bin/les-mpi-,$(skxsizes)): $(fortrandir)/bin/les-mpi-%: ../../code/performance/setup-skx.patch
	mkdir -p $(fortrandir)
	rm -rf $(fortrandir)/src
	cp -r ../../code/albertson-parlange-les/src $(fortrandir)
	cat $< | sed  "/^\+.*nproc =/ s/=.*/= $*/" | patch -d$(fortrandir) -p4
	ln -rsf ../../code/albertson-parlange-les/Makefile $(fortrandir)/
	ln -rsf ../../code/albertson-parlange-les/util $(fortrandir)/
	cd $(fortrandir) && make -j BUILDDIR=.
	mv $(fortrandir)/bin/les-mpi $@
$(addprefix $(fortrandir)/bin/les-mpi-,$(icxsizes)): $(fortrandir)/bin/les-mpi-%: ../../code/performance/setup-icx.patch
	mkdir -p $(fortrandir)
	rm -rf $(fortrandir)/src
	cp -r ../../code/albertson-parlange-les/src $(fortrandir)
	cat $< | sed  "/^\+.*nproc =/ s/=.*/= $*/" | patch -d$(fortrandir) -p4
	ln -rsf ../../code/albertson-parlange-les/Makefile $(fortrandir)/
	ln -rsf ../../code/albertson-parlange-les/util $(fortrandir)/
	cd $(fortrandir) && make -j BUILDDIR=.
	mv $(fortrandir)/bin/les-mpi $@

# launch julia in interactive mode
julia: $(juliabin)
	$(modulecmd) && srun --pty -N 1 -n 1 -t 30:00 -p skx-dev -A $(account) $(juliacmd)

# download julia
$(juliabin):
	@mkdir -p $(tooldir)
	wget -O $(tooldir)/julia.tar.gz https://julialang-s3.julialang.org/bin/linux/x64/$(basename $(juliaversion))/julia-$(juliaversion)-linux-x86_64.tar.gz
	tar -xzf $(tooldir)/julia.tar.gz -C $(tooldir)

# initialize the julia package depot with all dependencies
$(juliadepot): $(juliabin)
	$(modulecmd) && srun -N 1 -n 1 -t 30:00 -p skx-dev -A $(account) $(juliacmd) -e 'import Pkg; Pkg.instantiate()'
