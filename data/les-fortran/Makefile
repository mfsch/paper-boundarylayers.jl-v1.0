tmpdir := .cache
codedir := ../../code/albertson-parlange-les
simdir := ../../code/validation

.DEFAULT_GOAL := instantaneous-fields/grid.yaml
.PHONY: clean update-patch

clean:
	rm -rf $(tmpdir) instantaneous-fields

# create a new patch based on the current code in the temp directory
update-patch:
	diff -Naur $(codedir)/src $(tmpdir)/src | tee $(simdir)/les-fortran-setup.patch

# build the fortran code with simulation-specific patch
$(tmpdir)/bin/les-mpi: $(simdir)/les-fortran-setup.patch
	rm -rf $(tmpdir)/{src,bin}
	mkdir -p $(@D)
	cp -r $(codedir)/src $(@D)/..
	ln -rsf $(codedir)/Makefile $(codedir)/util $(@D)/../
	patch -d$(tmpdir) -p4 < $<
	cd $(@D)/.. && make -j BUILDDIR=.

# run simulation (overwrites data)
instantaneous-fields/grid.yaml: $(simdir)/les-fortran-config.yaml $(tmpdir)/bin/les-mpi
	mkdir -p "$(tmpdir)"/{input,output}
	ln -rsf "$<" "$(tmpdir)/input/config.yaml"
	cd $(tmpdir) && env LES_TOTAL_LAUNCHES=1 LES_CURRENT_LAUNCH=1 OMP_NUM_THREADS=1 \
	    mpirun -np $$(grep 'nproc =' "src/PARAM.f90" | cut -d' ' -f5) \
	    --oversubscribe bin/les-mpi 2>&1 | tee output.log
