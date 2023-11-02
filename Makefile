.DEFAULT_GOAL := help
.PHONY: manuscript archive dns-julia les-julia les-fortran perf
.INTERMEDIATE: git-commit-ids

help:
	@echo 'Run `make <goal>` with one of these goals:'
	@echo '- `manuscript`: Build PDF manuscript with LaTeX'
	@echo '- `archive`: Create ZIP archive of project files'
	@echo '- `dns-julia`: Run DNS with new code for validation'
	@echo '- `les-julia`: Run LES with new code for validation'
	@echo '- `les-fortran`: Run LES with Fortran code for validation'
	@echo '- `perf`: Launch simulations to gather performance data'

manuscript:
	cd paper && make default

archive: latest.zip

dns-julia:
	cd data/les-fortran && make

les-julia:
	cd data/les-fortran && make

les-fortran:
	cd data/les-fortran && make

perf:
	cd data/performance && make


# --- ARCHIVE ---------------------------------- #

git-commit-ids:
	git rev-parse HEAD > $@
	git submodule status >> $@

latest.zip: git-commit-ids
	git archive -o $@ $(foreach file,$^,--prefix=$(dir $(file)) --add-file=$(file)) --prefix= HEAD
