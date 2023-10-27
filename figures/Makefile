code := ../code/analysis
data := ../data
benchmarks = $(wildcard $(data)/performance/benchmarks/suite-icx-*.json)
integration = $(wildcard $(data)/performance/benchmarks/integration-icx-*.json)

.DEFAULT_GOAL := figures
.PHONY: figures

figures: fig01.pdf fig02.pdf fig03.pdf fig04.pdf

fig01.pdf: $(data)/errors-laminar.h5
fig02.pdf: $(data)/dns-profiles.h5
fig03.pdf: $(data)/les-profiles.h5
fig04.pdf: $(data)/performance.h5

$(data)/performance.h5: $(code)/figures.jl $(benchmarks) $(integration)
	cd "$(code)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(code)" && env WRITEH5=true julia --project figures.jl

fig01.pdf fig02.pdf fig03.pdf fig04.pdf &: $(code)/figures.jl
	cd "$(code)" && julia --project -e "import Pkg; Pkg.instantiate()"
	cd "$(code)" && env WRITE=true julia --project figures.jl
	mv validation-laminar.pdf fig01.pdf
	mv validation-dns.pdf fig02.pdf
	mv validation-les.pdf fig03.pdf
	mv performance.pdf fig04.pdf
