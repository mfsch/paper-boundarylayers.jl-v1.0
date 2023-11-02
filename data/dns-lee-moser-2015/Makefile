.DEFAULT_GOAL := all
.PHONY: all

baseurl := http://turbulence.oden.utexas.edu/channel2015/data

all: \
    LM_Channel_0180_mean_prof.dat \
    LM_Channel_0180_vel_fluc_prof.dat \
    LM_Channel_0180_RSTE_k_prof.dat \
    LM_Channel_0180_1d_energy_spectra.h5

LM_Channel_%:
	wget $(baseurl)/$@
