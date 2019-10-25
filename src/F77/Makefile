#!/bin/bash

#FF=gfortran
FF=ifort

SACLIB='/Users/yruan/src/sac-64/lib/libsacio.a'

findtransferfnsAuto:
	$(FF) -o findtransferfnsAuto findtransferfnsAuto.f -lm $(SACLIB)
findP2ZtransferfnsAuto:
	$(FF) -o findP2ZtransferfnsAuto findP2ZtransferfnsAuto.f -lm $(SACLIB)
remveTiltNoiseAuto:
	$(FF) -o remveTiltNoiseAuto remveTiltNoiseAuto.f -lm $(SACLIB)
remveWaterNoiseAuto:
	$(FF) -o remveWaterNoiseAuto remveWaterNoiseAuto.f -lm $(SACLIB)
rotate2tiltdir:
	$(FF) -o rotate2tiltdir rotate2tiltdir.f -lm $(SACLIB)
clean:
	\rm -f findP2ZtransferfnsAuto findtransferfnsAuto remveTiltNoiseAuto remveWaterNoiseAuto rotate2tiltdir
