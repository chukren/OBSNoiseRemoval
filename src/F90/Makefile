#!/bin/bash

#FF=gfortran
FF=ifort

SACLIB='/Users/yruan/src/sac-64/lib/libsacio.a'

findtransferfnsAuto:
	$(FF) -o findtransferfnsAuto findtransferfnsAuto.f90 -lm $(SACLIB)
clean:
	\rm -f transferp2z removeWaterNoise findtiltaxis findtiltBYtransfer
