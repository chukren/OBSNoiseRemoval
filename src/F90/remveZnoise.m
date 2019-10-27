# removes tilt effects and pressure/compliance,.ie. water gravity wave, noise
# from vertical given previously determined transfr functions

# assumes already have synched, decimated displacement and dpg files, but 
# dpg is not corrected for response, unless used response corrected dpg in 
# finding transfr function

echo on
# enter files names in order h1, h2, dpg, vert.  Then enter label for transfr
# function input, like OBS1day330
echo off

cut off

setbb bh1fn $bh1fn
setbb bh2fn $bh2fn
setbb bhdfn $bhdfn
setbb bhzfn $bhzfn
setbb outfn $outfn

setbb outhtrans %outfn%H2Ztransfer
setbb outptrans %outfn%P2Ztransfer
setbb outhcoeff %outfn%H2Zcoeff
setbb outpcoeff %outfn%P2Zcoeff
setbb outrecord %bhzfn%.tw

r %bh1fn %bh2fn %bhdfn %bhzfn

w temp1 temp2 tempd tempz

cp %outhcoeff tempH2Zcoeff
sc ./../rotate2tiltdir < tempH2Zcoeff


r temph 
#predict vertical noise from rotated horizontal record and transferfunction
fft amph
writesp amph temph
sc ./../remveTiltNoiseAuto
readsp amph tempf
ifft
rmean
w tempPred2

# subtract from vertical
r tempz
subf tempPred2
w tempztiltc

r tempd
rmean
taper
# tapering because very long-period tidal components
# cause discontinuities at ends of record
fft amph
writesp amph temp
cp %outpcoeff tempP2Zcoeff
sc ./../remveWaterNoiseAuto
readsp amph tempf
ifft
rmean
w tempPred
r tempztiltc
rmean
taper
subf tempPred
w %outrecord


