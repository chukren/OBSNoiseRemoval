# given list of signal free time windows in a day, find horizontal to vertical
# and dpg to vertical displacement transfer functions
#  specific for PLATE experiment data
echo on
# enter files names in order h1, h2, dpg, vert.  Then enter label for transfer
# function output, like OBS1day330
echo off
setbb bh1fn $bh1fn
setbb bh2fn $bh2fn
setbb bhdfn $bhdfn
setbb bhzfn $bhzfn
setbb outfn $outfn

setbb outhtrans %outfn%H2Ztransfer
setbb outptrans %outfn%P2Ztransfer
setbb outhcoeff %outfn%H2Zcoeff
setbb outpcoeff %outfn%P2Zcoeff
setbb outrecord %bhzfn%.tpzw

echo on
# enter polezero file for converting to displacement and removing 
# instrument response
echo off
setbb pzfile $pzfile

r %bh1fn %bh2fn %bhdfn %bhzfn
synch
w temp1 temp2 tempd tempz
cut 0 84500
r temp1 temp2 tempd tempz
rmean
dec 5
dec 2
rmean
w temp1 temp2 tempd tempz

# convert to displacement, removing instrument response 
cut off
r temp1 temp2 tempz
transfer from polezero subtype %pzfile to none freq .003 .005 1. 2. prew 2
int rect
w temp1 temp2 tempz
sc ./../findtransferfnsAuto < findtransferfnsAuto.inp .
cp tempH2Ztransfer %outhtrans
cp tempH2Zcoeff %outhcoeff

r temph 
#predict vertical noise from rotated horizontal record and transferfunction
fft amph
writesp amph temph
sc ./../remveTiltNoiseAuto
readsp amph tempf
ifft
rmean
w tempPred

# subtract from vertical
r tempz
subf tempPred
w tempztiltc

sc ./../findP2ZtransferfnsAuto < findtransferfnsAuto.inp .
cp tempP2Ztransfer %outptrans
cp tempP2Zcoeff %outpcoeff
cut off
r tempd
rmean
taper
# tapering because very long-period tidal components
# cause discontinuities at ends of record
fft amph
writesp amph temp
sc ./../remveWaterNoiseAuto
readsp amph tempf
ifft
rmean
w tempPred
# enter vertical filename
r tempztiltc
rmean
taper
subf tempPred
w %outrecord

r tempz tempztiltc %outrecord
rmean
bp co .01 .02 n 4 p 2
p1

echo off

