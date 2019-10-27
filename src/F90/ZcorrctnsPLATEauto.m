# given list of signal free time windows in a day, find horizontal to vertical
# and dpg to vertical displacement transfer functions
#  specific for PLATE experiment data
echo on
# enter day number for label for transfer
# function output, like 330
echo off

setbb day1 $day1

setbb list1 'OBS1 OBS16 OBS11 OBS7'
setbb list2 'OBS4 OBS15'

do STA list %list1%

 setbb bh1fn *$STA$.BH1.SAC
 setbb bh2fn *$STA$.BH2.SAC
 setbb bhdfn *$STA$.BHD.SAC
 setbb bhzfn *$STA$.BHZ.SAC

 setbb outhtrans $STA$day%day1%H2Ztransfer
 setbb outptrans $STA$day%day1%P2Ztransfer
 setbb outhcoeff $STA$day%day1%H2Zcoeff
 setbb outpcoeff $STA$day%day1%P2Zcoeff
 setbb outrecord %bhzfn%.tpzw


# enter polezero file for converting to displacement and removing 
# instrument response
 setbb pzfile ../responses/obs1.zp


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
enddo

do STA list %list2%

 setbb bh1fn *$STA$.BH1.SAC
 setbb bh2fn *$STA$.BH2.SAC
 setbb bhdfn *$STA$.BHD.SAC
 setbb bhzfn *$STA$.BHZ.SAC

 setbb outhtrans $STA$day%day1%H2Ztransfer
 setbb outptrans $STA$day%day1%P2Ztransfer
 setbb outhcoeff $STA$day%day1%H2Zcoeff
 setbb outpcoeff $STA$day%day1%P2Zcoeff
 setbb outrecord %bhzfn%.tpzw


# enter polezero file for converting to displacement and removing 
# instrument response
 setbb pzfile ../responses/obs4.zp


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
enddo

echo off
cut off
