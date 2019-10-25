cremveTiltNoiseAuto.f
c Use in conjunction with sac macro ZcorrctnsAuto.m to remove tilt noise from
c vertical using transfer function from rotated horizontal
c  This routine generates the predicted noise in the 
c  frequency domain to be subtracted from the original z in the time domain 
c 
c Assumes polynomial form of transfer function (quadratic) in frequency ranges

      character*80 f1a, f1p, f2a, f2p, transfer 
      character*10 OBSname
      parameter (nmax = 655360)
      dimension x1(nmax), x2(nmax), x3(nmax), x4(nmax)
      real*4 lfreq(10),ufreq(10)
      real*4 m0amp(10),m1amp(10),m2amp(10),m0ph(10),m1ph(10),m2ph(10)
      
c  f1a is amplitude spectrum of original DPG, f1p same for phase, generated 
c  in sac
c  f2a is amplitude spectrum of predicted noise, f2p same for phase
     
      f1a = 'temph.am'
      f2a = 'tempf.am'
      f1p = 'temph.ph'
      f2p = 'tempf.ph'
      transfer = 'tempH2Zcoeff'
      
      pi = 3.1415928
      twopi = 2.0*pi
      tenlg = 2.302585
      
c      open(11, file = 'dummy')
      
      call rsac1(f1a,x1,npts,beg,delf,nmax,nerr)
      call rsac1(f1p,x2,npts,beg,delf,nmax,nerr)
c      write(11,*) npts
      if (npts.gt.nmax) then
         write(*,*)  'time series too long'
	 go to 100
      end if

c  read in number of frequency bands, band limits, and polynomial transfer functions
c  of degree 2 (some coefficients may be zero) in each band
      open(10, file = transfer)
c      read(10,*) OBSname
c      write(*,*) 'correcting for water noise using coeff for ',OBSname
      read(10,*) dummyang
      read(10,*) nbands
      do i = 1, nbands
        read(10,*) dummyx
        read(10,*) lfreq(i),ufreq(i)
	read(10,*) m0amp(i),m1amp(i),m2amp(i)
	read(10,*) m0ph(i),m1ph(i),m2ph(i)
      enddo
	      
      do i = 2, npts       
        freq= delf*(i-1)
c  make no change to original fft if not within a correction band
	do j = 1, nbands
	  if ((freq.ge.lfreq(j)).and.(freq.le.ufreq(j))) then
	    x3(i) = x1(i)*(m0amp(j) + m1amp(j)*freq
     1                              + m2amp(j)*freq*freq)
            x4(i) = x2(i)-twopi*(m0ph(j) + m1ph(j)*freq
     1                              + m2ph(j)*freq*freq)
	    go to 110
	  endif
	enddo
	x3(i) = 0.0
	x4(i) = 0.0
110     continue
      enddo
      x3(1) = 0.0
      x4(1) = 0.0
	call wsac0(f2a,xdummy,x3,nerr)
	call wsac0(f2p,xdummy,x4,nerr)
100   continue
        close(unit=10)
c	close(unit=11)
      end      
