!==========================================================================
!  Author: Donald Forsyth 
!          Department of Geological Sciences, Brown Univesity
!  Date:   April 2013, modified August 2013
!
!  Modified: chukren August 2013
!
!
!  findtransferfnsAuto.f90
!
!  This program rotates horizontal components to find maximum coherence 
!  with the vertical component (using signal-free noise samples as input).
!  Uses coherence with complex transfer function to account for phase 
!  shifts between vertical and horizontal components.  Then tilts z-axis 
!  away from that direction to minimize horizontal noise on the vertical.  
!
!  For input, pick multiple windows (~20?) files as input, avoiding times 
!  with any earthquake signals (use maximum noise periods)  
!
!  Usage:  findtransferfnsAuto  < findtransferfnsAuto.inp
!==========================================================================

program findtransferfnsAuto
    implicit none

    double precision, parameter :: pi = 3.1415926535898
    double precision, parameter :: convdeg = pi/180.
    double precision, parameter :: circ  = 6371.*pi/180.
    double precision, parameter :: twopi = pi*2.

    integer, parameter :: maxnfreq = 400
    integer, parameter :: maxpts   = 1000000
    integer, parameter :: maxevnts = 100
    integer, parameter :: maxnobs  = 800
    integer, parameter :: nparam   = 6

    real :: bh1(maxpts), bh2(maxpts), bhz(maxpts) 
    real :: freq(maxnfreq)
    real :: timebeg(maxevnts,maxevnts)
    real :: sumcrossre(maxnfreq), sumcrossim(maxnfreq)
    real :: sumppower(maxnfreq)
    real :: sumzpower(maxnfreq), tlength(maxevnts)
    real :: begh1(maxevnts), begz(maxevnts), delth1(maxevnts)
    real :: deltz(maxevnts), begh2(maxevnts), delta(maxevnts)
    real :: winh1(maxevnts,maxpts),winz(maxevnts,maxpts)
    real :: winh2(maxevnts,maxpts)
    real :: gausst(maxpts), tdataz(maxpts), tdatahr(maxpts)
    real :: tdatah1(maxpts), tdatah2(maxpts), delth2(maxevnts)
    real :: qreal(maxnfreq), qimag(maxnfreq), coher(maxnfreq)
    real :: epsilon(maxnfreq)
    real :: qamp(maxnfreq), qphse(maxnfreq), ampdev(maxnfreq)
    real :: phsedev(maxnfreq)
    real :: preal,pimag
    real :: horz(maxpts),vert(maxpts)
    real :: freqlo(10), freqhi(10)
    real :: dumfreqlo, dumfreqhi
    real :: theta,costheta,avgcoher,cohermax,maxtheta
    real :: crossim,crossre,ppower,zpower,zimag,zreal
    real :: tg,sigma,sigma2,sintheta,tmiddle
    real :: xdummy

    real :: g(maxnobs,nparam)
    real :: stddevdata, gtd(nparam), stddev(nparam)
    real :: d(maxnobs),resid(maxnobs), dsig(maxnobs)
    real :: sumsqam,sumsqph,sigma2am,sigma2ph
    real :: stddevdataam,stddevdataph
    
    integer :: nfreq, nwindows(maxevnts), nptsw(maxevnts)
    integer :: nptsh1(maxevnts), nptsh2(maxevnts), nptsz(maxevnts)
    integer :: ifreq,i,j,jj,iev,nevents,iwin,ig,jang,iwind,iw,ijk
    integer :: npts,nptsall,nerr,iwstart,np
    integer :: nobs,nogauss,nbands,iband,nbandsp
    integer :: indx(nparam)
    
    double precision :: change(nparam)
    double precision :: gtg(nparam,nparam), gtginv(nparam,nparam),ddd

    character(len=70) :: foutput, foutput2, foutput3
    character(len=70) :: BH1fn(maxevnts), BH2fn(maxevnts), BHZfn(maxevnts)


!  read in number of sets of days/stations/and/or sampling periods for which 
!  you will find transfer functions
!       read(*,*) nsets
!       do iset = 1, nsets
!  initially, just do a few frequencies for which we will find average coherence
    nfreq = 7
    do ifreq = 1, nfreq
       freq(ifreq) = (ifreq-1)*.005 + .01
    enddo

!  foutput: file for output of transfer function and coherence 
    foutput = 'tempH2Ztransfer'
    open(11, file = foutput)

!  foutput2: file for summary output of angle and polynomial fit
    foutput2 = 'tempH2Zcoeff'
    open(12, file = foutput2)

    foutput3 = 'temph'

!  read list of files to be analyzed (events or time series)
!  nogauss = 0 means don't apply gaussian window - use record as is 
!	     (For already windowed signals)
!  nogauss = 1 means do apply gaussian window (For noise sample)
!  nevents here means number of seaparate SAC files from which to 
!  extract windows - typical nevents = 1 for a 1-day file
!      read(*,*) nevents, nogauss
    nevents = 1
    nogauss = 1
    nobs = 0

!  read number of frequency bands for fitting polynomial descriptions 
!  of transfer function.  Usually this is 1 unless there is a problem 
!  with instrument.
    read(*,*) nbands

!  read frequency range for fitting polynomial descriptions of transfer function
    do iband = 1,nbands
      read(*,*) freqlo(iband), freqhi(iband)
    enddo

! repeat for dummy info for P2Z
    read(*,*) nbandsp

!  read frequency range for fitting polynomial descriptions of transfer function
    do iband = 1,nbandsp
      read(*,*) dumfreqlo, dumfreqhi
    enddo
    
evloop: do iev = 1, nevents
    ! read in horizontal BH1 file name, then BH2, then associated vertical 
    ! file name - all should have same beginning time and same time increments
    BH1fn(iev) = 'temp1'
    BH2fn(iev) = 'temp2'
    BHZfn(iev) = 'tempz'
     
    ! read in number of windows selected and timelength of window for seismogram
    read(*,*) nwindows(iev), tlength(iev)

    ! timebeg is beginning of window relative to start of record
    do iwin = 1, nwindows(iev)
      read(*,*) timebeg(iev,iwin)
    enddo

    call rsac1(BH1fn(iev),tdatah1,nptsh1(iev),begh1(iev),delth1(iev),maxpts,nerr)
    call rsac1(BH2fn(iev),tdatah2,nptsh2(iev),begh2(iev),delth2(iev),maxpts,nerr)
    call rsac1(BHZfn(iev),tdataz,nptsz(iev),begz(iev),deltz(iev),maxpts,nerr)

    if ((begh1(iev).ne.begz(iev)).or.(begh2(iev).ne.begz(iev))) then
      write (*,*) iev,' beginning times not identical'
    endif
    if ((delth1(iev).ne.deltz(iev)).or.(delth2(iev).ne.deltz(iev))) then
      write(*,*) iev,' delt not identical'
    endif

    !  set up Gaussian window function - want 3 sigma at limits of time window
    !  and amplitude = 1 at center
    nptsall = nptsz(iev)
    npts = tlength(iev)/deltz(iev) + 1
    if (nogauss.eq.1) then
        sigma = tlength(iev)/6.0
        sigma2 = 2.0*sigma*sigma
        tmiddle = tlength(iev)/2.0
        do ig = 1, npts
            tg = (ig-1)*deltz(iev) - tmiddle
            gausst(ig) = exp(-(tg*tg)/sigma2)
        enddo
    endif
  
    ! now apply gaussian to windows
    do iwin = 1, nwindows(iev)
    nobs = nobs + 1
    iwstart = timebeg(iev,iwin)/deltz(iev)
    delta(nobs) = deltz(iev)
    nptsw(nobs) = npts
      do iw = 1, npts
        if (nogauss.eq.1) then
          winh1(nobs,iw) = tdatah1(iw+iwstart)*gausst(iw)
          winh2(nobs,iw) = tdatah2(iw+iwstart)*gausst(iw)
          winz(nobs,iw) = tdataz(iw+iwstart)*gausst(iw)
        else
          winh1(nobs,iw) = tdatah1(iw+iwstart)
          winh2(nobs,iw) = tdatah2(iw+iwstart)
          winz(nobs,iw) = tdataz(iw+iwstart)
        endif
      enddo
    enddo

enddo evloop
!  end of input loop over events


          
! rotate horizontal components every 5 degrees and calculate coherence
!  with vertical
    cohermax = 0.0
    angloop: do jang = 0, 175, 5
      theta = jang*pi/180.0
      costheta = cos(theta)
      sintheta = sin(theta)
      avgcoher = 0.0

      do ifreq = 1, nfreq
        sumcrossre(ifreq) = 0.0
        sumcrossim(ifreq) = 0.0
        sumppower(ifreq) = 0.0
        sumzpower(ifreq) = 0.0
      enddo

      ! cycle through all windows, all events to find coherence
      do iwind = 1, nobs
        do i = 1, nptsw(iwind)
          vert(i) = winz(iwind,i)
          horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
        enddo

        ! get real and imaginary (phases and amplitudes) at desired frequencies
        do ifreq = 1, nfreq
          call frt2(horz,freq(ifreq),preal,pimag,nptsw(iwind),delta(iwind))
          call frt2(vert,freq(ifreq),zreal,zimag,nptsw(iwind),delta(iwind))
          crossre = preal*zreal + pimag*zimag
          crossim = preal*zimag - pimag*zreal
          ppower = preal*preal + pimag*pimag
          zpower = zreal*zreal + zimag*zimag
          sumcrossre(ifreq) = sumcrossre(ifreq) + crossre
          sumcrossim(ifreq) = sumcrossim(ifreq) + crossim
          sumppower(ifreq) = sumppower(ifreq) + ppower
          sumzpower(ifreq) = sumzpower(ifreq) + zpower
        enddo
      enddo 
 
      !  calculate coherence as function of frequency and find average
      do ifreq = 1, nfreq	
        coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/ &
                       (sumppower(ifreq)*sumzpower(ifreq))
        avgcoher = avgcoher + coher(ifreq)
      enddo

      avgcoher = avgcoher/nfreq
      if (avgcoher.gt.cohermax) then
         maxtheta = jang
         cohermax = avgcoher
      endif

    enddo angloop
    
    !  now search every 1 degree to find better estimate of maximum
    !  coherence direction between horizontal and vertical	
    angloop2: do jang = maxtheta-3,maxtheta+3, 1
      theta = jang*pi/180.0
      costheta = cos(theta)
      sintheta = sin(theta)
      avgcoher = 0.0

      do ifreq = 1, nfreq
        sumcrossre(ifreq) = 0.0
        sumcrossim(ifreq) = 0.0
        sumppower(ifreq) = 0.0
        sumzpower(ifreq) = 0.0
      enddo

      ! cycle through all windows, all events to find coherence
      do iwind = 1, nobs
        do i = 1, nptsw(iwind)
          vert(i) = winz(iwind,i)
          horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
        enddo

        ! get real and imaginary (phases and amplitudes) at desired frequencies
        do ifreq = 1, nfreq
          call frt2(horz,freq(ifreq),preal,pimag,nptsw(iwind),delta(iwind))
          call frt2(vert,freq(ifreq),zreal,zimag,nptsw(iwind),delta(iwind))
          crossre = preal*zreal + pimag*zimag
          crossim = preal*zimag - pimag*zreal
          ppower  = preal*preal + pimag*pimag
          zpower  = zreal*zreal + zimag*zimag
          sumcrossre(ifreq) = sumcrossre(ifreq) + crossre
          sumcrossim(ifreq) = sumcrossim(ifreq) + crossim
          sumppower(ifreq)  = sumppower(ifreq) + ppower
          sumzpower(ifreq)  = sumzpower(ifreq) + zpower
        enddo
      enddo  

      ! calculate coherence as function of frequency and find average 
      do ifreq = 1, nfreq
        qreal(ifreq) = sumcrossre(ifreq)/sumppower(ifreq)
        qimag(ifreq) = sumcrossim(ifreq)/sumppower(ifreq)
        coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/ &
                       (sumppower(ifreq)*sumzpower(ifreq))
        avgcoher = avgcoher + coher(ifreq)
        call phase(qreal(ifreq),qimag(ifreq),qphse(ifreq))
        if (qphse(ifreq).gt.0.5) qphse(ifreq)=qphse(ifreq)-1.0
      enddo

      avgcoher = avgcoher/nfreq
      if (avgcoher.gt.cohermax) then
        maxtheta = jang
        cohermax = avgcoher
      endif

    enddo angloop2

    !  check to see if better to rotate by 180 degrees for optimum direction 
    !  to get phase close to zero instead of close to pi
    if (abs(qphse(4)).gt. 0.3) maxtheta = maxtheta + 180.
    ! phideg = phi * 180./pi
    write(*,*) maxtheta, cohermax 
    write(11,*) maxtheta, cohermax
    write(12,*) maxtheta, cohermax

    !at optimum direction, calculate transfer function
    theta = maxtheta*pi/180.
    costheta = cos(theta)
    sintheta = sin(theta)
    avgcoher = 0.0

    ! make rotated single horizontal for output
    do ijk = 1,nptsall 
      tdatahr(ijk) = tdatah1(ijk)*costheta+tdatah2(ijk)*sintheta
    enddo
    call wsac0(foutput3,xdummy,tdatahr,nerr)

    nfreq = 396
    do ifreq = 1, nfreq
      freq(ifreq) = (ifreq-1)*.001 + .005
      sumcrossre(ifreq) = 0.0
      sumcrossim(ifreq) = 0.0
      sumppower(ifreq) = 0.0
      sumzpower(ifreq) = 0.0
    enddo

    ! cycle through all windows, all events to find coherence
    do iwind = 1, nobs
       do i = 1, nptsw(iwind)
         vert(i) = winz(iwind,i)
         horz(i) = winh1(iwind,i)*costheta+winh2(iwind,i)*sintheta
       enddo

       ! get real and imaginary (phases and amplitudes) at desired frequencies
       do ifreq = 1, nfreq
          call frt2(horz,freq(ifreq),preal,pimag,nptsw(iwind),delta(iwind))
          call frt2(vert,freq(ifreq),zreal,zimag,nptsw(iwind),delta(iwind))
          crossre = preal*zreal + pimag*zimag
          crossim = preal*zimag - pimag*zreal
          ppower = preal*preal + pimag*pimag
          zpower = zreal*zreal + zimag*zimag

          sumcrossre(ifreq) = sumcrossre(ifreq) + crossre
          sumcrossim(ifreq) = sumcrossim(ifreq) + crossim
          sumppower(ifreq) = sumppower(ifreq) + ppower
          sumzpower(ifreq) = sumzpower(ifreq) + zpower
        enddo
    enddo
  
    write(11,*) 'freq   admittance std(admit)  phase   std(ph) coher      real        imag'  

    do ifreq = 1, nfreq
      qreal(ifreq) = sumcrossre(ifreq)/sumppower(ifreq)
      qimag(ifreq) = sumcrossim(ifreq)/sumppower(ifreq)
      coher(ifreq) = (sumcrossre(ifreq)**2 + sumcrossim(ifreq)**2)/ &
                     (sumppower(ifreq)*sumzpower(ifreq))

      ! epsilon is fractional uncertainty in amplitude of transfer function (1 stddev)
      ! and also approximately 1 stddev in phase (in radians) Bendat and Piersol page 317
      epsilon(ifreq) = sqrt((1.0 - coher(ifreq))/(2.0*coher(ifreq)*nobs))

      call phase(qreal(ifreq),qimag(ifreq),qphse(ifreq))

      if (qphse(ifreq).gt.0.5) qphse(ifreq)=qphse(ifreq)-1.0

      qamp(ifreq) = sqrt(qreal(ifreq)**2 + qimag(ifreq)**2)
      ampdev(ifreq) = qamp(ifreq)*epsilon(ifreq)
      !  phase deviation and phase in cycles
      phsedev(ifreq) = epsilon(ifreq)/twopi

      !  output 
      write(11,100) freq(ifreq),qamp(ifreq),ampdev(ifreq),&
          qphse(ifreq),phsedev(ifreq),coher(ifreq),& 
          qreal(ifreq), qimag(ifreq)
100   format(f7.4,2(1x,e11.4),3(1x,f7.4),2(1x,e11.4))      
    enddo

    
! ***************************************************************************
!  Now calculate polynomial (quadratic) fits to amplitude and phase of transfer
!  function over specified frequency range using weighted least squares

!  Parameter 1 is constant, parameter 2 is slope
!  and parameter 3 is quadratic term for amplitude.  4,5,6 same for phase

!  Loop through frequency bands for polynomial fits
bandloop:  do iband = 1, nbands

!  set up normalized data vector, i.e., divide each observation by the 
!  assigned standard deviation
    nobs = 0
    do ifreq = 1, nfreq
      !  check whether within desired frequency range
      if ((freq(ifreq).ge.freqlo(iband)).and.(freq(ifreq).le.freqhi(iband))) then
      nobs = nobs+2
      dsig(nobs-1) = ampdev(ifreq)
      d(nobs-1) = qamp(ifreq)/dsig(nobs-1)
      dsig(nobs) = phsedev(ifreq)
      d(nobs) = qphse(ifreq)/dsig(nobs)

      !  set up partial derivative matrix, again normalizing by dividing 
      !  by the standard deviation
      g(nobs-1,1) = 1./dsig(nobs-1)
      g(nobs-1,2) = freq(ifreq)/dsig(nobs-1)
      g(nobs-1,3) = freq(ifreq)*freq(ifreq)/dsig(nobs-1)
      g(nobs-1,4) = 0.0
      g(nobs-1,5) = 0.0
      g(nobs-1,6) = 0.0
      g(nobs,1) = 0.0
      g(nobs,2) = 0.0
      g(nobs,3) = 0.0
      g(nobs,4) = 1./dsig(nobs)
      g(nobs,5) = freq(ifreq)/dsig(nobs)
      g(nobs,6) = freq(ifreq)*freq(ifreq)/dsig(nobs)
      !	write(*,*) nobs, freq(ifreq), d(nobs-1),d(nobs)
      endif
    enddo

!  Calculate gtg and gtd

    np = 6
    do j = 1, np
      gtd(j) = 0.0
      do i = 1, nobs
          gtd(j) = gtd(j) + g(i,j)*d(i)
      enddo
      do jj = 1,j
          gtg(jj,j) = 0.0
          do i = 1, nobs
            gtg(jj,j)= gtg(jj,j) + g(i,jj)*g(i,j)
          enddo
          gtg(j,jj) = gtg(jj,j)
      enddo
    enddo

!  Invert gtg.  gtg will be destroyed.  
!  Not the most efficient approach because doesn't take advantage of 
!  symmetry of gtg.  Use LU decomposition from Press et al.
    do i= 1,np
        do j = 1, np
          gtginv(i,j)= 0.0D0
        enddo
        gtginv(i,i) =1.0D0
    enddo

    call dludcmp(gtg,np,nparam,indx,ddd)

    do j = 1,np
        call dlubksb(gtg,np,nparam,indx,gtginv(1,j))
    enddo

    !  Find solution
    do i= 1, np
        change(i)=0.0
        do j = 1,np
          change(i) =change(i) + gtd(j)*gtginv(i,j)
        enddo
        ! write(*,*), i, change(i)
    enddo


    !  Find normalized residuals  (linear problem with zero start 
    !  so just use partial derivatives) and sum of squares of errors
    sumsqam = 0.0
    sumsqph = 0.0
    do i = 2, nobs, 2
      resid(i) = d(i) - g(i,4)*change(4) - g(i,5)*change(5) - g(i,6)*change(6)
      sumsqph = sumsqph + resid(i)**2
      resid(i-1) = d(i-1) - g(i-1,1)*change(1) - g(i-1,2)*change(2) - g(i-1,3)*change(3)
      sumsqam = sumsqam + resid(i-1)**2
    enddo

    sigma2am = sumsqam/(nobs/2-3)
    stddevdataam = sqrt(sigma2am)
    sigma2ph = sumsqph/(nobs/2-3)
    stddevdataph = sqrt(sigma2ph)

    write(11,*) stddevdataam, stddevdataph, ' amp and phase normalized standard deviations'
    write(11,*) freqlo(iband), freqhi(iband)
    write(11,*) change(1),change(2),change(3)
    write(11,*) change(4),change(5),change(6)
    if (iband.eq.1) write(*,'(a)') foutput
!	write(*,*) stddevdataam, stddevdataph, ' amp and phase 
!     1 normalized standard deviations'
!	write(*,*) freqlo(iband), freqhi(iband)
!        write(*,*) change(1),change(2),change(3)
!	write(*,*) change(4),change(5),change(6)
    if (iband.eq.1) write(12,*) nbands
    write(12,*) stddevdataam, stddevdataph, ' amp and phase normalized standard deviations'
    write(12,*) freqlo(iband), freqhi(iband)
    write(12,*) change(1),change(2),change(3)
    write(12,*) change(4),change(5),change(6)
!  end loop over frequency bands
enddo bandloop

close(unit = 11)
close(unit = 12)
    
end program


!---positive fourier transform

SUBROUTINE FRT(UP,FR,ARZ,PRZ,NALL,DELT)
    DIMENSION UP(1000000)
    DIMENSION W(1000000)
    THETA=6.283185*FR*DELT
    C=COS(THETA)*2.0
    NR1=1
    NR2=NALL
    NDR1=NR2-1
    W(1)=UP(NR2)
    W(2)=C*W(1)+UP(NDR1)
    NDATA=NR2-NR1+1
    NTR1=NDATA-1
    NTR2=NDATA-2
    DO I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
    ENDDO
    ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
    ZIM=W(NTR1)*SIN(THETA)*DELT
    CALL PHASE(ZRE,ZIM,PRZ)
    ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
    RETURN
END


SUBROUTINE FRT2(UP,FR,ZRE,ZIM,NALL,DELT)
    DIMENSION UP(1000000)
    DIMENSION W(1000000)
    REAL :: ZRE,ZIM,PRZ
    THETA=6.283185*FR*DELT
    C=COS(THETA)*2.0
    NR1=1
    NR2=NALL
    NDR1=NR2-1
    W(1)=UP(NR2)
    W(2)=C*W(1)+UP(NDR1)
    NDATA=NR2-NR1+1
    NTR1=NDATA-1
    NTR2=NDATA-2
    DO I=3,NDATA
      I1=I-1
      I2=I-2
      NDRI=NR2-I+1
      W(I)=C*W(I1)-W(I2)+UP(NDRI)
    ENDDO
    ZRE=(W(NDATA)-W(NTR2)+UP(NR1))*DELT/2.0
    ZIM=W(NTR1)*SIN(THETA)*DELT
!      CALL PHASE(ZRE,ZIM,PRZ)
!      ARZ=SQRT(ZRE*ZRE+ZIM*ZIM)
    RETURN
END

SUBROUTINE PHASE(X,Y,PHI)
    REAL :: x,y,phi
    IF(X) 21,20,22
20    IF(Y) 23,24,25
23    PHI=1.5*3.141592
    GO TO 28
24    PHI=0.0
    GO TO 28
25    PHI=0.5*3.141592
    GO TO 28
21    PHI=ATAN(Y/X) +3.141592
    GO TO 28
22    IF(Y) 26,27,27
27    PHI=ATAN(Y/X)
    GO TO 28
26    PHI=ATAN(Y/X)+2.0*3.141592
    GO TO 28
28    CONTINUE
    PHI=PHI/6.283184
    PHI=PHI-AINT(PHI)
    RETURN
END

 
subroutine dlubksb(a,n,np,indx,b)
    implicit none

    integer :: n,np,indx(n)
    integer :: i,ii,j,ll
    double precision :: a(np,np),b(n)
    double precision :: sum

    ii=0
    do i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)

      if (ii.ne.0)then
        do j=ii,i-1
          sum=sum-a(i,j)*b(j)
        enddo
      else if (sum.ne.0.) then
        ii=i
      endif

      b(i)=sum
    enddo

    do i=n,1,-1
      sum=b(i)
      do j=i+1,n
        sum=sum-a(i,j)*b(j)
      enddo
      b(i)=sum/a(i,i)
    enddo

    return
end



subroutine dludcmp(a,n,np,indx,d)

    integer, parameter :: NMAX=500
    double precision, parameter :: TINY=1.0e-20

    integer :: n,np,indx(n)
    integer :: i,imax,j,k
    double precision :: d,a(np,np)
    double precision :: aamax,dum,sum,vv(NMAX)

    d=1.

    do i=1,n
        aamax=0.

        do j=1,n
            if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo

        if (aamax.eq.0.) pause 'singular matrix in ludcmp'

        vv(i)=1./aamax
    enddo

    do j=1,n
        do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
                sum=sum-a(i,k)*a(k,j)
            enddo
            a(i,j)=sum
        enddo

        aamax=0.

        do i=j,n
            sum=a(i,j)
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            enddo

            a(i,j)=sum
            dum=vv(i)*abs(sum)

            if (dum.ge.aamax) then
                imax=i
                aamax=dum
            endif
        enddo

        if (j.ne.imax)then
            do k=1,n
                dum=a(imax,k)
                a(imax,k)=a(j,k)
                a(j,k)=dum
            enddo

            d=-d
            vv(imax)=vv(j)
        endif

        indx(j)=imax

        if(a(j,j).eq.0.) a(j,j)=TINY
        if(j.ne.n) then
            dum=1./a(j,j)
            do i=j+1,n
                a(i,j)=a(i,j)*dum
            enddo
        endif

    enddo

    return

end

