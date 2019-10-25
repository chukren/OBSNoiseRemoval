c  combine two horizontals into one component in direction of vertical tilt
      parameter (maxpts = 1000000)
      real*4 tdatah1(maxpts), tdatah2(maxpts),tdatahr(maxpts)      
      character*70 foutput3
      
      character*70 BH1fn, BH2fn
        pi = 3.14159
        BH1fn = 'temp1'
	BH2fn = 'temp2'
        foutput3 = 'temph'
	call rsac1(BH1fn,tdatah1,nptsh1,begh1,
     1       delth1,maxpts,nerr)
        call rsac1(BH2fn,tdatah2,nptsh2,begh2,
     1       delth2,maxpts,nerr)
        nptsall = nptsh1
        read(*,*) maxtheta
        theta = maxtheta*pi/180.
	costheta = cos(theta)
	sintheta = sin(theta)
c  make rotated single horizontal for output
        do ijk = 1,nptsall 
	  tdatahr(ijk) = tdatah1(ijk)*costheta
     1                   +tdatah2(ijk)*sintheta
	enddo
	call wsac0(foutput3,xdummy,tdatahr,nerr)
        end
	
