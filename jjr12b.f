
!djh: top is maximum lowest principal quantum number
!     max is maximum highest principal quantum number
	parameter(itopfinal=18, maxfinal=4)
! maab: Any changes in the number of partial waves and other parameters
!       have to reflect here.
!		ncr == number of closed type orbitals
!		nor == number of open type orbitals
!		nljmax == number of different partial waves orbitals
        parameter(ncr=28,nor=2,nr=ncr+nor,mm=4321,nblk=90)
        parameter(nopmax=40,nop2=nopmax*2,nls=60,nkval=60,ntwo=9000)
        parameter(nljmax=13,norbs=nljmax*nr)
	implicit real*8(a-h,o-z)
	dimension lb(norbs),jb(norbs),n1(ntwo),n2(ntwo)
	dimension l1(ntwo),l2(ntwo),j1(ntwo),j2(ntwo)
	dimension kkstart(nblk),kkend(nblk),k1o(nop2),k2o(nop2)
	dimension l1blk(nblk),l2blk(nblk),j1blk(nblk),j2blk(nblk)
	dimension xrecp(nls)
	dimension dirang(nls,nkval),excang(nls,nkval)
	dimension direct(nkval),exch(nkval),r12(ntwo,ntwo)
	dimension y(mm,norbs)
	common/radial/r(mm),rp(mm),rpor(mm),xh
 
	call setup
	one=1.d0
	write(6,*) "START r12b"
	write(6,*) "Input Final state J"
	read(5,*) jin
	write(6,*) 'Input final state parity'
	read(5,*) iod
	open(unit=55, file='finalJ.bin', status='replace')
	write(55,*) jin
	write(55,*) iod
	close(55)
	open(unit=19,file='wfn19.dat',status='old',
     1		access='sequential',form='unformatted')
	rewind(19)
	read(19) m,r(1),xh
	read(19) ((y(j,i),j=1,m),i=1,norbs)
	close(19)
	h=dsqrt(r(1))
	write(6,*) m,r(1),xh,h
	do 10 i=2,m
	   r(i)=r(1)*i*i
	   rpor(i)=1.d0
	   rp(i)=r(i)
   10   continue
	mx=m
	write(6,*) 'r0=',r(m)
 
c ***** input open basis information
	open(unit=5,file='RaFinalStateJ'//char(jin+48)//'.dat'
     1       ,status='old')
        rewind(5)
	read(5,*) ipar,j
	jd=j*2
	read(5,*) nkopen,nchan
	do 110 i=1,nkopen
		read(5,*) k1o(i),k2o(i)
  110   continue
        close(5)
	nopper=nkopen/nchan
 
c ***** jj-coupling basis set key
	do 120 ilj=1,nljmax
	  do 118 i=1,nr
	    k=i+(ilj-1)*nr
	    if(ilj.eq.1) then
	      lb(k)=0
	      jb(k)=1
	    endif
	    if(ilj.gt.1) then
	      lb(k)=ilj/2
	      lb(k)=2*lb(k)
	      int=1-(-1)**ilj
	      jb(k)=lb(k)+int-1
	    endif
  118   continue
  120   continue
 
	write(2,*) 'one-electron basis set key: ',norbs
	write(2,4477) (lb(k),jb(k),k=1,norbs)
 4477   format(5(2i3,3x))
 
 
 
c ******************  AUTOMATIC two-electron basis set selector|
 
	write(2,*) '**'
	write(2,*) 'two-electron configurations:'
	kkopen=0
	kk=0
	iblocks=0
		l1prev=-1
		l2prev=-1
		j1prev=-1
		j2prev=-1
 
	do 140 lj1=1,nljmax
	do 139 lj2=1,nljmax
	do 138 ir1=1,nr
	do 137 ir2=1,nr
 
c
c
 
	  k1 = ir1 + (lj1-1)*nr
	  k2 = ir2 + (lj2-1)*nr
	  if(k2.lt.k1) go to 137

c *************************************************************
c
c ** Use the following lines to skip configurations with
c ******* both electrons very highly excited:
c

!djh
	  if ( (ir1.gt.itopfinal) .or. (ir2.gt.itopfinal) ) go to 137
c	  if ( (ir1.gt.maxfinal) .and. (ir2.gt.maxfinal) ) go to 137



!previous          if(ir1.gt.4 .and. (ir2-ir1).gt.3) go to 137
!          if(ir2.gt.4 .and. (ir1-ir2).gt.3) go to 137
!          if(lb(k1)+lb(k2).ge.2*5 .and.ir1+ir2.gt.12) go to 137
!          if(lb(k1).ge. 2*4 .and. ir1.gt.10) go to 137
!          if(lb(k2).ge. 2*4 .and. ir2.gt.10) go to 137
!          if(ir1.gt.18 .and. ir2.gt.3) go to 137
!          if(ir2.gt.18 .and. ir1.gt.3) go to 137
!          if(ir1.gt.13 .and. lb(k1).gt.2*3) go to 137
!          if(ir2.gt.13 .and. lb(k2).gt.2*3) go to 137
!          if(ir1.gt.15 .and. lb(k1).eq.2*2) go to 137
!          if(ir2.gt.15 .and. lb(k2).eq.2*2) go to 137


c	  if(ir1.gt.2 .and. (ir2-ir1).gt.14) go to 137
c  	  if(ir2.gt.2 .and. (ir1-ir2).gt.14) go to 137
c	  if(lb(k1).ge. 2*3 .and. ir1.gt.16) go to 137
c	  if(lb(k2).ge. 2*3 .and. ir2.gt.16) go to 137
c	  if(lb(k1).ge. 2*4 .and. ir1.gt.15) go to 137
c	  if(lb(k2).ge. 2*4 .and. ir2.gt.15) go to 137
c
c**************************************************************

 
	  ip=(-1)**( (lb(k1)+lb(k2))/2 )
	  jgo=1
	  if(iabs(jb(k1)-jb(k2)) .gt. jd) jgo=0
	  if(jb(k1)+jb(k2) .lt. jd) jgo=0
 
          iPauli=(-1)**((jb(k1)+jb(k2)+jd)/2)
          iSymmTest=0
          IF(ip.EQ.ipar) iSymmTest=1
          IF(ir1.EQ.ir2 .AND. lb(k1).EQ.lb(k2)) THEN
          IF(jb(k1).EQ.jb(k2) .AND. iPauli.EQ.1) THEN
            iSymmTest=0
          ENDIF
          ENDIF
 
	       iopcl1=k1 - ((k1-1)/nr)*nr
	       iopcl2=k2 - ((k2-1)/nr)*nr
	       iop1=0
	       if(iopcl1.gt.ncr) iop1=1
	       iop2=0
 
	       if(iopcl2.gt.ncr) iop2=1
	       iopgo=0
	       if(iop1+iop2 .eq.0) iopgo=1
	  if(iSymmTest.EQ.1 .and. jgo.eq.1 .and. iopgo.eq.1) then
 
	     kk=kk+1
	     n1(kk)=k1
	     n2(kk)=k2
	     l1(kk)=lb(k1)
	     l2(kk)=lb(k2)
	     j1(kk)=jb(k1)
	     j2(kk)=jb(k2)
	if(kk.ge.ntwo-1) write(6,*) 'kk,ntwo',kk,ntwo
 
	     itst=iabs(l1(kk)-l1prev)+iabs(l2(kk)-l2prev)
	     itst=itst+iabs(j1(kk)-j1prev)+iabs(j2(kk)-j2prev)
	     if(itst.ne.0) then
		if(iblocks .ne.0) kkend(iblocks)=kk-1
		iblocks=iblocks+1
		l1prev=l1(kk)
		l2prev=l2(kk)
		j1prev=j1(kk)
		j2prev=j2(kk)
		l1blk(iblocks)=l1(kk)
		l2blk(iblocks)=l2(kk)
		j1blk(iblocks)=j1(kk)
		j2blk(iblocks)=j2(kk)
		kkstart(iblocks)=kk
	     endif
		kkend(iblocks)=kk
	 endif
 
  137   continue
  138   continue
  139   continue
  140   continue
 
	ncc=kk
 
	write(2,*) ncc,' closed configurations'
	write(6,*) ncc,' closed configurations'
	do 4488 k=1,ncc
	   write(2,*) k,': ',n1(k),l1(k),j1(k),'- ',n2(k),l2(k),j2(k)
 4488   continue
 
c ***** Note that iblocks = number of distinct angular configurations
        IF(iblocks.GT.nblk) THEN
            WRITE(6,*) 'STOP:iblocks,nblk=',iblocks,nblk
            STOP
        ENDIF
 
c ***** Now add open-type configurations as read in above
 
	i=0
	do 160 ii=1,nchan
	do 158 ip=1,nopper
	    i=i+1
	    kk=kk+1
	    n1(kk)=k1o(i)
	    n2(kk)=k2o(i)
	    l1(kk)=lb(k1o(i))
	    l2(kk)=lb(k2o(i))
	    j1(kk)=jb(k1o(i))
	    j2(kk)=jb(k2o(i))

          IF(iabs(j1(kk)-j2(kk)).gt.jd) THEN
              WRITE(6,*) 'CHANNEL ',ii,' fails triang-.'
              STOP
          ENDIF
          IF((j1(kk)+j2(kk)).lt.jd) THEN
              WRITE(6,*) 'CHANNEL ',ii,' fails triang+.'
              STOP
          ENDIF
	    ipTEST=(-1)**( (l1(kk)+l2(kk))/2 )
          IF(ipTEST.NE.ipar) THEN
             WRITE(6,*) 'CHANNEL ',ii,' fails parity check'
             STOP
          ENDIF

  158   continue
	  iblocks=iblocks+1
	  l1blk(iblocks)=l1(kk)
	  l2blk(iblocks)=l2(kk)
	  j1blk(iblocks)=j1(kk)
	  j2blk(iblocks)=j2(kk)
	  kkstart(iblocks)=ncc + (ii-1)*nopper + 1
	  kkend(iblocks)=ncc + ii*nopper
  160   continue
 
	noo=nkopen
	ntot=ncc + noo
 
	write(2,*) noo,' open configurations'
	do 4499 k=ncc+1,ntot
	  write(2,*) k,': ',n1(k),l1(k),j1(k),'- ',n2(k),l2(k),j2(k)
 4499   continue
 
	write(2,*) ' '
	write(2,*) iblocks,' is number of distinct angular configs'
	do 3377 ib=1,iblocks
	   write(2,*) ib,l1blk(ib),j1blk(ib),'- ',l2blk(ib),j2blk(ib)
 3377   continue
	close(2)

CCC
	open(unit=8,file='config1.dat',status='unknown')
	write(8,181) ncc,noo,iblocks,nljmax,nr,ipar,j
  181	format(12i5)
	nnorbs=nr*nljmax
	ntot=ncc+noo
	write(8,181) (lb(k),jb(k),k=1,nnorbs)
	write(8,181) (n1(k),l1(k),j1(k),n2(k),l2(k),j2(k),k=1,ntot)
	write(8,181) (l1blk(ib),j1blk(ib),ib=1,iblocks)
      write(8,181) (l2blk(ib),j2blk(ib),ib=1,iblocks)
	write(8,181) (kkstart(ib),kkend(ib),ib=1,iblocks)
	close(8)
CCC

 
c ***** Calculate 1/r12 integrals in jj-coupling
	itotint=0
	do 300 ib1=1,iblocks
	do 300 ib2=ib1,iblocks
 
c ***** First set up the range of L-values to sum over
cc        in the r12-matrix calculation
 
	   lmin1=iabs(l1blk(ib1)-l2blk(ib1))
	   lmin2=iabs(l1blk(ib2)-l2blk(ib2))
	   llmin=max0(lmin1,lmin2)
	   lmax1=l1blk(ib1)+l2blk(ib1)
	   lmax2=l1blk(ib2)+l2blk(ib2)
	   llmax=min0(lmax1,lmax2)
 
	   iterm=0
	   do 220 is=0,2,2
	     do 219 ll=llmin,llmax,2
 
		jj1=j1blk(ib1)
		jj2=j2blk(ib1)
		ll1=l1blk(ib1)
		ll2=l2blk(ib1)
		front1=(is+1)*(ll+1)*(jj1+1)*(jj2+1)
		recoup1=xninej(1,ll1,jj1,1,ll2,jj2,is,ll,jd)
		recoup1=dsqrt(front1)*recoup1
 
		jj1p=j1blk(ib2)
		jj2p=j2blk(ib2)
		ll1p=l1blk(ib2)
		ll2p=l2blk(ib2)
		l1t=ll1/2
		l2t=ll2/2
		l1pt=ll1p/2
		l2pt=ll2p/2
		ltrue=ll/2
 
		front2=(is+1)*(ll+1)*(jj1p+1)*(jj2p+1)
		recoup2=xninej(1,ll1p,jj1p,1,ll2p,jj2p,is,ll,jd)
		recoup2=dsqrt(front2)*recoup2
 
		iterm=iterm+1
		xrecp(iterm)=recoup1*recoup2
 
c****   store the multipole integrals for each L,k
c********       Direct angular integrals, then exchange
 
		kmin=iabs(ll1-ll1p)
		kmax=ll1+ll1p
		ktest=iabs(ll2-ll2p)
		if(ktest.gt.kmin) kmin=ktest
		ktest=ll2+ll2p
		if(ktest.lt.kmax) kmax=ktest
 
		   do 212 k=kmin,kmax,4
		     ktrue=k/2
	   dirang(iterm,k/4+1)=coeff(l1t,l2t,l1pt,l2pt,ltrue,ktrue)
  212              continue
 
		kminx=iabs(ll1-ll2p)
		ktest=iabs(ll2-ll1p)
		kminx=max0(kminx,ktest)
 
		kmaxx=ll1+ll2p
		ktest=ll2+ll1p
		kmaxx=min0(kmaxx,ktest)
 
		do 216 k=kminx,kmaxx,4
		  ktrue=k/2
	   excang(iterm,k/4+1)=coeff(l2t,l1t,l1pt,l2pt,ltrue,ktrue)
  216           continue
  219   continue
  220   continue
 
c ****  Fold in the radial integrals and everything
 
	do 270 k1=kkstart(ib1),kkend(ib1)
	  kkgo=kkstart(ib2)
	  if(ib1.eq.ib2) kkgo=k1
	  do 268 k2=kkgo,kkend(ib2)
	    kr1=n1(k1)
	    kr2=n2(k1)
	    kr1p=n1(k2)
	    kr2p=n2(k2)
 
c ****  Direct term first, then exchange
 
	      do 250 k=kmin,kmax,4
		ktrue=k/2
		direct(k/4+1)=sl(y(1,kr1),y(1,kr2),y(1,kr1p),y(1,kr2p),
     x                       ktrue,l2t,l4t,one,one,m,mx)
  250         continue
 
	      do 252 k=kminx,kmaxx,4
		ktrue=k/2
		exch(k/4+1)=sl(y(1,kr2),y(1,kr1),y(1,kr1p),y(1,kr2p),
     x                        ktrue,l1t,l4t,one,one,m,mx)
  252         continue
 
	    iterm=0
	    totdir=0.d0
	    totexc=0.d0
	    do 260 is=0,2,2
	    do 259 ll=llmin,llmax,2
	       iterm=iterm+1
 
	       do 257 k=kmin,kmax,4
	 totdir=totdir+direct(k/4+1)*dirang(iterm,k/4+1)*xrecp(iterm)
  257          continue
 
	       totkkk=0.d0
	       do 258 k=kminx,kmaxx,4
	 totkkk=totkkk+exch(k/4+1)*excang(iterm,k/4+1)*xrecp(iterm)
  258          continue
	       totkkk=totkkk*ipar* (-1)**( (is+ll)/2 )
	       totexc=totexc+totkkk
  259       continue
  260       continue
 
	    xn=delt(kr1,kr2)
	    xn=1.d0/dsqrt(1.d0+xn)
	    xnp=delt(kr1p,kr2p)
	    xnp=1.d0/dsqrt(1.d0+xnp)
	    r12int=(totdir+totexc)*xn*xnp
	    r12(k1,k2)=r12int
	    r12(k2,k1)=r12int
 
	    itotint=itotint+1
	    ntotint=((itotint+1)*itotint)/2
cccccc	    write(6,*) 'ib: ',ib1,ib2,' k vals ',k1,k2,itotint,r12int
 
  268      continue
  270      continue
 
  300   continue
	open(unit=12,file='final12.dat',status='unknown')
	rewind(12)
	write(12,*) ntot
	write(12,4) ((r12(k1,k2),k1=1,k2),k2=1,ntot)
	close(12)
cc OLDformat-->    4   format(5(1pd16.9))
   4    format(1p3d22.15)
c *************************
	write(6,*) 'r12b calculation finished'
	end
