	parameter(nr=30,nljmax=13,norbs=nr*nljmax,nls=99)
	parameter(ntwo=9000,nopmax=40,nop2=nopmax*2)
	parameter(nblk=90,nEig0=20)
	implicit real*8(a-h,o-z)
	real*4 time(2)

        logical, parameter :: expth_flag = .false.
        integer, parameter :: num_expth = 5
c       Thresholds for the Ra atom in the range we want
c       From NIST table
c        real*8, parameter :: expth(5)=(/-0.4475768d0, 
c     1                                  -0.3242741d0, 
c     1                                  -0.3090996d0,
c     1                                  -0.3429597d0,
c     1                                  -0.3367089d0/)
c	    real*8, parameter :: expgnd = -0.67741295d0
c       From our fits
        real*8, parameter :: expth(5)=(/-0.3729019d0,   !S12
     1                                  -0.2756181d0,   !P12
     1                                  -0.2534856d0,   !P32
     1                                  -0.317842d0,   !D32
     1                                  -0.3102842d0/)  !D52
	    real*8, parameter :: expgnd = -0.566880d0     !GS
		character(len=6) :: ident_flag
		integer :: jint, jfin
        dimension nodes(norbs)
	dimension en(norbs),dy(norbs),y(norbs),ov(norbs,norbs)
	dimension rm(norbs,norbs),rin(norbs,norbs),dr(norbs,norbs)
	dimension n10(ntwo),n20(ntwo),l1b0(nblk),l2b0(nblk)
	dimension j1b0(nblk),j2b0(nblk),kkst0(nblk),kkend0(nblk)
      dimension vec00(ntwo,nEig0),eval0(nEig0)
	dimension vec0(ntwo),n1f(ntwo),n2f(ntwo),l1bf(nblk)
	dimension l2bf(nblk),j1bf(nblk),j2bf(nblk),kkstf(nblk)
	dimension kkendf(nblk),elam(ntwo),x(ntwo,ntwo)
	dimension hlco(ntwo,nop2),hloc(nop2,ntwo),olco(ntwo,nop2)
	dimension oloc(nop2,ntwo),hloo(nop2,nop2),oloo(nop2,nop2)
	dimension xrecp(nls),dlen(ntwo,ntwo),cc(4,nls),dd(4,nls)
	dimension omega(nop2,nop2),qlam(nop2,nop2),alfr(nop2)
	dimension alfi(nop2),bvec(nop2),vec(nop2,nop2),bvar(nop2)
	dimension vec1(nop2,nop2),zvec(ntwo,nopmax),dvel(ntwo,ntwo)
	dimension dbetal(nopmax),dbetav(nopmax),ps(nopmax,nopmax)
	dimension dps(nopmax,nopmax),eps(nopmax),lchan(nopmax)
	dimension f(nopmax),df(nopmax),g(nopmax),dg(nopmax)
	dimension qi(nopmax,nopmax),qj(nopmax,nopmax),cm(nopmax,nopmax)
	dimension dm(nopmax),xk(nopmax,nopmax),xmu(nopmax)
	dimension umat(nopmax,nopmax),dalphl(nopmax),dalphv(nopmax)
	dimension xgam(nopmax,nopmax),xlam(nopmax,nopmax)
	dimension avec(nopmax,nopmax),tao(nopmax),xnorm(nopmax)
	dimension drhol(nopmax),drhov(nopmax),ddlen(ntwo),ddvel(ntwo)

c!c! djh (see above)	dimension expth(nopmax),

        dimension Tmat(nopmax,nopmax),PartialL(nopmax),PartialV(nopmax)

	dimension epsexp(nopmax),beta(nopmax)
	dimension ipivot(nopmax),ind1(nopmax),ind2(nopmax)
	common/agqdt/aqdt,gqdt,f0,fp0,g0,gp0
c -- Thresholds for Yb:
c   J=0
c     data expth/-0.4475768d0,
c     1          -0.3429597d0,   
c     1          -0.3367089d0,
c     1          -0.3242741d0,
c     1          -0.3090996d0/
c   J=
c     data expth/-0.4475768d0,-0.4475768d0,
c     1          -0.3429597d0,-0.3429597d0,-0.3429597d0,   
c     1          -0.3367089d0,-0.3367089d0,-0.3367089d0,
c     1          -0.3242741d0,-0.3242741d0,
c     1          -0.3090996d0,-0.3090996d0,-0.3090996d0/
c	data expgnd/-0.67741295/
c -- Thresholds for Sr:
c	data expth/-0.4053499,-0.4053499,
c    1             -0.3390283,-0.3390283,-0.3390283,
c    1             -0.3377510,-0.3377510,-0.3377510,
c    1             -0.2972955,-0.2972955,
c    1             -0.2936438,-0.2936438,-0.2936438/
c	data expgnd/-0.614632353d0/

c -- Thresholds for Ba:
c	data expth/-0.36763789,
c    1             -0.34543090,
c    1             -0.34178138,
c    1		   -0.27531905,
c    1             -0.26761490/

cbg  1             -0.27543173,-0.27543173,
cbg  1             -0.26772747,-0.26772747,-0.26772747/
c -- above two lines are revised thresholds according to Brown & Ginter
ccc	data expgnd/-0.559162786d0/
c -- Ba 6s6p 1Po level:
c	data expgnd/-0.4768728993d0/

c -- Radium thresholds:
c -- Ra ground state:
c	data expgnd/-0.566897764d0/

c -- Calcium thresholds:
c	data expth/-0.43623565,-0.43623565,
c    1		-0.37408275,-0.37408275,-0.37408275,-0.37408275,
c    1		-0.3738062,-0.3738062,-0.3738062,-0.3738062,
c    1		-0.32035008744/
c	data expgnd/-0.660932093d0/
c -- Calcium thresholds:
c        data expth/-0.43623565,-0.43623565,
c     1          -0.97408275,-0.97408275,-0.97408275,-0.97408275,
c     1          -0.9738062,-0.9738062,-0.9738062,-0.9738062,
c     1          -0.32035008744/

                  
	data pi/3.14159265358979d0/
cgcg	call etime(time)
cgcg	write(6,*) 'time',time
 
	open(10, file='finalJ.bin', status='old')
	open(11, file='initJ.bin', status='old')
	read(11,*) jint
	read(10,*) jfin
	close(10)
	close(11)
	ident_flag='Ji'//char(jint+48)//'Jf'//char(jfin+48)
	open(35,file="Xsect.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(65, file="Xsect.nu.OnePhoton"//ident_flag//".dat",
     1 status="unknown")
	open(135,file="PartialsL.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(155,file="PartialsL.nu.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(136,file="PartialsV.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(156,file="PartialsV.nu.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(37,file="Nopen.OnePhoton"//ident_flag//".Dat",
     1 status="unknown")
	open(20,file='Chan_Elimination'//ident_flag//'.dat',
     1 status="unknown")


	kchn=nopmax
	ndim=nop2
 
	call setup
	open(15, file="Einputs.dat", status='old')
	read(15,*) e0,ef,ne,nchanNEW
	close(15)
      nooNEW=nchanNEW*2
c---Energy bounds set at the experimental p12 threshold
c   and study the photoionization cross section between
c   the 65th and 180th state, we should see the series
c   converging. 
	e6p12 = -0.3242741d0
c   Alternatively we can do the threshold obtained from
c   the potentials.
c	e6p12 = -3.240880722315326D-01
chg	e0 = e6p12-0.5d0/65.d0**2
chg	ef = e6p12-0.5d0/85.d0**2
c	e0+3.05040d-03
chg	ne = 1200
 
	open(unit=18,file='engy18.dat',status='old')
	rewind(18)
	read(18,*) zz,r0,nlj
	read(18,35) (en(i),dy(i),y(i),nodes(i),i=1,norbs)
	do 146 i=1,norbs
	do 146 ii=1,norbs
	   ov(i,ii)=0.d0
  146   continue
	iend=0
	do 151 kk=1,nlj
	   ibegin=iend+1
	   iend=iend+nr
	   read(18,4) ((ov(i,ii),ii=ibegin,i),i=ibegin,iend)
  151   continue
	   read(18,4) ((rm(i,ii),ii=1,i),i=1,norbs)
	   read(18,4) ((rin(i,ii),ii=1,i),i=1,norbs)
	   read(18,4) ((dr(i,ii),ii=1,norbs),i=1,norbs)
	do 155 i=1,norbs-1
	do 155 ii=i+1,norbs
	      rm(i,ii)=rm(ii,i)
	      rin(i,ii)=rin(ii,i)
	      ov(i,ii)=ov(ii,i)
  155   continue
 
cc OLDformat-->   35   format(1p3d18.10)
  35    format(1p3d22.15,i5)

	close(18)
	write(6,*) 'unit 18 input done ',en(1)
 
c -- Initial state data


	  open(unit=16,file='initvec.dat',status='unknown',
     1		access='sequential',form='unformatted')
	  rewind(16)
	  read(16) ipar0,j0t,n0,nbl0,NumEig0,itake
        write(6,*) 'ipar0,j0t,n0,nbl0,NumEig0,itake '
        write(6,*) ipar0,j0t,n0,nbl0,NumEig0,itake
        IF(NumEig0 .GT. nEig0) THEN
            WRITE(6,*) '**** ->In parameter statement'
            WRITE(6,*) 'increase nEig0 to:',NumEig0,' from ',nEig0
            STOP
        ENDIF
        IF(n0 .GT. ntwo) THEN
            WRITE(6,*) '**** ->In parameter statement'
            WRITE(6,*) 'increase ntwo to:',n0,' from ',ntwo
            STOP
        ENDIF
	  read(16) (n10(kk),n20(kk),kk=1,n0)
        write(6,*) 'n10,n20 ',n10(n0),n20(n0)
	  read(16) (l1b0(i),l2b0(i),j1b0(i),j2b0(i),i=1,nbl0)
        write(6,*) 'l1b0, l2b0,j1b0,j2b0 ',j2b0(nbl0)
	  read(16) (kkst0(i),kkend0(i),i=1,nbl0)
        write(6,*) 'kkst0,kkend0 ',kkst0(nbl0),kkend0(nbl0)
        do ia=1,NumEig0
	    read(16) eval0(ia)
          write(6,*) ia, eval0(ia)
	    read(16) (vec00(i,ia),i=1,n0)
          write(6,*) vec00(1,ia),vec00(2,ia),vec00(3,ia)
        enddo
	  close(16)

c -- The following (itake) specifies which initial eigenstate to use in the 
c      photoabsorption calculation:

      egnd=eval0(itake)
      write(6,*) 'itake, egnd ',itake, egnd
      do i=1,n0
        vec0(i)=vec00(i,itake)
      enddo
        write(6,*) itake,egnd,(vec0(i),i=1,10)
 
	j0=j0t*2
 
c -- Final state data
 
	open(unit=15,file='finvec.dat',status='old',
     1		access='sequential',form='unformatted')
	rewind(15)
	read(15) iparf,jt,ncc,noo,nchan,nblock,ntot
	write(6,*) iparf,jt,ncc,noo,nchan,nblock,ntot
	write(6,*) 'ntot,noo,nchan',ntot,noo,nchan
	read(15) (n1f(k),n2f(k),k=1,ntot)
	read(15) (l1bf(k),l2bf(k),j1bf(k),j2bf(k),k=1,nblock)
	read(15) (kkstf(k),kkendf(k),k=1,nblock)
	read(15) (elam(ia),ia=1,ncc)
	read(15) ((x(i,ia),i=1,ncc),ia=1,ncc)
cc OLDformat-->    4   format(1p5d16.9)
    4   format(1p3d22.15)
	read(15) ((hlco(ia,i),hloc(i,ia),ia=1,ncc),i=1,noo)
	read(15) ((olco(ia,i),oloc(i,ia),ia=1,ncc),i=1,noo)
	read(15) ((hloo(i,j),j=1,noo),i=1,noo)
	read(15) ((oloo(i,j),j=1,noo),i=1,noo)
	close(15)
      write(6,*) 'final state input finished ',iparf, jt, ncc, noo
      write(6,*) nchan, nblock, ntot
	jd=jt*2
 
c         write(6,*) 'elam='
c         write(6,*) (elam(ia),elam(ncc-ia+1),ia=1,ncc/2)
c         write(6,*) 'x vector='
c         write(6,*) (x(i,1),i=1,10)
c         write(6,*) (x(i,2),i=1,10)
c         write(6,*) (x(i,ncc),i=1,10)
 
c ***********************************************
 
c -- Calculate dipole matrix elements in the primitive basis
 
cgcg	call etime(time)
    	write(6,*) 'input done' ,nbl0, nblock
	do 3333 ibl0=1,nbl0
	do 3333 iblf=1,nblock
	   ll10=l1b0(ibl0)
	   ll20=l2b0(ibl0)
	   jj10=j1b0(ibl0)
	   jj20=j2b0(ibl0)
 
	   ll1f=l1bf(iblf)
	   ll2f=l2bf(iblf)
	   jj1f=j1bf(iblf)
	   jj2f=j2bf(iblf)
	   lmin0=iabs(ll10-ll20)
	   lmax0=ll10+ll20
	   lminf=iabs(ll1f-ll2f)
	   lmaxf=ll1f+ll2f
	   iterm=0
	   do 1500 is=0,2,2
	   do 1450 ll0=lmin0,lmax0,2
	   do 1400 llf=lminf,lmaxf,2
	      if(iabs(ll0-llf) .gt. 2) go to 1400
	      iterm=iterm+1
	      front0=(is+1)*(ll0+1)*(jj10+1)*(jj20+1)
	      recoup0=dsqrt(front0)*xninej(1,ll10,jj10,1,ll20,
     1                            jj20,is,ll0,j0)
	      front1=(is+1)*(llf+1)*(jj1f+1)*(jj2f+1)
	      recoup1=dsqrt(front1)*xninej(1,ll1f,jj1f,1,ll2f,
     1                            jj2f,is,llf,jd)
       front2=(-1)**((2+is+ll0+jd)/2)*dsqrt((llf+1.d0)*(j0+1.d0))
	      recoup2=front2*sixj(jd,is,llf,ll0,2,j0)
c!djh            write(549,*) iterm,is,ll0,llf,recoup0,recoup1,recoup2
	      xrecp(iterm)=recoup0*recoup1*recoup2
	      la=ll10/2
	      lb=ll20/2
	      lc=ll1f/2
	      ld=ll2f/2
	      l0=ll0/2
	      lf=llf/2
	      ip0=(-1)**(la+lb+l0+(is/2))
	      ipf=(-1)**(lc+ld+lf+(is/2))
	      cc(1,iterm)=tm1(la,lb,l0,lc,ld,lf)
	      cc(2,iterm)=tm1(lb,la,l0,ld,lc,lf)*ip0*ipf
	      cc(3,iterm)=tm1(lb,la,l0,lc,ld,lf)*ip0
	      cc(4,iterm)=tm1(la,lb,l0,ld,lc,lf)*ipf
	dd(1,iterm)=ang(la,lb,l0,lc,ld,lf)
	dd(2,iterm)=ang(lb,la,l0,ld,lc,lf)*ip0*ipf
	dd(3,iterm)=ang(lb,la,l0,lc,ld,lf)*ip0
	dd(4,iterm)=ang(la,lb,l0,ld,lc,lf)*ipf

c!djh      write(644,*)iterm,dd(1,iterm),dd(2,iterm),dd(3,iterm),dd(4,iterm)

 1400     continue
 1450     continue
 1500     continue
 
	iip0=(-1)**iabs((j0-jj10-jj20+2)/2)
	iipf=(-1)**iabs((jd-jj1f-jj2f+2)/2)
	cc1=angjj(ll1f,jj1f,ll2f,jj2f,jd,ll10,jj10,ll20,jj20,j0)
	cc2=angjj(ll2f,jj2f,ll1f,jj1f,jd,ll20,jj20,ll10,jj10,j0)
	cc2=cc2*iip0*iipf
	cc3=angjj(ll1f,jj1f,ll2f,jj2f,jd,ll20,jj20,ll10,jj10,j0)
	cc3=cc3*iip0
	cc4=angjj(ll2f,jj2f,ll1f,jj1f,jd,ll10,jj10,ll20,jj20,j0)
	cc4=cc4*iipf

	dd1=cc1*ratiov(ll1f,ll10)
	dd2=cc2*ratiov(ll2f,ll20)
	dd3=cc3*ratiov(ll1f,ll20)
	dd4=cc4*ratiov(ll2f,ll10)
 
	      do 2500 k0=kkst0(ibl0),kkend0(ibl0)
	      do 2499 kf=kkstf(iblf),kkendf(iblf)
		 ka=n10(k0)
		 kb=n20(k0)
		 kc=n1f(kf)
		 kd=n2f(kf)
 
c -- New jj-coupling length calculation:

	suml=rm(ka,kc)*ov(kb,kd)*cc1+rm(kb,kd)*ov(ka,kc)*cc2
	suml=suml+rm(kb,kc)*ov(ka,kd)*cc3+rm(kd,ka)*ov(kc,kb)*cc4

c -- New jj-coupling velocity calculation:

	sumvv=ov(kb,kd)*(dr(kc,ka)-(ll10/2)*rin(kc,ka))*cc1
	sumvv=sumvv+ov(kb,kd)*rin(kc,ka)*dd1*ll10/2.d0
	add=ov(ka,kc)*(dr(kd,kb)-(ll20/2)*rin(kd,kb))*cc2
	add=add + ov(ka,kc)*rin(kd,kb)*dd2*ll20/2.d0
	sumvv=sumvv+add
	add=ov(ka,kd)*(dr(kc,kb)-(ll20/2)*rin(kc,kb))*cc3
	add=add+ ov(ka,kd)*rin(kc,kb)*dd3*ll20/2.d0
	sumvv=sumvv+add
	add=ov(kc,kb)*(dr(kd,ka)-(ll10/2)*rin(kd,ka))*cc4
	add=add+ ov(kc,kb)*rin(kd,ka)*dd4*ll10/2.d0
	sumvv=sumvv+add

c -- The following is the old (incorrect length+velocity calculation)

 
	    xn=delt(ka,kb)
	    xn=1.d0/dsqrt(1.d0+xn)
	    xnf=delt(kc,kd)
	    xnf=1.d0/dsqrt(1.d0+xnf)
	    dlen(k0,kf)=xn*xnf*suml
	    dvel(k0,kf)= xn*xnf*sumvv

c!djh          write(645,*) k0,kf, dlen(k0,kf), dvel(k0,kf)

2499	continue
 2500       continue
 3333   continue
 
cgcg	call etime(time)
	write(6,*) 'primitive dipole calc done',time

c -- Sum primitive dipole matrix elements over initial eigenvector
c ------- and transform into eigenrepresentation

	do 4700 ia=1,ncc
	   ddlen(ia)=0.d0
	   ddvel(ia)=0.d0
	   do 4600 i=1,ncc
	   do 4500 k0=1,n0
		prod=vec0(k0)*x(i,ia)*dlen(k0,i)
		ddlen(ia)=ddlen(ia)+prod
		prodv=vec0(k0)*x(i,ia)*dvel(k0,i)
		ddvel(ia)=ddvel(ia)+prodv
 4500	   continue
 4600	   continue
 4700	continue
	kf=ncc
	do 4900 ii=1,noo
	   kf=kf+1
	   ddlen(kf)=0.d0
	   ddvel(kf)=0.d0
	   do 4850 k0=1,n0
	      ddlen(kf)=ddlen(kf)+dlen(k0,kf)*vec0(k0)
	      ddvel(kf)=ddvel(kf)+dvel(k0,kf)*vec0(k0)
 4850	   continue
 4900   continue

         nchan=nchanNEW
         noo=nooNEW
         ntot = ncc + noo
c -- Begin main energy loop
	xne=ne
	if(ne.eq.1) xne=0.d0
	he = (ef-e0)/(xne-1.d0)
	write(6,*) xne,e0,ef,he,'energy loop begins'
	open(unit=27, file="time_dels_jjstream_jf"
     1   //char(jt+48)//".dat",status="unknown")

	open(unit=28, file="kmatrix_jjstream_jf"
     1   //char(jt+48)//".dat", status="unknown")
 6038 format(1000(3x,e20.8))
 6138 format(1000(3x,e20.8))
	rewind(27)
	do 9999 ie=1,ne
	   e = e0 + (ie-1)*he
chg         xnuxnu = xnu00+(ie-1)*hnu
chg         e = ethresh-0.5d0/xnuxnu**2

	   irat= noo/nchan
c	   write(6,*) 'irat=',irat,noo,nchan,e
	   do 200 k1=1,noo
	   do 200 k2=1,noo
	      sum1=2.d0*(e*oloo(k1,k2)-hloo(k1,k2))
	      sum2=0.d0
	      do 150 ia=1,ncc
		prod1=e*oloc(k1,ia)-hloc(k1,ia)
		prod2=e*olco(ia,k2)-hlco(ia,k2)
		sum2=sum2 - 2.d0 * prod1*prod2/(e-elam(ia))
  150         continue
	      omega(k1,k2)=sum1 + sum2
	      i=k1 + ncc
	      j=k2 + ncc
	      nn2 = n2f(i)
	      nn2p= n2f(j)
	      ic= (k1+1)/irat
	      icp= (k2+1)/irat
            dlt=delt(ic,icp)
	      qlam(k1,k2)=delt(ic,icp)*y(nn2)*y(nn2p)
  200      continue
	   asym=0.d0
	   do 250 k1=1,noo
	   do 250 k2=1,k1
	      err=dabs(omega(k1,k2)-omega(k2,k1))
	      if(err.gt.asym) asym=err
		ave=(omega(k1,k2)+omega(k2,k1))/2.d0
		omega(k1,k2)=ave
		omega(k2,k1)=ave
  250      continue
 
c	call etime(time)
c	write(6,*) 'R-matrix calc',time
		write(6,*) 'dim and oo',ndim, noo
	   call rgg(ndim,noo,omega,qlam,alfr,alfi,bvec,1,vec,ierr)
c	call etime(time)
	write(6,*) 'rgg done w/ ierr',ierr
 
	   ib=0
	   do 500 kb=1,noo
c		  write(6,*) "GO kb ", kb, noo, bvec(kb)
	      if(dabs(bvec(kb)) .lt. 1.d-10) go to 500
	      ib=ib+1
	      bvar(ib)=alfr(kb)/bvec(kb)
	      do 310 i=1,noo
		 vec1(i,ib)=vec(i,kb)
  310         continue
 
c *** Construct the closed portion of the eigenvector now:
 
	      do 340 ia=1,ncc
		zvec(ia,ib)=0.d0
	      do 335 k=1,noo
		 prod1=(e*olco(ia,k)-hlco(ia,k))*vec1(k,ib)
		 zvec(ia,ib)=zvec(ia,ib)- prod1/(e-elam(ia))
  335         continue
  340         continue
	   k=ncc
	   do 410 i=1,noo
	      k=k+1
	      zvec(k,ib)=vec1(i,ib)
  410      continue
 
c	call etime(time)
c	write(6,*) 'dbeta calc starts',time
c	write(6,*) 'ib,n0,ntot,noo,nchan',ib,n0,ntot,noo,nchan
	   dbetal(ib)=0.d0
	   dbetav(ib)=0.d0
	   do 450 kf=1,ntot
	      term=ddlen(kf)*zvec(kf,ib)
	      dbetal(ib)=dbetal(ib) + term
	      termv=ddvel(kf)*zvec(kf,ib)
	      dbetav(ib)=dbetav(ib) + termv
  450	   continue
 
c!djh	write(546,*) 'dbetal,dbetav',dbetal(ib),dbetav(ib)

  500      continue
	   write(6,*) 'ib ',ib,' nchan ', nchan
	   if(ib.ne.nchan) stop 'error in number of solutions'
 
c	call etime(time)
c	write(6,*) 'dbeta calc finished',time
c *** Obtain wavefunctions and derivatives on the boundary
 
	   irat=noo/nchan
	   do 600 ib=1,nchan
	   do 600 i=1,nchan
	      ps(i,ib)=0.d0
	      do 580 j=1,irat
	       xx=y(n2f(ncc+irat*(i-1)+j))*vec1(irat*(i-1)+j,ib)
	       ps(i,ib)=ps(i,ib)+xx
  580         continue
	      dps(i,ib)= -ps(i,ib)*bvar(ib)
  600      continue
 
c -- Check orthogonality of ps(i,ib)
	do 617 ib=1,nchan
	do 617 ibp=1,ib-1
	   test=0.d0
	   do 616 i=1,nchan
		test=test+ps(i,ib)*ps(i,ibp)
  616	   continue
           if(dabs(test).gt.1.d-08) then 
		write(6,*) 'ib,ibp,test',ib,ibp,test
	   endif
  617	continue
c *** Obtain I and J matrices by matching to Coulomb functions
 
c	write(6,*) 'orthogonality test done'
 
	   ic=0
	   do 620 iblk=1,nblock
	      isum=0
	      do 618 k=kkstf(iblk),kkendf(iblk)
		 isum=isum+1
		 if(k.le.ncc .or. isum.ne.1 .OR. k.gt.ntot) go to 618

		 ic=ic+1
		    
c! djh  how about here, this way, instead of chg below
		if((expth_flag).and. (ll.le.num_expth)) then
			kk = mod(n1f(k)-1,nr)+1
				ll = (n1f(k)-kk)/nr + 1
				eps(ic) = e - expth(ll)

cc!		    	if (ie.eq.1) then
					print *, "For final J ", jt
					print *, "   Channel ", ic, " reset:"
					print *, "     old: ", en(n1f(k))
					print *, "     now: ", expth(ll)
cc!		    	endif
		else
						
			eps(ic) = e - en(n1f(k))
			
		endif
			
		 
		 epsexp(ic)=eps(ic)

chg
chg  -- The following commented-out option is occasionally used when
chg       it is desirable to use experimental thresholds, in which case
chg       they need to be entered in a data statement at the beginning
chg       of the program:
chg   		    epsexp(ic) = e - expth(ic)
chg
		    lchan(ic)=l2bf(iblk)/2
			eryd=2.d0*eps(ic)
	      call seaton(lchan(ic),eryd,r0,1.d0,f(ic),
     o             df(ic),g(ic),dg(ic))
c	write(6,*) 'seaton=',ic,f(ic)
c	write(6,*) 'eryd,lchan(ic)=',eryd,lchan(ic)
	xnu=1.d0/dsqrt(dabs(eryd))
	beta(ic)=pi*(xnu-lchan(ic))
	if(eryd.lt. -1.d0/(lchan(ic)**2+1.d-10)) then
		f(ic)=f0
		df(ic)=fp0
		g(ic)=g0
		dg(ic)=gp0
		xxxx=aqdt/dtan(pi*(xnu-lchan(ic))) + gqdt
		beta(ic)=datan(1.d0/xxxx)
	endif
	wron=f(ic)*dg(ic)-df(ic)*g(ic) - 2.d0/pi
c	write(6,*) 'wronskian=',wron
	if(dabs(wron).gt.1.d-04) write(6,*) ic,wron
3123	format(5(1x,1pd11.4))
  618         continue
  620      continue
 
 
	   if(ic.ne.nchan) THEN 
               write(6,*)   'error in block indexing'
               write(6,*) 'ic,isum,nchan=',ic,isum,nchan
           STOP
         endif
 
	   pi2=pi/2.d0
	   do 650 ib=1,nchan
	   do 650 i=1,nchan
	      qi(i,ib)= (ps(i,ib)*dg(i)-dps(i,ib)*g(i))*pi2
	      qj(i,ib)= (ps(i,ib)*df(i)-dps(i,ib)*f(i))*pi2
ccc		write(67,*) i,ib,ps(i,ib),dps(i,ib)
		omega(i,ib)=qi(i,ib)
  650      continue
 
c **   Invert qi
c		write(150,*) e,nchan,'qi-matrix before inversion'
c		write(150,6644) ((qi(i,ib),ib=1,nchan),i=1,nchan)
 6644		format(5(1x,1pd14.7))
c		write(151,*) e,nchan,'qj-matrix'
c		write(151,6644) ((qj(i,ib),ib=1,nchan),i=1,nchan)
	
ccc	   call gmis(nchan,qi,cm,dm,det,0,kchn)
c   invert matrix qi to compute the k matrix.
	   call inv(qi,nchan,ipivot,ind1,ind2,kchn)
c	write(6,*) 'inv finished ',qi(1,1)
 
c	write(152,*) e,nchan,'qi after inversion'
c	write(152,6644) ((qi(i,ib),ib=1,nchan),i=1,nchan)

c		write(153,*) e,nchan,'f'
c		write(153,6644) (f(ic),ic=1,nchan)
c		write(154,*) e,nchan,'g'
c		write(154,6644) (g(ic),ic=1,nchan)
c -- Test inverse
	do 670 i=1,nchan
	do 670 j=1,nchan
		test=0.d0
		do 665 k=1,nchan
		test=test+omega(i,k)*qi(k,j)
  665		continue
  		if(i.ne.j .and. dabs(test).gt.1.d-08) then
	write(6,*) 'i,j,inverse test=',i,j,test,omega(i,j)
		endif
  670	continue

c **   Get K-matrix
 
	   do 700 i=1,nchan
	   do 700 j=1,nchan
	      xk(i,j)=0.d0
	      do 680 k=1,nchan
		 xk(i,j)=xk(i,j)+qj(i,k)*qi(k,j)
  680         continue
c   	write(66,*) i,j,xk(i,j)
  700      continue
ccc	write(54,3123) ((xk(i,j),j=1,nchan),i=1,nchan)
	asym=0.d0
 	do 702 i=1,nchan
	do 702 j=1,i
		test=dabs(xk(i,j)-xk(j,i))
		if(test.gt.asym) asym=test
		xk(i,j)=(xk(i,j)+xk(j,i))/2.d0
		xk(j,i)=xk(i,j)
  702	continue
	if(dabs(asym).gt.1.d-06) then
	write(6,*) 'maximum asymmetry in K=',asym,e
	endif
c testing K-Matrix values

	open(unit=333, file='kmatrix.dat')
	do i=1,nchan
	   do j=1,nchan
	      write(333,103) xk(i,j), e
 103	      format(2d20.8)
	   end do
	end do
	close(333)
c   writing K matrix values.
	write(6, *)"Saving k mat for j ", jt
	write(28,6038) e, ((xk(i,j), j=i, nchan), i=1, nchan)
c **   Eigenchannel parameters
 
	   call vdiag(xk,xmu,umat,nchan,ierr,bvec,kchn)
	   do 720 ia=1,nchan
	      xmu(ia)=datan(xmu(ia))/pi
 
	   dalphl(ia)=0.d0
	   dalphv(ia)=0.d0
	   do 835 ib=1,nchan
	   do 835 i=1,nchan
	      term1=qi(ib,i)*umat(i,ia)*dcos(pi*xmu(ia))
	      dalphl(ia)=dalphl(ia) + dbetal(ib)*term1
	      dalphv(ia)=dalphv(ia) + dbetav(ib)*term1
  835      continue
ccc	write(6,*) 'dalphl,dalphv',dalphl(ia),dalphv(ia)
  720      continue
 
ccc	write(55,*) e
ccc	write(55,3123) (xmu(ia),ia=1,nchan)
 
 
c  ** MQDT calculation
chg -- Note that the lowest channel will be treated as open even if closed:
      IF(DABS(epsexp(2)-epsexp(1)).LT.1.d-06) epsexp(2)=10.d0
	epsexp(1)=10.d0
ccc      epsexp(2)=10.d0
chg
           eOpenCutoff=0.d0!-0.0325d0

	    iopen=0
	   do 800 i=1,nchan
	      xnu=0.d0
	      if(epsexp(i).ge.eOpenCutoff) iopen=iopen+1
            if(epsexp(i).lt.eOpenCutoff) xnu=beta(i)/pi
	   do 800 ia=1,nchan
		 xgam(i,ia)=umat(i,ia)*dsin(pi*(xnu+xmu(ia)))
	      if(epsexp(i).lt.eOpenCutoff) then
		 xlam(i,ia)=0.d0
	      else
		 xlam(i,ia)=umat(i,ia)*dcos(pi*xmu(ia))
	      endif
  800      continue
	   nopen=iopen
c	call etime(time)
c	write(6,*) 'MQDT calculation',time
 
	   call rgg(kchn,nchan,xgam,xlam,alfr,alfi,bvec,1,avec,ierr)
 
	   suml=0.d0
	   sumv=0.d0
	   irho=0
	   do 850 kk=1,nchan
	   if(dabs(bvec(kk)).lt. 1.d-08) go to 850
	      irho=irho+1
	      tao(irho)=datan(alfr(kk)/bvec(kk))/pi
	      xnorm(irho)=0.d0
	      do 815 i=1,nchan
	      if(epsexp(i).lt.eOpenCutoff) go to 815
	      xnm=0.d0
	      do 810 ia=1,nchan
	p=umat(i,ia)*avec(ia,kk)*dcos(pi*(-tao(irho)+xmu(ia)))
		 xnm=xnm+p
  810         continue
	      xnorm(irho)=xnorm(irho)+xnm*xnm
                Tmat(i,irho)=xnm
  815		continue

	      drhol(irho)=0.d0
	      drhov(irho)=0.d0
 
	      do 825 ia=1,nchan
		 drhol(irho)=drhol(irho)+dalphl(ia)*avec(ia,kk)
		 drhov(irho)=drhov(irho)+dalphv(ia)*avec(ia,kk)
  825      continue
 
	   suml=suml+drhol(irho)**2 / xnorm(irho)
	   sumv=sumv+drhov(irho)**2 / xnorm(irho)



ccc	write(6,*) 'drhol,drhov',irho,drhol(irho),drhov(irho) 
  850   continue
	if(nopen.ne.irho) then
		write(6,8892) (bvec(i),i=1,nchan)
 8892	format(5(1x,1pd12.5))
	endif


              WRITE(19,*) e,Nopen
              do 8815 i=1,nchan
                if(epsexp(i).lt.eOpenCutoff) go to 8815
              do 8814 irho=1,Nopen
                Tmat(i,irho)=Tmat(i,irho)/dsqrt(xnorm(irho))
 8814	      continue
              WRITE(19,*) (Tmat(i,iirho),iirho=1,Nopen)
 8815	      continue
              WRITE(41,*) e,tao(1)
              WRITE(42,*) e, tao(2)
              IF(ie. EQ. 1) THEN
                xnu=1.d0/dsqrt(-2.d0*epsexp(1))
                faclast=dsin(pi*(tao(1)+xnu))*dsin(pi*(tao(2)+xnu))
                xnulast=xnu
              ELSE
              IF(epsexp(1).LT. 0.d0) THEN
               xnu=1.d0/dsqrt(-2.d0*epsexp(1))
               fac1= dsin(pi*(tao(1)+xnu))*dsin(pi*(tao(2)+xnu))
               firstfac=dsin(pi*(tao(1)+xnu))
               secondfac=dsin(pi*(tao(2)+xnu))
               testprod=fac1*faclast
               abstest=dabs(testprod)
               WRITE(401,*) xnu,fac1
               WRITE(402,*) xnu,faclast
               WRITE(403,*) xnu,testprod
               IF(testprod .LT.0.d0 .AND. abstest.LT.1.d-03) THEN
                  write(400,*) xnu,testprod
                  xnuroot= (xnu*facLast-xnuLast*fac1)/(faclast-fac1)
                  write(500,*) xnuroot, firstfac, secondfac
               ENDIF
               WRITE(40,*) xnu,fac1
              ENDIF
               faclast=fac1
               xnulast=xnu
             ENDIF


	sumv=sumv/(e-egnd)**2
	omau= e-egnd
	const=4.d0*pi*pi*omau/(3.d0*(j0+1.d0))
	const=const*5.29d0**2 /137.036d0
	omcm=omau*2.d0*109737.3d0
	wlam=1.d+07 / omcm
c
c -- to refer energies to the ground state energy, comment out the following
c
c        expgnd=egnd
	eFinalCm = omcm + (egnd-expgnd)*2.d0*109737.3d0

	suml=suml*const
	sumv=sumv*const

c!djh    so suml and sumv are xsects in megabarnes


c!djh 		write(35,*) eFinalCm,suml,sumv
 		    write(35,*) (e-en(1))*2.1947e+5,suml,sumv
			write(65,*) dsqrt(-0.5/(e-e6p12)),suml,sumv
            write(31,*) e, suml
            write(32,*) e, sumv



c*** new section with partial cross sections

	do 9010 irh=1,nopen
		drhol(irh)=drhol(irh)/dsqrt(xnorm(irh))
		drhov(irh)=drhov(irh)/dsqrt(xnorm(irh))
		drhov(irh)=drhov(irh)/(e-egnd)
9010	continue

          write(23,*) e,(drhol(irh),irh=1,nopen)
          write(24,*) e,(drhov(irh),irh=1,nopen)
          write(43,*) e,drhol(1)
          write(44,*) e,drhol(2)
          write(53,*) e,drhol(1)*drhol(1)*2*(e-egnd)/3.d0
          write(54,*) e,drhol(2)*drhol(2)*2*(e-egnd)/3.d0
          write(55,*) e,drhov(1)*drhov(1)*2*(e-egnd)/3.d0
          write(56,*) e,drhov(2)*drhov(2)*2*(e-egnd)/3.d0

	ii=0	

        do 9040 i=1,nchan
         if(epsexp(i).lt.eOpenCutoff) go to 9040
	 
	 ii=ii+1

	partialL(ii)=0.d0
        partialV(ii)=0.d0

        do 9035 irho=1,nopen
        do 9034 irhop=1,nopen
         coscos=dcos(pi*(tao(irho)-tao(irhop)))
         prodL=Tmat(i,irho)*Tmat(i,irhop)*coscos
         prodL=prodL*drhol(irho)*drhol(irhop)
         prodV=Tmat(i,irho)*Tmat(i,irhop)*coscos
         prodV=prodV*drhov(irho)*drhov(irhop)
         partialL(ii)=partialL(ii)+prodL
         partialV(ii)=partialV(ii)+prodV
 9034	 continue
 9035	 continue
 9040	 continue


	 if (ii/=nopen) then
	    print *, "NOPEN ERR  ", ii,nopen,e,ie
	 else
	    write(135,185) (e-e0)*2.1947e+5,partiall(1:nopen)
		write(155,185) dsqrt(-0.5/(e-e6p12)), partiall(1:nopen)
	    write(136,185) (e-e0)*2.1947e+5,partialv(1:nopen)
		write(156,185) dsqrt(-0.5/(e-e6p12)), partialv(1:nopen)
	 endif
  185 format((f18.12, 100e20.10))



  45    format(2(1x,f40.18))


c MAAB rvise this calculation of the time delays, where is the
c S matrix calculated?
c where is the derivative to get Q?	
	write(37,*) eFinalCm,nopen
	write(20,8895) e, (tao(irh),irh=1,nopen)
 8895 format(25(1x,e18.12))
         sumtao=0.d0
         DO irh=1,nopen
         sumtao=sumtao+tao(irh)
         ENDDO
		 IF(e.eq.e0) THEN
			npis= 0
		 ELSE
			IF(dabs(sumtao-sumtaoprev).gt.0.5) THEN
				npis = npis+1
			ENDIF
		 ENDIF
		 sumtaoprev = sumtao
         WRITE(26,*) e,sumtao
		 WRITE(27,*) e, sumtao+npis


      IF(nopen.EQ.1) THEN
        write(38,91) eps(1),tao(1)
        write(39,91) beta(1)/pi,dsin(tao(1)*pi+beta(1))
  91    FORMAT(2(1x,f30.15))
      ENDIF
	


 9999   continue
	end
ccc
	double precision function tm1(l1,l2,l,l1p,l2p,lp)
	implicit double precision (a-h,o-z)
	tm1=0.d0
	if(l2.ne.l2p) return
	l1d=l1*2
	l2d=l2*2
	ld=l*2
cgcg	md=m*2
	l1pd=l1p*2
	l2pd=l2p*2
	lpd=lp*2
cgcg	mpd=mp*2
cgcg	iqd=iq*2
	tm1=(-1)**(   l1+l2+lp+1+l1)*dsqrt((2.d0*l+1)*(2*lp+1)
     1      *(2*l1+1)*(2*l1p+1))
     1      *sixj(l1d,ld,l2d,lpd,l1pd,2)
     1      *thrj(l1d,2,l1pd,0,0,0)*delt(l2,l2p)
	return
	end
******************************************************************************
	double precision function tm2(l1,l2,l,m,l1p,l2p,lp,mp,iq)
	implicit double precision (a-h,o-z)
	tm2=0.d0
	if(l1.ne.l1p) return
	l1d=l1*2
	l2d=l2*2
	ld=l*2
	md=m*2
	l1pd=l1p*2
	l2pd=l2p*2
	lpd=lp*2
	mpd=mp*2
	iqd=iq*2
	tm2=(-1)**(    l1+l2p+l+1+l2)*dsqrt((2.d0*l+1)*(2*lp+1)
     1      *(2*l2+1)*(2*l2p+1))
     1      *sixj(l2d,ld,l1d,lpd,l2pd,2)
     1      *thrj(l2d,2,l2pd,0,0,0)*delt(l1,l1p)
	return
	end
******************************************************************************
	double precision function ang(l1,l2,l,l1p,l2p,lp)
	implicit double precision(a-h,o-z)
	ang=0.d0
	if(l1.ne.l1p-1.and.l1.ne.l1p+1) return
	if(l2p.ne.l2) return
	if(l1p.eq.0) return
	l1d=2*l1
	l2d=2*l2
	ld=2*l
cgcg	md=2*m
	l1pd=2*l1p
	l2pd=2*l2p
	lpd=2*lp
cgcg	mpd=2*mp
cgcg	iqd=2*iq
 
 
	if(l1.eq.l1p-1)b=2.d0*l1p/(2*l1+1)-2.d0*l1p*l1p/(4*l1p*l1p-1)
	if(l1.eq.l1p+1) b=-2.d0*l1p*l1/(4.d0*l1*l1-1)
	a=dsqrt((2.d0*l1+1)*(2*l1p+1))*b/2.d0
	ang=(-1)**(    l1+l2+lp+1+l1)*dsqrt((2.d0*l+1)*(2*lp+1))
     1                              *sixj(l1d,ld,l2d,lpd,l1pd,2)
     1      /thrj(l1d,2,l1pd,0,0,0)*delt(l2,l2p)*a
	return
	end
*******************************************************************************
	double precision function r12(ii,jj)
	parameter (kbs=150,kbs0=200)
	implicit double precision(a-h,o-z)
	common /a/ndx(kbs0),rij(kbs0*(kbs0+1)/2)
	if(ii.gt.kbs.or.jj.gt.kbs) stop 'r12,i,j>kbs'
	i=ndx(ii)
	j=ndx(jj)
	if(i.ge.j) then
	k=(i*(i-1))/2+j
	else
	k=(j*(j-1))/2+i
	endif
	r12=rij(k)
	return
	end

	function angjj(l1d,j1d,l2d,j2d,jd,l10d,j10d,l20d,j20d,j0d)   
	implicit real*8(a-h,o-z)
	phase=(-1)**iabs((j0d+j2d-1-l10d+l1d)/2)
	front=(j1d+1)*(j10d+1)*(jd+1)*(j0d+1)*(l1d+1)*(l10d+1)
	s1=sixj(j1d,jd,j2d,j0d,j10d,2)
	s2=sixj(l1d,j1d,1,j10d,l10d,2)
	tj=thrj(l1d,2,l10d,0,0,0)
	angjj=phase*dsqrt(front)*s1*s2*tj
	return
	end
c
	function ratiov(ld,lpd)
	implicit real*8(a-h,o-z)
	ratiov=0.d0
	if( ld.eq. lpd-2) then

	tj=thrj(ld,2,lpd,0,0,0)
	term0=dsqrt((ld+1.d0)*(lpd+1.d0))*tj
	prod=dsqrt( (lpd+1.d0)/(lpd-1.d0) )
	   ratiov=prod/(tj*term0)
	endif
	return
	end

	function xnu(e,et)
		implicit real*8(a-h,o-z)
		xnu = SQRT(0.5/(et-e))
		return 
	end
