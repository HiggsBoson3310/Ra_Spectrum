    ! maab: Any changes in the number of partial waves and other parameters
!       have to reflect here.
!		ncr == number of closed type orbitals
!		nor == number of open type orbitals
!		nljmax == number of different partial waves orbitals if more orbitals are added.    
		parameter(ncr=28,nor=2,nr=ncr+nor,mm=4321,nblk=90)
        parameter(nopmax=40,nop2=nopmax*2,nls=40,nkval=60,ntwo=9000)
        parameter(nljmax=13,norbs=nljmax*nr)
	implicit real*8(a-h,o-z)
	dimension lb(norbs),jb(norbs),n1(ntwo),n2(ntwo)
	dimension l1(ntwo),l2(ntwo),j1(ntwo),j2(ntwo)
	dimension kkstart(nblk),kkend(nblk),k1o(nop2),k2o(nop2)
	dimension l1blk(nblk),l2blk(nblk),j1blk(nblk),j2blk(nblk)
	dimension xrecp(nls)
	dimension dirang(nls,nkval),excang(nls,nkval)
	dimension direct(nkval),exch(nkval),r12(ntwo,ntwo)
	dimension ov(norbs,norbs),rm(norbs,norbs),rin(norbs,norbs)
	dimension dr(norbs,norbs),en(norbs),y(norbs),dy(norbs)
	dimension subd(ntwo),eval(ntwo)
	dimension vec(ntwo,ntwo),olap(ntwo,ntwo)
	dimension hlco(ntwo,nop2),hloc(nop2,ntwo)
	dimension olco(ntwo,nop2),oloc(nop2,ntwo)
 
	one=1.d0
	open(unit=55, file='initJ.bin', status='old')
	read(55,*) jin
	read(55,*) iod
	close(55)
c ***** input open basis information

        open(unit=5,file='RaInitialStateJ'//char(jin+48)//'.dat'
     1       ,status='old')
	
	rewind(5)	
	read(5,*) ipar,jtrue,init,itake,Neig
c ****  ipar=+/- 1, jtrue is total J, init=1 for initial state only
	jd=jtrue*2
      nchan=0

      IF(init.NE.1) THEN
	read(5,*) nkopen,nchan
	do 110 i=1,nkopen
		read(5,*) k1o(i),k2o(i)
  110   continue
        close(5)
	nopper=nkopen/nchan
      ENDIF
 
CCC
	open(unit=8,file='config0.dat',status='unknown')
	read(8,181) ncc,noo,iblocks,nnljmax,nnr,ipar,j
	write(6,181) ncc,noo,iblocks,nnljmax,nnr,ipar,j
  181	format(12i5)
	nnorbs=nr*nljmax
      write(6,*) 'nnorbs = ',nnorbs
 
	ntot=ncc+noo
	read(8,181) (lb(k),jb(k),k=1,nnorbs)
      write(6,*) 'lb,jb ', jb(nnorbs)
	read(8,181) (n1(k),l1(k),j1(k),n2(k),l2(k),j2(k),k=1,ntot)
      write(6,*) 'n1 l1 j1 n2 l2 j2 ', j2(ntot)
	read(8,181) (l1blk(ib),j1blk(ib),ib=1,iblocks)
      write(6,*) 'l1blk j1blk iblocks ',j1blk(iblocks), iblocks
      read(8,181) (l2blk(ib),j2blk(ib),ib=1,iblocks)
      write(6,*) 'l2blk j2blk ', j2blk(iblocks)
	read(8,181) (kkstart(ib),kkend(ib),ib=1,iblocks)
      write(6,*) 'kkstart,kkend ',kkend(iblocks)
	close(8)
CCC


c	write(2,*) 'one-electron basis set key: ',norbs
c	write(2,4477) (lb(k),jb(k),k=1,norbs)
 4477   format(5(2i3,3x))
 
 
c	ncc=kk
 
c	write(2,*) ncc,' closed configurations'
c	do 4488 k=1,ncc
c	   write(2,*) k,': ',n1(k),l1(k),j1(k),'- ',n2(k),l2(k),j2(k)
c 4488   continue
 
	ntot=ncc + noo
 
c	write(2,*) noo,' open configurations'
c	do 4499 k=ncc+1,ntot
c	  write(2,*) k,': ',n1(k),l1(k),j1(k),'- ',n2(k),l2(k),j2(k)
c4499   continue
 
c	write(2,*) ' '
c	write(2,*) iblocks,' is number of distinct angular configs'
c	do 3377 ib=1,iblocks
c	   write(2,*) ib,l1blk(ib),j1blk(ib),'- ',l2blk(ib),j2blk(ib)
c	   write(2,*) '      kk-range',ib,kkstart(ib),kkend(ib)
c 3377   continue
c	close(2)
 
 
c ***** Input r12 integrals
 
	open(unit=12,file='init12.dat',status='old')
	rewind(12)
	read(12,*) ntot
	read(12,4) ((r12(k1,k2),k1=1,k2),k2=1,ntot)
	do 344 k2=1,ntot
	do 344 k1=1,k2
		r12(k2,k1)=r12(k1,k2)
  344   continue
cc OLDversion-->    4   format(5(1pd16.9))
   4    format(1p3d22.15)
c *************************
	write(6,*) 'r12 input ok',r12(1,1),r12(ntot,ntot)
	close(12)
 
 
c **************************************************
 
c **       UNIT 18 INPUT
 
	open(unit=18,file='engy18.dat',status='old')
	rewind(18)
	read(18,990) zz,r0,nl
	write(6,990) zz,r0,nl
 990    format(f4.0,f6.1,8i4)
        do 347 i=1,norbs
cc OLD	read(18,35) (en(i),dy(i),y(i),i=1,norbs)
	  read(18,35) en(i),dy(i),y(i)
 347    continue
cc        write(6,35)e0,derv,yr0,node
cc OLDversion-->  35    format(1p3d18.10)
  35    format(1p3d22.15,i5)

 
**** <i|j> matrix *****
	do 146 i=1,norbs
	do 146 ii=1,norbs
	ov(i,ii)=0.d0
  146   continue
	iend=0
	do 151 kk=1,nljmax
	ibegin=iend+1
	iend=iend+nr
	read(18,4) ((ov(i,ii),ii=ibegin,i),i=ibegin,iend)
 151    continue
 
**** <i|r|j> matrix(radial only) *******
	read(18,4) ((rm(i,ii),ii=1,i),i=1,norbs)
****<i|1/r|j> matrix (radial only)*******
	read(18,4) ((rin(i,ii),ii=1,i),i=1,norbs)
****<i|d/dr|j> matrix (radial only)*******
	read(18,4) ((dr(i,ii),ii=1,norbs),i=1,norbs)
 
	do 155 i=1,norbs-1
	do 155 ii=i+1,norbs
	  rm(i,ii)=rm(ii,i)
	  rin(i,ii)=rin(ii,i)
	  ov(i,ii)=ov(ii,i)
 155    continue
	close(18)
 
c ****************************************************************
 
	write(6,*) 'dr(norbs,norbs) ',dr(norbs,norbs)
	write(6,*) 'input finished'
 
c ***  SET up the hamiltonian matrix elements and diagonalize closed-part
 
	do 700 i=1,ntot
	   ll1=l1(i)
	   ll2=l2(i)
	   jj1=j1(i)
	   jj2=j2(i)
	   k1=n1(i)
	   k2=n2(i)
	   xn=delt(ll1,ll2)*delt(jj1,jj2)*delt(k1,k2)
	   xn=1.d0/dsqrt(1.d0 + xn)
 
	do 700 j=1,ntot
	   ll1p=l1(j)
	   ll2p=l2(j)
	   jj1p=j1(j)
	   jj2p=j2(j)
	   k1p=n1(j)
	   k2p=n2(j)
	   xnp=delt(ll1p,ll2p)*delt(jj1p,jj2p)*delt(k1p,k2p)
	   xnp=1.d0/dsqrt(1.d0+dsqrt(xnp))
	   phase=(-1)**( (jd+jj1p+jj2p)/2 )
	   angdir=delt(ll1,ll1p)*delt(jj1,jj1p)*delt(ll2,ll2p)
	   angdir=angdir*delt(jj2,jj2p)
	   angexc=delt(ll1,ll2p)*delt(ll2,ll1p)*delt(jj1,jj2p)
	   angexc=angexc*delt(jj2,jj1p)
 
	   olap(i,j)= angdir*ov(k1,k1p)*ov(k2,k2p)
	   olap(i,j)=olap(i,j) - phase*angexc*ov(k1,k2p)*ov(k2,k1p)
	   olap(i,j)=olap(i,j) * xn*xnp
ccc           write(6,*) i,j,olap(i,j)
 
	   r12(i,j)=r12(i,j) + olap(i,j)*(en(k1p)+en(k2p))
  700   continue
 
	nmax=ntwo
	write(6,*) 'diagonalization begins',ncc,nmax
	write(6,7703) ((r12(i,j),j=1,6),i=1,6)
 7703   format(6(1x,1pd12.5))
	nmax=ntwo
	call vdiag(r12,eval,vec,ncc,ierr,subd,nmax)


c        if (jtrue.eq.1) then
c
c           print *, "FIXING J=1 init state"
c           print *, "Calc - ", eval(2)
c           eval(2) = -0.47687010110d0
c           print *, "Expt - ", eval(2)

c        endif





	open(unit=3,file='eclosed.dat',status='unknown')
	rewind(3)
	write(3,7711) (eval(i),i=1,ncc)
	close(3)
 7711   format(5(1pd14.7))
	write(6,*) 'done ',eval(1),eval(2),eval(3)
 
	if(init.eq.1) then
	  ibclos=iblocks-nchan
	  open(unit=16,file='initvec.dat',status='unknown',
     1		access='sequential',form='unformatted')
	  rewind(16)
c
c -- Note that Neig, specified in the first line of input(unit 5) above, 
c      gives the number of initial states for which oscillator strengths
c      are desired.
c
	  write(16) ipar,jtrue,ncc,ibclos,Neig,itake
	  write(16) (n1(kk),n2(kk),kk=1,ncc)
	  write(16) (l1blk(i),l2blk(i),j1blk(i),j2blk(i),i=1,ibclos)
	  write(16) (kkstart(i),kkend(i),i=1,ibclos)
        do ia=1,Neig
	    write(16) eval(ia)
	    write(16) (vec(i,ia),i=1,ncc)
        enddo
	  close(16)
	  stop
	else


c -- if this is the final state, write out HL, olap, elam,X, etc
 
	open(unit=15,file='finvec.dat',status='unknown',
     1		access='sequential',form='unformatted')
	rewind(15)
	write(15) ipar,jtrue,ncc,noo,nchan,iblocks,ntot
	write(15) (n1(kk),n2(kk),kk=1,ntot)
       write(15) (l1blk(i),l2blk(i),j1blk(i),j2blk(i),i=1,iblocks)
	write(15) (kkstart(i),kkend(i),i=1,iblocks)
	write(15) (eval(ia),ia=1,ncc)
	write(15) ((vec(i,ia),i=1,ncc),ia=1,ncc)
 
c -- Next add Bloch Operator to relevant H matrix elements
c --------  and transform to the eigenrepresentation
 
	do 850 i=1,ntot
	do 850 j=1,ntot
	   bloch = 0.d0
	   if(i.le.ncc .and. j.le.ncc) go to 850
 
	   k1=n1(i)
	   k2=n2(i)
	   ll1=l1(i)
	   ll2=l2(i)
	   jj1=j1(i)
	   jj2=j2(i)
	   xn=delt(k1,k2)
	   xn=1.d0/dsqrt(1.d0+xn)
 
	   k1p=n1(j)
	   k2p=n2(j)
	   ll1p=l1(j)
	   ll2p=l2(j)
	   jj1p=j1(j)
	   jj2p=j2(j)
	   xnp=delt(k1p,k2p)
	   xnp=1.d0/dsqrt(1.d0+xnp)
 
	   term1=delt(ll1,ll1p)*delt(jj1,jj1p)*delt(ll2,ll2p)
	   term1=term1*delt(jj2,jj2p)*ov(k1,k1p)*y(k2)*dy(k2p)
	   term2=delt(ll1,ll2p)*delt(ll2,ll1p)*delt(jj1,jj2p)
	   term2=term2*delt(jj2,jj1p)*ov(k1,k2p)*y(k2)*dy(k1p)
	   phase=(-1)**( (jd+jj1+jj2)/2 )
	   bloch= 0.5d0*xn*xnp*(term1-phase*term2)
	   r12(i,j)=r12(i,j) + bloch
  850   continue
 
	diff=0.d0
	do 860 i=1,ntot
	do 860 j=1,i
	   err=dabs(r12(i,j)-r12(j,i))
	   if(err.gt.diff) diff=err
  860   continue
 
	write(6,*) 'maximum asymmetry in HL is ',diff
 
	do 1000 ia=1,ncc
	do 1000 i=1,noo
	   hlco(ia,i)=0.d0
	   hloc(i,ia)=0.d0
	   olco(ia,i)=0.d0
	   oloc(i,ia)=0.d0
	   do 910 j=1,ncc
	      hlco(ia,i)=hlco(ia,i)+vec(j,ia)*r12(j,i+ncc)
	      olco(ia,i)=olco(ia,i)+vec(j,ia)*olap(j,i+ncc)
	      hloc(i,ia)=hloc(i,ia)+vec(j,ia)*r12(i+ncc,j)
	      oloc(i,ia)=oloc(i,ia)+vec(j,ia)*olap(i+ncc,j)
  910      continue
 1000   continue
 
	write(15) ((hlco(ia,i),hloc(i,ia),ia=1,ncc),i=1,noo)
	write(15) ((olco(ia,i),oloc(i,ia),ia=1,ncc),i=1,noo)
	write(15) ((r12(i,j),j=ncc+1,ntot),i=ncc+1,ntot)
	write(15) ((olap(i,j),j=ncc+1,ntot),i=ncc+1,ntot)
 
	endif
 
	stop
	end
 
c
c	function delt(i,j)
c	implicit real*8(a-h,o-z)
c	delt=0.d0
c	if(i.eq.j) delt=1.d0
c	return
c	end
