	       PROGRAM  WAVEX
c   This program is the same as WAVE except that it uses different
c   radial mesh points. Here r=x**2 with x increasing linearly.
c   The potentials of the ions are read on unit 1,2,3 and 4.
c   This program generates numerical basis orbitals  in the ionic potl.
c
c   Link: gensub.f, gksub.f
c
c  nsc=no. of closed 's' orbitals ,npc= etc.
c  nso=no. of open   's' orbitals , npo =etc.
c
c   Parameters used in dimension statement:
c
c    kw=total number of one-electron wavefunctions
c    kbs=maximum number of two-electron base functions
c    ken=maximum number of energies
c    kchn=maximum number of channels
c    kms=maximum number of radial mesh
c
c   Variables:
c
c     nl=total mumber of different angular momenta
c     nwv(1,i)=number of wavefunctions (of l=lwv(i)) vanishing at r0
c     nwv(2,i)=number of wavefunctions (of l=lwv(i)) not vanishing at r0
c     en(i)=energy eigen values
c     y(i,j)=ith one-electron wavefunctions(multiplied by r) at jth radial point
c     yd(i,j)=derivative of y(i,j) with respect to x=sqrt(r)
c     rr(i)=radial distance=x(i)**2
c     vs(i,j)=ith potential(for the ith angular momentum) at rs(i,j)
c     zi(i)=ri(i)*vs(i)
c     eng(i)= estimated energy eigenvalue of the first wavefunction for ith l
c     ni, if r>ri(ni), use Coulomb potential instead of the numerical potential
c
 
 
 
 
	parameter(kw=390,kbs=150,ken=100,kchn=40,kms=4321,nlj=13)
	implicit real*8(a-h,o-z)
	dimension wk(kw,kw),en(25),yd(kw,kms),fdum(2,2)
	common/blok/y(kw,kms),ps(2,kms),
     &           rr(kms),qq(kms),x(kms),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid
	common/ncbcom/ncb(nlj)
cc -- old common block-> gave warnings!  common/sporb/ll,xjj,ijj,corr
	common/sporb/xjj,corr,ll,ijj
	common/rmc/vrmc
	external rwave
c -- The next call initializes the nuclear charge by calling fpot
	call fpot(0.1d0,fdum)
 
c   Read potential data
	call readin
      irt=0
c   Radial mesh
	xf=dsqrt(r0)
	xh=xf/m
	x0=xh
	do 1 i=1,m
	x(i)=xh*i
   1    rr(i)=x(i)*x(i)
	n=m
	ep=1.d-11
 
 
	write(6,*) nwv
	rmid=6.0d0
	open(unit=18,file='engy18.dat',status='unknown')
	rewind(18)
	write(18,990) zz,r0,nl,ntot
	write(6,990) zz,r0,nl,ntot,lwv,nwv
 990    format(f4.0,f6.1,2i4,12(i2,i2))
	ia=1
	do 200 iob=1,nl
CG        call pot(iob)
CG        l=lwv(iob)
	  l=iob/2
	  ld=l*2
	  int=1 - (-1)**iob
	  jd=ld + int - 1
	  xjj=jd/2.d0
	  ijj=1
	  ll=l
	  corr=1.d0
	  if(ll.eq.1) corr=1.d0/1.2210d0
	  if(ll.eq.2) corr=1.d0/1.2855d0
 
	hh=0.04d0-0.01*l
	if(l.gt.3) hh=0.004
	etry=eng(iob)
	do 10 ie=1,nwv(1,iob)
	if (ie.gt.2) hh=0.1d0*(en(ie-1)-en(ie-2))
	if (ie.eq.2) etry=en(ie-1)+0.10*dabs(en(ie-1))
	if (ie.gt.2) etry=en(ie-1)+0.1d0*(en(ie-1)-en(ie-2))
	nodeol=node
 3000    continue
c	write(6,*) 'etry,hh,',etry,hh
	if(hh.lt.1.0d-6) stop 'hh too small, infinite loop'
CG        write(6,*) 'll,xjj,ijj,etry,hh,irt'
CG        read(5,*) ll,xjj,ijj,etry,hh,irt
CG        l = ll
	call root(etry,hh,irt,ep,e0,err,rwave)
 
	do 20 i=1,m
  20    qq(i)=(ps(1,i)**2)*2.d0*x(i)
	xnorm=rint(qq,1,m,10,xh)
	xnorm=1.d0/dsqrt(xnorm)
	do 30 i=1,m
	y(ia,i)=ps(1,i)*xnorm
  30    yd(ia,i)=ps(2,i)*xnorm
chg
chg
	derv=ps(2,m)*xnorm/(2.d0*xf)
	yr0=0.d0
c nodes
	node=0
	do 32 i=3,m-1
	if(dsign(1.d0,y(ia,i))*dsign(1.d0,y(ia,i-1)).lt.0.)node=node+1
  32    continue
	write(18,35)e0,derv,yr0,node
	write(6,35)e0,derv,yr0,node
chg
chg -- remove previous statement
	if(ie.eq.1) goto 36
	if(node.ne.nodeol+1) then
	hh=hh/2.d0
	etry=(etry+en(ie-1))/2.d0
	goto 3000
	endif
  36    continue
  35    format(1p3d22.15,i5)
	ia=ia+1
	en(ie)=e0
  10    continue
 
 
	do 90 ie=1,nwv(2,iob)
	e1=(en(ie+ncb(iob))+en(ie+ncb(iob)+1))/2.d0
	wron = yfcnn(e1)
	do 100 i=1,m
 100    qq(i)=(ps(1,i)**2)*2.d0*x(i)
	xnorm=rint(qq,1,m,10,xh)
	xnorm=1.d0/dsqrt(xnorm)
	do 110  i=1,m
	y(ia,i)=ps(1,i)*xnorm
 110    yd(ia,i)=ps(2,i)*xnorm
	derv=ps(2,m)*xnorm/(2.d0*xf)
	yr0=ps(1,m)*xnorm
c nodes
	node=0
	do 132 i=3,m-1
	if(dsign(1.d0,y(ia,i))*dsign(1.d0,y(ia,i-1)).lt.0.)node=node+1
 132    continue
	write(18,35)e1,derv,yr0,node
	write(6,35)e1,derv,yr0,node
	ia=ia+1
  90    continue
 200    continue
 
 
**** <i|j> matrix *****
	iend=0
	do 151 ii=1,nl
	ibegin=iend+1
	iend=iend+nwv(1,ii)+nwv(2,ii)
	do 142 i=ibegin,iend
	do 142 j=ibegin,i
	do 143 k=1,m
 143    qq(k)=y(i,k)*y(j,k)*2.d0*x(k)
 142    wk(i,j)=rint(qq,1,m,10,xh)
	write(18,21) ((wk(i,j),j=ibegin,i),i=ibegin,iend)
 151    continue
**** <i|r|j> matrix(radial only) *******
	do 300 i=1,ntot
	do 300 j=1,i
	do 290 k=1,m
 290    qq(k)=y(i,k)*y(j,k)*rr(k)*2.d0*x(k)
 300    wk(i,j)=rint(qq,1,m,10,xh)
	write(18,21) ((wk(i,j),j=1,i),i=1,ntot)
****<i|1/r|j> matrix (radial only)*******
	do 500 i=1,ntot
	do 500 j=1,i
	do 490 k=1,m
 490    qq(k)=y(i,k)*y(j,k)*2.d0/x(k)
 500    wk(i,j)=rint(qq,1,m,10,xh)
	write(18,21) ((wk(i,j),j=1,i),i=1,ntot)
****<i|d/dr|j> matrix (radial only)*******
	do 400 i=1,ntot
	do 400 j=1,ntot
	do 390 k=1,m
 390    qq(k)=y(i,k)*(yd(j,k)-y(j,k)*2.d0/x(k))
 400    wk(i,j)=rint(qq,1,m,10,xh)
	write(18,21) ((wk(i,j),j=1,ntot),i=1,ntot)
****Write no. of points , r(1), step size and wavefunctions to unit 19
****  This version of wave writes this out in binary format
	h=h0*x(m)/(m-1)
	open(unit=19,file='wfn19.dat',status='unknown',
     1		access='sequential',form='unformatted')
	rewind(19)
	write(19)m,rr(1),xh,r0
	write(19)((y(i,j),j=1,m),i=1,ntot)
c       write(19)((yd(i,j),j=1,m),i=1,ntot)
        write(6,*) 'test:',y(1,1),y(1,m),y(ntot,1),y(ntot,m)
        write(6,*) 'test2:',y(1,10),y(1,m-10),y(ntot,10),y(ntot,m-10)
  21    format(1p3d22.15)
  37    format(i5,1p3d22.15)
	end
******************************************************************
	subroutine readin
	parameter(kw=390,kms=4321,nlj=13)
	implicit double precision(a-h,o-z)
	dimension ndevic(nlj)
	common/blok/y(kw,kms),ps(2,kms),
     &           rr(kms),qq(kms),x(kms),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid
	common/ncbcom/ncb(nlj)
 
 
**** readin r0, estimated energies for s, p, and d states.
  2     format(3i2)
	open(unit=25,file='OneElectron.dat',status='old')
	rewind(25)
	read(25,*) r0,m,nl,ntot
	write(6,*) r0,m,nl,ntot
	read(25,*) (lwv(i),i=1,nl)
	write(6,*) (lwv(i),i=1,nl)
	read(25,*) ((nwv(i,j),i=1,2),j=1,nl)
	write(6,*) ((nwv(i,j),i=1,2),j=1,nl)
	n=m
	read(25,*) (eng(i),i=1,nl)
	write(6,*) (eng(i),i=1,nl)
	read(25,*) (ncb(i),i=1,nl)
	write(6,*) (ncb(i),i=1,nl)
	read(25,*) (ndevic(i),i=1,4)
	write(6,*) (ndevic(i),i=1,4)
	close(25)
	h0 = 0.05d0
	rv0 = 0.263d-03
 
 
	ni=300
	zz=2.d0
c
c ** the following statement should be used only if an
c       analytical potential is used in fpot instead of
c       a corrected Hartree-Slater potential
c
		IF(zz.GT.0.d0) RETURN
c
	do 10 j=1,4
	read(ndevic(j),1)rv0,h0
	if(j.gt.1) goto 9
	do 8 i=1,ni
   8    ri(i)=rv0*dexp((i-1)*h0)
	write(6,*)rv0,h0
   1    format(1p5d14.7)
   9    read(ndevic(j),1)(zi(i),i=1,ni)
	do 100 i=1,ni
  100   vs(j,i)=-zi(i)*ri(i)
  10    continue
	write(6,*) 'charge at ri(ni)=', ri(ni), ' is ', vs(1,ni)
	return
	end
************************************************************************
c --The following subroutine is used only when numerical potentials are
c ------- used.
	subroutine pot(k)
	parameter(kw=390,kms=4321,nlj=13)
	implicit double precision(a-h,o-z)
	common/blok/y(kw,kms),ps(2,kms),
     &           rr(kms),qq(kms),x(kms),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid
	k1=k
	if(k.gt.4) k1=4
	do 10 j=1,ni
10      zi(j)=vs(k1,j)
	return
	end
