******************************************************************************
c	gksub.f
c
c	This is a library of subroutines to do two-electron radial integrals,
c		 and to solve differential equations numerically
c *****************************************************************************
c
c -- function sl evaluates numerically the radial portion of a
c ---- two-electron integral 1/r12 using an indep. eln. basis
c ---- for the l'th multipole.  Be sure to initialize common blocks.
c ---- This version is for a sqrt mesh of radial points.
c
        function sl(ra,rb,rc,rd,l,l2,l4,a2,a4,m,max)
        implicit real*8(a-h,o-z)
        common/radial/r(4321),rp(4321),rpor(4321),h
        dimension ra(4321),rb(4321),rc(4321),rd(4321)
        dimension x(4321),q(4321),f(4321),y(4321)


c   ra,rb,rc,rd are the radial wavefns. l=multipolarity , l2 is ang.mom.
c  of rb orbital , l4 ang. mom. of rd. a2,a4 can be set to 1 .
c     m=max =no.of radial points.

        a2=1.d0
        a4=1.d0
        ssl=0.d0
        do 10 i=1,m
        xr=dsqrt(r(i))
        x(i)=rb(i)*rd(i)  *2.d0*xr
  10    q(i)=ra(i)*rc(i)  *2.d0*xr
c       x(i)=rb(i)*rd(i)
c 10    q(i)=ra(i)*rc(i)
***yfun gives the hartree y-function
 
        call yfun(x,y,l,l2,l4,a2,a4,m,max)
        do 20 i=1,m
  20    f(i)=q(i)*y(i)
        ssl=rint(f,1,m,10,h) + q(1)*y(1)*r(1)/2
	sl=ssl
        return
        end
*******************************************************************************
        subroutine yfun(x,y,l,l1,l2,a1,a2,m,max)
 
***this program calculates the hartree y-function :y(l,r)/r
c  x : input function
c  y :  output hartree y-function
c  l : order of y-function
c  m : number of tabulation points for x
c
        implicit real*8(a-h,o-z)
        common/radial/r(4321),rp(4321),rpor(4321),h
     &       /temsp/u(4321),v(4321),w(4321),s(4321)
        dimension x(m),y(m)
        if(l.gt.0) go to 35
  20    do 22 i=1,m
        v(i)=x(i)
  22    w(i)=x(i)/r(i)
        call yint(v,w,y,u,r,l,l1,l2,a1,a2,m,h)
        ym=y(m)
        do 24 i=1,m
  24    y(i)=y(i)/r(i) +u(i)
        m1=m+1
        if(m1.gt.max) go to 900
        do 26 i=m1,max
  26    y(i)=ym/r(i)
        go to 900
  35    continue
        do 40 i=1,max
  40    s(i)=r(i)**l
        do 50 i=1,m
        v(i)=x(i)*s(i)
  50    w(i)=x(i)/(r(i)*s(i))
        call yint(v,w,y,u,r,l,l1,l2,a1,a2,m,h)
        ym=y(m)
        do 60 i=1,m
  60    y(i)=y(i)/(r(i)*s(i)) + u(i)*s(i)
        m1=m+1
        if(m1.gt.max) go to 900
        do 90 i=m1,max
  90    y(i)=ym/(r(i)*s(i))
  900   return
        end
*****************************************************************************
        subroutine yint(v,w,y,z,r,l,l1,l2,a1,a2,m,h)
c
c  this program calculates the indefinite integrals y and z
c        using the lagragian integration formula
c
c   y(r) = integral of v from 0 to r
c   z(r) = integral of w from r to infinity
c   m is the maximum tabulation point of the integrands
c   h is the step size
c
        implicit real*8(a-h,o-z)
        dimension v(m),w(m),y(m),z(m),r(m)
cg    data ho/0.d0/
c
c------------------------------------------------------------
c         lagrange 6 point integration formula
c               ------------------------------------
        dimension aa(3,6),a(6,3),b(3)
        equivalence (b(1),a(4,3))
	data ho/0.d0/
        data ia/3/, ja/6/, da/1440.d0/, a/18*0.d0/,
     1     aa/  475.d0,   -27.d0,  11.d0,
     2         1427.d0,   637.d0, -93.d0,
     3         -798.d0,  1022.d0, 802.d0,
     4          482.d0,  -258.d0, 802.d0,
     5         -173.d0,    77.d0, -93.d0,
     6           27.d0,   -11.d0,  11.d0/
c------------------------------------------------------------
        if(h.eq.ho) go to 20
        hd=h/da
        do 10 i=1,ia
        do 10 j=1,ja
  10    a(j,i)=aa(i,j)*hd
        ho=h
  20    z(m)=0.d0
        lsum=l1+l2+l+3
        y(1)=a1*a2*(r(1)**lsum)/lsum
        do 30 i=2,ia
        k=m-i+1
        y(i)=y(i-1)
        z(k)=z(k+1)
        ii=i-1
        do 30 j=1,ja
        y(i)=y(i)+a(j,ii)*v(j)
  30    z(k)=z(k)+a(j,ii)*w(m-j+1)
        im=ia + 1
        in=m -ia +1
        do 40 i=im,in
        k=m-i+1
        y(i)=y(i-1)+b(1)*(v(i)+v(i-1))+b(2)*(v(i+1)+v(i-2))
     &    +b(3)*(v(i+2)+v(i-3))
  40    z(k)=z(k+1)+b(1)*(w(k)+w(k+1))+b(2)*(w(k-1)+w(k+2))
     &    +b(3)*(w(k-2)+w(k+3))
        in=in+1
        do 50 i=in,m
        k=m-i+1
        y(i)=y(i-1)
        z(k)=z(k+1)
        do 50 j=1,ja
        y(i)=y(i)+a(j,k)*v(m-j+1)
  50    z(k)=z(k)+a(j,k)*w(j)
        return
        end
****************************************************************************
****************************************************************************
	function rwave(energy)
c
c -- function rwave is the main driver program used to solve a radial
c ---- differential equation using the subroutine prco, which uses
c ---- a 4-point predictor-corrector method.  rwave uses a sqrt mesh.
c
      ! maab: This must be changed for any change in the input files.  
       parameter(kw=390,nlj=13)
	implicit real*8(a-h,o-z)
	dimension yy(4,2), psa(2,4321),psb(2,4321), xa(4321),xb(4321)
        common/blok/y(kw,4321),ps(2,4321),
     &           rr(4321),qq(4321),x(4321),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid


	external fpot
c
c   calculate wavefunction usin prco pred-corr method
c
c      write(6,*) "Inside of rwave ",
        n=m
        r=x0*x0
	e=energy
	cc=1.d0
        z = dabs(zi(1))
        rho = 2.d0*z*r
        yy(1,1)=rho**(l+1)-(1.d0/(2*l+2.d0))*(rho**(l+2))
        yy(1,1)=yy(1,1)+(.5d0/((2*l+2.d0)*(2*l+3.d0)))*(rho**(l+3))
        yy(1,2)=(l+1)*(2*z)*rho**l
        yy(1,2)=yy(1,2)-(2*z)*((l+2.d0)/(2*l+2.d0))*rho**(l+1)
        yy(1,2)=yy(1,2)+z*((l+3.d0)/((2*l+2.d0)*(2*l+3.d0)))*rho**(l+2)
        yy(1,2)=yy(1,2)*(2.d0*x0)
        na=dsqrt(rmid)/xh
        nb=n-na+1
        xm=xh*na
c        write(6,120) rv0*dexp(h0*xm),na,nb
 120    format(2x,'xm,na,nb',f10.5,2i5)
        do 50 i=1,na
  50    xa(i)=x(i)
         call prco(xa,1,na,yy,psa,2,fpot)
         yy(1,1)=0.d0
         yy(1,2)=1.d-10
        do 55 i=1,nb
  55    xb(i)=x(m-i+1)
        call prco(xb,1,nb,yy,psb,2,fpot)
        wron=psa(1,na)*psb(2,nb)-psa(2,na)*psb(1,nb)
        rwave=wron
        do 30 i=1,na
        ps(1,i)=psa(1,i)
  30    ps(2,i)=psa(2,i)
        fit=psa(1,na)/psb(1,nb)
        do 31  i=1,nb-1
        ps(1,i+na)=psb(1,nb-i)*fit 
  31    ps(2,i+na)=psb(2,nb-i)*fit
c       write(6,350) energy,wron
 350    format(1p2e18.10)
	return
	end
****************************************************************************
	function yfcnn(energy)
        parameter(kw=390,nlj=13)
	implicit real*8(a-h,o-z)
	dimension yy(4,2)
        common/blok/y(kw,4321),ps(2,4321),
     &           rr(4321),qq(4321),x(4321),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid


	external fpot
c
c   calculate wavefunction using prco pred-corr method
c
        n=m
        r=x0*x0
	e=energy
	cc=1.d0
	       z = dabs(zi(1))
        rho = 2.d0*z*r
        yy(1,1)=rho**(l+1)-(1.d0/(2*l+2.d0))*(rho**(l+2))
        yy(1,1)=yy(1,1)+(.5d0/((2*l+2.d0)*(2*l+3.d0)))*(rho**(l+3))
        yy(1,2)=(l+1)*(2*z)*rho**l
        yy(1,2)=yy(1,2)-(2*z)*((l+2.d0)/(2*l+2.d0))*rho**(l+1)
        yy(1,2)=yy(1,2)+z*((l+3.d0)/((2*l+2.d0)*(2*l+3.d0)))*rho**(l+2)
        yy(1,2)=yy(1,2)*(2.d0*x0)


        call prco(x,1,n,yy,ps,2,fpot)
        qy=ps(1,n)
        qyp=ps(2,n)/(2.d0*xf)
        wron=qy
        yfcnn= wron
c       write(6,350) energy,wron
 350    format(2e18.10)
	return
        end
**********************************************************************
***********************************************************************
	subroutine cpot(xx,f)
	implicit double precision(a-h,o-z)
	dimension f(2,2)
        common /blk1/zi,ekk,le


c  evaluate potl. for pred-corr
c  ek=kinetic energy
        r=xx*xx
	z=zi
        l=le
        ek=ekk
	v=-z/r +0.5d0*l*(l+1)/(r*r)
	f(1,1)=0.d0
	f(1,2)=1.d0
        f(2,1)= -8.d0*r*(ek-v)
        f(2,2)=1.d0/xx
	return
	end
**********************************************************************
********************************************************************************
      subroutine prco(x,n0,n,yy,ps,m,fcn)
c
c -- general differential equation solver using a predictor-corrector
c ----- algorithm after starting the integration using runge-kutta
c
      implicit real*8(a-h,o-z)
      dimension x(m),yy(4,2),ps(2,4321),ym3(2),dyy(4,2)
      dimension ym2(2),ym1(2),ym0(2),dym3(2),dym2(2)
      dimension dym1(2),dym0(2),f(2,2),yp1(2),dyp1(2)
      dimension qyp1(2),qdyp1(2),wyp1(2),wdyp1(2)
      dimension qy(02),y(10),yp(10),dy(2)
 	external fcn,fpot
      ep = 1.d-10
      h = (x(n) - x(n0))/(n-n0)
	i0=n0
	i1=n0+1
	i2=n0+2
	i3=n0+3
	i4=n0+4
 
      y(1) = yy(1,1)
      yp(1) = yy(1,2)
      qy(1) = y(1)
      qy(2) = yp(1)
      call fcn(x(i0),f)
      call mmlt(dy,f,qy,2)
      dyy(1,1) = dy(1)
      dyy(1,2) = dy(2)
      x(i1) = x(i0) + h
      call runge2(x(i0),x(i1),qy,5,2,fpot)
c     qy(1) = y(1)
c     qy(2) = yp(1)
      call fcn(x(i1),f)
      call mmlt(dy,f,qy,2)
      do 7330 i = 1,2
      yy(2,i) = qy(i)
      dyy(2,i) = dy(i)
 7330 continue
      x(i2) = x(i1) + h
      call runge2(x(i1),x(i2),qy,5,2,fpot)
c     qy(1) = y(1)
c     qy(2) = yp(1)
      call fcn(x(i2),f)
      call mmlt(dy,f,qy,2)
      do 7332 i = 1,2
      yy(3,i) = qy(i)
      dyy(3,i) = dy(i)
 7332 continue
      x(i3) = x(i2) + h
      call runge2(x(i2),x(i3),qy,5,2,fpot)
c     qy(1) = y(1)
c     qy(2) = yp(1)
      call fcn(x(i3),f)
      call mmlt(dy,f,qy,2)
      do 7334 i = 1,2
      yy(4,i) = qy(i)
      dyy(4,i) = dy(i)
 7334 continue


      n04 = n0 + 4
      do 5 j = 1,m
      ym3(j) = yy(1,j)
      ym2(j) = yy(2,j)
      ym1(j) = yy(3,j)
      ym0(j) = yy(4,j)
      ps(j,n0) = ym3(j)
      ps(j,n0+1) = ym2(j)
      ps(j,n0+2) = ym1(j)
      ps(j,n0+3) = ym0(j)
    5 continue
      do 6 i = 1,3
      x(n0+i) = x(n0) + i*h
    6 continue
 
      call fcn(x(n0),f)
      call mmlt(dym3,f,ym3,m)
      xx = x(n0+1)
      call fcn(xx,f)
      call mmlt(dym2,f,ym2,m)
      xx = x(n0+2)
      call fcn(xx,f)
      call mmlt(dym1,f,ym1,m)
      xx = x(n0+3)
      call fcn(xx,f)
      call mmlt(dym0,f,ym0,m)
 
      do 500 i = n04,n
      x(i) = x(n0+3) + (i-n0-3)*h
      xx = x(i)
      ist = i
cccc  call fc2(xx,f)
      call fcn(xx,f)
      call ppred(qyp1,ym0,dym0,dym1,dym2,dym3,h,m)
      call mmlt(qdyp1,f,qyp1,m)
c	write(6,*) qyp1(1),qyp1(2),"qyp1"
c	write(6,*) qdyp1(1),qdyp1(2),"qdyp1"
      do 30 j = 1,m
      wyp1(j) = qyp1(j)
      wdyp1(j) = qdyp1(j)
   30 continue
      do 50 ii = 1,10
      call ccorr(yp1,ym0,qdyp1,dym0,dym1,dym2,h,m)
      err = 0.d0
      qpsi = 0.d0
      do 45 jj = 1,m
      err = err + (yp1(jj) - qyp1(jj)) ** 2
      qpsi = qpsi + yp1(jj) ** 2
   45 continue
      err = dsqrt(err/qpsi)
      call mmlt(dyp1,f,yp1,m)
c	write(6,*) yp1(1),yp1(2),"yp1"
c	write(6,*) dyp1(1),dyp1(2),"dyp1"
      if(err.lt.ep) go to 60
c     write(6,46) yp1(1),wyp1(1),yp1(2),wyp1(2)
   46 format(4(1x,1pd11.4))
      do 48 jj = 1,m
      qyp1(jj) = yp1(jj)
      qdyp1(jj) = dyp1(jj)
   48 continue
   50 continue
c     write(6,56) xx,err,ep
   60 do 53 jj = 1,m
      yp1(jj) = yp1(jj) + (wyp1(jj)-yp1(jj))*19.d0/270.d0
   53 continue
      call mmlt(dyp1,f,yp1,m)
 
c     write(6,56)xx,err,ep
   56 format(1x,'no convergence in prco, x=',3(1x,1pd12.5))
  960 continue
 
      do 75 jj = 1,m
      ps(jj,i) = yp1(jj)
      ym3(jj) = ym2(jj)
      ym2(jj) = ym1(jj)
      ym1(jj) = ym0(jj)
      ym0(jj) = yp1(jj)
      dym3(jj) = dym2(jj)
      dym2(jj) = dym1(jj)
      dym1(jj) = dym0(jj)
      dym0(jj) = dyp1(jj)
   75 continue
c     if( 20*(i/20) .ne. i) go to 500
c     write(6,77) xx,ps(1,i),ps(2,i)
c  77 format(3(1x,1pd12.5))
  500 continue
      return
      end
*******************************************************************************
      subroutine ppred(yp1,y0,dy0,dym1,dym2,dym3,h,m)
      implicit real*8(a-h,o-z)
      dimension yp1(2),y0(2),dy0(2),dym1(2),dym2(2),dym3(2)
      do 100 i = 1,m
      q = 55.d0 * dy0(i) - 59.d0*dym1(i)
      q = q + 37.d0*dym2(i) - 9.d0*dym3(i)
      yp1(i) = y0(i) + q * h/24.d0
  100 continue
      return
      end
*******************************************************************************
      subroutine mmlt(a,b,c,m)
      implicit real*8(a-h,o-z)
      dimension a(2),b(2,2),c(2)
      do 100 i = 1,m
      a(i) = 0.d0
      do 50 j = 1,m
      a(i) = a(i) + b(i,j) * c(j)
   50 continue
  100 continue
      return
      end
********************************************************************************
      subroutine ccorr(yp1,y0,dyp1,dy0,dym1,dym2,h,m)
      implicit real*8(a-h,o-z)
      dimension yp1(2),y0(2),dyp1(2),dy0(2),dym1(2),dym2(2)
      do 100 i = 1,m
      q = 9.d0 * dyp1(i) + 19.d0*dy0(i) - 5.d0 * dym1(i)
      q = q + dym2(i)
      yp1(i) = y0(i) + q * h / 24.d0
  100 continue
      return
      end
*******************************************************************************
      subroutine runge2(x0,xf,y,npts,m,fcn)                             
      implicit real*8(a-h,o-z)                                          
      dimension a(2),b(2),c(2),d(2),q(2),f(2,2)                         
      dimension y(2),z1(2),z2(2)                                        
	external fcn
                                                                        
      h = (xf-x0)/(npts-1)                                              
      nd1 = npts - 1                                                    
      x = x0                                                            
      call fcn(x,f)                                                     
      do 100 i = 1,nd1                                                  
      j = i + 1                                                         
      do 5 l = 1,m                                                      
      a(l) = 0.d0                                                       
      do 4 k = 1,m                                                      
      a(l) = a(l) + f(l,k) * y(k)                                       
    4 continue                                                          
      a(l) = a(l) * h                                                   
      q(l) = y(l) + 0.5d0 * a(l)                                        
    5 continue                                                          
      xint = x + h/2.d0                                                 
      call fcn(xint,f)                                                  
      do 10 l = 1,m                                                     
      b(l) = 0.d0                                                       
      do 8 k = 1,m                                                      
      b(l) = b(l) + f(l,k)*q(k)                                         
    8 continue                                                          
      b(l) = b(l) * h                                                   
      z1(l) = y(l) + 0.5d0 * b(l)                                       
   10 continue                                                          
      do 15 l = 1,m                                                     
      c(l) = 0.d0                                                       
      do 13 k = 1,m                                                     
      c(l) = c(l) + f(l,k)*z1(k)                                        
   13 continue                                                          
      c(l) = c(l) * h                                                   
      z2(l) = y(l) + c(l)                                               
   15 continue                                                          
      x = x + h                                                         
      call fcn(x,f)                                                     
      do 20 l = 1,m                                                     
      d(l) = 0.d0                                                       
      do 18 k = 1,m                                                     
      d(l) = d(l) + f(l,k) * z2(k)                                      
   18 continue                                                          
      d(l) = d(l) * h                                                   
      y(l) = y(l) + (a(l)+2.d0*b(l)+2.d0*c(l)+d(l))/6.d0                
   20 continue                                                          
  100 continue                                                          
      return                                                            
      end                                                            
	subroutine fpot(xx,f)
	parameter(kw=390,nlj=13)
	implicit real*8(a-h,o-z)
	dimension f(2,2),aal1(7),aal2(7),aal3(7),rrc(7)
	common/blok/y(kw,4321),ps(2,4321),
     &           rr(4321),qq(4321),x(4321),vs(nlj,351),
     &           ri(351),zi(351),eng(nlj),e,r0,etry,rv0,xh,x0,xf,
     &           h0,zz,l,m,ni,ntot,nl,nwv(2,nlj),lwv(nlj),rmid
	common/sporb/xjj,corr,ll,ijj
	common/rmc/vrmc
c -- the following set of parameters is from Aymar's paper I
cI        data al1,al2,al3,zion/3.4675d0,2.6357d0,1.0078d0,38.d0/
c -- the following set is from Aymar's paper II:
ccc	data al1,al2,al3,zion/3.55156d0,6.03726d0,1.43891d0,38.d0/
c -- New refitted potential parameters for Sr in jj-coupling:
c 	data aal1/3.4187d0,3.46432d0,3.25329d0,5.35402d0,5.35402/
c 	data aal2/4.7332d0,4.3960d0,3.2330d0,6.7517d0,6.7517/
c 	data aal3/1.59149d0,1.64505d0,1.59956d0,5.66239d0,5.66239/
c 	data rrc/1.796570d0,1.522933d0,1.681960d0,1.005666d0,1.0057/
c 	data zion,alpha/38.d0,7.5d0/
      data corrP,corrD/0.948042d0,0.924857d0/
c -- Aymar's new fitted potentials for barium:
c	data aal1/3.07514,3.23037,3.29609,3.62371,3.62371/
c	data aal2/2.6107,2.9561,3.0248,6.7416,6.7416/
c	data aal3/1.20259,1.19230,1.29432,2.03792,2.03792/
c	data rrc/2.600400,2.049703,1.894587,1.047344,1.047344/
c	data zion,alpha/56.d0,11.4d0/
c
c -- Miguel's fitted parameters for Yb in jj coupling
c      data aal1/4.08411932d0,  9.58745009d0,  3.67940204d0,  
c     1          4.92419126d0, 4.92419126d0, 4.92419126d0, 4.92419126d0/
c      data aal2/15.72608295d0, 60.15460544d0, 4.78566475d0,
c     1          27.54429042d0, 27.54429042d0, 27.54429042d0
c     1          , 27.54429042d0/
c      data aal3/2.33138144d0, 2.01238056d0, 1.80601090d0, 
c     1          1.84848583d0, 1.84848583d0, 1.84848583d0, 1.84848583d0/
c      data rrc/1.22728665d0, 4.27297600d0,  1.73219589d0, 
c     1         2.52242492d0,  2.52242492d0, 2.52242492d0, 2.52242492d0/
c      data zion,alpha/70.d0,7.35d0/
c -- Radium parameters
      data aal1/3.7702,3.9430,3.7008,3.8125,3.8125,3.8125,3.8125/
      data aal2/4.9928,5.0552,4.7748,5.0332,5.0332,5.0332,5.0332/
      data aal3/1.5179,3.6770,1.4956,2.1016,2.1016,2.1016,2.1016/
      data rrc/1.3691,1.0924,2.2784,1.2707,1.2707,1.2707,1.2707/
      data zion,alpha/88.0d0,18.0d0/
	data ccc/137.036d0/
c  evaluate potl. for pred-corr
c
	zi(1)=zion
	al1=aal1(ll+1)
	al2=aal2(ll+1)
	al3=aal3(ll+1)
	rc=rrc(ll+1)

      corr=1.d0
      IF(dabs(zion-38.d0).LT.0.001) THEN
          IF(ll.EQ.1) corr= corrP
          IF(ll.EQ.2) corr= corrD
      ENDIF

c
	n=m
	r=xx*xx
	z=zz
CG        if(r.gt.ri(ni)) goto 10
CG      z = yfix(r,ri,zi,ni,0)
  10    continue
c ** the following line contains Aymar's analytical potential
	f1=dexp(-al1*r)*(zion-2.d0)
	f2=dexp(-al3*r)*al2
CG      z = 2.d0 + (zion-2.d0)*dexp(-al1*r)+al2*r*dexp(-al3*r)
	z = 2.d0 + f1 + r*f2
	zp= -al1*f1 + f2*(1.d0-al3*r)
	expf=dexp(-(r/rc)**6)
	upol=(-0.5d0*alpha/r**4)*(1.d0-expf)
	upolp=(-5.d0*upol/r)- (3.d0*alpha*r/rc**6)*expf
	vp= -zp/r + z/r/r +upolp
	vp= 0.25d0*(xjj*(xjj+1)-ll*(ll+1)-0.75d0)*vp/r/ccc/ccc
	vp=vp/ ( 1.d0 - (upol-z/r)/(ccc*ccc*2.d0) )**2
	v= upol-z/r +0.5d0*ll*(ll+1)/(r*r) + ijj*vp  *corr
ccc	vrmc = - 0.5d0*(upol-z/r-e)**2 /ccc/ccc
ccc	v = v + vrmc

	f(1,1)=0.d0
	f(1,2)=1.d0
	f(2,1)= -8.d0*r*(e-v)
	f(2,2)=1.d0/xx
	return
	end
