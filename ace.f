C     real -> double precision conversion for R use
C     <TSL>
c     mortran 2.0     (version of 6/24/75)
      subroutine mace (p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
      implicit none
c
c   subroutine mace(p,n,x,y,w,l,delrsq,ns,tx,ty,rsq,ierr,m,z)
c------------------------------------------------------------------
c
c estimate multiple optimal transformations for regression and
c correlation by alternating conditional expectation estimates.
c
c version 3/28/85.
c
c breiman and friedman, journal of the american statistical
c association (september, 1985)
c
c coded  and copywrite (c) 1985 by:
c
c                        jerome h. friedman
c                     department of statistics
c                               and
c                stanford linear accelerator center
c                        stanford university
c
c all rights reserved.
c
c
c input:
c
c    n : number of observations.
c    p : number of predictor variables for each observation.
c    x(p,n) : predictor data matrix.
c    y(n) : response values for the observations.
c       missing values are signified by a value (response or
c       predictor) greater than or equal to big.
c       (see below - default, big = 1.0e20)
c    w(n) : weights for the observations.
c    l(p+1) : flag for each variable.
c       l(1) through l(p) : predictor variables.
c       l(p+1) : response variable.
c       l(i)=0 => ith variable not to be used.
c       l(i)=1 => ith variable assumes orderable values.
c       l(i)=2 => ith variable assumes circular (periodic) values
c                 in the range (0.0,1.0) with period 1.0.
c       l(i)=3 => ith variable transformation is to be monotone.
c       l(i)=4 => ith variable transformation is to be linear.
c       l(i)=5 => ith variable assumes categorical (unorderable) values.
c   delrsq : termination threshold. iteration stops when
c       rsq changes less than delrsq in nterm
c       consecutive iterations (see below - default, nterm=3).
c   ns : number of eigensolutions (sets of transformations).
c
c output:
c
c   tx(n,p,ns) : predictor transformations.
c      tx(j,i,k) = transformed value of ith predictor for jth obs
c                  for kth eigensolution.
c   ty(n,ns) = response transformations.
c      ty(j,k) = transformed response value for jth observation
c                for kth eigensolution.
c   rsq(ns) = fraction of variance(ty<y>)
c                       p
c         explained by sum tx(i)<x(i)>  for each eigensolution.
c                      i=1
c   ierr : error flag.
c      ierr = 0 : no errors detected.
c      ierr > 0 : error detected - see format statements below.
c
c scratch:
c
c    m(n,p+1), z(n,12) : internal working storage.
c
c note: mace uses an iterative procedure for solving the optimization
c    problem. default starting transformations are ty(j,k)=y(j),
c    tx(j,i,k)=x(i,j) : j=1,n, i=1,p, k=1,ns. other starting transformat
c    can be specified (if desired) for either the response and/or any of
c    the predictor variables. this is signaled by negating the
c    corresponding l(i) value and storing the starting transformed
c    values in the corresponding array (ty(j,k), tx(j,i,k)) before
c    calling mace.
c
c------------------------------------------------------------------
c
      integer n,p,pp1,m(n,p+1),l(p+1)
      integer ns,ierr,i,is,ism1,iter,j,js,k,maxit,nit,np,nt,nterm
      double precision rsqi,span
      double precision alpha, big, cmn, cmx
      double precision y(n),x(p,n),w(n),ty(n,ns),tx(n,p,ns)
      double precision z(n,12),ct(10),rsq(ns)
      double precision delrsq
      common /prams/ alpha,big,span,maxit,nterm
      double precision sm,sv,sw,sw1
      ierr=0
      pp1=p+1
      sm=0.0
      sv=sm
      sw=sv
      sw1=sw
      do 10 i=1,pp1
      if (l(i).ge.-5.and.l(i).le.5) go to 10
      ierr=6
 10   continue
      if (ierr.ne.0) return
      if (l(pp1).ne.0) go to 20
      ierr=4
      return
 20   np=0
      do 30 i=1,p
      if (l(i).ne.0) np=np+1
 30   continue
      if (np.gt.0) go to 40
      ierr=5
      return
 40   do 50 j=1,n
      sw=sw+w(j)
 50   continue
      if (sw.gt.0.0) go to 60
      ierr=1
      return
 60   do 580 is=1,ns
      do 70 j=1,n
      if (l(pp1).gt.0) ty(j,is)=y(j)
 70   continue
      do 170 i=1,p
      if (l(i).ne.0) go to 90
      do 80 j=1,n
      tx(j,i,is)=0.0
 80   continue
      go to 170
 90   if (l(i).le.0) go to 110
      do 100 j=1,n
      tx(j,i,is)=x(i,j)
 100  continue
 110  do 120 j=1,n
      if (tx(j,i,is).ge.big) go to 120
      sm=sm+w(j)*tx(j,i,is)
      sw1=sw1+w(j)
 120  continue
      if (sw1.gt.0.0) go to 140
      do 130 j=1,n
      tx(j,i,is)=0.0
 130  continue
      sm=0.0
      sw1=sm
      go to 170
 140  sm=sm/sw1
      do 160 j=1,n
      if (tx(j,i,is).ge.big) go to 150
      tx(j,i,is)=tx(j,i,is)-sm
      go to 160
 150  tx(j,i,is)=0.0
 160  continue
      sm=0.0
      sw1=sm
 170  continue
      do 180 j=1,n
      if (ty(j,is).ge.big) go to 180
      sm=sm+w(j)*ty(j,is)
      sw1=sw1+w(j)
 180  continue
      if (sw1.gt.0.0) go to 190
      ierr=1
      return
 190  sm=sm/sw1
      do 210 j=1,n
      if (ty(j,is).ge.big) go to 200
      ty(j,is)=ty(j,is)-sm
      go to 210
 200  ty(j,is)=0.0
 210  continue
      do 220 j=1,n
      sv=sv+w(j)*ty(j,is)**2
 220  continue
      sv=sv/sw
      if (sv.le.0.0) go to 230
      sv=1.0/dsqrt(sv)
      go to 260
 230  if (l(pp1).le.0) go to 240
      ierr=2
      go to 250
 240  ierr=3
 250  return
 260  do 270 j=1,n
      ty(j,is)=ty(j,is)*sv
 270  continue
      if (is.ne.1) go to 310
      do 280 j=1,n
      m(j,pp1)=j
      z(j,2)=y(j)
 280  continue
      call sort (z(1,2),m(1,pp1),1,n)
      do 300 i=1,p
      if (l(i).eq.0) go to 300
      do 290 j=1,n
      m(j,i)=j
      z(j,2)=x(i,j)
 290  continue
      call sort (z(1,2),m(1,i),1,n)
 300  continue
 310  call scail (p,n,w,sw,ty(1,is),tx(1,1,is),delrsq,p,z(1,5),z(1,6))
      rsq(is)=0.0
      iter=0
      nterm=min0(nterm,10)
      nt=0
      do 320 i=1,nterm
      ct(i)=100.0
 320  continue
 330  iter=iter+1
      nit=0
 340  rsqi=rsq(is)
      nit=nit+1
      do 360 j=1,n
      z(j,5)=ty(j,is)
      do 350 i=1,p
      if (l(i).ne.0) z(j,5)=z(j,5)-tx(j,i,is)
 350  continue
 360  continue
      do 420 i=1,p
      if (l(i).eq.0) go to 420
      do 370 j=1,n
      k=m(j,i)
      z(j,1)=z(k,5)+tx(k,i,is)
      z(j,2)=x(i,k)
      z(j,4)=w(k)
 370  continue
      call smothr (iabs(l(i)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
      sm=0.0
      do 380 j=1,n
      sm=sm+z(j,4)*z(j,3)
 380  continue
      sm=sm/sw
      do 390 j=1,n
      z(j,3)=z(j,3)-sm
 390  continue
      sv=0.0
      do 400 j=1,n
      sv=sv+z(j,4)*(z(j,1)-z(j,3))**2
 400  continue
      sv=1.0-sv/sw
      if (sv.le.rsq(is)) go to 420
      rsq(is)=sv
      do 410 j=1,n
      k=m(j,i)
      tx(k,i,is)=z(j,3)
      z(k,5)=z(j,1)-z(j,3)
 410  continue
 420  continue
      if (np.eq.1.or.rsq(is)-rsqi.le.delrsq.or.nit.ge.maxit) go to 430
      go to 340
 430  do 450 j=1,n
      k=m(j,pp1)
      z(j,2)=y(k)
      z(j,4)=w(k)
      z(j,1)=0.0
      do 440 i=1,p
      if (l(i).ne.0) z(j,1)=z(j,1)+tx(k,i,is)
 440  continue
 450  continue
      call smothr (iabs(l(pp1)),n,z(1,2),z,z(1,4),z(1,3),z(1,6))
      if (is.le.1) go to 490
      ism1=is-1
      do 480 js=1,ism1
      sm=0.0
      do 460 j=1,n
      k=m(j,pp1)
      sm=sm+w(k)*z(j,3)*ty(k,js)
 460  continue
      sm=sm/sw
      do 470 j=1,n
      k=m(j,pp1)
      z(j,3)=z(j,3)-sm*ty(k,js)
 470  continue
 480  continue
 490  sm=0.0
      sv=sm
      do 500 j=1,n
      k=m(j,pp1)
      sm=sm+w(k)*z(j,3)
      z(k,2)=z(j,1)
 500  continue
      sm=sm/sw
      do 510 j=1,n
      z(j,3)=z(j,3)-sm
      sv=sv+z(j,4)*z(j,3)**2
 510  continue
      sv=sv/sw
      if (sv.le.0.0) go to 520
      sv=1.0/dsqrt(sv)
      go to 530
 520  ierr=3
      return
 530  do 540 j=1,n
      k=m(j,pp1)
      ty(k,is)=z(j,3)*sv
 540  continue
      sv=0.0
      do 550 j=1,n
      sv=sv+w(j)*(ty(j,is)-z(j,2))**2
 550  continue
      rsq(is)=1.0-sv/sw
      nt=mod(nt,nterm)+1
      ct(nt)=rsq(is)
      cmn=100.0
      cmx=-100.0
      do 560 i=1,nterm
      cmn=min(cmn,ct(i))
      cmx=max(cmx,ct(i))
 560  continue
      if (cmx-cmn.le.delrsq.or.iter.ge.maxit) go to 570
      go to 330
 570  continue
 580  continue
      return
      end

      subroutine model (p,n,y,w,l,tx,ty,f,t,m,z)
      implicit none
c
c          subroutine model(p,n,y,w,l,tx,ty,f,t,m,z)
c--------------------------------------------------------------------
c
c computes response predictive  function f for the model yhat = f(t),
c where
c                                        p
c            f(t) = e(y : t),     t =   sum  tx<i> ( x<i> )
c                                       i=1
c using the x transformations tx constructed by subroutine ace.
c if y is a categorical variable (classification) then
c                                -1
c                       f(t) = ty  (t).
c input:
c
c    p,n,y,w,l : same input as for subroutine ace.
c    tx,ty,m,z : output from subroutine ace.
c
c output:
c
c    f(n),t(n) : input for subroutine acemod.
c
c note: this subroutine must be called before subroutine acemod.
c
c-------------------------------------------------------------------
c
      integer n,p,pp1,m(n,1),l(1)
      integer i,j,j1,j2,k,maxit,nterm
      double precision y(n),w(n),tx(n,p),ty(n),f(n),t(n),z(n,12)
      double precision alpha, big, s, span
      common /prams/ alpha,big,span,maxit,nterm
      pp1=p+1
      if (iabs(l(pp1)).ne.5) go to 20
      do 10 j=1,n
      t(j)=ty(j)
      m(j,pp1)=j
 10   continue
      go to 50
 20   do 40 j=1,n
      s=0.0
      do 30 i=1,p
      s=s+tx(j,i)
 30   continue
      t(j)=s
      m(j,pp1)=j
 40   continue
 50   call sort (t,m(1,pp1),1,n)
      do 140 j=1,n
      k=m(j,pp1)
      z(j,2)=w(k)
      if (y(k).ge.big) go to 60
      z(j,1)=y(k)
      go to 140
 60   j1=j
      j2=j1
 70   if (y(m(j1,pp1)).lt.big) go to 80
      j1=j1-1
      if (j1.lt.1) go to 80
      go to 70
 80   if (y(m(j2,pp1)).lt.big) go to 90
      j2=j2+1
      if (j2.gt.n) go to 90
      go to 80
 90   if (j1.ge.1) go to 100
      k=j2
      go to 130
 100  if (j2.le.n) go to 110
      k=j1
      go to 130
 110  if (t(j)-t(j1).ge.t(j2)-t(j)) go to 120
      k=j1
      go to 130
 120  k=j2
 130  z(j,1)=y(m(k,pp1))
      t(j)=t(k)
 140  continue
      if (iabs(l(pp1)).ne.5) go to 160
      do 150 j=1,n
      f(j)=z(j,1)
 150  continue
      go to 170
 160  call smothr (1,n,t,z,z(1,2),f,z(1,6))
 170  return
      end
      
      subroutine acemod (v,p,n,x,l,tx,f,t,m,yhat)
      implicit none
c          subroutine acemod(v,p,n,x,l,tx,f,t,m,yhat)
c--------------------------------------------------------------------
c
c computes response y estimates from the model
c
c                yhat =  f ( t( v ) )
c
c using the x transformations tx constructed by subroutine ace and
c the predictor function (f,t) constructed by subroutine model.
c
c input:
c
c       v(p) : vector of predictor values.
c    p,n,x,l : same input as for subroutine ace.
c       tx,m : output from subroutine ace.
c        f,t : output from subroutine model.
c
c output:
c
c    yhat : estimated response value for v.
c
c note: this subroutine must not be called before subroutine model.
c
c-------------------------------------------------------------------
c
      integer n,p,m(n,1),l(1),low,high,place
      integer maxit,nterm,i,jh,jl
      double precision alpha, big, span, th, vi, xt
      double precision  v(p),x(p,n),f(n),t(n),tx(n,p), yhat
      common /prams/ alpha,big,span,maxit,nterm
      th=0.0
      do 90 i=1,p
      if (l(i).eq.0) go to 90
      vi=v(i)
      if (vi.lt.big) go to 10
      if (x(i,m(n,i)).ge.big) th=th+tx(m(n,i),i)
      go to 90
 10   if (vi.gt.x(i,m(1,i))) go to 20
      place=1
      go to 80
 20   if (vi.lt.x(i,m(n,i))) go to 30
      place=n
      go to 80
 30   low=0
      high=n+1
 40   if (low+1.ge.high) go to 60
      place=(low+high)/2
      xt=x(i,m(place,i))
      if (vi.eq.xt) go to 80
      if (vi.ge.xt) go to 50
      high=place
      go to 40
 50   low=place
      go to 40
 60   if (iabs(l(i)).eq.5) go to 90
      jl=m(low,i)
      jh=m(high,i)
      if (x(i,jh).lt.big) go to 70
      th=th+tx(jl,i)
      go to 90
 70   th=th+tx(jl,i)+(tx(jh,i)-tx(jl,i))*(vi-x(i,jl))/(x(i,jh)-x(i,jl))
      go to 90
 80   th=th+tx(m(place,i),i)
 90   continue
      if (th.gt.t(1)) go to 100
      yhat=f(1)
      return
 100  if (th.lt.t(n)) go to 110
      yhat=f(n)
      return
 110  low=0
      high=n+1
 120  if (low+1.ge.high) go to 150
      place=(low+high)/2
      xt=t(place)
      if (th.ne.xt) go to 130
      yhat=f(place)
      return
 130  if (th.ge.xt) go to 140
      high=place
      go to 120
 140  low=place
      go to 120
 150  if (iabs(l(p+1)).ne.5) go to 170
      if (th-t(low).gt.t(high)-th) go to 160
      yhat=f(low)
      go to 180
 160  yhat=f(high)
      go to 180
 170  yhat=f(low)+(f(high)-f(low))*(th-t(low))/(t(high)-t(low))
 180  return
      end
      
c
c     block data
c     common /prams/ maxit,nterm,span,alpha,big
c
c------------------------------------------------------------------
c
c these procedure parameters can be changed in the calling routine
c by defining the above labeled common and resetting the values with
c executable statements.
c
c maxit : maximum number of iterations.
c nterm : number of consecutive iterations for which
c         rsq must change less than delcor for convergence.
c span, alpha : super smoother parameters (see below).
c big : a large representable floating point number.
c
c------------------------------------------------------------------
c
      block data acedata
        implicit double precision (A-H,O-Z)
        common /prams/ alpha,big,span,maxit,nterm
        data maxit,nterm,span,alpha,big /20,3,0.0,0.0,1.0e20/
      end

c
c 2016/10/7 Shawn Garbett Refactor to insure initialized variable h, and no division by zero
c
c Note: The original function was named "scale", but this is now part of the Fortran 95 namespace
c       So this was changed to "scail"
c
      subroutine scail (p,n,w,sw,ty,tx,eps,maxit,r,sc)
        
        implicit none
        integer p,n,maxit,i,iter,j,nit
        double precision w(n),ty(n),tx(n,p),r(n),sc(p,5)
        double precision s,h,t,u,gama,delta,sw,eps,v
        
      ! Initialization
      do i=1,p
         sc(i,1)=0.0
      end do
      nit=0
      
      do
        nit=nit+1   ! nit is iteration number
        do i=1,p
           sc(i,5)=sc(i,1)
        end do
        
        h = 1.0 ! Gets rid of unitialized warning
        do iter=1,p
          
          do j=1,n
            s=0.0
            do i=1,p
              s=s+sc(i,1)*tx(j,i)
            end do
            r(j)=(ty(j)-s)*w(j)
          end do
          
          do i=1,p
            s=0.0
            do j=1,n
              s=s+r(j)*tx(j,i)
            end do
            sc(i,2)=-2.0*s/sw
          end do
          
          s=0.0
          do i=1,p
             s=s+sc(i,2)**2
          end do
          
          ! Make sure that h gets initialized with s, and division by zero is not possible
          if (iter.eq.1.or.h.le.0.0) then 
            h = s
          end if
          
          ! Patch to ensure sum of sc(i,2)^2 is not zero
          if (s.le.0.0) exit

          if (iter.eq.1) then
            do i=1,p
              sc(i,3)=-sc(i,2)
            end do
          else
            gama=s/h
            do i=1,p
              sc(i,3)=-sc(i,2)+gama*sc(i,4)
            end do
          end if
          h=s
          
          s=0.0
          t=s
          do j=1,n
            u=0.0
            do i=1,p
              u=u+sc(i,3)*tx(j,i)
            end do
              s=s+u*r(j)
              t=t+w(j)*u**2
          end do
          delta=s/t
          do i=1,p
              sc(i,1)=sc(i,1)+delta*sc(i,3)
              sc(i,4)=sc(i,3)
          end do
        end do ! iter=1,p
        ! if all  sc(i,2) is zero exits to here

        ! Check for convergence, or maximum iteration
        v=0.0
        do i=1,p
          v=max(v,abs(sc(i,1)-sc(i,5)))
        end do
        if (v.lt.eps.or.nit.ge.maxit) exit
      end do ! Main iterator
      
      ! Compute final answer and return
      do i=1,p
         do j=1,n
            tx(j,i)=sc(i,1)*tx(j,i)
         end do
      end do
      return
      end ! scail
