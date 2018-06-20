      subroutine avas(p,n,x,y,w,l,delrsq,tx,ty,rsq,ierr,m,z,yspan,iter,
     1     iters)
      implicit none
      integer n,p,pp1,pp2,m(n,*),l(*),ierr,iter
      double precision y(n),x(n,p),w(n),ty(n),tx(n,p),z(n,17),ct(10)
      double precision iters(100,2), delrsq, rsq, yspan, rss
      double precision sumlog, tres, rr, rnew, cmn, cmx
      common /parms/ span,alpha,itape,maxit,nterm
      integer maxit,nterm,i,j,k,np,nt
      integer itape
      double precision dof, alpha, span
      double precision sm,sv,sw,svx
      ierr = 0
      pp1 = p+1
      pp2 = p + 2
      sm = 0.0
      sv = sm
      sw = sv
      np = 0
      do 23000 i = 1,p
         if(.not.(l(i).gt.0))goto 23002
         np = np+1
23002    continue
23000 continue
      do 23004 j = 1,n
         sm = sm+w(j)*y(j)
         sv = sv+w(j)*y(j)**2
         sw = sw+w(j)
         m(j,pp1) = j
         z(j,2) = y(j)
23004 continue
      sm = sm/sw
      sv = sv/sw-sm**2
      sv = 1.0/dsqrt(sv)
      do 23006 j = 1,n
         z(j,1) = (y(j)-sm)*sv
23006 continue
      call sort(z(1,2),m(1,pp1),1,n)
      do 23008 i = 1,p
         if(.not.(l(i) .gt. 0))goto 23010
         sm=0
         do 23012 j=1,n
            sm=sm+w(j)*x(j,i)
23012    continue
         sm=sm/sw
         do 23014 j = 1,n
            m(j,i) = j
            z(j,2) = x(j,i)
23014    continue
         call sort(z(1,2),m(1,i),1,n)
23010    continue
23008 continue
      rsq = 0.0
      iter = 0
      nterm = min0(nterm,10)
      nt = 0
      do 23016 i = 1,nterm
         ct(i) = 100.0
23016 continue
      do 23018 j=1,n
         ty(j)=z(j,1)
23018 continue
      do 23020 j = 1,n
         z(j,9)=ty(j)
23020 continue
      call bakfit(iter,delrsq,rsq,sw,l,z,m,x,z(1,9),tx,w,n,p,np)
      sumlog=0
23022 continue
      iter=iter +1
      if(.not.(l(pp1).eq.4))goto 23025
      go to 992
23025 continue
      call calcmu(n,p,l,z,tx)
      do 23027 j=1,n
         tres=(ty(j)-z(j,10))
         if(.not.(abs(tres).lt.1e-10))goto 23029
         tres=1e-10
23029    continue
         z(j,2)=log(sqrt(tres**2))
         m(j,pp2)=j
23027 continue
      call sort(z(1,10),m(1,pp2),1,n)
      do 23031 j=1,n
         k=m(j,pp2)
         z(j,4)=z(k,2)
         z(j,5)=w(k)
23031 continue
      call rlsmo(z(1,10),z(1,4),z(1,5),yspan,dof,n,z(1,6),rss,z(1,7))
      do 23033 j=1,n
         k=m(j,pp2)
         z(j,7)=exp(-z(j,6))
         sumlog=sumlog+n*(w(j)/sw)*2*z(j,6)
         z(j,8)=ty(k)
23033 continue
      call ctsub(n,z(1,10),z(1,7),z(1,8),z(1,9))
      sm=0
      do 23035 j=1,n
         sm=sm+w(j)*z(j,9)
23035 continue
      do 23037 j=1,n
         k=m(j,pp2)
         ty(k)=z(j,9)-sm/sw
23037 continue
      sv=0
      svx=0
      do 23039 j=1,n
         sv=sv+(w(j)/sw)*ty(j)*ty(j)
         svx=svx+(w(j)/sw)*z(j,10)*z(j,10)
23039 continue
      do 23041 j=1,n
         ty(j)=ty(j)/dsqrt(sv)
         do 23043 i=1,p
            if(.not.( l(i) .gt. 0))goto 23045
            tx(j,i)=tx(j,i)/dsqrt(svx)
23045       continue
23043    continue
23041 continue
 992  continue
      do 23047 j = 1,n
         z(j,9)=ty(j)
23047 continue
      call bakfit(iter,delrsq,rsq,sw,l,z,m,x,z(1,9),tx,w,n,p,np)
      sumlog=sumlog+n*dlog(sv)
      rr=0
      call calcmu(n,p,l,z,tx)
      do 23049 j=1,n
       rr=rr+(w(j)/sw)*(ty(j)-z(j,10))**2
23049 continue
      rsq=1-rr
      rnew=sumlog+rr
      iters(iter,1)=iter
      iters(iter,2)=rsq
      nt = mod(nt,nterm)+1
      ct(nt) = rsq
      cmn = 100.0
      cmx = -100.0
      do 23051 i = 1,nterm
         cmn = min(cmn,ct(i))
         cmx = max(cmx,ct(i))
23051 continue
      if(.not.(cmx-cmn.le.delrsq.or.iter.ge.maxit.or.l(pp1).eq.4))
     1     goto 23053
      return
23053 continue
      goto 23022
      return
      end


      subroutine calcmu(n,p,l,z,tx)
      implicit none
      integer p, l(*),j,n
      double precision z(n,17),tx(n,p)
      integer i
      do 23055 j=1,n
         z(j,10)=0
         do 23057 i=1,p
            if(.not.(l(i) .gt. 0))goto 23059
            z(j,10)=z(j,10)+tx(j,i)
23059       continue
23057    continue
23055 continue
      return
      end


      subroutine bakfit(iter,delrsq,rsq,sw,l,z,m,x,ty,tx,w,n,p,np)
      implicit none
      integer n,l(*),m(n,*),p,j,k,nit,i,np,iter
      double precision z(n,17),ty(n),tx(n,p),x(n,p),w(n)
      double precision sm,sv,sw, delrsq, rsq, rsqi
      double precision alpha,span
      integer itape,maxit,nterm
      common /parms/ span,alpha,itape,maxit,nterm
      call calcmu(n,p,l,z,tx)
      do 23061 j=1,n
         ty(j)=ty(j)-z(j,10)
23061 continue
      nit=0
23063 continue
      rsqi = rsq
      nit = nit+1
      do 23066 i = 1,p
         if(.not.(l(i).gt.0))goto 23068
         do 23070 j = 1,n
            k = m(j,i)
            z(j,1) = ty(k)+tx(k,i)
            z(j,2) = x(k,i)
            z(j,7) = w(k)
23070    continue
         call smothr(l(i),n,z(1,2),z,z(1,7),z(1,6),z(1,11))
         sm = 0.0
         do 23072 j = 1,n
            sm = sm+z(j,7)*z(j,6)
23072    continue
         sm = sm/sw
         do 23074 j = 1,n
            z(j,6) = z(j,6)-sm
23074    continue
         sv = 0.0
         do 23076 j = 1,n
            sv = sv+z(j,7)*(z(j,1)-z(j,6))**2
23076    continue
         sv = 1.0-sv/sw
         rsq = sv
         do 23078 j = 1,n
            k = m(j,i)
            tx(k,i) = z(j,6)
            ty(k) = z(j,1)-z(j,6)
23078    continue
23068    continue
23066 continue
      if(.not.(np.eq.1.or.abs(rsq-rsqi).le.delrsq.or.nit.ge.maxit))
     1     goto 23063
      if(.not.(rsq.eq.0.0.and.iter.eq.0))goto 23080
      do 23082 i = 1,p
         if(.not.(l(i).gt.0))goto 23084
         do 23086 j = 1,n
            tx(j,i) = x(j,i)
23086    continue
23084    continue
23082 continue
23080 continue
      return
      end


      subroutine ctsub(n,u,v,y,ty)
      implicit none
      double precision u(*),v(*),y(*),ty(*)
      integer i,j,n
      i=1
23088 if(.not.(i.le.n))goto 23090
      if(.not.(y(i).le.u(1)))goto 23091
      ty(i)=(y(i)-u(1))*v(1)
      goto 23092
23091 continue
      j=1
      ty(i)=0
23093 if(.not.((j.le.n) .and. (y(i).gt.u(j)) ))goto 23094
      if(.not.(j .gt. 1))goto 23095
      ty(i)=ty(i)+(u(j)-u(j-1))*(v(j)+v(j-1))/2
23095 continue
      j=j+1
      goto 23093
23094 continue
      if(.not.(y(i).le.u(n)))goto 23097
      ty(i)=ty(i)+.5*(y(i)-u(j-1))*(2*v(j-1)+(y(i)-u(j-1))*(v(j)-v(j-1))
     1     /(u(j)-u(j-1)))
      goto 23098
23097 continue
      ty(i)=ty(i)+(y(i)-u(n))*v(n)
23098 continue
23092 continue
      i=i+1
      goto 23088
23090 continue
      return
      end

      block data avasdata
      integer itape, maxit, nterm
      double precision span, alpha
      double precision big, sml, eps
      double precision spans(3)
      common /parms/ span,alpha,itape,maxit,nterm
      common /spans/ spans /consts/ big,sml,eps
c------------------------------------------------------------------
c
c these procedure parameters can be changed in the calling routine
c by defining the above labeled common and resetting the values with
c executable statements.
c
c itape : fortran file number for printer output.
c         (itape.le.0 => no printer output.)
c maxit : maximum number of iterations.
c nterm : number of consecutive iterations for which
c         rsq must change less than delcor for convergence.
c span, alpha : super smoother parameters.
c   (see - friedman and stuetzle, reference above.)
c
c------------------------------------------------------------------
      data itape,maxit,nterm,span,alpha /-6,20,3,0.0,5.0/
c---------------------------------------------------------------
c
c this sets the compile time (default) values for various
c internal parameters :
c
c spans : span values for the three running linear smoothers.
c spans(1) : tweeter span.
c spans(2) : midrange span.
c spans(3) : woofer span.
c (these span values should be changed only with care.)
c big : a large representable floating point number.
c sml : a small number. should be set so that (sml)**(10.0) does
c       not cause floating point underflow.
c eps : used to numerically stabilize slope calculations for
c       running linear fits.
c
c these parameter values can be changed by declaring the
c relevant labeled common in the main program and resetting
c them with executable statements.
c
c-----------------------------------------------------------------
      data spans,big,sml,eps /0.05,0.2,0.5,1.0e20,1.0e-4,1.0e-3/
      end

      subroutine smothr (l,n,x,y,w,smo,scr)
      implicit none
      integer n
      double precision x(n),y(n),w(n),smo(n),scr(n,7)
      common /parms/ span,alpha,itape,maxit,nterm
      double precision alpha,span
      integer itape,maxit,nterm
      double precision sm,sw,a,b,d
      integer i,j,j0,l
      if (l.lt.5) go to 50
      j=1
 10   j0=j
      sm=w(j)*y(j)
      sw=w(j)
      if (j.ge.n) go to 30
 20   if (x(j+1).gt.x(j)) go to 30
      j=j+1
      sm=sm+w(j)*y(j)
      sw=sw+w(j)
      if (j.ge.n) go to 30
      go to 20
 30   sm=sm/sw
      do 40 i=j0,j
         smo(i)=sm
 40   continue
      j=j+1
      if (j.gt.n) go to 250
      go to 10
 50   if (l.ne.4) go to 80
      sm=0.0
      sw=sm
      b=sw
      d=b
      do 60 j=1,n
         sm=sm+w(j)*x(j)*y(j)
         sw=sw+w(j)*x(j)**2
         b=b+w(j)*x(j)
         d=d+w(j)
 60   continue
      a=sm/(sw-(b**2)/d)
      b=b/d
      do 70 j=1,n
         smo(j)=a*(x(j)-b)
 70   continue
      go to 250
 80   call supsmu (n,x,y,w,l,span,alpha,smo,scr)
      if (l.ne.3) go to 250
      do 90 j=1,n
         scr(j,1)=smo(j)
         scr(n-j+1,2)=scr(j,1)
 90   continue
      call montne (scr,n)
      call montne (scr(1,2),n)
      sm=0.0
      sw=sm
      do 100 j=1,n
         sm=sm+(smo(j)-scr(j,1))**2
         sw=sw+(smo(j)-scr(n-j+1,2))**2
 100  continue
      if (sm.ge.sw) go to 120
      do 110 j=1,n
         smo(j)=scr(j,1)
 110  continue
      go to 140
 120  do 130 j=1,n
         smo(j)=scr(n-j+1,2)
 130  continue
 140  j=1
 150  j0=j
      if (j.ge.n) go to 170
 160  if (smo(j+1).ne.smo(j)) go to 170
      j=j+1
      if (j.ge.n) go to 170
      go to 160
 170  if (j.le.j0) go to 190
      a=0.0
      if (j0.gt.1) a=0.5*(smo(j0)-smo(j0-1))
      b=0.0
      if (j.lt.n) b=0.5*(smo(j+1)-smo(j))
      d=(a+b)/(j-j0)
      if (a.eq.0.0.or.b.eq.0.0) d=2.0*d
      if (a.eq.0.0) a=b
      do 180 i=j0,j
         smo(i)=smo(i)-a+d*(i-j0)
 180  continue
 190  j=j+1
      if (j.gt.n) go to 200
      go to 150
 200  j=1
 210  j0=j
      sm=smo(j)
      if (j.ge.n) go to 230
 220  if (x(j+1).gt.x(j)) go to 230
      j=j+1
      sm=sm+smo(j)
      if (j.ge.n) go to 230
      go to 220
 230  sm=sm/(j-j0+1)
      do 240 i=j0,j
         smo(i)=sm
 240  continue
      j=j+1
      if (j.gt.n) go to 250
      go to 210
 250  return
      end


      subroutine montne (x,n)
      double precision x(n)
      double precision pmn
      integer bb,eb,br,er,bl,el
      bb=0
      eb=bb
 10   if (eb.ge.n) go to 110
      bb=eb+1
      eb=bb
 20   if (eb.ge.n) go to 30
      if (x(bb).ne.x(eb+1)) go to 30
      eb=eb+1
      go to 20
 30   if (eb.ge.n) go to 70
      if (x(eb).le.x(eb+1)) go to 70
      br=eb+1
      er=br
 40   if (er.ge.n) go to 50
      if (x(er+1).ne.x(br)) go to 50
      er=er+1
      go to 40
 50   pmn=(x(bb)*(eb-bb+1)+x(br)*(er-br+1))/(er-bb+1)
      eb=er
      do 60 i=bb,eb
         x(i)=pmn
 60   continue
 70   if (bb.le.1) go to 10
      if (x(bb-1).le.x(bb)) go to 10
      bl=bb-1
      el=bl
 80   if (bl.le.1) go to 90
      if (x(bl-1).ne.x(el)) go to 90
      bl=bl-1
      go to 80
 90   pmn=(x(bb)*(eb-bb+1)+x(bl)*(el-bl+1))/(eb-bl+1)
      bb=bl
      do 100 i=bb,eb
         x(i)=pmn
 100  continue
      go to 30
 110  return
      end


      subroutine sort (v,a,ii,jj)
c     
c     puts into a the permutation vector which sorts v into
c     increasing order.  only elements from ii to jj are considered.
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements
c     
c     this is a modification of cacm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c     
      implicit none
      integer jj
      dimension a(jj),v(*)
      integer iu(20),il(20)
      integer t,tt,ij,j,k,l
      integer a,ii,m,i
      double precision v,vt,vtt
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=a(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   a(l)=a(k)
      a(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=a(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=a(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  a(k+1)=a(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      a(k+1)=t
      v(k+1)=vt
      go to 100
      end


      subroutine supsmu (n,x,y,w,iper,span,alpha,smo,sc)
c------------------------------------------------------------------
c     
c     super smoother (friedman and stuetzle, 1984).
c     
c     version 3/10/84
c     
c     coded by: j. h. friedman
c     department of statistics and
c     stanford linear accelerator center
c     stanford university
c     stanford ca. 94305
c     
c     input:
c     n : number of observations (x,y - pairs).
c     x(n) : ordered abscissa values.
c     y(n) : corresponding ordinate (response) values.
c     w(n) : weight for each (x,y) observation.
c     iper : periodic variable flag.
c     iper=1 => x is ordered interval variable.
c     iper=2 => x is a periodic variable with values
c     in the range (0.0,1.0) and peroid 1.0.
c     span : smoother span (fraction of observations in window).
c     span=0.0 => automatic (variable) span selection.
c     alpha : controles high frequency (small span) penality
c     used with automatic span selection (base tone control).
c     (alpha.le.0.0 or alpha.gt.10.0 => no effect.)
c     output:
c     smo(n) : smoothed ordinate (response) values.
c     scratch:
c     sc(n,7) : internal working storage.
c     
c     note:
c     for small samples (n < 40) or if there are substantial serial
c     correlations between obserations close in x - value, then
c     a prespecified fixed span smoother (span > 0) should be
c     used. reasonable span values are 0.3 to 0.5.
c     
c------------------------------------------------------------------
      implicit none
      integer n
      double precision x(n),y(n),w(n),smo(n),sc(n,7)
      double precision big,sml,eps
      double precision span,alpha
      integer i,j,jper,iper
      double precision spans(3)
      common /spans/ spans /consts/ big,sml,eps
      double precision h(1),sw,sy,a,scale,vsmlsq,resmin,f
      if (x(n).gt.x(1)) go to 30
      sy=0.0
      sw=sy
      do 10 j=1,n
         sy=sy+w(j)*y(j)
         sw=sw+w(j)
 10   continue
      a=sy/sw
      do 20 j=1,n
         smo(j)=a
 20   continue
      return
 30   i=n/4
      j=3*i
      scale=x(j)-x(i)
 40   if (scale.gt.0.0) go to 50
      if (j.lt.n) j=j+1
      if (i.gt.1) i=i-1
      scale=x(j)-x(i)
      go to 40
 50   vsmlsq=(eps*scale)**2
      jper=iper
      if (iper.eq.2.and.(x(1).lt.0.0.or.x(n).gt.1.0)) jper=1
      if (jper.lt.1.or.jper.gt.2) jper=1
      if (span.le.0.0) go to 60
      call smooth (n,x,y,w,span,jper,vsmlsq,smo,sc)
      return
 60   do 70 i=1,3
         call smooth (n,x,y,w,spans(i),jper,vsmlsq,sc(1,2*i-1),sc(1,7))
         call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2*i),h)
 70   continue
      do 90 j=1,n
         resmin=big
         do 80 i=1,3
            if (sc(j,2*i).ge.resmin) go to 80
            resmin=sc(j,2*i)
            sc(j,7)=spans(i)
 80      continue
         if (alpha.gt.0.0.and.alpha.le.10.0.and.resmin.lt.sc(j,6))
     1     sc(j,7) = sc(j,7) +
     2        (spans(3)-sc(j,7))*max(sml,resmin/sc(j,6))**
     3        (10.0-alpha)
 90   continue
      call smooth (n,x,sc(1,7),w,spans(2),-jper,vsmlsq,sc(1,2),h)
      do 110 j=1,n
         if (sc(j,2).le.spans(1)) sc(j,2)=spans(1)
         if (sc(j,2).ge.spans(3)) sc(j,2)=spans(3)
         f=sc(j,2)-spans(2)
         if (f.ge.0.0) go to 100
         f=-f/(spans(2)-spans(1))
         sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,1)
         go to 110
 100     f=f/(spans(3)-spans(2))
         sc(j,4)=(1.0-f)*sc(j,3)+f*sc(j,5)
 110  continue
      call smooth (n,x,sc(1,4),w,spans(1),-jper,vsmlsq,smo,h)
      return
      end

      subroutine smooth (n,x,y,w,span,iper,vsmlsq,smo,acvr)
      implicit none
      integer n
      double precision  x(n),y(n),w(n),smo(n),acvr(n),vsmlsq
      double precision  xti,wt,fbw,xm,ym,var,cvar,fbo,tmp,xto,a,h,sy
      double precision  span
      integer in,out,ibw
      integer i,j,j0,it,jper,iper
      xm=0.0
      ym=xm
      var=ym
      cvar=var
      fbw=cvar
      jper=iabs(iper)
      ibw=int(0.5*span*n+0.5)
      if (ibw.lt.2) ibw=2
      it=2*ibw+1
      do 20 i=1,it
         j=i
         if (jper.eq.2) j=i-ibw-1
         xti=x(j)
         if (j.ge.1) go to 10
         j=n+j
         xti=x(j)-1.0
 10      wt=w(j)
         fbo=fbw
         fbw=fbw+wt
         xm=(fbo*xm+wt*xti)/fbw
         ym=(fbo*ym+wt*y(j))/fbw
         tmp=0.0
         if (fbo.gt.0.0) tmp=fbw*wt*(xti-xm)/fbo
         var=var+tmp*(xti-xm)
         cvar=cvar+tmp*(y(j)-ym)
 20   continue
      do 70 j=1,n
         out=j-ibw-1
         in=j+ibw
         if ((jper.ne.2).and.(out.lt.1.or.in.gt.n)) go to 60
         if (out.ge.1) go to 30
         out=n+out
         xto=x(out)-1.0
         xti=x(in)
         go to 50
 30      if (in.le.n) go to 40
         in=in-n
         xti=x(in)+1.0
         xto=x(out)
         go to 50
 40      xto=x(out)
         xti=x(in)
 50      wt=w(out)
         fbo=fbw
         fbw=fbw-wt
         tmp=0.0
         if (fbw.gt.0.0) tmp=fbo*wt*(xto-xm)/fbw
         var=var-tmp*(xto-xm)
         cvar=cvar-tmp*(y(out)-ym)
         xm=(fbo*xm-wt*xto)/fbw
         ym=(fbo*ym-wt*y(out))/fbw
         wt=w(in)
         fbo=fbw
         fbw=fbw+wt
         xm=(fbo*xm+wt*xti)/fbw
         ym=(fbo*ym+wt*y(in))/fbw
         tmp=0.0
         if (fbo.gt.0.0) tmp=fbw*wt*(xti-xm)/fbo
         var=var+tmp*(xti-xm)
         cvar=cvar+tmp*(y(in)-ym)
 60      a=0.0
         if (var.gt.vsmlsq) a=cvar/var
         smo(j)=a*(x(j)-xm)+ym
         if (iper.le.0) go to 70
         h=1.0/fbw
         if (var.gt.vsmlsq) h=h+(x(j)-xm)**2/var
         acvr(j)=abs(y(j)-smo(j))/(1.0-w(j)*h)
 70   continue
      j=1
 80   j0=j
      sy=smo(j)*w(j)
      fbw=w(j)
      if (j.ge.n) go to 100
 90   if (x(j+1).gt.x(j)) go to 100
      j=j+1
      sy=sy+w(j)*smo(j)
      fbw=fbw+w(j)
      if (j.ge.n) go to 100
      go to 90
 100  if (j.le.j0) go to 120
      sy=sy/fbw
      do 110 i=j0,j
         smo(i)=sy
 110  continue
 120  j=j+1
      if (j.gt.n) go to 130
      go to 80
 130  return
      end
