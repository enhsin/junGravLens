 subroutine fastelldefl(x1in,x2in,q,gam,arat,s2,
     *     defl)
c
c     This routine calculates the deflection due to an elliptical
c     mass distribution, quickly and accurately.
c     The parameters are position (x1in,x2in), overall factor
c     (q), power (gam) which should be between -1 and 2, axis ratio
c     (arat) which is <=1, core radius squared (s2), and the output
c     two-component deflection (defl).
c     The projected mass density distribution, in units of the
c     critical density, is kappa(x1,x2)=q [u2+s2]^(-gam), where
c     u2=[x1^2+x2^2/(arat^2)].
c

      implicit double precision (a-h,o-z)
      double precision defl(2),ares(2)
      common /ellfirsttime/ ifirst1,ifirst2
      if ((ifirst1 .eq. 0) .and. (ifirst2 .eq. 0)) then
         write(*,*) 'Initializing FASTELL routines'
         call fastellprep()
         ifirst1=1
         ifirst2=1
      endif         
      if (arat .gt. (1.d0+1.d-10)) then
         write(*,*) 'Error: axis ratio is greater than 1'
         call setzerod(defl)
         return
      endif
      if (arat .lt. -1.d-10) then
         write(*,*) 'Error: axis ratio is less than 0'
         call setzerod(defl)
         return
      endif
      if (s2 .lt. -1.d-10) then
         write(*,*) 'Error: core radius squared is less than 0'
         call setzerod(defl)
         return
      endif
      if ((gam .lt. -1.01d0).or.(gam .gt. 2.01d0)) then
         write(*,*) 'Warning: power gam= ',gam
         write(*,*) 'is outside tested range (-1 to 2) in fastelldefl'
         inrange=0
      else
         inrange=1
      endif
      if (x1in .lt. 0.d0) then
         sgnx1=-1.d0
         x1=-x1in
      else
         sgnx1=1.d0
         x1=x1in
      endif
      if (x2in .lt. 0.d0) then
         sgnx2=-1.d0
         x2=-x2in
      else
         sgnx2=1.d0
         x2=x2in
      endif
      facx=1.d0
      facy=1.d0
c     At the center, we set deflection and magmx to zero. 
      if ((x1+x2) .le. 0.d0) then
         call setzerod(defl)
         return
      endif
c     Here we deal with the nearly-spherical case:
      if (arat .gt. .998d0) then
         r2=x1*x1+x2*x2
         rho2=x1*x1+x2*x2/(arat*arat)
         g1=1.d0-gam
         g2=2.d0-gam
         t6=rho2+s2
         t1=t6**g1
         t2=s2**g1
         gsm=3.d-4
         if (dabs(g1) .lt. gsm) then
            call helplog(s2,t6,1.d0,g1,gsm,t8)
         else
            t8=(t1-t2)/g1
         endif
         if (dabs(g2) .lt. gsm) then
            call helplog(s2,t6,1.d0,g2,gsm,t8p)
         else
            t8p=(t1*t6-t2*s2)/g2
         endif
c         t3=(t1*(rho2-s2/g1)+t2*s2/g1)/(1.d0+g1)
         t3=t8p-s2*t8
         t4=x1*x1-3.d0*x2*x2
         t5=3.d0*x1*x1-x2*x2
         a1=t8
         a2=t3*(1.d0-arat*arat)/(2.d0*r2*r2)
         defl(1)=x1*sgnx1*(a1+t4*a2)*arat*(q/r2)
         defl(2)=x2*sgnx2*(a1+t5*a2)*arat*(q/r2)
         return
      endif
c     If the power gam is a half integer, change it slightly 
c     to avoid division by zero in some cases below. But if
c     gam is in the tested range this is not needed.
      g2n=dnint(2.d0*gam)
      if ((inrange .eq. 0).and.(dabs(2.d0*gam-g2n).lt.2.d-8)) then
         gin=(g2n+2.d-8)/2.d0
      else
         gin=gam
      endif
c     Here we deal with extreme ratios of x1/x2:
      d1=8.d-6*x1
      d2=8.d-6*x2
      if (x2 .lt. d1) then
         facy=x2/d1
         x2=d1
      elseif (x1 .lt. d2) then
         facx=x1/d2
         x1=d2
      endif
      r=x1/x2
      sn=(1.d0-arat*arat)/2.d0
      e1=s2*sn/(x1*x2)
c     In case of an extremely small axis ratio 
c     (arat), we make a slight change:
      if (arat .lt. 1.d-4) then
         e2e1=sn*(r+1.d8/r)
      else
         e2e1=sn*(r+1.d0/(r*arat*arat))
      endif
      e2=e1+e2e1
      e1sb=(1.d0/r-r)/2.d0
      sb=e1-e1sb
      e1g1=e1**(1.d0-gin)
      e2g1=e2**(1.d0-gin)
      erat=e2e1/(e1+e2)
      if (erat.lt.1.d-7) then
         eave=e1+e2e1/2.d0
         epow=eave**(-gin)
         e1in=1.d-2
         e2in=e2e1+1.d-2
         sbin=1.d-2-e1sb
         gamin=0.d0
         e1g1in=e1in
         e2g1in=e2in
         ismall=1
      else
         e1in=e1
         e2in=e2
         sbin=sb
         e1g1in=e1g1
         e2g1in=e2g1
         gamin=gin
         ismall=0
      endif
      call approxintd(ares,gamin,sbin,e1in,e2in,e1g1in,e2g1in)
      if (ismall .eq. 1) then
         ares(1)=ares(1)*epow
         ares(2)=ares(2)*epow
      endif
      quan=sn/(x1*x2)
      xpr=x1*x2
      t1=q*(quan**gin)*dsqrt(xpr)*arat/(1.d0-arat*arat)
      defl(1)=facx*sgnx1*ares(1)*t1
      defl(2)=facy*sgnx2*ares(2)*t1
      return
      end
