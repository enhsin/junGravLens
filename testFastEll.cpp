#include <iostream> 
#define MaxN 50 
#define NN 49
using namespace std; 

/*      subroutine fastellprep()
c
c     Routine to precalculate expansion coefficients
c     which are needed in approxint() and approxintd(). The 
c     coefficients are passed on through common blocks. Time 
c     is saved by only calculating these once.
c
 */


/*
    SUBROUTINE CHEBCOMB(A,B,P,M,F)
c
c     Finds a polynomial P of degree M-1 which
c     approximates F on the interval [A,B],
c     using Chebyshev polynomials.
c     
c     Based on the free code 
c        www.netlib.org/textbook/mathews/chap4.f
c     distributed in the Netlib archive.
c
      implicit double precision (a-h,o-z)
      PARAMETER(MaxN=50,NN=MaxN-1)
      double precision C(MaxN),Y(MaxN),F,
     *     T(MaxN,MaxN),P(MaxN),PT(MAXN),
     *     CT(MaxN)
      EXTERNAL F
      N=M-1
      AL=2.d0/(B-A)
      BL=-(B+A)/(B-A)
      D=3.1415926535897932d0/(2.d0*NN+2.d0)
      DO K=0,NN
        X=DCOS((2.d0*K+1.d0)*D)
        Y(K+1)=F((X-BL)/AL)
        C(K+1)=0.d0
      ENDDO
      DO K=0,NN
        Z=(2.d0*K+1.d0)*D
        DO J=0,NN
          C(J+1)=C(J+1)+Y(K+1)*DCOS(J*Z)
        ENDDO
      ENDDO
      C(1)=C(1)/(NN+1.d0)
      DO J=1,NN
        C(J+1)=2.d0*C(J+1)/(NN+1.d0)
      ENDDO
      DO K=0,N
        DO J=0,N
          T(K+1,J+1)=0.d0
        ENDDO
      ENDDO
      T(1,1)=1.d0
      T(2,1)=0.d0
      T(2,2)=1.d0
      DO K=2,N
        DO J=1,K
          T(K+1,J+1)=2.d0*T(K,J)-T(K-1,J+1)
        ENDDO
        T(K+1,1)=-T(K-1,1)
      ENDDO
      DO J=0,N
        Sum=0.d0
        DO K=0,N
          Sum=Sum+C(K+1)*T(K+1,J+1)
        ENDDO
        PT(J+1)=Sum        
        P(J+1)=0.d0
        CT(J+1)=0.d0
      ENDDO
      CT(1)=1.d0
      P(1)=P(1)+PT(1)*CT(1)
      DO J=2,M
         DO I=J,2,-1
            CT(I)=CT(I)*BL+CT(I-1)*AL
         ENDDO
         CT(1)=CT(1)*BL
         DO I=1,J
            P(I)=P(I)+PT(J)*CT(I)
         ENDDO
      ENDDO
      RETURN
      END
*/
/*
void fastellprep() {

   double cp[50],ellipfn,ellprepf3,ellprepf4,ellprepf5,ellipdfn,cp2[50],
     ac[45],bc[45],ec[10],fc[10],rcp[45],rc[45],
          oc[28],pc2[70],acp[45],ocp[36],tc2p[45],uc[28],
          bcp[55],pc2p[70],qc[21],qcp[28],fc2[10],cc[10],
          sc[30],sc2[30],scp[35],pc[28],pcp[32],sc2p[35],
          gc[45],gc2[45],hc2[45],gcp[45],gc2p[45],ucp[28],
          hcp[50],hc2p[50],hc[45],tc[45],tc2[45],tcp[45],
          vc[35],vc2[35],vcp[35],vc2p[35],wc[45],wcp[45],
          xc[45],xc2[45],xcp[45],xc2p[45]; 

      h=3.0; 
      bnd=0.42; 
      d=2.50; 
      dbig=6.650 ; 
      dhuge=15.70 ; 
       x=-dhuge ; 

       chebcomb(x, -dbig, cp, 7, ellipfn); 


}
*/
/*
      implicit double precision (a-h,o-z)
      double precision cp(50),ellipfn,ellprepf3,ellprepf4,
     *     ellprepf5,ellipdfn,cp2(50),
     *     ac(45),bc(45),ec(10),fc(10),rcp(45),rc(45),
     *     oc(28),pc2(70),acp(45),ocp(36),tc2p(45),uc(28),
     *     bcp(55),pc2p(70),qc(21),qcp(28),fc2(10),cc(10),
     *     sc(30),sc2(30),scp(35),pc(28),pcp(32),sc2p(35),
     *     gc(45),gc2(45),hc2(45),gcp(45),gc2p(45),ucp(28),
     *     hcp(50),hc2p(50),hc(45),tc(45),tc2(45),tcp(45),
     *     vc(35),vc2(35),vcp(35),vc2p(35),wc(45),wcp(45),
     *     xc(45),xc2(45),xcp(45),xc2p(45)
      common /fastellpars/ d,h,bnd,dbig,dhuge
      common /fastellcom/ bc,ac,cc,ec,fc,gc,hc,oc,pc,gc2,hc2,
     *     pc2,acp,gcp,gc2p,ocp,pcp,bcp,hcp,hc2p,pc2p,qc,qcp,
     *     rc,rcp,fc2,sc,sc2,scp,sc2p,tc,tc2,tcp,tc2p,uc,ucp,
     *     vc,vc2,vcp,vc2p,wc,wcp,xc,xc2,xcp,xc2p
      external ellipfn,ellprepf3,ellprepf4,ellprepf5,ellipdfn
     
      call chebcomb(x,-dbig,cp,7,ellipfn)
      d2=(-dhuge-dbig)/2.d0
      f1=-dhuge-d2
      f2=-dbig-d2
      call hprep7(cp,uc,d2,vc,vc2,f1,f2,pc2,45,60)
      call chebcomb(x,-dbig,cp2,7,ellipdfn)
      call hprep7(cp2,ucp,d2,vcp,vc2p,f1,f2,pc2p,45,60)
      x=-dbig
      call chebcomb(x,-d,cp,6,ellipfn)
      d2=(-dbig-d)/2.d0
      f1=-dbig-d2
      f2=-d-d2
      call hprep6(cp,qc,d2,sc,sc2,f1,f2,pc2,25,40)
      call chebcomb(x,-d,cp2,7,ellipdfn)
      call hprep7(cp2,qcp,d2,scp,sc2p,f1,f2,pc2p,25,40)
      x=-d
      call chebcomb(x,-bnd,cp,9,ellipfn)
      d2=(-d-bnd)/2.d0
      f1=-d-d2
      f2=-bnd-d2
      call hprep9(cp,ac,d2,gc,gc2,f1,f2,pc2,5,20)
      call chebcomb(x,-bnd,cp2,9,ellipdfn)
      call hprep9(cp2,acp,d2,gcp,gc2p,f1,f2,pc2p,5,20)
      call chebcomb(-bnd,bnd,cp,7,ellipfn)
      f1=bnd
      f2=f1*f1
      call hprepm7(cp,oc,pc,f1,f2,pc2,1)
      call chebcomb(-bnd,bnd,cp2,8,ellipdfn)
      call hprepm8(cp2,ocp,pcp,f1,f2,pc2p,1)
      x=bnd
      call chebcomb(x,d,cp,9,ellipfn)
      d2=(bnd+d)/2.d0
      f1=-d+d2
      f2=-bnd+d2
      call hprep9(cp,bc,d2,hc,hc2,f1,f2,pc2,15,10)
      call chebcomb(x,d,cp2,10,ellipdfn)
      do i=1,10
         cp2(i)=-cp2(i)
      enddo
      call hprep10(cp2,bcp,d2,hcp,hc2p,f1,f2,pc2p,15,10)
      x=d
      call chebcomb(x,dbig,cp,9,ellipfn)
      d2=(d+dbig)/2.d0
      f1=-dbig+d2
      f2=-d+d2
      call hprep9(cp,rc,d2,tc,tc2,f1,f2,pc2,35,30)
      call chebcomb(x,dbig,cp2,9,ellipdfn)
      do i=1,9
         cp2(i)=-cp2(i)
      enddo
      call hprep9(cp2,rcp,d2,tcp,tc2p,f1,f2,pc2p,35,30)
      x=dbig
      call chebcomb(x,dhuge,cp,9,ellipfn)
      d2=(dbig+dhuge)/2.d0
      f1=-dhuge+d2
      f2=-dbig+d2
      call hprep9(cp,wc,d2,xc,xc2,f1,f2,pc2,55,50)
      call chebcomb(x,dhuge,cp2,9,ellipdfn)
      do i=1,9
         cp2(i)=-cp2(i)
      enddo
      call hprep9(cp2,wcp,d2,xcp,xc2p,f1,f2,pc2p,55,50)
      x=1.d0/h
      d2=.81d0-1.d0/h
      call chebcomb(x,x+d2,cp,9,ellprepf3)
      do i=1,9
         cc(i)=cp(i)
      enddo
      x=.4d0/h
      d2=3.5d0-.4d0/h
      call chebcomb(x,x+d2,cp,10,ellprepf4)
      do i=1,10
         ec(i)=cp(i)
      enddo
      x=1.2d0
      d2=.5d0
      call chebcomb(x,x+d2,cp,9,ellprepf5)
      do i=1,9
         fc(i)=cp(i)
      enddo
      x=1.7d0
      d2=4.d0-x
      call chebcomb(x,x+d2,cp,10,ellprepf5)
      do i=1,10
         fc2(i)=cp(i)
      enddo
      return
      end

*/


/*
void fastelldefl_(double *x1in, double *x2in, double *q, double *gam, double *arat, 
   double *s2, double defl[2]) {
   /*
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
*/

/*
   int inrange = 0; 
   double sgnx1, sgnx2, x1, x2, facx, facy; 


   if (arat > 1.0 || artt < 0 || s2 < 0   )  {
    
      cout << "Axis ratio is large than 1.0  or smaller than 0 or core square smaller than 0" << endl; 
      defl[0] = 0 ; 
      defl[1] = 0 ;   
      return ; 
   }

   if (gam < -1.0 || gam > 2.0) 
      inrange = 0; 
   else 
      inrange = 1; 

   if(x1in <0.0) {
      sgnx1 = -1.0; 
      x1 = -x1in  ; 
   }
   else {
      gnx1 = 1.0; 
      x1 = x1in  ; 
   }  

   if(x2in <0.0) {
      sgnx2 = -1.0; 
      x2 = -x2in  ; 
   }
   else {
      gnx2 = 1.0; 
      x2 = x2in  ; 
   }  

   facx = 1.0; 
   facy = 1.0; 
   if((x1+x2) < 0.0 ) {
      defl[0] = 0 ; 
      defl[1] = 0 ;   
      return ; 
   }




}





*/
extern void __stdcall fastelldefl_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2]);
int main() {



   cout << "hello world!" << endl; 
   return 0; 
}



/*
 subroutine fastelldefl(x1in,x2in,q,gam,arat,s2,
     *     defl)
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


*/