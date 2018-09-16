
      subroutine gen_zgamjj(r,p,wt,*) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h' 
      include 'mxdim.f' 
      include 'limits.f'
      
      real(dp):: r(mxdim),p(mxpart,4),wt
      real(dp):: wt0,p12(4),p34(4),p3(4),p4(4),p5(4),p6(4),p7(4)
      real(dp):: p567(4),p1(4),p2(4),p67(4)
      real(dp):: lntaum,tau,xx(2),xmin,xmax,xjac,lnxmin
      integer:: i,nu 
      real(dp):: Qsq,taumin
      parameter(wt0=one/twopi**2) 
      real(dp):: wt34567,wt567,wt67,wt34
      real(dp):: sqrts
      include 'energy.f'

      wt=wt0

!----- initialize p 
      do i=1,mxpart 
         do nu=1,4 
            p(i,nu)=0._dp 
         enddo
      enddo

      do nu=1,4 
         p12(nu)=0._dp 
         p1(nu)=0._dp 
         p2(nu)=0._dp 
         p3(nu)=0._dp 
         p4(nu)=0._dp 
         p5(nu)=0._dp 
         p6(nu)=0._dp 
         p7(nu)=0._dp 
         p34(nu)=0._dp 
   
      enddo

!-------- generate Qsq. 

      if(wsqmin>1.e-4_dp) then 
         taumin=(wsqmin/sqrts**2)
      else
         taumin=(1.e-4_dp/sqrts**2)
      endif

      lntaum=log(taumin)
      tau=exp(lntaum*(one-r(1)))         
      xjac=-lntaum*tau      
      Qsq=tau*sqrts**2
      xjac=xjac*sqrts**2
  
!------ generate x1 
      xmin=tau 
      xmax=1._dp 
      lnxmin=log(xmin/xmax)
      xx(1)=xmax*exp(lnxmin*(one-r(2)))
      xjac=xjac*(-lnxmin*xx(1))/(xx(1)*sqrts**2)
      xx(2)=Qsq/(xx(1)*sqrts**2) 
  
       
      if((xx(1)>0.999_dp).or.(xx(2)>0.999_dp)) return 1
         
      
      wt=wt

!---------- generate intial state 

      p1(3)=-sqrts/2._dp*xx(1)
      p1(4)=-sqrts/2._dp*xx(1)

      p2(3)=sqrts/2._dp*xx(2)
      p2(4)=-sqrts/2._dp*xx(2)

      do nu=1,4 
         p12(nu)=-p1(nu)-p2(nu) 
      enddo 

!----- decay to Z + (gamma+jet+jet)  
      
      call phi1_2(r(3),r(4),r(5),r(6),p12,p567,p34,wt34567,*999) 

!      write(6,*) p12     
!      write(6,*) p56
!      write(6,*) p34

  
!------- decay Z=>ll
      call phi3m0(r(7),r(8),p34,p3,p4,wt34,*999) 
!-----decay gamma+(jet+jet)  
      call phi1_2m(0._dp,r(9),r(10),r(11),1.e-4_dp,p567,p5,p67,wt567,*999)
      call phi3m0(r(12),r(13),p67,p6,p7,wt67,*999) 

      wt=wt*wt34567*wt34*wt567*wt67*xjac

!----- bulid p 
      do nu=1,4 
         p(1,nu)=p1(nu)
         p(2,nu)=p2(nu)
         p(3,nu)=p3(nu)
         p(4,nu)=p4(nu)
         p(5,nu)=p5(nu)
         p(6,nu)=p6(nu)
         p(7,nu)=p7(nu)
      enddo
      

      return 
 999  continue 
      return 1 
      end 
      

