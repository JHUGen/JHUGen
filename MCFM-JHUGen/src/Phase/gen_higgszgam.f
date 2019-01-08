      subroutine gen_higgszgam(r,p,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'mxdim.f' 
      include 'limits.f'
      include 'masses.f' 
      integer i,nu 
      double precision r(mxdim),p(mxpart,4),wt
      double precision wt0,p12(4),p34(4),p3(4),p4(4),p5(4)
      double precision  p1(4),p2(4)
      double precision tau,xx(2),xmin,xmax,xjac,lnxmin
      double precision Qsq,wt345,wt34
      include 'energy.f'
      parameter(wt0=one/twopi) 

      wt=wt0

!----- initialize p 
      do i=1,mxpart 
         do nu=1,4 
            p(i,nu)=0d0 
         enddo
      enddo

      do nu=1,4 
         p12(nu)=0d0 
         p1(nu)=0d0 
         p2(nu)=0d0 
         p3(nu)=0d0 
         p4(nu)=0d0 
         p5(nu)=0d0 
      enddo

!----- generate Qsq around Higgs BW
      call breitw(r(1),wsqmin,sqrts**2,hmass,hwidth,Qsq,xjac)
        
      tau=Qsq/sqrts**2
!------ generate x1 
      xmin=tau 
      xmax=1d0 
      lnxmin=dlog(xmin/xmax)
      xx(1)=xmax*dexp(lnxmin*(one-r(2)))
      xjac=xjac*(-lnxmin*xx(1))/(xx(1)*sqrts**2)
      xx(2)=Qsq/(xx(1)*sqrts**2) 
  
c--- check that xx(1) and xx(2) are in range
      if((xx(1).gt.1d0).or.(xx(2).gt.1d0)) return 1         
   
!---------- generate intial state 

      p1(3)=-sqrts/2d0*xx(1)
      p1(4)=-sqrts/2d0*xx(1)

      p2(3)=sqrts/2d0*xx(2)
      p2(4)=-sqrts/2d0*xx(2)

      do nu=1,4 
         p12(nu)=-p1(nu)-p2(nu) 
      enddo 

!----- decay to Z + (gamma)  
      
      call phi1_2m(0d0,r(3),r(4),r(5),wsqmin,p12,p5,p34,wt345,*999) 
  
!------- decay Z=>ll
      call phi3m0(r(6),r(7),p34,p3,p4,wt34,*999) 

      wt=wt*wt345*wt34*xjac

!----- build p 
      do nu=1,4 
         p(1,nu)=p1(nu)
         p(2,nu)=p2(nu)
         p(3,nu)=p3(nu)
         p(4,nu)=p4(nu)
         p(5,nu)=p5(nu)
      enddo
      
      
      return 
 999  continue 
      return 1 
      end 
      

