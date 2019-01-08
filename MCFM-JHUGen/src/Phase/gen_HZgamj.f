      subroutine gen_HZgamj(r,p,wt,*) 
      implicit none 
      include 'constants.f' 
      include 'mxdim.f' 
      include 'limits.f'
      include 'breit.f'
      include 'masses.f'
      include 'energy.f'
      
      double precision r(mxdim),p(mxpart,4),wt
      double precision wt0,p12(4),p34(4),p3(4),p4(4),p5(4),p6(4)
      double precision p345(4),p1(4),p2(4)
      double precision lntaum,tau,xx(2),xmin,xmax,xjac,lnxmin,cutoff
      integer i,nu 
      double precision Qsq,taumin
      parameter(wt0=one/twopi**2) 
      double precision wt3456,wt345,wt34

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
         p6(nu)=0d0 
         p34(nu)=0d0    
      enddo

!-------- generate Qsq. 
      cutoff=max(1d-4,(hmass-25d0*hwidth)**2)
      taumin=(cutoff/sqrts**2)

      lntaum=dlog(taumin)
      tau=dexp(lntaum*(one-r(1)))         
      xjac=-lntaum*tau      
      Qsq=tau*sqrts**2
      xjac=xjac*sqrts**2
  
!------ generate x1 
      xmin=tau 
      xmax=1d0 
      lnxmin=dlog(xmin/xmax)
      xx(1)=xmax*dexp(lnxmin*(one-r(2)))
      xjac=xjac*(-lnxmin*xx(1))/(xx(1)*sqrts**2)
      xx(2)=Qsq/(xx(1)*sqrts**2) 
                
!---------- generate intial state 
      p1(3)=-sqrts/2d0*xx(1)
      p1(4)=-sqrts/2d0*xx(1)

      p2(3)=sqrts/2d0*xx(2)
      p2(4)=-sqrts/2d0*xx(2)

      do nu=1,4 
         p12(nu)=-p1(nu)-p2(nu) 
      enddo 

!----- decay to H + jet        
      n2=0
      n3=1
      mass3=hmass
      width3=hwidth
      call phi1_2m(0d0,r(3),r(4),r(5),cutoff,p12,p6,p345,wt3456,*999)

!----- decay to Z+photon
      n2=0
      n3=1
      mass3=zmass
      width3=zwidth
      call phi1_2m(0d0,r(6),r(7),r(8),wsqmin,p345,p5,p34,wt345,*999)
  
!------- decay Z=>ll
      call phi3m0(r(9),r(10),p34,p3,p4,wt34,*999) 

      wt=wt0*xjac*wt3456*wt345*wt34

!----- build p 
      do nu=1,4 
         p(1,nu)=p1(nu)
         p(2,nu)=p2(nu)
         p(3,nu)=p3(nu)
         p(4,nu)=p4(nu)
         p(5,nu)=p5(nu)
         p(6,nu)=p6(nu)
      enddo
      
      return 
 999  continue 
      return 1 
      end 
      

