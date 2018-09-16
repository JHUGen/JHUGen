      subroutine gen_HZgamj(r,p,wt,*) 
      implicit none
      include 'types.f'
       
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'mxdim.f' 
      include 'limits.f'
      include 'breit.f'
      include 'masses.f'
      include 'energy.f'
      
      real(dp):: r(mxdim),p(mxpart,4),wt
      real(dp):: wt0,p12(4),p34(4),p3(4),p4(4),p5(4),p6(4)
      real(dp):: p345(4),p1(4),p2(4)
      real(dp):: lntaum,tau,xx(2),xmin,xmax,xjac,lnxmin,cutoff
      integer:: i,nu 
      real(dp):: Qsq,taumin
      parameter(wt0=one/twopi**2) 
      real(dp):: wt3456,wt345,wt34

!----- initialize p 
      p(:,:)=zip 

      p12(:)=zip 
      p1(:)=zip 
      p2(:)=zip 
      p3(:)=zip 
      p4(:)=zip 
      p5(:)=zip 
      p6(:)=zip 
      p34(:)=zip    

!-------- generate Qsq. 
      cutoff=max(1.e-4_dp,(hmass-25._dp*hwidth)**2)
      taumin=(cutoff/sqrts**2)

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
                
!---------- generate intial state 
      p1(3)=-sqrts/2._dp*xx(1)
      p1(4)=-sqrts/2._dp*xx(1)

      p2(3)=sqrts/2._dp*xx(2)
      p2(4)=-sqrts/2._dp*xx(2)

      do nu=1,4 
         p12(nu)=-p1(nu)-p2(nu) 
      enddo 

!----- decay to H + jet        
      n2=0
      n3=1
      mass3=hmass
      width3=hwidth
      call phi1_2m(zip,r(3),r(4),r(5),cutoff,p12,p6,p345,wt3456,*999)

!----- decay to Z+photon
      n2=0
      n3=1
      mass3=zmass
      width3=zwidth
      call phi1_2m(zip,r(6),r(7),r(8),wsqmin,p345,p5,p34,wt345,*999)
  
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
      

