      subroutine gen4handc(r,p,wt4,*)
c--- Generates 2->4 phase space with (3456) a Breit-Wigner around
c--- the Higgs mass;   width of the BW is artificially enhanced to
c--- better map the continuum contribution
      implicit none
      include 'constants.f'
      include 'mxdim.f'
      include 'debug.f'
      include 'masses.f'
      include 'phasemin.f'
      include 'breit.f'
      include 'interference.f'
      include 'x1x2.f'
      include 'first.f'
      integer nu
C      integer,save::icount=1
C!$omp threadprivate(icount)
      double precision r(mxdim)
      double precision wt4,p1(4),p2(4),p3(4),p4(4),p5(4),p6(4)
      double precision p(mxpart,4),rtshat
      double precision pswt,xjac,pswidth
      double precision s3456,wt3456,ymax,yave,lntaum,tau
      include 'energy.f'
      save pswidth
!$omp threadprivate(pswidth)

      wt4=0d0
      
      if (first) then
         first=.false.
c--- width to use in generation of PS: if too far from threshold, just
c--- use a width of 10 GeV in the B.W. to sample PS adequately
        if (hmass .lt. mass2+mass3-hwidth*5d0) then
          pswidth=10d0
        else
          pswidth=hwidth
        endif
      endif
      
      call breitw(r(9),0d0,sqrts**2,hmass,pswidth,s3456,wt3456)
            
      rtshat=dsqrt(s3456)
      ymax=dlog(sqrts/rtshat)
      yave=ymax*(two*r(10)-1d0)
      xjac=two*ymax*wt3456
           
      xx(1)=rtshat/sqrts*exp(+yave)
      xx(2)=rtshat/sqrts*exp(-yave)

      if   ((xx(1) .gt. 1d0) 
     & .or. (xx(2) .gt. 1d0)
     & .or. (xx(1) .lt. xmin)
     & .or. (xx(2) .lt. xmin)) then
c      write(6,*) 'problems with xx(1),xx(2) in gen4h',xx(1),xx(2)  
      return 1 
      endif

      p1(4)=-0.5d0*xx(1)*sqrts
      p1(1)=0d0
      p1(2)=0d0
      p1(3)=-0.5d0*xx(1)*sqrts
      
      p2(4)=-0.5d0*xx(2)*sqrts
      p2(1)=0d0
      p2(2)=0d0
      p2(3)=+0.5d0*xx(2)*sqrts

      call phase4(r,p1,p2,p3,p4,p5,p6,pswt,*999) 

      do nu=1,4
      p(1,nu)=p1(nu)
      p(2,nu)=p2(nu)
      p(3,nu)=p3(nu)
      p(4,nu)=p4(nu)
      p(5,nu)=p5(nu)
      p(6,nu)=p6(nu)
      p(7,nu)=0d0
      enddo 

      if (interference) then
!        if (icount .eq. 1) then
          bw34_56=.true.
!          icount=icount-1
!        else
!          bw34_56=.false.
!          do nu=1,4
!            p(4,nu)=p6(nu)
!            p(6,nu)=p4(nu)
!          enddo 
!          icount=icount+1
!        endif
      endif
      
      wt4=xjac*pswt/sqrts**2
      
      if (debug) write(6,*) 'wt4 in gen4handc',wt4
      
      return

 999  return 1
      end

