      subroutine computescalars(k1,k2,k3,k4,k5,k6,scints)
      implicit none
      include 'types.f'
c--- routine to compute all scalar integrals used in the Wbbm calculation
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'momwbbm.f'
      include 'scale.f'
      include 'Wbbmlabels.f'
      include 'first.f'
      logical:: writescalars
      integer:: k1,k2,k3,k4,k5,k6,nu,iep
      real(dp):: p2(4),p3(4),p23(4),p123(4),p234(4),p1234(4),
     & p12(4),p34(4),s23,s123,s234,s34,s12,s1234,msq
      complex(dp):: qlI4,qlI3,qlI2,qlI1
      common/writescalars/writescalars
!$omp threadprivate(/writescalars/)

c--- initialize QCDLoop, if necessary
      if (first) then
        call qlinit
        first=.false.
        if (writescalars) then
          open(unit=67,file='scalars.out',status='unknown')
        endif
      endif

      
      do nu=1,4
      p2(nu)=bp*mom(k2,nu)+bm*mom(k3,nu)
      p3(nu)=bp*mom(k3,nu)+bm*mom(k2,nu)
      p12(nu)=mom(k1,nu)+p2(nu)
      p34(nu)=p3(nu)+mom(k4,nu)
      p23(nu)=p2(nu)+p3(nu)
      p123(nu)=p23(nu)+mom(k1,nu)
      p234(nu)=p23(nu)+mom(k4,nu)
      p1234(nu)=p123(nu)+mom(k4,nu)
      enddo
      
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      s123=p123(4)**2-p123(1)**2-p123(2)**2-p123(3)**2
      s234=p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2
      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s1234=p1234(4)**2-p1234(1)**2-p1234(2)**2-p1234(3)**2
      msq=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
c      mb=sqrt(msq)
 
c      write(6,*)    
c      write(6,*) 'Kinematics:'    
c      write(6,*) 's123=',s123
c      write(6,*) 's234=',s234
c      write(6,*) 's1234=',s1234
c      write(6,*) 's23=',s23
c      write(6,*) 's12=',s12
c      write(6,*) 's13=',msq+two*(mom(k1,4)*p3(4)
c     & -mom(k1,1)*p3(1)-mom(k1,2)*p3(2)-mom(k1,3)*p3(3))
c      write(6,*) 's34=',s34
c      write(6,*) 's24=',msq+two*(mom(k4,4)*p2(4)
c     & -mom(k4,1)*p2(1)-mom(k4,2)*p2(2)-mom(k4,3)*p2(3))
c      write(6,*) 'bm= ',bm
c      write(6,*) 'bp= ',bp
c      write(6,*)          

c---  NB: only need finite pieces, poles will be handled separately
c      do iep=-2,0
      do iep=0,0

c--- box integrals
      scints(4,d2x3x4,iep)=qlI4(msq,msq,zip,s234,s23,s34,
     &                          zip,msq,zip,zip,musq,iep)
      scints(4,d1x2x3,iep)=qlI4(zip,msq,msq,s123,s12,s23,
     &                          zip,zip,msq,zip,musq,iep)
      scints(4,d1x23x4,iep)=qlI4(zip,s23,zip,s1234,s123,s234,
     &                          zip,zip,zip,zip,musq,iep)
      scints(4,d1x2x34,iep)=qlI4(zip,msq,s34,s1234,s12,s234,
     &                          zip,zip,msq,zip,musq,iep)
      scints(4,d12x3x4,iep)=qlI4(s12,msq,zip,s1234,s123,s34,
     &                          zip,msq,zip,zip,musq,iep)

c--- triangle integrals
      scints(3,c23x4,iep)=qlI3(s23,zip,s234,zip,zip,zip,musq,iep)
      scints(3,c1x23,iep)=qlI3(zip,s23,s123,zip,zip,zip,musq,iep)
      scints(3,c2x3,iep)=qlI3(msq,msq,s23,zip,msq,zip,musq,iep)
      scints(3,c2x34,iep)=qlI3(msq,s34,s234,zip,msq,zip,musq,iep)
      scints(3,c3x4,iep)=qlI3(msq,zip,s34,msq,zip,zip,musq,iep)
      scints(3,c12x3,iep)=qlI3(s12,msq,s123,zip,msq,zip,musq,iep)
      scints(3,c1x2,iep)=qlI3(zip,msq,s12,zip,zip,msq,musq,iep)
      scints(3,c1x234,iep)=qlI3(zip,s234,s1234,zip,zip,zip,musq,iep)
      scints(3,c123x4,iep)=qlI3(s123,zip,s1234,zip,zip,zip,musq,iep)
      scints(3,c12x34,iep)=qlI3(s12,s34,s1234,zip,msq,zip,musq,iep)

c--- bubble integrals
      scints(2,b23,iep)=qlI2(s23,zip,zip,musq,iep)
      scints(2,b234,iep)=qlI2(s234,zip,zip,musq,iep)
      scints(2,b123,iep)=qlI2(s123,zip,zip,musq,iep)
      scints(2,b2x1m,iep)=qlI2(msq,zip,msq,musq,iep)
      scints(2,b34,iep)=qlI2(s34,zip,msq,musq,iep)
      scints(2,b12,iep)=qlI2(s12,zip,msq,musq,iep)
      scints(2,b1234,iep)=qlI2(s1234,zip,zip,musq,iep)
 
c--- tadpole integrals
      scints(1,a0m,iep)=qlI1(msq,musq,iep)
      enddo
 
      if (writescalars) then
      write(67,*) 'scints(4,d2x3x4,0) ',scints(4,d2x3x4,0)
      write(67,*) 'scints(4,d1x2x3,0) ',scints(4,d1x2x3,0)
      write(67,*) 'scints(4,d1x23x4,0)',scints(4,d1x23x4,0)
      write(67,*) 'scints(4,d1x2x34,0)',scints(4,d1x2x34,0)
      write(67,*) 'scints(4,d12x3x4,0)',scints(4,d12x3x4,0)
      write(67,*) 'done boxes'
      write(67,*) 'scints(3,c23x4,0)  ',scints(3,c23x4,0)
      write(67,*) 'scints(3,c1x23,0)  ',scints(3,c1x23,0)
      write(67,*) 'scints(3,c2x3,0)   ',scints(3,c2x3,0)
      write(67,*) 'scints(3,c2x34,0)  ',scints(3,c2x34,0)
      write(67,*) 'scints(3,c3x4,0)   ',scints(3,c3x4,0)
      write(67,*) 'scints(3,c12x3,0)  ',scints(3,c12x3,0)
      write(67,*) 'scints(3,c1x2,0)   ',scints(3,c1x2,0)
      write(67,*) 'scints(3,c1x234,0) ',scints(3,c1x234,0)
      write(67,*) 'scints(3,c123x4,0) ',scints(3,c123x4,0)
      write(67,*) 'scints(3,c12x34,0) ',scints(3,c12x34,0)
      write(67,*) 'done triangles'
      write(67,*) 'scints(2,b23,0)    ',scints(2,b23,0)   
      write(67,*) 'scints(2,b234,0)   ',scints(2,b234,0)
      write(67,*) 'scints(2,b123,0)   ',scints(2,b123,0)
      write(67,*) 'scints(2,b2x1m,0)  ',scints(2,b2x1m,0)
      write(67,*) 'scints(2,b34,0)    ',scints(2,b34,0)
      write(67,*) 'scints(2,b12,0)    ',scints(2,b12,0)
      write(67,*) 'scints(2,b1234,0)  ',scints(2,b1234,0)
      write(67,*) 'done bubbles'
      write(67,*) 'scints(1,a0m,0)    ',scints(1,a0m,0)
      write(67,*) 'done tadpole'
      write(67,*)
      endif

      return
      end
      
