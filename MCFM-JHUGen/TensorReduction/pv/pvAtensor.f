      subroutine pvAtensor(m1s,FA0,FA1,FA2,FA3,FA4,FA5,FA6)
      implicit none
      include 'TRconstants.f'
      include 'pvAnames.f'
      include 'pvAv.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      double complex FA0(-2:0),FA1(y1max,-2:0),
     . FA2(y2max,-2:0),FA3(y3max,-2:0),FA4(y4max,-2:0),FA5(y5max,-2:0),
     . FA6(y6max,-2:0)
      double precision m1s
      double precision pvSDDDD,pvSDDDDDD
      integer n1,n2,n3,n4,n5,n6,ep,A0i,pvAcache
      logical,save:: first=.true.
!$omp threadprivate(first)
      if (first) then
      first=.false.
      call pvarraysetup
      endif


      A0i=pvAcache(m1s)
      do ep=-2,0
      FA0(ep)=Av(A0i+aa0,ep)
      enddo

      do ep=-2,0
      do n1=1,4
      FA1(n1,ep)=czip
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4 
      FA2(y2(n1,n2),ep)=g(n1,n2)*Av(A0i+aa00,ep)
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4 
      do n3=n2,4 
      FA3(y3(n1,n2,n3),ep)=czip
      enddo
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      FA4(y4(n1,n2,n3,n4),ep)=pvSDDDD(n1,n2,n3,n4)*Av(A0i+aa0000,ep)
      enddo
      enddo
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      FA5(y5(n1,n2,n3,n4,n5),ep)=czip
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      FA6(y6(n1,n2,n3,n4,n5,n6),ep)=
     . pvSDDDDDD(n1,n2,n3,n4,n5,n6)*Av(A0i+aa000000,ep)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      return
      end

