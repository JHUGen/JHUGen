      subroutine pvBtensor(q1,m1s,m2s,FB0,FB1,FB2,FB3,FB4,FB5,FB6)
      implicit none
C     q1 is the momentum in the loop = p1 the external momenta
C     m1s,m2s are the squares of the internal masses
      include 'pvBnames.f'
      include 'pvBv.f'
      include 'TRydef.f'
      double complex FB0(-2:0),FB1(y1max,-2:0),
     . FB2(y2max,-2:0),FB3(y3max,-2:0),FB4(y4max,-2:0),FB5(y5max,-2:0)
     . ,FB6(y6max,-2:0)
      double precision q1(4),q1Dq1,m1s,m2s
      double precision pvSDDP,pvSDDDD,pvSDDPP,pvSDDDDP,pvSDDPPP,
     . pvSDDDDPP,pvSDDPPPP,pvSDDDDDD
      integer n1,n2,n3,n4,n5,n6,ep,B0i,pvBcache
      include 'TRmetric.f'
      logical,save:: first=.true.
!$omp threadprivate(first)

      if (first) then
      first=.false.
      call pvarraysetup
      endif


      q1Dq1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      B0i=pvBcache(q1Dq1,m1s,m2s)
      do ep=-2,0
      FB0(ep)=Bv(B0i+bb0,ep)
      enddo

      do ep=-2,0
      do n1=1,4
      FB1(n1,ep)=Bv(B0i+bb1,ep)*q1(n1)
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4 
      FB2(y2(n1,n2),ep)=q1(n1)*q1(n2)*Bv(B0i+bb11,ep)
     . +g(n1,n2)*Bv(B0i+bb00,ep)
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4 
      do n3=n2,4 
      FB3(y3(n1,n2,n3),ep)=q1(n1)*q1(n2)*q1(n3)*Bv(B0i+bb111,ep)
     . +pvSDDP(n1,n2,n3,q1)*Bv(B0i+bb001,ep)
      enddo
      enddo
      enddo
      enddo


      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      FB4(y4(n1,n2,n3,n4),ep)=
     . +q1(n1)*q1(n2)*q1(n3)*q1(n4)*Bv(B0i+bb1111,ep)
     . +pvSDDPP(n1,n2,n3,n4,q1)*Bv(B0i+bb0011,ep)
     . +pvSDDDD(n1,n2,n3,n4)*Bv(B0i+bb0000,ep)


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

      FB5(y5(n1,n2,n3,n4,n5),ep)=
     . +Bv(B0i+bb00001,ep)*pvSDDDDP(n1,n2,n3,n4,n5,q1)
     . +Bv(B0i+bb00111,ep)*pvSDDPPP(n1,n2,n3,n4,n5,q1)
     . +Bv(B0i+bb11111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)


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

      FB6(y6(n1,n2,n3,n4,n5,n6),ep)=
     . +Bv(B0i+bb000000,ep)*pvSDDDDDD(n1,n2,n3,n4,n5,n6)
     . +Bv(B0i+bb000011,ep)*pvSDDDDPP(n1,n2,n3,n4,n5,n6,q1)
     . +Bv(B0i+bb001111,ep)*pvSDDPPPP(n1,n2,n3,n4,n5,n6,q1)
     . +Bv(B0i+bb111111,ep)*q1(n1)*q1(n2)*q1(n3)*q1(n4)*q1(n5)*q1(n6)


      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


      return
      end

