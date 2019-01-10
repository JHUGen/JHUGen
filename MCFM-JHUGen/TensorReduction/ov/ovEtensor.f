      subroutine ovEtensor(p1,p2,p3,p4,m1s,m2s,m3s,m4s,m5s,
     &  FE0,FE1,FE2,FE3,FE4,FE5)
C****NB The arguments of this routine are the external momenta "p",
C       and not the momentum offsets in the loop "q" ******
C****   m1s,m2s,m3s,m4s,m5s are the internal masses squared.

      implicit none
      include 'TRydef.f'
      include 'TRbadpoint.f'
      include 'TRmaxindex.f'
      include 'TRconstants.f'
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),
     & inGram(4,4),Gram(4,4),
     & s1Dp1,s1Dp2,s1Dp3,s1Dp4,s2Dp2,s2Dp3,s2Dp4,s3Dp3,s3Dp4,s4Dp4
      double complex FE0(-2:0),FE1(y1max,-2:0),FE2(y2max,-2:0),
     & FE3(y3max,-2:0),FE4(y4max,-2:0),FE5(y5max,-2:0),
     & FD0_1(-2:0),FD1_1(y1max,-2:0),FD2_1(y2max,-2:0),
     & FD3_1(y3max,-2:0),FD4_1(y4max,-2:0),
     & FD0_2(-2:0),FD1_2(y1max,-2:0),FD2_2(y2max,-2:0),
     & FD3_2(y3max,-2:0),FD4_2(y4max,-2:0),
     & FD0_3(-2:0),FD1_3(y1max,-2:0),FD2_3(y2max,-2:0),
     & FD3_3(y3max,-2:0),FD4_3(y4max,-2:0),
     & FD0_4(-2:0),FD1_4(y1max,-2:0),FD2_4(y2max,-2:0),
     & FD3_4(y3max,-2:0),FD4_4(y4max,-2:0),
     & FD0_5(-2:0),FD1_5(y1max,-2:0),FD2_5(y2max,-2:0),
     & FD3_5(y3max,-2:0),FD4_5(y4max,-2:0),
     & FD1_1a(y1max,-2:0),FD2_1a(y2max,-2:0),
     & FD3_1a(y3max,-2:0),FD4_1a(y4max,-2:0),RHS(4),inRHS(4)
      double precision p1Dp1,p2Dp2,p3Dp3,p4Dp4,p5Dp5,d,
     & p1Dp2,p1Dp3,p2Dp3,p1Dp4,p2Dp4,p3Dp4,
     & m1s,m2s,m3s,m4s,m5s,s12,s23,s34,s45,s51,
     & p12(4),p23(4),p34(4),p45(4),p51(4),p123(4),p1234(4)
      integer nu,n1,n2,n3,n4,n5,ep,indx(4)
      logical failed
      logical,save:: first=.true.
!$omp threadprivate(first)
      if (first) then
      first=.false.
      call ovarraysetup
      endif


      do nu=1,4
      p5(nu)=-p1(nu)-p2(nu)-p3(nu)-p4(nu)
      p12(nu)=p1(nu)+p2(nu)
      p23(nu)=p2(nu)+p3(nu)
      p34(nu)=p3(nu)+p4(nu)
      p45(nu)=p4(nu)+p5(nu)
      p51(nu)=p5(nu)+p1(nu)
      enddo


      p1Dp1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2
      p2Dp2=p2(4)**2-p2(1)**2-p2(2)**2-p2(3)**2
      p3Dp3=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4Dp4=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      p5Dp5=p5(4)**2-p5(1)**2-p5(2)**2-p5(3)**2

      p3Dp4=p3(4)*p4(4)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      p2Dp4=p2(4)*p4(4)-p2(1)*p4(1)-p2(2)*p4(2)-p2(3)*p4(3)
      p1Dp4=p1(4)*p4(4)-p1(1)*p4(1)-p1(2)*p4(2)-p1(3)*p4(3)

      p2Dp3=p2(4)*p3(4)-p2(1)*p3(1)-p2(2)*p3(2)-p2(3)*p3(3)
      p1Dp3=p1(4)*p3(4)-p1(1)*p3(1)-p1(2)*p3(2)-p1(3)*p3(3)

      p1Dp2=p1(4)*p2(4)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)

      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      s23=p23(4)**2-p23(1)**2-p23(2)**2-p23(3)**2
      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s45=p45(4)**2-p45(1)**2-p45(2)**2-p45(3)**2
      s51=p51(4)**2-p51(1)**2-p51(2)**2-p51(3)**2

      s1Dp1=0.5d0*(m2s-m1s-p1Dp1)
      s2Dp2=0.5d0*(m3s-m2s-p2Dp2)
      s3Dp3=0.5d0*(m4s-m3s-p3Dp3)
      s4Dp4=0.5d0*(m5s-m4s-p4Dp4)

      s3Dp4=s4Dp4-p3Dp4
      s2Dp4=s3Dp4-p2Dp4
      s1Dp4=s2Dp4-p1Dp4

      s2Dp3=s3Dp3-p2Dp3
      s1Dp3=s2Dp3-p1Dp3

      s1Dp2=s2Dp2-p1Dp2


      call ovE0scalar(FE0,p1Dp1,p2Dp2,p3Dp3,p4Dp4,p5Dp5,
     & s12,s23,s34,s45,s51,m1s,m2s,m3s,m4s,m5s)

      inGram(1,1)=p1Dp1
      inGram(2,2)=p2Dp2
      inGram(3,3)=p3Dp3
      inGram(4,4)=p4Dp4
      inGram(1,2)=p1Dp2
      inGram(1,3)=p1Dp3
      inGram(1,4)=p1Dp4
      inGram(2,3)=p2Dp3
      inGram(2,4)=p2Dp4
      inGram(3,4)=p3Dp4
      inGram(2,1)=inGram(1,2)
      inGram(3,1)=inGram(1,3)
      inGram(4,1)=inGram(1,4)
      inGram(3,2)=inGram(2,3)
      inGram(4,2)=inGram(2,4)
      inGram(4,3)=inGram(3,4)

      call ludcmp(inGram,4,indx,d,Gram)

      call ovDtensor(p2,p3,p4,m2s,m3s,m4s,m5s,
     & FD0_1,FD1_1,FD2_1,FD3_1,FD4_1)
      call ovDtensor(p12,p3,p4,m1s,m3s,m4s,m5s,
     & FD0_2,FD1_2,FD2_2,FD3_2,FD4_2)
      call ovDtensor(p1,p23,p4,m1s,m2s,m4s,m5s,
     & FD0_3,FD1_3,FD2_3,FD3_3,FD4_3)
      call ovDtensor(p1,p2,p34,m1s,m2s,m3s,m5s,
     & FD0_4,FD1_4,FD2_4,FD3_4,FD4_4)
      call ovDtensor(p1,p2,p3,m1s,m2s,m3s,m4s,
     & FD0_5,FD1_5,FD2_5,FD3_5,FD4_5)

      call ovswitch(p1,FD0_1,FD1_1,FD2_1,FD3_1,FD4_1,
     & FD1_1a,FD2_1a,FD3_1a,FD4_1a)


      do ep=-2,0
      do n1=1,4
      inRHS(1)=s1Dp1*FE0(ep)+half*(FD0_2(ep)-FD0_1(ep))
      inRHS(2)=s1Dp2*FE0(ep)+half*(FD0_3(ep)-FD0_2(ep))
      inRHS(3)=s1Dp3*FE0(ep)+half*(FD0_4(ep)-FD0_3(ep))
      inRHS(4)=s1Dp4*FE0(ep)+half*(FD0_5(ep)-FD0_4(ep))
      call zlubksb(Gram,4,indx,inRHS,RHS)
      FE1(n1,ep)=RHS(1)*p1(n1)+RHS(2)*p2(n1)+RHS(3)*p3(n1)+RHS(4)*p4(n1)
      enddo
      enddo
      if (maxeindex .eq. 1) goto 99


      do ep=-2,0
      do n1=1,4
      inRHS(1)=s1Dp1*FE1(n1,ep)+half*(FD1_2(n1,ep)-FD1_1a(n1,ep))
      inRHS(2)=s1Dp2*FE1(n1,ep)+half*(FD1_3(n1,ep)-FD1_2(n1,ep))
      inRHS(3)=s1Dp3*FE1(n1,ep)+half*(FD1_4(n1,ep)-FD1_3(n1,ep))
      inRHS(4)=s1Dp4*FE1(n1,ep)+half*(FD1_5(n1,ep)-FD1_4(n1,ep))
      call zlubksb(Gram,4,indx,inRHS,RHS)
      do n2=n1,4
      FE2(y2(n1,n2),ep)=
     & RHS(1)*p1(n2)+RHS(2)*p2(n2)+RHS(3)*p3(n2)+RHS(4)*p4(n2)
      enddo
      enddo
      enddo
      if (maxeindex .eq. 2) goto 99

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      inRHS(1)=s1Dp1*FE2(y2(n1,n2),ep)
     & +half*(FD2_2(y2(n1,n2),ep)-FD2_1a(y2(n1,n2),ep))
      inRHS(2)=s1Dp2*FE2(y2(n1,n2),ep)
     & +half*(FD2_3(y2(n1,n2),ep)-FD2_2(y2(n1,n2),ep))
      inRHS(3)=s1Dp3*FE2(y2(n1,n2),ep)
     & +half*(FD2_4(y2(n1,n2),ep)-FD2_3(y2(n1,n2),ep))
      inRHS(4)=s1Dp4*FE2(y2(n1,n2),ep)
     & +half*(FD2_5(y2(n1,n2),ep)-FD2_4(y2(n1,n2),ep))
      call zlubksb(Gram,4,indx,inRHS,RHS)
      do n3=n2,4
      FE3(y3(n1,n2,n3),ep)=
     & RHS(1)*p1(n3)+RHS(2)*p2(n3)+RHS(3)*p3(n3)+RHS(4)*p4(n3)
      enddo
      enddo
      enddo
      enddo
      if (maxeindex .eq. 3) goto 99

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      inRHS(1)=s1Dp1*FE3(y3(n1,n2,n3),ep)
     & +half*(FD3_2(y3(n1,n2,n3),ep)-FD3_1a(y3(n1,n2,n3),ep))
      inRHS(2)=s1Dp2*FE3(y3(n1,n2,n3),ep)
     & +half*(FD3_3(y3(n1,n2,n3),ep)-FD3_2(y3(n1,n2,n3),ep))
      inRHS(3)=s1Dp3*FE3(y3(n1,n2,n3),ep)
     & +half*(FD3_4(y3(n1,n2,n3),ep)-FD3_3(y3(n1,n2,n3),ep))
      inRHS(4)=s1Dp4*FE3(y3(n1,n2,n3),ep)
     &+half*(FD3_5(y3(n1,n2,n3),ep)-FD3_4(y3(n1,n2,n3),ep))
      call zlubksb(Gram,4,indx,inRHS,RHS)
      do n4=n3,4
      FE4(y4(n1,n2,n3,n4),ep)=
     & RHS(1)*p1(n4)+RHS(2)*p2(n4)+RHS(3)*p3(n4)+RHS(4)*p4(n4)
      enddo
      enddo
      enddo
      enddo
      enddo
      if (maxeindex .eq. 4) goto 99

      do ep=-2,0
      do n1=1,4
      do n2=n1,4
      do n3=n2,4
      do n4=n3,4
      inRHS(1)=s1Dp1*FE4(y4(n1,n2,n3,n4),ep)
     & +half*(FD4_2(y4(n1,n2,n3,n4),ep)-FD4_1a(y4(n1,n2,n3,n4),ep))
      inRHS(2)=s1Dp2*FE4(y4(n1,n2,n3,n4),ep)
     & +half*(FD4_3(y4(n1,n2,n3,n4),ep)-FD4_2(y4(n1,n2,n3,n4),ep))
      inRHS(3)=s1Dp3*FE4(y4(n1,n2,n3,n4),ep)
     & +half*(FD4_4(y4(n1,n2,n3,n4),ep)-FD4_3(y4(n1,n2,n3,n4),ep))
      inRHS(4)=s1Dp4*FE4(y4(n1,n2,n3,n4),ep)
     & +half*(FD4_5(y4(n1,n2,n3,n4),ep)-FD4_4(y4(n1,n2,n3,n4),ep))
      call zlubksb(Gram,4,indx,inRHS,RHS)
      do n5=n4,4
      FE5(y5(n1,n2,n3,n4,n5),ep)=
     & RHS(1)*p1(n5)+RHS(2)*p2(n5)+RHS(3)*p3(n5)+RHS(4)*p4(n5)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


c--- before returning, check tensor computed correctly        
   99 continue
      
      p123=p12+p3
      p1234=p123+p4
      call ovEcheck(maxeindex,p1,p12,p123,p1234,m1s,m2s,m3s,m4s,m5s,
     &  FE0,FE1,FE2,FE3,FE4,FE5,failed)
      if (failed) then
        write(6,*) 'badpoint set in ovEtensor'
        pvbadpoint=.true.
      endif


      return
      end

