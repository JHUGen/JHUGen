      subroutine pvEcheck(rank,q1,q2,q3,q4,m1s,m2s,m3s,m4s,m5s,
     & FE0,FE1,FE2,FE3,FE4,FE5,failed)
      implicit none
      include 'TRydef.f'
      include 'pvverbose.f'
      integer n2,n3,n4,ep,nu,rank,epmin
c     integer n5
      double precision q1(4),q2(4),q3(4),q4(4),p2(4),p23(4),p234(4),Eacc
      double precision q1Dq1,q2Dq2,q3Dq3,q4Dq4,m1s,m2s,m3s,m4s,m5s,
     & f1,f2,f3,f4
      double complex FE0(-2:0),FE1(y1max,-2:0),FE2(y2max,-2:0),
     & FE3(y3max,-2:0),FE4(y4max,-2:0),FE5(y5max,-2:0),
     & FD01(-2:0),FD11(y1max,-2:0),FD21(y2max,-2:0),
     & FD31(y3max,-2:0),FD41(y4max,-2:0),FD51(y5max,-2:0),
     & FD61(y6max,-2:0),
     & FD02(-2:0),FD12(y1max,-2:0),FD22(y2max,-2:0),
     & FD32(y3max,-2:0),FD42(y4max,-2:0),FD52(y5max,-2:0),
     & FD62(y6max,-2:0),
     & FD03(-2:0),FD13(y1max,-2:0),FD23(y2max,-2:0),
     & FD33(y3max,-2:0),FD43(y4max,-2:0),FD53(y5max,-2:0),
     & FD63(y6max,-2:0),
     & FD04(-2:0),FD14(y1max,-2:0),FD24(y2max,-2:0),
     & FD34(y3max,-2:0),FD44(y4max,-2:0),FD54(y5max,-2:0),
     & FD64(y6max,-2:0),
     & FD05(-2:0),FD15(y1max,-2:0),FD25(y2max,-2:0),
     & FD35(y3max,-2:0),FD45(y4max,-2:0),FD55(y5max,-2:0),
     & FD65(y6max,-2:0),
     & FD15a(y1max,-2:0),FD25a(y2max,-2:0),
     & FD35a(y3max,-2:0),FD45a(y4max,-2:0),
     & trhs,tq
      parameter(epmin=0) ! Only check finite pieces

      logical failed
      
      failed=.false.

      Eacc=1d-8
      
      do nu=1,4
      p2(nu)=q2(nu)-q1(nu)
      p23(nu)=q3(nu)-q1(nu)
      p234(nu)=q4(nu)-q1(nu)
      enddo

      q1Dq1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2Dq2=q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2
      q3Dq3=q3(4)**2-q3(1)**2-q3(2)**2-q3(3)**2
      q4Dq4=q4(4)**2-q4(1)**2-q4(2)**2-q4(3)**2

      f1=m2s-m1s-q1Dq1
      f2=m3s-m1s-q2Dq2
      f3=m4s-m1s-q3Dq3
      f4=m5s-m1s-q4Dq4

c      do ep=epmin,-1
c      sing4(ep)=zip
c      sing5(ep)=zip
c      sing6(ep)=zip
c      enddo

      call pvDtensor(q2,q3,q4,m1s,m3s,m4s,m5s,
     & FD01,FD11,FD21,FD31,FD41,FD51,FD61)
      call pvDtensor(q1,q3,q4,m1s,m2s,m4s,m5s,
     & FD02,FD12,FD22,FD32,FD42,FD52,FD62)
      call pvDtensor(q1,q2,q4,m1s,m2s,m3s,m5s,
     & FD03,FD13,FD23,FD33,FD43,FD53,FD63)
      call pvDtensor(q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD04,FD14,FD24,FD34,FD44,FD54,FD64)
      call pvDtensor(p2,p23,p234,m2s,m3s,m4s,m5s,
     & FD05,FD15,FD25,FD35,FD45,FD55,FD65)

      if ((rank .eq. 2) .or. (rank .eq. 3))
     &  call pvswitch1(q1,FD05,FD15,FD15a)
      if ((rank .eq. 3) .or. (rank .eq. 4))
     &  call pvswitch2(q1,FD05,FD15,FD25,FD25a)
      if ((rank .eq. 4) .or. (rank .eq. 5))
     &  call pvswitch3(q1,FD05,FD15,FD25,FD35,FD35a)
      if ((rank .eq. 5) .or. (rank .eq. 6))
     &  call pvswitch4(q1,FD05,FD15,FD25,FD35,FD45,FD45a)

      if (rank .eq. 1) then

      if (pvverbose) write(6,*) 'q1.FE1'
      do ep=epmin,0
      tq=q1(4)*FE1(4,ep)
     &  -q1(1)*FE1(1,ep)
     &  -q1(2)*FE1(2,ep)
     &  -q1(3)*FE1(3,ep)
      trhs=-0.5d0*(FD01(ep)-FD05(ep)+f1*FE0(ep))
      call checkaccuracy(trhs,tq,Eacc,failed)
c      if (pvverbose) write(6,*) tq
c      if (pvverbose) write(6,*) -0.5d0*(FD01(ep))
c      if (pvverbose) write(6,*) -0.5d0*(-FD04(ep))
c      if (pvverbose) write(6,*) -0.5d0*(+f1*FE0(ep))
      enddo

      if (pvverbose) write(6,*) 'q2.FE1'
      do ep=epmin,0
      tq=q2(4)*FE1(4,ep)
     &  -q2(1)*FE1(1,ep)
     &  -q2(2)*FE1(2,ep)
     &  -q2(3)*FE1(3,ep)
      trhs=-0.5d0*(FD02(ep)-FD05(ep)+f2*FE0(ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo

      if (pvverbose) write(6,*) 'q3.FE1'
      do ep=epmin,0
      tq=q3(4)*FE1(4,ep)
     &  -q3(1)*FE1(1,ep)
     &  -q3(2)*FE1(2,ep)
     &  -q3(3)*FE1(3,ep)
      trhs=-0.5d0*(FD03(ep)-FD05(ep)+f3*FE0(ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo

      if (pvverbose) write(6,*) 'q4.FE1'
      do ep=epmin,0
      tq=q4(4)*FE1(4,ep)
     &  -q4(1)*FE1(1,ep)
     &  -q4(2)*FE1(2,ep)
     &  -q4(3)*FE1(3,ep)
      trhs=-0.5d0*(FD04(ep)-FD05(ep)+f4*FE0(ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      
      
      elseif (rank .eq. 2) then

      if (pvverbose) write(6,*) 'q1.FE2'
      do ep=epmin,0
      do n2=1,4
      tq =q1(4)*FE2(y2(4,n2),ep)
     &   -q1(1)*FE2(y2(1,n2),ep)
     &   -q1(2)*FE2(y2(2,n2),ep)
     &   -q1(3)*FE2(y2(3,n2),ep)   
      trhs= 
     &   -0.5d0*(FD11(n2,ep)-FD15a(n2,ep)+f1*FE1(n2,ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FE2'
      do ep=epmin,0
      do n2=1,4
      tq =q2(4)*FE2(y2(4,n2),ep)
     &   -q2(1)*FE2(y2(1,n2),ep)
     &   -q2(2)*FE2(y2(2,n2),ep)
     &   -q2(3)*FE2(y2(3,n2),ep)
      trhs=
     &   -0.5d0*(FD12(n2,ep)-FD15a(n2,ep)+f2*FE1(n2,ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'q3.FE2'
      do ep=epmin,0
      do n2=1,4
      tq = q3(4)*FE2(y2(4,n2),ep)
     &    -q3(1)*FE2(y2(1,n2),ep)
     &    -q3(2)*FE2(y2(2,n2),ep)
     &    -q3(3)*FE2(y2(3,n2),ep)
      trhs=
     & -0.5d0*(FD13(n2,ep)-FD15a(n2,ep)+f3*FE1(n2,ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'q4.FE2'
      do ep=epmin,0
      do n2=1,4
      tq = q4(4)*FE2(y2(4,n2),ep)
     &    -q4(1)*FE2(y2(1,n2),ep)
     &    -q4(2)*FE2(y2(2,n2),ep)
     &    -q4(3)*FE2(y2(3,n2),ep)
      trhs=
     & -0.5d0*(FD14(n2,ep)-FD15a(n2,ep)+f4*FE1(n2,ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FE2'
      do ep=epmin,0
      tq = 
     & +FE2(y2(4,4),ep)
     & -FE2(y2(1,1),ep)
     & -FE2(y2(2,2),ep)
     & -FE2(y2(3,3),ep)
     & -m1s*FE0(ep)
      trhs=-FD05(ep)
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo


      elseif (rank .eq. 3) then
      if (pvverbose) write(6,*) 'q1.FE3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq = q1(4)*FE3(y3(4,n2,n3),ep)
     &    -q1(1)*FE3(y3(1,n2,n3),ep)
     &    -q1(2)*FE3(y3(2,n2,n3),ep)
     &    -q1(3)*FE3(y3(3,n2,n3),ep)
      trhs=
     &   -0.5d0*(FD21(y2(n2,n3),ep)
     & -FD25a(y2(n2,n3),ep)+f1*FE2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q2.FE3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq =q2(4)*FE3(y3(4,n2,n3),ep)
     &   -q2(1)*FE3(y3(1,n2,n3),ep)
     &   -q2(2)*FE3(y3(2,n2,n3),ep)
     &   -q2(3)*FE3(y3(3,n2,n3),ep)   
      trhs=
     & -0.5d0*(FD22(y2(n2,n3),ep)
     & -FD25a(y2(n2,n3),ep)+f2*FE2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q3.FE3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq =q3(4)*FE3(y3(4,n2,n3),ep)
     &   -q3(1)*FE3(y3(1,n2,n3),ep)
     &   -q3(2)*FE3(y3(2,n2,n3),ep)
     &   -q3(3)*FE3(y3(3,n2,n3),ep)   
      trhs=
     & -0.5d0*(FD23(y2(n2,n3),ep)
     & -FD25a(y2(n2,n3),ep)+f3*FE2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q4.FE3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq =q4(4)*FE3(y3(4,n2,n3),ep)
     &   -q4(1)*FE3(y3(1,n2,n3),ep)
     &   -q4(2)*FE3(y3(2,n2,n3),ep)
     &   -q4(3)*FE3(y3(3,n2,n3),ep)   
      trhs=
     & -0.5d0*(FD24(y2(n2,n3),ep)
     & -FD25a(y2(n2,n3),ep)+f4*FE2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FE3'
c      do ep=epmin,0
      ep=0
      do n3=1,4
      tq =    
     & +FE3(y3(4,4,n3),ep)
     & -FE3(y3(1,1,n3),ep)
     & -FE3(y3(2,2,n3),ep)
     & -FE3(y3(3,3,n3),ep)
     & -m1s*FE1(n3,ep)
      trhs=
     & -FD15a(n3,ep)
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
c      enddo


      elseif (rank .eq. 4) then
      if (pvverbose) write(6,*) 'q1.FE4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q1(4)*FE4(y4(4,n2,n3,n4),ep)
     &   -q1(1)*FE4(y4(1,n2,n3,n4),ep)
     &   -q1(2)*FE4(y4(2,n2,n3,n4),ep)
     &   -q1(3)*FE4(y4(3,n2,n3,n4),ep)
      trhs=
     & -0.5d0*(FD31(y3(n2,n3,n4),ep)
     & -FD35a(y3(n2,n3,n4),ep)+f1*FE3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
      enddo
      enddo
     
      if (pvverbose) write(6,*) 'q2.FE4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q2(4)*FE4(y4(4,n2,n3,n4),ep)
     &   -q2(1)*FE4(y4(1,n2,n3,n4),ep)
     &   -q2(2)*FE4(y4(2,n2,n3,n4),ep)
     &   -q2(3)*FE4(y4(3,n2,n3,n4),ep)
      trhs=
     & -0.5d0*(FD32(y3(n2,n3,n4),ep)
     & -FD35a(y3(n2,n3,n4),ep)+f2*FE3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
      enddo
      enddo
     
      if (pvverbose) write(6,*) 'q3.FE4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q3(4)*FE4(y4(4,n2,n3,n4),ep)
     &   -q3(1)*FE4(y4(1,n2,n3,n4),ep)
     &   -q3(2)*FE4(y4(2,n2,n3,n4),ep)
     &   -q3(3)*FE4(y4(3,n2,n3,n4),ep)
      trhs=
     &   -0.5d0*(FD33(y3(n2,n3,n4),ep)
     & -FD35a(y3(n2,n3,n4),ep)+f3*FE3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
      enddo
      enddo
    
      if (pvverbose) write(6,*) 'q4.FE4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q4(4)*FE4(y4(4,n2,n3,n4),ep)
     &   -q4(1)*FE4(y4(1,n2,n3,n4),ep)
     &   -q4(2)*FE4(y4(2,n2,n3,n4),ep)
     &   -q4(3)*FE4(y4(3,n2,n3,n4),ep)
      trhs=
     &   -0.5d0*(FD34(y3(n2,n3,n4),ep)
     & -FD35a(y3(n2,n3,n4),ep)+f4*FE3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Eacc,failed) 
      enddo
      enddo
      enddo
      enddo
 
c--- This test needs to be thought about some more    
c      if (pvverbose) write(6,*) 'g_(mu,nu)*FE4'
c      do ep=epmin,0
c      do n3=1,4
c      do n4=n3,4
c      sing4(0)=-1d0/12d0*g(n3,n4)
c      tq = +FE4(y4(4,4,n3,n4),ep)
c     & -FE4(y4(1,1,n3,n4),ep)
c     & -FE4(y4(2,2,n3,n4),ep)
c     & -FE4(y4(3,3,n3,n4),ep)
c     & -m1s*FE2(y2(n3,n4),ep)
c      trhs= tq
c     & -FD25a(y2(n3,n4),ep)
c     & +dcmplx(sing4(ep))
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      enddo
c      enddo
c      enddo


      elseif (rank .eq. 5) then
      
      write(6,*) 'No check available for rank 5 pentagon.'
      stop
      
c      if (pvverbose) write(6,*) 'q1.FE5'
c      do ep=epmin,0
c      do n2=1,4
c      do n3=n2,4
c      do n4=n3,4
c      do n5=n4,4
c      tq=+q1(4)*FE5(y5(4,n2,n3,n4,n5),ep)
c     &   -q1(1)*FE5(y5(1,n2,n3,n4,n5),ep)
c     &   -q1(2)*FE5(y5(2,n2,n3,n4,n5),ep)
c     &   -q1(3)*FE5(y5(3,n2,n3,n4,n5),ep)
c      trhs=
c     & -0.5d0*(FD41(y4(n2,n3,n4,n5),ep)
c     & -FD45a(y4(n2,n3,n4,n5),ep)+f1*FE4(y4(n2,n3,n4,n5),ep))
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c      if (pvverbose) write(6,*) 'q2.FE5'
c      do ep=epmin,0
c      do n2=1,4
c      do n3=n2,4
c      do n4=n3,4
c      do n5=n4,4
c      tq=+q2(4)*FE5(y5(4,n2,n3,n4,n5),ep)
c     &   -q2(1)*FE5(y5(1,n2,n3,n4,n5),ep)
c     &   -q2(2)*FE5(y5(2,n2,n3,n4,n5),ep)
c     &   -q2(3)*FE5(y5(3,n2,n3,n4,n5),ep)
c      trhs=
c     & -0.5d0*(FD42(y4(n2,n3,n4,n5),ep)
c     & -FD45a(y4(n2,n3,n4,n5),ep)+f2*FE4(y4(n2,n3,n4,n5),ep))
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c      if (pvverbose) write(6,*) 'q3.FE5'
c      do ep=epmin,0
c      do n2=1,4
c      do n3=n2,4
c      do n4=n3,4
c      do n5=n4,4
c      tq=+q3(4)*FE5(y5(4,n2,n3,n4,n5),ep)
c     &   -q3(1)*FE5(y5(1,n2,n3,n4,n5),ep)
c     &   -q3(2)*FE5(y5(2,n2,n3,n4,n5),ep)
c     &   -q3(3)*FE5(y5(3,n2,n3,n4,n5),ep)
c      trhs=
c     & -0.5d0*(FD43(y4(n2,n3,n4,n5),ep)
c     & -FD45a(y4(n2,n3,n4,n5),ep)+f3*FE4(y4(n2,n3,n4,n5),ep))
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c      if (pvverbose) write(6,*) 'q4.FE5'
c      do ep=epmin,0
c      do n2=1,4
c      do n3=n2,4
c      do n4=n3,4
c      do n5=n4,4
c      tq=+q4(4)*FE5(y5(4,n2,n3,n4,n5),ep)
c     &   -q4(1)*FE5(y5(1,n2,n3,n4,n5),ep)
c     &   -q4(2)*FE5(y5(2,n2,n3,n4,n5),ep)
c     &   -q4(3)*FE5(y5(3,n2,n3,n4,n5),ep)
c      trhs=
c     & -0.5d0*(FD44(y4(n2,n3,n4,n5),ep)
c     & -FD45a(y4(n2,n3,n4,n5),ep)+f4*FE4(y4(n2,n3,n4,n5),ep))
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      write(6,*) tq,tq+trhs
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c--- This test needs to be thought about some more    
c      if (pvverbose) write(6,*) 'g_(mu,nu)*FE5'
c      do ep=epmin,0
c      do n3=1,4
c      do n4=n3,4
c      do n5=n4,4
c      sing5(0)=
c     & +(g(n3,n4)*q1(n5)+g(n4,n5)*q1(n3)+g(n5,n3)*q1(n4))/48d0
c     & +(g(n3,n4)*q2(n5)+g(n4,n5)*q2(n3)+g(n5,n3)*q2(n4))/48d0
c     & +(g(n3,n4)*q3(n5)+g(n4,n5)*q3(n3)+g(n5,n3)*q3(n4))/48d0
c      tq = 
c     & +FE5(y5(4,4,n3,n4,n5),ep)
c     & -FE5(y5(1,1,n3,n4,n5),ep)
c     & -FE5(y5(2,2,n3,n4,n5),ep)
c     & -FE5(y5(3,3,n3,n4,n5),ep)
c     & -m1s*FE3(y3(n3,n4,n5),ep)
c      trhs= tq 
c     & -FD34a(y3(n3,n4,n5),ep)
c     & +dcmplx(sing5(ep))  
c      call checkaccuracy(trhs,tq,Eacc,failed) 
c      enddo
c      enddo
c      enddo
c      enddo
     
      endif
      return
      
      end
      
      
