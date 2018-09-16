      subroutine pvDcheck(rank,q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4,FD5,FD6,failed)
      implicit none
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'pvverbose.f'
      integer n2,n3,n4,n5,n6,ep,nu,rank,epmin
      double precision q1(4),q2(4),q3(4),p2(4),p23(4),Dacc
      double precision q1Dq1,q2Dq2,q3Dq3,q1Dq2,q2Dq3,q1Dq3,pvSDDDD,
     & pvSDDPP,pvSDDPK,s12,s13,s23,m1s,m2s,m3s,m4s,
     & sing4(-2:0),sing5(-2:0),sing6(-2:0),f1,f2,f3
      double complex 
     & FC01(-2:0),FC11(y1max,-2:0),FC21(y2max,-2:0),FC31(y3max,-2:0),
     & FC41(y4max,-2:0),FC51(y5max,-2:0),FC61(y6max,-2:0)
      double complex 
     & FC02(-2:0),FC12(y1max,-2:0),FC22(y2max,-2:0),FC32(y3max,-2:0),
     & FC42(y4max,-2:0),FC52(y5max,-2:0),FC62(y6max,-2:0)
      double complex 
     & FC03(-2:0),FC13(y1max,-2:0),FC23(y2max,-2:0),FC33(y3max,-2:0),
     & FC43(y4max,-2:0),FC53(y5max,-2:0),FC63(y6max,-2:0)
      double complex 
     & FC04(-2:0),FC14(y1max,-2:0),FC24(y2max,-2:0),FC34(y3max,-2:0),
     & FC44(y4max,-2:0),FC54(y5max,-2:0),FC64(y6max,-2:0)
      double complex FC14a(y1max,-2:0),FC24a(y2max,-2:0),
     & FC34a(y3max,-2:0),FC44a(y4max,-2:0),FC54a(y5max,-2:0)
      double complex FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     & FD3(y3max,-2:0),FD4(y4max,-2:0),FD5(y5max,-2:0),
     & FD6(y6max,-2:0),trhs,tq
      logical failed
      include 'TRmetric.f'
      parameter(epmin=0) ! Only check finite pieces
      
      failed=.false.

      Dacc=1d-8
      
      q1Dq1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2Dq2=q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2
      q3Dq3=q3(4)**2-q3(1)**2-q3(2)**2-q3(3)**2
      q1Dq2=q1(4)*q2(4)-q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)
      q1Dq3=q1(4)*q3(4)-q1(1)*q3(1)-q1(2)*q3(2)-q1(3)*q3(3)
      q2Dq3=q2(4)*q3(4)-q2(1)*q3(1)-q2(2)*q3(2)-q2(3)*q3(3)

c      if (pvverbose) write(6,*)
c      if (pvverbose) write(6,*) '(p1sq, p2sq, p3sq, m1sq, m2sq, m3sq, m4sq) = '
c      if (pvverbose) write(6,'(a2,7(e12.5,a2))'),
c     . '( ',q1Dq1,', ',q2Dq2+q1Dq1-2d0*q1Dq2,', ',q3Dq3+q2Dq2-2d0*q2Dq3,
c     . ', ',m1s,', ',m2s,', ',m3s,', ',m4s,' )' 
            
      do nu=1,4
      p2(nu)=q2(nu)-q1(nu)
      p23(nu)=q3(nu)-q1(nu)
      enddo  

      do ep=-2,-1
      sing4(ep)=zip
      sing5(ep)=zip
      sing6(ep)=zip
      enddo
      call pvCtensor(p2,p23,m2s,m3s,m4s,
     & FC04,FC14,FC24,FC34,FC44,FC54,FC64)
      call pvCtensor(q2,q3,m1s,m3s,m4s,
     & FC01,FC11,FC21,FC31,FC41,FC51,FC61)
      call pvCtensor(q1,q3,m1s,m2s,m4s,
     & FC02,FC12,FC22,FC32,FC42,FC52,FC62)
      call pvCtensor(q1,q2,m1s,m2s,m3s,
     & FC03,FC13,FC23,FC33,FC43,FC53,FC63)

      if ((rank .eq. 2) .or. (rank .eq. 3))
     &  call pvswitch1(q1,FC04,FC14,FC14a)
      if ((rank .eq. 3) .or. (rank .eq. 4))
     &  call pvswitch2(q1,FC04,FC14,FC24,FC24a)
      if ((rank .eq. 4) .or. (rank .eq. 5))
     &  call pvswitch3(q1,FC04,FC14,FC24,FC34,FC34a)
      if ((rank .eq. 5) .or. (rank .eq. 6))
     &  call pvswitch4(q1,FC04,FC14,FC24,FC34,FC44,FC44a)
      if (rank .eq. 6)
     &  call pvswitch5(q1,FC04,FC14,FC24,FC34,FC44,FC54,FC54a)

      f1=m2s-m1s-q1Dq1
      f2=m3s-m1s-q2Dq2
      f3=m4s-m1s-q3Dq3

      if (rank .eq. 1) then

      if (pvverbose) write(6,*) 'q1.FD1'
      do ep=epmin,0
      tq=q1(4)*FD1(4,ep)
     &  -q1(1)*FD1(1,ep)
     &  -q1(2)*FD1(2,ep)
     &  -q1(3)*FD1(3,ep)
      trhs=-0.5d0*(FC01(ep)-FC04(ep)+f1*FD0(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
c      if (pvverbose) write(6,*) tq
c      if (pvverbose) write(6,*) -0.5d0*(FC01(ep))
c      if (pvverbose) write(6,*) -0.5d0*(-FC04(ep))
c      if (pvverbose) write(6,*) -0.5d0*(+f1*FD0(ep))
      enddo

      if (pvverbose) write(6,*) 'q2.FD1'
      do ep=epmin,0
      tq=q2(4)*FD1(4,ep)
     &  -q2(1)*FD1(1,ep)
     &  -q2(2)*FD1(2,ep)
     &  -q2(3)*FD1(3,ep)
      trhs=-0.5d0*(FC02(ep)-FC04(ep)+f2*FD0(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo

      if (pvverbose) write(6,*) 'q3.FD1'
      do ep=epmin,0
      tq=q3(4)*FD1(4,ep)
     &    -q3(1)*FD1(1,ep)
     &    -q3(2)*FD1(2,ep)
     &    -q3(3)*FD1(3,ep)
      trhs=-0.5d0*(FC03(ep)-FC04(ep)+f3*FD0(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo


      elseif (rank .eq. 2) then


      if (pvverbose) write(6,*) 'q1.FD2'
      do ep=epmin,0
      do n2=1,4
      tq =q1(4)*FD2(y2(4,n2),ep)
     &   -q1(1)*FD2(y2(1,n2),ep)
     &   -q1(2)*FD2(y2(2,n2),ep)
     &   -q1(3)*FD2(y2(3,n2),ep)   
      trhs= 
     &   -0.5d0*(FC11(n2,ep)-FC14a(n2,ep)+f1*FD1(n2,ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FD2'
      do ep=epmin,0
      do n2=1,4
      tq =q2(4)*FD2(y2(4,n2),ep)
     &   -q2(1)*FD2(y2(1,n2),ep)
     &   -q2(2)*FD2(y2(2,n2),ep)
     &   -q2(3)*FD2(y2(3,n2),ep)
      trhs=
     &   -0.5d0*(FC12(n2,ep)-FC14a(n2,ep)+f2*FD1(n2,ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'q3.FD2'
      do ep=epmin,0
      do n2=1,4
      tq = q3(4)*FD2(y2(4,n2),ep)
     &    -q3(1)*FD2(y2(1,n2),ep)
     &    -q3(2)*FD2(y2(2,n2),ep)
     &    -q3(3)*FD2(y2(3,n2),ep)
      trhs=
     & -0.5d0*(FC13(n2,ep)-FC14a(n2,ep)+f3*FD1(n2,ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FD2'
      do ep=epmin,0
      tq = 
     & +FD2(y2(4,4),ep)
     & -FD2(y2(1,1),ep)
     & -FD2(y2(2,2),ep)
     & -FD2(y2(3,3),ep)
     & -m1s*FD0(ep)
      trhs=-FC04(ep)
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo



      elseif (rank .eq. 3) then
      if (pvverbose) write(6,*) 'q1.FD3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq = q1(4)*FD3(y3(4,n2,n3),ep)
     &    -q1(1)*FD3(y3(1,n2,n3),ep)
     &    -q1(2)*FD3(y3(2,n2,n3),ep)
     &    -q1(3)*FD3(y3(3,n2,n3),ep)
      trhs=
     &   -0.5d0*(FC21(y2(n2,n3),ep)
     & -FC24a(y2(n2,n3),ep)+f1*FD2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q2.FD3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq =q2(4)*FD3(y3(4,n2,n3),ep)
     &   -q2(1)*FD3(y3(1,n2,n3),ep)
     &   -q2(2)*FD3(y3(2,n2,n3),ep)
     &   -q2(3)*FD3(y3(3,n2,n3),ep)   
      trhs=
     & -0.5d0*(FC22(y2(n2,n3),ep)
     & -FC24a(y2(n2,n3),ep)+f2*FD2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'q3.FD3'
c      do ep=epmin,0
      ep=0
      do n2=1,4
      do n3=n2,4
      tq =q3(4)*FD3(y3(4,n2,n3),ep)
     &   -q3(1)*FD3(y3(1,n2,n3),ep)
     &   -q3(2)*FD3(y3(2,n2,n3),ep)
     &   -q3(3)*FD3(y3(3,n2,n3),ep)   
      trhs=
     & -0.5d0*(FC23(y2(n2,n3),ep)
     & -FC24a(y2(n2,n3),ep)+f3*FD2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
c      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FD3'
c      do ep=epmin,0
      ep=0
      do n3=1,4
      tq =    
     & +FD3(y3(4,4,n3),ep)
     & -FD3(y3(1,1,n3),ep)
     & -FD3(y3(2,2,n3),ep)
     & -FD3(y3(3,3,n3),ep)
     & -m1s*FD1(n3,ep)
      trhs=
     & -FC14a(n3,ep)
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
c      enddo

      elseif (rank .eq. 4) then
      if (pvverbose) write(6,*) 'q1.FD4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q1(4)*FD4(y4(4,n2,n3,n4),ep)
     &   -q1(1)*FD4(y4(1,n2,n3,n4),ep)
     &   -q1(2)*FD4(y4(2,n2,n3,n4),ep)
     &   -q1(3)*FD4(y4(3,n2,n3,n4),ep)
      trhs=
     & -0.5d0*(FC31(y3(n2,n3,n4),ep)
     & -FC34a(y3(n2,n3,n4),ep)+f1*FD3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
     
      if (pvverbose) write(6,*) 'q2.FD4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q2(4)*FD4(y4(4,n2,n3,n4),ep)
     &   -q2(1)*FD4(y4(1,n2,n3,n4),ep)
     &   -q2(2)*FD4(y4(2,n2,n3,n4),ep)
     &   -q2(3)*FD4(y4(3,n2,n3,n4),ep)
      trhs=
     & -0.5d0*(FC32(y3(n2,n3,n4),ep)
     & -FC34a(y3(n2,n3,n4),ep)+f2*FD3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
     
      if (pvverbose) write(6,*) 'q3.FD4'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      tq= q3(4)*FD4(y4(4,n2,n3,n4),ep)
     &   -q3(1)*FD4(y4(1,n2,n3,n4),ep)
     &   -q3(2)*FD4(y4(2,n2,n3,n4),ep)
     &   -q3(3)*FD4(y4(3,n2,n3,n4),ep)
      trhs=
     &   -0.5d0*(FC33(y3(n2,n3,n4),ep)
     & -FC34a(y3(n2,n3,n4),ep)+f3*FD3(y3(n2,n3,n4),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
    
      if (pvverbose) write(6,*) 'g_(mu,nu)*FD4'
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      sing4(0)=-1d0/12d0*g(n3,n4)
      tq = +FD4(y4(4,4,n3,n4),ep)
     & -FD4(y4(1,1,n3,n4),ep)
     & -FD4(y4(2,2,n3,n4),ep)
     & -FD4(y4(3,3,n3,n4),ep)
     & -m1s*FD2(y2(n3,n4),ep)
      trhs=
     & -FC24a(y2(n3,n4),ep)
     & +dcmplx(sing4(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo

      elseif (rank .eq. 5) then
      if (pvverbose) write(6,*) 'q1.FD5'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      tq=+q1(4)*FD5(y5(4,n2,n3,n4,n5),ep)
     &   -q1(1)*FD5(y5(1,n2,n3,n4,n5),ep)
     &   -q1(2)*FD5(y5(2,n2,n3,n4,n5),ep)
     &   -q1(3)*FD5(y5(3,n2,n3,n4,n5),ep)
      trhs=
     & -0.5d0*(FC41(y4(n2,n3,n4,n5),ep)
     & -FC44a(y4(n2,n3,n4,n5),ep)+f1*FD4(y4(n2,n3,n4,n5),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FD5'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      tq=+q2(4)*FD5(y5(4,n2,n3,n4,n5),ep)
     &   -q2(1)*FD5(y5(1,n2,n3,n4,n5),ep)
     &   -q2(2)*FD5(y5(2,n2,n3,n4,n5),ep)
     &   -q2(3)*FD5(y5(3,n2,n3,n4,n5),ep)
      trhs=
     & -0.5d0*(FC42(y4(n2,n3,n4,n5),ep)
     & -FC44a(y4(n2,n3,n4,n5),ep)+f2*FD4(y4(n2,n3,n4,n5),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q3.FD5'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      tq=+q3(4)*FD5(y5(4,n2,n3,n4,n5),ep)
     &   -q3(1)*FD5(y5(1,n2,n3,n4,n5),ep)
     &   -q3(2)*FD5(y5(2,n2,n3,n4,n5),ep)
     &   -q3(3)*FD5(y5(3,n2,n3,n4,n5),ep)
      trhs=
     & -0.5d0*(FC43(y4(n2,n3,n4,n5),ep)
     & -FC44a(y4(n2,n3,n4,n5),ep)+f3*FD4(y4(n2,n3,n4,n5),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FD5'
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      do n5=n4,4
      sing5(0)=
     & +(g(n3,n4)*q1(n5)+g(n4,n5)*q1(n3)+g(n5,n3)*q1(n4))/48d0
     & +(g(n3,n4)*q2(n5)+g(n4,n5)*q2(n3)+g(n5,n3)*q2(n4))/48d0
     & +(g(n3,n4)*q3(n5)+g(n4,n5)*q3(n3)+g(n5,n3)*q3(n4))/48d0
      tq = 
     & +FD5(y5(4,4,n3,n4,n5),ep)
     & -FD5(y5(1,1,n3,n4,n5),ep)
     & -FD5(y5(2,2,n3,n4,n5),ep)
     & -FD5(y5(3,3,n3,n4,n5),ep)
     & -m1s*FD3(y3(n3,n4,n5),ep)
      trhs=
     & -FC34a(y3(n3,n4,n5),ep)
     & +dcmplx(sing5(ep))  
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
     

      elseif (rank .eq. 6) then
      if (pvverbose) write(6,*) 'q1.FD6'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      tq=q1(4)*FD6(y6(4,n2,n3,n4,n5,n6),ep)
     &    -q1(1)*FD6(y6(1,n2,n3,n4,n5,n6),ep)
     &    -q1(2)*FD6(y6(2,n2,n3,n4,n5,n6),ep)
     &    -q1(3)*FD6(y6(3,n2,n3,n4,n5,n6),ep)

      trhs=
     &   -0.5d0*(FC51(y5(n2,n3,n4,n5,n6),ep)
     & -FC54a(y5(n2,n3,n4,n5,n6),ep)+f1*FD5(y5(n2,n3,n4,n5,n6),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      
      if (pvverbose) write(6,*) 'q2.FD6'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      tq =q2(4)*FD6(y6(4,n2,n3,n4,n5,n6),ep)
     &   -q2(1)*FD6(y6(1,n2,n3,n4,n5,n6),ep)
     &   -q2(2)*FD6(y6(2,n2,n3,n4,n5,n6),ep)
     &   -q2(3)*FD6(y6(3,n2,n3,n4,n5,n6),ep)
   
      trhs=
     & -0.5d0*(FC52(y5(n2,n3,n4,n5,n6),ep)
     & -FC54a(y5(n2,n3,n4,n5,n6),ep)+f2*FD5(y5(n2,n3,n4,n5,n6),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      
      if (pvverbose) write(6,*) 'q3.FD6'
      do ep=epmin,0
      do n2=1,4
      do n3=n2,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4
      tq = q3(4)*FD6(y6(4,n2,n3,n4,n5,n6),ep)
     &   -q3(1)*FD6(y6(1,n2,n3,n4,n5,n6),ep)
     &   -q3(2)*FD6(y6(2,n2,n3,n4,n5,n6),ep)
     &   -q3(3)*FD6(y6(3,n2,n3,n4,n5,n6),ep) 
      trhs=
     & -0.5d0*(FC53(y5(n2,n3,n4,n5,n6),ep)
     & -FC54a(y5(n2,n3,n4,n5,n6),ep)+f3*FD5(y5(n2,n3,n4,n5,n6),ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


      if (pvverbose) write(6,*) 'g_(mu,nu)*FD6'
      do ep=epmin,0
      do n3=1,4
      do n4=n3,4
      do n5=n4,4
      do n6=n5,4


      s12=q1Dq1+q2Dq2-2d0*q1Dq2
      s13=q1Dq1+q3Dq3-2d0*q1Dq3
      s23=q2Dq2+q3Dq3-2d0*q2Dq3
      sing6(0)=
     & +pvSDDDD(n3,n4,n5,n6)
     & *((s12+s23+s13+q1Dq1+q2Dq2+q3Dq3)/480d0-(m1s+m2s+m3s+m4s)/96d0)
     & -pvSDDPP(n3,n4,n5,n6,q1)/120D0
     & -pvSDDPP(n3,n4,n5,n6,q2)/120D0
     & -pvSDDPP(n3,n4,n5,n6,q3)/120D0
     & -pvSDDPK(n3,n4,n5,n6,q1,q2)/240D0
     & -pvSDDPK(n3,n4,n5,n6,q2,q3)/240D0
     & -pvSDDPK(n3,n4,n5,n6,q3,q1)/240D0
      tq = 
     & +FD6(y6(4,4,n3,n4,n5,n6),ep)
     & -FD6(y6(1,1,n3,n4,n5,n6),ep)
     & -FD6(y6(2,2,n3,n4,n5,n6),ep)
     & -FD6(y6(3,3,n3,n4,n5,n6),ep)
     & -m1s*FD4(y4(n3,n4,n5,n6),ep)
     & -FC44a(y4(n3,n4,n5,n6),ep)


      trhs=+dcmplx(sing6(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      enddo
      enddo
      enddo
      enddo
      enddo

      endif

      return
      end


