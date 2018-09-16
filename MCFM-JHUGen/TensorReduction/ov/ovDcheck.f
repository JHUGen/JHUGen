      subroutine ovDcheck(rank,q1,q2,q3,m1s,m2s,m3s,m4s,
     & FD0,FD1,FD2,FD3,FD4,failed)
      implicit none
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'pvverbose.f'
      include 'TRmetric.f'
      integer n2,n3,n4,ep,nu,rank,epmin
      double precision q1(4),q2(4),q3(4),p2(4),p3(4),p23(4),Dacc
      double precision q1Dq1,q2Dq2,q3Dq3,q1Dq2,q2Dq3,q1Dq3,
     & m1s,m2s,m3s,m4s,sing4(-2:0),f1,f2,f3
      double complex 
     & FC01(-2:0),FC11(y1max,-2:0),FC21(y2max,-2:0),FC31(y3max,-2:0),
     & FC02(-2:0),FC12(y1max,-2:0),FC22(y2max,-2:0),FC32(y3max,-2:0),
     & FC03(-2:0),FC13(y1max,-2:0),FC23(y2max,-2:0),FC33(y3max,-2:0),
     & FC04(-2:0),FC14(y1max,-2:0),FC24(y2max,-2:0),FC34(y3max,-2:0)
      double complex FC14a(y1max,-2:0),FC24a(y2max,-2:0),
     & FC34a(y3max,-2:0)
      double complex FD0(-2:0),FD1(y1max,-2:0),FD2(y2max,-2:0),
     & FD3(y3max,-2:0),FD4(y4max,-2:0),trhs,tq,C00(-2:0),tau3(4,-2:0)
      logical failed
      integer ierr
      parameter(epmin=0) ! Only check finite pieces
      
      failed=.false.
      ierr=0

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
      p3(nu)=q3(nu)-q2(nu)
      enddo  

      do ep=-2,0
      sing4(ep)=zip
      enddo
      call ovCtensor(p2,p3,m2s,m3s,m4s,
     & FC04,FC14,FC24,FC34,C00,tau3)
      call ovCtensor(q2,p3,m1s,m3s,m4s,
     & FC01,FC11,FC21,FC31,C00,tau3)
      call ovCtensor(q1,p23,m1s,m2s,m4s,
     & FC02,FC12,FC22,FC32,C00,tau3)
      call ovCtensor(q1,p2,m1s,m2s,m3s,
     & FC03,FC13,FC23,FC33,C00,tau3)

      if ((rank .eq. 2) .or. (rank .eq. 3))
     &  call pvswitch1(q1,FC04,FC14,FC14a)
      if ((rank .eq. 3) .or. (rank .eq. 4))
     &  call pvswitch2(q1,FC04,FC14,FC24,FC24a)
      if (rank .eq. 4)
     &  call pvswitch3(q1,FC04,FC14,FC24,FC34,FC34a)

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
      if (failed) then
        ierr=11
      goto 77
      endif
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
      if (failed) then
        ierr=12
      goto 77
      endif
      enddo

      if (pvverbose) write(6,*) 'q3.FD1'
      do ep=epmin,0
      tq=q3(4)*FD1(4,ep)
     &    -q3(1)*FD1(1,ep)
     &    -q3(2)*FD1(2,ep)
     &    -q3(3)*FD1(3,ep)
      trhs=-0.5d0*(FC03(ep)-FC04(ep)+f3*FD0(ep))
      call checkaccuracy(trhs,tq,Dacc,failed) 
      if (failed) then
        ierr=13
      goto 77
      endif
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
      if (failed) then
        ierr=21
      goto 77
      endif
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
      if (failed) then
        ierr=22
      goto 77
      endif
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
      if (failed) then
        ierr=23
      goto 77
      endif
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
      if (failed) then
        ierr=20
      goto 77
      endif
      enddo



      elseif (rank .eq. 3) then
      if (pvverbose) write(6,*) 'q1.FD3'
      do ep=epmin,0
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
      if (failed) then
        ierr=31
      goto 77
      endif
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FD3'
      do ep=epmin,0
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
      if (failed) then
        ierr=32
      goto 77
      endif
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'q3.FD3'
      do ep=epmin,0
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
      if (failed) then
        ierr=33
      goto 77
      endif
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FD3'
      do ep=epmin,0
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
      if (failed) then
        ierr=30
      goto 77
      endif
      enddo
      enddo

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
      if (failed) then
        ierr=41
      goto 77
      endif
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
      if (failed) then
        ierr=42
      goto 77
      endif
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
      if (failed) then
        ierr=43
      goto 77
      endif
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
      if (failed) then
        ierr=40
      goto 77
      endif
      enddo
      enddo
      enddo

      endif

   77 continue
c--- print out error message if necessary and return
      if (failed) then
c        write(6,*) 'ovDcheck: error code',ierr
c        write(6,*) 'tq,trhs,tq+trhs',tq,trhs,
c     &   (tq+trhs)/(tq-trhs)
      endif
   
      return
      end


