      subroutine ovCcheck(rank,q1,q2,m0s,m1s,m2s,
     & FC0,FC1,FC2,FC3,failed)
      implicit none
      include 'TRconstants.f'
      include 'TRydef.f'
      include 'TRmetric.f'
      include 'pvverbose.f'
      integer n2,n3,ep,nu,rank,epmin
      double precision q1(4),q2(4),p2(4),f1,f2,Cacc
      double precision q1Dq1,q2Dq2,q1Dq2,s12,m0s,m1s,m2s
      double precision sing2(-2:0),sing3(-2:0)
      double complex FB01(-2:0),FB11(y1max,-2:0),FB21(y2max,-2:0)
      double complex FB02(-2:0),FB12(y1max,-2:0),FB22(y2max,-2:0)
      double complex FB03(-2:0),FB13(y1max,-2:0),FB23(y2max,-2:0)
      double complex FB13a(y1max,-2:0),FB23a(y2max,-2:0)
      double complex B00(-2:0)
      double complex FC0(-2:0),FC1(y1max,-2:0),FC2(y2max,-2:0),
     & FC3(y3max,-2:0),trhs,tq
      logical failed
      integer ierr
      parameter(epmin=0) ! Only check finite pieces
      
      failed=.false.
      ierr=0
      
      Cacc=1d-8

      q1Dq1=q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2
      q2Dq2=q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2
      q1Dq2=q1(4)*q2(4)-q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)
      s12=q1Dq1+q2Dq2-2d0*q1Dq2

c      write(6,'(a35,5(e12.5,a2))') 
c     . '(p1sq, p2sq, m1sq, m2sq, m3sq) = ( ',
c     .  q1Dq1,', ',s12,', ',m0s,', ',m1s,', ',m2s,' )'
      
      do ep=epmin,0
      sing2(ep)=zip
      sing3(ep)=zip
      enddo
      do nu=1,4
      p2(nu)=q2(nu)-q1(nu)
      enddo  

      call ovBtensor(q2,m0s,m2s,FB01,FB11,FB21,B00)
      call ovBtensor(q1,m0s,m1s,FB02,FB12,FB22,B00)
      call ovBtensor(p2,m1s,m2s,FB03,FB13,FB23,B00)

      if ((rank .eq. 2) .or. (rank .eq. 3))
     &  call pvswitch1(q1,FB03,FB13,FB13a)
      if (rank .eq. 3)
     &  call pvswitch2(q1,FB03,FB13,FB23,FB23a)

      f1=m1s-m0s-q1Dq1
      f2=m2s-m0s-q2Dq2

c--- check rank 1
      if (rank .eq. 1) then
      if (pvverbose) write(6,*) 'q1.FC1'
      do ep=epmin,0
      tq=  q1(4)*FC1(4,ep)
     &    -q1(1)*FC1(1,ep)
     &    -q1(2)*FC1(2,ep)
     &    -q1(3)*FC1(3,ep)
      trhs=
     & -0.5d0*(FB01(ep)-FB03(ep)+f1*FC0(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=11
      goto 77
      endif
      enddo

      if (pvverbose) write(6,*) 'q2.FC1'
      do ep=epmin,0
      tq=  q2(4)*FC1(4,ep)
     &    -q2(1)*FC1(1,ep)
     &    -q2(2)*FC1(2,ep)
     &    -q2(3)*FC1(3,ep)
      trhs=
     & -0.5d0*(FB02(ep)-FB03(ep)+f2*FC0(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=12
      goto 77
      endif
      enddo
      
      endif
      

c--- check for rank 2
      if (rank .eq. 2) then
      if (pvverbose) write(6,*) 'q1.FC2'
      do ep=epmin,0
      do n2=1,4
      tq=  q1(4)*FC2(y2(4,n2),ep)
     &    -q1(1)*FC2(y2(1,n2),ep)
     &    -q1(2)*FC2(y2(2,n2),ep)
     &    -q1(3)*FC2(y2(3,n2),ep)
      trhs=
     & -0.5d0*(FB11(n2,ep)-FB13a(n2,ep)+f1*FC1(n2,ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=21
      goto 77
      endif
      enddo
      enddo

      if (pvverbose) write(6,*) 'q2.FC2'
      do ep=epmin,0
      do n2=1,4
      tq= q2(4)*FC2(y2(4,n2),ep)
     &   -q2(1)*FC2(y2(1,n2),ep)
     &   -q2(2)*FC2(y2(2,n2),ep)
     &   -q2(3)*FC2(y2(3,n2),ep)
      trhs=
     & -0.5d0*(FB12(n2,ep)-FB13a(n2,ep)+f2*FC1(n2,ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=22
      goto 77
      endif
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC2'
      sing2(0)=-0.5d0 
      do ep=epmin,0 
      tq=FC2(y2(4,4),ep)
     & -FC2(y2(1,1),ep)
     & -FC2(y2(2,2),ep)
     & -FC2(y2(3,3),ep)
     & -m0s*FC0(ep)-FB03(ep)
      trhs=
     & +dcmplx(sing2(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=20
      goto 77
      endif
      enddo

      endif
      

c--- check for rank 2
      if (rank .eq. 3) then
      if (pvverbose) write(6,*) 'q1.FC3'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      tq=  +q1(4)*FC3(y3(4,n2,n3),ep)
     &     -q1(1)*FC3(y3(1,n2,n3),ep)
     &     -q1(2)*FC3(y3(2,n2,n3),ep)
     &     -q1(3)*FC3(y3(3,n2,n3),ep)
      trhs=
     &    -0.5d0*(FB21(y2(n2,n3),ep)
     &           -FB23a(y2(n2,n3),ep)+f1*FC2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=31
      goto 77
      endif
      enddo
      enddo
      enddo
      
      if (pvverbose) write(6,*) 'q2.FC3'
      do ep=epmin,0
      do n2= 1,4
      do n3=n2,4
      tq=+q2(4)*FC3(y3(4,n2,n3),ep)
     &     -q2(1)*FC3(y3(1,n2,n3),ep)
     &     -q2(2)*FC3(y3(2,n2,n3),ep)
     &     -q2(3)*FC3(y3(3,n2,n3),ep)
      trhs=
     &    -0.5d0*(FB22(y2(n2,n3),ep)
     &           -FB23a(y2(n2,n3),ep)+f2*FC2(y2(n2,n3),ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=32
      goto 77
      endif
      enddo
      enddo
      enddo

      if (pvverbose) write(6,*) 'g_(mu,nu)*FC3'
      do ep=epmin,0
      do n3=1,4
      sing3(0)=+1d0/6d0*(q1(n3)+q2(n3))
 
      tq=FC3(y3(4,4,n3),ep)
     & -FC3(y3(1,1,n3),ep)
     & -FC3(y3(2,2,n3),ep)
     & -FC3(y3(3,3,n3),ep)
     & -m0s*FC1(n3,ep)-FB13a(n3,ep)
      trhs=
     & +dcmplx(sing3(ep))
      call checkaccuracy(trhs,tq,Cacc,failed) 
      if (failed) then
        ierr=30
      goto 77
      endif
      enddo
      enddo
  
      endif
      
   77 continue
c--- print out error message if necessary and return
      if (failed) then
c        write(6,*) 'ovCcheck: error code',ierr
c        write(6,*) 'tq,trhs,tq+trhs',tq,trhs,
c     &   (tq+trhs)/(tq-trhs)
      endif
   
      return
      end
