      subroutine qqb_z2jetx(p,msq,mqq,msqx,msqx_cs)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) --> Z +g(p5) +g(p6)
c                          |
c                          --> l(p3)+a(p4)
c                            
c--all momenta incoming
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'flags.f'
      include 'lc.f'
      include 'part.f'
      integer i,j,k,f,pq,pl,nquark,nup,ndo,
     .   j1,j2,j3,icol
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,faclo,
     .   qqbZgg2_cs(0:2,2,2),qbqZgg2_cs(0:2,2,2),
     .   qgZqg2_cs(0:2,2,2),gqZqg2_cs(0:2,2,2),
     .   qgZgq2_cs(0:2,2,2),gqZgq2_cs(0:2,2,2),
     .   qbgZqbg2_cs(0:2,2,2),gqbZqbg2_cs(0:2,2,2),
     .   qbgZgqb2_cs(0:2,2,2),gqbZgqb2_cs(0:2,2,2),
     .   ggZqbq2_cs(0:2,2,2),ggZqqb2_cs(0:2,2,2),ggtemp(0:2),
     .   qqbZgg2(2,2),
     .   qgZqg2(2,2),
c     .   gqbZgqb2(2,2),
     .   gqZgq2(2,2),
     .   qgZgq2(2,2),
c     .   gqbZqbg2(2,2),
     .   gqZqg2(2,2),
     .   ggZqbq2(2,2)
c     .   qbgZgqb2(2,2),
c     .   ggZqqb2(2,2),
c     .   qbgZqbg2(2,2),
c     .   qbqZgg2(2,2)
      double precision tup,tdo

      double complex qRb_a(2,2,2),qRb_b(2,2,2)
      double complex qqb_a(2,2,2),qqb_b(2,2,2),prop

      double complex qbq_a(2,2,2),qbq_b(2,2,2)
      double complex qbR_a(2,2,2),qbR_b(2,2,2)

      double complex qq_a(2,2,2),qq_b(2,2,2)
      double complex qR_a(2,2,2),qR_b(2,2,2)

      double complex qbRb_a(2,2,2),qbRb_b(2,2,2)
      double complex qbqb_a(2,2,2),qbqb_b(2,2,2)

      double complex qRb_ax(2,2,2),qRb_bx(2,2,2)
      double complex qqb_ax(2,2,2),qqb_bx(2,2,2)

      double complex qbR_ax(2,2,2),qbR_bx(2,2,2)
      double complex qbq_ax(2,2,2),qbq_bx(2,2,2)

      double precision mqq(0:2,fn:nf,fn:nf),msq0,msq1,msq2
      double precision msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      double precision msqx_cs(0:2,-nf:nf,-nf:nf)

      integer rcolourchoice

      logical rGflag

      integer,parameter::swap(2)=(/2,1/),swap1(0:2)=(/0,2,1/)

c--- if we're calculating the REAL or VIRT matrix elements, we
c--- need all the colour structures, but want to preserve
c--- the actual value of colourchoice
      if ((part .eq. 'real') .or. (part .eq. 'virt')) then
        rcolourchoice=colourchoice
        colourchoice=0
      endif     
c--- if we're calculating the REAL matrix elements with Qflag=TRUE,
c    the subtraction terms involve the (Gflag=TRUE) matrix elements
      if ((part .eq. 'real') .and. (Qflag .eqv. .true.)) then
        rGflag=Gflag
        Gflag=.true.
      endif     


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call spinoru(6,p,za,zb)

      prop=s(3,4)/dcmplx((s(3,4)-zmass**2),zmass*zwidth)

c--- calculate 2-quark, 2-gluon amplitudes
      if (Gflag) then
        call z2jetsq(1,2,3,4,5,6,za,zb,qqbZgg2)
        call storecsz(qqbZgg2_cs)
        call z2jetsq(1,5,3,4,2,6,za,zb,qgZqg2)
        call storecsz(qgZqg2_cs)
        call z2jetsq(2,5,3,4,1,6,za,zb,gqZqg2)
        call storecsz(gqZqg2_cs)
        call z2jetsq(1,6,3,4,2,5,za,zb,qgZgq2)
        call storecsz(qgZgq2_cs)
        call z2jetsq(2,6,3,4,1,5,za,zb,gqZgq2)
        call storecsz(gqZgq2_cs)
C --NB this is the matrix element for gg->Z qb(5) q(6)
        call z2jetsq(5,6,3,4,1,2,za,zb,ggZqbq2)
        call storecsz(ggZqbq2_cs)        

        do j=1,2
        do k=1,2
c        qbqZgg2(j,k)=qqbZgg2(swap(j),k)
c        qbgZqbg2(j,k)=qgZqg2(swap(j),k)
c        gqbZqbg2(j,k)=gqZqg2(swap(j),k)
c        qbgZgqb2(j,k)=qgZgq2(swap(j),k)
c        gqbZgqb2(j,k)=gqZgq2(swap(j),k)        
c        ggZqqb2(j,k)=ggZqbq2(swap(j),k)        
        do i=0,2
        qbqZgg2_cs(i,j,k)=qqbZgg2_cs(swap1(i),swap(j),k)
        qbgZqbg2_cs(i,j,k)=qgZqg2_cs(swap1(i),swap(j),k)
        gqbZqbg2_cs(i,j,k)=gqZqg2_cs(swap1(i),swap(j),k)
        qbgZgqb2_cs(i,j,k)=qgZgq2_cs(swap1(i),swap(j),k)
        gqbZgqb2_cs(i,j,k)=gqZgq2_cs(swap1(i),swap(j),k)        
        ggZqqb2_cs(i,j,k)=ggZqbq2_cs(swap1(i),swap(j),k)        
        enddo
        enddo
        enddo

c        call z2jetsq(2,1,3,4,5,6,za,zb,qbqZgg2)
c        call storecsz(qbqZgg2_cs)
c        call z2jetsq(5,1,3,4,2,6,za,zb,qbgZqbg2)
c        call storecsz(qbgZqbg2_cs)  
c        call z2jetsq(5,2,3,4,1,6,za,zb,gqbZqbg2)
c        call storecsz(gqbZqbg2_cs)
c        call z2jetsq(6,1,3,4,2,5,za,zb,qbgZgqb2)
c        call storecsz(qbgZgqb2_cs)  
c        call z2jetsq(6,2,3,4,1,5,za,zb,gqbZgqb2)
c        call storecsz(gqbZgqb2_cs)
C --NB this is the matrix element for gg->Z q(5) qb(6)
c        call z2jetsq(6,5,3,4,1,2,za,zb,ggZqqb2)
c        call storecsz(ggZqqb2_cs)        

        fac=v*xn/four*(esq*gsq)**2
        do pq=1,2
        do pl=1,2
        do i=0,2        
          qqbZgg2_cs(i,pq,pl) = half*aveqq*fac*qqbZgg2_cs(i,pq,pl)
          qbqZgg2_cs(i,pq,pl) = half*aveqq*fac*qbqZgg2_cs(i,pq,pl)
          gqZqg2_cs(i,pq,pl)  = aveqg*fac*gqZqg2_cs(i,pq,pl)
          qgZqg2_cs(i,pq,pl)  = aveqg*fac*qgZqg2_cs(i,pq,pl)
          gqbZqbg2_cs(i,pq,pl)= aveqg*fac*gqbZqbg2_cs(i,pq,pl)
          qbgZqbg2_cs(i,pq,pl)= aveqg*fac*qbgZqbg2_cs(i,pq,pl)

          gqZgq2_cs(i,pq,pl)  = aveqg*fac*gqZgq2_cs(i,pq,pl)
          qgZgq2_cs(i,pq,pl)  = aveqg*fac*qgZgq2_cs(i,pq,pl)
          gqbZgqb2_cs(i,pq,pl)= aveqg*fac*gqbZgqb2_cs(i,pq,pl)
          qbgZgqb2_cs(i,pq,pl)= aveqg*fac*qbgZgqb2_cs(i,pq,pl)

          ggZqbq2_cs(i,pq,pl) = avegg*fac*ggZqbq2_cs(i,pq,pl) 
          ggZqqb2_cs(i,pq,pl) = avegg*fac*ggZqqb2_cs(i,pq,pl) 
       enddo

        qqbZgg2(pq,pl) = qqbZgg2_cs(1,pq,pl)+qqbZgg2_cs(2,pq,pl)
     .                  +qqbZgg2_cs(0,pq,pl) 
        gqZqg2(pq,pl)  = gqZqg2_cs(1,pq,pl) +gqZqg2_cs(2,pq,pl)
     .                  +gqZqg2_cs(0,pq,pl)  
        qgZqg2(pq,pl)  = qgZqg2_cs(1,pq,pl)  +qgZqg2_cs(2,pq,pl)
     .                  +qgZqg2_cs(0,pq,pl)  
        gqZgq2(pq,pl)  = gqZgq2_cs(1,pq,pl) +gqZgq2_cs(2,pq,pl)
     .                  +gqZgq2_cs(0,pq,pl)  
        qgZgq2(pq,pl)  = qgZgq2_cs(1,pq,pl)  +qgZgq2_cs(2,pq,pl)
     .                  +qgZgq2_cs(0,pq,pl)  
        ggZqbq2(pq,pl) = ggZqbq2_cs(1,pq,pl) +ggZqbq2_cs(2,pq,pl)
     .                  +ggZqbq2_cs(0,pq,pl) 
c        gqbZqbg2(pq,pl)= gqbZqbg2_cs(1,pq,pl)+gqbZqbg2_cs(2,pq,pl)
c     .                  +gqbZqbg2_cs(0,pq,pl)
c        qbqZgg2(pq,pl) = qbqZgg2_cs(1,pq,pl)+qbqZgg2_cs(2,pq,pl)
c     .                  +qbqZgg2_cs(0,pq,pl) 
c        qbgZqbg2(pq,pl)= qbgZqbg2_cs(1,pq,pl)+qbgZqbg2_cs(2,pq,pl)
c     .                  +qbgZqbg2_cs(0,pq,pl)
c        gqbZgqb2(pq,pl)= gqbZgqb2_cs(1,pq,pl)+gqbZgqb2_cs(2,pq,pl)
c     .                  +gqbZgqb2_cs(0,pq,pl)
c        qbgZgqb2(pq,pl)= qbgZgqb2_cs(1,pq,pl)+qbgZgqb2_cs(2,pq,pl)
c     .                  +qbgZgqb2_cs(0,pq,pl)
c        ggZqqb2(pq,pl) = ggZqqb2_cs(1,pq,pl) +ggZqqb2_cs(2,pq,pl)
c     .                  +ggZqqb2_cs(0,pq,pl) 
        enddo
        enddo
      endif


      if (Qflag) then
      call spinoru(6,p,za,zb)

c--- qRb->qRb
      call ampqqb_qqb(1,5,2,6,qRb_a,qRb_b)
c--- qR->qR
c instead of calling ampqqb_qqb(1,5,6,2,qR_a,qR_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qR_a(j1,j2,j3)=+qRb_a(j1,swap(j2),j3)
      qR_b(j1,j2,j3)=-qRb_b(j1,swap(j2),j3)
      enddo
      enddo
      enddo
c--- qbR->qbR
      call ampqqb_qqb(6,1,5,2,qbR_a,qbR_b)

c--- qbRb->qbRb
c instead of calling ampqqb_qqb(5,1,2,6,qbRb_a,qbRb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbRb_a(j1,j2,j3)=-qRb_a(swap(j1),j2,j3)
      qbRb_b(j1,j2,j3)=+qRb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qqb->qqb
      call ampqqb_qqb(1,2,5,6,qqb_a,qqb_b)
c--- qbq->qqb
c instead of calling ampqqb_qqb(2,1,5,6,qbq_a,qbq_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_a(j1,j2,j3)=-qqb_a(swap(j1),j2,j3)
      qbq_b(j1,j2,j3)=+qqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo

c--- qq->qq
      call ampqqb_qqb(1,6,5,2,qq_a,qq_b)
c--- qbqb->qbqb
c instead of calling ampqqb_qqb(6,1,2,5,qbqb_a,qbqb_b)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbqb_a(j1,j2,j3)=-qq_a(swap(j1),swap(j2),j3)
      qbqb_b(j1,j2,j3)=-qq_b(swap(j1),swap(j2),j3)
      enddo
      enddo
      enddo

c--- four extra amplitudes added for 4-label matrix elements

c--- qRb->Rbq
c instead of calling ampqqb_qqb(1,6,2,5,qRb_ax,qRb_bx)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qRb_ax(j1,j2,j3)=-qbqb_a(swap(j1),j2,j3)
      qRb_bx(j1,j2,j3)=+qbqb_b(swap(j1),j2,j3)
      enddo
      enddo
      enddo
      
c--- qqb->qbq
c instead of calling ampqqb_qqb(1,2,6,5,qqb_ax,qqb_bx)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qqb_ax(j1,j2,j3)=+qqb_a(j1,swap(j2),j3)
      qqb_bx(j1,j2,j3)=-qqb_b(j1,swap(j2),j3)
      enddo
      enddo
      enddo

c--- qbR->qbR    
c instead of calling ampqqb_qqb(5,1,6,2,qbR_ax,qbR_bx)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbR_ax(j1,j2,j3)=-qRb_a(swap(j1),swap(j2),j3)
      qbR_bx(j1,j2,j3)=-qRb_b(swap(j1),swap(j2),j3)
      enddo
      enddo
      enddo

c--- qbq->qbq    
c instead of calling ampqqb_qqb(2,1,6,5,qbq_ax,qbq_bx)
      do j1=1,2
      do j2=1,2
      do j3=1,2
      qbq_ax(j1,j2,j3)=+qbq_a(j1,swap(j2),j3)
      qbq_bx(j1,j2,j3)=-qbq_b(j1,swap(j2),j3)
      enddo
      enddo
      enddo

      faclo=4d0*V*gsq**2*esq**2*aveqq 
      endif

      if (Gflag) then
      do j=-nf,nf
      do k=-nf,nf
      
      do icol=0,2
        msqx_cs(icol,j,k)=0d0
      enddo
      
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j .eq. 0) .and. (k .eq. 0)) then

           do icol=0,2
           ggtemp(icol)=0d0
           do nquark=1,nf
           msqx(icol,j,k,-nquark,nquark)=
     .      +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqbq2_cs(icol,1,1)
     .      +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqbq2_cs(icol,2,2)
     .      +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqbq2_cs(icol,1,2)
     .      +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqbq2_cs(icol,2,1)
           msqx(icol,j,k,nquark,-nquark)=
     .      +abs(Q(nquark)*q1+L(nquark)*l1*prop)**2*ggZqqb2_cs(icol,1,1)
     .      +abs(Q(nquark)*q1+R(nquark)*r1*prop)**2*ggZqqb2_cs(icol,2,2)
     .      +abs(Q(nquark)*q1+L(nquark)*r1*prop)**2*ggZqqb2_cs(icol,1,2)
     .      +abs(Q(nquark)*q1+R(nquark)*l1*prop)**2*ggZqqb2_cs(icol,2,1)
           ggtemp(icol)=ggtemp(icol)+msqx(icol,j,k,-nquark,nquark)
           enddo
           msqx_cs(icol,j,k)=ggtemp(icol)

           enddo
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZgg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZgg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZgg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZgg2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo

c---Statistical factor already included above
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZgg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZgg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZgg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZgg2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          do icol=0,2
c--- Madloop check: full Z in Born (no photon)
c             msqx_cs(icol,j,k)=
c     .       +abs(Q(j)*q1*0+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
c     .       +abs(Q(j)*q1*0+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
c     .       +abs(Q(j)*q1*0+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
c     .       +abs(Q(j)*q1*0+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1)
c--- Madloop check: axial coupling of Z only
c             msqx_cs(icol,j,k)=
c     .       +abs(Q(j)*q1*0+0.25d0/sin2w**2*prop)**2*qgZqg2_cs(icol,1,1)
c     .       +abs(Q(j)*q1*0+0.25d0/sin2w**2*prop)**2*qgZqg2_cs(icol,2,2)
c     .       +abs(Q(j)*q1*0+0.25d0/sin2w**2*prop)**2*qgZqg2_cs(icol,1,2)
c     .       +abs(Q(j)*q1*0+0.25d0/sin2w**2*prop)**2*qgZqg2_cs(icol,2,1)

c--- normal case
             msqx_cs(icol,j,k)=
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZqg2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZqg2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZqg2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZqg2_cs(icol,2,1)
             msqx(icol,j,k,j,k)=msqx_cs(icol,j,k)
             msqx(icol,j,k,k,j)=
     .       +abs(Q(j)*q1+L(j)*l1*prop)**2*qgZgq2_cs(icol,1,1)
     .       +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZgq2_cs(icol,2,2)
     .       +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZgq2_cs(icol,1,2)
     .       +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZgq2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqbg2_cs(icol,1,1)
     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqbg2_cs(icol,2,2)
     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqbg2_cs(icol,1,2)
     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqbg2_cs(icol,2,1)
             msqx(icol,j,k,j,k)=msqx_cs(icol,j,k)
             msqx(icol,j,k,k,j)=
     .       +abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZgqb2_cs(icol,1,1)
     .       +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZgqb2_cs(icol,2,2)
     .       +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZgqb2_cs(icol,1,2)
     .       +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZgqb2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZqg2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZqg2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZqg2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZqg2_cs(icol,2,1)
             msqx(icol,j,k,k,j)=msqx_cs(icol,j,k)
             msqx(icol,j,k,j,k)=
     .       +abs(Q(k)*q1+L(k)*l1*prop)**2*gqZgq2_cs(icol,1,1)
     .       +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZgq2_cs(icol,2,2)
     .       +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZgq2_cs(icol,1,2)
     .       +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZgq2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          do icol=0,2
             msqx_cs(icol,j,k)=
     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqbg2_cs(icol,1,1)
     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqbg2_cs(icol,2,2)
     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqbg2_cs(icol,1,2)
     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqbg2_cs(icol,2,1)
             msqx(icol,j,k,k,j)=msqx_cs(icol,j,k)
             msqx(icol,j,k,j,k)=
     .       +abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZgqb2_cs(icol,1,1)
     .       +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZgqb2_cs(icol,2,2)
     .       +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZgqb2_cs(icol,1,2)
     .       +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZgqb2_cs(icol,2,1)
             msqx_cs(icol,j,k)=msqx_cs(icol,j,k)
          enddo
      endif
      msq(j,k)=msqx_cs(0,j,k)+msqx_cs(1,j,k)+msqx_cs(2,j,k)
   19 continue
      enddo
      enddo
      endif

      if (Qflag) then

        do j=-nf,nf
        do k=-nf,nf

        do i=0,2
        mqq(i,j,k)=0d0
        enddo

        if ((j .gt. 0) .and. (k .gt. 0)) then
c----QQ case
          call msq_qq(j,k,qR_a,qR_b,qq_a,qq_b,prop,msq0,msq1,msq2)
                      
          msqx(0,j,k,j,k)=faclo*msq0
          msqx(1,j,k,j,k)=faclo*msq1
          msqx(2,j,k,j,k)=faclo*msq2
                      
          mqq(0,j,k)=msqx(0,j,k,j,k)
          mqq(1,j,k)=msqx(1,j,k,j,k)
          mqq(2,j,k)=msqx(2,j,k,j,k)

          if (j .ne. k) then
            call msq_qq(j,k,qq_a,qq_b,qR_a,qR_b,prop,msq0,msq1,msq2)
            msqx(0,j,k,k,j)=faclo*msq0
            msqx(1,j,k,k,j)=faclo*msq1
            msqx(2,j,k,k,j)=faclo*msq2
         endif
            
          elseif ((j .lt. 0) .and. (k .lt. 0)) then
c----QbQb case
          call msq_qq(-j,-k,qbRb_a,qbRb_b,qbqb_a,qbqb_b,prop,
     .                msq0,msq1,msq2)
                      
          msqx(0,j,k,j,k)=faclo*msq0
          msqx(1,j,k,j,k)=faclo*msq1
          msqx(2,j,k,j,k)=faclo*msq2
            
          mqq(0,j,k)=msqx(0,j,k,j,k)
          mqq(1,j,k)=msqx(1,j,k,j,k)
          mqq(2,j,k)=msqx(2,j,k,j,k)

          if (j .ne. k) then
            call msq_qq(-j,-k,qbqb_a,qbqb_b,qbRb_a,qbRb_b,prop,
     .                  msq0,msq1,msq2)
                      
            msqx(0,j,k,k,j)=faclo*msq0
            msqx(1,j,k,k,j)=faclo*msq1
            msqx(2,j,k,k,j)=faclo*msq2
          endif
            
        elseif ((j .gt. 0) .and. (k .lt. 0)) then
C---Q-Qb case
          call msq_qqb(j,k,qRb_a,qRb_b,qqb_a,qqb_b,prop,
     .                msq0,msq1,msq2,tup,tdo)
                      
          msqx(0,j,k,j,k)=faclo*msq0
          msqx(1,j,k,j,k)=faclo*msq1
          msqx(2,j,k,j,k)=faclo*msq2
            
          mqq(0,j,k)=msqx(0,j,k,j,k)
          mqq(1,j,k)=msqx(1,j,k,j,k)
          mqq(2,j,k)=msqx(2,j,k,j,k)
          
          if (j .eq. -k) then
            do f=1,5,2
              if (f .ne. j) then
                msqx(0,j,k,f,-f)=zip            
                msqx(1,j,k,f,-f)=faclo*tdo            
                msqx(2,j,k,f,-f)=zip 
              endif           
            enddo
            do f=2,4,2
              if (f .ne. j) then
                msqx(0,j,k,f,-f)=zip            
                msqx(1,j,k,f,-f)=faclo*tup            
                msqx(2,j,k,f,-f)=zip
              endif           
           enddo                              
            if ((j.eq.1).or.(j.eq.3).or.(j.eq.5)) then
              nup=2
              ndo=nf-3
            else
              nup=1
              ndo=nf-2
            endif
           mqq(1,j,k)=mqq(1,j,k)+faclo*(dfloat(nup)*tup+dfloat(ndo)*tdo)
         endif

         call msq_qqb(j,k,qRb_ax,qRb_bx,qqb_ax,qqb_bx,prop,
     .             msq0,msq1,msq2,tup,tdo)
                   
         msqx(0,j,k,k,j)=faclo*msq0
         msqx(1,j,k,k,j)=faclo*msq1
         msqx(2,j,k,k,j)=faclo*msq2
         
          if (j .eq. -k) then
            do f=1,5,2
              if (f .ne. j) then
                msqx(0,j,k,-f,f)=zip            
                msqx(1,j,k,-f,f)=faclo*tdo            
                msqx(2,j,k,-f,f)=zip
              endif            
            enddo
            do f=2,4,2
              if (f .ne. j) then
                msqx(0,j,k,-f,f)=zip            
                msqx(1,j,k,-f,f)=faclo*tup            
                msqx(2,j,k,-f,f)=zip
              endif
            enddo                              
         endif
                           
        elseif ((j .lt. 0) .and. (k .gt. 0)) then
C---Qb-Q case
          call msq_qqb(-j,-k,qbR_a,qbR_b,qbq_a,qbq_b,prop,
     .                msq0,msq1,msq2,tup,tdo)
                      
          msqx(0,j,k,k,j)=faclo*msq0
          msqx(1,j,k,k,j)=faclo*msq1
          msqx(2,j,k,k,j)=faclo*msq2
            
          mqq(0,j,k)=msqx(0,j,k,k,j)
          mqq(1,j,k)=msqx(1,j,k,k,j)
          mqq(2,j,k)=msqx(2,j,k,k,j)
          
          if (j .eq. -k) then
            do f=1,5,2
              if (f .ne. k) then
                msqx(0,j,k,f,-f)=zip            
                msqx(1,j,k,f,-f)=faclo*tdo            
                msqx(2,j,k,f,-f)=zip 
              endif           
            enddo
            do f=2,4,2
              if (f .ne. k) then
                msqx(0,j,k,f,-f)=zip            
                msqx(1,j,k,f,-f)=faclo*tup            
                msqx(2,j,k,f,-f)=zip
              endif           
           enddo                              
            if ((k.eq.1).or.(k.eq.3).or.(k.eq.5)) then
              nup=2
              ndo=nf-3
            else
              nup=1
              ndo=nf-2
            endif
           mqq(1,j,k)=mqq(1,j,k)+faclo*(dfloat(nup)*tup+dfloat(ndo)*tdo)
         endif

         call msq_qqb(-j,-k,qbR_ax,qbR_bx,qbq_ax,qbq_bx,prop,
     .             msq0,msq1,msq2,tup,tdo)
                   
         msqx(0,j,k,j,k)=faclo*msq0
         msqx(1,j,k,j,k)=faclo*msq1
         msqx(2,j,k,j,k)=faclo*msq2
         
         
          if (j .eq. -k) then
            do f=1,5,2
              if (f .ne. k) then
                msqx(0,j,k,-f,f)=zip            
                msqx(1,j,k,-f,f)=faclo*tdo            
                msqx(2,j,k,-f,f)=zip
              endif            
            enddo
            do f=2,4,2
              if (f .ne. k) then
                msqx(0,j,k,-f,f)=zip            
                msqx(1,j,k,-f,f)=faclo*tup            
                msqx(2,j,k,-f,f)=zip
              endif
            enddo                              
         endif
                           
      endif

      msq(j,k)=msq(j,k)+mqq(0,j,k)+mqq(1,j,k)+mqq(2,j,k)

      enddo
      enddo

      endif

  999 continue
c--- restore proper colourchoice if necessary
      if ((part .eq. 'real') .or. (part .eq. 'virt')) then
        colourchoice=rcolourchoice
      endif
c--- restore proper parton sub-process selection, if necessary
      if ((part .eq. 'real') .and. (Qflag .eqv. .true.)) then
        Gflag=rGflag
      endif     

      return
      end
          
    

      
     
