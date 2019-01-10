      subroutine qqb_w_cjet_massless_g(p,msq)
      implicit none
c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + c(p5) + f(p6)
c                            |
c                            --> e^-(p3) + nubar(p4)
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'msq_cs.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     .                 facgg,prop,Vfac
      double precision qgWqg2,qbgWqbg2,gqbWqbg2,gqWqg2,ggWqqb2,ggWqbq2
      double precision qgWqg2_cs(0:2),qbgWqbg2_cs(0:2),
     .                 gqbWqbg2_cs(0:2),gqWqg2_cs(0:2),
     .                 ggWqqb2_cs(0:2),ggWqbq2_cs(0:2)
      double precision sbqWcbq,sqWcq,qsWcq,qsbWcbq,qqbWcsb,qqbWcbs,s34
      common/facgg/facgg
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer i,j,k,n1,n2
!$omp threadprivate(/facgg/)

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

c--- calculate 4-quark contributions
c--- note that these are symmetric under interchange of q and qb
      call wqq_sc(2,5,3,4,6,1,p,qsWcq)
      call wqq_sc(1,5,3,4,6,2,p,sqWcq)

      call wqq_sc(1,5,4,3,6,2,p,sbqWcbq)
      call wqq_sc(2,5,4,3,6,1,p,qsbWcbq)

      call wqq_sc(6,5,3,4,2,1,p,qqbWcsb)
      call wqq_sc(6,5,4,3,1,2,p,qqbWcbs)

      s34=2d0*(p(3,4)*p(4,4)-p(3,1)*p(4,1)-p(3,2)*p(4,2)-p(3,3)*p(4,3))
      prop=s34**2/((s34-wmass**2)**2+wmass**2*wwidth**2)
      facgg=V*xn/four*(gwsq/2d0*gsq)**2*prop
 
c--- calculate 2-quark, 2-gluon amplitudes
  
        call w2jetsq_mass(1,5,4,3,2,6,p,qbgWqbg2)
        call storecs(qbgWqbg2_cs)

        call w2jetsq_mass(2,5,4,3,1,6,p,gqbWqbg2)
        call storecs(gqbWqbg2_cs)

        call w2jetsq_mass(6,5,4,3,1,2,p,ggWqbq2)
        call storecs(ggWqbq2_cs)        

        call w2jetsq_mass(6,5,3,4,1,2,p,ggWqqb2)
        call storecs(ggWqqb2_cs)        

        call w2jetsq_mass(2,5,3,4,1,6,p,gqWqg2)
        call storecs(gqWqg2_cs)

        call w2jetsq_mass(1,5,3,4,2,6,p,qgWqg2)
        call storecs(qgWqg2_cs)


        do i=0,2        
          gqWqg2_cs(i)  = aveqg*facgg*gqWqg2_cs(i)
          qgWqg2_cs(i)  = aveqg*facgg*qgWqg2_cs(i)
          gqbWqbg2_cs(i)= aveqg*facgg*gqbWqbg2_cs(i)
          qbgWqbg2_cs(i)= aveqg*facgg*qbgWqbg2_cs(i)
          ggWqqb2_cs(i) = avegg*facgg*ggWqqb2_cs(i)
          ggWqbq2_cs(i) = avegg*facgg*ggWqbq2_cs(i)
        enddo
        gqWqg2  = gqWqg2_cs(1)  +gqWqg2_cs(2)  +gqWqg2_cs(0)  
        qgWqg2  = qgWqg2_cs(1)  +qgWqg2_cs(2)  +qgWqg2_cs(0)  
        gqbWqbg2= gqbWqbg2_cs(1)+gqbWqbg2_cs(2)+gqbWqbg2_cs(0)
        qbgWqbg2= qbgWqbg2_cs(1)+qbgWqbg2_cs(2)+qbgWqbg2_cs(0)
        ggWqqb2 = ggWqqb2_cs(1) +ggWqqb2_cs(2) +ggWqqb2_cs(0) 
        ggWqbq2 = ggWqbq2_cs(1) +ggWqbq2_cs(2) +ggWqbq2_cs(0) 
          

      
      do j=-(nf-2),(nf-2)
      do k=-(nf-2),(nf-2)
      
c--- 2-quark, 2-gluon contribution to matrix elements      
      if ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,-4)*qgWqg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,-4)*qgWqg2_cs(i)
          enddo
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=msq(j,k)+Vsq(j,+4)*qbgWqbg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(j,+4)*qbgWqbg2_cs(i)
          enddo
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=msq(j,k)+Vsq(-4,k)*gqWqg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(-4,k)*gqWqg2_cs(i)
          enddo
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=msq(j,k)+Vsq(+4,k)*gqbWqbg2
          do i=0,2
            msq_cs(i,j,k)=Vsq(+4,k)*gqbWqbg2_cs(i)
          enddo
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          Vfac=0d0
            do n1=1,nf
            Vfac=Vfac+Vsq(n1,-4)
            enddo
          msq(j,k)=msq(j,k)+Vfac*ggWqqb2
          do i=0,2
            msq_cs(i,j,k)=Vfac*ggWqqb2_cs(i)
          enddo
          Vfac=0d0
            do n2=-nf,-1
            Vfac=Vfac+Vsq(4,n2)
            enddo
          msq(j,k)=msq(j,k)+Vfac*ggWqbq2
          do i=0,2
            msq_cs(i,j,k)=msq_cs(i,j,k)+Vfac*ggWqbq2_cs(i)
          enddo
      endif

c--- 4-quark contribution to matrix elements      
      if     ((j .gt. 0) .and. (k .gt. 0)) then
        msq(j,k)=Vsq(j,-4)*sqWcq+Vsq(k,-4)*qsWcq
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
        msq(j,k)=Vsq(j,+4)*sbqWcbq+Vsq(k,+4)*qsbWcbq
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
        msq(j,k)=Vsq(j,-4)*sqWcq+Vsq(k,+4)*qsbWcbq
        if (j .eq. -k) then
          msq(j,k)=msq(j,k)+Vsum(-4)*qqbWcsb+Vsum(+4)*qqbWcbs
        endif
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
        msq(j,k)=Vsq(j,+4)*sbqWcbq+Vsq(k,-4)*qsWcq
        if (j .eq. -k) then
          msq(j,k)=msq(j,k)+Vsum(-4)*qqbWcsb+Vsum(+4)*qqbWcbs
        endif
      endif

      enddo
      enddo
      
      return
      end
     
