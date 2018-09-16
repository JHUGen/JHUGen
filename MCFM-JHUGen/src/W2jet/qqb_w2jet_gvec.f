      subroutine qqb_w2jet_gvec(p,n,in,msq)
      implicit none
      include 'types.f'

c----Matrix element for W+2jet production
C----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
C For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p5)+ g(p6)+W^+(n(p3)+e^+(p4))
C For nwz=-1
c     d(-p1)+ubar(-p2)--> g(p5)+ g(p6)+W^-(e^-(p3)+nbar(p4))
c---It has been checked that this gives the right matrix element
c---squared when n is replaced by two physical polarizations,
c---but this routine contains ne factor of 1/2 for identical gluons
c   in the final state.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ckm.f'
      include 'msqv_cs.f'
C ip is the label of the emitter
      integer:: j,k,in,i,n1,n2
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),tmp
      real(dp):: fac,prop,p1p2(0:2,-1:1,-1:1),n(4),Vfac
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/p1p2/p1p2
!$omp threadprivate(/p1p2/)
c--- note that we will use the first index of p1p2 to label
c--- the colour structure of the squared matrix element:
c---     0 --> -ninth*(qed piece)
c---     1 --> +(qcd ordering 1)
c---     2 --> +(qcd ordering 1)
c---     3 --> Total (0+1+2)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      do i=0,2
        msqv_cs(i,j,k)=0._dp
      enddo
      enddo
      enddo

      do i=0,2
      do j=-1,1
      do k=-1,1
      p1p2(i,j,k)=0._dp
      enddo
      enddo
      enddo

      prop=s(3,4)**2/((s(3,4)-wmass**2)**2+wmass**2*wwidth**2)
      fac=V*xn/eight*(gsq*gwsq)**2*prop
      call spinoru(6,p,za,zb)
C   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Arbitrary conventions of Bern, Dixon, Kosower, Weinzierl,
c---za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)

      if (in == 1) then
C--initial-initial
        call w2jetn(2,5,3,4,6,1,p,n,za,zb,zab,zba)
        call storecsv(0,+1)
        call w2jetn(5,2,3,4,6,1,p,n,za,zb,zab,zba)
        call storecsv(0,-1)
        call w2jetn(5,6,3,4,2,1,p,n,za,zb,zab,zba)
        call storecsv(0,0)
        do i=0,2
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
          p1p2(i,0, 0)=avegg*fac*p1p2(i,0,0)
        enddo
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,6,3,4,5,1,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(6,2,3,4,5,1,p,n,za,zb,zab,zba)
      elseif (in == 2) then
        call w2jetn(1,5,3,4,6,2,p,n,za,zb,zab,zba)
        call storecsv(+1,0)
        call w2jetn(5,1,3,4,6,2,p,n,za,zb,zab,zba)
        call storecsv(-1,0)
        call w2jetn(5,6,3,4,1,2,p,n,za,zb,zab,zba)
        call storecsv(0,0)
c--- since we have interchanged 1 and 2 to get the gg matrix element,
c---  the colour structures should be interchanged too
        tmp=p1p2(1,0,0)
        p1p2(1,0,0)=p1p2(2,0,0)
        p1p2(2,0,0)=tmp
        do i=0,2
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0, 0)=avegg*fac*p1p2(i,0,0)
        enddo
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,6,3,4,5,2,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(6,1,3,4,5,2,p,n,za,zb,zab,zba)
      elseif (in == 5) then
        call w2jetn(1,2,3,4,6,5,p,n,za,zb,zab,zba)
        call storecsv(1,-1)
        call w2jetn(2,1,3,4,6,5,p,n,za,zb,zab,zba)
        call storecsv(-1,1)
        call w2jetn(1,6,3,4,2,5,p,n,za,zb,zab,zba)
        call storecsv(+1,0)
        call w2jetn(6,1,3,4,2,5,p,n,za,zb,zab,zba)
        call storecsv(-1,0)
        call w2jetn(2,6,3,4,1,5,p,n,za,zb,zab,zba)
        call storecsv(0,+1)
        call w2jetn(6,2,3,4,1,5,p,n,za,zb,zab,zba)
        call storecsv(0,-1)
        do i=0,2
          p1p2(i,1,-1)=half*aveqq*fac*p1p2(i,1,-1)
          p1p2(i,-1,1)=half*aveqq*fac*p1p2(i,-1,1)
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
        enddo
c      p1p2(1,-1)=aveqq*fac*w2jetn(1,2,3,4,6,5,p,n,za,zb,zab,zba)
c      p1p2(-1,1)=aveqq*fac*w2jetn(2,1,3,4,6,5,p,n,za,zb,zab,zba)
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,6,3,4,2,5,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(6,1,3,4,2,5,p,n,za,zb,zab,zba)
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,6,3,4,1,5,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(6,2,3,4,1,5,p,n,za,zb,zab,zba)
      elseif (in == 6) then
        call w2jetn(1,2,3,4,5,6,p,n,za,zb,zab,zba)
        call storecsv(1,-1)
        call w2jetn(2,1,3,4,5,6,p,n,za,zb,zab,zba)
        call storecsv(-1,1)
        call w2jetn(1,5,3,4,2,6,p,n,za,zb,zab,zba)
        call storecsv(+1,0)
        call w2jetn(5,1,3,4,2,6,p,n,za,zb,zab,zba)
        call storecsv(-1,0)
        call w2jetn(2,5,3,4,1,6,p,n,za,zb,zab,zba)
        call storecsv(0,+1)
        call w2jetn(5,2,3,4,1,6,p,n,za,zb,zab,zba)
        call storecsv(0,-1)
        do i=0,2
          p1p2(i,1,-1)=half*aveqq*fac*p1p2(i,1,-1)
          p1p2(i,-1,1)=half*aveqq*fac*p1p2(i,-1,1)
          p1p2(i,+1,0)=aveqg*fac*p1p2(i,+1,0)
          p1p2(i,-1,0)=aveqg*fac*p1p2(i,-1,0)
          p1p2(i,0,+1)=aveqg*fac*p1p2(i,0,+1)
          p1p2(i,0,-1)=aveqg*fac*p1p2(i,0,-1)
        enddo
c      p1p2(1,-1)=aveqq*fac*w2jetn(1,2,3,4,5,6,p,n,za,zb,zab,zba)
c      p1p2(-1,1)=aveqq*fac*w2jetn(2,1,3,4,5,6,p,n,za,zb,zab,zba)
c      p1p2(+1,0)=aveqg*fac*w2jetn(1,5,3,4,2,6,p,n,za,zb,zab,zba)
c      p1p2(-1,0)=aveqg*fac*w2jetn(5,1,3,4,2,6,p,n,za,zb,zab,zba)
c      p1p2(0,+1)=aveqg*fac*w2jetn(2,5,3,4,1,6,p,n,za,zb,zab,zba)
c      p1p2(0,-1)=aveqg*fac*w2jetn(5,2,3,4,1,6,p,n,za,zb,zab,zba)
      endif

      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
          do i=0,2
            msqv_cs(i,j,k)=Vsq(j,k)*p1p2(i,1,-1)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          do i=0,2
            msqv_cs(i,j,k)=Vsq(j,k)*p1p2(i,-1,1)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*p1p2(i,+1,0)
          enddo
      elseif ((j < 0) .and. (k == 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*p1p2(i,-1,0)
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*p1p2(i,0,+1)
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do i=0,2
            msqv_cs(i,j,k)=
     &  (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*p1p2(i,0,-1)
          enddo
      elseif ((j == 0) .and. (k == 0)) then
          Vfac=0._dp
          do n1=1,nf
            do n2=-nf,-1
              Vfac=Vfac+Vsq(n1,n2)
            enddo
          enddo
          do i=0,2
            msqv_cs(i,j,k)=Vfac*p1p2(i,0,0)
          enddo
      endif
      msq(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
      enddo
      enddo

      return
      end

      subroutine w2jetn(i1,i2,i3,i4,i5,i6,p,n,za,zb,zab,zba)
      implicit none
      include 'types.f'
C----matrix element squared with p5 line contracted with n(mu)
C----nDp6 should be equal to zero

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'mmsqv_cs.f'
      complex(dp):: qcdabn(2,2,2),qcdban(2,2,2),qedn(2,2,2)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msq1,msq2,msqq,n(4),p(mxpart,4)
      real(dp):: nDp5
      integer:: i1,i2,i3,i4,i5,i6

      nDp5=n(4)*p(i5,4)-n(3)*p(i5,3)-n(2)*p(i5,2)-n(1)*p(i5,1)

      call checkndotp(p,n,i6)

      call subqcdn(i1,i2,i3,i4,i5,i6,nDp5,za,zb,zab,zba,qcdabn,qcdban)

C--first argument is quark line
C--second argument is polarization of i5 line
C  1=L,2=R
      qedn(1,1,1)=qcdabn(1,1,1)+qcdban(1,1,1)
      qedn(2,1,1)=qcdabn(2,1,1)+qcdban(2,1,1)

      msq1= abs(qcdabn(1,1,1))**2+abs(qcdabn(2,1,1))**2
      msq2= abs(qcdban(1,1,1))**2+abs(qcdban(2,1,1))**2
      msqq= abs(qedn(1,1,1))**2+abs(qedn(2,1,1))**2

      mmsqv_cs(0,+1,+1)=-ninth*msqq
      mmsqv_cs(1,+1,+1)=msq1
      mmsqv_cs(2,+1,+1)=msq2

      return
      end

      subroutine storecsv(i,j)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet_gvec matrix elements into elements (..,i,j) of p1p2

      include 'mmsqv_cs.f'
      integer:: i,j,k
      real(dp):: p1p2(0:2,-1:1,-1:1)
      common/p1p2/p1p2
!$omp threadprivate(/p1p2/)

      do k=0,2
        p1p2(k,i,j)=mmsqv_cs(k,+1,+1)
      enddo

      return
      end



