      subroutine bq_tpq_v(p,msqv)
      implicit none
      include 'types.f'

c     Virtual matrix element for t-bbar production
C     (nwz=+1)
c      b(-p1)+u(-p2)-->n(p3)+e^+(p4)+b(p5)+d(p6)
C     or for
C     (nwz=-1)
c      b~(-p1)+d(-p2)-->e^-(p3)+n~(p4)+b~(p5)+u(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'scheme.f'
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,bu,ub,bubar,ubarb,virtqqb
      integer:: j,k,ib

      scheme='dred'

      call spinoru(6,p,za,zb)

      fac=ason2pi*cf
      fac=aveqq*gw**8*xn**2*fac
C      virtqqb(ju,jb,jn,je,jc,jd)

      if     (nwz == +1) then
        ub=+fac*virtqqb(1,2,3,4,5,6)
        bu=+fac*virtqqb(2,1,3,4,5,6)
        bubar=+fac*virtqqb(6,1,3,4,5,2)
        ubarb=+fac*virtqqb(6,2,3,4,5,1)
      elseif (nwz == -1) then
        ub=fac*virtqqb(6,2,4,3,5,1)
        bu=fac*virtqqb(6,1,4,3,5,2)
        bubar=fac*virtqqb(2,1,4,3,5,6)
        ubarb=fac*virtqqb(1,2,4,3,5,6)
      endif

c--- for nwz=+1, initial state is b, for nwz=-1 it is b~
      ib=5*nwz

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp

      if     ((j == ib) .and. (k > 0)) then
      msqv(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))
     & *bu
      elseif ((j == ib) .and. (k < 0)) then
      msqv(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))
     & *bubar
      elseif ((j > 0) .and. (k == ib)) then
      msqv(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))
     & *ub
      elseif ((j < 0) .and. (k == ib)) then
      msqv(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))
     & *ubarb
      endif

      enddo
      enddo

      return
      end

