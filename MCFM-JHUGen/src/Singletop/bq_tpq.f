      subroutine bq_tpq(p,msq)
      implicit none
      include 'types.f'

c     Matrix element for t-bbar production
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
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qqbtbbar
      real(dp):: fac,bu,ub,bubar,ubarb
      integer:: j,k,ib

      call spinoru(6,p,za,zb)

      fac=aveqq*gw**8*xn**2

      if     (nwz == +1) then
        ub=fac*qqbtbbar(1,2,3,4,5,6)
        bu=fac*qqbtbbar(2,1,3,4,5,6)
        ubarb=fac*qqbtbbar(6,2,3,4,5,1)
        bubar=fac*qqbtbbar(6,1,3,4,5,2)
      elseif (nwz == -1) then
        ub=fac*qqbtbbar(6,2,4,3,5,1)
        bu=fac*qqbtbbar(6,1,4,3,5,2)
        ubarb=fac*qqbtbbar(1,2,4,3,5,6)
        bubar=fac*qqbtbbar(2,1,4,3,5,6)
      endif

c--- for nwz=+1, initial state is b, for nwz=-1 it is b~
      ib=5*nwz

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if     ((j == ib) .and. (k > 0)) then
      msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*bu
      elseif ((j == ib) .and. (k < 0)) then
      msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*bubar
      elseif ((j > 0) .and. (k == ib)) then
      msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*ub
      elseif ((j < 0) .and. (k == ib)) then
      msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*ubarb
      endif

      enddo
      enddo

      return
      end

