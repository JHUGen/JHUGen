      subroutine bq_tpq_gdk(p,msq)
      implicit none
c     Matrix element for t-bbar production
c      b(-p1)+u(-p2)-->t(n(p3)+e^+(p4)+b(p5)+g(p7))+d(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     . qqbtbbargd
      double precision fac,ub,bu,bubar,ubarb
      integer j,k,ib

      call spinoru(7,p,za,zb)
      fac=2d0*gsq*cf*aveqq*gw**8*xn**2
      if     (nwz .eq. +1) then
        ub=fac*qqbtbbargd(1,2,3,4,5,6,7,p)
        bu=fac*qqbtbbargd(2,1,3,4,5,6,7,p)
        ubarb=fac*qqbtbbargd(6,2,3,4,5,1,7,p)
        bubar=fac*qqbtbbargd(6,1,3,4,5,2,7,p)
      elseif (nwz .eq. -1) then
        ub=fac*qqbtbbargd(6,2,4,3,5,1,7,p)
        bu=fac*qqbtbbargd(6,1,4,3,5,2,7,p)
        ubarb=fac*qqbtbbargd(1,2,4,3,5,6,7,p)
        bubar=fac*qqbtbbargd(2,1,4,3,5,6,7,p)
      endif

c--- for nwz=+1, initial state is b, for nwz=-1 it is b~
      ib=5*nwz
      
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      if     ((j .eq. ib) .and. (k .gt. 0)) then
      msq(j,k)=(Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*bu
      elseif ((j .eq. ib) .and. (k .lt. 0)) then
      msq(j,k)=(Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*bubar
      elseif ((j .gt. 0) .and. (k .eq. ib)) then
      msq(j,k)=(Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*ub
      elseif ((j .lt. 0) .and. (k .eq. ib)) then
      msq(j,k)=(Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*ubarb
      endif
      enddo
      enddo

      return
      end

