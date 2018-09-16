      subroutine qqb_tbb_vdk(p,msqv)
      implicit none
      include 'types.f'

c     Virtual matrix element for t-bbar decay
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->n(p3)+e^+(p4)+b(p5)+bbar(p6)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->e^-(p3)+n(p4)+bbar(p5)+b(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      integer:: j,k
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'scheme.f'
      real(dp):: p(mxpart,4),msqv(-nf:nf,-nf:nf),
     & virtqqbdk,qqb,qbq,fac

      scheme='dred'

      call spinoru(6,p,za,zb)
      fac=ason2pi*cf
      fac=aveqq*xn**2*gw**8*fac

      if     (nwz == +1) then
        qqb=fac*virtqqbdk(1,6,3,4,5,2)
        qbq=fac*virtqqbdk(2,6,3,4,5,1)
      elseif (nwz == -1) then
        qqb=fac*virtqqbdk(2,6,4,3,5,1)
        qbq=fac*virtqqbdk(1,6,4,3,5,2)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      if     ((j > 0) .and. (k < 0)) then
      msqv(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
      msqv(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end


      function virtqqbdk(ju,jb,jn,je,jc,jd)
      implicit none
      include 'types.f'
      real(dp):: virtqqbdk
C----virtual correction for decay of top quark.
C----nlo correction to the width is included.


      integer:: ju,jd,jn,je,jc,jb
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      real(dp):: snec,prop,mtsq,cv,ct,c1
      complex(dp):: amp,ampho


      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)

      call coefsdk(s(jn,je),mtsq,ct,cv,c1)
      if (s(ju,jd) < 0._dp) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      prop=prop*((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)

      amp=za(jc,jn)*zb(ju,jb)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))

      ampho=za(jc,jn)*zb(ju,jb)
     & *(cplx1(ct+cv)*(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd))
     & +cplx1(0.5_dp*c1)*zb(je,jc)*za(jc,jd))

      virtqqbdk=real(amp*conjg(ampho))/prop
      return
      end


