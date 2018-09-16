      subroutine qqb_tbb(p,msq)
      implicit none
      include 'types.f'

c     Matrix element for t-bbar production
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->n(p3)+e^+(p4)+b(p5)+bbar(p6)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->e^-(p3)+n(p4)+bbar(p5)+b(p6)
C     averaged(summed) over initial(final) colours and spins
c--NB average over spins only -- colour factors cancel
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'sprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qqbtbbar
      real(dp):: fac,qqb,qbq
      integer:: j,k

      call dotem(6,p,s)
      fac=aveqq*gw**8*xn**2

      if     (nwz == +1) then
        qqb=fac*qqbtbbar(1,6,3,4,5,2)
        qbq=fac*qqbtbbar(2,6,3,4,5,1)
      elseif (nwz == -1) then
        qqb=fac*qqbtbbar(2,6,4,3,5,1)
        qbq=fac*qqbtbbar(1,6,4,3,5,2)
      endif

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=Vsq(j,k)*qqb
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=Vsq(j,k)*qbq
      endif
      enddo
      enddo

      return
      end

      function qqbtbbar(ju,jb,jn,je,jc,jd)
      implicit none
      include 'types.f'
      real(dp):: qqbtbbar

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
c      u(ju)+b(jb)-->n(jn)+e^+(je)+c(jc)+d(jd)
C     This form is valid also if the mass of the b-quark is non-zero
C     or if the top is not on shell
      include 'sprods_com.f'
      include 'masses.f'
      integer:: ju,jb,jn,je,jc,jd
      real(dp):: st,prop

      st=s(jn,je)+s(je,jc)+s(jn,jc)
      if (s(ju,jd) < 0._dp) then
      prop=(s(ju,jd)-wmass**2)**2
      else
      prop=(s(ju,jd)-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=prop
     &    *((s(jn,je)-wmass**2)**2+(wmass*wwidth)**2)
     &    *((st-mt**2)**2+(mt*twidth)**2)
      qqbtbbar=s(jc,jn)*s(ju,jb)
     & *(-(s(ju,jd)+s(jb,jd))*(s(jn,je)+s(jc,je))-st*s(je,jd))/prop
      return
      end
