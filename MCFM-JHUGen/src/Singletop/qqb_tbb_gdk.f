      subroutine qqb_tbb_gdk(p,msq)
      implicit none
      include 'types.f'

c     Matrix element for real corrections to single top production
C     (nwz=+1)
c      u(-p1)+dbar(-p2)-->t(=> n(p3)+e^+(p4)+b(p5)+g(p7))+F(p6)
C     or for
C     (nwz=-1)
c      ubar(-p1)+d(-p2)-->t~(=> e^-(p3)+n(p4)+bbar(p5)+g(p7))+F(p6)
C     averaged(summed) over initial(final) colours and spins

c--- g(p7) represents a gluon


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'nwz.f'
      include 'zprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qqbtbbargd
      real(dp):: fac,qqb,qbq
      integer:: j,k

      call spinoru(7,p,za,zb)
      fac=2._dp*gsq*cf*gw**8*xn**2


      if (nwz == +1) then
c--- t production
        qqb=aveqq*fac*qqbtbbargd(1,6,3,4,5,2,7,p)
        qbq=aveqq*fac*qqbtbbargd(2,6,3,4,5,1,7,p)
      elseif (nwz == -1) then
c--- tbar production
        qqb=aveqq*fac*qqbtbbargd(2,6,4,3,5,1,7,p)
        qbq=aveqq*fac*qqbtbbargd(1,6,4,3,5,2,7,p)
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


      function qqbtbbargd(ju,jb,jn,je,jc,jd,jg,p)
      implicit none
      include 'types.f'
      real(dp):: qqbtbbargd
C     Matrix element squared for single top production with gluon
C     radiation in decay (radiation from final line)
C      u(ju) b(jb) -> t(n~(jn)+e+(je)+c(jc)+g(jg))+d(jd)
C     masses of b quarks c.c=b.b=0

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      integer:: ju,jb,jn,je,jc,jd,jg,nu
      real(dp):: p(mxpart,4),pt(4),ptDpt
      real(dp):: sne,sdu,prop,twoptg
      complex(dp):: ampf(2),amp0

      do nu=1,4
      pt(nu)=p(je,nu)+p(jn,nu)+p(jc,nu)+p(jg,nu)
      enddo
      ptDpt=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2
      twoptg=
     & 2._dp*(pt(4)*p(jg,4)-pt(1)*p(jg,1)-pt(2)*p(jg,2)-pt(3)*p(jg,3))

      sne=s(jn,je)
      sdu=s(jd,ju)
      if (sdu < 0._dp) then
      prop=(sdu-wmass**2)**2
      else
      prop=(sdu-wmass**2)**2+(wmass*wwidth)**2
      endif
      prop=((sne-wmass**2)**2+(wmass*wwidth)**2)
     &    *((ptDpt-mt**2)**2+(mt*twidth)**2)*prop

C  -Lefthanded gluon tb-line
c      ampf(1)=(zb(je,jc)*za(jc,jg)+zb(je,jn)*za(jn,jg))
c     & *(+zb(jc,je)*za(je,jd)+zb(jc,jn)*za(jn,jd)+zb(jc,jg)*za(jg,jd))
c     & +mt**2*zb(je,jc)*za(jg,jd)
c      ampf(1)=-za(jc,jn)*zb(ju,jb)/(twoptg*zb(jg,jc))*ampf(1)
c      ampf(1)=ampf(1)-za(jg,jn)*zb(ju,jb)/zb(jg,jc)
c     & *(+zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))
c      write(6,*) 'ampf(1)',ampf(1)


      amp0=(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))
     &  *za(jc,jn)*zb(ju,jb)
c---eikonal form
      ampf(1)=-amp0
     & *(za(jg,je)*zb(je,jc)+za(jg,jn)*zb(jn,jc))/zb(jg,jc)/twoptg
      ampf(1)=ampf(1)-za(jg,jn)*zb(ju,jb)/zb(jg,jc)
     & *(zb(je,jc)*za(jc,jd)+zb(je,jn)*za(jn,jd)+zb(je,jg)*za(jg,jd))

C  -Righthanded gluon tb-line
c      ampf(2)=-(zb(je,jn)*za(jn,jc)
c     & *(+zb(jg,jc)*za(jc,jd)+zb(jg,je)*za(je,jd)+zb(jg,jn)*za(jn,jd))
c     & +mt**2*zb(je,jg)*za(jc,jd))
c      ampf(2)=za(jc,jn)*zb(ju,jb)/za(jc,jg)/twoptg*ampf(2)
c      write(6,*) 'ampf(2)',ampf(2)

c---eikonal form
      ampf(2)=-amp0
     & *(zb(jg,je)*za(je,jc)+zb(jg,jn)*za(jn,jc))/za(jc,jg)/twoptg
     & -za(jc,jn)*zb(ju,jb)/twoptg*zb(je,jg)
     & *(zb(jg,jc)*za(jc,jd)+zb(jg,je)*za(je,jd)+zb(jg,jn)*za(jn,jd))


      qqbtbbargd=(abs(ampf(1))**2+abs(ampf(2))**2)/prop

      return
      end
