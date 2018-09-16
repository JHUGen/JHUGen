      subroutine WW_Hqq(p,msq)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion by W-W exchange only
c---Matrix element squared averaged over initial colors and spins
c
c     q(-p1)+q(-p2) -->  H(p3,p4)+q(p5)+q(p6)
c                           |
c                           |
c                           |
c                           ---> b(p3)+bbar(p4)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: p(mxpart,4),fac,s34
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_du,uub_ddb,msqgamgam

      integer,parameter::pn(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

      call dotem(6,p,s)

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      fac=0.25_dp*gwsq**3*hdecay
C Color cancels, 0.25_dp is spin average

C q-q and qbar-qbar
c--- u(1)+d(2) -> d(5)+u(6)
c--- ub(1)+db(2) -> db(5)+ub(6)
      call msqpieces_ww(1,2,6,5,ud_du)

C q-qbar and qbar-q
c--- u(1)+ub(2) -> d(5)+db(6)
c--- ub(1)+d(2) -> db(5)+u(6)
      call msqpieces_ww(1,6,2,5,uub_ddb)

c--- Only loop up to (nf-1) to avoid b->t transitions
      do j=-(nf-1),nf-1
      do k=-(nf-1),nf-1
      msq(j,k)=0._dp
        if     ((j > 0) .and. (k < 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j < 0) .and. (k > 0)) then
          if (pn(j) == -pn(k)) msq(j,k)=fac*uub_ddb
        elseif ((j > 0) .and. (k > 0)) then
          if (pn(j)+pn(k) == +3) msq(j,k)=fac*ud_du
        elseif ((j < 0) .and. (k < 0)) then
          if (pn(j)+pn(k) == -3) msq(j,k)=fac*ud_du
        endif
      enddo
      enddo

      return
      end


      subroutine msqpieces_ww(i1,i2,i5,i6,wll)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      real(dp):: wll,htheta
      real(dp):: propw,x
      integer:: i1,i2,i5,i6
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)
      propw(i1,i2)=sign(one,(s(i1,i2)-wmass**2))*sqrt(
     .((s(i1,i2)-wmass**2)**2+htheta(s(i1,i2))*(wmass*wwidth)**2)/wmass)
      wll=s(i1,i2)*s(i5,i6)/(propw(i1,i6)*propw(i2,i5))**2
      return
      end
