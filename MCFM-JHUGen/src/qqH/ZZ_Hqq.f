      subroutine ZZ_Hqq(p,msq)
      implicit none
      include 'types.f'

c--- Weak Bosion Fusion by Z-Z exchange only
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
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'hdecaymode.f'
      integer:: j,k
      real(dp):: p(mxpart,4),fac,s34
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,
     & ud_ud_LL,udb_udb_LL,ud_ud_LR,udb_udb_LR,msqgamgam

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
c--- u(1)+d(2) -> u(5)+d(6)
c--- ub(1)+db(2) -> ub(5)+db(6)
      call msqpieces_zz(1,2,5,6,ud_ud_LL,ud_ud_LR)

C q-qbar and qbar-q
c--- u(1)+db(2) -> u(5)+db(6)
c--- ub(1)+d(2) -> ub(5)+d(6)
      call msqpieces_zz(1,6,5,2,udb_udb_LL,udb_udb_LR)

      do j=-nf,nf
      do k=-nf,nf
        if     ((j > 0) .and. (k < 0)) then
          msq(j,k)=fac*(
     &     +udb_udb_LL*((L(+j)*L(-k))**2+(R(+j)*R(-k))**2)
     &     +udb_udb_LR*((L(+j)*R(-k))**2+(R(+j)*L(-k))**2))
        elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=fac*(
     &     +udb_udb_LL*((L(-j)*L(k))**2+(R(-j)*R(k))**2)
     &     +udb_udb_LR*((L(-j)*R(k))**2+(R(-j)*L(k))**2))
        elseif ((j > 0) .and. (k > 0)) then
          msq(j,k)=fac*(
     &     +ud_ud_LL*((L(+j)*L(+k))**2+(R(+j)*R(+k))**2)
     &     +ud_ud_LR*((L(+j)*R(+k))**2+(R(+j)*L(+k))**2))
        elseif ((j < 0) .and. (k < 0)) then
          msq(j,k)=fac*(
     &     +ud_ud_LL*((L(-j)*L(-k))**2+(R(-j)*R(-k))**2)
     &     +ud_ud_LR*((L(-j)*R(-k))**2+(R(-j)*L(-k))**2))
        endif
      enddo
      enddo

      return
      end


      subroutine msqpieces_zz(i1,i2,i5,i6,zll,zlr)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      real(dp):: zll,zLR,htheta
      real(dp):: propz,x
      integer:: i1,i2,i5,i6
C--- define Heaviside theta function (=1 for x>0) and (0 for x < 0)
      htheta(x)=half+sign(half,x)

      propz(i1,i2)=sign(one,(s(i1,i2)-zmass**2))
     & *sqrt(sqrt(1._dp-xw)/xw/2._dp/zmass
     & *((s(i1,i2)-zmass**2)**2+htheta(s(i1,i2))*(zmass*zwidth)**2))

      zll=s(i1,i2)*s(i5,i6)/(propz(i1,i5)*propz(i2,i6))**2
      zlr=s(i1,i6)*s(i2,i5)/(propz(i1,i5)*propz(i2,i6))**2
      return
      end
