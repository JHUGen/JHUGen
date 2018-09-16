      subroutine qqb_ZH_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
C***********************************************************************
c     Author: R.K. Ellis                                               *
c     September, 1999.                                                 *
c     Matrix element for Z production                                  *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     u(-p1)+dbar(-p2)--> g(p7)+ Z^+(l(p3)+a(p4)) +H(p5,p6)            *
C***********************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'hdecaymode.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'nflav.f'
      include 'hbbparams.f'
      integer:: j,k,in
      integer:: ng
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: z1jetn,fac,p1p2(-1:1,-1:1),n(4),facdecay,hdecay,
     & msqhbb,msqhtautau,msqhgamgam,s56,s127
      complex(dp):: prop
      real(dp):: sH
      integer ig
      msq(:,:)=zip

      if(hdecaymode.ne.'wpwm') then
         ig=7
      else
         ig=9
      endif
      ng=ig
      call dotem(ig,p,s)

      fac=16._dp*cf*xn*esq**2*gsq
      s127=s(1,2)+s(1,ig)+s(2,ig)
      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      prop=prop/cplx2((s127-zmass**2),zmass*zwidth)

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          sH=s(5,6)+two*mtau**2
          hdecay=msqhtautau(sH)
      elseif (hdecaymode == 'bqba') then
          sH=s(5,6)+two*mb**2
          hdecay=msqhbb(sH)
c---  adjust for fixed H->bb BR if necessary
          if (FixBrHbb) then
             hdecay=hdecay*GamHbb/GamHbb0
          endif
       elseif (hdecaymode == 'gaga') then
          sH=s(5,6)
          hdecay=msqhgamgam(sH)
       elseif (hdecaymode =='wpwm') then
          sH=s(5,6)+s(5,7)+s(5,8)+s(6,7)+s(6,8)+s(7,8)
          call hwwdecay(p,5,6,7,8,hdecay)

       else
          write(6,*) 'Unimplemented process in qqb_higgs'
          stop
      endif
      hdecay=hdecay/((sH-hmass**2)**2+(hmass*hwidth)**2)
      facdecay=gwsq*wmass**2/(one-xw)**2*hdecay

      fac=fac*facdecay

      p1p2(:,:)=zip

      if (in == 1) then
      p1p2(0,-1)=-aveqg*fac*z1jetn(ng,2,1,p,n)
      p1p2(0,+1)=-aveqg*fac*z1jetn(2,ng,1,p,n)
      elseif (in == 2) then
      p1p2(+1,0)=-aveqg*fac*z1jetn(1,ng,2,p,n)
      p1p2(-1,0)=-aveqg*fac*z1jetn(ng,1,2,p,n)
      elseif (in == ng) then
      p1p2(-1,1)=+aveqq*fac*z1jetn(2,1,ng,p,n)
      p1p2(1,-1)=+aveqq*fac*z1jetn(1,2,ng,p,n)
      endif

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=zip
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+(abs(Q(j)*q1+L(j)*l1*prop)**2
     &              +abs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(1,-1)
     &             +(abs(Q(j)*q1+L(j)*r1*prop)**2
     &              +abs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,1)
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+(abs(Q(k)*q1+L(k)*l1*prop)**2
     &              +abs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(-1,1)
     &             +(abs(Q(k)*q1+L(k)*r1*prop)**2
     &              +abs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(1,-1)
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+(abs(Q(j)*q1+L(j)*l1*prop)**2
     &              +abs(Q(j)*q1+R(j)*r1*prop)**2)*p1p2(+1,0)
     &             +(abs(Q(j)*q1+L(j)*r1*prop)**2
     &              +abs(Q(j)*q1+R(j)*l1*prop)**2)*p1p2(-1,0)
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+(abs(Q(-j)*q1+L(-j)*l1*prop)**2
     &              +abs(Q(-j)*q1+R(-j)*r1*prop)**2)*p1p2(-1,0)
     &             +(abs(Q(-j)*q1+L(-j)*r1*prop)**2
     &              +abs(Q(-j)*q1+R(-j)*l1*prop)**2)*p1p2(+1,0)
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+(abs(Q(k)*q1+L(k)*l1*prop)**2
     &              +abs(Q(k)*q1+R(k)*r1*prop)**2)*p1p2(0,+1)
     &             +(abs(Q(k)*q1+L(k)*r1*prop)**2
     &              +abs(Q(k)*q1+R(k)*l1*prop)**2)*p1p2(0,-1)
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+(abs(Q(-k)*q1+L(-k)*l1*prop)**2
     &              +abs(Q(-k)*q1+R(-k)*r1*prop)**2)*p1p2(0,-1)
     &             +(abs(Q(-k)*q1+L(-k)*r1*prop)**2
     &              +abs(Q(-k)*q1+R(-k)*l1*prop)**2)*p1p2(0,+1)
      endif

   19 continue
      enddo
      enddo

      return
      end



