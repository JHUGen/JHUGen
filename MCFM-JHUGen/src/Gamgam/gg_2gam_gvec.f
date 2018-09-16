      subroutine gg_2gam_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c--- For a description of the procedure implemented in this routine,
c--- see Bern, Dixon and Schmidt, hep-ph/0206194, Eqs.(9)-(11).
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zprods_com.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in,h1,h2,h3,h4,i1,i2
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),facgg,gg,Qsum
      complex(dp):: amp(2,2,2,2),phase,phasec
      real(dp),parameter::statfac=0.5_dp

      msq(:,:)=0._dp

c--- this amplitude doesn't contain any n dependence
c--- (it has been explicitly written out), and the emitted momentum
c--- (label 5) has been passed in from dipolesub via "n"
      do j=1,4
      p(5,j)=n(j)
      enddo

      call spinoru(5,p,za,zb)      

      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2 
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      if     (in == 1) then
        i1=1
        i2=2
      elseif (in == 2) then
        i1=2
        i2=1
      else
        write(6,*) 'Unexpected value of in in qqb_gamgam_gvec.f: in=',in
      endif
      phase=za(i1,5)*zb(5,i2)*za(i2,i1)
     &    /(zb(i1,5)*za(5,i2)*zb(i2,i1))
      phasec=conjg(phase)
            
      call gg_gaga_amps(1,2,3,4,amp)
      gg=0._dp
      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
c--- first, LO term
      gg=gg+abs(amp(h1,h2,h3,h4))**2
c--- now subtract correlated term
      if     (in == 1) then
        if (h1 == 2) then
          gg=gg-real(phase*amp(h1,h2,h3,h4)*conjg(amp(3-h1,h2,h3,h4)))
        endif
        if (h1 == 1) then
          gg=gg-real(phasec*amp(h1,h2,h3,h4)*conjg(amp(3-h1,h2,h3,h4)))
        endif
      elseif (in == 2) then
        if (h2 == 2) then
          gg=gg-real(phase*amp(h1,h2,h3,h4)*conjg(amp(h1,3-h2,h3,h4)))
        endif
        if (h2 == 1) then
          gg=gg-real(phasec*amp(h1,h2,h3,h4)*conjg(amp(h1,3-h2,h3,h4)))
        endif
      endif
      enddo
      enddo
      enddo
      enddo
      
      msq(0,0)=avegg*V*facgg**2*statfac*gg

c--- divide by extra factor of two to compensate for normalization
c--- of "subv" in dipolesub.f
      msq(0,0)=msq(0,0)/2._dp

      return
      end

      subroutine gg_gaga_amps(j1,j2,j3,j4,amp)
      implicit none
      include 'types.f'
c--- Implementation of amplitudes for g g -> gamma gamma
c--- from Bern, De Freitas and Dixon, arXiv:hep-ph/0109078
c--- These differ from those already coded in msqgggaga (qqb_gamgam.f)
c--- in that overall phases are included
c--- Additional factor added: sum over quark charges
            
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      real(dp):: ss,tt,uu
      complex(dp):: amp(2,2,2,2),phase
      integer:: j1,j2,j3,j4,h1,h2,h3,h4

      ss=s(j1,j2)
      tt=s(j1,j3)
      uu=s(j2,j3)

      do h1=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      amp(h1,h2,h3,h4)=cone
      enddo
      enddo
      enddo
      enddo

c--- definitions appearing in Eq. (3.16) [no overall phases]
      amp(1,1,2,2)=
     &  (-0.5_dp*(tt**2+uu**2)/ss**2*(log(tt/uu)**2+pisq)
     & -(tt-uu)/ss*log(tt/uu)-cone)
      amp(1,2,1,2)=
     &  (-0.5_dp*(tt**2+ss**2)/uu**2*log(-tt/ss)**2
     & -(tt-ss)/uu*log(-tt/ss)-cone
     & -impi*((tt**2+ss**2)/uu**2*log(-tt/ss)+(tt-ss)/uu))
      amp(2,1,1,2)=
     &  (-0.5_dp*(uu**2+ss**2)/tt**2*log(-uu/ss)**2
     & -(uu-ss)/tt*log(-uu/ss)-cone
     & -impi*((uu**2+ss**2)/tt**2*log(-uu/ss)+(uu-ss)/tt))

c--- now apply phases from Eq. (3.5)
      phase=
     &  im*zb(j1,j2)*zb(j3,j4)/(za(j1,j2)*za(j3,j4))
      amp(1,1,1,1)=amp(2,2,2,2)*conjg(phase)
      amp(2,2,2,2)=amp(2,2,2,2)*phase
      
      phase=
     &  im*za(j1,j2)*za(j1,j4)*zb(j2,j4)/(za(j3,j4)*za(j2,j3)*za(j2,j4))
      amp(2,1,1,1)=amp(1,2,2,2)*conjg(phase)
      amp(1,2,2,2)=amp(1,2,2,2)*phase
      
      phase=
     &  im*za(j2,j3)*za(j2,j4)*zb(j3,j4)/(za(j1,j4)*za(j3,j1)*za(j3,j4))
      amp(1,2,1,1)=amp(2,1,2,2)*conjg(phase)
      amp(2,1,2,2)=amp(2,1,2,2)*phase
      
      phase=
     &  im*za(j3,j2)*za(j3,j4)*zb(j2,j4)/(za(j1,j4)*za(j2,j1)*za(j2,j4))
      amp(1,1,2,1)=amp(2,2,1,2)*conjg(phase)
      amp(2,2,1,2)=amp(2,2,1,2)*phase

      phase=
     &  im*za(j1,j2)*zb(j3,j4)/(zb(j1,j2)*za(j3,j4))
      amp(2,2,1,1)=amp(1,1,2,2)*conjg(phase)
      amp(1,1,2,2)=amp(1,1,2,2)*phase
      
      phase=
     &  im*za(j1,j3)*zb(j2,j4)/(zb(j1,j3)*za(j2,j4))
      amp(2,1,2,1)=amp(1,2,1,2)*conjg(phase)
      amp(1,2,1,2)=amp(1,2,1,2)*phase
      
      phase=
     &  im*za(j2,j3)*zb(j1,j4)/(zb(j2,j3)*za(j1,j4))
      amp(1,2,2,1)=amp(2,1,1,2)*conjg(phase)
      amp(2,1,1,2)=amp(2,1,1,2)*phase
      
c--- the following amplitude has been added by me (permutation of above)
      phase=
     &  im*za(j4,j2)*za(j4,j3)*zb(j2,j3)/(za(j1,j3)*za(j2,j1)*za(j2,j3))
      amp(1,1,1,2)=amp(2,2,2,1)*conjg(phase)
      amp(2,2,2,1)=amp(2,2,2,1)*phase
      
      return
      end
      
