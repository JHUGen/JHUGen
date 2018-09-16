      subroutine epem3j_g(p,msq)
      implicit none
      include 'types.f'

c--- simple modification of qqb_w1jet_gs.f: permuted 1 and 4, 2 and 3
c--- to switch leptons with quarks and added a factor of Nc

c--- matrix element squared and averaged over initial colours and spins
c     q(-p1) + qbar(-p2) --> W + f(p5) + f(p6)
c                            |
c                            --> nu(p3) + e^+(p4)
c     where the fermions are either q(p5) and qbar(p6) [Qflag = .true.]
c                                or g(p5) and g(p6)    [Gflag = .true.]
c--- all momenta are incoming
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'nflav.f'
      include 'mmsq_cs.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & facgg,facqq,prop,qqbWgg2_slc,qqbWgg2_lc,qqbWgg2,qqb_ijkk(0:2)
      complex(dp):: qqb1(3),qqb2(3),qqb3(3),qqb4(3)
c--- we label the amplitudes by helicity (qqb1 ... qqb4)
c--- and by type of contribution qqb(1) ... qqb(n)
      integer:: j,k

c--- initialize matrix elements
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- set up spinors
      call spinoru(6,p,za,zb)
      prop=s(1,2)**2/((s(1,2)-wmass**2)**2+wmass**2*wwidth**2)
      facqq=4._dp*V*xn*gsq**2*(gwsq/2._dp)**2*aveqq*prop
      facgg=V*xn**2/four*(gwsq/2._dp)**2*gsq**2*prop

      call w2jetsq(4,3,2,1,5,6,za,zb,qqbWgg2)
      qqbWgg2_lc=mmsq_cs(1,1,1)+mmsq_cs(2,1,1)
      qqbWgg2_slc=mmsq_cs(0,1,1)

      qqbWgg2 = half*aveqq*facgg*qqbWgg2
      qqbWgg2_lc = half*aveqq*facgg*qqbWgg2_lc
      qqbWgg2_slc = half*aveqq*facgg*qqbWgg2_slc

c--- calculate four-quark amplitudes
c--- basic amplitudes - q qb --> W + g* (--> q qb) (amps 1 and 3)
c---                and q qb --> g* --> q (--> W q) qb (amps 2 and 4)
      call amp_epem3j_4q(4,3,2,1,5,6,qqb1(1),qqb2(1),qqb3(1),qqb4(1))

c--- now square these amplitudes separating into color structures
c   1) Amplitude
c   2) Amplitude with (5<-->6)
c   0) Interference between above
c
c--- q(i) qb(j) --> W + g* (--> q(k) qb(k)) with k != i,j
      qqb_ijkk(1)=abs(qqb1(1))**2+abs(qqb3(1))**2
      qqb_ijkk(2)=zip
      qqb_ijkk(0)=zip


c--- leading colour (Ca) goes in (0,0)
c--- subleading colour (2*Cf-Ca) goes in (0,1)
c--- Tr terms go in (1,0)
      msq(0,0)=qqbWgg2_lc
      msq(0,1)=qqbWgg2_slc
      msq(1,0)=2._dp*tr*real(nflav,dp)*facqq*qqb_ijkk(1)

      return
      end


      subroutine amp_epem3j_4q(i1,i2,i3,i4,i5,i6,amp1,amp2,amp3,amp4)
      implicit none
      include 'types.f'
c--- modification of this subroutine in qqb_w2jet.f:
c      subroutine amp_q_QbQ_qb(i1,i2,i5,i6,amp1,amp2,amp3,amp4)

      integer:: i1,i2,i3,i4,i5,i6
      complex(dp):: aqqb_zbb_new,amp1,amp2,amp3,amp4
c--- Amplitudes for q(i1) + qb(i2) --> qb(i6) + q(i5) + W (-> 3+4)
c
c--- the form of this function is taken from the subroutine msqzbb
c--- in qqb_zbb.f. See there for comments.
c--- We return four amplitudes, in pairs corresponding to diagrams
c--- where the W couples to both quark lines

c--- quark i5 is left-handed
      amp1=+aqqb_zbb_new(i1,i6,i5,i2,i3,i4)
      amp2=-conjg(aqqb_zbb_new(i5,i2,i1,i6,i4,i3))

c--- quark i5 is right-handed
      amp3=-aqqb_zbb_new(i1,i5,i6,i2,i3,i4)
      amp4=-conjg(aqqb_zbb_new(i5,i1,i2,i6,i4,i3))

      return
      end

