      subroutine qqb_gamgam_gvec(p,n,in,msq)
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
      include 'noglue.f'
C  in is the label of the momentum contracted with n
      integer:: j,k,in,h1,h2,h3,h4,i1,i2
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),facgg,gg,Qsum
      complex(dp):: amp(2,2,2,2),phase,phasec
      real(dp),parameter::statfac=0.5_dp

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      enddo
      enddo

c--- no contribution if we are omitting gg piece
      if (omitgg) return

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
