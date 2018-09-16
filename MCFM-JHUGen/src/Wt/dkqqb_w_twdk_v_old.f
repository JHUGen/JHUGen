      subroutine dkqqb_w_twdk_v_old(p,msqv)
      implicit none
      include 'types.f'
************************************************************************
*     Author: Francesco Tramontano                                     *
*     February, 2005.                                                  *
*     Virtual corrections in the decay, averaged over initial          * 
*      colours and spins                                               *
*--- W+t production (nwz=-1)                                           *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> nu(p5) + e^+(p6) + b(p7)          *
*                            |                                         *
*                            --> e^-(p3) + nubar(p4)                   *
*--- W+tbar production (nwz=+1)                                        *
*     q(-p1) + qbar(-p2) --> W + t(p567)                               *
*                            |   |                                     *
*                            |   --> e^-(p5) + nubar(p6) + b(p7)       *
*                            |                                         *
*                            --> nu(p3) + e^+(p4)                      *
************************************************************************
      
      integer:: i,j,k
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'scheme.f'
      include 'masses.f'
      include 'nwz.f'
      integer:: i3,i4,i5,i6,iq
      real(dp):: p(mxpart,4),msqv(-nf:nf,-nf:nf),
     & msq(-nf:nf,-nf:nf),virtgqdk,gq,qg,t(4),fac,corr,nloratiotopdecay
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

      scheme='dred'
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=zero
      enddo
      enddo

c--- set up lepton variables depending on whether it's t or tbar
      if     (nwz == -1) then
        i3=3
	i4=4
	i5=5
	i6=6
	iq=1 ! top quark
      elseif (nwz == +1) then
        i3=4
	i4=3
	i5=6
	i6=5
	iq=-1 ! antitop quark
      else
        write(6,*) 'Error in qqb_w_twdk_vdk, nwz is not +1 or -1 : ',nwz
	stop
      endif

      call qqb_w_twdk(p,msq) 
      corr=nloratiotopdecay(mt,mb,wmass,wwidth)
c---fill matrices of spinor products

      do i=1,4
        t(i)=p(5,i)+p(6,i)+p(7,i)
      enddo

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,t)

      fac=ason2pi*cf
      fac=aveqg*gwsq**4*gsq*V*fac
      gq=virtgqdk(1,2,i3,i4,i5,i6,7)
      qg=virtgqdk(2,1,i3,i4,i5,i6,7)

      do j=-nf,nf,nf
      do k=-nf,nf,nf
      msqv(j,k)=zero
      if     ((j == 5*iq) .and. (k == 0)) then
          msqv(j,k)=fac*qg-corr*msq(j,k)
      elseif ((j == 0) .and. (k == 5*iq)) then
          msqv(j,k)=fac*gq-corr*msq(j,k)
      endif
      enddo
      enddo
      
      return
      end


      function virtgqdk(ig,is,ie,in,jn,je,jc)
      implicit none
      include 'types.f'
      real(dp):: virtgqdk
      

      integer:: ig,is,ie,in,je,jn,jc
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      real(dp):: snec,taugt,prop,mtsq,cv,ct,c1
      complex(dp):: amp(2),amp1(2),ampho(2)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba

      mtsq=mt**2
      snec=+s(jn,je)+s(je,jc)+s(jc,jn)
      taugt=s(ig,je)+s(ig,jn)+s(ig,jc)

      call coefsdk(s(jn,je),mtsq,ct,cv,c1)
      prop=(s(ie,in)-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((s(je,jn)-wmass**2)**2+(wmass*wwidth)**2)
      prop=prop*((snec-mt**2)**2+(mt*twidth)**2)
      
c--- tree level amplitudes for the decay t->Wb, "c0 contribution"
      amp(1) =
     &  - za(ie,in)*za(jc,jn)*zb(is,in)**2/zb(ig,is)*zab(ig,je)*
     & taugt**(-1)
     &  + za(ig,ie)*za(jc,jn)*zb(is,in)*zb(is,je)/zb(ig,is)*mtsq*
     & taugt**(-1)
     &

      amp(2) =
     &  - za(is,ie)*za(jc,jn)*zb(is,in)*zb(ig,je)/za(ig,is)*mtsq*
     & taugt**(-1)
     &  + za(jc,jn)*zb(is,in)/za(ig,is)*zab(is,je)*zab(ie,ig)*
     & taugt**(-1)
     &  + za(jc,jn)*zb(ig,in)/za(ig,is)*zab(ie,je)
     &

c--- virtual amplitudes for the decay t->Wb, "c1 contribution"
      amp1(1) =
     &  - za(ie,in)*za(ig,jc)*za(jc,jn)*zb(is,in)**2*zb(jc,je)/zb(ig,
     & is)*taugt**(-1)
     &  + za(ig,ie)*za(jc,jn)*zb(is,in)*zb(jc,je)/zb(ig,is)*zab(jc,is)
     & *taugt**(-1)
     &

      amp1(2) =
     &  - za(is,ie)*za(jc,jn)*zb(is,in)*zb(jc,je)/za(ig,is)*zab(jc,ig)
     & *taugt**(-1)
     &  + za(is,jc)*za(jc,jn)*zb(is,in)*zb(jc,je)/za(ig,is)*zab(ie,ig)
     & *taugt**(-1)
     &  + za(ie,jc)*za(jc,jn)*zb(ig,in)*zb(jc,je)/za(ig,is)
     &

      ampho(1)=cplx2(ct+cv)*amp(1)+cplx2(half*c1)*amp1(1)
      ampho(2)=cplx2(ct+cv)*amp(2)+cplx2(half*c1)*amp1(2)

      virtgqdk=real(amp(1)*conjg(ampho(1)))/prop
     &        +real(amp(2)*conjg(ampho(2)))/prop
      return
      end


