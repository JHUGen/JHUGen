      subroutine qqb_zz_g(P,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  q'(p4)+bar{q'}(p5) + n(p6)+ebar(p7)+ g(p3)
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      include 'interference.f'
      include 'pchoice.f'
      integer:: jk,hq,h34,h56,hg,jp,kp,ii,nmax
      real(dp):: fac,fac1,q34,q56,s127
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: ave,v34(2),v56(2),rescale1,rescale2,oprat
      complex(dp):: aq12(2,2,2,2),aq34(2,2,2,2),aq56(2,2,2,2)
      complex(dp):: qa12(2,2,2,2),qa34(2,2,2,2),qa56(2,2,2,2)
      complex(dp):: gq12(2,2,2,2),gq34(2,2,2,2),gq56(2,2,2,2)
      complex(dp):: ag12(2,2,2,2),ag34(2,2,2,2),ag56(2,2,2,2)
      complex(dp):: ga12(2,2,2,2),ga34(2,2,2,2),ga56(2,2,2,2)
      complex(dp):: qg12(2,2,2,2),qg34(2,2,2,2),qg56(2,2,2,2)
      complex(dp):: Uncrossed(-nf:nf,-nf:nf,2,2,2,2)
      complex(dp):: amp,prop34,prop56,prop127
      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)


C-set Uncrossed array to zero
      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=zip
      Uncrossed(-nf:nf,-nf:nf,:,:,:,:)=czip
      enddo
      enddo

      fac=-four*esq**2
      fac1=two*cf*xn*gsq

      v34(1)=l1
      v34(2)=r1
      q34=q1
      v56(1)=l2
      v56(2)=r2
      q56=q2

C----setup factor to avoid summing over too many neutrinos
C----if coupling enters twice
      if (q34 == zip) then
      rescale1=one/sqrt(three)
      else
      rescale1=one
      endif
      if (q56 == zip) then
      rescale2=one/sqrt(three)
      else
      rescale2=one
      endif

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,p,za,zb)

      if (interference) then
         nmax=2
      else
         nmax=1
      endif

      do ii=1,nmax

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(5,i6(ii))/cplx2(s(5,i6(ii))-zmass**2,zmass*zwidth)

c--- Amplitude returned with arguments (hq,h34,h56,h7)
c---case qbar-q
      call zzgamp(1,2,3,i4(ii),5,i6(ii),7,za,zb,aq12,aq34,aq56)
c---case q-qbar
      call zzgamp(2,1,3,i4(ii),5,i6(ii),7,za,zb,qa12,qa34,qa56)
c---case qbar-g
      call zzgamp(1,7,3,i4(ii),5,i6(ii),2,za,zb,ag12,ag34,ag56)
c---case q-g
      call zzgamp(7,1,3,i4(ii),5,i6(ii),2,za,zb,qg12,qg34,qg56)
c---case g-q
      call zzgamp(7,2,3,i4(ii),5,i6(ii),1,za,zb,gq12,gq34,gq56)
c---case g-qbar
      call zzgamp(2,7,3,i4(ii),5,i6(ii),1,za,zb,ga12,ga34,ga56)


C---calculate over a limited flavour range (-2:2)
      do j=-2,2
      do k=-2,2
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

c---determine averaging factor for different channels
c     vsymfact=symmetry factor
      if ((j == 0) .or. (k == 0)) then
        jk=j+k
        ave=aveqg*vsymfact
      else
        jk=max(j,k)
        ave=aveqq*vsymfact
      endif

      if (jk == 0) goto 19

      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

      amp=zip

c---case qbar-q
      if    ((j < 0).and.(k > 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(k)+q56*q(k))
     & *(prop34*v34(h34)*l(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(k)+q56*q(k))
     & *(prop34*v34(h34)*r(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      endif

c---case q-qbar
      elseif((j > 0).and.(k < 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(j)+q56*q(j))
     & *(prop34*v34(h34)*l(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(j)+q56*q(j))
     & *(prop34*v34(h34)*r(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      endif

c---case qbar-g
      elseif((j < 0).and.(k == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*l(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*r(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      endif

c---case g-qbar
      elseif((k < 0).and.(j == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*l(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*r(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      endif

c---case q-g
      elseif((j > 0).and.(k == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*l(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*r(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      endif

      elseif((k > 0).and.(j == 0)) then
c---case g-q
      if (hq == 1) then
      amp=(prop56*v56(h56)*l(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*l(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*l(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*l(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*r(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*r(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*r(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*r(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      endif

      endif

      amp=amp*fac

c      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2
c      !set-up interference terms
c      if ((interference).and.(ii==1)) then
c      Uncrossed(j,k,hq,h34,h56,hg)=amp
c      elseif (ii==2) then
c      if (h34==h56) then
c      msq(j,k)=msq(j,k)
c     &    -fac1*two*ave*real(conjg(amp)*Uncrossed(j,k,hq,h34,h56,hg))
c      endif
c      endif


      if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          Uncrossed(j,k,hq,h34,h56,hg)=amp
        else
          if (h34 == h56) then
            oprat=one
     &           -two*real(conjg(amp)*Uncrossed(j,k,hq,h34,h56,hg))
     &            /(abs(amp)**2+abs(Uncrossed(j,k,hq,h34,h56,hg))**2)
          else
            oprat=one
          endif
          if (bw34_56) then
            msq(j,k)=msq(j,k)
     &        +fac1*ave*two*abs(Uncrossed(j,k,hq,h34,h56,hg))**2*oprat
          else
            msq(j,k)=msq(j,k)+fac1*ave*two*abs(amp)**2*oprat
          endif
        endif
      endif

      enddo  ! endloop hg
      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq



   19 continue
      enddo  !endloop j
      enddo  !endloop k
      enddo  !endloop ii


C---extend to full flavour range
      do j=-nf,nf
      do k=-nf,nf
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 20

      msq(j,k)=msq(jkswitch(j),jkswitch(k))

 20   continue
      enddo
      enddo

      return
      end







