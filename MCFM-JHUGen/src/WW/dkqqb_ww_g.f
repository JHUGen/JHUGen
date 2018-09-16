      subroutine dkqqb_ww_g(p,msq)
      implicit none
      include 'types.f'

C----Author R.K.Ellis June 2012
c----Matrix element for WW production with radiation in hadronic decay
C----averaged over initial colours and spins
C----massless final state particles
c     q(-p1)+qbar(-p2)-->n(p3)+ebar(p4)+q'(p5)+bar{q'}(p6)+g(p7)
c--- note that non-leptonic W decays do not include scattering diagrams
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
c      include 'anomcoup.f'
      include 'plabel.f'
      include 'srdiags.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & fac,s567,tanw,l3,l4,l5,l6,q3,q4,q5,q6
      real(dp), parameter :: mp(nf)=(/-1._dp,+1._dp,-1._dp,+1._dp,-1._dp/)
      complex(dp):: prop12,prop34,prop567,AWW(2,2)
      complex(dp):: propwp,propwm,propzg,cprop
      complex(dp):: Fa12(2,2),Fb12(2,2),Fc12(2,2),Fd12(2,2)
      complex(dp):: Fa21(2,2),Fb21(2,2),Fc21(2,2)
      complex(dp):: qqb(4,2,2),qbq(4,2,2)
      complex(dp):: cs_z(2,2),cs_g(2,2)
      integer:: j,k,jk,tjk,hg,hq

      fac=aveqq*xn*gw**8*gsq*2._dp*CF
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
      enddo
      enddo

c   We have --- f(p1) + f'(p2)--> nu(p3)+e^+(p4)+q(p5)+qbar(p6)+g(p7)
c   We have --- f(p1) + f'(p2)--> q(p3)+qbar(p4)+e^-(p5)+nu~(p6)+g(p7)

      call spinoru(7,p,za,zb)

      if ((plabel(3) == 'nl').and.(plabel(5) == 'qj')) then
c--- this is the normal case
         l4=le
         q4=qe
         l3=ln
         q3=0._dp
         l5=l(1)
         q5=q(1)
         l6=l(2)
         q6=q(2)
         fac=2._dp*xn*fac
      else if ((plabel(3) == 'el').and.(plabel(5) == 'qj')) then
c--- this is the swapped case
         l4=ln
         q4=0._dp
         l3=le
         q3=qe
         l5=l(2)
         q5=q(2)
         l6=l(1)
         q6=q(1)
         fac=2._dp*xn*fac
      else
         return
      endif

      if     (zerowidth  .eqv. .true.) then
      write(6,*) 'dkqqb_ww_g.f:zerowidth should be false'
      stop
      elseif (zerowidth .neqv. .true.) then
      s567=s(5,6)+s(5,7)+s(6,7)
      prop12=cplx1(s(1,2)/(s(1,2)-zmass**2))
      prop34=cplx1(s(3,4)/(s(3,4)-wmass**2))
      prop567=cplx1(s567/(s567-wmass**2))
      propwm=(s(3,4)-wmass**2)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      propwp=(s567-wmass**2)/cplx2(s567-wmass**2,wmass*wwidth)
      propzg=(s(1,2)-zmass**2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
      cprop=propwp*propwm*propzg
c--- debug
c----for Madgraph
c      prop12=s(1,2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
c      prop34=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
c      prop567=s567/cplx2(s567-wmass**2,wmass*wwidth)
c      cprop=cone
      endif

c-- couplings with or without photon pole
c---first index is quark helicity, second is weak isospin
      tanw=2._dp*xw/sin2w
      do j=1,2
      cs_z(1,j)=+mp(j)*l(j)*sin2w*prop12
      cs_z(2,j)=-mp(j)*2._dp*Q(j)*xw*prop12
      cs_g(1,j)=+mp(j)*2._dp*Q(j)*xw
      cs_g(2,j)=+mp(j)*2._dp*Q(j)*xw
      enddo


      if ((plabel(3) == 'nl').and.(plabel(5) == 'qj')) then
      call dkWWamps(1,2,3,4,5,6,7,Fa12,Fb12,Fc12,Fd12)
      elseif ((plabel(3) == 'el').and.(plabel(5) == 'qj')) then
      call dkWWamps(1,2,3,4,5,6,7,Fb12,Fa12,Fd12,Fc12)
      endif

C--calculate qbq from qqb
      do j=1,2
      do hg=1,2
      Fa21(j,hg)=Fa12(3-j,hg)
      Fb21(j,hg)=Fb12(3-j,hg)
      Fc21(j,hg)=Fc12(3-j,hg)
      enddo
      enddo
C----zero out right handed 12 quark line pieces that do not contribute
C----because W's attach to 12 line
      Fa12(2,:)=czip
      Fa21(2,:)=czip
      Fb12(2,:)=czip
      Fb21(2,:)=czip

      do j=-nf,nf
         k=-j
C---as we begin every different j,k flavour, zero out AWW
         AWW(:,:)=czip
c--   Exclude gluon-gluon initial state
         if (j==0) go to 20
         jk=max(j,k)

         do hg=1,2
            do hq=1,2
               if (j > 0) then
                  if  (tau(jk) == +1._dp) then
                     AWW(hq,hg)=prop567*prop34*(mp(jk)*Fb12(hq,hg)
     &                    +(cs_z(hq,2)+cs_g(hq,2))*Fc12(hq,hg))
                  elseif     (tau(jk) == -1._dp) then
                     AWW(hq,hg)=prop567*prop34*(mp(jk)*Fa12(hq,hg)
     &                    +(cs_z(hq,1)+cs_g(hq,1))*Fc12(hq,hg))
                  endif
               elseif (j < 0) then
                  if  (tau(jk) == +1._dp) then
                     AWW(hq,hg)=prop567*prop34*(mp(jk)*Fb21(hq,hg)
     &                    +(cs_z(hq,2)+cs_g(hq,2))*Fc21(hq,hg))
                  elseif     (tau(jk) == -1._dp) then
                     AWW(hq,hg)=prop567*prop34*(mp(jk)*Fa21(hq,hg)
     &                    +(cs_z(hq,1)+cs_g(hq,1))*Fc21(hq,hg))
                  endif
               endif
            enddo

            if (srdiags) then
c---we need supplementary diagrams for gauge invariance.

               call c7tree(1,2,3,4,5,6,7,qqb) ! qqb
               call c7tree(2,1,3,4,5,6,7,qbq) ! qbq

C---tjk is equal to 2 (u,c) or 1 (d,s,b)
               tjk=2-mod(abs(jk),2)
               do hq=1,2
                  if (j > 0) then
                     AWW(hq,hg)=AWW(hq,hg)
     &      +prop567*qqb(1,hq,hg)*(tanw*cs_z(hq,tjk)*l4+cs_g(hq,tjk)*q4)
     &      +prop567*qqb(2,hq,hg)*(tanw*cs_z(hq,tjk)*l3+cs_g(hq,tjk)*q3)
     &      +prop34 *qqb(3,hq,hg)*(tanw*cs_z(hq,tjk)*l5+cs_g(hq,tjk)*q5)
     &      +prop34 *qqb(4,hq,hg)*(tanw*cs_z(hq,tjk)*l6+cs_g(hq,tjk)*q6)
                  elseif (j < 0) then
                     AWW(hq,hg)=AWW(hq,hg)
     &      +prop567*qbq(1,hq,hg)*(tanw*cs_z(hq,tjk)*l4+cs_g(hq,tjk)*q4)
     &      +prop567*qbq(2,hq,hg)*(tanw*cs_z(hq,tjk)*l3+cs_g(hq,tjk)*q3)
     &      +prop34 *qbq(3,hq,hg)*(tanw*cs_z(hq,tjk)*l5+cs_g(hq,tjk)*q5)
     &      +prop34 *qbq(4,hq,hg)*(tanw*cs_z(hq,tjk)*l6+cs_g(hq,tjk)*q6)
                  endif
               enddo

            endif

C-- Inclusion of width for W's a la Baur and Zeppenfeld with cprop.
            AWW(1,hg)=cprop*AWW(1,hg)
            AWW(2,hg)=cprop*AWW(2,hg)
         enddo
         msq(j,k)=fac*
     &        (+abs(AWW(1,1))**2+abs(AWW(1,2))**2
     &        +abs(AWW(2,1))**2+abs(AWW(2,2))**2)

 20      continue
      enddo

      return
      end


