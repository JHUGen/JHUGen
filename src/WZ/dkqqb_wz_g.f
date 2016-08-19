      subroutine dkqqb_wz_g(p,msq)
      implicit none
C----Author R.K.Ellis June 2012
c----Matrix element for WZ production with radiation in the 
C----hadronic decay of the Z
C----averaged over initial colours and spins
C----massless final state particles
c     q(-p1)+qbar(-p2)-->n(p3)+ebar(p4)+q'(p5)+bar{q'}(p6)+g(p7)
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'zerowidth.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'ckm.f'
c      include 'anomcoup.f'
      include 'plabel.f'
      include 'srdiags.f'
      include 'nwz.f'
      include 'flqqb.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     & fac,s567,cotw,qqbsq,qbqsq
      double complex prop12,prop34,prop567,AWZqbq(2,2),AWZqqb(2,2)
      double complex cprop,props,propw,propz
      double complex Fa12(2,2),Fb12(2,2),Fc12(2,2),Fd12(2,2)
      double complex Fa21(2,2),Fb21(2,2),Fc21(2,2),Fd21(2,2)
      double complex qqb(4,2,2),qbq(4,2,2)
      double complex cs_z(2,2),cs_g(2,2)
      integer j,k,hg,hq

      fac=aveqq*xn*gw**8*gsq*2d0*CF
      flqqb=2
      srdiags=.true.      
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
      enddo
      enddo

c   We have --- f(p1) + f'(p2)--> nu(p3)+e^+(p4)+q(p5)+qbar(p6)+g(p7)
c   We have --- f(p1) + f'(p2)--> q(p3)+qbar(p4)+e^-(p5)+nu~(p6)+g(p7)

      call spinoru(7,p,za,zb)

  
      if     (zerowidth  .eqv. .true.) then
      write(6,*) 'dkqqb_ww_g.f:zerowidth should be false'
      stop
      elseif (zerowidth .neqv. .true.) then
c--   calculate propagators
      cotw=dsqrt((one-xw)/xw)
      s567=s(5,6)+s(5,7)+s(6,7)
      if     (zerowidth  .eqv. .true.) then
      prop12=s(1,2)/dcmplx(s(1,2)-wmass**2,wmass*wwidth)  
      prop34=s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      prop567=s567/dcmplx(s567-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
      elseif (zerowidth .neqv. .true.) then      
      prop12=dcmplx(s(1,2)/(s(1,2)-wmass**2))
      prop34=dcmplx(s(3,4)/(s(3,4)-wmass**2))
      prop567=dcmplx(s567/(s567-zmass**2))
      props=(s(1,2)-wmass**2)/dcmplx(s(1,2)-wmass**2,wmass*wwidth)
      propw=(s(3,4)-wmass**2)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      propz=(s(5,6)-zmass**2)/dcmplx(s(5,6)-zmass**2,zmass*zwidth)
      cprop=props*propw*propz
      endif

c--- DEBUG to compare with Madgraph
      prop12=s(1,2)/dcmplx(s(1,2)-wmass**2,wmass*wwidth)
      prop34=s(3,4)/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
      prop567=s567/dcmplx(s567-zmass**2,zmass*zwidth)
      cprop=dcmplx(1d0)
c--- DEBUG to compare with Madgraph

      endif
      
c-- couplings with or without photon pole
c---first index is quark helicity, second is weak isospin
      do j=1,2
      cs_z(1,j)=2d0*xw*l(j)*prop567
      cs_z(2,j)=2d0*xw*r(j)*prop567
      cs_g(1,j)=2d0*Q(j)*xw
      cs_g(2,j)=2d0*Q(j)*xw
      enddo


      if ((plabel(3) .eq. 'nl').and.(plabel(5) .eq. 'qj')) then
      call dkdrwz(1,2,3,4,5,6,7,Fa12,Fb12,Fc12,Fd12)
      call dkdrwz(2,1,3,4,5,6,7,Fa21,Fb21,Fc21,Fd21)
      elseif ((plabel(3) .eq. 'el').and.(plabel(5) .eq. 'qj')) then
      write(6,*) 'nwz',nwz
      call dkdrwz(1,2,3,4,5,6,7,Fb12,Fa12,Fd12,Fc12)
      call dkdrwz(2,1,3,4,5,6,7,Fb21,Fa21,Fd21,Fc21)
      endif

C----hq is the helicity of the 56 line
C----hg is the helicity of the gluon

      do hg=1,2
      do hq=1,2
      AWZqqb(hq,hg)=
     &  +prop34*(cs_z(hq,flqqb)*L(2)+cs_g(hq,flqqb)*Q(2))*Fa12(hq,hg) ! Coupling of strange line to up
     &  +prop34*(cs_z(hq,flqqb)*L(1)+cs_g(hq,flqqb)*Q(1))*Fb12(hq,hg) ! Coupling of strange line to down
     &  +prop12*prop34*(cotw*cs_z(hq,flqqb)+cs_g(hq,flqqb))*Fc12(hq,hg)

      AWZqbq(hq,hg)=
     &  +prop34*(cs_z(hq,flqqb)*L(2)+cs_g(hq,flqqb)*Q(2))*Fa21(hq,hg)
     &  +prop34*(cs_z(hq,flqqb)*L(1)+cs_g(hq,flqqb)*Q(1))*Fb21(hq,hg)
     &  +prop12*prop34*(cotw*cs_z(hq,flqqb)+cs_g(hq,flqqb))*Fc21(hq,hg)
      enddo
      enddo


      if (srdiags) then
c---we need supplementary diagrams for gauge invariance.

      if (nwz .eq. +1) then
      call dksrwz(1,2,3,4,5,6,7,qqb)   ! qqb
      call dksrwz(2,1,3,4,5,6,7,qbq)   ! qbq

      do hg=1,2
      do hq=1,2

      AWZqqb(hq,hg)=AWZqqb(hq,hg)
     & +prop12*qqb(1,hq,hg)*(cs_z(hq,flqqb)*ln+0d0*cs_g(hq,flqqb))
     & +prop12*qqb(2,hq,hg)*(cs_z(hq,flqqb)*le+qe*cs_g(hq,flqqb))
     & +prop12*prop34*qqb(2+flqqb,hq,hg)
      AWZqbq(hq,hg)=AWZqbq(hq,hg)
     & +prop12*qbq(1,hq,hg)*(cs_z(hq,flqqb)*ln+0d0*cs_g(hq,flqqb))
     & +prop12*qbq(2,hq,hg)*(cs_z(hq,flqqb)*le+qe*cs_g(hq,flqqb))
     & +prop12*prop34*qbq(2+flqqb,hq,hg)
      enddo
      enddo

      elseif (nwz .eq. -1) then
      call dksrwz(1,2,3,4,5,6,7,qqb)   ! qqb
      call dksrwz(2,1,3,4,5,6,7,qbq)   ! qbq
      do hg=1,2
      do hq=1,2
      AWZqqb(hq,hg)=AWZqqb(hq,hg)
     & +prop12*qqb(2,hq,hg)*(cs_z(hq,flqqb)*ln+0d0*cs_g(hq,flqqb))
     & +prop12*qqb(1,hq,hg)*(cs_z(hq,flqqb)*le+qe*cs_g(hq,flqqb))
     & +prop12*prop34*qqb(2+flqqb,hq,hg)
      AWZqbq(hq,hg)=AWZqbq(hq,hg)
     & +prop12*qbq(2,hq,hg)*(cs_z(hq,flqqb)*ln+0d0*cs_g(hq,flqqb))
     & +prop12*qbq(1,hq,hg)*(cs_z(hq,flqqb)*le+qe*cs_g(hq,flqqb))
     & +prop12*prop34*qbq(2+flqqb,hq,hg)
      enddo
      enddo
      endif
      endif

C-- Inclusion of width for W's a la Baur and Zeppenfeld with cprop.
      do hg=1,2
      do hq=1,2
      AWZqqb(hq,hg)=cprop*AWZqqb(hq,hg)
      AWZqbq(hq,hg)=cprop*AWZqbq(hq,hg)
      enddo
      enddo

      qqbsq=fac*
     & (+abs(AWZqqb(1,1))**2+abs(AWZqqb(1,2))**2
     &  +abs(AWZqqb(2,1))**2+abs(AWZqqb(2,2))**2)

      qbqsq=fac*
     & (+abs(AWZqbq(1,1))**2+abs(AWZqbq(1,2))**2
     &  +abs(AWZqbq(2,1))**2+abs(AWZqbq(2,2))**2)

      
      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
          if ((j .gt. 0) .and. (k .lt. 0)) then
            msq(j,k)=Vsq(j,k)*qqbsq
          elseif ((j .lt. 0) .and. (k .gt. 0)) then
            msq(j,k)=Vsq(j,k)*qbqsq
          endif
      enddo
      enddo
      return
      end


