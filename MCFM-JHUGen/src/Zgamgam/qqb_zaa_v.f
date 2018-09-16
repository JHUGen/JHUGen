****************************************************************
*   Virtual matrix element for
*   0 -> q(-p1) + qb(-p4) + a(p2) + a(p3) + lb(p5) + l(p6)
*   l = charged lepton/neutrino
****************************************************************
      subroutine qqb_zaa_v(p,msqv)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scheme.f'

      include 'scale.f'
      include 'epinv2.f'
      include 'epinv.f'

      real(dp):: p(mxpart,4),msqv(-nf:nf,-nf:nf)
      complex(dp):: a60h(4,16),a6vh(4,16)
      real(dp):: qqb(2),qbq(2),pbdk(mxpart,4)
      integer:: i,j,k
c-----set scheme
      scheme='dred'
c-----swap momenta
c-----implemented (ala bdk)
c----- 0 -> q(p1) + qb(p4) + gam(p2) + gam(p3) + lb(p5) + l(p6)
c-----mcfm input:
c----- 0 -> q(p1) + qb(p2) + l(p3) + lb(p4) + gam(p5) + gam(p6)
      do i=1,4
         pbdk(1,i)=p(2,i)
         pbdk(4,i)=p(1,i)
         pbdk(2,i)=p(6,i)
         pbdk(3,i)=p(5,i)
         pbdk(5,i)=p(4,i)
         pbdk(6,i)=p(3,i)
      enddo
c-----evaluate spinor products,invariants
      call spinoru(6,pbdk,za,zb)
c-----initialize matelem
      do i=1,2
         qqb(i)=zip
         qbq(i)=zip
      enddo
c-----compute color ordered amplitudes
c-----then compute squared matrix element
c-----qqbar subprocess
      call zaa_a60h(1,2,3,4,5,6,za,zb,a60h)
      call zaa_a6vh(1,2,3,4,5,6,za,zb,a6vh)
      call zaa_m6vsq(2,a60h,a6vh,qqb(2))
      call zaa_m6vsq(1,a60h,a6vh,qqb(1))
c-----qbarq subprocess
      call zaa_a60h(4,2,3,1,5,6,za,zb,a60h)
      call zaa_a6vh(4,2,3,1,5,6,za,zb,a6vh)
      call zaa_m6vsq(2,a60h,a6vh,qbq(2))
      call zaa_m6vsq(1,a60h,a6vh,qbq(1))
c-----initialize msq
      do i=-nf,nf
      do j=-nf,nf
         msqv(i,j)=zip
      enddo
      enddo
c-----fill msq
      msqv(1,-1)=qqb(1)
      msqv(3,-3)=qqb(1)
      msqv(5,-5)=qqb(1)
      msqv(2,-2)=qqb(2)
      msqv(4,-4)=qqb(2)
      msqv(-1,1)=qbq(1)
      msqv(-3,3)=qbq(1)
      msqv(-5,5)=qbq(1)
      msqv(-2,2)=qbq(2)
      msqv(-4,4)=qbq(2)
c-----donehere
      return
      end

***********************************************
* virtual matrix element squared
* qqb_zaa, for initial state flavor qi
***********************************************
      subroutine zaa_m6vsq(qi,a60h,a6vh,msqv)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'ipsgen.f'
      real(dp):: t
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a60h(4,16),m60hA(16),m60hB(16),m60hC(16),
     & m60hD(16),a6vh(4,16),m6vhA(16),m6vhB(16),m6vhC(16),m6vhD(16)
      complex(dp):: propzQQ,propzLL,propzQL(3:4)
      real(dp):: m6sqvh(16),msqv
      real(dp):: m6sqvhAA(16),m6sqvhBB(16),m6sqvhCC(16),
     .m6sqvhDD(16),m6sqvhAB(16),m6sqvhAC(16),m6sqvhAD(16),m6sqvhBC(16),
     .m6sqvhBD(16),m6sqvhCD(16)
      real(dp):: iwdth,iph,izz,qq1
      integer:: qi,i,j,k
c-----Z propagator
      iwdth = one
c-----QQ propagator
      propzQQ=s(5,6)/cplx2(s(5,6)-zmass**2,iwdth*zwidth*zmass)
c-----LL propagator
      propzLL=s(1,4)/cplx2(s(1,4)-zmass**2,iwdth*zwidth*zmass)
c-----QL propagator (2,3)
      propzQL(3)=t(1,4,2)/cplx2(t(1,4,2)-zmass**2,iwdth*zwidth*zmass)
c-----QL propagator (3,2)
      propzQL(4)=t(1,4,3)/cplx2(t(1,4,3)-zmass**2,iwdth*zwidth*zmass)
c-----dress the amplitude with couplings etc.
c-----iph,izz -> switches for photon/z contribution
      iph=one
      izz=one
c-----qq1 -> choose electron/neutrino contributions
      qq1=abs(q1)
c-----LO amplitude
c-----M(+-xx-+)
      do i=1,4
         if (ipsgen .ne. 2) then
         m60hA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+r1*r(qi)*propzQQ)*a60h(1,i))
         endif
         if (ipsgen .ne. 3) then
         m60hB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)*iph+izz*r1*r(qi)*propzLL)*a60h(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m60hC(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQL(3))*a60h(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m60hD(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQL(4))*a60h(4,i))
         endif
      enddo
c-----M(+-xx+-)
      do i=5,8
         if (ipsgen .ne. 2) then
         m60hA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+l1*r(qi)*propzQQ)*a60h(1,i))
         endif
         if (ipsgen .ne. 3) then
         m60hB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+l1*r(qi)*propzLL)*a60h(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m60hC(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQL(3))*a60h(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m60hD(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQL(4))*a60h(4,i))
         endif
      enddo
c-----(-+xx+-)
      do i=9,12
         if (ipsgen .ne. 2) then
         m60hA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+l1*l(qi)*propzQQ)*a60h(1,i))
         endif
         if (ipsgen .ne. 3) then
         m60hB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+l1*l(qi)*propzLL)*a60h(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m60hC(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQL(3))*a60h(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m60hD(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQL(4))*a60h(4,i))
         endif
      enddo
c-----(-+xx-+)
      do i=13,16
         if (ipsgen .ne. 2) then
         m60hA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+r1*l(qi)*propzQQ)*a60h(1,i))
         endif
         if (ipsgen .ne. 3) then
         m60hB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+r1*l(qi)*propzLL)*a60h(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m60hC(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQL(3))*a60h(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m60hD(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQL(4))*a60h(4,i))
         endif
      enddo
c-----Virtual amplitude
c-----M(+-xx-+)
      do i=1,4
         if (ipsgen .ne. 2) then
         m6vhA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+r1*r(qi)*propzQQ)*a6vh(1,i))
         endif
         if (ipsgen .ne. 3) then
         m6vhB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)*iph+izz*r1*r(qi)*propzLL)*a6vh(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m6vhC(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQL(3))*a6vh(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m6vhD(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQL(4))*a6vh(4,i))
         endif
      enddo
c-----M(+-xx+-)
      do i=5,8
         if (ipsgen .ne. 2) then
         m6vhA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+l1*r(qi)*propzQQ)*a6vh(1,i))
         endif
         if (ipsgen .ne. 3) then
         m6vhB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+l1*r(qi)*propzLL)*a6vh(2,i))
          endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
        m6vhC(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQL(3))*a6vh(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m6vhD(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQL(4))*a6vh(4,i))
         endif
      enddo
c-----(-+xx+-)
      do i=9,12
         if (ipsgen .ne. 2) then
         m6vhA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+l1*l(qi)*propzQQ)*a6vh(1,i))
         endif
         if (ipsgen .ne. 3) then
         m6vhB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+l1*l(qi)*propzLL)*a6vh(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m6vhC(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQL(3))*a6vh(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m6vhD(i)=four*esq**2*
     &   (-qq1*Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQL(4))*a6vh(4,i))
         endif
      enddo
c-----(-+xx-+)
      do i=13,16
         if (ipsgen .ne. 2) then
         m6vhA(i)=four*esq**2*
     &   (+Q(qi)**2*(q1*Q(qi)+r1*l(qi)*propzQQ)*a6vh(1,i))
         endif
         if (ipsgen .ne. 3) then
         m6vhB(i)=four*esq**2*
     &   (+qq1*(q1*Q(qi)+r1*l(qi)*propzLL)*a6vh(2,i))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m6vhC(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQL(3))*a6vh(3,i))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m6vhD(i)=four*esq**2*
     &   (+qq1*Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQL(4))*a6vh(4,i))
         endif
      enddo
c-----square the helicity amplitudes
      do i=1,16
         if (ipsgen == 1) then
         m6sqvhAA(i)=two*real(conjg(m60hA(i))*m6vhA(i))
         m6sqvhAB(i)=two*real( conjg(m60hA(i))*m6vhB(i)
     &                        +conjg(m60hB(i))*m6vhA(i) )
         endif
         if (ipsgen == 2) then
         m6sqvhBB(i)=two*real(conjg(m60hB(i))*m6vhB(i))
         m6sqvhBC(i)=two*real( conjg(m60hB(i))*m6vhC(i)
     &                        +conjg(m60hC(i))*m6vhB(i) )
         endif
         if (ipsgen == 3) then
         m6sqvhCC(i)=two*real(conjg(m60hC(i))*m6vhC(i))
         m6sqvhAC(i)=two*real( conjg(m60hA(i))*m6vhC(i)
     &                        +conjg(m60hC(i))*m6vhA(i) )
         m6sqvhCD(i)=two*real( conjg(m60hC(i))*m6vhD(i)
     &                        +conjg(m60hD(i))*m6vhC(i) )
         endif
         if (ipsgen == 4) then
         m6sqvhDD(i)=two*real(conjg(m60hD(i))*m6vhD(i))
         m6sqvhAD(i)=two*real( conjg(m60hA(i))*m6vhD(i)
     &                        +conjg(m60hD(i))*m6vhA(i) )
         m6sqvhBD(i)=two*real( conjg(m60hB(i))*m6vhD(i)
     &                        +conjg(m60hD(i))*m6vhB(i) )
         endif
      enddo
      do i=1,16
c--------commented pieces below only for cross checking with ols PSgen
c         m6sqvh(i)= m6sqvhAA(i)+m6sqvhBB(i)+m6sqvhCC(i)+m6sqvhDD(i)
c     &             +m6sqvhAB(i)+m6sqvhAC(i)+m6sqvhAD(i)+m6sqvhBC(i)
c     &             +m6sqvhBD(i)+m6sqvhCD(i)
         if (ipsgen==1) then
             m6sqvh(i)= m6sqvhAA(i)+m6sqvhAB(i) !AA+AB
         elseif (ipsgen==2) then
             m6sqvh(i)= m6sqvhBB(i)+m6sqvhBC(i) !BB+BC
         elseif (ipsgen==3) then
             m6sqvh(i)= m6sqvhCC(i)+m6sqvhAC(i)+m6sqvhCD(i) !CC+AC+CD
         elseif (ipsgen==4) then
             m6sqvh(i)= m6sqvhDD(i)+m6sqvhAD(i)+m6sqvhBD(i) !D.e+_dpA.e+_dpBD
         endif
      enddo
c-----sum them up
      msqv = zip
      do i=1,16
         msqv=msqv+m6sqvh(i)
      enddo
c-----multiply by color factors,identical particle factors
c-----average over colors,spins
      msqv=msqv*half*aveqq*(Nc-one/Nc)*three*gsq/(sixteen*pi**2)
c-----donehere
      return
      end

