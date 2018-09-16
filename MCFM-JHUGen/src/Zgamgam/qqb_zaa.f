****************************************************************
*   0 -> q(-p1) + qb(-p4) + a(p2) + a(p3) + lb(p5) + l(p6)
*   l = charged lepton/neutrino
****************************************************************
      subroutine qqb_zaa(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: a60h(4,16)
      real(dp):: qqb(2),qbq(2),pbdk(mxpart,4)
      integer:: i,j,k
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
      call zaa_m60sq(2,a60h,qqb(2))
      call zaa_m60sq(1,a60h,qqb(1))
c-----qbarq subprocess
      call zaa_a60h(4,2,3,1,5,6,za,zb,a60h)
      call zaa_m60sq(2,a60h,qbq(2))
      call zaa_m60sq(1,a60h,qbq(1))
c-----initialize msq
      do i=-nf,nf
         do j=-nf,nf
            msq(i,j)=zip
         enddo
      enddo
c-----fill msq
      msq(1,-1)=qqb(1)
      msq(3,-3)=qqb(1)
      msq(5,-5)=qqb(1)
      msq(2,-2)=qqb(2)
      msq(4,-4)=qqb(2)
      msq(-1,1)=qbq(1)
      msq(-3,3)=qbq(1)
      msq(-5,5)=qbq(1)
      msq(-2,2)=qbq(2)
      msq(-4,4)=qbq(2)
c-----donehere
      return
      end
***********************************************
* squared matrix element
* qqb_zaa, for initial state flavor qi
***********************************************
      subroutine zaa_m60sq(qi,a60h,msq)
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
      complex(dp)::a60h(4,16),m60hA(16),m60hB(16),m60hC(16),
     & m60hD(16),propzQQ,propzLL,propzQL(3:4)
      real(dp):: m6sq0h(16),msq
      real(dp):: m6sq0hAA(16),m6sq0hBB(16),m6sq0hCC(16),
     .m6sq0hDD(16),m6sq0hAB(16),m6sq0hAC(16),m6sq0hAD(16),m6sq0hBC(16),
     .m6sq0hBD(16),m6sq0hCD(16)
      real(dp):: iwdth,iph,izz,qq1
      integer:: qi,i,j,k
c-----Z propagator
      iwdth = 1._dp
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
c-----square the helicity amplitudes
      do i=1,16
         if (ipsgen == 1) then
         m6sq0hAA(i)=abs(m60hA(i))**2
         m6sq0hAB(i)=2._dp*real(m60hA(i)*conjg(m60hB(i)))
         endif
         if (ipsgen == 2) then
         m6sq0hBB(i)=abs(m60hB(i))**2
         m6sq0hBC(i)=2._dp*real(m60hB(i)*conjg(m60hC(i)))
         endif
         if (ipsgen == 3) then
         m6sq0hCC(i)=abs(m60hC(i))**2
         m6sq0hAC(i)=2._dp*real(m60hA(i)*conjg(m60hC(i)))
         m6sq0hCD(i)=2._dp*real(m60hC(i)*conjg(m60hD(i)))
         endif
         if (ipsgen == 4) then
         m6sq0hDD(i)=abs(m60hD(i))**2
         m6sq0hAD(i)=2._dp*real(m60hA(i)*conjg(m60hD(i)))
         m6sq0hBD(i)=2._dp*real(m60hB(i)*conjg(m60hD(i)))
         endif
      enddo
      do i=1,16
c---------commented pieces below only for crosscheck with old PSgen
c         m6sq0h(i)= m6sq0hAA(i)+m6sq0hBB(i)+m6sq0hCC(i)+m6sq0hDD(i)
c     &             +m6sq0hAB(i)+m6sq0hAC(i)+m6sq0hAD(i)+m6sq0hBC(i)
c     &             +m6sq0hBD(i)+m6sq0hCD(i)
c         m6sq0h(i)= m6sq0hAB(i)+m6sq0hAC(i)+m6sq0hAD(i)+m6sq0hBC(i)
c     &             +m6sq0hBD(i)+m6sq0hCD(i)
         if (ipsgen==1) then
             m6sq0h(i)= m6sq0hAA(i)+m6sq0hAB(i) !AA+AB
         elseif (ipsgen==2) then
             m6sq0h(i)= m6sq0hBB(i)+m6sq0hBC(i) !BB+BC
         elseif (ipsgen==3) then
             m6sq0h(i)= m6sq0hCC(i)+m6sq0hAC(i)+m6sq0hCD(i) !CC+AC+CD
         elseif (ipsgen==4) then
             m6sq0h(i)= m6sq0hDD(i)+m6sq0hAD(i)+m6sq0hBD(i) !D.e+_dpA.e+_dpBD
         endif
      enddo
c-----sum them up
      msq = zip
      do i=1,16
         msq=msq+m6sq0h(i)
      enddo
c-----multiply by color factors,identical particle factors
c-----average over colors,spins
      msq=half*spinave*(one/three)*msq
c-----donehere
      return
      end
