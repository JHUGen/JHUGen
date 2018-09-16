
      subroutine qqb_zaj_v(p,msqv)
      implicit none
      include 'types.f'
****************************************************************
*   Virtual matrix element for
*   0 -> q(-p1) + qb(-p4) + g(p2) + a(p3) + lb(p5) + l(p6)
*   l = charged lepton/neutrino
****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'scale.f'
      include 'epinv2.f'
      include 'epinv.f'
      include 'nflav.f'
      real(dp):: p(mxpart,4),msqv(-nf:nf,-nf:nf)
      real(dp):: msq0(-nf:nf,-nf:nf),subuv
      complex(dp):: a60h(4,16),a6vh(5,16)
      real(dp):: pbdk(mxpart,4)
      real(dp):: qqb(2),qbq(2),qg(2),qbg(2),gq(2),gqb(2)
      integer:: i,j,k
      integer,parameter:: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter:: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
c-----set scheme
      scheme='dred'
c-----calculate lowest order for UV renormalization
      call qqb_zaj(p,msq0)
c-----UV counterterm contains the finite renormalization to arrive
c-----at MS bar scheme. 
      subuv=ason2pi*xn
     & *(epinv*(11._dp-2._dp*real(nflav,dp)/xn)-1._dp)/6._dp
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
      do j=1,2
         qqb(j)=zip
         qbq(j)=zip
         qg(j) =zip
         qbg(j)=zip
         gq(j) =zip
         gqb(j)=zip
      enddo
c-----compute color ordered amplitudes
c-----then compute squared matrix element
c-----qqbar subprocess
      call zaj_a60h(1,2,3,4,5,6,za,zb,a60h)
      call zaj_a6vh(1,2,3,4,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,qqb(2))
      call zaj_m6vsq(1,a60h,a6vh,qqb(1))
c-----qbarq subprocess
      call zaj_a60h(4,2,3,1,5,6,za,zb,a60h)
      call zaj_a6vh(4,2,3,1,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,qbq(2))
      call zaj_m6vsq(1,a60h,a6vh,qbq(1))
c-----qg subprocess
      call zaj_a60h(2,1,3,4,5,6,za,zb,a60h)
      call zaj_a6vh(2,1,3,4,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,qg(2))
      call zaj_m6vsq(1,a60h,a6vh,qg(1))
c-----qbg subprocess
      call zaj_a60h(4,1,3,2,5,6,za,zb,a60h)
      call zaj_a6vh(4,1,3,2,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,qbg(2))
      call zaj_m6vsq(1,a60h,a6vh,qbg(1))
c-----gq subprocess
      call zaj_a60h(2,4,3,1,5,6,za,zb,a60h)
      call zaj_a6vh(2,4,3,1,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,gq(2))
      call zaj_m6vsq(1,a60h,a6vh,gq(1))
c-----gqb subprocess
      call zaj_a60h(1,4,3,2,5,6,za,zb,a60h)
      call zaj_a6vh(1,4,3,2,5,6,za,zb,a6vh)
      call zaj_m6vsq(2,a60h,a6vh,gqb(2))
      call zaj_m6vsq(1,a60h,a6vh,gqb(1))
c-----initialize msq
      do i=-nf,nf
      do j=-nf,nf
         msqv(i,j)=zip
      enddo
      enddo
c-----fill msq
      do j=-nf,nf
      do k=-nf,nf
          if ((j == 0) .and. (k == 0)) then
            msqv(j,k)=0._dp
          elseif ((j == 0) .and. (k < 0)) then
            msqv(j,k)=aveqg*gqb(-kk(k))-subuv*msq0(j,k)
          elseif ((j == 0) .and. (k > 0)) then
            msqv(j,k)=aveqg*gq(kk(k))-subuv*msq0(j,k)
          elseif ((j > 0) .and. (k == -j)) then
            msqv(j,k)=aveqq*qqb(jj(j))-subuv*msq0(j,k)
          elseif ((j < 0) .and. (k == -j)) then
            msqv(j,k)=aveqq*qbq(kk(k))-subuv*msq0(j,k)
          elseif ((j > 0) .and. (k == 0)) then
            msqv(j,k)=aveqg*qg(jj(j))-subuv*msq0(j,k)
          elseif ((j < 0) .and. (k == 0)) then
            msqv(j,k)=aveqg*qbg(-jj(j))-subuv*msq0(j,k)
          else
            msqv(j,k)=0._dp
          endif
      enddo
      enddo
c-----donehere
      return
      end

      subroutine zaj_m6vsq(qi,a60h,a6vh,msqv)
      implicit none
      include 'types.f'
***********************************************
* virtual matrix element squared
* qqb_zaj, for initial state flavor qi
***********************************************
      
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
      include 'new_pspace.f'
      include 'ipsgen.f'
      real(dp):: t
      integer:: j1,j2,j3,j4,j5,j6
      complex(dp):: a60h(4,16),m60hA(16),m60hB(16)
      complex(dp):: a6vh(5,16),m6vhA(16),m6vhB(16)
      complex(dp):: propzQQ,propzLL,propzQL(3:4)
      real(dp):: m6sqvhAA(16),m6sqvhBB(16),m6sqvhAB(16)
      real(dp):: m6sqvh(16),msqv
      real(dp):: iwdth,iph,izz,qq1,ee,gg
      real(dp):: SumQiNf,SumViNf,sw,cw
      integer:: qi,i,j,k
c-----electric and strong coupling
      ee=sqrt(esq)
      gg=sqrt(gsq)
      sw=sqrt(xw)
      cw=sqrt(one-xw)
c-----Z propagator
      iwdth = one
c-----QQ propagator
      propzQQ=s(5,6)/cplx2(s(5,6)-zmass**2,iwdth*zwidth*zmass)
c-----QL propagator (2,3)
      propzQL(3)=t(1,4,2)/cplx2(t(1,4,2)-zmass**2,iwdth*zwidth*zmass)
c-----dress the amplitude with couplings etc.
c-----iph,izz -> switches for photon/z contribution
      iph=one
      izz=one
c-----qq1 -> choose electron/neutrino contributions
      qq1=abs(q1)
c-----LO amplitude
c-----M(+-xx-+)
      do i=1,4
         m60hA(i)=twort2*ee**3*gg*
     &   (+Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQQ)*a60h(1,i) )
         m60hB(i)=twort2*ee**3*gg*
     &   (-qq1*(q1*Q(qi)+r1*r(qi)*propzQL(3))*a60h(3,i) )
      enddo
c-----M(+-xx+-)
      do i=5,8
         m60hA(i)=twort2*ee**3*gg*
     &   (+Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQQ)*a60h(1,i) )
         m60hB(i)=twort2*ee**3*gg*
     &   (+qq1*(q1*Q(qi)+l1*r(qi)*propzQL(3))*a60h(3,i) )
      enddo
c-----(-+xx+-)
      do i=9,12
         m60hA(i)=twort2*ee**3*gg*
     &   (+Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQQ)*a60h(1,i)) 
         m60hB(i)=twort2*ee**3*gg*
     &   (-qq1*(q1*Q(qi)+l1*l(qi)*propzQL(3))*a60h(3,i)) 
      enddo
c-----(-+xx-+)
      do i=13,16
         m60hA(i)=twort2*ee**3*gg*
     &   (+Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQQ)*a60h(1,i)) 
         m60hB(i)=twort2*ee**3*gg*
     &   (+qq1*(q1*Q(qi)+r1*l(qi)*propzQL(3))*a60h(3,i)) 
      enddo
c-----Virtual amplitude
c-----Sum of nf light quark electric charge
      SumQiNf=zip
      do i=1,nf
         SumQiNf=SumQiNf+Q(i)
      enddo
c-----Sum of nf light quark z-coupling: vL(i)+vR(i)
      SumViNf=zip
      do i=1,nf
         SumViNf=SumViNf + l(i)+r(i)
      enddo
c-----M(+-xx-+)
      do i=1,4
         m6vhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*( 
     &             +(q1*Q(qi)+r1*r(qi)*propzQQ)*a6vh(1,i)
     &             +(q1*SumQiNf + half*r1*SumViNf*propzQQ)*a6vh(2,i)
     &             +(r1/two/sw/cw)*propzQQ*a6vh(3,i) 
     &             ) )
         m6vhB(i)=twort2*ee**3*gg*
     &   ( - qq1*( 
     &           +(q1*Q(qi)+r1*r(qi)*propzQL(3))*a6vh(4,i)
     &           +(r1/two/sw/cw)*propzQL(3)*a6vh(5,i) 
     &           ))
      enddo
c-----M(+-xx+-)
      do i=5,8
         m6vhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(
     &             +(q1*Q(qi)+l1*r(qi)*propzQQ)*a6vh(1,i)
     &             +(q1*SumQiNf + half*l1*SumViNf*propzQQ)*a6vh(2,i)
     &             +(l1/two/sw/cw)*propzQQ*a6vh(3,i) 
     &             ))
         m6vhB(i)=twort2*ee**3*gg*
     &   ( + qq1*(
     &           +(q1*Q(qi)+l1*r(qi)*propzQL(3))*a6vh(4,i)
     &           +(l1/two/sw/cw)*propzQL(3)*a6vh(5,i) 
     &           ))
      enddo
c-----(-+xx+-)
      do i=9,12
         m6vhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(
     &             +(q1*Q(qi)+l1*l(qi)*propzQQ)*a6vh(1,i)
     &             +(q1*SumQiNf + half*l1*SumViNf*propzQQ)*a6vh(2,i)
     &             +(l1/two/sw/cw)*propzQQ*a6vh(3,i)
     &             ))
         m6vhB(i)=twort2*ee**3*gg*
     &   ( - qq1*(
     &           +(q1*Q(qi)+l1*l(qi)*propzQL(3))*a6vh(4,i)
     &           +(l1/two/sw/cw)*propzQL(3)*a6vh(5,i)
     &           ))
      enddo
c-----(-+xx-+)
      do i=13,16
         m6vhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(
     &             +(q1*Q(qi)+r1*l(qi)*propzQQ)*a6vh(1,i)
     &             +(q1*SumQiNf + half*r1*SumViNf*propzQQ)*a6vh(2,i)
     &             +(r1/two/sw/cw)*propzQQ*a6vh(3,i)
     &             ))
         m6vhB(i)=twort2*ee**3*gg*
     &   ( + qq1*(
     &           +(q1*Q(qi)+r1*l(qi)*propzQL(3))*a6vh(4,i)
     &           +(r1/two/sw/cw)*propzQL(3)*a6vh(5,i)
     &           ))
      enddo
c-----square the helicity amplitudes
      do i=1,16
         m6sqvhAA(i)=two*real(conjg(m60hA(i))*m6vhA(i))
         m6sqvhBB(i)=two*real(conjg(m60hB(i))*m6vhB(i))
         m6sqvhAB(i)= two*real(conjg(m60hA(i))*m6vhB(i))
     &               +two*real(conjg(m60hB(i))*m6vhA(i))
      enddo
      do i=1,16
c         m6sqvh(i)=m6sqvhAA(i)+m6sqvhBB(i)+m6sqvhAB(i)
         if     (ipsgen == 1) then
            m6sqvh(i)=m6sqvhAA(i)
         elseif (ipsgen == 2) then
            m6sqvh(i)=m6sqvhBB(i)+m6sqvhAB(i)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
c         if (new_pspace) then
c            m6sqvh(i)=m6sqvhBB(i)+m6sqvhAB(i)
c         else
c            m6sqvh(i)=m6sqvhAA(i)
c         endif
      enddo
c-----sum them up
      msqv = zip
      do i=1,16
         msqv=msqv+m6sqvh(i)
      enddo
c-----multiply by color factors
c-----average over spins,no average over colors
      msqv=8._dp*msqv*gsq/16._dp/pi**2
c-----donehere
      return
      end

