****************************************************************
*  Author: H. Hartanto                                         *
*                                                              *
*   0 -> q(-p1) + qb(-p4) + g(p2) + a(p3) + lb(p5) + l(p6)     *
*                                                              *
****************************************************************
      subroutine qqb_zaj(p,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: pbdk(mxpart,4)
      complex(dp):: aqqb(4,16),aqbq(4,16),aqg(4,16)
      complex(dp):: aqbg(4,16),agq(4,16),agqb(4,16)
      real(dp):: qqb(2),qbq(2),qg(2),qbg(2),gq(2),gqb(2)
      integer:: i,j,k
      integer, parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer, parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
c-----swap momenta
c-----implemented (ala bdk)
c----- 0 -> q(p1) + qb(p4) + g(p2) + gam(p3) + lb(p5) + l(p6)
c-----mcfm input:
c----- 0 -> q(p2) + qb(p1) + l(p3) + lb(p4) + gam(p5) + g(p6)
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
c-----call color ordered amplitudes
      call zaj_a60h(1,2,3,4,5,6,za,zb,aqqb)
      call zaj_a60h(4,2,3,1,5,6,za,zb,aqbq)
      call zaj_a60h(2,1,3,4,5,6,za,zb,aqg)
      call zaj_a60h(4,1,3,2,5,6,za,zb,aqbg)
      call zaj_a60h(2,4,3,1,5,6,za,zb,agq)
      call zaj_a60h(1,4,3,2,5,6,za,zb,agqb)
c-----initialized matelem2 for each subprocess
      do j=1,2
         qqb(j)=zip
         qbq(j)=zip
         qg(j) =zip
         qbg(j)=zip
         gq(j) =zip
         gqb(j)=zip
      enddo
c-----call squared matrix element
      do j=1,2
         call zaj_m60sq(j,aqqb,qqb(j))
         call zaj_m60sq(j,aqbq,qbq(j))
         call zaj_m60sq(j,aqg,qg(j))
         call zaj_m60sq(j,aqbg,qbg(j))
         call zaj_m60sq(j,agq,gq(j))
         call zaj_m60sq(j,agqb,gqb(j))
      enddo
c-----initialize msq
      do i=-nf,nf
      do j=-nf,nf
         msq(i,j)=zip
      enddo
      enddo
c-----fill msq
      do j=-nf,nf
      do k=-nf,nf
          if ((j == 0) .and. (k == 0)) then
            msq(j,k)=zip
          elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=aveqg*gqb(-kk(k))
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=aveqg*gq(kk(k))
          elseif ((j > 0) .and. (k == -j)) then
            msq(j,k)=aveqq*qqb(jj(j))
          elseif ((j < 0) .and. (k == -j)) then
            msq(j,k)=aveqq*qbq(kk(k))
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=aveqg*qg(jj(j))
          elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=aveqg*qbg(-jj(j))
          else
            msq(j,k)=zip
          endif
      enddo
      enddo
c-----donehere
      return
      end


      subroutine zaj_m60sq(qi,a60h,msq)
      implicit none
      include 'types.f'
***********************************************
* squared matrix element
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
      complex(dp):: propzQQ,propzLL,propzQL(3:4)
      real(dp):: m6sq0h(16),msq
      real(dp):: m6sq0hAA(16),m6sq0hBB(16),m6sq0hAB(16)
      real(dp):: ee,gg,iwdth,qq
      integer:: qi,i,j,k
c-----electric and strong coupling
      ee=sqrt(esq)
      gg=sqrt(gsq)
c-----Z propagator
      iwdth = one
      qq=abs(q1)
c-----QQ propagator
      propzQQ=s(5,6)/cplx2(s(5,6)-zmass**2,iwdth*zwidth*zmass)
c-----QL propagator (2,3)
      propzQL(3)=t(1,4,2)/cplx2(t(1,4,2)-zmass**2,iwdth*zwidth*zmass)
c-----dress the amplitude with couplings etc.
c-----M(+-xx-+)
      do i=1,4
         m60hA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQQ)*a60h(1,i) )
         m60hB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+r1*r(qi)*propzQL(3))*a60h(3,i) )
      enddo
c-----M(+-xx+-)
      do i=5,8
         m60hA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQQ)*a60h(1,i) )
         m60hB(i)=twort2*ee**3*gg*
     &   ( + qq*(q1*Q(qi)+l1*r(qi)*propzQL(3))*a60h(3,i) )
      enddo
c-----(-+xx+-)
      do i=9,12
         m60hA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQQ)*conjg(a60h(1,i-8)) )
         m60hB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+l1*l(qi)*propzQL(3))*conjg(a60h(3,i-8)) )
      enddo
c-----(-+xx-+)
      do i=13,16
         m60hA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQQ)*conjg(a60h(1,i-8)) )
         m60hB(i)=twort2*ee**3*gg*
     &   ( + qq*(q1*Q(qi)+r1*l(qi)*propzQL(3))*conjg(a60h(3,i-8)) )
      enddo

      if (.false.) then
      write(*,*) ' '
      do i=1,8
         write(*,*) 'A6q(',i,') = ',a60h(1,i)
      enddo
      do i=9,16
         write(*,*) 'A6q(',i,') = ',conjg(a60h(1,i-8))
         write(*,*) '                   = ',a60h(1,i)
      enddo
      write(*,*) ' '
      do i=1,8
         write(*,*) 'A6l(',i,') = ',a60h(3,i)
      enddo
      do i=9,16
         write(*,*) 'A6l(',i,') = ',conjg(a60h(3,i-8))
         write(*,*) '                   = ',a60h(3,i)
      enddo
      pause
      endif

c-----square the helicity amplitudes
      do i=1,16
         m6sq0hAA(i)=real(m60hA(i)*conjg(m60hA(i)))
         m6sq0hBB(i)=real(m60hB(i)*conjg(m60hB(i)))
         m6sq0hAB(i)=2._dp*real(m60hA(i)*conjg(m60hB(i)))
      enddo
      do i=1,16
c         m6sq0h(i)=m6sq0hAA(i)+m6sq0hBB(i)+m6sq0hAB(i)
c         m6sq0h(i)=m6sq0hBB(i)

c--- ipsgen may now be equal to 3 or 4 when this routine is
c--- called by fragmentation dipoles in qqb_zaa_gs.f
         if     (ipsgen == 1) then
            m6sq0h(i)=m6sq0hAA(i)
         elseif (ipsgen > 1) then
            m6sq0h(i)=m6sq0hBB(i)+m6sq0hAB(i)
         endif

c         if     (ipsgen == 1) then
c            m6sq0h(i)=m6sq0hAA(i)
c         elseif (ipsgen == 2) then
c            m6sq0h(i)=m6sq0hBB(i)+m6sq0hAB(i)
c         else
c         write(6,*) 'Parameter ipsgen should be 1 or 2'
c         write(6,*) 'ipsgen = ',ipsgen
c         stop
c         endif

c         if (new_pspace) then
c            m6sq0h(i)=m6sq0hBB(i)+m6sq0hAB(i)
c         else
c            m6sq0h(i)=m6sq0hAA(i)
c         endif
      enddo
c-----sum them up
      msq = zip
      do i=1,16
         msq=msq+m6sq0h(i)
      enddo
c-----multiply by color factors
c-----average over spins,no average over colors
      msq=8._dp*msq
c-----donehere
      return
      end

