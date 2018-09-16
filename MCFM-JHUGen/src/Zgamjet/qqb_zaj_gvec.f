      subroutine qqb_zaj_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
C*********************************************************************** 
c     Matrix element for Z+gamma+jet production                        *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p6)              *
c     0 -> q(p1)+qbar(p2)+l(p3)+lb(p4)+gam(p5)+glu(p6)                 *
C*********************************************************************** 
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_decl.f'
C--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4)
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      complex(dp):: aqqbn(2,8),aqbqn(2,8),aqgn(2,8),aqbgn(2,8)
      complex(dp):: agqn(2,8),agqbn(2,8)
      real(dp):: qqbn(2),qbqn(2),qgn(2),qbgn(2),gqn(2),gqbn(2)
      integer:: in,i,j,k
      integer, parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer, parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c-----initialized matelem2 for each subprocess
      do j=1,2
         qqbn(j)=zip
         qbqn(j)=zip
         qgn(j) =zip
         qbgn(j)=zip
         gqn(j) =zip
         gqbn(j)=zip
      enddo
c-----evaluate spinor products
      call spinoru(6,p,za,zb)
C-----zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector 
c-----conventions of Bern, Dixon, Kosower, Weinzierl, 
c-----ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)
c-----call color ordered amplitudes,then
c-----call squared matrix element
      if (in==6) then
         call zajn_a60h(1,2,3,4,5,6,p,n,za,zb,zab,zba,aqbqn) !qbq
         call zajn_a60h(2,1,3,4,5,6,p,n,za,zb,zab,zba,aqqbn) !qqb
         do j=1,2
            call zajn_m60sq(j,aqqbn,qqbn(j))
            call zajn_m60sq(j,aqbqn,qbqn(j))
         enddo
      elseif (in==2) then
         call zajn_a60h(1,6,3,4,5,2,p,n,za,zb,zab,zba,aqbgn)  !qbg
         call zajn_a60h(6,1,3,4,5,2,p,n,za,zb,zab,zba,aqgn) !qg
         do j=1,2
            call zajn_m60sq(j,aqgn,qgn(j))
            call zajn_m60sq(j,aqbgn,qbgn(j))
         enddo
      elseif (in==1) then
         call zajn_a60h(2,6,3,4,5,1,p,n,za,zb,zab,zba,agqbn) !gqb
         call zajn_a60h(6,2,3,4,5,1,p,n,za,zb,zab,zba,agqn)  !gq
         do j=1,2
            call zajn_m60sq(j,agqn,gqn(j))
            call zajn_m60sq(j,agqbn,gqbn(j))
         enddo
      endif
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
            msq(j,k)=aveqg*gqbn(-kk(k))
          elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=aveqg*gqn(kk(k))
          elseif ((j > 0) .and. (k == -j)) then
            msq(j,k)=aveqq*qqbn(jj(j))
          elseif ((j < 0) .and. (k == -j)) then
            msq(j,k)=aveqq*qbqn(kk(k))
          elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=aveqg*qgn(jj(j))
          elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=aveqg*qbgn(-jj(j))
          else
            msq(j,k)=zip
          endif
      enddo
      enddo
c-----donehere
      return
      end


      subroutine zajn_m60sq(qi,a6nh,msqn)
      implicit none
      include 'types.f'
***********************************************
* squared matrix element
* qqb_zag, for initial state flavor qi
* gluon line is contracted with a vector n(mu)
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
      complex(dp):: a6nh(2,8),m6nhA(8),m6nhB(8)
      complex(dp):: propzQ,propzL
      real(dp):: m6sqnhAA(16),m6sqnhBB(16),m6sqnhAB(16)
      real(dp):: m6sqnh(16),msqn
      real(dp):: ee,gg,qq,isgn
      integer:: qi,i,j,k
c-----electric and strong coupling
      ee=sqrt(esq)
      gg=sqrt(gsq)
c-----Z propagator
      qq=abs(q1)
c-----QQ propagator
      propzQ=s(3,4)/cplx2(s(3,4)-zmass**2,zwidth*zmass)
c-----QL propagator (2,3)
      propzL=t(3,4,5)/cplx2(t(3,4,5)-zmass**2,zwidth*zmass)
c-----dress the amplitude with couplings etc.
c-----M(+-xx-+)
      do i=1,2
         m6nhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+r1*r(qi)*propzL)*a6nh(2,i) )
      enddo
c-----M(+-xx+-)
      do i=3,4
         m6nhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+l1*r(qi)*propzL)*a6nh(2,i) )
      enddo
c-----(-+xx+-)
      do i=5,6
         m6nhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+l1*l(qi)*propzL)*a6nh(2,i) )
      enddo
c-----(-+xx-+)
      do i=7,8
         m6nhA(i)=twort2*ee**3*gg*
     &   ( + Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     &   ( - qq*(q1*Q(qi)+r1*l(qi)*propzL)*a6nh(2,i) )
      enddo
c-----square the helicity amplitudes
      do i=1,8
         m6sqnhAA(i)=real(m6nhA(i)*conjg(m6nhA(i)))
         m6sqnhBB(i)=real(m6nhB(i)*conjg(m6nhB(i)))
         m6sqnhAB(i)=two*real(m6nhA(i)*conjg(m6nhB(i)))
      enddo
      do i=1,8
c         m6sqnh(i)=m6sqnhAA(i)+m6sqnhBB(i)+m6sqnhAB(i)
         if     (ipsgen == 1) then
            m6sqnh(i)=m6sqnhAA(i)
         elseif (ipsgen == 2) then
            m6sqnh(i)=m6sqnhBB(i)+m6sqnhAB(i)
         else
            write(6,*) 'Parameter ipsgen should be 1 or 2'
            write(6,*) 'ipsgen = ',ipsgen
            stop
         endif
c         if (new_pspace) then
c            m6sqnh(i)=m6sqnhBB(i)+m6sqnhAB(i)
c         else
c            m6sqnh(i)=m6sqnhAA(i)
c         endif
      enddo
c-----sum them up
      msqn = zip
      do i=1,8
         msqn=msqn+m6sqnh(i)
      enddo
c-----multiply by color factors
c-----no average over spins,no average over colors
      msqn=msqn
c-----donehere
      return
      end

