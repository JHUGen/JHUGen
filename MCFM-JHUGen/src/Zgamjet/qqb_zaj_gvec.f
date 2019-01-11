      subroutine qqb_zaj_gvec(p,n,in,msq)
C*********************************************************************** 
c     Matrix element for Z+gamma+jet production                        *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p6)              *
c     0 -> q(p1)+qbar(p2)+l(p3)+lb(p4)+gam(p5)+glu(p6)                 *
C*********************************************************************** 
      implicit none
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),n(4)
      double complex zab(mxpart,mxpart),zba(mxpart,mxpart)
      double complex aqqbn(2,8),aqbqn(2,8),aqgn(2,8),aqbgn(2,8)
      double complex agqn(2,8),agqbn(2,8)
      double precision qqbn(2),qbqn(2),qgn(2),qbgn(2),gqn(2),gqbn(2)
      integer in,i,j,k
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
      if (in.eq.6) then
         call zajn_a60h(1,2,3,4,5,6,p,n,za,zb,zab,zba,aqbqn) !qbq
         call zajn_a60h(2,1,3,4,5,6,p,n,za,zb,zab,zba,aqqbn) !qqb
         do j=1,2
            call zajn_m60sq(j,aqqbn,qqbn(j))
            call zajn_m60sq(j,aqbqn,qbqn(j))
         enddo
      elseif (in.eq.2) then
         call zajn_a60h(1,6,3,4,5,2,p,n,za,zb,zab,zba,aqbgn)  !qbg
         call zajn_a60h(6,1,3,4,5,2,p,n,za,zb,zab,zba,aqgn) !qg
         do j=1,2
            call zajn_m60sq(j,aqgn,qgn(j))
            call zajn_m60sq(j,aqbgn,qbgn(j))
         enddo
      elseif (in.eq.1) then
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
          if ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=0d0
          elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=aveqg*gqbn(-kk(k))
          elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=aveqg*gqn(kk(k))
          elseif ((j .gt. 0) .and. (k .eq. -j)) then
            msq(j,k)=aveqq*qqbn(jj(j))
          elseif ((j .lt. 0) .and. (k .eq. -j)) then
            msq(j,k)=aveqq*qbqn(kk(k))
          elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=aveqg*qgn(jj(j))
          elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=aveqg*qbgn(-jj(j))
          else
            msq(j,k)=0d0
          endif
      enddo
      enddo
c-----donehere
      return
      end


      subroutine zajn_m60sq(qi,a6nh,msqn)
***********************************************
* squared matrix element
* qqb_zag, for initial state flavor qi
* gluon line is contracted with a vector n(mu)
***********************************************
      implicit none
      include 'constants.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'new_pspace.f'
      include 'ipsgen.f'
      double precision t
      integer j1,j2,j3,j4,j5,j6
      double complex a6nh(2,8),m6nhA(8),m6nhB(8)
      double complex propzQ,propzL
      double precision m6sqnhAA(16),m6sqnhBB(16),m6sqnhAB(16)
      double precision m6sqnh(16),msqn
      double precision ee,gg,qq,isgn
      integer qi,i,j,k
c-----electric and strong coupling
      ee=dsqrt(esq)
      gg=dsqrt(gsq)
c-----Z propagator
      qq=dabs(q1)
c-----QQ propagator
      propzQ=s(3,4)/Dcmplx(s(3,4)-zmass**2,zwidth*zmass)
c-----QL propagator (2,3)
      propzL=t(3,4,5)/Dcmplx(t(3,4,5)-zmass**2,zwidth*zmass)
c-----dress the amplitude with couplings etc.
c-----M(+-xx-+)
      do i=1,2
         m6nhA(i)=twort2*ee**3*gg*
     .   ( + Q(qi)*(q1*Q(qi)+r1*r(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     .   ( - qq*(q1*Q(qi)+r1*r(qi)*propzL)*a6nh(2,i) )
      enddo
c-----M(+-xx+-)
      do i=3,4
         m6nhA(i)=twort2*ee**3*gg*
     .   ( + Q(qi)*(q1*Q(qi)+l1*r(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     .   ( - qq*(q1*Q(qi)+l1*r(qi)*propzL)*a6nh(2,i) )
      enddo
c-----(-+xx+-)
      do i=5,6
         m6nhA(i)=twort2*ee**3*gg*
     .   ( + Q(qi)*(q1*Q(qi)+l1*l(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     .   ( - qq*(q1*Q(qi)+l1*l(qi)*propzL)*a6nh(2,i) )
      enddo
c-----(-+xx-+)
      do i=7,8
         m6nhA(i)=twort2*ee**3*gg*
     .   ( + Q(qi)*(q1*Q(qi)+r1*l(qi)*propzQ)*a6nh(1,i) )
         m6nhB(i)=twort2*ee**3*gg*
     .   ( - qq*(q1*Q(qi)+r1*l(qi)*propzL)*a6nh(2,i) )
      enddo
c-----square the helicity amplitudes
      do i=1,8
         m6sqnhAA(i)=dreal(m6nhA(i)*dconjg(m6nhA(i)))
         m6sqnhBB(i)=dreal(m6nhB(i)*dconjg(m6nhB(i)))
         m6sqnhAB(i)=2D0*dreal(m6nhA(i)*dconjg(m6nhB(i)))
      enddo
      do i=1,8
c         m6sqnh(i)=m6sqnhAA(i)+m6sqnhBB(i)+m6sqnhAB(i)
         if     (ipsgen .eq. 1) then
            m6sqnh(i)=m6sqnhAA(i)
         elseif (ipsgen .eq. 2) then
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

