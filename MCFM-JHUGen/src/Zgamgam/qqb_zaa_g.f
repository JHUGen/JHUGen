*******************************************************************
*  0 -> q(-p1) + qb(-p5) + a(p2) + a(p3) + g(p4) + lb(p6) + l(p7) *
*******************************************************************
      subroutine qqb_zaa_g(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: pnt(mxpart,4)
      integer:: i,j,k
      integer:: hq,h2,h3,h4,lh
      complex(dp):: a70h1(2,2,2,2,2),a70h2(2,2,2,2,2)
      complex(dp):: a70h3(2,2,2,2,2),a70h4(2,2,2,2,2)
      real(dp):: qqb(2),qbq(2),qg(2),qbg(2),gq(2),gqb(2)
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c-----convert to Nagy-Trocsanyi momentum convention
      do i=1,4
         pnt(1,i)=p(2,i)
         pnt(2,i)=p(5,i)
         pnt(3,i)=p(6,i)
         pnt(4,i)=p(7,i)
         pnt(5,i)=p(1,i)
         pnt(6,i)=p(4,i)
         pnt(7,i)=p(3,i)
      enddo
c-----initialize matelem
      do j=1,2
         qqb(j)=zip
         qbq(j)=zip
      enddo
c-----calculate sipnor products, invariants
      call spinoru(7,pnt,za,zb)
c-----call color ordered amplitudes and compute msq for qqb
      call zaag_a70h(1,2,3,4,5,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,qqb(j))
      enddo
c-----call color ordered amplitudes and compute msq for qbq
      call zaag_a70h(5,2,3,4,1,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,qbq(j))
      enddo
c-----call color ordered amplitudes and compute msq for qg
      call zaag_a70h(4,2,3,1,5,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,qg(j))
      enddo
c-----call color ordered amplitudes and compute msq for qbg
      call zaag_a70h(5,2,3,1,4,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,qbg(j))
      enddo
c-----call color ordered amplitudes and compute msq for gq
      call zaag_a70h(4,2,3,5,1,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,gq(j))
      enddo
c-----call color ordered amplitudes and compute msq for gqb
      call zaag_a70h(1,2,3,5,4,6,7,a70h1,a70h2,a70h3,a70h4)
      do j=1,2
         call zaag_m70sq(j,a70h1,a70h2,a70h3,a70h4,gqb(j))
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
            msq(j,k)=0._dp
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
            msq(j,k)=0._dp
          endif
      enddo
      enddo
c-----done
      return
      end

*******************************************************************
* 0 -> q(-p1) + qb(-p5) + a(p2) + a(p3) + g(p4) + lb(p6) + l(p7)
* return helicity amplitudes for each channel
********************************************************************
      subroutine zaag_a70h(j1,j2,j3,j4,j5,j6,j7,
     & a70h1,a70h2,a70h3,a70h4)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ipsgen.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      integer:: hq,h2,h3,h4,lh
      complex(dp):: a70h1(2,2,2,2,2),a70h2(2,2,2,2,2)
      complex(dp):: a70h3(2,2,2,2,2),a70h4(2,2,2,2,2)
c-----initialize
      do hq=1,2
      do h2=1,2
      do h3=1,2
      do h4=1,2
      do lh=1,2
      a70h1(hq,h2,h3,h4,lh)=czip
      a70h2(hq,h2,h3,h4,lh)=czip
      a70h3(hq,h2,h3,h4,lh)=czip
      a70h4(hq,h2,h3,h4,lh)=czip
      enddo
      enddo
      enddo
      enddo
      enddo
c----compute amplitude
      if (ipsgen .ne. 2) then
      call xzqqaag_qq(j1,j2,j3,j4,j5,j6,j7,a70h1)
      endif
      if (ipsgen .ne. 3) then
      call xzqqaag_ll(j7,j4,j3,j2,j6,j5,j1,a70h2)
      endif
      if (ipsgen .ne. 1) then
      call xzqqaag_ql(j1,j2,j3,j4,j5,j6,j7,a70h3,a70h4)
      endif
c-----done
      return 
      end

********************************************************************
* 0 -> q(-p1) + qb(-p5) + a(p2) + a(p3) + g(p4) + lb(p6) + l(p7)
* return matrix element squared, for initial flavor qi
* given the helicity amplitudes from each channel
* no averaging
********************************************************************
      subroutine zaag_m70sq(qi,a70h1,a70h2,a70h3,a70h4,msq)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'ipsgen.f'
      real(dp):: msq,msqsum(32),t,qq,gg
      real(dp):: msqAA(32),msqBB(32),msqCC(32),msqDD(32),msqAB(32)
      real(dp):: msqAC(32),msqAD(32),msqBC(32),msqBD(32),msqCD(32)
      integer:: i,j,k,ihel
      integer:: h2,h3,h4,qi
      complex(dp):: a70h1(2,2,2,2,2),a70h2(2,2,2,2,2)
      complex(dp):: a70h3(2,2,2,2,2),a70h4(2,2,2,2,2)
      complex(dp):: propQQ,propLL,propQL3,propQL4
      complex(dp):: m70hA(32),m70hB(32),m70hC(32),m70hD(32)
c-----
      gg=sqrt(gsq)
      qq=abs(q1)
c-----
      propQQ =s(6,7)/cplx2(s(6,7)-zmass**2,zwidth*zmass)
      propQL3=t(2,6,7)/cplx2(t(2,6,7)-zmass**2,zwidth*zmass)
      propQL4=t(3,6,7)/cplx2(t(3,6,7)-zmass**2,zwidth*zmass)
      propLL =t(1,4,5)/cplx2(t(1,4,5)-zmass**2,zwidth*zmass)
c-----      
      ihel=1
      do h2=1,2
      do h3=1,2
      do h4=1,2
         if (ipsgen .ne. 2) then
         m70hA(ihel)=four*esq**2*gg*
     &   (+(q1*Q(qi)+l(qi)*l1*propQQ)*Q(qi)**2*a70h1(1,h2,h3,h4,1))
         endif
         if (ipsgen .ne. 3) then
         m70hB(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+l(qi)*l1*propLL)*a70h2(1,h2,h3,h4,1))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m70hC(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+l(qi)*l1*propQL3)*(-Q(qi))*a70h3(1,h2,h3,h4,1))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m70hD(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+l(qi)*l1*propQL4)*(-Q(qi))*a70h4(1,h2,h3,h4,1))
         endif
         ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
         if (ipsgen .ne. 2) then
         m70hA(ihel)=four*esq**2*gg*
     &   (+(q1*Q(qi)+r(qi)*l1*propQQ)*Q(qi)**2*a70h1(2,h2,h3,h4,1))
         endif
         if (ipsgen .ne. 3) then
         m70hB(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+r(qi)*l1*propLL)*a70h2(2,h2,h3,h4,1))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m70hC(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+r(qi)*l1*propQL3)*(-Q(qi))*a70h3(2,h2,h3,h4,1))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m70hD(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+r(qi)*l1*propQL4)*(-Q(qi))*a70h4(2,h2,h3,h4,1))
         endif
         ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
         if (ipsgen .ne. 2) then
         m70hA(ihel)=four*esq**2*gg*
     &   (+(q1*Q(qi)+r(qi)*r1*propQQ)*Q(qi)**2*a70h1(2,h2,h3,h4,2))
         endif
         if (ipsgen .ne. 3) then
         m70hB(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+r(qi)*r1*propLL)*a70h2(2,h2,h3,h4,2))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m70hC(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+r(qi)*r1*propQL3)*(-Q(qi))*a70h3(2,h2,h3,h4,2))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m70hD(ihel)=four*esq**2*gg*
     &   (+qq*(q1*Q(qi)+r(qi)*r1*propQL4)*(-Q(qi))*a70h4(2,h2,h3,h4,2))
         endif
         ihel=ihel+1
      enddo
      enddo
      enddo
c-----
      do h2=1,2
      do h3=1,2
      do h4=1,2
         if (ipsgen .ne. 2) then
         m70hA(ihel)=four*esq**2*gg*
     &   (+(q1*Q(qi)+l(qi)*r1*propQQ)*Q(qi)**2*a70h1(1,h2,h3,h4,2))
         endif
         if (ipsgen .ne. 3) then
         m70hB(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+l(qi)*r1*propLL)*a70h2(1,h2,h3,h4,2))
         endif
         if ((ipsgen == 2) .or. (ipsgen == 3)) then
         m70hC(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+l(qi)*r1*propQL3)*(-Q(qi))*a70h3(1,h2,h3,h4,2))
         endif
         if ((ipsgen == 3) .or. (ipsgen == 4)) then
         m70hD(ihel)=four*esq**2*gg*
     &   (-qq*(q1*Q(qi)+l(qi)*r1*propQL4)*(-Q(qi))*a70h4(1,h2,h3,h4,2))
         endif
         ihel=ihel+1
c----
      enddo
      enddo
      enddo
c-----overcount ihel
      ihel=ihel-1
c-----choose part
      do i=1,32
         if (ipsgen==1) then
         msqAA(i)=abs(m70hA(i))**2
         msqAB(i)=2._dp*real(conjg(m70hA(i))*m70hB(i))
         endif
         if (ipsgen==2) then
         msqBB(i)=abs(m70hB(i))**2
         msqBC(i)=2._dp*real(conjg(m70hB(i))*m70hC(i))
         endif
         if (ipsgen==3) then
         msqCC(i)=abs(m70hC(i))**2
         msqAC(i)=2._dp*real(conjg(m70hA(i))*m70hC(i))
         msqCD(i)=2._dp*real(conjg(m70hC(i))*m70hD(i))
         endif
         if (ipsgen==4) then
         msqDD(i)=abs(m70hD(i))**2
         msqAD(i)=2._dp*real(conjg(m70hA(i))*m70hD(i))
         msqBD(i)=2._dp*real(conjg(m70hB(i))*m70hD(i))
         endif
      enddo
      do i=1,32
c---------commented pieces below only for crosscheck with old PSgen
c         msqsum(i)= msqAA(i)+msqBB(i)+msqCC(i)+msqDD(i)+msqAB(i)
c     &             +msqAC(i)+msqAD(i)+msqBC(i)+msqBD(i)+msqCD(i)
         if (ipsgen==1) then
            msqsum(i)=msqAA(i)+msqAB(i) !AA+AB
         elseif (ipsgen==2) then
            msqsum(i)=msqBB(i)+msqBC(i) !BB+BC
         elseif (ipsgen==3) then
            msqsum(i)=msqCC(i)+msqAC(i)+msqCD(i) !CC+AC+CD
         elseif (ipsgen==4) then
            msqsum(i)=msqDD(i)+msqAD(i)+msqBD(i) !D.e+_dpA.e+_dpBD
         endif
      enddo
c-----square them up and sum them up
      msq=zip
      do i=1,32
         msq=msq+msqsum(i)
      enddo
c-----include color factors and identical particle factors
      msq=8._dp/2._dp*msq
c-----done
      return
      end

