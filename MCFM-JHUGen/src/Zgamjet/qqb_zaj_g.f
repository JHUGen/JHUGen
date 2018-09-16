      subroutine qqb_zaj_g(p,msq)
      implicit none
      include 'types.f'
*******************************************************************
*  Author: H. Hartanto                                            *
*                                                                 *
*  return matrix element squared for                              *
*  0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7)  *
*******************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      real(dp):: msq3(-nf:nf,-nf:nf)
      real(dp):: qqb_gg(2),qbq_gg(2),gg_qqb(2)
      real(dp):: qg_qg(2),qbg_qbg(2),gq_gq(2),gqb_gqb(2)
      real(dp):: iib_jjb(4),ibi_jjb(4)
      real(dp):: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
      real(dp):: qiqj(2,2),qjqi(2,2),qbiqbj(2,2),qbjqbi(2,2)
      real(dp):: qiqbj(2,2),qjqbi(2,2),qbiqj(2,2),qbjqi(2,2)
      integer:: i,j,k
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
c-----call squared matrix elements from each subprocesses
      call msq_zqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     .qg_qg,qbg_qbg,gq_gq,gqb_gqb)
      call msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
     .iib_iib,ibi_iib,ii_ii,ibib_ibib,
     .qiqj,qbiqbj,qiqbj,qbiqj,qjqi,qbjqbi,qjqbi,qbjqi)
c-----initialize msq
      do i=-nf,nf
      do j=-nf,nf
         msq(i,j)=zip
         msq1(i,j)=zip
         msq2(i,j)=zip
         msq3(i,j)=zip
      enddo
      enddo
c-----fill msq1 and msq3
      do j=-nf,nf
      do k=-nf,nf
c---------fill msq1
          if ((j == 0) .and. (k == 0)) then
            msq1(j,k)=(two*gg_qqb(2)+three*gg_qqb(1))
          elseif ((j == 0) .and. (k < 0)) then
            msq1(j,k)=gqb_gqb(-kk(k))
          elseif ((j == 0) .and. (k > 0)) then
            msq1(j,k)=gq_gq(kk(k))
          elseif ((j > 0) .and. (k == -j)) then
            msq1(j,k)=qqb_gg(jj(j))
          elseif ((j < 0) .and. (k == -j)) then
            msq1(j,k)=qbq_gg(kk(k))
          elseif ((j > 0) .and. (k == 0)) then
            msq1(j,k)=qg_qg(jj(j))
          elseif ((j < 0) .and. (k == 0)) then
            msq1(j,k)=qbg_qbg(-jj(j))
          else
            msq1(j,k)=0._dp
          endif
c---------fill msq3
          if ((j>0).and.(k>0)) then
            if (j<k) then
               msq3(j,k)=qiqj(jj(j),kk(k))
            elseif (k<j) then
               msq3(j,k)=qjqi(kk(k),jj(j))
            endif
          elseif ((j<0).and.(k<0)) then
            if (j<k) then
               msq3(j,k)=qbjqbi(-kk(k),-jj(j))
            elseif (k<j) then
               msq3(j,k)=qbiqbj(-jj(j),-kk(k))
            endif
          elseif ((j>0).and.(k<0)) then
            if (j<abs(k)) then
               msq3(j,k)=qiqbj(jj(j),-kk(k))
            elseif (abs(k)<j) then
               msq3(j,k)=qjqbi(-kk(k),jj(j))
            endif
          elseif ((j<0).and.(k>0)) then
            if (abs(j)<k) then
               msq3(j,k)=qbiqj(-jj(j),kk(k))
            elseif (k<abs(j)) then
               msq3(j,k)=qbjqi(kk(k),-jj(j))
            endif
          endif
c---------
      enddo
      enddo
c-----fill msq2
      msq2(2,-2)=iib_iib(2)+iib_jjb(3)+3._dp*iib_jjb(1)
      msq2(4,-4)=msq2(2,-2)
      msq2(1,-1)=iib_iib(1)+2._dp*iib_jjb(2)+2._dp*iib_jjb(4)
      msq2(3,-3)=msq2(1,-1)
      msq2(5,-5)=msq2(1,-1)
c-----
      msq2(-2,2)=ibi_iib(2)+ibi_jjb(3)+3._dp*ibi_jjb(1)
      msq2(-4,4)=msq2(-2,2)
      msq2(-1,1)=ibi_iib(1)+2._dp*ibi_jjb(2)+2._dp*ibi_jjb(4)
      msq2(-3,3)=msq2(-1,1)
      msq2(-5,5)=msq2(-1,1)
c-----
      msq2(2,2)=ii_ii(2)
      msq2(4,4)=msq2(2,2)
      msq2(1,1)=ii_ii(1)
      msq2(3,3)=msq2(1,1)
      msq2(5,5)=msq2(1,1)
c-----
      msq2(-2,-2)=ibib_ibib(2)
      msq2(-4,-4)=msq2(-2,-2)
      msq2(-1,-1)=ibib_ibib(1)
      msq2(-3,-3)=msq2(-1,-1)
      msq2(-5,-5)=msq2(-1,-1)
c-----sum msq1,2,3
      do i=-nf,nf
      do j=-nf,nf
         msq(i,j)=msq1(i,j)+msq2(i,j)+msq3(i,j)
      enddo
      enddo
c-----done
      return
      end


      subroutine msq_zqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     & qg_qg,qbg_qbg,gq_gq,gqb_gqb)
      implicit none
      include 'types.f'
*****************************************************************
* return averaged matrix element squared for                    *
* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+gam+g+g+lb+l) amplitude                  *
*****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4)
      real(dp):: pnt(mxpart,4)
      integer:: i,j,k
      complex(dp):: a70h1(2,2,2,2,2,2),a70h3(2,2,2,2,2,2)
      real(dp):: qqb_gg(2),qbq_gg(2),gg_qqb(2)
      real(dp):: qg_qg(2),qbg_qbg(2),gq_gq(2),gqb_gqb(2)
      integer:: i1,i2,i3,i4,i5,i6
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
         qqb_gg(j)=zip
         qbq_gg(j)=zip
         gg_qqb(j)=zip
         qg_qg(j)=zip
         qbg_qbg(j)=zip
         gq_gq(j)=zip
         gqb_gqb(j)=zip
      enddo
c-----calculate sipnor products, invariants
      call spinoru(7,pnt,za,zb)
c-----call color ordered amplitudes and compute msq for qqb_gg
      call zagg_a70h(1,2,3,4,5,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,qqb_gg(j))
      enddo
c-----call color ordered amplitudes and compute msq for qbq_gg
      call zagg_a70h(5,2,3,4,1,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,qbq_gg(j))
      enddo
c-----call color ordered amplitudes and compute msq for gg_qqb
      call zagg_a70h(3,2,1,5,4,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,gg_qqb(j))
      enddo
c-----call color ordered amplitudes and compute msq for qg_qg
      call zagg_a70h(3,2,1,4,5,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,qg_qg(j))
      enddo
c-----call color ordered amplitudes and compute msq for qbg_qbg
      call zagg_a70h(5,2,1,4,3,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,qbg_qbg(j))
      enddo
c-----call color ordered amplitudes and compute msq for gq_gq
      call zagg_a70h(3,2,5,4,1,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,gq_gq(j))
      enddo
c-----call color ordered amplitudes and compute msq for gqb_gqb
      call zagg_a70h(1,2,5,4,3,6,7,a70h1,a70h3)
      do j=1,2
         call zagg_m70sq(j,a70h1,a70h3,gqb_gqb(j))
      enddo
c-----averaging and put identical particle factor
      do j=1,2
         qqb_gg(j)=qqb_gg(j)*half*aveqq
         qbq_gg(j)=qbq_gg(j)*half*aveqq
         gg_qqb(j)=gg_qqb(j)*avegg
         qg_qg(j)=qg_qg(j)*aveqg
         qbg_qbg(j)=qbg_qbg(j)*aveqg
         gq_gq(j)=gq_gq(j)*aveqg
         gqb_gqb(j)=gqb_gqb(j)*aveqg
      enddo
c-----done
      return
      end


      subroutine msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
     & iib_iib,ibi_iib,ii_ii,ibib_ibib,
     & qiqj,qbiqbj,qiqbj,qbiqj,qjqi,qbjqbi,qjqbi,qbjqi)
      implicit none
      include 'types.f'
*****************************************************************
* return matrix element squared for                             *
* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+Q+Qb+gam+lb+l) amplitude                 *
*****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'ewcharge.f'
      real(dp):: p(mxpart,4),pnt(mxpart,4)
      real(dp):: iib_jjb(4),ibi_jjb(4)
      real(dp):: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
      real(dp):: qiqj(2,2),qjqi(2,2),qbiqbj(2,2),qbjqbi(2,2)
      real(dp):: qiqbj(2,2),qjqbi(2,2),qbiqj(2,2),qbjqi(2,2)
      real(dp):: xqiqj(4),xqjqi(4),xqbiqbj(4),xqbjqbi(4)
      real(dp):: xqiqbj(4),xqjqbi(4),xqbiqj(4),xqbjqi(4)
      integer:: i,j,k
c-----convert to Nagy-Trocsanyi momentum convention
      do i=1,4
         pnt(1,i)=p(2,i)
         pnt(2,i)=p(1,i)
         pnt(3,i)=p(6,i)
         pnt(4,i)=p(7,i)
         pnt(5,i)=p(5,i)
         pnt(6,i)=p(4,i)
         pnt(7,i)=p(3,i)
      enddo
c-----initialize msq
      do i=1,4
         iib_jjb(i)=zip
         ibi_jjb(i)=zip
      enddo
      do i=1,2
         iib_iib(i)  =zip
         ibi_iib(i)  =zip
         ii_ii(i)    =zip
         ibib_ibib(i)=zip
      enddo
      do i=1,2
      do j=1,2
         qiqj(i,j)=zip
         qbiqbj(i,j)=zip
         qiqbj(i,j)=zip
         qbiqj(i,j)=zip
         qjqi(i,j)=zip
         qbjqbi(i,j)=zip
         qjqbi(i,j)=zip
         qbjqi(i,j)=zip
      enddo
      enddo
c-----evaluate spinor products, invariants
      call spinoru(7,pnt,za,zb)
c-----
      call zqqbQQba_msqij(1,2,3,4,5,6,7,iib_jjb) 
c-----
      call zqqbQQba_msqij(2,1,3,4,5,6,7,ibi_jjb) 
c-----
      call zqqbQQba_msqii(1,2,3,4,5,6,7,iib_iib) 
c-----
      call zqqbQQba_msqii(2,1,3,4,5,6,7,ibi_iib) 
c-----
      call zqqbQQba_msqii(3,1,4,2,5,6,7,ii_ii) 
c-----
      call zqqbQQba_msqii(1,3,2,4,5,6,7,ibib_ibib)
c-----
      call zqqbQQba_msqij(3,2,4,1,5,6,7,xqiqj)
      call zqqbQQba_msqij(4,1,3,2,5,6,7,xqjqi)
      qiqj(1,1)=xqiqj(4)
      qiqj(1,2)=xqjqi(1)
      qiqj(2,1)=xqiqj(1)
      qiqj(2,2)=xqiqj(3)
      qjqi(1,1)=xqjqi(4)
      qjqi(1,2)=xqiqj(1)
      qjqi(2,1)=xqjqi(1)
      qjqi(2,2)=xqjqi(3)
c-----
      call zqqbQQba_msqij(2,3,1,4,5,6,7,xqbiqbj)
      call zqqbQQba_msqij(1,4,2,3,5,6,7,xqbjqbi)
      qbiqbj(1,1)=xqbiqbj(4)
      qbiqbj(1,2)=xqbjqbi(1)
      qbiqbj(2,1)=xqbiqbj(1)
      qbiqbj(2,2)=xqbiqbj(3)
      qbjqbi(1,1)=xqbjqbi(4)
      qbjqbi(1,2)=xqbiqbj(1)
      qbjqbi(2,1)=xqbjqbi(1)
      qbjqbi(2,2)=xqbjqbi(3)
c-----
      call zqqbQQba_msqij(3,2,1,4,5,6,7,xqiqbj) 
      call zqqbQQba_msqij(1,4,3,2,5,6,7,xqjqbi) 
      qiqbj(1,1)=xqiqbj(4)
      qiqbj(1,2)=xqjqbi(1)
      qiqbj(2,1)=xqiqbj(1)
      qiqbj(2,2)=xqiqbj(3)
      qjqbi(1,1)=xqjqbi(4)
      qjqbi(1,2)=xqiqbj(1)
      qjqbi(2,1)=xqjqbi(1)
      qjqbi(2,2)=xqjqbi(3)
c-----
      call zqqbQQba_msqij(2,4,3,1,5,6,7,xqbiqj) 
      call zqqbQQba_msqij(3,1,2,4,5,6,7,xqbjqi)
      qbiqj(1,1)=xqbiqj(4)
      qbiqj(1,2)=xqbjqi(1)
      qbiqj(2,1)=xqbiqj(1)
      qbiqj(2,2)=xqbiqj(3)
      qbjqi(1,1)=xqbjqi(4)
      qbjqi(1,2)=xqbiqj(1)
      qbjqi(2,1)=xqbjqi(1)
      qbjqi(2,2)=xqbjqi(3)
c-----averaging and put identical particle factors
      do i=1,4
         iib_jjb(i)=iib_jjb(i)*aveqq
         ibi_jjb(i)=ibi_jjb(i)*aveqq
      enddo
      do i=1,2
         iib_iib(i)  =iib_iib(i)*aveqq
         ibi_iib(i)  =ibi_iib(i)*aveqq
         ii_ii(i)    =ii_ii(i)*aveqq*half
         ibib_ibib(i)=ibib_ibib(i)*aveqq*half
      enddo
      do i=1,2
      do j=1,2
         qiqj(i,j)=qiqj(i,j)*aveqq
         qjqi(i,j)=qjqi(i,j)*aveqq
         qbiqbj(i,j)=qbiqbj(i,j)*aveqq
         qbjqbi(i,j)=qbjqbi(i,j)*aveqq
         qiqbj(i,j)=qiqbj(i,j)*aveqq
         qbiqj(i,j)=qbiqj(i,j)*aveqq
         qjqbi(i,j)=qjqbi(i,j)*aveqq
         qbjqi(i,j)=qbjqi(i,j)*aveqq
      enddo
      enddo
c-----done
      return
      end

