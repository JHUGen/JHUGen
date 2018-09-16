      subroutine qqb_dirgam_g(p,msq)
      implicit none
      include 'types.f'
      
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
*                                                                      *
*     Matrix element for gamma + 2 parton production                   *
*     in order alpha_s, averaged over initial colours and spins        *
*                                                                      *
*     q(-p1)+qbar(-p2) --> gamma(p3)+f(p4)+f(p5)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'msqbits.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     &  Bigagam,Bigbgam,Bigcgam,
     & qr_qr(2,2),ab_ab(2,2),ra_ra(2,2),ar_ra(2,2),
     & qa_rb(2,2),aq_rb(2,2),
     & qq_qq(2),aa_aa(2),
     & qa_qa(2),aq_qa(2),
     & aq_gg(2),qa_gg(2),
     & ag_ag(2),qg_qg(2),ga_ag(2),gq_qg(2),
     & gg_qa(2)
      integer:: j,k
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      call dotem(5,p,s)

      do j=1,2
      do k=1,2
      qr_qr(j,k)=+aveqq*Bigagam(1,2,4,5,3,j,k)
      ab_ab(j,k)=+aveqq*Bigagam(4,5,1,2,3,j,k)
      ra_ra(j,k)=+aveqq*Bigagam(1,5,4,2,3,j,k)
      ar_ra(j,k)=+aveqq*Bigagam(5,2,1,4,3,j,k)
      qa_rb(j,k)=+aveqq*Bigagam(1,5,2,4,3,j,k)
      aq_rb(j,k)=+aveqq*Bigagam(2,5,1,4,3,j,k)
      enddo
      qq_qq(j)=+aveqq*Bigbgam(1,2,4,5,3,j)*0.5_dp
      aa_aa(j)=+aveqq*Bigbgam(4,5,1,2,3,j)*0.5_dp
      qa_qa(j)=+aveqq*Bigbgam(1,5,4,2,3,j)
      aq_qa(j)=+aveqq*Bigbgam(2,5,4,1,3,j)

      aq_gg(j)=+aveqq*Bigcgam(1,2,4,5,3,j)*0.5_dp
      qa_gg(j)=+aveqq*Bigcgam(2,1,4,5,3,j)*0.5_dp

      ag_ag(j)=-aveqg*Bigcgam(1,4,5,2,3,j)
      qg_qg(j)=-aveqg*Bigcgam(4,1,5,2,3,j)
      ga_ag(j)=-aveqg*Bigcgam(2,4,5,1,3,j)
      gq_qg(j)=-aveqg*Bigcgam(4,2,5,1,3,j)

      gg_qa(j)=+avegg*Bigcgam(4,5,2,1,3,j)
      enddo

c--- fill ancillary array that is used for fragmentation
c--- contributions in gmgmjt process
      msqbits(ddb_ddb)=qa_qa(1)
      msqbits(ddb_ssb)=qa_rb(1,1)
      msqbits(ddb_uub)=qa_rb(1,2)
      msqbits(uub_uub)=qa_qa(2)
      msqbits(uub_ccb)=qa_rb(2,2)
      msqbits(uub_ddb)=qa_rb(2,1)
      msqbits(dbd_ddb)=aq_qa(1)
      msqbits(dbd_ssb)=aq_rb(1,1)
      msqbits(dbd_uub)=aq_rb(1,2)
      msqbits(ubu_uub)=aq_qa(2)
      msqbits(ubu_ccb)=aq_rb(2,2)
      msqbits(ubu_ddb)=aq_rb(2,1)

      do j=-nf,nf
      do k=-nf,nf

c--set msq=0 to initalizec--set msq=0 to initalize
      msq(j,k)=0._dp
C--qq      
      if ((j > 0) .and. (k > 0)) then
          if (j == k) then
c---- u u -> gamma u u
            msq(j,k)=qq_qq(jj(j))
          else
c---- u d -> gamma u d
            msq(j,k)=qr_qr(jj(j),kk(k))
          endif

C--aa      
      elseif ((j < 0) .and. (k < 0)) then
          if (j == k) then
c---- ub ub -> gamma ub ub
            msq(j,k)=aa_aa(-jj(j))
          else
c---- ub db -> gamma ub db
            msq(j,k)=ab_ab(-jj(j),-kk(k))
          endif

C--qa      
      elseif ((j > 0) .and. (k < 0)) then
          if (j == -k) then             
            msq(j,k)=qa_qa(jj(j))+qa_gg(jj(j))
            if (jj(j)==1) then
c---- (d db -> gamma d db) + 2*(d db -> gamma s sb) + 2*(d db -> gamma u ub) + (d db -> gamma g g) 
               msq(j,k)=msq(j,k)
     &              +2._dp*qa_rb(jj(j),1)+2._dp*qa_rb(jj(j),2)
            elseif (jj(j)==2) then
c---- (u ub -> gamma u ub) + (u ub -> gamma c cb) + 3*(u ub -> gamma d db) + (u ub -> gamma g g) 
                msq(j,k)=msq(j,k)+qa_rb(jj(j),2)+3._dp*qa_rb(jj(j),1)
            endif
          else
c---- u db -> gamma u db
            msq(j,k)=ra_ra(jj(j),-kk(k))
          endif

C--aq      
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
            msq(j,k)=aq_qa(kk(k))+aq_gg(kk(k))
            if (kk(k)==1) then
c---- (db d -> gamma d db) + 2*(db d -> gamma s sb) + 2*(db d -> gamma u ub) + (db d -> gamma g g) 
               msq(j,k)=msq(j,k)+2._dp*aq_rb(-jj(j),1)+2._dp*aq_rb(-jj(j),2)
            elseif (kk(k)==2) then
c---- (ub u -> gamma u ub) + (ub u -> gamma c cb) + 3*(ub u -> gamma d db) + (ub u -> gamma g g) 
               msq(j,k)=msq(j,k)+3._dp*aq_rb(-jj(j),1)+aq_rb(-jj(j),2)
            endif
          else
c---- db u -> gamma u db
            msq(j,k)=ar_ra(-jj(j),kk(k))
          endif

C--qg      
      elseif ((j > 0) .and. (k == 0)) then
c---- u g -> gamma u g
            msq(j,k)=qg_qg(jj(j))

C--ag      
      elseif ((j < 0) .and. (k == 0)) then
c---- ub g -> gamma ub g
            msq(j,k)=ag_ag(-jj(j))

C--gq      
      elseif ((j == 0) .and. (k > 0)) then
c---- g u -> gamma u g
            msq(j,k)=gq_qg(kk(k))

C--ga      
      elseif ((j == 0) .and. (k < 0)) then
c---- g ub -> gamma ub g
            msq(j,k)=ga_ag(-kk(k))

C--gg      
      elseif ((j == 0) .and. (k == 0)) then
c---- 2*(g g -> gamma u ub) + (nf-2)*(g g -> gamma d db)
            msq(j,k)=+2._dp*gg_qa(2)+real(nf-2,dp)*gg_qa(1)
      endif

      enddo
      enddo



      return
      end
