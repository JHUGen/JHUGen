      subroutine qqb_hflgam_g(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     January, 2013.                                                   *
************************************************************************
*                                                                      *
*     Matrix element for gamma + 2 parton production                   *
*     in order alpha_s, averaged over initial colours and spins        *
*                                                                      *
*     q(-p1)+qbar(-p2) --> gamma(p3)+f(p4)+f(p5)                       *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'sprods_com.f'
      include 'heavyflav.f'
      include 'masses.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     .  Bigagam,Bigbgam,Bigcgam,
     . qr_qr(2,2),rq_qr(2,2),ra_ra(2,2),ar_ra(2,2),
     . qa_rb(2,2),aq_rb(2,2),
     . qq_qq(2),qa_qa(2),aq_qa(2),
     . qg_qg(2),gq_qg(2),gg_qa(2)
      integer j,k
      
      call dotem(5,p,s)

      do j=1,2
      do k=1,2
      qr_qr(j,k)=+aveqq*Bigagam(1,2,4,5,3,j,k)
      rq_qr(j,k)=+aveqq*Bigagam(2,1,4,5,3,j,k)
c      ab_ab(j,k)=+aveqq*Bigagam(4,5,1,2,3,j,k)

      ra_ra(j,k)=+aveqq*Bigagam(5,1,2,4,3,j,k)
      ar_ra(j,k)=+aveqq*Bigagam(5,2,1,4,3,j,k)

      qa_rb(j,k)=+aveqq*Bigagam(1,5,2,4,3,j,k)
      aq_rb(j,k)=+aveqq*Bigagam(2,5,1,4,3,j,k)

      enddo
      qq_qq(j)=+aveqq*Bigbgam(1,2,4,5,3,j)*0.5d0
c      aa_aa(j)=+aveqq*Bigbgam(4,5,1,2,3,j)*0.5d0
      qa_qa(j)=+aveqq*Bigbgam(1,5,4,2,3,j)
      aq_qa(j)=+aveqq*Bigbgam(2,5,4,1,3,j)

c      ag_ag(j)=-aveqg*Bigcgam(1,4,5,2,3,j)
      qg_qg(j)=-aveqg*Bigcgam(4,1,5,2,3,j)
c      ga_ag(j)=-aveqg*Bigcgam(2,4,5,1,3,j)
      gq_qg(j)=-aveqg*Bigcgam(4,2,5,1,3,j)

      gg_qa(j)=+avegg*Bigcgam(4,5,2,1,3,j)
      enddo

c      do j=-flav,flav
c      do k=-flav,flav
c      if ((j.eq.flav).or.(k.eq.flav)
c     & .or.(j.eq.-k).or.((j.eq.0).and.(k.eq.0)))
c     & msq(j,k)=10d0
c      enddo
c      enddo
      
      msq(:,:)=zip
      
      if (flav .eq. 4) then 
        msq(+1,+4)=rq_qr(2,1)
        msq(+2,+4)=rq_qr(2,2)
        msq(+3,+4)=rq_qr(2,1)
        msq(-4,+4)=aq_qa(2)
        msq(-3,+3)=aq_rb(1,2)
        msq(-3,+4)=ar_ra(1,2)
        msq(-2,+2)=aq_rb(2,2)
        msq(-2,+4)=ar_ra(2,2)
        msq(-1,+4)=ar_ra(1,2)
        msq(-1,+1)=aq_rb(1,2)
        msq(+1,-1)=qa_rb(1,2)
        msq(+2,-2)=qa_rb(2,2)
        msq(+3,-3)=qa_rb(1,2)
        msq(+4,-4)=qa_qa(2)
        msq(+4,-3)=ra_ra(1,2)
        msq(+4,-2)=ra_ra(2,2)
        msq(+4,-1)=ra_ra(1,2)
        msq(+4,+0)=qg_qg(2)
        msq(+4,+1)=qr_qr(2,1)
        msq(+4,+2)=qr_qr(2,2)
        msq(+4,+3)=qr_qr(2,1)
        msq(+4,+4)=qq_qq(2)
        msq(+0,+0)=gg_qa(2)
        msq(+0,+4)=gq_qg(2)

      elseif (flav .eq. 5) then 
        msq(+1,+5)=rq_qr(1,1)
        msq(+2,+5)=rq_qr(1,2)
        msq(+4,+5)=rq_qr(1,2)
        msq(+3,+5)=rq_qr(1,1)
        msq(-5,+5)=aq_qa(1)
        msq(-3,+3)=aq_rb(1,1)
        msq(-3,+5)=ar_ra(1,1)
        
        msq(-2,+2)=aq_rb(2,1)
        msq(-4,+4)=aq_rb(2,1)
        msq(+2,-2)=qa_rb(2,1)
        msq(+4,-4)=qa_rb(2,1)
        
        msq(-2,+5)=ar_ra(2,1)
        msq(-4,+5)=ar_ra(2,1)
        msq(-1,+5)=ar_ra(1,1)
        msq(-1,+1)=aq_rb(1,1)
        msq(+1,-1)=qa_rb(1,1)
                
        msq(+3,-3)=qa_rb(1,1)
        msq(+5,-5)=qa_qa(1)
        msq(+5,-4)=ra_ra(2,1)
        msq(+5,-3)=ra_ra(1,1)
        msq(+5,-2)=ra_ra(2,1)
        msq(+5,-1)=ra_ra(1,1)
        msq(+5,+0)=qg_qg(1)
        msq(+5,+1)=qr_qr(1,1)
        msq(+5,+2)=qr_qr(1,2)
        msq(+5,+3)=qr_qr(1,1)
        msq(+5,+4)=qr_qr(1,2)
        msq(+5,+5)=qq_qq(1)
        msq(+0,+0)=gg_qa(1)
        msq(+0,+5)=gq_qg(1)

      else      
        write(6,*) 'qqb_hflgam:unimplemented value of flav'
        stop
        
      endif
      
c--- remove contributions for qq~ -> photon+QQ~ with m(QQ~)<4*mQ^2  
c--- (to screen collinear divergence in g->QQ~)    
      if     (flav .eq. 4) then
        if (s(4,5) .lt. 4d0*mcsq) then
          do j=1,4
          msq(j,-j)=0d0
          msq(-j,j)=0d0
          enddo
        endif
      elseif (flav. eq. 5) then
        if (s(4,5) .lt. 4d0*mbsq) then
          do j=1,5
          msq(j,-j)=0d0
          msq(-j,j)=0d0
          enddo
        endif
      endif
      
      return
      end
