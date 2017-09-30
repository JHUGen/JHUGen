      subroutine qqb_gamgam_g(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John M. Campbell                         *
*     December, 2010.                                                  *
************************************************************************
*                                                                      *
*     Matrix element for gamma + gamma + 1 parton production           *
*     in order alpha_s, averaged over initial colours and spins        *
*                                                                      *
*     q(-p1)+qbar(-p2) --> gamma(p3)+gamma(p4)+gluon(p5)               *
*                                                                      *
*     (and all crossings)                                              *
*                                                                      *
************************************************************************
      include 'constants.f'
      include 'zprods_com.f'
      include 'noglue.f'
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     & Cgamgam,ggtogagag,qa_gagag(2),
     & ag_gagaa(2),qg_gagaq(2),ga_gagaa(2),gq_gagaq(2),gg_gaga
      integer j,k
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      
      call spinoru(5,p,za,zb)
c      call dotem(5,p,s)

      do j=1,2
      ag_gagaa(j)=-aveqg*Cgamgam(5,1,3,4,2,j)
      qg_gagaq(j)=-aveqg*Cgamgam(1,5,3,4,2,j)
      ga_gagaa(j)=-aveqg*Cgamgam(5,2,3,4,1,j)
      gq_gagaq(j)=-aveqg*Cgamgam(2,5,3,4,1,j)

      qa_gagag(j)=+aveqq*Cgamgam(1,2,3,4,5,j)
      enddo

c--- averaging factor already included in ggtogagag 

      if (omitgg) then
        gg_gaga=0d0
      else
        gg_gaga=ggtogagag()
      endif
      
      do j=-nf,nf
      do k=-nf,nf

c--set msq=0 to initalize
      msq(j,k)=0d0
C--qa      
      if ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
          msq(j,k)=qa_gagag(jj(j))
          endif
C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
          msq(j,k)=qa_gagag(kk(k))
          endif
C--qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_gagaq(jj(j))
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_gagaa(-jj(j))
C--gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gagaq(kk(k))
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_gagaa(-kk(k))

C--gg     
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_gaga
      endif

      enddo
      enddo

      return
      end
