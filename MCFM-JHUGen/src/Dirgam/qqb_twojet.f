      subroutine qqb_twojet(p,msq)
********************************************************************************
*                                                                              *
*   Author: J. Campbell, September 2014                                        *
*   Matrix elements for two-jet production                                     *
*                                                                              *
*   Note that this differs from routine qqb_2jet.f (used for direct photon     *
*   fragmentation contribution) since gg->gg matrix element is included here;  *
*   in addition, ordering of final state partons differs for some channels     *
*                                                                              *
********************************************************************************
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'nflav.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  QjQk_QjQk,AjAk_AjAk,
     &  QjAk_QjAk,AjQk_AjQk,AjQj_AkQk,QjAj_QkAk,
     &  QjQj_QjQj,AjAj_AjAj,AjQj_AjQj,QjAj_QjAj,
     &  AQ_GG,GQ_GQ,GA_GA,QG_QG,AG_AG,GG_GG,GG_QA,QA_GG,ss,tt,uu,
     &  smalla,smallb,smallc,smalld 

      call dotem(4,p,s)
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      fac=gsq**2
      
c--- four quarks
      QjQk_QjQk=fac*aveqq*smalla(ss,tt,uu)
      AjAk_AjAk=fac*aveqq*smalla(ss,tt,uu)
      
      QjAj_QkAk=fac*aveqq*smalla(tt,ss,uu)
      QjAk_QjAk=fac*aveqq*smalla(uu,tt,ss)
      
      AjQj_AkQk=fac*aveqq*smalla(tt,ss,uu)
      AjQk_AjQk=fac*aveqq*smalla(uu,tt,ss)

      QjQj_QjQj=fac*aveqq*smallb(ss,tt,uu)*half
      AjAj_AjAj=fac*aveqq*smallb(ss,tt,uu)*half
      QjAj_QjAj=fac*aveqq*smallb(uu,tt,ss)
      AjQj_AjQj=fac*aveqq*smallb(uu,tt,ss)

c--- two quarks, two gluons
      GG_QA=fac*avegg*smallc(ss,tt,uu)
      QA_GG=fac*aveqq*smallc(ss,tt,uu)*half
      AQ_GG=fac*aveqq*smallc(ss,uu,tt)*half

      QG_QG=-fac*aveqg*smallc(tt,ss,uu)
      AG_AG=-fac*aveqg*smallc(tt,uu,ss)
      GQ_GQ=-fac*aveqg*smallc(tt,uu,ss)
      GA_GA=-fac*aveqg*smallc(tt,ss,uu)

c--- four gluons
      GG_GG=fac*avegg*smalld(ss,tt,uu)*half

c--- initalize matrix elements to zero
      msq(:,:)=0._dp

c      nflav=4

      do j=-nflav,nflav
      do k=-nflav,nflav

C-- QQ
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=QjQj_QjQj
          else
            msq(j,k)=QjQk_QjQk
          endif

C-- QA
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=QjAj_QjAj+dfloat(nflav-1)*QjAj_QkAk+QA_GG
          else
            msq(j,k)=QjAk_QjAk
          endif

C-- AA
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=AjAj_AjAj
          else
            msq(j,k)=AjAk_AjAk
          endif

C-- AQ
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=AjQj_AjQj+dfloat(nflav-1)*AjQj_AkQk+AQ_GG
          else
            msq(j,k)=AjQk_AjQk
          endif

C-- QG
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=QG_QG

C-- AG
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=AG_AG

C-- GQ
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=GQ_GQ

C-- GA
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=GA_GA

C-- GG
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=GG_GG+dfloat(nflav)*GG_QA

      endif

      enddo
      enddo

      return
      end



