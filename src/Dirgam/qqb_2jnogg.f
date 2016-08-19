      subroutine qqb_2jnogg(p,msq)
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qqij_ij,aaij_ij,
     .  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     .  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     .  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     .  smalla,smallb,smallc 
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(ss,uu,tt)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

     

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)
      
    
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0d0
      

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=qaii_ii+dfloat(nf-1)*qaii_jj+qa_gg
          else
            msq(j,k)=qaij_ij
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aqii_ii+dfloat(nf-1)*aqii_jj+aq_gg
          else
            msq(j,k)=aqij_ij
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end

      subroutine qqb_2jnoggswap(pin,msq)
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision pin(mxpart,4)
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     .  qqij_ij,aaij_ij,
     .  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     .  qqii_ii,aaii_ii,
     .  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     .  smalla,smallb,smallc 
      
      p=0d0

      do j=1,4 
         p(1,j)=pin(1,j)
         p(2,j)=pin(2,j)
         p(4,j)=pin(3,j) 
         p(3,j)=pin(4,j) 
      enddo
      


      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(ss,uu,tt)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
     

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)
      
    
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0d0
      

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=qaij_ij+0d0*(dfloat(nf-1)*qaii_jj+qa_gg)
          else
            msq(j,k)=qaij_ij
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aqij_ij+0d0*(dfloat(nf-1)*aqii_jj+aq_gg)
          else
            msq(j,k)=aqij_ij
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo

      
      return
      end


      subroutine qqb_2j_t(p,msq)
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
!----- This routine contains t-channel type diagrams which will contain an intial-final 
!---- photon singularity other terms are set to zero  
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,
     &  qqii_ii,aaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,ss,tt,uu,
     &  smalla,smallb,smallc
C    & ,gq_qg,qaii_jj,aqii_ii,qaii_ii,
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

     
      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
c      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
!      aqij_ij=fac*aveqq*smalla(uu,tt,ss)
      aqij_ij=fac*aveqq*smalla(ss,uu,tt)
      
      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
c      aqii_ii=fac*aveqq*smallb(tt,uu,ss)
c      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

     

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
c      gq_qg=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)
      
c      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)
      
    
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0d0
      

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
             
            msq(j,k)=qaij_ij !--- want only t-channel scattering
          else
            msq(j,k)=qaij_ij
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=aqij_ij+0d0*(dfloat(nf-1)*aqii_jj+aq_gg) !-- want only t-channel scattering
          else
            msq(j,k)=aqij_ij
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
 !        msq(j,k)=gq_qg
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end



      subroutine qqb_2j_s(p,msq)
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc
      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)

      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(uu,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

     

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
c      gq_qg=-fac*aveqg*smallc(tt,ss,uu)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)
      
    
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0d0
      

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=0d0*qaii_ii+qaii_jj+0d0*qa_gg  !-- want 1 s-channel scattering
          else
            msq(j,k)=qaij_ij*0d0
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=0d0*aqii_ii+aqii_jj+0d0*aq_gg !-- want 1 s-channel scattering 
          else
            msq(j,k)=aqij_ij*0d0
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
 !           pause
 !        msq(j,k)=gq_qg
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end


      subroutine qqb_2j_sswap(pin,msq)
c---- Calcuated g g -> q q without p p -> g g needed for photon dipoles
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      double precision pin(mxpart,4) 
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac,
     &  qqij_ij,aaij_ij,
     &  qaij_ij,aqij_ij,aqii_jj,qaii_jj,
     &  qqii_ii,aaii_ii,aqii_ii,qaii_ii,
     &  aq_gg,gq_gq,ga_ga,qg_qg,ag_ag,gg_qa,qa_gg,ss,tt,uu,
     &  smalla,smallb,smallc 

      p=0d0
      
      do j=1,4 
         p(1,j)=pin(1,j)
         p(2,j)=pin(2,j)
         p(4,j)=pin(3,j) 
         p(3,j)=pin(4,j) 
      enddo
      



      call dotem(4,p,s)
      fac=gsq**2
      ss=s(1,2)
      tt=s(1,3)
      uu=s(2,3)



      qqij_ij=fac*aveqq*smalla(ss,tt,uu)
      aaij_ij=fac*aveqq*smalla(ss,tt,uu)
      qaii_jj=fac*aveqq*smalla(tt,ss,uu)
      qaij_ij=fac*aveqq*smalla(ss,tt,uu)
      aqii_jj=fac*aveqq*smalla(tt,ss,uu)
      aqij_ij=fac*aveqq*smalla(ss,tt,ss)

      qqii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aaii_ii=fac*aveqq*smallb(ss,tt,uu)*half
      aqii_ii=fac*aveqq*smallb(uu,tt,ss)
      qaii_ii=fac*aveqq*smallb(uu,tt,ss)

     

      qg_qg=-fac*aveqg*smallc(tt,ss,uu)
      ag_ag=-fac*aveqg*smallc(tt,uu,ss)
      gq_gq=-fac*aveqg*smallc(uu,ss,tt)
      ga_ga=-fac*aveqg*smallc(uu,tt,ss)

      qa_gg=+fac*aveqq*smallc(ss,tt,uu)*half
      gg_qa=+fac*avegg*smallc(ss,tt,uu)
      
    
      aq_gg=+fac*aveqq*smallc(ss,tt,uu)*half*0d0
      

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0d0
C--qq      
      if ((j .gt. 0) .and. (k .gt. 0)) then
          if (j .eq. k) then
            msq(j,k)=qqii_ii
          else
            msq(j,k)=qqij_ij
          endif

C--qa      
      elseif ((j .gt. 0) .and. (k .lt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=0d0*qaii_ii+qaii_jj+0d0*qa_gg  !-- want 1 s-channel scattering
          else
            msq(j,k)=qaij_ij*0d0
          endif

C--aa      
      elseif ((j .lt. 0) .and. (k .lt. 0)) then
          if (j .eq. k) then
            msq(j,k)=aaii_ii
          else
            msq(j,k)=aaij_ij
          endif

C--aq      
      elseif ((j .lt. 0) .and. (k .gt. 0)) then
          if (j .eq. -k) then
            msq(j,k)=0d0*aqii_ii+aqii_jj+0d0*aq_gg !-- want 1 s-channel scattering 
          else
            msq(j,k)=aqij_ij*0d0
          endif

C--qg_qg      
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
            msq(j,k)=qg_qg
C--ag      
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
            msq(j,k)=ag_ag
C--gq_gq      
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
            msq(j,k)=gq_gq
C--ga      
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
            msq(j,k)=ga_ga
C--gg      
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
            msq(j,k)=gg_qa
      endif

      enddo
      enddo


      return
      end


