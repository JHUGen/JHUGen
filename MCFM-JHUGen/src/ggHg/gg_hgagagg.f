      subroutine gg_hgagagg(p,msq)
      implicit none
      include 'types.f'
      
c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2)-->H -->  gamma(p3)+gamma(p4)) +g(p_i5=5)+g(p_i6=6) 
c

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'msq_struc.f'
      integer:: j,k,i5,i6
      real(dp):: p(mxpart,4),Asq,fac
      real(dp):: Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625
c     &                     ,Hgggg_1652,Hgggg_1562,Hgggg_1526
      real(dp):: Hqagg,Haqgg,Hgqqg,Hgaag,Hqgqg,Hagag,Hggqa
      real(dp):: Hggqa_ab,Hggqa_ba,Hggqa_sym
      real(dp):: Hqgqg_ab,Hqgqg_ba,Hqgqg_sym
      real(dp):: Hgqqg_ab,Hgqqg_ba,Hgqqg_sym
      real(dp):: Hagag_ab,Hagag_ba,Hagag_sym
      real(dp):: Hgaag_ab,Hgaag_ba,Hgaag_sym
      real(dp):: Hqagg_ab,Hqagg_ba,Hqagg_sym
      real(dp):: Haqgg_ab,Haqgg_ba,Haqgg_sym
      real(dp):: Hqqqq_a,Hqqqq_b,Hqqqq_i
      real(dp):: Hqaqa_a,Hqaqa_b,Hqaqa_i
      real(dp):: Haqaq_a,Haqaq_b,Haqaq_i
      real(dp):: Hqaaq_a,Hqaaq_b,Hqaaq_i
      real(dp):: 
     & Hqrqr,Hqqqq,
     & Habab,Haaaa,
     & Hqarb,Hqaqa,Hqbqb,
     & Haqbr,Haqaq,Hbqbq,
     & Hqaaq,msqgamgam
      real(dp):: msq(-nf:nf,-nf:nf),hdecay
      parameter(i5=5,i6=6)


C---fill spinor products up to maximum number
      call spinoru(i6,p,za,zb)  

C   Deal with Higgs decay to b-bbar
      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

C--four gluon terms
      call h4g(1,2,i5,i6,Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
      
C--two quark two gluon terms
      call hqqggdfm(1,2,i5,i6,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)
      call hqqggdfm(2,1,i5,i6,Haqgg,Haqgg_ab,Haqgg_ba,Haqgg_sym)
C====symmetric in first two arguments, but not the ab, ba terms
c      Haqgg=Hqagg

      call hqqggdfm(1,i5,2,i6,Hqgqg,Hqgqg_ab,Hqgqg_ba,Hqgqg_sym)
      call hqqggdfm(i5,1,2,i6,Hagag,Hagag_ab,Hagag_ba,Hagag_sym)
C====symmetric in first two arguments
c      Hagag=Hqgqg

      call hqqggdfm(2,i5,1,i6,Hgqqg,Hgqqg_ab,Hgqqg_ba,Hgqqg_sym)
      call hqqggdfm(i5,2,1,i6,Hgaag,Hgaag_ab,Hgaag_ba,Hgaag_sym)
C====symmetric in first two arguments
c      Hgaag=Hgqqg

      call hqqggdfm(i6,i5,1,2,Hggqa,Hggqa_ab,Hggqa_ba,Hggqa_sym)
      
C---four quark terms
      call H4qn(1,2,i5,i6,Hqrqr)
      call H4qi(1,2,i5,i6,Hqqqq,Hqqqq_a,Hqqqq_b,Hqqqq_i)
C---four anti-quark terms
c      call H4qn(i5,i6,1,2,Habab)
c      call H4qi(i5,i6,1,2,Haaaa)
      Habab=Hqrqr
      Haaaa=Hqqqq

C-qqb
      call H4qn(1,i6,2,i5,Hqarb)
      call H4qi(1,i6,i5,2,Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i)
      call H4qn(1,i6,i5,2,Hqbqb)
c      write(6,*) 'Hqaqa',Hqaqa_a,Hqaqa_b,Hqaqa_i
c      write(6,*) 'Hqbqb',Hqbqb
C-qbq
      Haqbr=Hqarb
      
      Haqaq=Hqaqa
      Haqaq_a=Hqaqa_a
      Haqaq_b=Hqaqa_b
      Haqaq_i=Hqaqa_i
      Hbqbq=Hqbqb

      do j=fn,nf
      do k=fn,nf
      msq(j,k)=0._dp
      msq_struc(iqr,j,k)=0._dp

      if ((j>0).and.(k>0)) then 
        if (j==k) then
          msq(j,k)=0.5_dp*aveqq*fac*Hqqqq
          msq_struc(iqq_a,j,k)=0.5_dp*aveqq*fac*Hqqqq_a
          msq_struc(iqq_b,j,k)=0.5_dp*aveqq*fac*Hqqqq_b
          msq_struc(iqq_i,j,k)=0.5_dp*aveqq*fac*Hqqqq_i
        else
          msq(j,k)=aveqq*fac*Hqrqr
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif
      
      if ((j<0).and.(k<0)) then 
        if (j==k) then
          msq(j,k)=0.5_dp*aveqq*fac*Haaaa
        else
          msq(j,k)=aveqq*fac*Habab
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j>0).and.(k<0)) then
        if (j==-k) then
          msq(j,k)=aveqq*fac*(0.5_dp*Hqagg+Hqaqa+real(nf-1,dp)*Hqarb)
          msq_struc(iqr,j,k)=aveqq*fac*real(nf-1,dp)*Hqarb
          msq_struc(iqq_a,j,k)=aveqq*fac*Hqaqa_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Hqaqa_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Hqaqa_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5_dp*Hqagg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5_dp*Hqagg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5_dp*Hqagg_sym
        else
          msq(j,k)=aveqq*fac*Hqbqb
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j<0).and.(k>0)) then
        if (j==-k) then
          msq(j,k)=aveqq*fac*(0.5_dp*Haqgg+Haqaq+real(nf-1,dp)*Haqbr)
          msq_struc(iqr,j,k)=aveqq*fac*real(nf-1,dp)*Haqbr
          msq_struc(iqq_a,j,k)=aveqq*fac*Haqaq_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Haqaq_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Haqaq_i
          msq_struc(igg_ab,j,k)=aveqq*fac*0.5_dp*Haqgg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*0.5_dp*Haqgg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*0.5_dp*Haqgg_sym
        else
          msq(j,k)=aveqq*fac*Hbqbq
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j>0).and.(k==0)) then
        msq(j,0)=aveqg*fac*Hqgqg
        msq_struc(igg_ab,j,0)=aveqg*fac*Hqgqg_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hqgqg_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hqgqg_sym
      endif
      
      if ((j<0).and.(k==0)) then
        msq(j,0)=aveqg*fac*Hagag
        msq_struc(igg_ab,j,0)=aveqg*fac*Hagag_ab
        msq_struc(igg_ba,j,0)=aveqg*fac*Hagag_ba
        msq_struc(igg_sym,j,0)=aveqg*fac*Hagag_sym
      endif

      if ((j==0).and.(k>0)) then
        msq(0,k)=aveqg*fac*Hgqqg
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgqqg_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgqqg_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgqqg_sym
      endif

      if ((j==0).and.(k<0)) then
        msq(0,k)=aveqg*fac*Hgaag
        msq_struc(igg_ab,0,k)=aveqg*fac*Hgaag_ab
        msq_struc(igg_ba,0,k)=aveqg*fac*Hgaag_ba
        msq_struc(igg_sym,0,k)=aveqg*fac*Hgaag_sym
      endif

      if ((j==0).and.(k==0)) then
        msq(0,0)=avegg*fac*(0.5_dp*Hgggg+real(nf,dp)*Hggqa)
        msq_struc(igg_ab,0,0)=avegg*fac*real(nf,dp)*Hggqa_ab
        msq_struc(igg_ba,0,0)=avegg*fac*real(nf,dp)*Hggqa_ba
        msq_struc(igg_sym,0,0)=avegg*fac*real(nf,dp)*Hggqa_sym
        msq_struc(igggg_a,0,0)=avegg*fac*0.5_dp*Hgggg_1256
        msq_struc(igggg_b,0,0)=avegg*fac*0.5_dp*Hgggg_1625
        msq_struc(igggg_c,0,0)=avegg*fac*0.5_dp*Hgggg_1265
      endif
      
      enddo
      enddo

c--- subtraction matrix elements use qa->aq; calculate this and
c--- artificially store it in msq_struc(iqr,0,0), which is not
c--- used for anything else
      call H4qi(1,i5,i6,2,Hqaaq,Hqaaq_a,Hqaaq_b,Hqaaq_i)
      msq_struc(iqr,0,0)=aveqq*fac*Hqaaq
      
      return
      end

 
