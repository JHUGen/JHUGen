      subroutine gg_hgg(p,msq)
      implicit none
      include 'types.f'

c---Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p3)+g(p_iglue1=5)+g(p_iglue2=6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'msq_struc.f'
      include 'nflav.f'
      include 'hdecaymode.f'
      integer:: j,k,iglue1,iglue2
      real(dp):: p(mxpart,4),Asq,fac,msqhgamgam
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
     & Hqaaq
      real(dp):: msq(-nf:nf,-nf:nf),hdecay,s34

      parameter(iglue1=5,iglue2=6)


C---fill spinor products up to maximum number
      call spinoru(iglue2,p,za,zb)

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqhgamgam(s34)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(three*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

c--- four gluon terms
      call HggggLO(1,2,iglue1,iglue2,
     &             Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)

c--- two quark two gluon terms
      call HQAggLO(1,2,iglue1,iglue2,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)
      call HQAggLO(2,1,iglue1,iglue2,Haqgg,Haqgg_ab,Haqgg_ba,Haqgg_sym)
c---   note: symmetric in first two arguments, but not the ab, ba terms

      call HQAggLO(1,iglue1,2,iglue2,Hqgqg,Hqgqg_ab,Hqgqg_ba,Hqgqg_sym)
      call HQAggLO(iglue1,1,2,iglue2,Hagag,Hagag_ab,Hagag_ba,Hagag_sym)

      call HQAggLO(2,iglue1,1,iglue2,Hgqqg,Hgqqg_ab,Hgqqg_ba,Hgqqg_sym)
      call HQAggLO(iglue1,2,1,iglue2,Hgaag,Hgaag_ab,Hgaag_ba,Hgaag_sym)

      call HQAggLO(iglue2,iglue1,1,2,Hggqa,Hggqa_ab,Hggqa_ba,Hggqa_sym)

c--- four quark terms
      call HqarbLO(1,2,iglue1,iglue2,Hqrqr)
      call HqaqaLO(1,2,iglue1,iglue2,Hqqqq,Hqqqq_a,Hqqqq_b,Hqqqq_i)

C---four anti-quark terms
c      call H4qn(iglue1,iglue2,1,2,Habab)
c      call H4qi(iglue1,iglue2,1,2,Haaaa)
      Habab=Hqrqr
      Haaaa=Hqqqq

C-qqb
      call HqarbLO(1,iglue2,2,iglue1,Hqarb)
      call HqaqaLO(1,iglue2,iglue1,2,Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i)
      call HqarbLO(1,iglue2,iglue1,2,Hqbqb)

C-qbq
      Haqbr=Hqarb

      Haqaq=Hqaqa
      Haqaq_a=Hqaqa_a
      Haqaq_b=Hqaqa_b
      Haqaq_i=Hqaqa_i
      Hbqbq=Hqbqb

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      msq_struc(:,j,k)=zip

      if ((j>0).and.(k>0)) then
        if (j==k) then
          msq(j,k)=half*aveqq*fac*Hqqqq
          msq_struc(iqq_a,j,k)=half*aveqq*fac*Hqqqq_a
          msq_struc(iqq_b,j,k)=half*aveqq*fac*Hqqqq_b
          msq_struc(iqq_i,j,k)=half*aveqq*fac*Hqqqq_i
        else
          msq(j,k)=aveqq*fac*Hqrqr
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j<0).and.(k<0)) then
        if (j==k) then
          msq(j,k)=half*aveqq*fac*Haaaa
        else
          msq(j,k)=aveqq*fac*Habab
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j>0).and.(k<0)) then
        if (j==-k) then
          msq(j,k)=aveqq*fac*(half*Hqagg+Hqaqa+real(nflav-1,dp)*Hqarb)
          msq_struc(iqr,j,k)=aveqq*fac*real(nflav-1,dp)*Hqarb
          msq_struc(iqq_a,j,k)=aveqq*fac*Hqaqa_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Hqaqa_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Hqaqa_i
          msq_struc(igg_ab,j,k)=aveqq*fac*half*Hqagg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*half*Hqagg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*half*Hqagg_sym
        else
          msq(j,k)=aveqq*fac*Hqbqb
          msq_struc(iqq_a,j,k)=msq(j,k)
          msq_struc(iqq_b,j,k)=0._dp
          msq_struc(iqq_i,j,k)=0._dp
        endif
      endif

      if ((j<0).and.(k>0)) then
        if (j==-k) then
          msq(j,k)=aveqq*fac*(half*Haqgg+Haqaq+real(nflav-1,dp)*Haqbr)
          msq_struc(iqr,j,k)=aveqq*fac*real(nflav-1,dp)*Haqbr
          msq_struc(iqq_a,j,k)=aveqq*fac*Haqaq_a
          msq_struc(iqq_b,j,k)=aveqq*fac*Haqaq_b
          msq_struc(iqq_i,j,k)=aveqq*fac*Haqaq_i
          msq_struc(igg_ab,j,k)=aveqq*fac*half*Haqgg_ab
          msq_struc(igg_ba,j,k)=aveqq*fac*half*Haqgg_ba
          msq_struc(igg_sym,j,k)=aveqq*fac*half*Haqgg_sym
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
        msq(0,0)=avegg*fac*(half*Hgggg+real(nflav,dp)*Hggqa)
        msq_struc(igg_ab,0,0)=avegg*fac*real(nflav,dp)*Hggqa_ab
        msq_struc(igg_ba,0,0)=avegg*fac*real(nflav,dp)*Hggqa_ba
        msq_struc(igg_sym,0,0)=avegg*fac*real(nflav,dp)*Hggqa_sym
        msq_struc(igggg_a,0,0)=avegg*fac*half*Hgggg_1256
        msq_struc(igggg_b,0,0)=avegg*fac*half*Hgggg_1625
        msq_struc(igggg_c,0,0)=avegg*fac*half*Hgggg_1265
      endif

      enddo
      enddo

c--- subtraction matrix elements use qa->aq; calculate this and
c--- artificially store it in msq_struc(iqr,0,0), which is not
c--- used for anything else
      call H4qi(1,iglue1,iglue2,2,Hqaaq,Hqaaq_a,Hqaaq_b,Hqaaq_i)
      msq_struc(iqr,0,0)=aveqq*fac*Hqaaq

      return
      end


