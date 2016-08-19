      subroutine gg_hggg(p,msq)
      implicit none
c--- Matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2) -->  H(p34)+g(p_iglue1=5)+g(p_iglue2=6)+g(p_iglue2=7) 
c
c--- Using the results of Frizzo and Company
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'nflav.f'
      include 'hdecaymode.f'
      include 'bitflags.f'
      integer j,k,nu
      double precision p(mxpart,4),Asq,fac,q(mxpart,4)
      double precision Hggggg,msqgamgam,
     . Hqaggg,Haqggg,Hgqqgg,Hgaagg,Hqgqgg,Hagagg,Hggqag
      double precision qr_qrg,ar_arg,ab_abg,qa_rbg,
     .                 gr_rqa,gb_baq,rg_rqa,bg_baq,aq_brg

      double precision qq_qqg,aq_aqg,aa_aag,
     .                 gq_qqa,ga_aaq,qg_qqa,ag_aaq
      double precision ra_rag,qa_qag

      double precision dummy
       
      double precision msq(-nf:nf,-nf:nf),hdecay,s34
      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in gg_hgg_v'
      stop
      endif

      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(3d0*pi))**2/vevsq

C---swap momenta so that Higgs decay products are last
      do nu=1,4
      q(1,nu)=p(1,nu)
      q(2,nu)=p(2,nu)
      q(3,nu)=p(5,nu)
      q(4,nu)=p(6,nu)
      q(5,nu)=p(7,nu)
      q(6,nu)=p(3,nu)
      q(7,nu)=p(4,nu)
      enddo

C---fill spinor products up to maximum number
      call spinoru(5,q,za,zb)  


C--five gluon terms
      call h5g(Hggggg)

C--two quark three gluon terms
      call h2q3g(1,2,3,4,5,Hqaggg)
      call h2q3g(2,1,3,4,5,Haqggg)

      call h2q3g(1,3,2,4,5,Hqgqgg)
      call h2q3g(2,3,1,4,5,Hgqqgg)

      call h2q3g(3,1,2,4,5,Hagagg)
      call h2q3g(3,2,1,4,5,Hgaagg)
      call h2q3g(4,3,1,2,5,Hggqag)

C--four quark one gluon terms
C-----q r-->q r g
      call h4qg(3,1,4,2,5,qr_qrg,qq_qqg)

C---  q~ r --> q~ r g
      call h4qg(4,2,1,3,5,ar_arg,aq_aqg)

C---  q r~ --> q r~ g
      call h4qg(3,1,2,4,5,ra_rag,qa_qag)
      
C---  q~ r~--> q~ r~ g
      call h4qg(4,2,3,1,5,ab_abg,aa_aag)

C---  q q~ -> r r~ g (note that dummy is the same as qa_qag)
      call h4qg(2,1,3,4,5,qa_rbg,dummy)

C---  q~ q -> r~ r g (note that dummy is the same as aq_aqg)
      call h4qg(1,2,4,3,5,aq_brg,dummy)
      
C---  g r --> r q q~
      call h4qg(3,2,4,5,1,gr_rqa,gq_qqa)

C---  g r~ --> r~ q~ q
      call h4qg(2,3,5,4,1,gb_baq,ga_aaq)

C---  r g --> r q q~
      call h4qg(3,1,4,5,2,rg_rqa,qg_qqa)

C---  r~ g --> r~ q~ q
      call h4qg(1,3,5,4,2,bg_baq,ag_aaq)

      fac=gsq**3*Asq*hdecay

c--- apply flags
      Hggggg=f0q*Hggggg
      
      Hqaggg=f2q*Hqaggg
      Haqggg=f2q*Haqggg
      Hqgqgg=f2q*Hqgqgg
      Hgqqgg=f2q*Hgqqgg

      Hagagg=f2q*Hagagg
      Hgaagg=f2q*Hgaagg
      Hggqag=f2q*Hggqag

      qr_qrg=f4q*qr_qrg
      ar_arg=f4q*ar_arg
      ra_rag=f4q*ra_rag
      ab_abg=f4q*ab_abg
      gr_rqa=f4q*gr_rqa
      gb_baq=f4q*gb_baq
      rg_rqa=f4q*rg_rqa
      bg_baq=f4q*bg_baq

      qq_qqg=f4q*qq_qqg
      aq_aqg=f4q*aq_aqg
      qa_qag=f4q*qa_qag
      aa_aag=f4q*aa_aag
      gq_qqa=f4q*gq_qqa
      ga_aaq=f4q*ga_aaq
      qg_qqa=f4q*qg_qqa
      ag_aaq=f4q*ag_aaq

      qa_qag=f4q*qa_qag
      ra_rag=f4q*ra_rag
      qa_rbg=f4q*qa_rbg
      aq_brg=f4q*aq_brg

C----Fill up array with values;
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

C ---qq
      if ((j.gt.0).and.(k.gt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*qq_qqg
        else
          msq(j,k)=aveqq*fac*qr_qrg
        endif
      endif

C ---aa
      if ((j.lt.0).and.(k.lt.0)) then 
        if (j.eq.k) then
          msq(j,k)=0.5d0*aveqq*fac*aa_aag
        else
          msq(j,k)=aveqq*fac*ab_abg
        endif
      endif

c ---qa
      if ((j.gt.0).and.(k.lt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(Hqaggg/6d0+qa_qag+dfloat(nflav-1)*qa_rbg)
        else
          msq(j,k)=aveqq*fac*ra_rag
        endif
      endif

c ---aq
      if ((j.lt.0).and.(k.gt.0)) then
        if (j.eq.-k) then
          msq(j,k)=aveqq*fac*(Haqggg/6d0+aq_aqg+dfloat(nflav-1)*aq_brg)
        else
          msq(j,k)=aveqq*fac*ar_arg
        endif
      endif

c--- qg
      if ((j.gt.0).and.(k.eq.0)) then
       msq(j,0)=aveqg*fac*((Hqgqgg+qg_qqa)*0.5d0+dfloat(nflav-1)*rg_rqa)
      endif
      
c--- ag
      if ((j.lt.0).and.(k.eq.0)) then
       msq(j,0)=aveqg*fac*((Hagagg+ag_aaq)*0.5d0+dfloat(nflav-1)*bg_baq)
      endif

c--- gq
      if ((j.eq.0).and.(k.gt.0)) then
       msq(0,k)=aveqg*fac*((Hgqqgg+gq_qqa)*0.5d0+dfloat(nflav-1)*gr_rqa)
      endif

c--- ga
      if ((j.eq.0).and.(k.lt.0)) then
       msq(0,k)=aveqg*fac*((Hgaagg+ga_aaq)*0.5d0+dfloat(nflav-1)*gb_baq)
      endif

c--- gg
      if ((j.eq.0).and.(k.eq.0)) then
        msq(0,0)=avegg*fac*(Hggggg/6d0+dfloat(nflav)*Hggqag)
      endif
      
      enddo
      enddo

      return
      end

 
