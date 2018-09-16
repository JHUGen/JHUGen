      subroutine qq_tchan_htqg_dk(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->h(p3,p4)+t(nu(p5)+e(p6)+b(p7))+d(p8)+g(p9)
      
      include 'ewcouple.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'nwz.f'
      include 'zprods_com.f'
      integer:: nu,eta,e5
      real(dp):: p(mxpart,4),q(mxpart,4),q5Deta,s34,hdecay,
     & msq(-nf:nf,-nf:nf),msqgamgam,fac,ampsqL,ampsqH,
     & b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub,
     & db_g,g_b,u_g,g_db,b_g,g_u,
     & d_g,g_bb,ub_g,g_d,bb_g,g_ub,tpropsq
      complex(dp):: mdecaymb(2,2),mdecay

c -- mb set to zero for rest of code, needs to be zero here too, so this is to be sure
c -- tdecay assumes a massive b, hence 4 polarizations -- we only need one
      mb=0d0
      if (nwz == 1) then
      call tdecay(p,5,6,7,mdecaymb)
      mdecay=mdecaymb(1,1)
      elseif (nwz == -1) then
      call adecay(p,5,6,7,mdecaymb)
      mdecay=mdecaymb(2,2)
      endif

C   Deal with Higgs decay
      s34=(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &   -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

c -- now assign momenta q for undecayed top 
      q(1,:)=p(1,:)
      q(2,:)=p(2,:)
      q(3,:)=p(3,:)
      q(4,:)=p(4,:)
      q(5,:)=p(5,:)+p(6,:)+p(7,:)
      q(6,:)=p(8,:)
      q(7,:)=p(9,:)

C   Construct demassified momentum for q5 wrt electron momentum and store in position e5=8
      e5=8
      eta=6
      q5Deta=q(5,4)*p(eta,4)
     & -q(5,1)*p(eta,1)-q(5,2)*p(eta,2)-q(5,3)*p(eta,3)

      do nu=1,4
      q(e5,nu)=mt**2*p(eta,nu)/(2d0*q5Deta)
      q(5,nu)=q(5,nu)-q(e5,nu)
      enddo

C -- zero all components of msq
      msq(:,:)=0d0

      tpropsq=1d0/(mt*twidth)**2
      fac=xn**2*(2d0*CF)*gwsq**5*gsq*hdecay*tpropsq
      if (nwz == 1) then
      call spinoru(8,q,za,zb) 

      call ubhtdgsqdk(1,2,3,4,5,6,7,e5,mdecay,ampsqL,ampsqH)
      u_b=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(2,1,3,4,5,6,7,e5,mdecay,ampsqL,ampsqH)
      b_u=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(6,2,3,4,5,1,7,e5,mdecay,ampsqL,ampsqH)
      db_b=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(6,1,3,4,5,2,7,e5,mdecay,ampsqL,ampsqH)
      b_db=aveqq*fac*(ampsqL+ampsqH)
c--- In g-b diagrams, remove corrections on heavy line
      call ubhtdgsqdk(7,2,3,4,5,6,1,e5,mdecay,ampsqL,ampsqH)
      g_b=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
      call ubhtdgsqdk(7,1,3,4,5,6,2,e5,mdecay,ampsqL,ampsqH)
      b_g=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
c--- In g-q diagrams, remove corrections on light line
      call ubhtdgsqdk(2,7,3,4,5,6,1,e5,mdecay,ampsqL,ampsqH)
      g_u=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(1,7,3,4,5,6,2,e5,mdecay,ampsqL,ampsqH)
      u_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(6,7,3,4,5,1,2,e5,mdecay,ampsqL,ampsqH)
      db_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(6,7,3,4,5,2,1,e5,mdecay,ampsqL,ampsqH)
      g_db=aveqg*fac*(zip*ampsqL+ampsqH)

      msq(-1,+5)=db_b
      msq(-3,+5)=db_b
      msq(+2,+5)=u_b
      msq(+4,+5)=u_b

      msq(+5,-1)=b_db
      msq(+5,-3)=b_db
      msq(+5,+2)=b_u
      msq(+5,+4)=b_u

      msq(+2,0)=u_g
      msq(+4,0)=u_g
      msq(-1,0)=db_g
      msq(-3,0)=db_g
      msq( 0,5)=g_b

      msq(0,+2)=g_u
      msq(0,+4)=g_u
      msq(0,-1)=g_db
      msq(0,-3)=g_db
      msq(5,0)=b_g

      elseif(nwz == -1) then
      call spinoru(8,q,zb,za) 

      call ubhtdgsqdk(1,2,4,3,5,6,7,e5,mdecay,ampsqL,ampsqH)
      d_bb=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(2,1,4,3,5,6,7,e5,mdecay,ampsqL,ampsqH)
      bb_d=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(6,2,4,3,5,1,7,e5,mdecay,ampsqL,ampsqH)
      ub_bb=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsqdk(6,1,4,3,5,2,7,e5,mdecay,ampsqL,ampsqH)
      bb_ub=aveqq*fac*(ampsqL+ampsqH)

c--- In g-b diagrams, remove corrections on heavy line
      call ubhtdgsqdk(7,2,4,3,5,6,1,e5,mdecay,ampsqL,ampsqH)
      g_bb=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
      call ubhtdgsqdk(7,1,4,3,5,6,2,e5,mdecay,ampsqL,ampsqH)
      bb_g=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
c--- In g-q diagrams, remove corrections on light line
      call ubhtdgsqdk(2,7,4,3,5,6,1,e5,mdecay,ampsqL,ampsqH)
      g_d=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(1,7,4,3,5,6,2,e5,mdecay,ampsqL,ampsqH)
      d_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(6,7,4,3,5,1,2,e5,mdecay,ampsqL,ampsqH)
      ub_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsqdk(6,7,4,3,5,2,1,e5,mdecay,ampsqL,ampsqH)
      g_ub=aveqg*fac*(zip*ampsqL+ampsqH)


      msq(+1,-5)=d_bb
      msq(+3,-5)=d_bb
      msq(-2,-5)=ub_bb
      msq(-4,-5)=ub_bb

      msq(-5,+1)=bb_d
      msq(-5,+3)=bb_d
      msq(-5,-2)=bb_ub
      msq(-5,-4)=bb_ub
      
      msq(+1,0)=d_g
      msq(+3,0)=d_g
      msq(-2,0)=ub_g
      msq(-4,0)=ub_g
      msq(0,-5)=g_bb

      msq(0,+1)=g_d
      msq(0,+3)=g_d
      msq(0,-2)=g_ub
      msq(0,-4)=g_ub
      msq(-5,0)=bb_g

      endif
      return
      end
