      subroutine qq_tchan_htqg(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->h(p3,p4)+t(p5)+d(p6)+g(p7)
      
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
      integer:: j,nu,eta,e5
      real(dp):: p(mxpart,4),q(mxpart,4),p5Deta,s34,hdecay,
     & msq(-nf:nf,-nf:nf),dot,msqgamgam,fac,ampsqL,ampsqH,
     & b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub,
     & db_g,g_b,u_g,g_db,b_g,g_u,
     & d_g,g_bb,ub_g,g_d,bb_g,g_ub

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

C   Construct demassified momentum for p5 wrt eta and store in position e5
C   Set eta=3 for now
      eta=3
      p5Deta=dot(p,5,eta)

      e5=8
      do nu=1,4
      do j=1,7
      q(j,nu)=p(j,nu)
      enddo
      q(e5,nu)=mt**2*p(eta,nu)/(2d0*p5Deta)
      q(5,nu)=p(5,nu)-q(e5,nu)
      enddo

C---zero all matrix elements
      msq(:,:)=0d0
      fac=xn**2*(2d0*CF)*gwsq**3*gsq*hdecay

      if (nwz == 1) then
      call spinoru(8,q,za,zb) 

      call ubhtdgsq(1,2,3,4,5,6,7,e5,ampsqL,ampsqH)
      u_b=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(2,1,3,4,5,6,7,e5,ampsqL,ampsqH)
      b_u=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(6,2,3,4,5,1,7,e5,ampsqL,ampsqH)
      db_b=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(6,1,3,4,5,2,7,e5,ampsqL,ampsqH)
      b_db=aveqq*fac*(ampsqL+ampsqH)
c--- In g-b diagrams, remove corrections on heavy line
      call ubhtdgsq(7,2,3,4,5,6,1,e5,ampsqL,ampsqH)
      g_b=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
      call ubhtdgsq(7,1,3,4,5,6,2,e5,ampsqL,ampsqH)
      b_g=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
c--- In g-q diagrams, remove corrections on light line
      call ubhtdgsq(2,7,3,4,5,6,1,e5,ampsqL,ampsqH)
      g_u=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(1,7,3,4,5,6,2,e5,ampsqL,ampsqH)
      u_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(6,7,3,4,5,1,2,e5,ampsqL,ampsqH)
      db_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(6,7,3,4,5,2,1,e5,ampsqL,ampsqH)
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
      
      call ubhtdgsq(1,2,4,3,5,6,7,e5,ampsqL,ampsqH)
      d_bb=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(2,1,4,3,5,6,7,e5,ampsqL,ampsqH)
      bb_d=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(6,2,4,3,5,1,7,e5,ampsqL,ampsqH)
      ub_bb=aveqq*fac*(ampsqL+ampsqH)
      call ubhtdgsq(6,1,4,3,5,2,7,e5,ampsqL,ampsqH)
      bb_ub=aveqq*fac*(ampsqL+ampsqH)

c--- In g-b diagrams, remove corrections on heavy line
      call ubhtdgsq(7,2,4,3,5,6,1,e5,ampsqL,ampsqH)
      g_bb=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
      call ubhtdgsq(7,1,4,3,5,6,2,e5,ampsqL,ampsqH)
      bb_g=2d0*aveqg*fac*(ampsqL+zip*ampsqH)
c--- In g-q diagrams, remove corrections on light line
      call ubhtdgsq(2,7,4,3,5,6,1,e5,ampsqL,ampsqH)
      g_d=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(1,7,4,3,5,6,2,e5,ampsqL,ampsqH)
      d_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(6,7,4,3,5,1,2,e5,ampsqL,ampsqH)
      ub_g=aveqg*fac*(zip*ampsqL+ampsqH)
      call ubhtdgsq(6,7,4,3,5,2,1,e5,ampsqL,ampsqH)
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
