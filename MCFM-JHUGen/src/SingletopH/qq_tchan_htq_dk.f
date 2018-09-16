      subroutine qq_tchan_htq_dk(p,msq)
      implicit none
      include 'types.f'
      
c---Matrix element squared averaged over initial colors and spins
c     u(-p1)+b(p2)->h(p3,p4)+t(nu(p5)+e(p6)+b(p7))+d(p8)
      include 'ewcouple.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'nwz.f'
      include 'zprods_com.f'
      integer:: nu,eta,e5
      real(dp):: p(mxpart,4),q(mxpart,4),q5Deta,s34,hdecay,
     & msq(-nf:nf,-nf:nf),msqgamgam,fac,tpropsq,
     & b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub,ubhtdsqdk
       complex(dp):: mdecaymb(2,2),mdecay


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

c --  Deal with top decay
c --  tdecay assumes a massive b, hence 4 polarizations -- we only need one
      mb=0d0
      if (nwz == 1) then
      call tdecay(p,5,6,7,mdecaymb)
      mdecay=mdecaymb(1,1)
      elseif (nwz == -1) then
      call adecay(p,5,6,7,mdecaymb)
      mdecay=mdecaymb(2,2)
      endif
c -- now assign momenta for undecayed top 
      q(1,:)=p(1,:)
      q(2,:)=p(2,:)
      q(3,:)=p(3,:)
      q(4,:)=p(4,:)
      q(5,:)=p(5,:)+p(6,:)+p(7,:)
      q(6,:)=p(8,:)


      eta=6
      q5Deta=q(5,4)*p(eta,4)
     & -q(5,1)*p(eta,1)-q(5,2)*p(eta,2)-q(5,3)*p(eta,3)

      e5=7
      do nu=1,4
C   store rescaled electron momentum in position e5
      q(e5,nu)=mt**2*p(eta,nu)/(2d0*q5Deta)
C   Construct demassified momentum for q(5,:) 
C   using rescaled electron momentum q(e5,:) 
      q(5,nu)=q(5,nu)-q(e5,nu)
      enddo
      tpropsq=1d0/(mt*twidth)**2


      fac=xn**2*aveqq*gwsq**5*hdecay*tpropsq

C -- zero all components of msq
      msq(:,:)=0d0

      if (nwz == +1) then
      call spinoru(7,q,za,zb) 

      u_b=ubhtdsqdk(1,2,3,4,5,6,e5,mdecay)*fac
      db_b=ubhtdsqdk(6,2,3,4,5,1,e5,mdecay)*fac
      b_u=ubhtdsqdk(2,1,3,4,5,6,e5,mdecay)*fac
      b_db=ubhtdsqdk(6,1,3,4,5,2,e5,mdecay)*fac

      msq(-1,+5)=db_b
      msq(-3,+5)=db_b
      msq(+2,+5)=u_b
      msq(+4,+5)=u_b
      msq(+5,-1)=b_db
      msq(+5,-3)=b_db
      msq(+5,+2)=b_u
      msq(+5,+4)=b_u

      
      elseif(nwz == -1) then
      call spinoru(7,q,zb,za)

      d_bb=ubhtdsqdk(1,2,4,3,5,6,e5,mdecay)*fac
      ub_bb=ubhtdsqdk(6,2,4,3,5,1,e5,mdecay)*fac
      bb_d=ubhtdsqdk(2,1,4,3,5,6,e5,mdecay)*fac
      bb_ub=ubhtdsqdk(6,1,4,3,5,2,e5,mdecay)*fac


      msq(+1,-5)=d_bb
      msq(+3,-5)=d_bb
      msq(-2,-5)=ub_bb
      msq(-4,-5)=ub_bb

      msq(-5,+1)=bb_d
      msq(-5,+3)=bb_d
      msq(-5,-2)=bb_ub
      msq(-5,-4)=bb_ub
      
      endif
      return
      end
