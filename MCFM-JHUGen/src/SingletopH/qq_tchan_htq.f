      subroutine qq_tchan_htq(p,msq)
      implicit none
      include 'types.f'
      
c---  Matrix element squared averaged over initial colors and spins
c---  u(-p1)+b(p2)->H(p3,p4)+t(p5)+d(p6) etc
      include 'ewcouple.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'nwz.f'
      include 'zprods_com.f'
      integer:: j,nu,eta,e5
      real(dp):: p(mxpart,4),q(mxpart,4),p5Deta,s34,hdecay,
     & msq(-nf:nf,-nf:nf),dot,msqgamgam,fac,
     & b_u,u_b,db_b,b_db,d_bb,ub_bb,bb_d,bb_ub,ubhtdsq


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

      e5=7
      do nu=1,4
      do j=1,6
      q(j,nu)=p(j,nu)
      enddo
      q(e5,nu)=mt**2*p(eta,nu)/(2d0*p5Deta)
      q(5,nu)=p(5,nu)-q(e5,nu)
      enddo

C -- zero all matrix elements
      msq(:,:)=0d0
      
      fac=xn**2*aveqq*gwsq**3*hdecay

      if (nwz == 1) then
        call spinoru(7,q,za,zb) 

        u_b=ubhtdsq(1,2,3,4,5,6,e5)*fac
        db_b=ubhtdsq(6,2,3,4,5,1,e5)*fac
        b_u=ubhtdsq(2,1,3,4,5,6,e5)*fac
        b_db=ubhtdsq(6,1,3,4,5,2,e5)*fac

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

        d_bb=ubhtdsq(1,2,4,3,5,6,e5)*fac
        ub_bb=ubhtdsq(6,2,4,3,5,1,e5)*fac
        bb_d=ubhtdsq(2,1,4,3,5,6,e5)*fac       
        bb_ub=ubhtdsq(6,1,4,3,5,2,e5)*fac
               
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
