      subroutine gg_ww(pin,msqgg)
      implicit none
      include 'types.f'

c--- Author: J. M. Campbell, March 2011
c--- For now, work in the approximation of two massless isodoublets
c--- Box contributions are then complete
c--- Triangle (vector) pieces always vanish
c--- Triangle (axial) pieces cancel for massless isodoublets

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'Higgsint.f'
      include 'noglue.f'
      include 'plabel.f'
      include 'first.f'
      integer:: j,h1,h2,nu,del1,del2,k12h,k34h,k56h11,k34h11,e
      real(dp):: p(mxpart,4),pin(mxpart,4),msqgg,fac
      complex(dp):: Avec(2,2),Agen3(2,2),sum(2,2,-2:0),
     & box(2,2,-2:0),triang(2,2,-2:0),bub(2,2,-2:0),faccont
      complex(dp):: a64v
      real(dp):: dot,s12,s34,s56,dot1256,afac,bfac,gden,delta,
     & dot1234,dot3456,pttwo,ptWsafetycut_massive,ptWsafetycut_massless
      logical:: includegens1and2,includegen3
      parameter(del1=7,del2=8)
      parameter(k12h=9,k34h=10,k56h11=11,k34h11=12)

c--- if noglue or omitgg, set to zero and return
      if (noglue .or. omitgg) then
         msqgg=0._dp
       return
      endif

      if ((plabel(3) == 'el') .and. (plabel(5) == 'qj')) then
C----swapped case for hadronic decay of Wp
      do j=1,4
      p(1,j)=pin(1,j)
      p(2,j)=pin(2,j)
      p(3,j)=pin(5,j)
      p(4,j)=pin(6,j)
      p(5,j)=pin(3,j)
      p(6,j)=pin(4,j)
      enddo
      else
      do j=1,4
      p(1,j)=pin(1,j)
      p(2,j)=pin(2,j)
      p(3,j)=pin(3,j)
      p(4,j)=pin(4,j)
      p(5,j)=pin(5,j)
      p(6,j)=pin(6,j)
      enddo
      endif

c--- set this to true to include generations 1 and 2 of (light) quarks
      includegens1and2=.true.
c--- set this to true to include 3rd generation of quarks (t,b)
      includegen3=.true.

c--- omit massive loops for pt(W) < "ptWsafetycut_massive" (for num. stability)
      ptWsafetycut_massive=2._dp

c--- omit massless loops for pt(W) < "ptWsafetycut_massless" (for num. stability)
      ptWsafetycut_massless=0.05_dp

      if (first) then
        write(6,*)'****************************************************'
        write(6,*)'*                                                  *'
        if (includegens1and2) then
        write(6,*)'*  gg->WW box loop includes gens. 1 and 2          *'
        else
        write(6,*)'*  gg->WW box loop does not include gens. 1 and 2  *'
        endif
        if (includegen3) then
        write(6,*)'*  gg->WW box loop includes 3rd generation         *'
        else
        write(6,*)'*  gg->WW box loop does not include 3rd gen.       *'
        endif
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        first=.false.
      if (includegen3) call qlinit
      endif

c--- if neither contribution is included, set to zero and return
      if ((includegens1and2 .eqv. .false.) .and.
     &    (includegen3      .eqv. .false.)) then
         msqgg=0._dp
       return
      endif

c--- if computing 3rd generation, set up extra flattened vectors
      if (includegen3) then
C--- define flattened vectors (k12 and k56)
      s12=2._dp*dot(p,1,2)
      s56=2._dp*dot(p,5,6)
      s34=2._dp*dot(p,3,4)
      dot1256=0.5_dp*(s34-s12-s56)
      delta=dot1256**2-s12*s56
      gden=dot1256+sqrt(delta)
      afac=s12/gden
      bfac=s56/gden
      do nu=1,4
      p(del2,nu)=one/(one-afac*bfac)
     & *(p(1,nu)+p(2,nu)-afac*(p(5,nu)+p(6,nu)))
      p(del1,nu)=one/(one-afac*bfac)
     & *(p(5,nu)+p(6,nu)-bfac*(p(1,nu)+p(2,nu)))
      enddo

C--- define flattened vectors (k12 and k34)
      dot1234=0.5_dp*(s56-s12-s34)
      delta=dot1234**2-s12*s34
      gden=dot1234+sqrt(delta)
      afac=s12/gden
      bfac=s34/gden

      do nu=1,4
      p(k12h,nu)=one/(one-afac*bfac)
     & *(p(1,nu)+p(2,nu)-afac*(p(3,nu)+p(4,nu)))
      p(k34h,nu)=one/(one-afac*bfac)
     & *(p(3,nu)+p(4,nu)-bfac*(p(1,nu)+p(2,nu)))
      enddo

C--- define flattened vectors (k56 and k34)
      dot3456=0.5_dp*(s12-s56-s34)
      delta=dot3456**2-s56*s34
      gden=dot3456+sqrt(delta)
      afac=s56/gden
      bfac=s34/gden

      do nu=1,4
      p(k56h11,nu)=one/(one-afac*bfac)
     & *(p(5,nu)+p(6,nu)-afac*(p(3,nu)+p(4,nu)))
      p(k34h11,nu)=one/(one-afac*bfac)
     & *(p(3,nu)+p(4,nu)-bfac*(p(5,nu)+p(6,nu)))
      enddo

      endif
c--- end of 3rd generation initialization

c--- set up spinor products (including for flat vectors, for 3rd gen)
      if (includegen3) then
        call spinoru(12,p,za,zb)
      else
        call spinoru(6,p,za,zb)
      endif

c--- fill amplitudes used for generations 1 and 2
      if (includegens1and2) then
c--- ensure numerical stability: set massless loops to zero
c--- for pt(W) < "ptWsafetycut_massless" GeV
        if (pttwo(3,4,p) < ptWsafetycut_massless) then
          do h1=1,2
          do h2=1,2
            Avec(h1,h2)=czip
          enddo
          enddo
      else
          Avec(2,2)=a64v('q+qb-g-g-',3,4,1,2,6,5,zb,za)*(-im)
          Avec(2,1)=a64v('q+qb-g-g+',3,4,1,2,6,5,zb,za)*(-im)
          Avec(1,2)=a64v('q+qb-g+g-',3,4,1,2,6,5,zb,za)*(-im)
          Avec(1,1)=a64v('q+qb-g+g+',3,4,1,2,6,5,zb,za)*(-im)
        endif
      else
        do h1=1,2
      do h2=1,2
      Avec(h1,h2)=czip
      enddo
      enddo
      endif

c--- fill amplitudes used for 3rd generation
      if (includegen3) then
c--- ensure numerical stability: set massive loops to zero
c--- for pt(W) < "ptWsafetycut_massive" GeV
        if ((pttwo(3,4,p)/sqrt(s(1,2)) < 1.e-2_dp)
     &  .or.(pttwo(3,4,p) < ptWsafetycut_massive)) then
          do h1=1,2
          do h2=1,2
            Agen3(h1,h2)=czip
          enddo
          enddo
        else
          Higgsint=.false.
          do h1=1,2
          do h2=1,2
          do e=-2,0
          box(h1,h2,e)=czip
          triang(h1,h2,e)=czip
          bub(h1,h2,e)=czip
          enddo
          enddo
          enddo
c---    compute integrals and their coefficients
          call massivebox6(1,2,3,4,5,6,za,zb,box)
          call massivetri6(1,2,3,4,5,6,za,zb,triang)
          call massivebub(1,2,3,4,5,6,za,zb,bub)
c---    this contribution is finite, so we only retain "0" piece
          e=0
          do h1=1,2
          do h2=1,2
           sum(h1,h2,e)=box(h1,h2,e)+triang(h1,h2,e)+bub(h1,h2,e)
           Agen3(h1,h2)=sum(h1,h2,0)
          enddo
          enddo
      endif
      else
        do h1=1,2
      do h2=1,2
      Agen3(h1,h2)=czip
      enddo
      enddo
      endif

c--- factor two for complete massless isodoublets
      faccont=ctwo

c--- sum and square amplitudes
      msqgg=0._dp
      do h1=1,2
      do h2=1,2
      msqgg=msqgg+abs(faccont*Avec(h1,h2)+Agen3(h1,h2))**2
      enddo
      enddo

c--- overall factor from diagrams
      fac=avegg*V*(2._dp*gwsq*gsq/(16._dp*pisq)*gwsq/2._dp)**2
     & *s(3,4)**2/((s(3,4)-wmass**2)**2+(wwidth*wmass)**2)
     & *s(5,6)**2/((s(5,6)-wmass**2)**2+(wwidth*wmass)**2)

C-----add in factor for hadronic decays
      if (plabel(5) == 'qj') fac=2._dp*xn*fac

      msqgg=msqgg*fac

      return
      end


