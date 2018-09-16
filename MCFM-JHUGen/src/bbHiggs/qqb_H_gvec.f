      subroutine qqb_H_gvec(p,n,in,msq)
C*********************************************************************** 
c     Author: R.K. Ellis                                               *
c     September, 2001.                                                 *
c     Matrix element for H production                                  *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     f(-p1)+f(-p2)--> H(b(p3)+b~(p4))+f(p5)                           * 
C*********************************************************************** 
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'susycoup.f'
      include 'scale.f'
      include 'part.f'
      include 'couple.f'
      integer j,k,in
C--in is the label of the parton dotted with n
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4)
      double precision h1jetn,fac,n(4),propsq,hdecay
      double precision coupsq_eff,ghbb_eff
      double precision mb_eff,massfrun

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      call dotem(5,p,s)

      if (s(3,4) .lt. 4d0*mbsq) return

c--- run mb to appropriate scale
      if (part .eq. 'lord') then
        mb_eff=massfrun(mb_msbar,scale,amz,1)
      else
        mb_eff=massfrun(mb_msbar,scale,amz,2)
      endif
c       mb_eff=mb_msbar
      
      call hbbdecay(p,3,4,hdecay)
      hdecay=hdecay*susycoup**2
      propsq=1d0/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)
c--- The _eff couplings include the running mass
c--- We need to separate these from the factors associated with the
c--- Higgs decay, because the Br. Ratio does not include running mb
      ghbb_eff=dsqrt(esq/xw)*mb_eff/2d0/wmass
      coupsq_eff=susycoup**2*ghbb_eff**2

      fac=CF*xn*gsq*coupsq_eff*propsq*hdecay

      if (in .eq. 1) then
      msq(0,+5)=-fac*aveqg*h1jetn(2,5,1,p,n)
      msq(0,-5)=-fac*aveqg*h1jetn(5,2,1,p,n)
      elseif (in .eq. 2) then
      msq(+5,0)=-fac*aveqg*h1jetn(1,5,2,p,n)
      msq(-5,0)=-fac*aveqg*h1jetn(5,1,2,p,n)
      endif
      return
      end
 
      double precision function h1jetn(j1,j2,j5,p,n)
      implicit none 
C---calculates the amplitude squared for the process 
c   b(p1)+bbar(p2) --> H(b(p3)+b~(p4))+g(p5)
c   contracted with the vector n(mu) 
c   before spin/color average
      include 'constants.f'
      include 'sprods_com.f'

      integer j1,j2,j3,j4,j5
      double precision n(4),p(mxpart,4),nDn,nDp1,nDp2
      j3=3
      j4=4

      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

      call checkndotp(p,n,j5)

      h1jetn=4d0*(  
     . +2d0*nDp2**2*(s(j1,j2)+s(j1,j5))/s(j2,j5)**2
     . +2d0*nDp1**2*(s(j1,j2)+s(j2,j5))/s(j1,j5)**2 
     . +(2d0*s(j2,j5)*nDp1**2+2d0*s(j1,j5)*nDp2**2
     .  -nDn/2d0*s(j1,j5)**2-nDn/2d0*s(j2,j5)**2 
     . -4d0*nDp1*nDp2
     . *(s(j1,j2)+s(j1,j5)+s(j2,j5)))/s(j1,j5)/s(j2,j5)- nDn)

      return
      end




