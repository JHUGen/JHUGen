      subroutine qqb_tth(pin,msq)
      implicit none
      include 'types.f'

************************************************************************
*     Author: R.K. Ellis                                               *
*     December, 1999.                                                  *
*     calculate the element squared and subtraction terms              *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=tbar(bbar(p6)+e-(p7)+nubar(p8))
C                      +t(b(p5)+nu(p3)+e+(p4))
C                      +H(b(p9)+bbar(p10))
C
************************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'zprods_com.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'masses.f'
      include 'plabel.f'
      include 'msbarmasses.f'
      include 'hdecaymode.f'
      include 'first.f'

      integer:: j,k,nu
      real(dp):: msq(-nf:nf,-nf:nf),pin(mxpart,4),p(mxpart,4)
      real(dp):: pw1(4),pw2(4),q(4),a(4),r(4),b(4),h(4),
     & q1(4),q2(4),a1(4),a2(4)
      real(dp):: sh,sw1,sw2,qDq,aDa,rDr,bDb,densq,
     & p3Dp5,p6Dp8,q1Dq1,q2Dq2,a1Da1,a2Da2
      real(dp):: wtqqb,wtgg,hdecay,mt_eff,massfrun
      real(dp):: gamr1,gamr2,gamb1,gamb2,gamq4,gama7,dot
      real(dp):: fac,facqq,facgg,denr,denb,denq1,denq2,dena1,dena2
      real(dp):: p1Dr,p2Dr,p1Db,p2Db,p4Dq,p7Da,msqgamgam
      integer:: q4,a7,r1,r2,b1,b2
      parameter(q4=3,a7=5,r1=6,r2=8,b1=9,b2=10)
      save mt_eff
!$omp threadprivate(mt_eff)
      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif


      do nu=1,4
      do j=1,mxpart
      p(j,nu)=pin(j,nu)
      enddo

      h(nu)=p(9,nu)+p(10,nu)+p(11,nu)+p(12,nu)
      pw1(nu)=p(3,nu)+p(4,nu)
      pw2(nu)=p(7,nu)+p(8,nu)

      q(nu)=+p(3,nu)+p(4,nu)+p(5,nu)
      r(nu)=q(nu)+h(nu)
      q1(nu)=q(nu)+p(1,nu)
      q2(nu)=q(nu)+p(2,nu)

      a(nu)=-p(6,nu)-p(7,nu)-p(8,nu)
      b(nu)=a(nu)-h(nu)
      a1(nu)=a(nu)-p(1,nu)
      a2(nu)=a(nu)-p(2,nu)
      enddo


      sh=(h(4)**2-h(1)**2-h(2)**2-h(3)**2)
      sw1=(pw1(4)**2-pw1(1)**2-pw1(2)**2-pw1(3)**2)
      sw2=(pw2(4)**2-pw2(1)**2-pw2(2)**2-pw2(3)**2)
      qDq=(q(4)**2-q(1)**2-q(2)**2-q(3)**2)
      q1Dq1=(q1(4)**2-q1(1)**2-q1(2)**2-q1(3)**2)
      q2Dq2=(q2(4)**2-q2(1)**2-q2(2)**2-q2(3)**2)
      aDa=(a(4)**2-a(1)**2-a(2)**2-a(3)**2)
      a1Da1=(a1(4)**2-a1(1)**2-a1(2)**2-a1(3)**2)
      a2Da2=(a2(4)**2-a2(1)**2-a2(2)**2-a2(3)**2)
      rDr=(r(4)**2-r(1)**2-r(2)**2-r(3)**2)
      bDb=(b(4)**2-b(1)**2-b(2)**2-b(3)**2)

      p3Dp5=dot(p,3,5)
      p6Dp8=dot(p,6,8)

      densq=      ((sw1-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sw2-wmass**2)**2+wmass**2*wwidth**2)
      densq=densq*((sh-hmass**2)**2+hmass**2*hwidth**2)
      densq=densq*((qDq-mt**2)**2+mt**2*twidth**2)
      densq=densq*((aDa-mt**2)**2+mt**2*twidth**2)

      denr=rDr-mt**2
      denb=bDb-mt**2
      denq1=q1Dq1-mt**2
      denq2=q2Dq2-mt**2
      dena1=a1Da1-mt**2
      dena2=a2Da2-mt**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,9,10,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,9,10,hdecay)
      elseif (hdecaymode == 'wpwm') then
          call hwwdecay(p,9,10,11,12,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in qqb_tth'
      stop
      endif

      fac=2d0*p3Dp5*2d0*p6Dp8/densq
     &   *gwsq**5*gsq**2*mt_eff**2/wmass**2*hdecay
C---include factor for hadronic decays of top
      if (plabel(3) == 'pp') fac=2d0*xn*fac
      if (plabel(7) == 'pp') fac=2d0*xn*fac
      facqq=aveqq*V/4d0*fac
      facgg=avegg*V*xn/4d0*fac

      p4Dq=p(4,4)*q(4)-p(4,1)*q(1)-p(4,2)*q(2)-p(4,3)*q(3)
      p7Da=p(7,4)*a(4)-p(7,1)*a(1)-p(7,2)*a(2)-p(7,3)*a(3)

      p1Dr=p(1,4)*r(4)-p(1,1)*r(1)-p(1,2)*r(2)-p(1,3)*r(3)
      p2Dr=p(2,4)*r(4)-p(2,1)*r(1)-p(2,2)*r(2)-p(2,3)*r(3)

      p1Db=p(1,4)*b(4)-p(1,1)*b(1)-p(1,2)*b(2)-p(1,3)*b(3)
      p2Db=p(2,4)*b(4)-p(2,1)*b(1)-p(2,2)*b(2)-p(2,3)*b(3)

      gamq4=qDq/(2d0*p4Dq)
      gama7=aDa/(2d0*p7Da)

      gamb1=bDb/(2d0*p1Db)
      gamb2=bDb/(2d0*p2Db)
      gamr1=rDr/(2d0*p1Dr)
      gamr2=rDr/(2d0*p2Dr)

C     now the momenta 3,5,6,8,9,10 are no longer needed
C     so set
      do nu=1,4
      p(q4,nu)=q(nu)-gamq4*p(4,nu)
      p(a7,nu)=a(nu)-gama7*p(7,nu)
      p(r1,nu)=r(nu)-gamr1*p(1,nu)
      p(r2,nu)=r(nu)-gamr2*p(2,nu)
      p(b1,nu)=b(nu)-gamb1*p(1,nu)
      p(b2,nu)=b(nu)-gamb2*p(2,nu)
      enddo

      call spinoru(10,p,za,zb)
      call qqbtth(denr,denb,wtqqb)
      call ggtth(denr,denb,denq1,denq2,dena1,dena2,wtgg)

C----set all elements to zero
      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=facqq*wtqqb
      elseif (j == 0) then
          msq(j,j)=facgg*wtgg
      elseif (j > 0) then
          msq(j,-j)=facqq*wtqqb
      endif
      enddo

      return
      end
