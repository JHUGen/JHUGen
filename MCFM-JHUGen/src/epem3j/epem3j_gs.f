      subroutine epem3j_gs(p,msq)
      implicit none
      include 'types.f'
************************************************************************
*     Author: J.M. Campbell                                            *
*     November, 2008.                                                  *
************************************************************************
c--- simple modification of qqb_w1jet_gs.f: permuted 1 and 4, 2 and 3
c--- to switch leptons with quarks and added a factor of Nc

c---Matrix element SUBTRACTION squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  W + parton(p5) + parton(p6)
c                           |
c                            -->l(p3)+a(p4)
c   positively charged W only


      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'nflav.f'
      integer:: j,k,nd
c --- remember: nd will count the dipoles

      real(dp):: p(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
      real(dp)::
     & msq45_3(-nf:nf,-nf:nf),msq35_4(-nf:nf,-nf:nf),
     & msq46_3(-nf:nf,-nf:nf),msq36_4(-nf:nf,-nf:nf),
     & msq45_6(-nf:nf,-nf:nf),msq46_5(-nf:nf,-nf:nf),
     & msq35_6(-nf:nf,-nf:nf),msq36_5(-nf:nf,-nf:nf),
     & msq56_3(-nf:nf,-nf:nf),msq56_4(-nf:nf,-nf:nf),
     & msq56_4v(-nf:nf,-nf:nf),msq56_3v(-nf:nf,-nf:nf),
     & msq36_5v(-nf:nf,-nf:nf),msq36_4v(-nf:nf,-nf:nf),
     & msq45_6v(-nf:nf,-nf:nf),msq46_2v(-nf:nf,-nf:nf),
     & msq45_3v(-nf:nf,-nf:nf),msq46_3v(-nf:nf,-nf:nf),
     & msq35_4v(-nf:nf,-nf:nf),
     & dummy(-nf:nf,-nf:nf),
     & sub45_3(4),sub35_4(4),sub46_3(4),sub36_4(4),
     & sub45_6(4),sub46_5(4),sub35_6(4),sub36_5(4),
     & sub56_4(4),sub56_3(4),sub56_4v,sub56_3v,
     & sub36_5v,sub36_4v,sub46_5v,sub46_3v,sub45_3v,sub45_6v,sub35_6v,
     & sub35_4v
      external epem3j,epem3j_gvec,donothing_gvec

      ndmax=10

c--- calculate all the dipoles
      call dips( 1,p,4,5,3,sub45_3,sub45_3v,msq45_3,dummy,
     & epem3j,donothing_gvec)
      call dips( 2,p,3,5,4,sub35_4,sub35_4v,msq35_4,dummy,
     & epem3j,donothing_gvec)
      call dips( 3,p,4,6,3,sub46_3,sub46_3v,msq46_3,dummy,
     & epem3j,donothing_gvec)
      call dips( 4,p,3,6,4,sub36_4,sub36_4v,msq36_4,dummy,
     & epem3j,donothing_gvec)

      call dips( 7,p,4,5,6,sub45_6,sub45_6v,msq45_6,dummy,
     & epem3j,donothing_gvec)
      call dips( 8,p,3,5,6,sub35_6,sub35_6v,msq35_6,dummy,
     & epem3j,donothing_gvec)
      call dips( 9,p,4,6,5,sub46_5,sub46_5v,msq46_5,dummy,
     & epem3j,donothing_gvec)
      call dips(10,p,3,6,5,sub36_5,sub36_5v,msq36_5,dummy,
     & epem3j,donothing_gvec)

      call dips(5,p,5,6,4,sub56_4,sub56_4v,msq56_4,msq56_4v,
     & epem3j,epem3j_gvec)
      call dips(6,p,5,6,3,sub56_3,sub56_3v,msq56_3,msq56_3v,
     & epem3j,epem3j_gvec)

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msq(nd,j,k)=zip
      enddo
      enddo
      enddo

      msq(1,0,1)=-msq45_3(0,0)*sub45_3(qq)/xn
      msq(2,0,1)=-msq35_4(0,0)*sub35_4(qq)/xn
      msq(3,0,1)=-msq46_3(0,0)*sub46_3(qq)/xn
      msq(4,0,1)=-msq36_4(0,0)*sub36_4(qq)/xn
      msq(5,0,0)=xn*(
     &   msq56_4(0,0)*sub56_4(gg)+msq56_4v(0,0)*sub56_4v)
      msq(5,1,0)=real(nflav,dp)*(
     &   msq56_4(0,0)*sub56_4(gq)-msq56_4v(0,0)*sub56_4v)
      msq(6,0,0)=xn*(
     &   msq56_3(0,0)*sub56_3(gg)+msq56_3v(0,0)*sub56_3v)
      msq(6,1,0)=real(nflav,dp)*(
     &   msq56_3(0,0)*sub56_3(gq)-msq56_3v(0,0)*sub56_3v)

      msq( 7,0,0)=xn*msq45_6(0,0)*sub45_6(qq)
      msq(10,0,0)=xn*msq36_5(0,0)*sub36_5(qq)

      msq( 9,0,0)=xn*msq46_5(0,0)*sub46_5(qq)
      msq( 8,0,0)=xn*msq35_6(0,0)*sub35_6(qq)

c--- note statistical factor of one half for two gluons in the final state
      do nd=1,ndmax
        msq(nd,0,0)=half*msq(nd,0,0)
        msq(nd,0,1)=half*msq(nd,0,1)
        msq(nd,1,0)=half*msq(nd,1,0)
      enddo

      return
      end

