      subroutine gg_hgg_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'msq_struc.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'kprocess.f'
      include 'hdecaymode.f'
C  in is the label of the momentum contracted with n
      integer::j,k,in
      real(dp)::msq(-nf:nf,-nf:nf),msqgamgam
      real(dp)::n(4),p(mxpart,4),hdecay,s34,fac,
     & qqgghn_ab,qqgghn_ba,qqgghn_sym,
     & c1234,c1243,c1423,Asq

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=zip
      enddo
      enddo

C---fill dot products
      call spinoru(6,p,za,zb)

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
      write(6,*) 'Unimplemented process in gg_hgg_gvec'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(three*pi))**2/vevsq
      fac=gsq**2*Asq*hdecay

c      do j=-1,+1
c      do k=-1,+1
c      p1p2(j,k)=zip
c      enddo
c      enddo

c--- NOTE: there seems to be some redundancy in the function calls, e.g.
c--- the entries for (0,-1) = (0,+1). Need to check if this is true
c--- and eliminate where possible

C     The function qqgghn calculates q(p1)+qbar(p2)-->H+g(p3)+g(p4)
C     with p4 contracted with n
C     The function gggghn calculates g(p1)+g(p2)-->H+g(p3)+g(p4)
C     with p1 contracted with n

c--- Note that I have removed all references to p1p2(j,k) since
c--- the appropriate terms are actually used only in their
c--- colour-separated forms

      if     (in == 1) then
        call qqgghn(2,5,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,-1)=msq_strucv(igg_ab,0,-1)+msq_strucv(igg_ba,0,-1)
c     &            +msq_strucv(igg_sym,0,-1)
c        p1p2(0,+1)=-aveqg*fac*qqgghn_old(5,2,6,1,p,n)
        call qqgghn(5,2,6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,-1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,0,-1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,0,-1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     &            +msq_strucv(igg_sym,0,+1)
        call qqgghn(5,6,2,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,0)=avegg*fac*real(nf,dp)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*real(nf,dp)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*real(nf,dp)*qqgghn_sym
c        p1p2(0,0)=half*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     &                  +msq_strucv(igggg_c,0,0))/two
c     &      +real(nf,dp)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     &                  +msq_strucv(igg_sym,0,0))
        call gggghn_amp(1,2,5,6,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1243
c        p1p2(0,0)=avegg*fac*(half*(c1234+c1243+c1423)
c     &                 +real(nf,dp)*qqgghn_old(5,6,2,1,p,n))

c        p1p2(0,-1)=-aveqg*fac*qqgghn_old(2,5,6,1,p,n)

      elseif (in == 2) then
c        p1p2(+1,0)=-aveqg*fac*qqgghn_old(1,5,6,2,p,n)
        call qqgghn(1,5,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     &            +msq_strucv(igg_sym,+1,0)
c        p1p2(-1,0)=-aveqg*fac*qqgghn_old(5,1,6,2,p,n)
        call qqgghn(5,1,6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,-1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,-1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(-1,0)=msq_strucv(igg_ab,-1,0)+msq_strucv(igg_ba,-1,0)
c     &            +msq_strucv(igg_sym,-1,0)
        call qqgghn(6,5,1,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,0)=avegg*fac*real(nf,dp)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*real(nf,dp)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*real(nf,dp)*qqgghn_sym
c        p1p2(0,0)=half*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     &                  +msq_strucv(igggg_c,0,0))/two
c     &      +real(nf,dp)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     &                  +msq_strucv(igg_sym,0,0))
        call gggghn_amp(2,1,5,6,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1243
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1234
c        p1p2(0,0)=+avegg*fac*(half*(c1234+c1243+c1324)
c     &                 +real(nf,dp)*qqgghn_old(5,6,1,2,p,n))

      elseif (in == 5) then
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,6,5,p,n)
        call qqgghn(1,2,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*half*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     &             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,6,5,p,n)
        call qqgghn(2,1,6,5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*half*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     &             +msq_strucv(igg_sym,-1,+1)
       call gggghn_amp(5,6,1,2,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1243
c        p1p2(0,0)=+half*avegg*fac*(c1234+c1243+c1423)

      elseif (in == 6) then
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,5,6,p,n)
        call qqgghn(1,2,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*half*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     &             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,5,6,p,n)
        call qqgghn(2,1,5,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*half*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     &             +msq_strucv(igg_sym,-1,+1)
c--- for the qg, gq pieces, note that qbar-g and g-qbar are never used
        call qqgghn(1,5,2,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     &            +msq_strucv(igg_sym,+1,0)
        call qqgghn(2,5,1,6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     &            +msq_strucv(igg_sym,0,+1)
        call gggghn_amp(6,1,2,5,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1243
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1423
c        p1p2(0,0)=+half*avegg*fac*(c1234+c1243+c1423)
      endif

      return
      end


