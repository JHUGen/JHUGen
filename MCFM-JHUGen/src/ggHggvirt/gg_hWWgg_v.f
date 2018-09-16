      subroutine gg_hWWgg_v(p,msq)
      implicit none
      include 'types.f'
c--- Virtual matrix element squared averaged over initial colors and spins
c
c     g(-p1)+g(-p2)-->H -->  W^- (e^-(p5)+nubar(p6))
c                          + W^+ (nu(p3)+e^+(p4))+g(p_iglue1=7)+g(p_iglue2=8)
c
c    Calculation is fully analytic

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'deltar.f'
      integer:: j,k,i5,i6
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),s3456
      real(dp):: hdecay,Asq,fac
      real(dp):: qrqr,qarb,aqbr,abab,qbra,bqar
      real(dp):: qaqa,aqaq,qqqq,aaaa
      real(dp):: qagg,aqgg,qgqg,gqqg,agag,gaag,ggqa
      real(dp):: gggg
      real(dp):: Hqarbvsqanal
      real(dp):: Hqaqavsqanal
      real(dp):: HAQggvsqanal
      real(dp):: Hggggvsqanal
      logical:: CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)
      parameter(i5=7,i6=8)
C***************************************************
      scheme='dred'
C***************************************************

      if     (scheme == 'dred') then
        deltar=0._dp
      elseif (scheme == 'tH-V') then
        deltar=1._dp
      else
        write(6,*) 'Invalid scheme in gg_hgg_v.f'
        stop
      endif

c--- Set this to true to check squared matrix elements against
c--- hep-ph/0506196 using the point specified in Eq. (51)
      CheckEGZ=.false.

c--- Set up spinor products
      call spinoru(i6,p,za,zb)

      Asq=(as/(3._dp*pi))**2/vevsq

C   Deal with Higgs decay to WW
      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*wmass**2*s(3,5)*s(6,4)
      hdecay=hdecay/(((s3456-hmass**2)**2+(hmass*hwidth)**2)
     &   *((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
     &   *((s(5,6)-wmass**2)**2+(wmass*wwidth)**2))
      Asq=(as/(3._dp*pi))**2/vevsq
      fac=ason2pi*Asq*gsq**2*hdecay

c--- for checking EGZ
      if (CheckEGZ) then
        call CheckEGZres
      endif

c--- for checking scheme dependence of amplitudes
c      call CheckScheme(1,2,i5,i6)

C--- Note that Hqarbvsqanal(1,2,i5,i6)=Hqarbvsqanal(i6,i5,2,1)
C--- and the basic process is q(-ki6)+r(-k2)-->q(-k5)+r(-k1)

c--- FOUR-QUARK PROCESSES WITH NON-IDENTICAL QUARKS

C---quark-quark
C     q(1)+r(2)->q(i5)+r(i6)
      qrqr=Hqarbvsqanal(i6,2,i5,1)

C----quark-antiquark annihilation (i6-->i5-->2-->i6) wrt q(1)+r(2)->q(i5)+r(i6)
c     q(1)+a(2)->r(i5)+b(i6)
      qarb=Hqarbvsqanal(i5,i6,2,1)

C----antiquark-quark annihilation (1<-->2, i5<-->6) wrt to the above
c     a(1)+q(2)->b(i5)+r(i6)
c      aqbr=Hqarbvsqanal(i6,i5,1,2)
      aqbr=qarb

C----quark-antiquark scattering (i6<-->2) wrt q(1)+r(2)->q(i5)+r(i6)
c     q(1)+b(2)->r(i5)+a(i6)
      qbra=Hqarbvsqanal(2,i6,i5,1)

C----antiquark-quark scattering
c     b(1)+q(2)->a(i5)+r(i6) (1<-->2, i5<-->i6) wrt to the above
c      bqar=Hqarbvsqanal(1,i5,i6,2)
      bqar=qbra

C---antiquark-antiquark scattering (1<-->i5,2<-->i6) wrt q(1)+r(2)->q(i5)+r(i6)
C     a(1)+b(2)->a(i5)+b(i6)
      abab=Hqarbvsqanal(2,i6,1,i5)

C--- FOUR-QUARK PROCESSES WITH IDENTICAL QUARKS

C     q(1)+q(2)->q(i5)+q(i6)
      qqqq=qrqr+Hqarbvsqanal(i5,2,i6,1)+Hqaqavsqanal(i6,2,i5,1)

C     a(1)+a(2)->a(i5)+a(i6) (1<-->i5,2<-->i6) wrt q(1)+q(2)->q(i5)+q(i6)
      aaaa=abab+Hqarbvsqanal(2,i5,1,i6)+Hqaqavsqanal(2,i6,1,i5)

C     q(1)+a(2)->q(i5)+a(i6) (2<-->i6) wrt q(1)+q(2)->q(i5)+q(i6)
      qaqa=qbra+qarb+Hqaqavsqanal(2,i6,i5,1)

C     a(1)+q(2)->a(i5)+q(i6) (1<-->2, i5<-->i6) wrt the above
C      aqqa=qbra+qarb+Hqaqavsqanal(1,i5,i6,2)
      aqaq=qaqa

c--- TWO-QUARK, TWO GLUON PROCESSES

C     a(1)+q(2)->g(3)+g(4)
      aqgg=+HAQggvsqanal(2,1,i5,i6)

C     q(1)+g(2)->q(i5)+g(i6)
      qgqg=+HAQggvsqanal(1,i5,2,i6)

C     g(1)+q(2)->q(i5)+g(i6)
      gqqg=+HAQggvsqanal(2,i5,1,i6)

C     a(1)+g(2)->a(i5)+g(i6)
c      agag=+HAQggvsqanal(i5,1,2,i6)
      agag=qgqg

C     g(1)+a(2)->a(i5)+g(i6)
c      gaag=+HAQggvsqanal(i5,2,1,i6)
      gaag=gqqg

C     g(1)+g(2)->q(i5)+a(i6)
      ggqa=+HAQggvsqanal(i6,i5,1,2)

C     q(1)+a(2)->g(i5)+g(i6)
c      qagg=+HAQggvsqanal(1,2,i5,i6)
      qagg=aqgg

c--- FOUR GLUON PROCESS
      gggg=+Hggggvsqanal(1,2,i5,i6)


C--- DEBUGGING OUTPUT
C      write(6,*) 'qrqr',qrqr
C      write(6,*) 'qarb',qarb
C      write(6,*) 'aqrb',aqrb
C      write(6,*) 'abab',abab
C      write(6,*) 'qbra',qbra
C      write(6,*) 'bqra',bqra

C      write(6,*) 'Identical'
C      write(6,*) 'qaqa',qaqa
C      write(6,*) 'aqqa',aqqa
C      write(6,*) 'qqqq',qqqq
C      write(6,*) 'aaaa',aaaa


      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp

      if ((j==0).and.(k==0)) then
C---gg - all poles cancelled
             msq(j,k)=fac*avegg*(half*gggg+real(nflav,dp)*ggqa)

      elseif ((j>0).and.(k>0)) then
C---qq - all poles cancelled
             if (j==k) then
             msq(j,k)=aveqq*fac*half*qqqq
             else
             msq(j,k)=aveqq*fac*qrqr
             endif

      elseif ((j<0).and.(k<0)) then
C---aa - all poles cancelled
             if (j==k) then
             msq(j,k)=aveqq*fac*half*aaaa
             else
             msq(j,k)=aveqq*fac*abab
             endif

      elseif ((j>0).and.(k<0)) then
C----qa scattering - all poles cancelled
         if (j==-k) then
         msq(j,k)=aveqq*fac*(real(nflav-1,dp)*qarb+qaqa+half*qagg)
             else
         msq(j,k)=aveqq*fac*qbra
         endif

      elseif ((j<0).and.(k>0)) then
C----aq scattering - all poles cancelled
         if (j==-k) then
         msq(j,k)=aveqq*fac*(real(nflav-1,dp)*aqbr+aqaq+half*aqgg)
             else
         msq(j,k)=aveqq*fac*bqar
         endif

      elseif ((j==0).and.(k>0)) then
C----gq scattering - all poles cancelled
         msq(j,k)=aveqg*fac*gqqg

      elseif ((j==0).and.(k<0)) then
C----ga scattering - all poles cancelled
         msq(j,k)=aveqg*fac*gaag

      elseif ((j>0).and.(k==0)) then
C----qg scattering - all poles cancelled
         msq(j,k)=aveqg*fac*qgqg

      elseif ((j<0).and.(k==0)) then
C----ag scattering - all poles cancelled
         msq(j,k)=aveqg*fac*agag
      endif

      enddo
      enddo

      return
      end




