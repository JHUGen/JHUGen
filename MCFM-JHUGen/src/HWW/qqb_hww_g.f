      subroutine qqb_hww_g(p,msq)
      implicit none
c----NLO matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
C----averaged over initial colours and spins
c     g(-p1)+g(-p2)-->H -->  W^- (e^-(p5)+nubar(p6)) 
c                          + W^+ (nu(p3)+e^+(p4))+g(p7)
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'process.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),fac
      double precision sh,ss,tt,uu,decay,s(mxpart,mxpart)
      double precision aw,qqb,qg,gq,gg

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      aw=gwsq/(4d0*pi)
      call dotem(7,p,s)
      decay=gwsq**3*wmass**2*s(5,3)*s(6,4)
      if (case .eq. 'HWW2lq') decay=2d0*xn*decay

c--   calculate propagators
      ss=s(1,2)
      tt=s(1,7)
      uu=s(2,7)
      sh=s(1,2)+s(1,7)+s(2,7)

      gg=aw*as**3*4d0*V/9d0*xn*(sh**4+ss**4+tt**4+uu**4)
     . /(ss*tt*uu*wmass**2)
      qqb=aw*as**3*2d0*V/9d0*(tt**2+uu**2)/(ss*wmass**2)
      gq=-aw*as**3*2d0*V/9d0*(ss**2+tt**2)/(uu*wmass**2)
      qg=-aw*as**3*2d0*V/9d0*(ss**2+uu**2)/(tt*wmass**2)

      fac=one/((s(3,4)-wmass**2)**2+(wmass*wwidth)**2)
      fac=fac/((s(5,6)-wmass**2)**2+(wmass*wwidth)**2)
      fac=fac/((sh-hmass**2)**2+(hmass*hwidth)**2)


      gg=avegg*fac*gg*decay
      gq=aveqg*fac*gq*decay
      qg=aveqg*fac*qg*decay
      qqb=aveqq*fac*qqb*decay

c--set msq=0 to initialize
      do j=-nf,nf
      do k=-nf,nf
      if ((k .eq. -j) .and. (j .ne. 0)) then
      msq(j,k)=qqb
      elseif ((j .eq. 0) .and. (k .ne. 0)) then
      msq(j,k)=gq
      elseif ((j .ne. 0) .and. (k .eq. 0)) then
      msq(j,k)=qg
      elseif ((k .eq. 0) .and. (j .eq. 0)) then
      msq(j,k)=gg
      endif
      enddo
      enddo
      return
      end
