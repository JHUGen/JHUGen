      subroutine gg_HH(p,msq)
      implicit none
      include 'types.f'
C-----Matrix element squared for double Higgs production
C-----g(-p1)+g(-p2) --> H(p3,p4)+H(p5,p6)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'hdecaymode.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & s12,s34,s56,alphaw,p1(4),p2(4),p3(4),p4(4)
      real(dp):: msqgamgam,hdecay,hdecay2,fac,facPSZ
      complex(dp):: gauge(1,2),amp(2)
      do j=1,4
      p1(j)=p(1,j)
      p2(j)=p(2,j)
      p3(j)=p(3,j)+p(4,j)
      p4(j)=p(5,j)+p(6,j)
      enddo
C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          call hgamgamdecay(p,3,4,hdecay)
      else
      write(6,*) 'Unimplemented process in gg_HH'
      stop
      endif
      if (hdecaymode2 == 'tlta') then
          call htautaudecay(p,5,6,hdecay2)
      elseif (hdecaymode2 == 'bqba') then
          call hbbdecay(p,5,6,hdecay2)
      elseif (hdecaymode2 == 'gaga') then
          call hgamgamdecay(p,5,6,hdecay2)
      else
      write(6,*) 'Unimplemented process in gg_HH'
      stop
      endif
      s12=2._dp*(p(1,4)*p(2,4)-p(1,1)*p(2,1)-p(1,2)*p(2,2)-p(1,3)*p(2,3))
      s34=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      s56=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)
      hdecay2=hdecay2/((s56-hmass**2)**2+(hmass*hwidth)**2)
      alphaw=gwsq/(4._dp*pi)
      as=gsq/(4._dp*pi)
C     1024=2^10
C     dsig/dt=1/(2s)*(1/(8*pi*s)) |M|^2  (1/2)
C             (Flux)*(Phase space)       (Statistical Factor)
C     Thus matrix element without Statistical factor has factor:-
      fac=as**2*alphaw**2/(1024._dp*wmass**4)*hdecay*hdecay2
      call HHamps(p1,p2,p3,p4,gauge)
      msq(:,:)=0._dp
      msq(0,0)=fac
     & *(real(gauge(1,1)*conjg(gauge(1,1)))
     &  +real(gauge(1,2)*conjg(gauge(1,2))))

c      write(6,*) 'msq:GlB*2',2._dp*msq(0,0)
c      msq(0,0)=facPSZ
c     & *(real(amp(1)*conjg(amp(1)))
c     &  +real(amp(2)*conjg(amp(2))))
c      call PSZHHamps(p1,p2,p3,p4,amp)
c      GF=gwsq/8._dp/wmass**2*sqrt(2._dp)
c      facPSZ=GF**2*as**2*s12**2/(128._dp*pi**2)*hdecay*hdecay2
c      write(6,*) 'msq:PSZ',msq(0,0)
c      pause

      return
      end
