      subroutine middleh(p,first,middlebit)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'metric0.f'
      include 'alpha1.f'
      include 'masses.f'
      include 'TRydef.f'
      include 'currentdecl.f'
      include 'tensordecl.f'
      integer:: fi,nu,ro,ep
      complex(dp):: prW,string,middlebit(-2:0),part3,part1,bit
      real(dp):: p(mxpart,4),vec,mtsq
      real(dp):: p1(4),p16(4),p2(4),p5(4),p6(4),p25(4),p34(4),
     & s16,s25,epbit,omal
      integer:: epmin
      logical:: first
c
c     statement function
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
c     end statement function

      omal=1d0-alpha1

      mtsq=mt**2
      p1(:)= p(1,:)
      p2(:)= p(2,:)
      p5(:)= p(5,:)
      p6(:)= p(6,:)
      p16(:)=p(1,:)+p(6,:)
      p34(:)=p(3,:)+p(4,:)
      p25(:)=p(2,:)+p(5,:)
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s25=p25(4)**2-p25(1)**2-p25(2)**2-p25(3)**2

      call doBtensor(p5,0d0,mtsq,FB0x5,FB1x5,
     &  FB2x5,FB3x5,FB4x5,FB5x5,FB6x5)

      call doCtensor(p1,p16,0d0,0d0,0d0,FC0x1x6,FC1x1x6,
     &  FC2x1x6,FC3x1x6,FC4x1x6,FC5x1x6,FC6x1x6)
      call doCtensor(p2,p25,0d0,0d0,mtsq,FC0x2x5,FC1x2x5,
     &  FC2x2x5,FC3x2x5,FC4x2x5,FC5x2x5,FC6x2x5)

c--- only compute poles for checking on first call
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif
c      write(*,*) 'epmin in middleh', epmin

      middlebit(:)=czip

      do ep=epmin,0
      epbit=0d0
      if (ep == 0) epbit=1d0

      part3=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      string=g0(fi)*g0(nu)*g0(ro)*j61x3(fi,nu,ro)
      bit= - prW(s16)*prW(s25)*J52x0*vec(p1,ro)*vec(p34,nu)*FC1x1x6(fi,
     & ep)*omal*mt*wmass**(-1) - prW(s16)*prW(s25)*J52x0*vec(p6,ro)*
     & vec(p34,nu)*FC1x1x6(fi,ep)*omal*mt*wmass**(-1) - prW(s16)*prW(
     & s25)*J52x0*vec(p34,nu)*FC2x1x6(y2(fi,ro),ep)*omal*mt*wmass**(-1)
     &  - prW(s16)*prW(s25)*J52x1(nu)*vec(p1,ro)*FC1x1x6(fi,ep)*wmass
     &  - prW(s16)*prW(s25)*J52x1(nu)*vec(p6,ro)*FC1x1x6(fi,ep)*wmass
     &  - prW(s16)*prW(s25)*J52x1(nu)*FC2x1x6(y2(fi,ro),ep)*wmass

      part3=part3+bit*string
      enddo
      enddo
      enddo


      part1=czip
      do fi=1,4
      string=g0(fi)*j61x1(fi)
      bit= - 1.D0/2.D0*prW(s16)*prW(s25)*J52x0*vec(p34,fi)*omal*mt*
     & epbit*wmass**(-1) + 2.D0*prW(s16)*prW(s25)*J52x0*FC1x2x5(fi,ep)*
     & mt*wmass - prW(s16)*prW(s25)*J52x1(fi)*epbit*wmass

      part1=part1+bit*string
      do nu=1,4
      bit=2.D0*prW(s16)*prW(s25)*J52x0*vec(p25,nu)*vec(p34,fi)*FC1x2x5(
     & nu,ep)*omal*mt*wmass**(-1) + prW(s16)*prW(s25)*J52x1(nu)*vec(p5,
     & nu)*vec(p34,fi)*FB0x5(ep)*omal*wmass**(-1) - prW(s16)*prW(s25)*
     & J52x1(nu)*vec(p34,fi)*FC1x2x5(nu,ep)*omal*mtsq*wmass**(-1) +
     & prW(s16)*prW(s25)*J52x1(nu)*vec(p34,fi)*FB1x5(nu,ep)*omal*
     & wmass**(-1)

      part1=part1+g0(nu)*bit*string
      do ro=1,4
      bit= - prW(s16)*prW(s25)*vec(p25,nu)*FC1x2x5(ro,ep)*J52x3(ro,fi,
     & nu)*wmass - prW(s16)*prW(s25)*FC2x2x5(y2(nu,ro),ep)*J52x3(nu,fi,
     & ro)*wmass

      part1=part1+g0(nu)*g0(ro)*bit*string
      enddo
      enddo
      enddo

      middlebit(ep)=part1+part3

      enddo


      return
      end
