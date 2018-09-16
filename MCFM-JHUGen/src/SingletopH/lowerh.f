      subroutine lowerh(p,first,lowerbit)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'metric0.f'
      include 'masses.f'
      include 'scale.f'
      include 'TRydef.f' 
      include 'currentdecl.f'
      include 'tensordecl.f'
      integer:: fi,nu,ro,si,ep
      complex(dp):: prW,prt,string,lowerbit(-2:0),bit,
     & part4,part3,part2,part1,part0
      real(dp):: p(mxpart,4),vec,mtsq,epbit,Zm(-2:0),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     & p25(4),p34(4),p16(4),p345(4),p2345(4),s16,s34,s345
      integer:: epmin
      logical:: first
c     
c     statement functions
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
      prt(s16)=cone/cplx2(s16-mt**2,zip)
c     end statement functions

      mtsq=mt**2
      p1(:)= p(1,:)
      p2(:)= p(2,:)
      p3(:)= p(3,:)
      p4(:)= p(4,:)
      p5(:)= p(5,:)
      p6(:)= p(6,:)
      p16(:)=p1(:)+p6(:)
      p25(:)=p2(:)+p5(:)
      p34(:)=p3(:)+p4(:)
      p345(:)=p3(:)+p4(:)+p5(:)
      p2345(:)=p2(:)+p3(:)+p4(:)+p5(:)
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s345=p345(4)**2-p345(1)**2-p345(2)**2-p345(3)**2
       
      call doBtensor(p345,0d0,mtsq,FB0x345,FB1x345,
     &  FB2x345,FB3x345,FB4x345,FB5x345,FB6x345)

      call doCtensor(p1,p16,0d0,0d0,0d0,FC0x1x6,FC1x1x6,
     &  FC2x1x6,FC3x1x6,FC4x1x6,FC5x1x6,FC6x1x6)
      call doCtensor(p5,p345,0d0,mtsq,mtsq,FC0x5x34,FC1x5x34,
     &  FC2x5x34,FC3x5x34,FC4x5x34,FC5x5x34,FC6x5x34)
      call doCtensor(p2,p2345,0d0,0d0,mtsq,FC0x2x345,FC1x2x345,
     &  FC2x2x345,FC3x2x345,FC4x2x345,FC5x2x345,FC6x2x345)
       
      call doDtensor(p2,p25,p2345,0d0,0d0,mtsq,mtsq,FD0x2x5x34,
     &FD1x2x5x34,FD2x2x5x34,FD3x2x5x34,FD4x2x5x34,FD5x2x5x34,FD6x2x5x34)
       
      Zm(-2)=0d0
      Zm(-1)=3d0
      Zm( 0)=3d0*log(musq/mtsq)+5d0

c--- only compute poles for checking on first call
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif
      
      lowerbit(:)=czip
     
      do ep=epmin,0
     
      epbit=0d0
      if (ep == 0) epbit=1d0

      part4=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      string=g0(fi)*g0(nu)*g0(ro)*g0(si)*J52x4(fi,nu,ro,si)
      bit= + prt(s345)**2*J61x1(si)*mt*wmass**(-1) * (  - 1.D0/2.D0*
     &    vec(p345,ro)*vec(p345,fi)*FB1x345(nu,ep) )
      bit = bit + prt(s345)*J61x1(ro)*mt*wmass**(-1) * (  - 1.D0/2.D0*
     &    vec(p345,si)*vec(p345,fi)*FC1x2x345(nu,ep) - 1.D0/2.D0*vec(
     &    p345,fi)*FC2x2x345(y2(nu,si),ep) )
      bit = bit + J61x1(si)*mt*wmass**(-1) * ( 1.D0/2.D0*vec(p25,nu)*
     &    vec(p34,ro)*FD1x2x5x34(fi,ep) + 1.D0/2.D0*vec(p34,ro)*
     &    FD2x2x5x34(y2(fi,nu),ep) )
      bit = bit + J61x1(fi)*mt*wmass**(-1) * ( 1.D0/2.D0*vec(p25,ro)*
     &    vec(p34,nu)*FD1x2x5x34(si,ep) + 1.D0/2.D0*vec(p34,nu)*
     &    FD2x2x5x34(y2(ro,si),ep) )

      part4=part4+bit*string
      enddo
      enddo
      enddo
      enddo


      part3=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      string=g0(fi)*g0(nu)*g0(ro)*J52x3(fi,nu,ro)
      bit= + prt(s345)**2*J61x1(ro)*mt**2*wmass**(-1) * (  - 1.D0/2.D0*
     &    vec(p345,nu)*FB1x345(fi,ep) - 1.D0/2.D0*vec(p345,fi)*FB1x345(
     &    nu,ep) )
      bit = bit + prt(s345)*J61x1(nu)*mt**2*wmass**(-1) * (  - 1.D0/2.D0
     &    *vec(p345,ro)*FC1x2x345(fi,ep) - 1.D0/2.D0*FC2x2x345(y2(fi,ro
     &    ),ep) )
      bit = bit + prt(s345)*J61x1(ro)*mt**2*wmass**(-1) * (  - 1.D0/2.D0
     &    *vec(p5,fi)*vec(p345,nu)*FC0x5x34(ep) - vec(p345,nu)*
     &    FC1x5x34(fi,ep) )
      bit = bit + J61x1(nu)*mt**2*wmass**(-1) * (  - 1.D0/2.D0*vec(p5,
     &    ro)*FD1x2x5x34(fi,ep) - 1.D0/2.D0*vec(p345,ro)*FD1x2x5x34(fi,
     &    ep) - FD2x2x5x34(y2(fi,ro),ep) )

      part3=part3+bit*string
      enddo
      enddo
      enddo


      part2=czip
      do fi=1,4
      do nu=1,4
      string=g0(fi)*g0(nu)*J52x2(fi,nu)
      bit= + prt(s345)**2*J61x1(nu)*mt**3*wmass**(-1) * (  - 1.D0/2.D0*
     &    vec(p345,fi)*Zm(ep) + vec(p345,fi)*FB0x345(ep) - 1.D0/2.D0*
     &    FB1x345(fi,ep) )
      bit = bit + prt(s345)*J61x1(nu)*mt*wmass**(-1) * (  - 1.D0/2.D0*
     &    vec(p345,fi)*epbit - 1.D0/4.D0*vec(p345,fi)*Zm(ep) + 1.D0/2.D0
     &    *vec(p345,fi)*FC0x5x34(ep)*s345 - 1.D0/2.D0*vec(p345,fi)*
     &    FC0x5x34(ep)*s34 + 1.D0/2.D0*vec(p345,fi)*FB0x345(ep) )
      bit = bit + prt(s345)*J61x1(nu)*mt**3*wmass**(-1) * (  - 1.D0/2.D0
     &    *vec(p5,fi)*FC0x5x34(ep) + vec(p345,fi)*FC0x5x34(ep) - 
     &    FC1x5x34(fi,ep) )

      part2=part2+bit*string
      do ro=1,4
      bit= + prt(s345)*J61x1(nu)*mt*wmass**(-1) * ( vec(p34,ro)*vec(
     &    p345,fi)*FC1x5x34(ro,ep) )

      part2=part2+g0(ro)*bit*string
      do si=1,4
      bit= + prt(s345)*mt*wmass**(-1) * (  - 1.D0/2.D0*vec(p6,ro)*vec(
     &    p345,fi)*FC1x1x6(si,ep)*J61x3(si,nu,ro) - 1.D0/2.D0*vec(p345,
     &    fi)*FC2x1x6(y2(ro,si),ep)*J61x3(ro,nu,si) )

      part2=part2+g0(ro)*g0(si)*bit*string
      enddo
      enddo
      enddo
      enddo


      part1=czip
      do fi=1,4
      string=g0(fi)*J52x1(fi)
      bit= + prt(s345)**2*J61x1(fi)*mt**4*wmass**(-1) * (  - 1.D0/2.D0*
     &    Zm(ep) + FB0x345(ep) )
      bit = bit + prt(s345)*J61x1(fi)*mt**2*wmass**(-1) * (  - 1.D0/2.D0
     &    *epbit - 1.D0/2.D0*Zm(ep) + 1.D0/2.D0*FC0x5x34(ep)*s345 - 1.D0
     &    /2.D0*FC0x5x34(ep)*s34 + FB0x345(ep) )
      bit = bit + prt(s345)*J61x1(fi)*mt**4*wmass**(-1) * ( FC0x5x34(ep
     &    ) )
      bit = bit + J61x1(fi)*mt**2*wmass**(-1) * (  - 1.D0/2.D0*
     &    FC0x5x34(ep) )

      part1=part1+bit*string
      do nu=1,4
      bit= + prt(s345)*J61x1(nu)*mt**2*wmass**(-1) * ( vec(p345,fi)*
     &    FC1x2x345(nu,ep) )
      bit = bit + prt(s345)*J61x1(fi)*mt**2*wmass**(-1) * ( vec(p34,nu)
     &    *FC1x5x34(nu,ep) )

      part1=part1+g0(nu)*bit*string
      do ro=1,4
      bit= + prt(s345)*mt**2*wmass**(-1) * (  - 1.D0/2.D0*vec(p6,nu)*
     &    FC1x1x6(ro,ep)*J61x3(ro,fi,nu) - 1.D0/2.D0*FC2x1x6(y2(nu,ro),
     &    ep)*J61x3(nu,fi,ro) )

      part1=part1+g0(nu)*g0(ro)*bit*string
      enddo
      enddo
      enddo

      part0=czip
      string=J52x0
      do fi=1,4
      bit= + prt(s345)*J61x1(fi)*mt**3*wmass**(-1) * ( FC1x2x345(fi,ep)
     &     )
      bit = bit + J61x1(fi)*mt*wmass**(-1) * ( FC1x2x345(fi,ep) )
      bit = bit + J61x1(fi)*mt**3*wmass**(-1) * ( 2.D0*FD1x2x5x34(fi,ep
     &    ) )

      part0=part0+g0(fi)*bit*string
      enddo

      lowerbit(ep)=prW(s16)*(part0+part1+part2+part3+part4)

      enddo

      return
      end
