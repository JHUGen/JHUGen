      subroutine scalarh(p,first,scalarbit)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'alpha1.f'
      include 'metric0.f'
      include 'masses.f'
      include 'scale.f'
      include 'TRydef.f' 
      include 'currentdecl.f'
      include 'tensordecl.f'
      integer:: fi,nu,ro,ep
      complex(dp):: prW,string,scalarbit(-2:0),bit,
     & part3,part1
      real(dp):: p(mxpart,4),vec,mtsq,epbit,Zm(-2:0),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),s25,
     & p25(4),p34(4),p16(4),s16
      integer:: epmin
      logical:: first
c     
c     statement functions
      prW(s16)=cone/cplx2(s16-wmass**2,zip)
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
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s25=p25(4)**2-p25(1)**2-p25(2)**2-p25(3)**2
       
      call doBtensor(p5,0d0,mtsq,FB0x5,FB1x5,
     &  FB2x5,FB3x5,FB4x5,FB5x5,FB6x5)
       
      call doCtensor(p1,p16,0d0,0d0,0d0,FC0x1x6,FC1x1x6,
     &  FC2x1x6,FC3x1x6,FC4x1x6,FC5x1x6,FC6x1x6)
      call doCtensor(p2,p25,0d0,0d0,mtsq,FC0x2x5,FC1x2x5,
     &  FC2x2x5,FC3x2x5,FC4x2x5,FC5x2x5,FC6x2x5)
      Zm(-2)=0d0
      Zm(-1)=3d0
      Zm( 0)=3d0*log(musq/mtsq)+5d0

c--- only compute poles for checking on first call
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif
c      write(*,*) 'epmin in scalarh', epmin
     
      scalarbit(:)=czip
     
      do ep=epmin,0
      epbit=0d0
      if (ep == 0) epbit=1d0

      part3=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      string=g0(fi)*g0(nu)*g0(ro)*j61x3(fi,nu,ro)
      bit= - prW(s16)*prW(s25)*J52x0*vec(p1,ro)*vec(p34,nu)*FC1x1x6(fi,
     & ep)*alpha1*mt*wmass**(-1) - prW(s16)*prW(s25)*J52x0*vec(p6,ro)*
     & vec(p34,nu)*FC1x1x6(fi,ep)*alpha1*mt*wmass**(-1) - prW(s16)*prW(
     & s25)*J52x0*vec(p34,nu)*FC2x1x6(y2(fi,ro),ep)*alpha1*mt*
     & wmass**(-1)

      part3=part3+bit*string
      enddo
      enddo
      enddo


      part1=czip
      do fi=1,4
      string=g0(fi)*j61x1(fi)
      bit= - 1.D0/2.D0*prW(s16)*prW(s25)*J52x0*vec(p34,fi)*alpha1*mt*
     & epbit*wmass**(-1) - 1.D0/2.D0*prW(s16)*prW(s25)*J52x0*vec(p34,fi
     & )*Zm(ep)*alpha1*mt*wmass**(-1) + 2.D0*prW(s16)*prW(s25)*J52x0*
     & vec(p34,fi)*FB0x5(ep)*alpha1*mt*wmass**(-1)

      part1=part1+bit*string
      do nu=1,4
      bit= - prW(s16)*prW(s25)*J52x1(nu)*vec(p34,fi)*FC1x2x5(nu,ep)*
     & alpha1*mtsq*wmass**(-1)

      part1=part1+g0(nu)*bit*string
      do ro=1,4
      bit=prW(s16)*prW(s25)*vec(p25,nu)*vec(p34,fi)*FC1x2x5(ro,ep)*
     & J52x2(nu,ro)*alpha1*mt*wmass**(-1) + prW(s16)*prW(s25)*vec(p25,
     & nu)*vec(p34,fi)*FC1x2x5(ro,ep)*J52x2(ro,nu)*alpha1*mt*
     & wmass**(-1)

      part1=part1+g0(nu)*g0(ro)*bit*string
      enddo
      enddo
      enddo

      scalarbit(ep)=part1+part3

      enddo

      return
      end
