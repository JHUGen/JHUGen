      subroutine lowerdk_partbox(p,h3,lowerbox,first)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zcouple.f'
      include 'metric0.f'
      include 'masses.f'
      include 'scale.f'
      include 'TRydef.f' 
      include 'currentdecl.f'
      include 'tensordecl.f'
      integer:: fi,nu,ro,si,om,ep,h3
      complex(dp):: prW,prt,prq,string,lowerbox(-2:0),bit,iprZ,
     & part5,part4,part3,part2,part1,cprop
      complex(dp):: facuRl,facuLl,facdLl
      real(dp):: p(mxpart,4),mtsq,Zm(-2:0),
     & p1(4),p2(4),p3(4),p4(4),p5(4),p6(4),
     & p25(4),p34(4),p16(4),p345(4),p234(4),p2345(4),s16,s34,s234,s345
      integer:: epmin
      logical:: first
c     
c     statement functions
      prW(s34)=cone/cplx2(s34-wmass**2,zip)
      prt(s34)=cone/cplx2(s34-mt**2,zip)
      prq(s34)=cone/cplx1(s34)
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
      p234(:)=p2(:)+p3(:)+p4(:)
      p2345(:)=p2(:)+p3(:)+p4(:)+p5(:)
      s34=p34(4)**2-p34(1)**2-p34(2)**2-p34(3)**2
      s16=p16(4)**2-p16(1)**2-p16(2)**2-p16(3)**2
      s234=p234(4)**2-p234(1)**2-p234(2)**2-p234(3)**2
      s345=p345(4)**2-p345(1)**2-p345(2)**2-p345(3)**2
       
c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=cplx1(1d0/sqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      iprZ=cplx1(s34-zmass**2)

       cprop=cprop/cplx2(zip,mt*twidth)

      if (h3 == -1) then
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*le)
        facuRl=cplx1(Qu*q1)*iprZ/s34+cplx1(R(2)*le)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*le)
c        facdRl=cplx1(Qd*q1)*iprZ/s34+cplx1(R(1)*le)
      else
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*re)
        facuRl=cplx1(Qu*q1)*iprZ/s34+cplx1(R(2)*re)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*re)
c        facdRl=cplx1(Qd*q1)*iprZ/s34+cplx1(R(1)*re)
      endif
       
      call doBtensor(p234,0d0,0d0,FB0x234,FB1x234,
     &  FB2x234,FB3x234,FB4x234,FB5x234,FB6x234)
      call doBtensor(p345,0d0,mtsq,FB0x345,FB1x345,
     &  FB2x345,FB3x345,FB4x345,FB5x345,FB6x345)

      call doCtensor(p1,p16,0d0,0d0,0d0,FC0x1x6,FC1x1x6,
     &  FC2x1x6,FC3x1x6,FC4x1x6,FC5x1x6,FC6x1x6)
      call doCtensor(p5,p345,0d0,mtsq,mtsq,FC0x5x34,FC1x5x34,
     &  FC2x5x34,FC3x5x34,FC4x5x34,FC5x5x34,FC6x5x34)
      call doCtensor(p2,p234,0d0,0d0,0d0,FC0x2x34,FC1x2x34,
     &  FC2x2x34,FC3x2x34,FC4x2x34,FC5x2x34,FC6x2x34)
      call doCtensor(p2,p2345,0d0,0d0,mtsq,FC0x2x345,FC1x2x345,
     &  FC2x2x345,FC3x2x345,FC4x2x345,FC5x2x345,FC6x2x345)
      call doCtensor(p234,p2345,0d0,0d0,mtsq,FC0x234x5,FC1x234x5,
     &  FC2x234x5,FC3x234x5,FC4x234x5,FC5x234x5,FC6x234x5)
       
      call doDtensor(p34,p234,p2345,0d0,0d0,0d0,mtsq,FD0x34x2x5,
     &FD1x34x2x5,FD2x34x2x5,FD3x34x2x5,FD4x34x2x5,FD5x34x2x5,FD6x34x2x5)
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

c      write(*,*) 'epmin in lowerdk_box', epmin
      if (epmin == -2) then
c        write(*,*) ' calculating all poles' 
      endif
      do ep=epmin,0
     

      part5=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do om=1,4
      string=g0(fi)*g0(nu)*g0(ro)*g0(si)*g0(om)*J52x5(fi,nu,ro,si,om)
      bit= + prt(s345)**2*J34x1(fi)*J61x1(om)*facuLl * ( FB1x345(ro,ep)
     &    *p345(nu)*p345(si) )
      bit = bit + J34x1(nu)*J61x1(ro)*facuLl * (  - FD3x2x5x34(y3(fi,si
     &    ,om),ep) - FD2x2x5x34(y2(si,om),ep)*p25(fi) - FD2x2x5x34(y2(
     &    fi,si),ep)*p345(om) - FD1x2x5x34(si,ep)*p25(fi)*p345(om) )
      bit = bit + J34x1(nu)*J61x1(si)*facuLl * ( FD3x2x5x34(y3(fi,ro,om
     &    ),ep) + FD2x2x5x34(y2(ro,om),ep)*p25(fi) + FD2x2x5x34(y2(fi,
     &    om),ep)*p2345(ro) + FD1x2x5x34(om,ep)*p25(fi)*p2345(ro) )
      bit = bit + J34x1(nu)*J61x1(om)*facuLl * ( FD3x2x5x34(y3(fi,ro,si
     &    ),ep) + FD2x2x5x34(y2(ro,si),ep)*p25(fi) + FD2x2x5x34(y2(fi,
     &    si),ep)*p2345(ro) + FD1x2x5x34(si,ep)*p25(fi)*p2345(ro) )
      bit = bit + J34x1(ro)*J61x1(nu)*facdLl * (  - FD3x34x2x5(y3(fi,si
     &    ,om),ep) - FD2x34x2x5(y2(si,om),ep)*p2345(fi) - FD2x34x2x5(
     &    y2(fi,om),ep)*p34(si) - FD1x34x2x5(om,ep)*p34(si)*p2345(fi) )
      bit = bit + J34x1(ro)*J61x1(fi)*facdLl * (  - FD3x34x2x5(y3(nu,si
     &    ,om),ep) - FD2x34x2x5(y2(nu,si),ep)*p345(om) - FD2x34x2x5(y2(
     &    nu,om),ep)*p34(si) - FD1x34x2x5(nu,ep)*p34(si)*p345(om) )
      bit = bit + J34x1(ro)*J61x1(om)*facdLl * ( FD3x34x2x5(y3(fi,nu,si
     &    ),ep) + FD2x34x2x5(y2(nu,si),ep)*p2345(fi) + FD2x34x2x5(y2(fi
     &    ,nu),ep)*p34(si) + FD1x34x2x5(nu,ep)*p34(si)*p2345(fi) )
      bit = bit + J34x1(si)*J61x1(nu)*facdLl * ( FD3x34x2x5(y3(fi,ro,om
     &    ),ep) + FD2x34x2x5(y2(ro,om),ep)*p2345(fi) + FD2x34x2x5(y2(fi
     &    ,ro),ep)*p34(om) + FD1x34x2x5(ro,ep)*p34(om)*p2345(fi) )
      bit = bit + J34x1(fi)*J61x1(ro)*facuLl * (  - FD3x2x5x34(y3(nu,si
     &    ,om),ep) - FD2x2x5x34(y2(nu,si),ep)*p5(om) - FD2x2x5x34(y2(si
     &    ,om),ep)*p2345(nu) - FD1x2x5x34(si,ep)*p5(om)*p2345(nu) )
      bit = bit + J34x1(om)*J61x1(nu)*facdLl * ( FD3x34x2x5(y3(fi,ro,si
     &    ),ep) + FD2x34x2x5(y2(ro,si),ep)*p2345(fi) + FD2x34x2x5(y2(fi
     &    ,ro),ep)*p34(si) + FD1x34x2x5(ro,ep)*p34(si)*p2345(fi) )
      bit = bit + J34x1(om)*J61x1(ro)*facuLl * ( FD3x2x5x34(y3(fi,nu,si
     &    ),ep) + FD2x2x5x34(y2(nu,si),ep)*p25(fi) + FD2x2x5x34(y2(fi,
     &    si),ep)*p2345(nu) + FD1x2x5x34(si,ep)*p25(fi)*p2345(nu) )
      bit = bit + J34x1(om)*J61x1(fi)*facdLl * (  - prq(s234)**2*
     &    FB1x234(ro,ep)*p234(nu)*p234(si) )

      part5=part5+bit*string
      enddo
      enddo
      enddo
      enddo
      enddo


      part4=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      string=g0(fi)*g0(nu)*g0(ro)*g0(si)*J52x4(fi,nu,ro,si)
      bit= + prt(s345)**2*J34x1(fi)*J61x1(si)*facuRl*mt * ( FB1x345(nu,
     &    ep)*p345(ro) + FB1x345(ro,ep)*p345(nu) )
      bit = bit + J34x1(nu)*J61x1(si)*facuLl*mt * (  - FD2x2x5x34(y2(fi
     &    ,ro),ep) - FD1x2x5x34(fi,ep)*p2345(ro) )
      bit = bit + J34x1(nu)*J61x1(fi)*facuRl*mt * (  - FD2x2x5x34(y2(ro
     &    ,si),ep) - FD1x2x5x34(si,ep)*p25(ro) )
      bit = bit + J34x1(ro)*J61x1(si)*facuRl*mt * (  - FD2x2x5x34(y2(fi
     &    ,nu),ep) - FD1x2x5x34(fi,ep)*p25(nu) )
      bit = bit + J34x1(ro)*J61x1(fi)*facuLl*mt * (  - FD2x2x5x34(y2(nu
     &    ,si),ep) - FD1x2x5x34(si,ep)*p2345(nu) )
      bit = bit + J34x1(si)*J61x1(nu)*facdLl*mt * (  - FD2x34x2x5(y2(fi
     &    ,ro),ep) - FD1x34x2x5(ro,ep)*p34(fi) )
      bit = bit + J34x1(fi)*J61x1(ro)*facdLl*mt * (  - FD2x34x2x5(y2(nu
     &    ,si),ep) - FD1x34x2x5(nu,ep)*p34(si) )

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
      bit= + prt(s345)**2*J34x1(fi)*J61x1(ro)*facuLl*mtsq * ( FB1x345(
     &    nu,ep) + Zm(ep)*p345(nu) - 2.D0*FB0x345(ep)*p345(nu) )
      bit = bit + prt(s345)*J34x1(fi)*J61x1(ro)*facuLl * ( FB0x345(ep)*
     &    p345(nu) )
      bit = bit + J34x1(ro)*J61x1(nu)*facuRl*mtsq * ( FD1x2x5x34(fi,ep)
     &     )
      bit = bit + J34x1(ro)*J61x1(fi)*facdLl * (  - prq(s234)*FB0x234(
     &    ep)*p234(nu) )

      part3=part3+bit*string
      enddo
      enddo
      enddo


      part2=czip
      do fi=1,4
      do nu=1,4
      string=g0(fi)*g0(nu)*J52x2(fi,nu)
      bit= + prt(s345)**2*J34x1(fi)*J61x1(nu)*facuRl*mt*mtsq * ( Zm(ep)
     &     - 2.D0*FB0x345(ep) )
      bit = bit + prt(s345)*J34x1(fi)*J61x1(nu)*facuRl*mt * ( 1.D0/2.D0
     &    *Zm(ep) )

      part2=part2+bit*string
      enddo
      enddo


      part1=czip
      do fi=1,4
      string=g0(fi)*J52x1(fi)
      enddo

      lowerbox(ep)=prW(s16)*(part1+part2+part3+part4+part5)*cprop

      enddo


      return
      end
