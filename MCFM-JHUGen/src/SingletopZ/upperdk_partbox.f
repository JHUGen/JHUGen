      subroutine upperdk_partbox(p,h3,upperbox,first)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'metric0.f'
      include 'alpha1.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'TRydef.f' 
      include 'currentdecl.f'
      include 'tensordecl.f'
      include 'decl_kininv.f'
      integer:: fi,nu,ro,si,om,ep,h3
      complex(dp):: prW,prq,string,upperbox(-2:0),bit,iprZ,
     & part5,part3,part1,cprop
      complex(dp):: facuLl,facdLl,facLdiff
       real(dp):: p(mxpart,4),mwsq,
     & omal,p1(4),p6(4)
      integer:: epmin
      logical:: first

c     statement function
      prW(s34)=cone/cplx1(s34-wmass**2)
      prq(s34)=cone/cplx1(s34)
c     end statement function

      omal=1d0-alpha1

      mwsq=wmass**2
      p1(:)= p(1,:)
      p6(:)= p(6,:)

c--- Implementation of Baur-Zeppenfeld treatment of Z width
      cprop=cplx1(1d0/sqrt((s34-zmass**2)**2+(zmass*zwidth)**2))
      iprZ=cplx1(s34-zmass**2)

      cprop=cprop/cplx2(zip,mt*twidth)

      if (h3 == -1) then
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*le)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*le)
      else
        facuLl=cplx1(Qu*q1)*iprZ/s34+cplx1(L(2)*re)
        facdLl=cplx1(Qd*q1)*iprZ/s34+cplx1(L(1)*re)
      endif
      facLdiff=facuLl-facdLl

      call doBtensor(p134,0d0,0d0,FB0x134,FB1x134,
     &  FB2x134,FB3x134,FB4x134,FB5x134,FB6x134)
      call doBtensor(p346,0d0,0d0,FB0x346,FB1x346,
     &  FB2x346,FB3x346,FB4x346,FB5x346,FB6x346)
       
      call doCtensor(p6,p346,0d0,0d0,0d0,FC0x6x34,FC1x6x34,
     &  FC2x6x34,FC3x6x34,FC4x6x34,FC5x6x34,FC6x6x34)
      call doCtensor(p34,p134,0d0,0d0,0d0,FC0x34x1,FC1x34x1,
     &  FC2x34x1,FC3x34x1,FC4x34x1,FC5x34x1,FC6x34x1)
      call doCtensor(p1,p16,0d0,0d0,0d0,FC0x1x6,FC1x1x6,
     &  FC2x1x6,FC3x1x6,FC4x1x6,FC5x1x6,FC6x1x6)
       
      call doDtensor(p1,p16,p1346,0d0,0d0,0d0,0d0,FD0x1x6x34,
     &FD1x1x6x34,FD2x1x6x34,FD3x1x6x34,FD4x1x6x34,FD5x1x6x34,FD6x1x6x34)
      call doDtensor(p34,p134,p1346,0d0,0d0,0d0,0d0,FD0x34x1x6,
     &FD1x34x1x6,FD2x34x1x6,FD3x34x1x6,FD4x34x1x6,FD5x34x1x6,FD6x34x1x6)
       
c--- only compute poles for checking on first call
      if (first) then
         epmin=-2
      else
         epmin=-1
      endif

c      write(*,*) 'epmin in upperdk_box', epmin

      do ep=epmin,0


      part5=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      do si=1,4
      do om=1,4
      string=g0(fi)*g0(nu)*g0(ro)*g0(si)*g0(om)*J61x5(fi,nu,ro,si,om)
      bit= + prq(s134)**2*J34x1(om)*J52x1(fi)*facuLl * ( FB1x134(ro,ep)
     &    *p134(nu)*p134(si) )
      bit = bit + prq(s346)**2*J34x1(fi)*J52x1(om)*facdLl * ( FB1x346(
     &    ro,ep)*p346(nu)*p346(si) )
      bit = bit + J34x1(nu)*J52x1(ro)*facdLl * (  - FD3x1x6x34(y3(fi,si
     &    ,om),ep) - FD2x1x6x34(y2(si,om),ep)*p16(fi) - FD2x1x6x34(y2(
     &    fi,si),ep)*p1346(om) - FD1x1x6x34(si,ep)*p16(fi)*p1346(om) )
      bit = bit + J34x1(nu)*J52x1(si)*facdLl * ( FD3x1x6x34(y3(fi,ro,om
     &    ),ep) + FD2x1x6x34(y2(ro,om),ep)*p16(fi) + FD2x1x6x34(y2(fi,
     &    om),ep)*p1346(ro) + FD1x1x6x34(om,ep)*p16(fi)*p1346(ro) )
      bit = bit + J34x1(nu)*J52x1(om)*facdLl * ( FD3x1x6x34(y3(fi,ro,si
     &    ),ep) + FD2x1x6x34(y2(ro,si),ep)*p16(fi) + FD2x1x6x34(y2(fi,
     &    si),ep)*p1346(ro) + FD1x1x6x34(si,ep)*p16(fi)*p1346(ro) )
      bit = bit + J34x1(ro)*J52x1(nu)*facuLl * (  - FD3x34x1x6(y3(fi,si
     &    ,om),ep) - FD2x34x1x6(y2(si,om),ep)*p1346(fi) - FD2x34x1x6(
     &    y2(fi,om),ep)*p34(si) - FD1x34x1x6(om,ep)*p34(si)*p1346(fi) )
      bit = bit + J34x1(ro)*J52x1(fi)*facuLl * (  - FD3x34x1x6(y3(nu,si
     &    ,om),ep) - FD2x34x1x6(y2(nu,si),ep)*p1346(om) - FD2x34x1x6(
     &    y2(nu,om),ep)*p34(si) - FD1x34x1x6(nu,ep)*p34(si)*p1346(om) )
      bit = bit + J34x1(ro)*J52x1(om)*facuLl * ( FD3x34x1x6(y3(fi,nu,si
     &    ),ep) + FD2x34x1x6(y2(nu,si),ep)*p1346(fi) + FD2x34x1x6(y2(fi
     &    ,nu),ep)*p34(si) + FD1x34x1x6(nu,ep)*p34(si)*p1346(fi) )
      bit = bit + J34x1(si)*J52x1(nu)*facuLl * ( FD3x34x1x6(y3(fi,ro,om
     &    ),ep) + FD2x34x1x6(y2(ro,om),ep)*p1346(fi) + FD2x34x1x6(y2(fi
     &    ,ro),ep)*p34(om) + FD1x34x1x6(ro,ep)*p34(om)*p1346(fi) )
      bit = bit + J34x1(fi)*J52x1(ro)*facdLl * (  - FD3x1x6x34(y3(nu,si
     &    ,om),ep) - FD2x1x6x34(y2(nu,si),ep)*p16(om) - FD2x1x6x34(y2(
     &    si,om),ep)*p1346(nu) - FD1x1x6x34(si,ep)*p16(om)*p1346(nu) )
      bit = bit + J34x1(om)*J52x1(nu)*facuLl * ( FD3x34x1x6(y3(fi,ro,si
     &    ),ep) + FD2x34x1x6(y2(ro,si),ep)*p1346(fi) + FD2x34x1x6(y2(fi
     &    ,ro),ep)*p34(si) + FD1x34x1x6(ro,ep)*p34(si)*p1346(fi) )
      bit = bit + J34x1(om)*J52x1(ro)*facdLl * ( FD3x1x6x34(y3(fi,nu,si
     &    ),ep) + FD2x1x6x34(y2(nu,si),ep)*p16(fi) + FD2x1x6x34(y2(fi,
     &    si),ep)*p1346(nu) + FD1x1x6x34(si,ep)*p16(fi)*p1346(nu) )

      part5=part5+bit*string
      enddo
      enddo
      enddo
      enddo
      enddo


      part3=czip
      do fi=1,4
      do nu=1,4
      do ro=1,4
      string=g0(fi)*g0(nu)*g0(ro)*J61x3(fi,nu,ro)
      bit= + prq(s134)*J34x1(ro)*J52x0*omal*mwsq**(-1)*mt * ( FB1x134(
     &    fi,ep)*p134(nu)*facLdiff )
      bit = bit + prq(s134)*J34x1(ro)*J52x0*omal*mwsq**(-1)*facdLl*mt
     &  * ( FB1x134(fi,ep)*p134(nu) )
      bit = bit + prq(s346)*J34x1(fi)*J52x0*omal*mwsq**(-1)*facdLl*mt
     &  * ( FB1x346(ro,ep)*p346(nu) )
      bit = bit + prq(s346)*J34x1(fi)*J52x1(ro)*facdLl * ( FB0x346(ep)*
     &    p346(nu) )
      bit = bit + J34x1(nu)*J52x0*omal*mwsq**(-1)*mt * ( FC2x34x1(y2(fi
     &    ,ro),ep)*facLdiff - FC2x1x6(y2(fi,ro),ep)*facLdiff + 
     &    FC1x34x1(ro,ep)*p34(fi)*facLdiff - FC1x1x6(fi,ep)*p16(ro)*
     &    facLdiff )
      bit = bit + J34x1(nu)*J52x0*omal*mwsq**(-1)*facdLl*mt * (  - 
     &    FC2x6x34(y2(fi,ro),ep) + FC2x34x1(y2(fi,ro),ep) - FC1x6x34(ro
     &    ,ep)*p346(fi) - FC1x6x34(fi,ep)*p6(ro) + FC1x34x1(ro,ep)*
     &    p34(fi) - FC0x6x34(ep)*p6(ro)*p346(fi) )

      part3=part3+bit*string
      enddo
      enddo
      enddo


      part1=czip
      do fi=1,4
      string=g0(fi)*J61x1(fi)
      bit= + J34x1(fi)*J52x0*omal*mwsq**(-1)*facdLl*mt * ( FB0x346(ep)
     &     )

      part1=part1+bit*string

      enddo

      upperbox(ep)=prW(s25)*(part1+part3+part5)*cprop

      enddo


      return
      end
