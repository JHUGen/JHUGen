MODULE ModTopDecay
implicit none

public :: TopDecay
private

 CONTAINS





SUBROUTINE TopDecay(Flavor,Mom,Spinor,TopHel)
use ModMisc
use ModParameters
implicit none
real(8) :: Mom(1:4,1:3)
integer :: flavor
integer,optional :: TopHel
complex(8) :: Spinor(1:4)
real(8) :: zeros(1:6),TopMom(1:4),NWAFactor_Top
complex(8) :: Spi(1:4),BarSpi(1:4),BotSpi(1:4),WCurr(1:4)
real(8) :: NWAFactor_W
complex(8) :: WProp
real(8),parameter :: Nc=xn,NFlav=2

NWAFactor_W   = 1d0/dsqrt(2d0*Ga_W*m_W)
WProp = (0d0,-1d0)*NWAFactor_W



!     zeros(1:4) = dble(TopMom(1:4)) - Mom(1:4,1)-Mom(1:4,2)-Mom(1:4,3)
!     if( any(abs(zeros(1:4)/dble(TopMom(1))).gt.1d-8) ) then
!         print *, "ERROR: energy-momentum violation in SUBROUTINE TopDecay(): ",zeros(1:4)
!         print *, "momentum dump:"
!         print *, Mom(1:4,1:12)
!     endif
!
!     zeros(1) = dble(TopMom(1:4).dot.TopMom(1:4)) - m_Top**2
!     zeros(2) = Mom(1:4,1).dot.Mom(1:4,1)
!     zeros(3) = Mom(1:4,2).dot.Mom(1:4,2)
!     zeros(4) = Mom(1:4,3).dot.Mom(1:4,3)
!     if( any(abs(zeros(1:6)/dble(TopMom(1))**2).gt.1d-8) ) then
!         print *, "ERROR: onshell-ness violation in SUBROUTINE TopDecay(): ",zeros(1:6)
!         print *, "momentum dump:"
!         print *, Mom(1:4,1:12)
!     endif



    TopMom(1:4) = Mom(1:4,1)+Mom(1:4,2)+Mom(1:4,3)
    if( TOPDECAYS.eq.0 ) then
!         if( Flavor.eq.Top_  ) call ubarSpi_Dirac(dcmplx(TopMom(1:4)),M_Top,TopHel,Spinor(1:4))
!         if( Flavor.eq.ATop_ ) call    vSpi_Dirac(dcmplx(TopMom(1:4)),M_Top,TopHel,Spinor(1:4))
        Spinor(1:4) = 0d0
        RETURN
    endif
!     elseif( TOPDECAYS.eq.1 ) then
    NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)
!     elseif( TOPDECAYS.eq.2 ) then
!         NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)
!         NWAFactor_Top = NWAFactor_Top * dsqrt(dsqrt(Nc*NFlav)**2)
!     elseif( TOPDECAYS.eq.3 .or. TOPDECAYS.eq.4 ) then
!         NWAFactor_Top = 1d0/dsqrt(2d0*Ga_Top*m_Top)
!         NWAFactor_Top = NWAFactor_Top * dsqrt(dsqrt(Nc*NFlav))
!     endif


    if( Flavor.eq.Top_ ) then ! Top quark decay
!       assemble lepton current
        call ubarSpi_Weyl(dcmplx(Mom(1:4,1)),-1,BotSpi(1:4))  ! bot
        call    vSpi_Weyl(dcmplx(Mom(1:4,2)),+1,Spi(1:4))     ! l+ or dn_bar
        call ubarSpi_Weyl(dcmplx(Mom(1:4,3)),-1,BarSpi(1:4))  ! nu or up
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp * gwsq ! vbqq introduces -i/Sqrt(2)

!       connect to quark current
        BarSpi(1:4) = BotSpi(1:4)
        Spinor(1:4) = vgq_Weyl( WCurr(1:4),BarSpi(1:4) ) ! vgq introduces -i/Sqrt(2)
        Spinor(1:4) =( spb2_Weyl(Spinor(1:4),dcmplx(TopMom(1:4))) + m_Top*Spinor(1:4) ) * NWAFactor_Top
        Spinor(1:4) = WeylToDirac(Spinor(1:4))
    elseif( Flavor.eq.ATop_ ) then ! Anti-Top quark decay
!       assemble lepton current
        call    vSpi_Weyl(dcmplx(Mom(1:4,1)),+1,BotSpi(1:4))  ! Abot
        call ubarSpi_Weyl(dcmplx(Mom(1:4,2)),-1,BarSpi(1:4))  ! l- or dn
        call    vSpi_Weyl(dcmplx(Mom(1:4,3)),+1,Spi(1:4))     ! nubar or up_bar
        WCurr(1:4)  = vbqq_Weyl(BarSpi(1:4),Spi(1:4)) * WProp * gwsq ! vbqq introduces -i/Sqrt(2)

!       connect to quark current:
        Spi(1:4) = BotSpi(1:4)
        Spinor(1:4) = vbqg_Weyl( Spi(1:4),WCurr(1:4) )! vbqg introduces -i/Sqrt(2)
        Spinor(1:4) = ( spi2_Weyl(dcmplx(TopMom(1:4)),Spinor(1:4)) - m_Top*Spinor(1:4) ) * NWAFactor_Top
        Spinor(1:4) = WeylToDirac(Spinor(1:4))
    endif


RETURN
END SUBROUTINE





          SUBROUTINE ubarSpi_Weyl(p,i,ubarSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          use modMisc
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: ubarSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          complex(8) :: fc, fc2

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            ubarSpi(1)=(0d0,0d0)
            ubarSpi(2)=(0d0,0d0)
            ubarSpi(3)=fc
            ubarSpi(4)=(px-(0d0,1d0)*py)/fc
         elseif (i.eq.-1) then
            ubarSpi(1)=(px+(0d0,1d0)*py)/fc
            ubarSpi(2)=-fc
            ubarSpi(3)=(0d0,0d0)
            ubarSpi(4)=(0d0,0d0)
         else
          call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            ubarSpi(1) = (0d0,0d0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = dsqrt(2d0*p0)
         elseif (i.eq.-1) then
            ubarSpi(1) = dsqrt(2d0*p0)
            ubarSpi(2) = (0d0,0d0)
            ubarSpi(3) = (0d0,0d0)
            ubarSpi(4) = (0d0,0d0)
         else
            call Error("wrong helicity setting in ubarSpi_Weyl")
         endif

         endif
        return
        END SUBROUTINE





          SUBROUTINE vSpi_Weyl(p,i,vSpi)  ! i=+1 is ES to Chir_Weyl(.false.), i=-1 is ES to Chir_Weyl(.true.)
          use ModMisc
          implicit none
          integer, intent(in):: i
          complex(8), intent(in) :: p(4)
          complex(8) :: vSpi(4)
          complex(8) :: ephi
          real(8) :: p0,px,py,pz
          real(8) :: nx,ny,nz,theta,phi
          real(8) :: ct,ct2,st,st2,cphi,sphi
          complex(8) :: fc2, fc

          p0=dreal(p(1))
          px=dreal(p(2))
          py=dreal(p(3))
          pz=dreal(p(4))

         fc2 = p0 + pz
         fc=cdsqrt(fc2)

         if (cdabs(fc2).gt.1D-15) then

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=(px-(0d0,1d0)*py)/fc
            vSpi(4)=-fc
         elseif (i.eq.-1) then
            vSpi(1)=fc
            vSpi(2)=(px+(0d0,1d0)*py)/fc
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         else

         if (i.eq.1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=(0d0,0d0)
            vSpi(3)=dsqrt(2d0*p0)
            vSpi(4)=(0d0,0d0)
         elseif (i.eq.-1) then
            vSpi(1)=(0d0,0d0)
            vSpi(2)=dsqrt(2d0*p0)
            vSpi(3)=(0d0,0d0)
            vSpi(4)=(0d0,0d0)
         else
            call Error("wrong helicity setting in vSpi_Weyl")
         endif

         endif

         RETURN
         END SUBROUTINE








        FUNCTION vbqg_Weyl(sp,e1)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vbqg_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vbqg_Weyl = (0d0,-1d0)/sqrt2*spi2_Weyl(e1,sp)

        END FUNCTION




        FUNCTION vbqq_Weyl(sp1,sp2)
        implicit none
        complex(8), intent(in) :: sp1(:), sp2(:)
        integer, parameter ::  Dv=4
        integer :: i
        complex(8) :: vbqq_Weyl(Dv)
        complex(8) :: rr, va(Dv),sp1a(size(sp1))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            va=(0d0,0d0)
            vbqq_Weyl=(0d0,0d0)

            do i=1,Dv
              if (i.eq.1) then
                va(1)=(1d0,0d0)
              else
                va(i)=(-1d0,0d0)
              endif
              sp1a=spb2_Weyl(sp1,va)

              rr=(0d0,-1d0)/sqrt2*psp1_(sp1a,sp2)
              if (i.eq.1) then
                    vbqq_Weyl = vbqq_Weyl + rr*va
                else
                    vbqq_Weyl = vbqq_Weyl - rr*va
              endif
              va(i)=(0d0,0d0)
            enddo

        END FUNCTION



        function vgq_Weyl(e1,sp)
        implicit none
        complex(8), intent(in) :: e1(:)
        complex(8), intent(in) :: sp(:)
        complex(8) :: vgq_Weyl(size(sp))
        real(8), parameter :: sqrt2 = 1.4142135623730950488016887242096980786d0

            vgq_Weyl = (0d0,-1d0)/sqrt2*spb2_Weyl(sp,e1)

        end function





      function spb2_Weyl(sp,v)
      implicit none
      complex(8), intent(in) :: sp(:),v(:)
      integer, parameter ::  Dv=4, Ds=4
      complex (8) :: spb2_Weyl(size(sp))
      complex(8) :: x0(4,4),xx(4,4),xy(4,4)
      complex(8) :: xz(4,4),x5(4,4)
      complex(8) :: y1,y2,y3,y4,bp,bm,cp,cm
      integer :: i,i1,i2,i3,imax



      imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2

           xx(1,i) = y4
           xx(2,i) = y3
           xx(3,i) = -y2
           xx(4,i) = -y1

           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=y3
           xz(2,i)=-y4
           xz(3,i)=-y1
           xz(4,i)=y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo



           do i=1,4

           spb2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1)-v(3)*xy(i,1)-v(4)*xz(i,1)

           enddo


           end function





         function spi2_Weyl(v,sp)
         implicit none
         complex(8), intent(in) :: sp(:),v(:)
         complex(8) :: spi2_Weyl(size(sp))
         integer, parameter ::  Dv=4,Ds=4
         complex(8) :: x0(4,4),xx(4,4),xy(4,4)
         complex(8) :: xz(4,4),x5(4,4)
         complex(8) ::  y1,y2,y3,y4,bp,bm,cp,cm
         integer :: i,i1,i2,i3,imax


         imax = Ds/4

           do i=1,imax
           i1= 1+4*(i-1)
           i2=i1+3

           y1=sp(i1)
           y2=sp(i1+1)
           y3=sp(i1+2)
           y4=sp(i1+3)

           x0(1,i)=y3
           x0(2,i)=y4
           x0(3,i)=y1
           x0(4,i)=y2


           xx(1,i) = -y4
           xx(2,i) = -y3
           xx(3,i) = y2
           xx(4,i) = y1


           xy(1,i)=(0d0,1d0)*y4
           xy(2,i)=-(0d0,1d0)*y3
           xy(3,i)=-(0d0,1d0)*y2
           xy(4,i)=(0d0,1d0)*y1

           xz(1,i)=-y3
           xz(2,i)=y4
           xz(3,i)=y1
           xz(4,i)=-y2

           x5(1,i)=y1
           x5(2,i)=y2
           x5(3,i)=-y3
           x5(4,i)=-y4

           enddo


           do i=1,4

           spi2_Weyl(i)=v(1)*x0(i,1)-v(2)*xx(i,1) -v(3)*xy(i,1)-v(4)*xz(i,1)
           enddo


           end function



          function WeylToDirac(sp)   ! unitary transformation U to convert a Weyl spinor into the Dirac representation
          implicit none              ! sp can be spinor or bar-spinor, i.e. U^dagger.sp = barsp.U
          double complex :: sp(1:4)
          double complex :: WeylToDirac(1:4)
          double precision,parameter :: SqrtFac=1d0/dsqrt(2d0)

              WeylToDirac(1) = SqrtFac*(sp(1)+sp(3))
              WeylToDirac(2) = SqrtFac*(sp(2)+sp(4))
              WeylToDirac(3) = SqrtFac*(sp(1)-sp(3))
              WeylToDirac(4) = SqrtFac*(sp(2)-sp(4))
          return
          end function




          function psp1_(sp1,sp2) result(res)
          implicit none
          complex(8), intent(in) :: sp1(:)
          complex(8), intent(in) :: sp2(:)
          complex(8) :: res

            res = sum(sp1(1:)*sp2(1:))

           end function





END MODULE
