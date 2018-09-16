*
* Note added: July 2015, J. Campbell
* The files in this directory contain the subroutines present in the
* TDHPL package (described below);  they have been split over several
* files to assist compilation on some machines.  In addition, type
* declarations and OMP directives have been added in order to function
* properly with the multi-thread version of MCFM.  A single-file version
* has been included here to aid comparison with the original package.
*
******************************************************************************
**  tdhpl:  a subroutine for the evaluation of
**          two-dimensional harmonic polylogarithms
**          Version 1.0         09/11/2001
**  described in:
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of Two-Dimensional
**                            Harmonic Polylogarithms
**                            (hep-ph/0111255; CERN-TH/2001/326)
**  the harmonic polylogarithms are defined in:
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  and:
**  T.Gehrmann and E.Remiddi: Two-Loop Master Integrals for gamma*->3jets:
**                            the Planar Topologies
**                            (hep-ph/0008287; Nucl.Phys. B601(2001) 248)
**  this subroutine invokes the subroutine hplog, which is described in:
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of
**                            Harmonic Polylogarithms
**                            (hep-ph/0107173; Comp.Phys.Comm. 141(2001) 296)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
******************************************************************************
      subroutine tdhpl(y,z,nmax,GYZ1,GYZ2,GYZ3,GYZ4,
     $                          HZ1,HZ2,HZ3,HZ4)
*********************************************************************
***** Input:                                                    *****
***** y and z are the arguments of the 2dHPL to be evaluated    *****
***** nmax    is the maximum weight to be evaluated             *****
***** Output:                                                   *****
***** GYZ1,GYZ2,GYZ3,GYZ4 are the double precision 2dHPL of y,z *****
***** HZ1,HZ2,HZ3,HZ4     are the double precision 1dHPL of z   *****
*********************************************************************
      implicit none
      include 'types.f'
      integer i1,i2,i3,i4
      integer sgi1,sgi2,sgi3,sgi4
      integer iflag,nmax
      real(dp):: x,y,z
      real(dp):: zold
      real(dp):: HZ1,HZ2,HZ3,HZ4,GYZ1,GYZ2,GYZ3,GYZ4
      real(dp):: HZ1a,HZ2a,HZ3a,HZ4a,HYZ1a,HYZ2a,HYZ3a,HYZ4a
      dimension HZ1(0:1),HZ2(0:1,0:1),HZ3(0:1,0:1,0:1),
     $          HZ4(0:1,0:1,0:1,0:1)
      dimension GYZ1(0:3),GYZ2(0:3,0:3),GYZ3(0:3,0:3,0:3),
     $          GYZ4(0:3,0:3,0:3,0:3)
      common/HPL2OUT/
     $ HZ1a(0:1),HZ2a(0:1,0:1),HZ3a(0:1,0:1,0:1),
     $           HZ4a(0:1,0:1,0:1,0:1),
     $ HYZ1a(0:3),HYZ2a(0:3,0:3),HYZ3a(0:3,0:3,0:3),
     $            HYZ4a(0:3,0:3,0:3,0:3)
      data zold/ -999d0/
      save zold
!$omp threadprivate(zold,/hpl2out/)


      x=1d0-y-z
      if (x.lt.0d0.or.x.gt.1.d0) then
         write(6,991) y,z
         stop
      endif
      if (y.lt.0d0.or.y.gt.1.d0) then
         write(6,991) y,z
         stop
      endif
      if (z.lt.0d0.or.z.gt.1.d0) then
         write(6,991) y,z
         stop
      endif
      if (nmax.lt.1.or.nmax.gt.4) then
         write(6,992) nmax
         stop
      endif

      iflag=1
      if (z.eq.zold) iflag=2
      zold = z
*      write(6,*) iflag
      call fill2dhpl(iflag,nmax,y,z)

      do i1=0,1
         if (nmax.ge.2) then
         do i2=0,1
            if (nmax.ge.3) then
            do i3=0,1
               if (nmax.eq.4) then
               do i4=0,1
                  HZ4(i1,i2,i3,i4) = HZ4a(i1,i2,i3,i4)
               enddo
               endif
               HZ3(i1,i2,i3) = HZ3a(i1,i2,i3)
            enddo
            endif
            HZ2(i1,i2) = HZ2a(i1,i2)
         enddo
         endif
         HZ1(i1) = HZ1a(i1)
      enddo
      do i1=0,3
         if (i1.eq.1.or.i1.eq.2) then
            sgi1=-1
         else
            sgi1=1
         endif
         if (nmax.ge.2) then
         do i2=0,3
            if (i2.eq.1.or.i2.eq.2) then
               sgi2=-1
            else
               sgi2=1
            endif
            if (nmax.ge.3) then
            do i3=0,3
               if (i3.eq.1.or.i3.eq.2) then
                  sgi3=-1
               else
                  sgi3=1
               endif
               if (nmax.eq.4) then
               do i4=0,3
                  if (i4.eq.1.or.i4.eq.2) then
                     sgi4=-1
                  else
                     sgi4=1
                  endif
                  GYZ4(i1,i2,i3,i4) = sgi1*sgi2*sgi3*sgi4*
     $                                    HYZ4a(i1,i2,i3,i4)
               enddo
               endif
               GYZ3(i1,i2,i3) = sgi1*sgi2*sgi3*HYZ3a(i1,i2,i3)
            enddo
            endif
            GYZ2(i1,i2) = sgi1*sgi2*HYZ2a(i1,i2)
         enddo
         endif
         GYZ1(i1) = sgi1*HYZ1a(i1)
      enddo
      return
 991  format('Parameters out of range, y=',f9.5,', z=',f9.5)
 992  format('Weight out of range, nmax=',i3)
      end


      subroutine fill2dhpl(iflag,nmax,y,z)
*********************************************************************
***** fill2dhpl evaluates all 2dhpl writes all                  *****
***** 2dhpl(y) and 1dhpl(z)                                     *****
***** up to weight nmax into HPL2OUT                            *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      integer i1,i2,i3,i4
      real(dp):: y,z,xi
      real(dp):: HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HZE1,HZE2,HZE3,HZE4
      dimension HZE1(-1:1),HZE2(-1:1,-1:1),HZE3(-1:1,-1:1,-1:1),
     $          HZE4(-1:1,-1:1,-1:1,-1:1)
      common/HPL2OUT/
     $ HZ1(0:1),HZ2(0:1,0:1),HZ3(0:1,0:1,0:1),
     $          HZ4(0:1,0:1,0:1,0:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      save HZE1,HZE2,HZE3,HZE4
!$omp threadprivate(HZE1,HZE2,HZE3,HZE4,/HPL2OUT/)
      if (iflag.eq.1) then
      call get1dhplog0m11(iflag,nmax,z,HZE1,HZE2,HZE3,HZE4)
      endif
      do i1=0,1
         if (nmax.ge.2) then
         do i2=0,1
            if (nmax.ge.3) then
            do i3=0,1
               if (nmax.eq.4) then
               do i4=0,1
                  HZ4(i1,i2,i3,i4) = HZE4(i1,i2,i3,i4)
               enddo
               endif
               HZ3(i1,i2,i3) = HZE3(i1,i2,i3)
            enddo
            endif
            HZ2(i1,i2) = HZE2(i1,i2)
         enddo
         endif
         HZ1(i1) = HZE1(i1)
      enddo
      xi = y/(1d0-z)
      if (xi.le.0.5d0) then
         call get2dhplat0(iflag,nmax,y,z,HYZ1,HYZ2,HYZ3,HYZ4,
     $                                   HZE1,HZE2,HZE3,HZE4)
*         write(6,*) 'call get2dhplat0'
      else
         call get2dhplat1(iflag,nmax,y,z,HYZ1,HYZ2,HYZ3,HYZ4,
     $                                   HZE1,HZE2,HZE3,HZE4)
*         write(6,*) 'call get2dhplat1'
      endif

      return
      end

      subroutine get1dhplog01(iflag,nw,z,HZ1,HZ2,HZ3,HZ4)
*********************************************************************
***** get1dhplog01 returns the 1dhpl (0,1) of z                 *****
***** up to weight nw                                           *****
***** evaluation using hplog                                    *****
*********************************************************************
      implicit none
      include 'types.f'
      integer n1,n2,nw,i1,i2,iflag
      parameter (n1=0)
      parameter (n2=1)
      complex(dp):: Hc1,Hc2,Hc3,Hc4
      real(dp)::     HZ1,HZ2,HZ3,HZ4
      real(dp)::     Hi1,Hi2,Hi3,Hi4
      real(dp)::     z
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2),
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension HZ1(n1:n2),HZ2(n1:n2,n1:n2),HZ3(n1:n2,n1:n2,n1:n2),
     $          HZ4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      call hplog(z,nw,Hc1,Hc2,Hc3,Hc4,
     $                HZ1,HZ2,HZ3,HZ4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      return
      end


      subroutine get1dhplog0m11(iflag,nw,z,HZ1,HZ2,HZ3,HZ4)
*********************************************************************
***** get1dhplog0m11 returns the 1dhpl (-1,0,1) of z            *****
***** up to weight nw                                           *****
***** evaluation using hplog                                    *****
*********************************************************************
      implicit none
      include 'types.f'
      integer n1,n2,nw,i1,i2,iflag
      parameter (n1=-1)
      parameter (n2=1)
      complex(dp):: Hc1,Hc2,Hc3,Hc4
      real(dp)::     HZ1,HZ2,HZ3,HZ4
      real(dp)::     Hi1,Hi2,Hi3,Hi4
      real(dp)::     z
      dimension Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2),
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension HZ1(n1:n2),HZ2(n1:n2,n1:n2),HZ3(n1:n2,n1:n2,n1:n2),
     $          HZ4(n1:n2,n1:n2,n1:n2,n1:n2)
      dimension Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2),
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      call hplog(z,nw,Hc1,Hc2,Hc3,Hc4,
     $                HZ1,HZ2,HZ3,HZ4,Hi1,Hi2,Hi3,Hi4,n1,n2)
      return
      end


      subroutine get2dhplat1(iflag,nmax,y,z,HYZ11,HYZ12,HYZ13,HYZ14
     $                       ,HZE1,HZE2,HZE3,HZE4)
*********************************************************************
***** get2dhplat0 evaluates all 2dhpl from expansions           *****
***** for small 1-y-z=x <(1-z)/2                                *****
***** and subsequent transformation                             *****
***** up to weight nmax                                         *****
***** HZEn must be supplied (needed in cut-separation)          *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      real(dp):: x,y,z
      real(dp):: HZ1,HZ2,HZ3,HZ4,HYZ11,HYZ12,HYZ13,HYZ14
      real(dp):: HZE1,HZE2,HZE3,HZE4,HXZ01,HXZ02,HXZ03,HXZ04
      dimension HYZ11(0:3),HYZ12(0:3,0:3),HYZ13(0:3,0:3,0:3),
     $          HYZ14(0:3,0:3,0:3,0:3)
      dimension HXZ01(0:3),HXZ02(0:3,0:3),HXZ03(0:3,0:3,0:3),
     $          HXZ04(0:3,0:3,0:3,0:3)
      dimension HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $          HZ4(-1:1,-1:1,-1:1,-1:1)
      dimension HZE1(-1:1),HZE2(-1:1,-1:1),HZE3(-1:1,-1:1,-1:1),
     $          HZE4(-1:1,-1:1,-1:1,-1:1)

      x=1d0-y-z
      if (x.ge.y) then
         write(6,*) 'Error  1-y-z> y for eval. at 1'
         stop
      endif
      call get2dhplat0(iflag,nmax,x,z,HXZ01,HXZ02,HXZ03,HXZ04,
     $                                HZE1,HZE2,HZE3,HZE4)
      call swap2dhplxy(iflag,nmax,HXZ01,HXZ02,HXZ03,HXZ04,
     $                 HZE1,HZE2,HZE3,HZE4,HYZ11,HYZ12,HYZ13,HYZ14)
      return
      end


      subroutine get2dhplat0(iflag,nmax,y,z,HYZ01,HYZ02,HYZ03,HYZ04,
     $                                      HZE1,HZE2,HZE3,HZE4)
*********************************************************************
***** get2dhplat0 evaluates all 2dhpl from expansions           *****
***** for small y<(1-z)/2                                       *****
***** up to weight nmax                                         *****
***** HZEn must be supplied (needed in cut-separation)          *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax,n
      integer i1,i2,i3,i4
      real(dp):: y,z,xi
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HYZ01,HYZ02,HYZ03,HYZ04,HZE1,HZE2,HZE3,HZE4
      dimension HYZ01(0:3),HYZ02(0:3,0:3),HYZ03(0:3,0:3,0:3),
     $          HYZ04(0:3,0:3,0:3,0:3)
      dimension HZE1(-1:1),HZE2(-1:1,-1:1),HZE3(-1:1,-1:1,-1:1),
     $          HZE4(-1:1,-1:1,-1:1,-1:1)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
!$omp threadprivate(/HPL2/)

      do i1=-1,1
         if (nmax.ge.2) then
         do i2=-1,1
            if (nmax.ge.3) then
            do i3=-1,1
               if (nmax.eq.4) then
               do i4=-1,1
                  HZ4(i1,i2,i3,i4) = HZE4(i1,i2,i3,i4)
               enddo
               endif
               HZ3(i1,i2,i3) = HZE3(i1,i2,i3)
            enddo
            endif
            HZ2(i1,i2) = HZE2(i1,i2)
         enddo
         endif
         HZ1(i1) = HZE1(i1)
      enddo

      xi = y/(1d0-z)
      if (xi.gt.0.5d0) then
         write(6,*) 'Error y/(1-z) > 0.5 for eval. at 0'
         stop
      endif

      call fillirr2dhpl10(iflag,nmax,y,z)
      call fillirr2dhpl20(iflag,nmax,y,z)
      call fillirr2dhpl30(iflag,nmax,y,z)
      call makeexpar(y,z)

      do n=2,nmax
         call fillirr2dhpl210(iflag,n)
         call fillirr2dhpl310(iflag,n)
         if (z.le.0.5d0) then
            call fillirr2dhpl320(iflag,n)
            call fillirr2dhpl321(iflag,n)
         else
            call fillirr2dhpl320e(iflag,n)
            call fillirr2dhpl321e(iflag,n)
         endif
         call fillred2dhpl(iflag,n)
      enddo
      do i1=0,3
         if (nmax.ge.2) then
         do i2=0,3
            if (nmax.ge.3) then
            do i3=0,3
               if (nmax.eq.4) then
               do i4=0,3
                  HYZ04(i1,i2,i3,i4) = HYZ4(i1,i2,i3,i4)
               enddo
               endif
               HYZ03(i1,i2,i3) = HYZ3(i1,i2,i3)
            enddo
            endif
            HYZ02(i1,i2) = HYZ2(i1,i2)
         enddo
         endif
         HYZ01(i1) = HYZ1(i1)
      enddo
      return
      end


      subroutine fillirr2dhpl10(iflag,nmax,y,z)
*********************************************************************
***** fillirr2dhpl10 fills the irreducible 2dHPL                *****
***** with indices (1,0) up to weight n                         *****
***** by evaluating 1dHPL in y                                  *****
***** applicable for all z                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      integer i1,i2,i3,i4
      real(dp):: y,z
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HY1R,HY2R,HY3R,HY4R
      dimension HY1R(0:1),HY2R(0:1,0:1),HY3R(0:1,0:1,0:1),
     $          HY4R(0:1,0:1,0:1,0:1)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
!$omp threadprivate(/HPL2/)

      call get1dhplog01(iflag,nmax,y,HY1R,HY2R,HY3R,HY4R)
      do i1=0,1
         if (nmax.ge.2) then
         do i2=0,1
            if (nmax.ge.3) then
            do i3=0,1
               if (nmax.eq.4) then
               do i4=0,1
                  HYZ4(i1,i2,i3,i4) = HY4R(i1,i2,i3,i4)
                  HY4(i1,i2,i3,i4) = HY4R(i1,i2,i3,i4)
               enddo
               endif
               HYZ3(i1,i2,i3) = HY3R(i1,i2,i3)
               HY3(i1,i2,i3) = HY3R(i1,i2,i3)
            enddo
            endif
            HYZ2(i1,i2) = HY2R(i1,i2)
            HY2(i1,i2) = HY2R(i1,i2)
         enddo
         endif
         HYZ1(i1) = HY1R(i1)
         HY1(i1) = HY1R(i1)
      enddo
      return
      end


      subroutine fillirr2dhpl20(iflag,nmax,y,z)
*********************************************************************
***** fillirr2dhpl20 fills the irreducible 2dHPL                *****
***** with indices (2,0) up to weight n                         *****
***** by evaluating 1dHPL in xi=y/(1-z)                         *****
***** applicable for all z                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      integer i1,i2,i3,i4
      real(dp):: y,z,xi
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HXI1,HXI2,HXI3,HXI4
      dimension HXI1(0:1),HXI2(0:1,0:1),HXI3(0:1,0:1,0:1),
     $          HXI4(0:1,0:1,0:1,0:1)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
!$omp threadprivate(/HPL2/)

      xi = y/(1-z)
      call get1dhplog01(iflag,nmax,xi,HXI1,HXI2,HXI3,HXI4)
      i1=1
      if (nmax.ge.2) then
      do i2=0,1
         if (nmax.ge.3) then
         do i3=0,1
            if (nmax.eq.4) then
            do i4=0,1
               HYZ4(2*i4,2*i3,2*i2,2*i1) = HXI4(i4,i3,i2,i1)
            enddo
            endif
            HYZ3(2*i3,2*i2,2*i1) = HXI3(i3,i2,i1)
         enddo
         endif
         HYZ2(2*i2,2*i1) = HXI2(i2,i1)
      enddo
      endif
      HYZ1(2*i1) = HXI1(i1)
      return
      end


      subroutine fillirr2dhpl30(iflag,nmax,y,z)
*********************************************************************
***** fillirr2dhpl30 fills the irreducible 2dHPL                *****
***** with indices (3,0) up to weight n                         *****
***** by evaluating 1dHPL in xi=-y/z                            *****
***** applicable for all z                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      integer i1,i2,i3,i4
      real(dp):: y,z,xi
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HXI1,HXI2,HXI3,HXI4
      dimension HXI1(0:1),HXI2(0:1,0:1),HXI3(0:1,0:1,0:1),
     $          HXI4(0:1,0:1,0:1,0:1)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
!$omp threadprivate(/HPL2/)

      xi = -y/z
      call get1dhplog01(iflag,nmax,xi,HXI1,HXI2,HXI3,HXI4)
      i1=1
      if (nmax.ge.2) then
      do i2=0,1
         if (nmax.ge.3) then
         do i3=0,1
            if (nmax.eq.4) then
            do i4=0,1
                  HYZ4(3*i4,3*i3,3*i2,3*i1) = (-1)**(i4+i3+i2+i1)
     $                                         *HXI4(i4,i3,i2,i1)
            enddo
            endif
            HYZ3(3*i3,3*i2,3*i1) = (-1)**(i3+i2+i1)
     $                              *HXI3(i3,i2,i1)
         enddo
         endif
         HYZ2(3*i2,3*i1) = (-1)**(i2+i1)*HXI2(i2,i1)
      enddo
      endif
      HYZ1(3*i1) = (-1)**i1*HXI1(i1)
      return
      end


      subroutine fillirr2dhpl210(iflag,n)
*********************************************************************
***** fillirr2dhpl210(iflag,n) fills the irreducible 2dHPL      *****
***** with indices (2,1,0)                                      *****
***** up to weight n using the EXPANDED expressions for the     *****
***** z-dependent expansion coefficients                        *****
***** applicable for all z                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: CsHYZ2,CsHYZ3,CsHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      dimension CsHYZ2(0:3,0:3,22),CsHYZ3(0:3,0:3,0:3,22),
     $          CsHYZ4(0:3,0:3,0:3,0:3,22)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common /s/
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      save CsHYZ2,CsHYZ3,CsHYZ4
!$omp threadprivate(CsHYZ2,CsHYZ3,CsHYZ4,/s/,/HPL2/)

      if (iflag.eq.1) then
      call fillcoeff2dhpl210(iflag,n,CsHYZ2,CsHYZ3,CsHYZ4)
      endif
* 2001-04-03:16:14:02.hpl
* <- HPLexpand.210t
* produced by form-to-fortr for gehrt@pcth62

      if (n.eq.2) then
      HYZ2(2,1) =
     $  + s02 * CsHYZ2(2,1,02)
     $  + s03 * CsHYZ2(2,1,03)
     $  + s04 * CsHYZ2(2,1,04)
     $  + s05 * CsHYZ2(2,1,05)
     $  + s06 * CsHYZ2(2,1,06)
     $  + s07 * CsHYZ2(2,1,07)
     $  + s08 * CsHYZ2(2,1,08)
     $  + s09 * CsHYZ2(2,1,09)
     $  + s10 * CsHYZ2(2,1,10)
     $  + s11 * CsHYZ2(2,1,11)
     $  + s12 * CsHYZ2(2,1,12)
     $  + s13 * CsHYZ2(2,1,13)
     $  + s14 * CsHYZ2(2,1,14)
     $  + s15 * CsHYZ2(2,1,15)
     $  + s16 * CsHYZ2(2,1,16)
     $  + s17 * CsHYZ2(2,1,17)
     $  + s18 * CsHYZ2(2,1,18)
     $  + s19 * CsHYZ2(2,1,19)
     $  + s20 * CsHYZ2(2,1,20)
     $  + s21 * CsHYZ2(2,1,21)
      endif
      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(0,2,1) =
     $  + s02 * CsHYZ3(0,2,1,02)
     $  + s03 * CsHYZ3(0,2,1,03)
     $  + s04 * CsHYZ3(0,2,1,04)
     $  + s05 * CsHYZ3(0,2,1,05)
     $  + s06 * CsHYZ3(0,2,1,06)
     $  + s07 * CsHYZ3(0,2,1,07)
     $  + s08 * CsHYZ3(0,2,1,08)
     $  + s09 * CsHYZ3(0,2,1,09)
     $  + s10 * CsHYZ3(0,2,1,10)
     $  + s11 * CsHYZ3(0,2,1,11)
     $  + s12 * CsHYZ3(0,2,1,12)
     $  + s13 * CsHYZ3(0,2,1,13)
     $  + s14 * CsHYZ3(0,2,1,14)
     $  + s15 * CsHYZ3(0,2,1,15)
     $  + s16 * CsHYZ3(0,2,1,16)
     $  + s17 * CsHYZ3(0,2,1,17)
     $  + s18 * CsHYZ3(0,2,1,18)
     $  + s19 * CsHYZ3(0,2,1,19)
     $  + s20 * CsHYZ3(0,2,1,20)
      HYZ3(2,0,1) =
     $  + s02 * CsHYZ3(2,0,1,02)
     $  + s03 * CsHYZ3(2,0,1,03)
     $  + s04 * CsHYZ3(2,0,1,04)
     $  + s05 * CsHYZ3(2,0,1,05)
     $  + s06 * CsHYZ3(2,0,1,06)
     $  + s07 * CsHYZ3(2,0,1,07)
     $  + s08 * CsHYZ3(2,0,1,08)
     $  + s09 * CsHYZ3(2,0,1,09)
     $  + s10 * CsHYZ3(2,0,1,10)
     $  + s11 * CsHYZ3(2,0,1,11)
     $  + s12 * CsHYZ3(2,0,1,12)
     $  + s13 * CsHYZ3(2,0,1,13)
     $  + s14 * CsHYZ3(2,0,1,14)
     $  + s15 * CsHYZ3(2,0,1,15)
     $  + s16 * CsHYZ3(2,0,1,16)
     $  + s17 * CsHYZ3(2,0,1,17)
     $  + s18 * CsHYZ3(2,0,1,18)
     $  + s19 * CsHYZ3(2,0,1,19)
     $  + s20 * CsHYZ3(2,0,1,20)
      HYZ3(2,1,1) =
     $  + s03 * CsHYZ3(2,1,1,03)
     $  + s04 * CsHYZ3(2,1,1,04)
     $  + s05 * CsHYZ3(2,1,1,05)
     $  + s06 * CsHYZ3(2,1,1,06)
     $  + s07 * CsHYZ3(2,1,1,07)
     $  + s08 * CsHYZ3(2,1,1,08)
     $  + s09 * CsHYZ3(2,1,1,09)
     $  + s10 * CsHYZ3(2,1,1,10)
     $  + s11 * CsHYZ3(2,1,1,11)
     $  + s12 * CsHYZ3(2,1,1,12)
     $  + s13 * CsHYZ3(2,1,1,13)
     $  + s14 * CsHYZ3(2,1,1,14)
     $  + s15 * CsHYZ3(2,1,1,15)
     $  + s16 * CsHYZ3(2,1,1,16)
     $  + s17 * CsHYZ3(2,1,1,17)
     $  + s18 * CsHYZ3(2,1,1,18)
     $  + s19 * CsHYZ3(2,1,1,19)
     $  + s20 * CsHYZ3(2,1,1,20)
     $  + s21 * CsHYZ3(2,1,1,21)
     $  + s22 * CsHYZ3(2,1,1,22)
      HYZ3(2,2,1) =
     $  + s03 * CsHYZ3(2,2,1,03)
     $  + s04 * CsHYZ3(2,2,1,04)
     $  + s05 * CsHYZ3(2,2,1,05)
     $  + s06 * CsHYZ3(2,2,1,06)
     $  + s07 * CsHYZ3(2,2,1,07)
     $  + s08 * CsHYZ3(2,2,1,08)
     $  + s09 * CsHYZ3(2,2,1,09)
     $  + s10 * CsHYZ3(2,2,1,10)
     $  + s11 * CsHYZ3(2,2,1,11)
     $  + s12 * CsHYZ3(2,2,1,12)
     $  + s13 * CsHYZ3(2,2,1,13)
     $  + s14 * CsHYZ3(2,2,1,14)
     $  + s15 * CsHYZ3(2,2,1,15)
     $  + s16 * CsHYZ3(2,2,1,16)
     $  + s17 * CsHYZ3(2,2,1,17)
     $  + s18 * CsHYZ3(2,2,1,18)
     $  + s19 * CsHYZ3(2,2,1,19)
     $  + s20 * CsHYZ3(2,2,1,20)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,0,2,1) =
     $  + s02 * CsHYZ4(0,0,2,1,02)
     $  + s03 * CsHYZ4(0,0,2,1,03)
     $  + s04 * CsHYZ4(0,0,2,1,04)
     $  + s05 * CsHYZ4(0,0,2,1,05)
     $  + s06 * CsHYZ4(0,0,2,1,06)
     $  + s07 * CsHYZ4(0,0,2,1,07)
     $  + s08 * CsHYZ4(0,0,2,1,08)
     $  + s09 * CsHYZ4(0,0,2,1,09)
     $  + s10 * CsHYZ4(0,0,2,1,10)
     $  + s11 * CsHYZ4(0,0,2,1,11)
     $  + s12 * CsHYZ4(0,0,2,1,12)
     $  + s13 * CsHYZ4(0,0,2,1,13)
     $  + s14 * CsHYZ4(0,0,2,1,14)
     $  + s15 * CsHYZ4(0,0,2,1,15)
     $  + s16 * CsHYZ4(0,0,2,1,16)
     $  + s17 * CsHYZ4(0,0,2,1,17)
     $  + s18 * CsHYZ4(0,0,2,1,18)
      HYZ4(0,1,2,1) =
     $  + s03 * CsHYZ4(0,1,2,1,03)
     $  + s04 * CsHYZ4(0,1,2,1,04)
     $  + s05 * CsHYZ4(0,1,2,1,05)
     $  + s06 * CsHYZ4(0,1,2,1,06)
     $  + s07 * CsHYZ4(0,1,2,1,07)
     $  + s08 * CsHYZ4(0,1,2,1,08)
     $  + s09 * CsHYZ4(0,1,2,1,09)
     $  + s10 * CsHYZ4(0,1,2,1,10)
     $  + s11 * CsHYZ4(0,1,2,1,11)
     $  + s12 * CsHYZ4(0,1,2,1,12)
     $  + s13 * CsHYZ4(0,1,2,1,13)
     $  + s14 * CsHYZ4(0,1,2,1,14)
     $  + s15 * CsHYZ4(0,1,2,1,15)
     $  + s16 * CsHYZ4(0,1,2,1,16)
     $  + s17 * CsHYZ4(0,1,2,1,17)
     $  + s18 * CsHYZ4(0,1,2,1,18)
     $  + s19 * CsHYZ4(0,1,2,1,19)
     $  + s20 * CsHYZ4(0,1,2,1,20)
     $  + s21 * CsHYZ4(0,1,2,1,21)
      HYZ4(0,2,1,1) =
     $  + s03 * CsHYZ4(0,2,1,1,03)
     $  + s04 * CsHYZ4(0,2,1,1,04)
     $  + s05 * CsHYZ4(0,2,1,1,05)
     $  + s06 * CsHYZ4(0,2,1,1,06)
     $  + s07 * CsHYZ4(0,2,1,1,07)
     $  + s08 * CsHYZ4(0,2,1,1,08)
     $  + s09 * CsHYZ4(0,2,1,1,09)
     $  + s10 * CsHYZ4(0,2,1,1,10)
     $  + s11 * CsHYZ4(0,2,1,1,11)
     $  + s12 * CsHYZ4(0,2,1,1,12)
     $  + s13 * CsHYZ4(0,2,1,1,13)
     $  + s14 * CsHYZ4(0,2,1,1,14)
     $  + s15 * CsHYZ4(0,2,1,1,15)
     $  + s16 * CsHYZ4(0,2,1,1,16)
     $  + s17 * CsHYZ4(0,2,1,1,17)
     $  + s18 * CsHYZ4(0,2,1,1,18)
     $  + s19 * CsHYZ4(0,2,1,1,19)
     $  + s20 * CsHYZ4(0,2,1,1,20)
      HYZ4(0,2,0,1) =
     $  + s02 * CsHYZ4(0,2,0,1,02)
     $  + s03 * CsHYZ4(0,2,0,1,03)
     $  + s04 * CsHYZ4(0,2,0,1,04)
     $  + s05 * CsHYZ4(0,2,0,1,05)
     $  + s06 * CsHYZ4(0,2,0,1,06)
     $  + s07 * CsHYZ4(0,2,0,1,07)
     $  + s08 * CsHYZ4(0,2,0,1,08)
     $  + s09 * CsHYZ4(0,2,0,1,09)
     $  + s10 * CsHYZ4(0,2,0,1,10)
     $  + s11 * CsHYZ4(0,2,0,1,11)
     $  + s12 * CsHYZ4(0,2,0,1,12)
     $  + s13 * CsHYZ4(0,2,0,1,13)
     $  + s14 * CsHYZ4(0,2,0,1,14)
     $  + s15 * CsHYZ4(0,2,0,1,15)
     $  + s16 * CsHYZ4(0,2,0,1,16)
     $  + s17 * CsHYZ4(0,2,0,1,17)
     $  + s18 * CsHYZ4(0,2,0,1,18)
      HYZ4(0,2,2,1) =
     $  + s03 * CsHYZ4(0,2,2,1,03)
     $  + s04 * CsHYZ4(0,2,2,1,04)
     $  + s05 * CsHYZ4(0,2,2,1,05)
     $  + s06 * CsHYZ4(0,2,2,1,06)
     $  + s07 * CsHYZ4(0,2,2,1,07)
     $  + s08 * CsHYZ4(0,2,2,1,08)
     $  + s09 * CsHYZ4(0,2,2,1,09)
     $  + s10 * CsHYZ4(0,2,2,1,10)
     $  + s11 * CsHYZ4(0,2,2,1,11)
     $  + s12 * CsHYZ4(0,2,2,1,12)
     $  + s13 * CsHYZ4(0,2,2,1,13)
     $  + s14 * CsHYZ4(0,2,2,1,14)
     $  + s15 * CsHYZ4(0,2,2,1,15)
     $  + s16 * CsHYZ4(0,2,2,1,16)
     $  + s17 * CsHYZ4(0,2,2,1,17)
     $  + s18 * CsHYZ4(0,2,2,1,18)
     $  + s19 * CsHYZ4(0,2,2,1,19)
      HYZ4(2,0,0,1) =
     $  + s02 * CsHYZ4(2,0,0,1,02)
     $  + s03 * CsHYZ4(2,0,0,1,03)
     $  + s04 * CsHYZ4(2,0,0,1,04)
     $  + s05 * CsHYZ4(2,0,0,1,05)
     $  + s06 * CsHYZ4(2,0,0,1,06)
     $  + s07 * CsHYZ4(2,0,0,1,07)
     $  + s08 * CsHYZ4(2,0,0,1,08)
     $  + s09 * CsHYZ4(2,0,0,1,09)
     $  + s10 * CsHYZ4(2,0,0,1,10)
     $  + s11 * CsHYZ4(2,0,0,1,11)
     $  + s12 * CsHYZ4(2,0,0,1,12)
     $  + s13 * CsHYZ4(2,0,0,1,13)
     $  + s14 * CsHYZ4(2,0,0,1,14)
     $  + s15 * CsHYZ4(2,0,0,1,15)
     $  + s16 * CsHYZ4(2,0,0,1,16)
     $  + s17 * CsHYZ4(2,0,0,1,17)
     $  + s18 * CsHYZ4(2,0,0,1,18)
      HYZ4(2,0,1,1) =
     $  + s03 * CsHYZ4(2,0,1,1,03)
     $  + s04 * CsHYZ4(2,0,1,1,04)
     $  + s05 * CsHYZ4(2,0,1,1,05)
     $  + s06 * CsHYZ4(2,0,1,1,06)
     $  + s07 * CsHYZ4(2,0,1,1,07)
     $  + s08 * CsHYZ4(2,0,1,1,08)
     $  + s09 * CsHYZ4(2,0,1,1,09)
     $  + s10 * CsHYZ4(2,0,1,1,10)
     $  + s11 * CsHYZ4(2,0,1,1,11)
     $  + s12 * CsHYZ4(2,0,1,1,12)
     $  + s13 * CsHYZ4(2,0,1,1,13)
     $  + s14 * CsHYZ4(2,0,1,1,14)
     $  + s15 * CsHYZ4(2,0,1,1,15)
     $  + s16 * CsHYZ4(2,0,1,1,16)
     $  + s17 * CsHYZ4(2,0,1,1,17)
     $  + s18 * CsHYZ4(2,0,1,1,18)
     $  + s19 * CsHYZ4(2,0,1,1,19)
     $  + s20 * CsHYZ4(2,0,1,1,20)
      HYZ4(2,0,2,1) =
     $  + s03 * CsHYZ4(2,0,2,1,03)
     $  + s04 * CsHYZ4(2,0,2,1,04)
     $  + s05 * CsHYZ4(2,0,2,1,05)
     $  + s06 * CsHYZ4(2,0,2,1,06)
     $  + s07 * CsHYZ4(2,0,2,1,07)
     $  + s08 * CsHYZ4(2,0,2,1,08)
     $  + s09 * CsHYZ4(2,0,2,1,09)
     $  + s10 * CsHYZ4(2,0,2,1,10)
     $  + s11 * CsHYZ4(2,0,2,1,11)
     $  + s12 * CsHYZ4(2,0,2,1,12)
     $  + s13 * CsHYZ4(2,0,2,1,13)
     $  + s14 * CsHYZ4(2,0,2,1,14)
     $  + s15 * CsHYZ4(2,0,2,1,15)
     $  + s16 * CsHYZ4(2,0,2,1,16)
     $  + s17 * CsHYZ4(2,0,2,1,17)
     $  + s18 * CsHYZ4(2,0,2,1,18)
     $  + s19 * CsHYZ4(2,0,2,1,19)
      HYZ4(2,1,1,1) =
     $  + s04 * CsHYZ4(2,1,1,1,04)
     $  + s05 * CsHYZ4(2,1,1,1,05)
     $  + s06 * CsHYZ4(2,1,1,1,06)
     $  + s07 * CsHYZ4(2,1,1,1,07)
     $  + s08 * CsHYZ4(2,1,1,1,08)
     $  + s09 * CsHYZ4(2,1,1,1,09)
     $  + s10 * CsHYZ4(2,1,1,1,10)
     $  + s11 * CsHYZ4(2,1,1,1,11)
     $  + s12 * CsHYZ4(2,1,1,1,12)
     $  + s13 * CsHYZ4(2,1,1,1,13)
     $  + s14 * CsHYZ4(2,1,1,1,14)
     $  + s15 * CsHYZ4(2,1,1,1,15)
     $  + s16 * CsHYZ4(2,1,1,1,16)
     $  + s17 * CsHYZ4(2,1,1,1,17)
     $  + s18 * CsHYZ4(2,1,1,1,18)
     $  + s19 * CsHYZ4(2,1,1,1,19)
     $  + s20 * CsHYZ4(2,1,1,1,20)
     $  + s21 * CsHYZ4(2,1,1,1,21)
     $  + s22 * CsHYZ4(2,1,1,1,22)
      HYZ4(2,2,0,1) =
     $  + s03 * CsHYZ4(2,2,0,1,03)
     $  + s04 * CsHYZ4(2,2,0,1,04)
     $  + s05 * CsHYZ4(2,2,0,1,05)
     $  + s06 * CsHYZ4(2,2,0,1,06)
     $  + s07 * CsHYZ4(2,2,0,1,07)
     $  + s08 * CsHYZ4(2,2,0,1,08)
     $  + s09 * CsHYZ4(2,2,0,1,09)
     $  + s10 * CsHYZ4(2,2,0,1,10)
     $  + s11 * CsHYZ4(2,2,0,1,11)
     $  + s12 * CsHYZ4(2,2,0,1,12)
     $  + s13 * CsHYZ4(2,2,0,1,13)
     $  + s14 * CsHYZ4(2,2,0,1,14)
     $  + s15 * CsHYZ4(2,2,0,1,15)
     $  + s16 * CsHYZ4(2,2,0,1,16)
     $  + s17 * CsHYZ4(2,2,0,1,17)
     $  + s18 * CsHYZ4(2,2,0,1,18)
     $  + s19 * CsHYZ4(2,2,0,1,19)
      HYZ4(2,2,1,1) =
     $  + s04 * CsHYZ4(2,2,1,1,04)
     $  + s05 * CsHYZ4(2,2,1,1,05)
     $  + s06 * CsHYZ4(2,2,1,1,06)
     $  + s07 * CsHYZ4(2,2,1,1,07)
     $  + s08 * CsHYZ4(2,2,1,1,08)
     $  + s09 * CsHYZ4(2,2,1,1,09)
     $  + s10 * CsHYZ4(2,2,1,1,10)
     $  + s11 * CsHYZ4(2,2,1,1,11)
     $  + s12 * CsHYZ4(2,2,1,1,12)
     $  + s13 * CsHYZ4(2,2,1,1,13)
     $  + s14 * CsHYZ4(2,2,1,1,14)
     $  + s15 * CsHYZ4(2,2,1,1,15)
     $  + s16 * CsHYZ4(2,2,1,1,16)
     $  + s17 * CsHYZ4(2,2,1,1,17)
     $  + s18 * CsHYZ4(2,2,1,1,18)
     $  + s19 * CsHYZ4(2,2,1,1,19)
     $  + s20 * CsHYZ4(2,2,1,1,20)
     $  + s21 * CsHYZ4(2,2,1,1,21)
      HYZ4(2,2,2,1) =
     $  + s04 * CsHYZ4(2,2,2,1,04)
     $  + s05 * CsHYZ4(2,2,2,1,05)
     $  + s06 * CsHYZ4(2,2,2,1,06)
     $  + s07 * CsHYZ4(2,2,2,1,07)
     $  + s08 * CsHYZ4(2,2,2,1,08)
     $  + s09 * CsHYZ4(2,2,2,1,09)
     $  + s10 * CsHYZ4(2,2,2,1,10)
     $  + s11 * CsHYZ4(2,2,2,1,11)
     $  + s12 * CsHYZ4(2,2,2,1,12)
     $  + s13 * CsHYZ4(2,2,2,1,13)
     $  + s14 * CsHYZ4(2,2,2,1,14)
     $  + s15 * CsHYZ4(2,2,2,1,15)
     $  + s16 * CsHYZ4(2,2,2,1,16)
     $  + s17 * CsHYZ4(2,2,2,1,17)
     $  + s18 * CsHYZ4(2,2,2,1,18)
     $  + s19 * CsHYZ4(2,2,2,1,19)
      endif

      return
      end

      subroutine fillirr2dhpl310(iflag,n)
*********************************************************************
***** fillirr2dhpl310(iflag,n) fills the irreducible 2dHPL      *****
***** with indices (3,1,0)                                      *****
***** up to weight n using the EXPANDED expressions for the     *****
***** z-dependent expansion coefficients                        *****
***** applicable for all z                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      integer nmax
      real(dp):: CrHYZ2,CrHYZ3,CrHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,
     $     r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22
      dimension CrHYZ2(0:3,0:3,18),CrHYZ3(0:3,0:3,0:3,18),
     $          CrHYZ4(0:3,0:3,0:3,0:3,18)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common /rtdhpl/
     $     r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,
     $     r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22
      save CrHYZ2,CrHYZ3,CrHYZ4
!$omp threadprivate(CrHYZ2,CrHYZ3,CrHYZ4,/rtdhpl/,/HPL2/)

      if (iflag.eq.1) then
      call fillcoeff2dhpl310(iflag,n,CrHYZ2,CrHYZ3,CrHYZ4)
      endif
* 2001-04-03:17:41:05.hpl
* <- HPLexpand.310t
* produced by form-to-fortr for gehrt@pcth62

      if (n.eq.2) then
      HYZ2(3,1) =
     $  + r01 * CrHYZ2(3,1,01)
     $  + r02 * CrHYZ2(3,1,02)
     $  + r03 * CrHYZ2(3,1,03)
     $  + r04 * CrHYZ2(3,1,04)
     $  + r05 * CrHYZ2(3,1,05)
     $  + r06 * CrHYZ2(3,1,06)
     $  + r07 * CrHYZ2(3,1,07)
     $  + r08 * CrHYZ2(3,1,08)
     $  + r09 * CrHYZ2(3,1,09)
     $  + r10 * CrHYZ2(3,1,10)
     $  + r11 * CrHYZ2(3,1,11)
     $  + r12 * CrHYZ2(3,1,12)
     $  + r13 * CrHYZ2(3,1,13)
     $  + r14 * CrHYZ2(3,1,14)
     $  + r15 * CrHYZ2(3,1,15)
     $  + r16 * CrHYZ2(3,1,16)
     $  + r17 * CrHYZ2(3,1,17)
     $  - HZ1( -1)*HYZ1(3)
      endif
      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(0,3,1) =
     $  + r01 * CrHYZ3(0,3,1,01)
     $  + r02 * CrHYZ3(0,3,1,02)
     $  + r03 * CrHYZ3(0,3,1,03)
     $  + r04 * CrHYZ3(0,3,1,04)
     $  + r05 * CrHYZ3(0,3,1,05)
     $  + r06 * CrHYZ3(0,3,1,06)
     $  + r07 * CrHYZ3(0,3,1,07)
     $  + r08 * CrHYZ3(0,3,1,08)
     $  + r09 * CrHYZ3(0,3,1,09)
     $  + r10 * CrHYZ3(0,3,1,10)
     $  + r11 * CrHYZ3(0,3,1,11)
     $  + r12 * CrHYZ3(0,3,1,12)
     $  + r13 * CrHYZ3(0,3,1,13)
     $  + r14 * CrHYZ3(0,3,1,14)
     $  + r15 * CrHYZ3(0,3,1,15)
     $  + r16 * CrHYZ3(0,3,1,16)
     $  + r17 * CrHYZ3(0,3,1,17)
     $  - HZ1( -1)*HYZ2(0,3)
      HYZ3(3,0,1) =
     $  + r01 * CrHYZ3(3,0,1,01)
     $  + r02 * CrHYZ3(3,0,1,02)
     $  + r03 * CrHYZ3(3,0,1,03)
     $  + r04 * CrHYZ3(3,0,1,04)
     $  + r05 * CrHYZ3(3,0,1,05)
     $  + r06 * CrHYZ3(3,0,1,06)
     $  + r07 * CrHYZ3(3,0,1,07)
     $  + r08 * CrHYZ3(3,0,1,08)
     $  + r09 * CrHYZ3(3,0,1,09)
     $  + r10 * CrHYZ3(3,0,1,10)
     $  + r11 * CrHYZ3(3,0,1,11)
     $  + r12 * CrHYZ3(3,0,1,12)
     $  + r13 * CrHYZ3(3,0,1,13)
     $  + r14 * CrHYZ3(3,0,1,14)
     $  + r15 * CrHYZ3(3,0,1,15)
     $  + r16 * CrHYZ3(3,0,1,16)
     $  + r17 * CrHYZ3(3,0,1,17)
     $  - HZ2(0, -1)*HYZ1(3)
      HYZ3(3,1,1) =
     $  + r01 * CrHYZ3(3,1,1,01)
     $  + r02 * CrHYZ3(3,1,1,02)
     $  + r03 * CrHYZ3(3,1,1,03)
     $  + r04 * CrHYZ3(3,1,1,04)
     $  + r05 * CrHYZ3(3,1,1,05)
     $  + r06 * CrHYZ3(3,1,1,06)
     $  + r07 * CrHYZ3(3,1,1,07)
     $  + r08 * CrHYZ3(3,1,1,08)
     $  + r09 * CrHYZ3(3,1,1,09)
     $  + r10 * CrHYZ3(3,1,1,10)
     $  + r11 * CrHYZ3(3,1,1,11)
     $  + r12 * CrHYZ3(3,1,1,12)
     $  + r13 * CrHYZ3(3,1,1,13)
     $  + r14 * CrHYZ3(3,1,1,14)
     $  + r15 * CrHYZ3(3,1,1,15)
     $  + r16 * CrHYZ3(3,1,1,16)
     $  + r17 * CrHYZ3(3,1,1,17)
     $  + HZ2( -1,-1)*HYZ1(3)
      HYZ3(3,3,1) =
     $  + r01 * CrHYZ3(3,3,1,01)
     $  + r02 * CrHYZ3(3,3,1,02)
     $  + r03 * CrHYZ3(3,3,1,03)
     $  + r04 * CrHYZ3(3,3,1,04)
     $  + r05 * CrHYZ3(3,3,1,05)
     $  + r06 * CrHYZ3(3,3,1,06)
     $  + r07 * CrHYZ3(3,3,1,07)
     $  + r08 * CrHYZ3(3,3,1,08)
     $  + r09 * CrHYZ3(3,3,1,09)
     $  + r10 * CrHYZ3(3,3,1,10)
     $  + r11 * CrHYZ3(3,3,1,11)
     $  + r12 * CrHYZ3(3,3,1,12)
     $  + r13 * CrHYZ3(3,3,1,13)
     $  + r14 * CrHYZ3(3,3,1,14)
     $  + r15 * CrHYZ3(3,3,1,15)
     $  + r16 * CrHYZ3(3,3,1,16)
     $  + r17 * CrHYZ3(3,3,1,17)
     $  - HZ1( -1)*HYZ2(3,3)
     $  + HZ2( -1,-1)*HYZ1(3)
     $  - HZ2(0, -1)*HYZ1(3)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,0,3,1) =
     $  + r01 * CrHYZ4(0,0,3,1,01)
     $  + r02 * CrHYZ4(0,0,3,1,02)
     $  + r03 * CrHYZ4(0,0,3,1,03)
     $  + r04 * CrHYZ4(0,0,3,1,04)
     $  + r05 * CrHYZ4(0,0,3,1,05)
     $  + r06 * CrHYZ4(0,0,3,1,06)
     $  + r07 * CrHYZ4(0,0,3,1,07)
     $  + r08 * CrHYZ4(0,0,3,1,08)
     $  + r09 * CrHYZ4(0,0,3,1,09)
     $  + r10 * CrHYZ4(0,0,3,1,10)
     $  + r11 * CrHYZ4(0,0,3,1,11)
     $  + r12 * CrHYZ4(0,0,3,1,12)
     $  + r13 * CrHYZ4(0,0,3,1,13)
     $  + r14 * CrHYZ4(0,0,3,1,14)
     $  + r15 * CrHYZ4(0,0,3,1,15)
     $  + r16 * CrHYZ4(0,0,3,1,16)
     $  + r17 * CrHYZ4(0,0,3,1,17)
     $  - HZ1( -1)*HYZ3(0,0,3)
      HYZ4(0,1,3,1) =
     $  + r01 * CrHYZ4(0,1,3,1,01)
     $  + r02 * CrHYZ4(0,1,3,1,02)
     $  + r03 * CrHYZ4(0,1,3,1,03)
     $  + r04 * CrHYZ4(0,1,3,1,04)
     $  + r05 * CrHYZ4(0,1,3,1,05)
     $  + r06 * CrHYZ4(0,1,3,1,06)
     $  + r07 * CrHYZ4(0,1,3,1,07)
     $  + r08 * CrHYZ4(0,1,3,1,08)
     $  + r09 * CrHYZ4(0,1,3,1,09)
     $  + r10 * CrHYZ4(0,1,3,1,10)
     $  + r11 * CrHYZ4(0,1,3,1,11)
     $  + r12 * CrHYZ4(0,1,3,1,12)
     $  + r13 * CrHYZ4(0,1,3,1,13)
     $  + r14 * CrHYZ4(0,1,3,1,14)
     $  + r15 * CrHYZ4(0,1,3,1,15)
     $  + r16 * CrHYZ4(0,1,3,1,16)
     $  + r17 * CrHYZ4(0,1,3,1,17)
     $  - HZ1( -1)*HYZ3(0,1,3)
     $  + 2.000000000000000d+00*HZ2(-1,-1)*HYZ2(0,1)
      HYZ4(0,3,0,1) =
     $  + r01 * CrHYZ4(0,3,0,1,01)
     $  + r02 * CrHYZ4(0,3,0,1,02)
     $  + r03 * CrHYZ4(0,3,0,1,03)
     $  + r04 * CrHYZ4(0,3,0,1,04)
     $  + r05 * CrHYZ4(0,3,0,1,05)
     $  + r06 * CrHYZ4(0,3,0,1,06)
     $  + r07 * CrHYZ4(0,3,0,1,07)
     $  + r08 * CrHYZ4(0,3,0,1,08)
     $  + r09 * CrHYZ4(0,3,0,1,09)
     $  + r10 * CrHYZ4(0,3,0,1,10)
     $  + r11 * CrHYZ4(0,3,0,1,11)
     $  + r12 * CrHYZ4(0,3,0,1,12)
     $  + r13 * CrHYZ4(0,3,0,1,13)
     $  + r14 * CrHYZ4(0,3,0,1,14)
     $  + r15 * CrHYZ4(0,3,0,1,15)
     $  + r16 * CrHYZ4(0,3,0,1,16)
     $  + r17 * CrHYZ4(0,3,0,1,17)
     $  - HZ2(0, -1)*HYZ2(0,3)
      HYZ4(0,3,1,1) =
     $  + r01 * CrHYZ4(0,3,1,1,01)
     $  + r02 * CrHYZ4(0,3,1,1,02)
     $  + r03 * CrHYZ4(0,3,1,1,03)
     $  + r04 * CrHYZ4(0,3,1,1,04)
     $  + r05 * CrHYZ4(0,3,1,1,05)
     $  + r06 * CrHYZ4(0,3,1,1,06)
     $  + r07 * CrHYZ4(0,3,1,1,07)
     $  + r08 * CrHYZ4(0,3,1,1,08)
     $  + r09 * CrHYZ4(0,3,1,1,09)
     $  + r10 * CrHYZ4(0,3,1,1,10)
     $  + r11 * CrHYZ4(0,3,1,1,11)
     $  + r12 * CrHYZ4(0,3,1,1,12)
     $  + r13 * CrHYZ4(0,3,1,1,13)
     $  + r14 * CrHYZ4(0,3,1,1,14)
     $  + r15 * CrHYZ4(0,3,1,1,15)
     $  + r16 * CrHYZ4(0,3,1,1,16)
     $  + r17 * CrHYZ4(0,3,1,1,17)
     $  + HZ2( -1,-1)*HYZ2(0,3)
      HYZ4(0,3,3,1) =
     $  + r01 * CrHYZ4(0,3,3,1,01)
     $  + r02 * CrHYZ4(0,3,3,1,02)
     $  + r03 * CrHYZ4(0,3,3,1,03)
     $  + r04 * CrHYZ4(0,3,3,1,04)
     $  + r05 * CrHYZ4(0,3,3,1,05)
     $  + r06 * CrHYZ4(0,3,3,1,06)
     $  + r07 * CrHYZ4(0,3,3,1,07)
     $  + r08 * CrHYZ4(0,3,3,1,08)
     $  + r09 * CrHYZ4(0,3,3,1,09)
     $  + r10 * CrHYZ4(0,3,3,1,10)
     $  + r11 * CrHYZ4(0,3,3,1,11)
     $  + r12 * CrHYZ4(0,3,3,1,12)
     $  + r13 * CrHYZ4(0,3,3,1,13)
     $  + r14 * CrHYZ4(0,3,3,1,14)
     $  + r15 * CrHYZ4(0,3,3,1,15)
     $  + r16 * CrHYZ4(0,3,3,1,16)
     $  + r17 * CrHYZ4(0,3,3,1,17)
     $  - HZ1( -1)*HYZ3(0,3,3)
     $  + HZ2( -1,-1)*HYZ2(0,3)
     $  - HZ2(0, -1)*HYZ2(0,3)
      HYZ4(3,0,0,1) =
     $  + r01 * CrHYZ4(3,0,0,1,01)
     $  + r02 * CrHYZ4(3,0,0,1,02)
     $  + r03 * CrHYZ4(3,0,0,1,03)
     $  + r04 * CrHYZ4(3,0,0,1,04)
     $  + r05 * CrHYZ4(3,0,0,1,05)
     $  + r06 * CrHYZ4(3,0,0,1,06)
     $  + r07 * CrHYZ4(3,0,0,1,07)
     $  + r08 * CrHYZ4(3,0,0,1,08)
     $  + r09 * CrHYZ4(3,0,0,1,09)
     $  + r10 * CrHYZ4(3,0,0,1,10)
     $  + r11 * CrHYZ4(3,0,0,1,11)
     $  + r12 * CrHYZ4(3,0,0,1,12)
     $  + r13 * CrHYZ4(3,0,0,1,13)
     $  + r14 * CrHYZ4(3,0,0,1,14)
     $  + r15 * CrHYZ4(3,0,0,1,15)
     $  + r16 * CrHYZ4(3,0,0,1,16)
     $  + r17 * CrHYZ4(3,0,0,1,17)
     $  - HZ3(0,0, -1)*HYZ1(3)
      HYZ4(3,0,1,1) =
     $  + r01 * CrHYZ4(3,0,1,1,01)
     $  + r02 * CrHYZ4(3,0,1,1,02)
     $  + r03 * CrHYZ4(3,0,1,1,03)
     $  + r04 * CrHYZ4(3,0,1,1,04)
     $  + r05 * CrHYZ4(3,0,1,1,05)
     $  + r06 * CrHYZ4(3,0,1,1,06)
     $  + r07 * CrHYZ4(3,0,1,1,07)
     $  + r08 * CrHYZ4(3,0,1,1,08)
     $  + r09 * CrHYZ4(3,0,1,1,09)
     $  + r10 * CrHYZ4(3,0,1,1,10)
     $  + r11 * CrHYZ4(3,0,1,1,11)
     $  + r12 * CrHYZ4(3,0,1,1,12)
     $  + r13 * CrHYZ4(3,0,1,1,13)
     $  + r14 * CrHYZ4(3,0,1,1,14)
     $  + r15 * CrHYZ4(3,0,1,1,15)
     $  + r16 * CrHYZ4(3,0,1,1,16)
     $  + r17 * CrHYZ4(3,0,1,1,17)
     $  + r18 * CrHYZ4(3,0,1,1,18)
     $  + HZ3(0, -1,-1)*HYZ1(3)
      HYZ4(3,0,3,1) =
     $  + r01 * CrHYZ4(3,0,3,1,01)
     $  + r02 * CrHYZ4(3,0,3,1,02)
     $  + r03 * CrHYZ4(3,0,3,1,03)
     $  + r04 * CrHYZ4(3,0,3,1,04)
     $  + r05 * CrHYZ4(3,0,3,1,05)
     $  + r06 * CrHYZ4(3,0,3,1,06)
     $  + r07 * CrHYZ4(3,0,3,1,07)
     $  + r08 * CrHYZ4(3,0,3,1,08)
     $  + r09 * CrHYZ4(3,0,3,1,09)
     $  + r10 * CrHYZ4(3,0,3,1,10)
     $  + r11 * CrHYZ4(3,0,3,1,11)
     $  + r12 * CrHYZ4(3,0,3,1,12)
     $  + r13 * CrHYZ4(3,0,3,1,13)
     $  + r14 * CrHYZ4(3,0,3,1,14)
     $  + r15 * CrHYZ4(3,0,3,1,15)
     $  + r16 * CrHYZ4(3,0,3,1,16)
     $  + r17 * CrHYZ4(3,0,3,1,17)
     $  - HZ1( -1)*HYZ3(3,0,3)
     $  + HZ3( -1,0,-1)*HYZ1(3)
     $  - HZ3(0,0, -1)*HYZ1(3)
      HYZ4(3,1,1,1) =
     $  + r01 * CrHYZ4(3,1,1,1,01)
     $  + r02 * CrHYZ4(3,1,1,1,02)
     $  + r03 * CrHYZ4(3,1,1,1,03)
     $  + r04 * CrHYZ4(3,1,1,1,04)
     $  + r05 * CrHYZ4(3,1,1,1,05)
     $  + r06 * CrHYZ4(3,1,1,1,06)
     $  + r07 * CrHYZ4(3,1,1,1,07)
     $  + r08 * CrHYZ4(3,1,1,1,08)
     $  + r09 * CrHYZ4(3,1,1,1,09)
     $  + r10 * CrHYZ4(3,1,1,1,10)
     $  + r11 * CrHYZ4(3,1,1,1,11)
     $  + r12 * CrHYZ4(3,1,1,1,12)
     $  + r13 * CrHYZ4(3,1,1,1,13)
     $  + r14 * CrHYZ4(3,1,1,1,14)
     $  + r15 * CrHYZ4(3,1,1,1,15)
     $  + r16 * CrHYZ4(3,1,1,1,16)
     $  + r17 * CrHYZ4(3,1,1,1,17)
     $  - HZ3( -1,-1,-1)*HYZ1(3)
      HYZ4(3,3,0,1) =
     $  + r01 * CrHYZ4(3,3,0,1,01)
     $  + r02 * CrHYZ4(3,3,0,1,02)
     $  + r03 * CrHYZ4(3,3,0,1,03)
     $  + r04 * CrHYZ4(3,3,0,1,04)
     $  + r05 * CrHYZ4(3,3,0,1,05)
     $  + r06 * CrHYZ4(3,3,0,1,06)
     $  + r07 * CrHYZ4(3,3,0,1,07)
     $  + r08 * CrHYZ4(3,3,0,1,08)
     $  + r09 * CrHYZ4(3,3,0,1,09)
     $  + r10 * CrHYZ4(3,3,0,1,10)
     $  + r11 * CrHYZ4(3,3,0,1,11)
     $  + r12 * CrHYZ4(3,3,0,1,12)
     $  + r13 * CrHYZ4(3,3,0,1,13)
     $  + r14 * CrHYZ4(3,3,0,1,14)
     $  + r15 * CrHYZ4(3,3,0,1,15)
     $  + r16 * CrHYZ4(3,3,0,1,16)
     $  + r17 * CrHYZ4(3,3,0,1,17)
     $  - HZ2(0, -1)*HYZ2(3,3)
     $  + HZ3(0, -1,-1)*HYZ1(3)
     $  - HZ3(0,0, -1)*HYZ1(3)
      HYZ4(3,3,1,1) =
     $  + r01 * CrHYZ4(3,3,1,1,01)
     $  + r02 * CrHYZ4(3,3,1,1,02)
     $  + r03 * CrHYZ4(3,3,1,1,03)
     $  + r04 * CrHYZ4(3,3,1,1,04)
     $  + r05 * CrHYZ4(3,3,1,1,05)
     $  + r06 * CrHYZ4(3,3,1,1,06)
     $  + r07 * CrHYZ4(3,3,1,1,07)
     $  + r08 * CrHYZ4(3,3,1,1,08)
     $  + r09 * CrHYZ4(3,3,1,1,09)
     $  + r10 * CrHYZ4(3,3,1,1,10)
     $  + r11 * CrHYZ4(3,3,1,1,11)
     $  + r12 * CrHYZ4(3,3,1,1,12)
     $  + r13 * CrHYZ4(3,3,1,1,13)
     $  + r14 * CrHYZ4(3,3,1,1,14)
     $  + r15 * CrHYZ4(3,3,1,1,15)
     $  + r16 * CrHYZ4(3,3,1,1,16)
     $  + r17 * CrHYZ4(3,3,1,1,17)
     $  + HZ2( -1,-1)*HYZ2(3,3)
     $  - 2.000000000000000d+00*HZ3(-1,-1,-1)*HYZ1(3)
     $  + HZ3( -1,0,-1)*HYZ1(3)
     $  + HZ3(0, -1,-1)*HYZ1(3)
      HYZ4(3,3,3,1) =
     $  + r01 * CrHYZ4(3,3,3,1,01)
     $  + r02 * CrHYZ4(3,3,3,1,02)
     $  + r03 * CrHYZ4(3,3,3,1,03)
     $  + r04 * CrHYZ4(3,3,3,1,04)
     $  + r05 * CrHYZ4(3,3,3,1,05)
     $  + r06 * CrHYZ4(3,3,3,1,06)
     $  + r07 * CrHYZ4(3,3,3,1,07)
     $  + r08 * CrHYZ4(3,3,3,1,08)
     $  + r09 * CrHYZ4(3,3,3,1,09)
     $  + r10 * CrHYZ4(3,3,3,1,10)
     $  + r11 * CrHYZ4(3,3,3,1,11)
     $  + r12 * CrHYZ4(3,3,3,1,12)
     $  + r13 * CrHYZ4(3,3,3,1,13)
     $  + r14 * CrHYZ4(3,3,3,1,14)
     $  + r15 * CrHYZ4(3,3,3,1,15)
     $  + r16 * CrHYZ4(3,3,3,1,16)
     $  + r17 * CrHYZ4(3,3,3,1,17)
     $  - HZ1( -1)*HYZ3(3,3,3)
     $  + HZ2( -1,-1)*HYZ2(3,3)
     $  - HZ2(0, -1)*HYZ2(3,3)
     $  - HZ3( -1,-1,-1)*HYZ1(3)
     $  + HZ3( -1,0,-1)*HYZ1(3)
     $  + HZ3(0, -1,-1)*HYZ1(3)
     $  - HZ3(0,0, -1)*HYZ1(3)
      endif

      return
      end


      subroutine fillirr2dhpl320(iflag,n)
*********************************************************************
***** fillirr2dhpl320(iflag,n) fills the irreducible 2dHPL      *****
***** with indices (3,2,0)                                      *****
***** up to weight n using the EXPANDED expressions for the     *****
***** z-dependent expansion coefficients                        *****
***** applicable for z<0.5                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: CsHYZ2,CsHYZ3,CsHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      dimension CsHYZ2(0:3,0:3,22),CsHYZ3(0:3,0:3,0:3,22),
     $          CsHYZ4(0:3,0:3,0:3,0:3,22)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common /s/
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      save CsHYZ2,CsHYZ3,CsHYZ4
!$omp threadprivate(CsHYZ2,CsHYZ3,CsHYZ4,/s/,/HPL2/)

      if (iflag.eq.1) then
      call fillcoeff2dhpl320(iflag,n,CsHYZ2,CsHYZ3,CsHYZ4)
      endif

* 2001-04-20:10:54:20.hpl
* <- HPLexpand.320t
* produced by form-to-fortr for gehrt@pcth62

      if (n.eq.2) then
      HYZ2(3,2) =
     $  + s01 * CsHYZ2(3,2,01)
     $  + s02 * CsHYZ2(3,2,02)
     $  + s03 * CsHYZ2(3,2,03)
     $  + s04 * CsHYZ2(3,2,04)
     $  + s05 * CsHYZ2(3,2,05)
     $  + s06 * CsHYZ2(3,2,06)
     $  + s07 * CsHYZ2(3,2,07)
     $  + s08 * CsHYZ2(3,2,08)
     $  + s09 * CsHYZ2(3,2,09)
     $  + s10 * CsHYZ2(3,2,10)
     $  + s11 * CsHYZ2(3,2,11)
     $  + s12 * CsHYZ2(3,2,12)
     $  + s13 * CsHYZ2(3,2,13)
     $  + s14 * CsHYZ2(3,2,14)
     $  + s15 * CsHYZ2(3,2,15)
     $  + s16 * CsHYZ2(3,2,16)
     $  + s17 * CsHYZ2(3,2,17)
     $  - HZ1(1) *HYZ1(3)
      endif
      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(0,3,2) =
     $  + s01 * CsHYZ3(0,3,2,01)
     $  + s02 * CsHYZ3(0,3,2,02)
     $  + s03 * CsHYZ3(0,3,2,03)
     $  + s04 * CsHYZ3(0,3,2,04)
     $  + s05 * CsHYZ3(0,3,2,05)
     $  + s06 * CsHYZ3(0,3,2,06)
     $  + s07 * CsHYZ3(0,3,2,07)
     $  + s08 * CsHYZ3(0,3,2,08)
     $  + s09 * CsHYZ3(0,3,2,09)
     $  + s10 * CsHYZ3(0,3,2,10)
     $  + s11 * CsHYZ3(0,3,2,11)
     $  + s12 * CsHYZ3(0,3,2,12)
     $  + s13 * CsHYZ3(0,3,2,13)
     $  + s14 * CsHYZ3(0,3,2,14)
     $  + s15 * CsHYZ3(0,3,2,15)
     $  + s16 * CsHYZ3(0,3,2,16)
     $  + s17 * CsHYZ3(0,3,2,17)
     $  - HZ1(1) *HYZ2(0,3)
      HYZ3(3,0,2) =
     $  + s01 * CsHYZ3(3,0,2,01)
     $  + s02 * CsHYZ3(3,0,2,02)
     $  + s03 * CsHYZ3(3,0,2,03)
     $  + s04 * CsHYZ3(3,0,2,04)
     $  + s05 * CsHYZ3(3,0,2,05)
     $  + s06 * CsHYZ3(3,0,2,06)
     $  + s07 * CsHYZ3(3,0,2,07)
     $  + s08 * CsHYZ3(3,0,2,08)
     $  + s09 * CsHYZ3(3,0,2,09)
     $  + s10 * CsHYZ3(3,0,2,10)
     $  + s11 * CsHYZ3(3,0,2,11)
     $  + s12 * CsHYZ3(3,0,2,12)
     $  + s13 * CsHYZ3(3,0,2,13)
     $  + s14 * CsHYZ3(3,0,2,14)
     $  + s15 * CsHYZ3(3,0,2,15)
     $  + s16 * CsHYZ3(3,0,2,16)
     $  + s17 * CsHYZ3(3,0,2,17)
     $  - HZ2(0,1) *HYZ1(3)
     $  - HZ2(1,1) *HYZ1(3)
      HYZ3(3,2,2) =
     $  + s01 * CsHYZ3(3,2,2,01)
     $  + s02 * CsHYZ3(3,2,2,02)
     $  + s03 * CsHYZ3(3,2,2,03)
     $  + s04 * CsHYZ3(3,2,2,04)
     $  + s05 * CsHYZ3(3,2,2,05)
     $  + s06 * CsHYZ3(3,2,2,06)
     $  + s07 * CsHYZ3(3,2,2,07)
     $  + s08 * CsHYZ3(3,2,2,08)
     $  + s09 * CsHYZ3(3,2,2,09)
     $  + s10 * CsHYZ3(3,2,2,10)
     $  + s11 * CsHYZ3(3,2,2,11)
     $  + s12 * CsHYZ3(3,2,2,12)
     $  + s13 * CsHYZ3(3,2,2,13)
     $  + s14 * CsHYZ3(3,2,2,14)
     $  + s15 * CsHYZ3(3,2,2,15)
     $  + s16 * CsHYZ3(3,2,2,16)
     $  + s17 * CsHYZ3(3,2,2,17)
     $  + HZ2(1,1) *HYZ1(3)
      HYZ3(3,3,2) =
     $  + s01 * CsHYZ3(3,3,2,01)
     $  + s02 * CsHYZ3(3,3,2,02)
     $  + s03 * CsHYZ3(3,3,2,03)
     $  + s04 * CsHYZ3(3,3,2,04)
     $  + s05 * CsHYZ3(3,3,2,05)
     $  + s06 * CsHYZ3(3,3,2,06)
     $  + s07 * CsHYZ3(3,3,2,07)
     $  + s08 * CsHYZ3(3,3,2,08)
     $  + s09 * CsHYZ3(3,3,2,09)
     $  + s10 * CsHYZ3(3,3,2,10)
     $  + s11 * CsHYZ3(3,3,2,11)
     $  + s12 * CsHYZ3(3,3,2,12)
     $  + s13 * CsHYZ3(3,3,2,13)
     $  + s14 * CsHYZ3(3,3,2,14)
     $  + s15 * CsHYZ3(3,3,2,15)
     $  + s16 * CsHYZ3(3,3,2,16)
     $  + s17 * CsHYZ3(3,3,2,17)
     $  - HZ1(1) *HYZ2(3,3)
     $  - HZ2(0,1) *HYZ1(3)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,0,3,2) =
     $  + s01 * CsHYZ4(0,0,3,2,01)
     $  + s02 * CsHYZ4(0,0,3,2,02)
     $  + s03 * CsHYZ4(0,0,3,2,03)
     $  + s04 * CsHYZ4(0,0,3,2,04)
     $  + s05 * CsHYZ4(0,0,3,2,05)
     $  + s06 * CsHYZ4(0,0,3,2,06)
     $  + s07 * CsHYZ4(0,0,3,2,07)
     $  + s08 * CsHYZ4(0,0,3,2,08)
     $  + s09 * CsHYZ4(0,0,3,2,09)
     $  + s10 * CsHYZ4(0,0,3,2,10)
     $  + s11 * CsHYZ4(0,0,3,2,11)
     $  + s12 * CsHYZ4(0,0,3,2,12)
     $  + s13 * CsHYZ4(0,0,3,2,13)
     $  + s14 * CsHYZ4(0,0,3,2,14)
     $  + s15 * CsHYZ4(0,0,3,2,15)
     $  + s16 * CsHYZ4(0,0,3,2,16)
     $  + s17 * CsHYZ4(0,0,3,2,17)
     $  - HZ1(1) *HYZ3(0,0,3)
      HYZ4(0,2,3,2) =
     $  + s02 * CsHYZ4(0,2,3,2,02)
     $  + s03 * CsHYZ4(0,2,3,2,03)
     $  + s04 * CsHYZ4(0,2,3,2,04)
     $  + s05 * CsHYZ4(0,2,3,2,05)
     $  + s06 * CsHYZ4(0,2,3,2,06)
     $  + s07 * CsHYZ4(0,2,3,2,07)
     $  + s08 * CsHYZ4(0,2,3,2,08)
     $  + s09 * CsHYZ4(0,2,3,2,09)
     $  + s10 * CsHYZ4(0,2,3,2,10)
     $  + s11 * CsHYZ4(0,2,3,2,11)
     $  + s12 * CsHYZ4(0,2,3,2,12)
     $  + s13 * CsHYZ4(0,2,3,2,13)
     $  + s14 * CsHYZ4(0,2,3,2,14)
     $  + s15 * CsHYZ4(0,2,3,2,15)
     $  + s16 * CsHYZ4(0,2,3,2,16)
     $  + s17 * CsHYZ4(0,2,3,2,17)
     $  - HZ1(1) *HYZ3(0,2,3)
      HYZ4(0,3,0,2) =
     $  + s01 * CsHYZ4(0,3,0,2,01)
     $  + s02 * CsHYZ4(0,3,0,2,02)
     $  + s03 * CsHYZ4(0,3,0,2,03)
     $  + s04 * CsHYZ4(0,3,0,2,04)
     $  + s05 * CsHYZ4(0,3,0,2,05)
     $  + s06 * CsHYZ4(0,3,0,2,06)
     $  + s07 * CsHYZ4(0,3,0,2,07)
     $  + s08 * CsHYZ4(0,3,0,2,08)
     $  + s09 * CsHYZ4(0,3,0,2,09)
     $  + s10 * CsHYZ4(0,3,0,2,10)
     $  + s11 * CsHYZ4(0,3,0,2,11)
     $  + s12 * CsHYZ4(0,3,0,2,12)
     $  + s13 * CsHYZ4(0,3,0,2,13)
     $  + s14 * CsHYZ4(0,3,0,2,14)
     $  + s15 * CsHYZ4(0,3,0,2,15)
     $  + s16 * CsHYZ4(0,3,0,2,16)
     $  + s17 * CsHYZ4(0,3,0,2,17)
     $  - HZ2(0,1) *HYZ2(0,3)
     $  - HZ2(1,1) *HYZ2(0,3)
      HYZ4(0,3,2,2) =
     $  + s01 * CsHYZ4(0,3,2,2,01)
     $  + s02 * CsHYZ4(0,3,2,2,02)
     $  + s03 * CsHYZ4(0,3,2,2,03)
     $  + s04 * CsHYZ4(0,3,2,2,04)
     $  + s05 * CsHYZ4(0,3,2,2,05)
     $  + s06 * CsHYZ4(0,3,2,2,06)
     $  + s07 * CsHYZ4(0,3,2,2,07)
     $  + s08 * CsHYZ4(0,3,2,2,08)
     $  + s09 * CsHYZ4(0,3,2,2,09)
     $  + s10 * CsHYZ4(0,3,2,2,10)
     $  + s11 * CsHYZ4(0,3,2,2,11)
     $  + s12 * CsHYZ4(0,3,2,2,12)
     $  + s13 * CsHYZ4(0,3,2,2,13)
     $  + s14 * CsHYZ4(0,3,2,2,14)
     $  + s15 * CsHYZ4(0,3,2,2,15)
     $  + s16 * CsHYZ4(0,3,2,2,16)
     $  + s17 * CsHYZ4(0,3,2,2,17)
     $  + HZ2(1,1) *HYZ2(0,3)
      HYZ4(0,3,3,2) =
     $  + s01 * CsHYZ4(0,3,3,2,01)
     $  + s02 * CsHYZ4(0,3,3,2,02)
     $  + s03 * CsHYZ4(0,3,3,2,03)
     $  + s04 * CsHYZ4(0,3,3,2,04)
     $  + s05 * CsHYZ4(0,3,3,2,05)
     $  + s06 * CsHYZ4(0,3,3,2,06)
     $  + s07 * CsHYZ4(0,3,3,2,07)
     $  + s08 * CsHYZ4(0,3,3,2,08)
     $  + s09 * CsHYZ4(0,3,3,2,09)
     $  + s10 * CsHYZ4(0,3,3,2,10)
     $  + s11 * CsHYZ4(0,3,3,2,11)
     $  + s12 * CsHYZ4(0,3,3,2,12)
     $  + s13 * CsHYZ4(0,3,3,2,13)
     $  + s14 * CsHYZ4(0,3,3,2,14)
     $  + s15 * CsHYZ4(0,3,3,2,15)
     $  + s16 * CsHYZ4(0,3,3,2,16)
     $  + s17 * CsHYZ4(0,3,3,2,17)
     $  - HZ1(1) *HYZ3(0,3,3)
     $  - HZ2(0,1) *HYZ2(0,3)
      HYZ4(3,0,0,2) =
     $  + s01 * CsHYZ4(3,0,0,2,01)
     $  + s02 * CsHYZ4(3,0,0,2,02)
     $  + s03 * CsHYZ4(3,0,0,2,03)
     $  + s04 * CsHYZ4(3,0,0,2,04)
     $  + s05 * CsHYZ4(3,0,0,2,05)
     $  + s06 * CsHYZ4(3,0,0,2,06)
     $  + s07 * CsHYZ4(3,0,0,2,07)
     $  + s08 * CsHYZ4(3,0,0,2,08)
     $  + s09 * CsHYZ4(3,0,0,2,09)
     $  + s10 * CsHYZ4(3,0,0,2,10)
     $  + s11 * CsHYZ4(3,0,0,2,11)
     $  + s12 * CsHYZ4(3,0,0,2,12)
     $  + s13 * CsHYZ4(3,0,0,2,13)
     $  + s14 * CsHYZ4(3,0,0,2,14)
     $  + s15 * CsHYZ4(3,0,0,2,15)
     $  + s16 * CsHYZ4(3,0,0,2,16)
     $  + s17 * CsHYZ4(3,0,0,2,17)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(0,1,1) *HYZ1(3)
     $  - HZ3(1,0,1) *HYZ1(3)
     $  - HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,0,2,2) =
     $  + s01 * CsHYZ4(3,0,2,2,01)
     $  + s02 * CsHYZ4(3,0,2,2,02)
     $  + s03 * CsHYZ4(3,0,2,2,03)
     $  + s04 * CsHYZ4(3,0,2,2,04)
     $  + s05 * CsHYZ4(3,0,2,2,05)
     $  + s06 * CsHYZ4(3,0,2,2,06)
     $  + s07 * CsHYZ4(3,0,2,2,07)
     $  + s08 * CsHYZ4(3,0,2,2,08)
     $  + s09 * CsHYZ4(3,0,2,2,09)
     $  + s10 * CsHYZ4(3,0,2,2,10)
     $  + s11 * CsHYZ4(3,0,2,2,11)
     $  + s12 * CsHYZ4(3,0,2,2,12)
     $  + s13 * CsHYZ4(3,0,2,2,13)
     $  + s14 * CsHYZ4(3,0,2,2,14)
     $  + s15 * CsHYZ4(3,0,2,2,15)
     $  + s16 * CsHYZ4(3,0,2,2,16)
     $  + s17 * CsHYZ4(3,0,2,2,17)
     $  + s18 * CsHYZ4(3,0,2,2,18)
     $  + HZ3(0,1,1) *HYZ1(3)
     $  + HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,0,3,2) =
     $  + s01 * CsHYZ4(3,0,3,2,01)
     $  + s02 * CsHYZ4(3,0,3,2,02)
     $  + s03 * CsHYZ4(3,0,3,2,03)
     $  + s04 * CsHYZ4(3,0,3,2,04)
     $  + s05 * CsHYZ4(3,0,3,2,05)
     $  + s06 * CsHYZ4(3,0,3,2,06)
     $  + s07 * CsHYZ4(3,0,3,2,07)
     $  + s08 * CsHYZ4(3,0,3,2,08)
     $  + s09 * CsHYZ4(3,0,3,2,09)
     $  + s10 * CsHYZ4(3,0,3,2,10)
     $  + s11 * CsHYZ4(3,0,3,2,11)
     $  + s12 * CsHYZ4(3,0,3,2,12)
     $  + s13 * CsHYZ4(3,0,3,2,13)
     $  + s14 * CsHYZ4(3,0,3,2,14)
     $  + s15 * CsHYZ4(3,0,3,2,15)
     $  + s16 * CsHYZ4(3,0,3,2,16)
     $  + s17 * CsHYZ4(3,0,3,2,17)
     $  - HZ1(1) *HYZ3(3,0,3)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(0,1,1) *HYZ1(3)
      HYZ4(3,2,2,2) =
     $  + s01 * CsHYZ4(3,2,2,2,01)
     $  + s02 * CsHYZ4(3,2,2,2,02)
     $  + s03 * CsHYZ4(3,2,2,2,03)
     $  + s04 * CsHYZ4(3,2,2,2,04)
     $  + s05 * CsHYZ4(3,2,2,2,05)
     $  + s06 * CsHYZ4(3,2,2,2,06)
     $  + s07 * CsHYZ4(3,2,2,2,07)
     $  + s08 * CsHYZ4(3,2,2,2,08)
     $  + s09 * CsHYZ4(3,2,2,2,09)
     $  + s10 * CsHYZ4(3,2,2,2,10)
     $  + s11 * CsHYZ4(3,2,2,2,11)
     $  + s12 * CsHYZ4(3,2,2,2,12)
     $  + s13 * CsHYZ4(3,2,2,2,13)
     $  + s14 * CsHYZ4(3,2,2,2,14)
     $  + s15 * CsHYZ4(3,2,2,2,15)
     $  + s16 * CsHYZ4(3,2,2,2,16)
     $  + s17 * CsHYZ4(3,2,2,2,17)
     $  - HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,3,0,2) =
     $  + s01 * CsHYZ4(3,3,0,2,01)
     $  + s02 * CsHYZ4(3,3,0,2,02)
     $  + s03 * CsHYZ4(3,3,0,2,03)
     $  + s04 * CsHYZ4(3,3,0,2,04)
     $  + s05 * CsHYZ4(3,3,0,2,05)
     $  + s06 * CsHYZ4(3,3,0,2,06)
     $  + s07 * CsHYZ4(3,3,0,2,07)
     $  + s08 * CsHYZ4(3,3,0,2,08)
     $  + s09 * CsHYZ4(3,3,0,2,09)
     $  + s10 * CsHYZ4(3,3,0,2,10)
     $  + s11 * CsHYZ4(3,3,0,2,11)
     $  + s12 * CsHYZ4(3,3,0,2,12)
     $  + s13 * CsHYZ4(3,3,0,2,13)
     $  + s14 * CsHYZ4(3,3,0,2,14)
     $  + s15 * CsHYZ4(3,3,0,2,15)
     $  + s16 * CsHYZ4(3,3,0,2,16)
     $  + s17 * CsHYZ4(3,3,0,2,17)
     $  - HZ2(0,1) *HYZ2(3,3)
     $  - HZ2(1,1) *HYZ2(3,3)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,3,2,2) =
     $  + s01 * CsHYZ4(3,3,2,2,01)
     $  + s02 * CsHYZ4(3,3,2,2,02)
     $  + s03 * CsHYZ4(3,3,2,2,03)
     $  + s04 * CsHYZ4(3,3,2,2,04)
     $  + s05 * CsHYZ4(3,3,2,2,05)
     $  + s06 * CsHYZ4(3,3,2,2,06)
     $  + s07 * CsHYZ4(3,3,2,2,07)
     $  + s08 * CsHYZ4(3,3,2,2,08)
     $  + s09 * CsHYZ4(3,3,2,2,09)
     $  + s10 * CsHYZ4(3,3,2,2,10)
     $  + s11 * CsHYZ4(3,3,2,2,11)
     $  + s12 * CsHYZ4(3,3,2,2,12)
     $  + s13 * CsHYZ4(3,3,2,2,13)
     $  + s14 * CsHYZ4(3,3,2,2,14)
     $  + s15 * CsHYZ4(3,3,2,2,15)
     $  + s16 * CsHYZ4(3,3,2,2,16)
     $  + s17 * CsHYZ4(3,3,2,2,17)
     $  + HZ2(1,1) *HYZ2(3,3)
     $  + HZ3(0,1,1) *HYZ1(3)
     $  + HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,3,3,2) =
     $  + s01 * CsHYZ4(3,3,3,2,01)
     $  + s02 * CsHYZ4(3,3,3,2,02)
     $  + s03 * CsHYZ4(3,3,3,2,03)
     $  + s04 * CsHYZ4(3,3,3,2,04)
     $  + s05 * CsHYZ4(3,3,3,2,05)
     $  + s06 * CsHYZ4(3,3,3,2,06)
     $  + s07 * CsHYZ4(3,3,3,2,07)
     $  + s08 * CsHYZ4(3,3,3,2,08)
     $  + s09 * CsHYZ4(3,3,3,2,09)
     $  + s10 * CsHYZ4(3,3,3,2,10)
     $  + s11 * CsHYZ4(3,3,3,2,11)
     $  + s12 * CsHYZ4(3,3,3,2,12)
     $  + s13 * CsHYZ4(3,3,3,2,13)
     $  + s14 * CsHYZ4(3,3,3,2,14)
     $  + s15 * CsHYZ4(3,3,3,2,15)
     $  + s16 * CsHYZ4(3,3,3,2,16)
     $  + s17 * CsHYZ4(3,3,3,2,17)
     $  - HZ1(1) *HYZ3(3,3,3)
     $  - HZ2(0,1) *HYZ2(3,3)
     $  - HZ3(0,0,1) *HYZ1(3)
      endif

      return
      end


      subroutine fillirr2dhpl320e(iflag,n)
*********************************************************************
***** fillirr2dhpl320e(iflag,n) fills the irreducible 2dHPL     *****
***** with indices (3,2,0)                                      *****
***** up to weight n using the EXACT expressions for the        *****
***** z-dependent expansion coefficients                        *****
***** applicable for z>0.5                                      *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22
      dimension CbHYZ1(3:3,22),CbHYZ2(0:3,0:3,22),
     $          CbHYZ3(0:3,0:3,0:3,22),CbHYZ4(0:3,0:3,0:3,0:3,22)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common /b/
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22
      common/aux/CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4
!$omp threadprivate(/HPL2/,/b/,/aux/)

      if (iflag.eq.1) then
      call fillcoeff2dhplaux(iflag,n,CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4)
      call fillcoeff2dhpl320e(iflag,n,CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4)
      endif
* 2001-04-20:13:08:47.hpl
* <- HPLexpand.320e
* produced by form-to-fortr for gehrt@pcth62

      if (n.eq.2) then
      HYZ2(3,2) =
     $  + b01 * CbHYZ2(3,2,01)
     $  + b02 * CbHYZ2(3,2,02)
     $  + b03 * CbHYZ2(3,2,03)
     $  + b04 * CbHYZ2(3,2,04)
     $  + b05 * CbHYZ2(3,2,05)
     $  + b06 * CbHYZ2(3,2,06)
     $  + b07 * CbHYZ2(3,2,07)
     $  + b08 * CbHYZ2(3,2,08)
     $  + b09 * CbHYZ2(3,2,09)
     $  + b10 * CbHYZ2(3,2,10)
     $  + b11 * CbHYZ2(3,2,11)
     $  + b12 * CbHYZ2(3,2,12)
     $  + b13 * CbHYZ2(3,2,13)
     $  + b14 * CbHYZ2(3,2,14)
     $  + b15 * CbHYZ2(3,2,15)
     $  + b16 * CbHYZ2(3,2,16)
     $  + b17 * CbHYZ2(3,2,17)
     $  - HZ1(1) *HYZ1(3)
      endif
      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(0,3,2) =
     $  + b01 * CbHYZ3(0,3,2,01)
     $  + b02 * CbHYZ3(0,3,2,02)
     $  + b03 * CbHYZ3(0,3,2,03)
     $  + b04 * CbHYZ3(0,3,2,04)
     $  + b05 * CbHYZ3(0,3,2,05)
     $  + b06 * CbHYZ3(0,3,2,06)
     $  + b07 * CbHYZ3(0,3,2,07)
     $  + b08 * CbHYZ3(0,3,2,08)
     $  + b09 * CbHYZ3(0,3,2,09)
     $  + b10 * CbHYZ3(0,3,2,10)
     $  + b11 * CbHYZ3(0,3,2,11)
     $  + b12 * CbHYZ3(0,3,2,12)
     $  + b13 * CbHYZ3(0,3,2,13)
     $  + b14 * CbHYZ3(0,3,2,14)
     $  + b15 * CbHYZ3(0,3,2,15)
     $  + b16 * CbHYZ3(0,3,2,16)
     $  + b17 * CbHYZ3(0,3,2,17)
     $  - HZ1(1) *HYZ2(0,3)
      HYZ3(3,0,2) =
     $  + b01 * CbHYZ3(3,0,2,01)
     $  + b02 * CbHYZ3(3,0,2,02)
     $  + b03 * CbHYZ3(3,0,2,03)
     $  + b04 * CbHYZ3(3,0,2,04)
     $  + b05 * CbHYZ3(3,0,2,05)
     $  + b06 * CbHYZ3(3,0,2,06)
     $  + b07 * CbHYZ3(3,0,2,07)
     $  + b08 * CbHYZ3(3,0,2,08)
     $  + b09 * CbHYZ3(3,0,2,09)
     $  + b10 * CbHYZ3(3,0,2,10)
     $  + b11 * CbHYZ3(3,0,2,11)
     $  + b12 * CbHYZ3(3,0,2,12)
     $  + b13 * CbHYZ3(3,0,2,13)
     $  + b14 * CbHYZ3(3,0,2,14)
     $  + b15 * CbHYZ3(3,0,2,15)
     $  + b16 * CbHYZ3(3,0,2,16)
     $  + b17 * CbHYZ3(3,0,2,17)
     $  - HZ2(0,1) *HYZ1(3)
     $  - HZ2(1,1) *HYZ1(3)
      HYZ3(3,2,2) =
     $  + b01 * CbHYZ3(3,2,2,01)
     $  + b02 * CbHYZ3(3,2,2,02)
     $  + b03 * CbHYZ3(3,2,2,03)
     $  + b04 * CbHYZ3(3,2,2,04)
     $  + b05 * CbHYZ3(3,2,2,05)
     $  + b06 * CbHYZ3(3,2,2,06)
     $  + b07 * CbHYZ3(3,2,2,07)
     $  + b08 * CbHYZ3(3,2,2,08)
     $  + b09 * CbHYZ3(3,2,2,09)
     $  + b10 * CbHYZ3(3,2,2,10)
     $  + b11 * CbHYZ3(3,2,2,11)
     $  + b12 * CbHYZ3(3,2,2,12)
     $  + b13 * CbHYZ3(3,2,2,13)
     $  + b14 * CbHYZ3(3,2,2,14)
     $  + b15 * CbHYZ3(3,2,2,15)
     $  + b16 * CbHYZ3(3,2,2,16)
     $  + b17 * CbHYZ3(3,2,2,17)
     $  + b18 * CbHYZ3(3,2,2,18)
     $  + HZ2(1,1) *HYZ1(3)
      HYZ3(3,3,2) =
     $  + b01 * CbHYZ3(3,3,2,01)
     $  + b02 * CbHYZ3(3,3,2,02)
     $  + b03 * CbHYZ3(3,3,2,03)
     $  + b04 * CbHYZ3(3,3,2,04)
     $  + b05 * CbHYZ3(3,3,2,05)
     $  + b06 * CbHYZ3(3,3,2,06)
     $  + b07 * CbHYZ3(3,3,2,07)
     $  + b08 * CbHYZ3(3,3,2,08)
     $  + b09 * CbHYZ3(3,3,2,09)
     $  + b10 * CbHYZ3(3,3,2,10)
     $  + b11 * CbHYZ3(3,3,2,11)
     $  + b12 * CbHYZ3(3,3,2,12)
     $  + b13 * CbHYZ3(3,3,2,13)
     $  + b14 * CbHYZ3(3,3,2,14)
     $  + b15 * CbHYZ3(3,3,2,15)
     $  + b16 * CbHYZ3(3,3,2,16)
     $  + b17 * CbHYZ3(3,3,2,17)
     $  - HZ1(1) *HYZ2(3,3)
     $  - HZ2(0,1) *HYZ1(3)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,0,3,2) =
     $  + b01 * CbHYZ4(0,0,3,2,01)
     $  + b02 * CbHYZ4(0,0,3,2,02)
     $  + b03 * CbHYZ4(0,0,3,2,03)
     $  + b04 * CbHYZ4(0,0,3,2,04)
     $  + b05 * CbHYZ4(0,0,3,2,05)
     $  + b06 * CbHYZ4(0,0,3,2,06)
     $  + b07 * CbHYZ4(0,0,3,2,07)
     $  + b08 * CbHYZ4(0,0,3,2,08)
     $  + b09 * CbHYZ4(0,0,3,2,09)
     $  + b10 * CbHYZ4(0,0,3,2,10)
     $  + b11 * CbHYZ4(0,0,3,2,11)
     $  + b12 * CbHYZ4(0,0,3,2,12)
     $  + b13 * CbHYZ4(0,0,3,2,13)
     $  + b14 * CbHYZ4(0,0,3,2,14)
     $  + b15 * CbHYZ4(0,0,3,2,15)
     $  + b16 * CbHYZ4(0,0,3,2,16)
     $  + b17 * CbHYZ4(0,0,3,2,17)
     $  - HZ1(1) *HYZ3(0,0,3)
      HYZ4(0,2,3,2) =
     $  + b01 * CbHYZ4(0,2,3,2,01)
     $  + b02 * CbHYZ4(0,2,3,2,02)
     $  + b03 * CbHYZ4(0,2,3,2,03)
     $  + b04 * CbHYZ4(0,2,3,2,04)
     $  + b05 * CbHYZ4(0,2,3,2,05)
     $  + b06 * CbHYZ4(0,2,3,2,06)
     $  + b07 * CbHYZ4(0,2,3,2,07)
     $  + b08 * CbHYZ4(0,2,3,2,08)
     $  + b09 * CbHYZ4(0,2,3,2,09)
     $  + b10 * CbHYZ4(0,2,3,2,10)
     $  + b11 * CbHYZ4(0,2,3,2,11)
     $  + b12 * CbHYZ4(0,2,3,2,12)
     $  + b13 * CbHYZ4(0,2,3,2,13)
     $  + b14 * CbHYZ4(0,2,3,2,14)
     $  + b15 * CbHYZ4(0,2,3,2,15)
     $  + b16 * CbHYZ4(0,2,3,2,16)
     $  + b17 * CbHYZ4(0,2,3,2,17)
     $  - HZ1(1) *HYZ3(0,2,3)
      HYZ4(0,3,0,2) =
     $  + b01 * CbHYZ4(0,3,0,2,01)
     $  + b02 * CbHYZ4(0,3,0,2,02)
     $  + b03 * CbHYZ4(0,3,0,2,03)
     $  + b04 * CbHYZ4(0,3,0,2,04)
     $  + b05 * CbHYZ4(0,3,0,2,05)
     $  + b06 * CbHYZ4(0,3,0,2,06)
     $  + b07 * CbHYZ4(0,3,0,2,07)
     $  + b08 * CbHYZ4(0,3,0,2,08)
     $  + b09 * CbHYZ4(0,3,0,2,09)
     $  + b10 * CbHYZ4(0,3,0,2,10)
     $  + b11 * CbHYZ4(0,3,0,2,11)
     $  + b12 * CbHYZ4(0,3,0,2,12)
     $  + b13 * CbHYZ4(0,3,0,2,13)
     $  + b14 * CbHYZ4(0,3,0,2,14)
     $  + b15 * CbHYZ4(0,3,0,2,15)
     $  + b16 * CbHYZ4(0,3,0,2,16)
     $  + b17 * CbHYZ4(0,3,0,2,17)
     $  - HZ2(0,1) *HYZ2(0,3)
     $  - HZ2(1,1) *HYZ2(0,3)
      HYZ4(0,3,2,2) =
     $  + b01 * CbHYZ4(0,3,2,2,01)
     $  + b02 * CbHYZ4(0,3,2,2,02)
     $  + b03 * CbHYZ4(0,3,2,2,03)
     $  + b04 * CbHYZ4(0,3,2,2,04)
     $  + b05 * CbHYZ4(0,3,2,2,05)
     $  + b06 * CbHYZ4(0,3,2,2,06)
     $  + b07 * CbHYZ4(0,3,2,2,07)
     $  + b08 * CbHYZ4(0,3,2,2,08)
     $  + b09 * CbHYZ4(0,3,2,2,09)
     $  + b10 * CbHYZ4(0,3,2,2,10)
     $  + b11 * CbHYZ4(0,3,2,2,11)
     $  + b12 * CbHYZ4(0,3,2,2,12)
     $  + b13 * CbHYZ4(0,3,2,2,13)
     $  + b14 * CbHYZ4(0,3,2,2,14)
     $  + b15 * CbHYZ4(0,3,2,2,15)
     $  + b16 * CbHYZ4(0,3,2,2,16)
     $  + b17 * CbHYZ4(0,3,2,2,17)
     $  + HZ2(1,1) *HYZ2(0,3)
      HYZ4(0,3,3,2) =
     $  + b01 * CbHYZ4(0,3,3,2,01)
     $  + b02 * CbHYZ4(0,3,3,2,02)
     $  + b03 * CbHYZ4(0,3,3,2,03)
     $  + b04 * CbHYZ4(0,3,3,2,04)
     $  + b05 * CbHYZ4(0,3,3,2,05)
     $  + b06 * CbHYZ4(0,3,3,2,06)
     $  + b07 * CbHYZ4(0,3,3,2,07)
     $  + b08 * CbHYZ4(0,3,3,2,08)
     $  + b09 * CbHYZ4(0,3,3,2,09)
     $  + b10 * CbHYZ4(0,3,3,2,10)
     $  + b11 * CbHYZ4(0,3,3,2,11)
     $  + b12 * CbHYZ4(0,3,3,2,12)
     $  + b13 * CbHYZ4(0,3,3,2,13)
     $  + b14 * CbHYZ4(0,3,3,2,14)
     $  + b15 * CbHYZ4(0,3,3,2,15)
     $  + b16 * CbHYZ4(0,3,3,2,16)
     $  + b17 * CbHYZ4(0,3,3,2,17)
     $  - HZ1(1) *HYZ3(0,3,3)
     $  - HZ2(0,1) *HYZ2(0,3)
      HYZ4(3,0,0,2) =
     $  + b01 * CbHYZ4(3,0,0,2,01)
     $  + b02 * CbHYZ4(3,0,0,2,02)
     $  + b03 * CbHYZ4(3,0,0,2,03)
     $  + b04 * CbHYZ4(3,0,0,2,04)
     $  + b05 * CbHYZ4(3,0,0,2,05)
     $  + b06 * CbHYZ4(3,0,0,2,06)
     $  + b07 * CbHYZ4(3,0,0,2,07)
     $  + b08 * CbHYZ4(3,0,0,2,08)
     $  + b09 * CbHYZ4(3,0,0,2,09)
     $  + b10 * CbHYZ4(3,0,0,2,10)
     $  + b11 * CbHYZ4(3,0,0,2,11)
     $  + b12 * CbHYZ4(3,0,0,2,12)
     $  + b13 * CbHYZ4(3,0,0,2,13)
     $  + b14 * CbHYZ4(3,0,0,2,14)
     $  + b15 * CbHYZ4(3,0,0,2,15)
     $  + b16 * CbHYZ4(3,0,0,2,16)
     $  + b17 * CbHYZ4(3,0,0,2,17)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(0,1,1) *HYZ1(3)
     $  - HZ3(1,0,1) *HYZ1(3)
     $  - HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,0,2,2) =
     $  + b01 * CbHYZ4(3,0,2,2,01)
     $  + b02 * CbHYZ4(3,0,2,2,02)
     $  + b03 * CbHYZ4(3,0,2,2,03)
     $  + b04 * CbHYZ4(3,0,2,2,04)
     $  + b05 * CbHYZ4(3,0,2,2,05)
     $  + b06 * CbHYZ4(3,0,2,2,06)
     $  + b07 * CbHYZ4(3,0,2,2,07)
     $  + b08 * CbHYZ4(3,0,2,2,08)
     $  + b09 * CbHYZ4(3,0,2,2,09)
     $  + b10 * CbHYZ4(3,0,2,2,10)
     $  + b11 * CbHYZ4(3,0,2,2,11)
     $  + b12 * CbHYZ4(3,0,2,2,12)
     $  + b13 * CbHYZ4(3,0,2,2,13)
     $  + b14 * CbHYZ4(3,0,2,2,14)
     $  + b15 * CbHYZ4(3,0,2,2,15)
     $  + b16 * CbHYZ4(3,0,2,2,16)
     $  + b17 * CbHYZ4(3,0,2,2,17)
     $  + HZ3(0,1,1) *HYZ1(3)
     $  + HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,0,3,2) =
     $  + b01 * CbHYZ4(3,0,3,2,01)
     $  + b02 * CbHYZ4(3,0,3,2,02)
     $  + b03 * CbHYZ4(3,0,3,2,03)
     $  + b04 * CbHYZ4(3,0,3,2,04)
     $  + b05 * CbHYZ4(3,0,3,2,05)
     $  + b06 * CbHYZ4(3,0,3,2,06)
     $  + b07 * CbHYZ4(3,0,3,2,07)
     $  + b08 * CbHYZ4(3,0,3,2,08)
     $  + b09 * CbHYZ4(3,0,3,2,09)
     $  + b10 * CbHYZ4(3,0,3,2,10)
     $  + b11 * CbHYZ4(3,0,3,2,11)
     $  + b12 * CbHYZ4(3,0,3,2,12)
     $  + b13 * CbHYZ4(3,0,3,2,13)
     $  + b14 * CbHYZ4(3,0,3,2,14)
     $  + b15 * CbHYZ4(3,0,3,2,15)
     $  + b16 * CbHYZ4(3,0,3,2,16)
     $  + b17 * CbHYZ4(3,0,3,2,17)
     $  - HZ1(1) *HYZ3(3,0,3)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(0,1,1) *HYZ1(3)
      HYZ4(3,2,2,2) =
     $  + b01 * CbHYZ4(3,2,2,2,01)
     $  + b02 * CbHYZ4(3,2,2,2,02)
     $  + b03 * CbHYZ4(3,2,2,2,03)
     $  + b04 * CbHYZ4(3,2,2,2,04)
     $  + b05 * CbHYZ4(3,2,2,2,05)
     $  + b06 * CbHYZ4(3,2,2,2,06)
     $  + b07 * CbHYZ4(3,2,2,2,07)
     $  + b08 * CbHYZ4(3,2,2,2,08)
     $  + b09 * CbHYZ4(3,2,2,2,09)
     $  + b10 * CbHYZ4(3,2,2,2,10)
     $  + b11 * CbHYZ4(3,2,2,2,11)
     $  + b12 * CbHYZ4(3,2,2,2,12)
     $  + b13 * CbHYZ4(3,2,2,2,13)
     $  + b14 * CbHYZ4(3,2,2,2,14)
     $  + b15 * CbHYZ4(3,2,2,2,15)
     $  + b16 * CbHYZ4(3,2,2,2,16)
     $  + b17 * CbHYZ4(3,2,2,2,17)
     $  + b18 * CbHYZ4(3,2,2,2,18)
     $  - HZ3(1,1,1) *HYZ1(3)
      HYZ4(3,3,0,2) =
     $  + b01 * CbHYZ4(3,3,0,2,01)
     $  + b02 * CbHYZ4(3,3,0,2,02)
     $  + b03 * CbHYZ4(3,3,0,2,03)
     $  + b04 * CbHYZ4(3,3,0,2,04)
     $  + b05 * CbHYZ4(3,3,0,2,05)
     $  + b06 * CbHYZ4(3,3,0,2,06)
     $  + b07 * CbHYZ4(3,3,0,2,07)
     $  + b08 * CbHYZ4(3,3,0,2,08)
     $  + b09 * CbHYZ4(3,3,0,2,09)
     $  + b10 * CbHYZ4(3,3,0,2,10)
     $  + b11 * CbHYZ4(3,3,0,2,11)
     $  + b12 * CbHYZ4(3,3,0,2,12)
     $  + b13 * CbHYZ4(3,3,0,2,13)
     $  + b14 * CbHYZ4(3,3,0,2,14)
     $  + b15 * CbHYZ4(3,3,0,2,15)
     $  + b16 * CbHYZ4(3,3,0,2,16)
     $  + b17 * CbHYZ4(3,3,0,2,17)
     $  - HZ2(0,1) *HYZ2(3,3)
     $  - HZ2(1,1) *HYZ2(3,3)
     $  - HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,3,2,2) =
     $  + b01 * CbHYZ4(3,3,2,2,01)
     $  + b02 * CbHYZ4(3,3,2,2,02)
     $  + b03 * CbHYZ4(3,3,2,2,03)
     $  + b04 * CbHYZ4(3,3,2,2,04)
     $  + b05 * CbHYZ4(3,3,2,2,05)
     $  + b06 * CbHYZ4(3,3,2,2,06)
     $  + b07 * CbHYZ4(3,3,2,2,07)
     $  + b08 * CbHYZ4(3,3,2,2,08)
     $  + b09 * CbHYZ4(3,3,2,2,09)
     $  + b10 * CbHYZ4(3,3,2,2,10)
     $  + b11 * CbHYZ4(3,3,2,2,11)
     $  + b12 * CbHYZ4(3,3,2,2,12)
     $  + b13 * CbHYZ4(3,3,2,2,13)
     $  + b14 * CbHYZ4(3,3,2,2,14)
     $  + b15 * CbHYZ4(3,3,2,2,15)
     $  + b16 * CbHYZ4(3,3,2,2,16)
     $  + b17 * CbHYZ4(3,3,2,2,17)
     $  + b18 * CbHYZ4(3,3,2,2,18)
     $  + HZ2(1,1) *HYZ2(3,3)
     $  + HZ3(0,1,1) *HYZ1(3)
     $  + HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,3,3,2) =
     $  + b01 * CbHYZ4(3,3,3,2,01)
     $  + b02 * CbHYZ4(3,3,3,2,02)
     $  + b03 * CbHYZ4(3,3,3,2,03)
     $  + b04 * CbHYZ4(3,3,3,2,04)
     $  + b05 * CbHYZ4(3,3,3,2,05)
     $  + b06 * CbHYZ4(3,3,3,2,06)
     $  + b07 * CbHYZ4(3,3,3,2,07)
     $  + b08 * CbHYZ4(3,3,3,2,08)
     $  + b09 * CbHYZ4(3,3,3,2,09)
     $  + b10 * CbHYZ4(3,3,3,2,10)
     $  + b11 * CbHYZ4(3,3,3,2,11)
     $  + b12 * CbHYZ4(3,3,3,2,12)
     $  + b13 * CbHYZ4(3,3,3,2,13)
     $  + b14 * CbHYZ4(3,3,3,2,14)
     $  + b15 * CbHYZ4(3,3,3,2,15)
     $  + b16 * CbHYZ4(3,3,3,2,16)
     $  + b17 * CbHYZ4(3,3,3,2,17)
     $  + b18 * CbHYZ4(3,3,3,2,18)
     $  - HZ1(1) *HYZ3(3,3,3)
     $  - HZ2(0,1) *HYZ2(3,3)
     $  - HZ3(0,0,1) *HYZ1(3)
      endif
      return

      end

      subroutine fillirr2dhpl321(iflag,n)
*********************************************************************
***** fillirr2dhpl321(iflag,n) fills the irreducible 2dHPL      *****
***** with indices (3,2,1)                                      *****
***** up to weight n using the ONLY EXPANDED contributions      *****
***** to the z-dependent expansion coefficients                 *****
*****   applicable for z<0.5                                    *****
***** treatment of terms non-analytic in z=1:                   *****
*****        in expanded form (z<0.5: call fillcoeff2dhpl321u)  *****
***** this routine is also invoked by                           *****
*****    fillirr2dhpl321e, which handles the exact coefficients *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: CsHYZ2,CsHYZ3,CsHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      dimension CsHYZ2(0:3,0:3,22),CsHYZ3(0:3,0:3,0:3,22),
     $          CsHYZ4(0:3,0:3,0:3,0:3,22)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common /s/
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      save CsHYZ2,CsHYZ3,CsHYZ4
!$omp threadprivate(CsHYZ2,CsHYZ3,CsHYZ4,/s/,/HPL2/)

      if (iflag.eq.1.or.iflag.eq.-1) then
      call fillcoeff2dhpl321(iflag,n,CsHYZ2,CsHYZ3,CsHYZ4)
      if (iflag.gt.0)
     $    call fillcoeff2dhpl321u(iflag,n,CsHYZ2,CsHYZ3,CsHYZ4)
      endif

* 2001-04-20:16:55:01.hpl
* <- HPLexpand.321tcut
* produced by form-to-fortr for gehrt@pcth78

      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(2,3,1) =
     $  + s02 * CsHYZ3(2,3,1,02)
     $  + s03 * CsHYZ3(2,3,1,03)
     $  + s04 * CsHYZ3(2,3,1,04)
     $  + s05 * CsHYZ3(2,3,1,05)
     $  + s06 * CsHYZ3(2,3,1,06)
     $  + s07 * CsHYZ3(2,3,1,07)
     $  + s08 * CsHYZ3(2,3,1,08)
     $  + s09 * CsHYZ3(2,3,1,09)
     $  + s10 * CsHYZ3(2,3,1,10)
     $  + s11 * CsHYZ3(2,3,1,11)
     $  + s12 * CsHYZ3(2,3,1,12)
     $  + s13 * CsHYZ3(2,3,1,13)
     $  + s14 * CsHYZ3(2,3,1,14)
     $  + s15 * CsHYZ3(2,3,1,15)
     $  + s16 * CsHYZ3(2,3,1,16)
     $  + s17 * CsHYZ3(2,3,1,17)
     $  + s18 * CsHYZ3(2,3,1,18)
     $  + s19 * CsHYZ3(2,3,1,19)
     $  - HZ1( -1)*HYZ2(2,3)
      HYZ3(3,2,1) =
     $  + s01 * CsHYZ3(3,2,1,01)
     $  + s02 * CsHYZ3(3,2,1,02)
     $  + s03 * CsHYZ3(3,2,1,03)
     $  + s04 * CsHYZ3(3,2,1,04)
     $  + s05 * CsHYZ3(3,2,1,05)
     $  + s06 * CsHYZ3(3,2,1,06)
     $  + s07 * CsHYZ3(3,2,1,07)
     $  + s08 * CsHYZ3(3,2,1,08)
     $  + s09 * CsHYZ3(3,2,1,09)
     $  + s10 * CsHYZ3(3,2,1,10)
     $  + s11 * CsHYZ3(3,2,1,11)
     $  + s12 * CsHYZ3(3,2,1,12)
     $  + s13 * CsHYZ3(3,2,1,13)
     $  + s14 * CsHYZ3(3,2,1,14)
     $  + s15 * CsHYZ3(3,2,1,15)
     $  + s16 * CsHYZ3(3,2,1,16)
     $  + s17 * CsHYZ3(3,2,1,17)
     $  + s18 * CsHYZ3(3,2,1,18)
     $  + s19 * CsHYZ3(3,2,1,19)
     $  - HZ2(0, -1)*HYZ1(3)
     $  + HZ2(0,1) *HYZ1(3)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,2,3,1) =
     $  + s02 * CsHYZ4(0,2,3,1,02)
     $  + s03 * CsHYZ4(0,2,3,1,03)
     $  + s04 * CsHYZ4(0,2,3,1,04)
     $  + s05 * CsHYZ4(0,2,3,1,05)
     $  + s06 * CsHYZ4(0,2,3,1,06)
     $  + s07 * CsHYZ4(0,2,3,1,07)
     $  + s08 * CsHYZ4(0,2,3,1,08)
     $  + s09 * CsHYZ4(0,2,3,1,09)
     $  + s10 * CsHYZ4(0,2,3,1,10)
     $  + s11 * CsHYZ4(0,2,3,1,11)
     $  + s12 * CsHYZ4(0,2,3,1,12)
     $  + s13 * CsHYZ4(0,2,3,1,13)
     $  + s14 * CsHYZ4(0,2,3,1,14)
     $  + s15 * CsHYZ4(0,2,3,1,15)
     $  + s16 * CsHYZ4(0,2,3,1,16)
     $  + s17 * CsHYZ4(0,2,3,1,17)
     $  + s18 * CsHYZ4(0,2,3,1,18)
     $  - HZ1( -1)*HYZ3(0,2,3)
      HYZ4(0,3,2,1) =
     $  + s01 * CsHYZ4(0,3,2,1,01)
     $  + s02 * CsHYZ4(0,3,2,1,02)
     $  + s03 * CsHYZ4(0,3,2,1,03)
     $  + s04 * CsHYZ4(0,3,2,1,04)
     $  + s05 * CsHYZ4(0,3,2,1,05)
     $  + s06 * CsHYZ4(0,3,2,1,06)
     $  + s07 * CsHYZ4(0,3,2,1,07)
     $  + s08 * CsHYZ4(0,3,2,1,08)
     $  + s09 * CsHYZ4(0,3,2,1,09)
     $  + s10 * CsHYZ4(0,3,2,1,10)
     $  + s11 * CsHYZ4(0,3,2,1,11)
     $  + s12 * CsHYZ4(0,3,2,1,12)
     $  + s13 * CsHYZ4(0,3,2,1,13)
     $  + s14 * CsHYZ4(0,3,2,1,14)
     $  + s15 * CsHYZ4(0,3,2,1,15)
     $  + s16 * CsHYZ4(0,3,2,1,16)
     $  + s17 * CsHYZ4(0,3,2,1,17)
     $  + s18 * CsHYZ4(0,3,2,1,18)
     $  - HZ2(0, -1)*HYZ2(0,3)
     $  + HZ2(0,1) *HYZ2(0,3)
      HYZ4(2,0,3,1) =
     $  + s02 * CsHYZ4(2,0,3,1,02)
     $  + s03 * CsHYZ4(2,0,3,1,03)
     $  + s04 * CsHYZ4(2,0,3,1,04)
     $  + s05 * CsHYZ4(2,0,3,1,05)
     $  + s06 * CsHYZ4(2,0,3,1,06)
     $  + s07 * CsHYZ4(2,0,3,1,07)
     $  + s08 * CsHYZ4(2,0,3,1,08)
     $  + s09 * CsHYZ4(2,0,3,1,09)
     $  + s10 * CsHYZ4(2,0,3,1,10)
     $  + s11 * CsHYZ4(2,0,3,1,11)
     $  + s12 * CsHYZ4(2,0,3,1,12)
     $  + s13 * CsHYZ4(2,0,3,1,13)
     $  + s14 * CsHYZ4(2,0,3,1,14)
     $  + s15 * CsHYZ4(2,0,3,1,15)
     $  + s16 * CsHYZ4(2,0,3,1,16)
     $  + s17 * CsHYZ4(2,0,3,1,17)
     $  + s18 * CsHYZ4(2,0,3,1,18)
     $  - HZ1( -1)*HYZ3(2,0,3)
      HYZ4(2,2,3,1) =
     $  + s03 * CsHYZ4(2,2,3,1,03)
     $  + s04 * CsHYZ4(2,2,3,1,04)
     $  + s05 * CsHYZ4(2,2,3,1,05)
     $  + s06 * CsHYZ4(2,2,3,1,06)
     $  + s07 * CsHYZ4(2,2,3,1,07)
     $  + s08 * CsHYZ4(2,2,3,1,08)
     $  + s09 * CsHYZ4(2,2,3,1,09)
     $  + s10 * CsHYZ4(2,2,3,1,10)
     $  + s11 * CsHYZ4(2,2,3,1,11)
     $  + s12 * CsHYZ4(2,2,3,1,12)
     $  + s13 * CsHYZ4(2,2,3,1,13)
     $  + s14 * CsHYZ4(2,2,3,1,14)
     $  + s15 * CsHYZ4(2,2,3,1,15)
     $  + s16 * CsHYZ4(2,2,3,1,16)
     $  + s17 * CsHYZ4(2,2,3,1,17)
     $  + s18 * CsHYZ4(2,2,3,1,18)
     $  - HZ1( -1)*HYZ3(2,2,3)
      HYZ4(2,3,0,1) =
     $  + s02 * CsHYZ4(2,3,0,1,02)
     $  + s03 * CsHYZ4(2,3,0,1,03)
     $  + s04 * CsHYZ4(2,3,0,1,04)
     $  + s05 * CsHYZ4(2,3,0,1,05)
     $  + s06 * CsHYZ4(2,3,0,1,06)
     $  + s07 * CsHYZ4(2,3,0,1,07)
     $  + s08 * CsHYZ4(2,3,0,1,08)
     $  + s09 * CsHYZ4(2,3,0,1,09)
     $  + s10 * CsHYZ4(2,3,0,1,10)
     $  + s11 * CsHYZ4(2,3,0,1,11)
     $  + s12 * CsHYZ4(2,3,0,1,12)
     $  + s13 * CsHYZ4(2,3,0,1,13)
     $  + s14 * CsHYZ4(2,3,0,1,14)
     $  + s15 * CsHYZ4(2,3,0,1,15)
     $  + s16 * CsHYZ4(2,3,0,1,16)
     $  + s17 * CsHYZ4(2,3,0,1,17)
     $  + s18 * CsHYZ4(2,3,0,1,18)
     $  - HZ2(0, -1)*HYZ2(2,3)
      HYZ4(2,3,1,1) =
     $  + s02 * CsHYZ4(2,3,1,1,02)
     $  + s03 * CsHYZ4(2,3,1,1,03)
     $  + s04 * CsHYZ4(2,3,1,1,04)
     $  + s05 * CsHYZ4(2,3,1,1,05)
     $  + s06 * CsHYZ4(2,3,1,1,06)
     $  + s07 * CsHYZ4(2,3,1,1,07)
     $  + s08 * CsHYZ4(2,3,1,1,08)
     $  + s09 * CsHYZ4(2,3,1,1,09)
     $  + s10 * CsHYZ4(2,3,1,1,10)
     $  + s11 * CsHYZ4(2,3,1,1,11)
     $  + s12 * CsHYZ4(2,3,1,1,12)
     $  + s13 * CsHYZ4(2,3,1,1,13)
     $  + s14 * CsHYZ4(2,3,1,1,14)
     $  + s15 * CsHYZ4(2,3,1,1,15)
     $  + s16 * CsHYZ4(2,3,1,1,16)
     $  + s17 * CsHYZ4(2,3,1,1,17)
     $  + s18 * CsHYZ4(2,3,1,1,18)
     $  + s19 * CsHYZ4(2,3,1,1,19)
     $  + s20 * CsHYZ4(2,3,1,1,20)
     $  + HZ2( -1,-1)*HYZ2(2,3)
      HYZ4(2,3,2,1) =
     $  + s02 * CsHYZ4(2,3,2,1,02)
     $  + s03 * CsHYZ4(2,3,2,1,03)
     $  + s04 * CsHYZ4(2,3,2,1,04)
     $  + s05 * CsHYZ4(2,3,2,1,05)
     $  + s06 * CsHYZ4(2,3,2,1,06)
     $  + s07 * CsHYZ4(2,3,2,1,07)
     $  + s08 * CsHYZ4(2,3,2,1,08)
     $  + s09 * CsHYZ4(2,3,2,1,09)
     $  + s10 * CsHYZ4(2,3,2,1,10)
     $  + s11 * CsHYZ4(2,3,2,1,11)
     $  + s12 * CsHYZ4(2,3,2,1,12)
     $  + s13 * CsHYZ4(2,3,2,1,13)
     $  + s14 * CsHYZ4(2,3,2,1,14)
     $  + s15 * CsHYZ4(2,3,2,1,15)
     $  + s16 * CsHYZ4(2,3,2,1,16)
     $  + s17 * CsHYZ4(2,3,2,1,17)
     $  + s18 * CsHYZ4(2,3,2,1,18)
     $  - HZ2(0, -1)*HYZ2(2,3)
     $  + HZ2(0,1) *HYZ2(2,3)
      HYZ4(2,3,3,1) =
     $  + s02 * CsHYZ4(2,3,3,1,02)
     $  + s03 * CsHYZ4(2,3,3,1,03)
     $  + s04 * CsHYZ4(2,3,3,1,04)
     $  + s05 * CsHYZ4(2,3,3,1,05)
     $  + s06 * CsHYZ4(2,3,3,1,06)
     $  + s07 * CsHYZ4(2,3,3,1,07)
     $  + s08 * CsHYZ4(2,3,3,1,08)
     $  + s09 * CsHYZ4(2,3,3,1,09)
     $  + s10 * CsHYZ4(2,3,3,1,10)
     $  + s11 * CsHYZ4(2,3,3,1,11)
     $  + s12 * CsHYZ4(2,3,3,1,12)
     $  + s13 * CsHYZ4(2,3,3,1,13)
     $  + s14 * CsHYZ4(2,3,3,1,14)
     $  + s15 * CsHYZ4(2,3,3,1,15)
     $  + s16 * CsHYZ4(2,3,3,1,16)
     $  + s17 * CsHYZ4(2,3,3,1,17)
     $  + s18 * CsHYZ4(2,3,3,1,18)
     $  - HZ1( -1)*HYZ3(2,3,3)
     $  + HZ2( -1,-1)*HYZ2(2,3)
     $  - HZ2(0, -1)*HYZ2(2,3)
      HYZ4(3,0,2,1) =
     $  + s01 * CsHYZ4(3,0,2,1,01)
     $  + s02 * CsHYZ4(3,0,2,1,02)
     $  + s03 * CsHYZ4(3,0,2,1,03)
     $  + s04 * CsHYZ4(3,0,2,1,04)
     $  + s05 * CsHYZ4(3,0,2,1,05)
     $  + s06 * CsHYZ4(3,0,2,1,06)
     $  + s07 * CsHYZ4(3,0,2,1,07)
     $  + s08 * CsHYZ4(3,0,2,1,08)
     $  + s09 * CsHYZ4(3,0,2,1,09)
     $  + s10 * CsHYZ4(3,0,2,1,10)
     $  + s11 * CsHYZ4(3,0,2,1,11)
     $  + s12 * CsHYZ4(3,0,2,1,12)
     $  + s13 * CsHYZ4(3,0,2,1,13)
     $  + s14 * CsHYZ4(3,0,2,1,14)
     $  + s15 * CsHYZ4(3,0,2,1,15)
     $  + s16 * CsHYZ4(3,0,2,1,16)
     $  + s17 * CsHYZ4(3,0,2,1,17)
     $  + s18 * CsHYZ4(3,0,2,1,18)
     $  - 2.000000000000000d+00*HZ3(0,0,-1)*HYZ1(3)
     $  + 2.000000000000000d+00*HZ3(0,0,1)*HYZ1(3)
     $  + HZ3(0,1,1) *HYZ1(3)
     $  - 2.000000000000000d+00*HZ3(1,0,-1)*HYZ1(3)
     $  + HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,1,2,1) =
     $  + s01 * CsHYZ4(3,1,2,1,01)
     $  + s02 * CsHYZ4(3,1,2,1,02)
     $  + s03 * CsHYZ4(3,1,2,1,03)
     $  + s04 * CsHYZ4(3,1,2,1,04)
     $  + s05 * CsHYZ4(3,1,2,1,05)
     $  + s06 * CsHYZ4(3,1,2,1,06)
     $  + s07 * CsHYZ4(3,1,2,1,07)
     $  + s08 * CsHYZ4(3,1,2,1,08)
     $  + s09 * CsHYZ4(3,1,2,1,09)
     $  + s10 * CsHYZ4(3,1,2,1,10)
     $  + s11 * CsHYZ4(3,1,2,1,11)
     $  + s12 * CsHYZ4(3,1,2,1,12)
     $  + s13 * CsHYZ4(3,1,2,1,13)
     $  + s14 * CsHYZ4(3,1,2,1,14)
     $  + s15 * CsHYZ4(3,1,2,1,15)
     $  + s16 * CsHYZ4(3,1,2,1,16)
     $  + s17 * CsHYZ4(3,1,2,1,17)
     $  + s18 * CsHYZ4(3,1,2,1,18)
     $  + s19 * CsHYZ4(3,1,2,1,19)
     $  + s20 * CsHYZ4(3,1,2,1,20)
     $  + s21 * CsHYZ4(3,1,2,1,21)
     $  + HZ3( -1,0,-1)*HYZ1(3)
     $  - HZ3( -1,0,1)*HYZ1(3)
     $  - HZ3(0, -1,1)*HYZ1(3)
     $  - 2.000000000000000d+00*HZ3(0,0,-1)*HYZ1(3)
     $  + 2.000000000000000d+00*HZ3(0,0,1)*HYZ1(3)
     $  - HZ3(0,1, -1)*HYZ1(3)
      HYZ4(3,2,0,1) =
     $  + s01 * CsHYZ4(3,2,0,1,01)
     $  + s02 * CsHYZ4(3,2,0,1,02)
     $  + s03 * CsHYZ4(3,2,0,1,03)
     $  + s04 * CsHYZ4(3,2,0,1,04)
     $  + s05 * CsHYZ4(3,2,0,1,05)
     $  + s06 * CsHYZ4(3,2,0,1,06)
     $  + s07 * CsHYZ4(3,2,0,1,07)
     $  + s08 * CsHYZ4(3,2,0,1,08)
     $  + s09 * CsHYZ4(3,2,0,1,09)
     $  + s10 * CsHYZ4(3,2,0,1,10)
     $  + s11 * CsHYZ4(3,2,0,1,11)
     $  + s12 * CsHYZ4(3,2,0,1,12)
     $  + s13 * CsHYZ4(3,2,0,1,13)
     $  + s14 * CsHYZ4(3,2,0,1,14)
     $  + s15 * CsHYZ4(3,2,0,1,15)
     $  + s16 * CsHYZ4(3,2,0,1,16)
     $  + s17 * CsHYZ4(3,2,0,1,17)
     $  + s18 * CsHYZ4(3,2,0,1,18)
     $  + 2.000000000000000d+00*HZ3(1,0,-1)*HYZ1(3)
     $  - HZ3(1,0,1) *HYZ1(3)
      HYZ4(3,2,1,1) =
     $  + s01 * CsHYZ4(3,2,1,1,01)
     $  + s02 * CsHYZ4(3,2,1,1,02)
     $  + s03 * CsHYZ4(3,2,1,1,03)
     $  + s04 * CsHYZ4(3,2,1,1,04)
     $  + s05 * CsHYZ4(3,2,1,1,05)
     $  + s06 * CsHYZ4(3,2,1,1,06)
     $  + s07 * CsHYZ4(3,2,1,1,07)
     $  + s08 * CsHYZ4(3,2,1,1,08)
     $  + s09 * CsHYZ4(3,2,1,1,09)
     $  + s10 * CsHYZ4(3,2,1,1,10)
     $  + s11 * CsHYZ4(3,2,1,1,11)
     $  + s12 * CsHYZ4(3,2,1,1,12)
     $  + s13 * CsHYZ4(3,2,1,1,13)
     $  + s14 * CsHYZ4(3,2,1,1,14)
     $  + s15 * CsHYZ4(3,2,1,1,15)
     $  + s16 * CsHYZ4(3,2,1,1,16)
     $  + s17 * CsHYZ4(3,2,1,1,17)
     $  + s18 * CsHYZ4(3,2,1,1,18)
     $  + s19 * CsHYZ4(3,2,1,1,19)
     $  + s20 * CsHYZ4(3,2,1,1,20)
     $  + HZ3(0, -1,-1)*HYZ1(3)
     $  + HZ3(0,0, -1)*HYZ1(3)
     $  - HZ3(0,0,1) *HYZ1(3)
      HYZ4(3,2,2,1) =
     $  + s01 * CsHYZ4(3,2,2,1,01)
     $  + s02 * CsHYZ4(3,2,2,1,02)
     $  + s03 * CsHYZ4(3,2,2,1,03)
     $  + s04 * CsHYZ4(3,2,2,1,04)
     $  + s05 * CsHYZ4(3,2,2,1,05)
     $  + s06 * CsHYZ4(3,2,2,1,06)
     $  + s07 * CsHYZ4(3,2,2,1,07)
     $  + s08 * CsHYZ4(3,2,2,1,08)
     $  + s09 * CsHYZ4(3,2,2,1,09)
     $  + s10 * CsHYZ4(3,2,2,1,10)
     $  + s11 * CsHYZ4(3,2,2,1,11)
     $  + s12 * CsHYZ4(3,2,2,1,12)
     $  + s13 * CsHYZ4(3,2,2,1,13)
     $  + s14 * CsHYZ4(3,2,2,1,14)
     $  + s15 * CsHYZ4(3,2,2,1,15)
     $  + s16 * CsHYZ4(3,2,2,1,16)
     $  + s17 * CsHYZ4(3,2,2,1,17)
     $  + s18 * CsHYZ4(3,2,2,1,18)
     $  - HZ3(0,0, -1)*HYZ1(3)
     $  + HZ3(0,0,1) *HYZ1(3)
     $  - HZ3(0,1,1) *HYZ1(3)
      HYZ4(3,2,3,1) =
     $  + s01 * CsHYZ4(3,2,3,1,01)
     $  + s02 * CsHYZ4(3,2,3,1,02)
     $  + s03 * CsHYZ4(3,2,3,1,03)
     $  + s04 * CsHYZ4(3,2,3,1,04)
     $  + s05 * CsHYZ4(3,2,3,1,05)
     $  + s06 * CsHYZ4(3,2,3,1,06)
     $  + s07 * CsHYZ4(3,2,3,1,07)
     $  + s08 * CsHYZ4(3,2,3,1,08)
     $  + s09 * CsHYZ4(3,2,3,1,09)
     $  + s10 * CsHYZ4(3,2,3,1,10)
     $  + s11 * CsHYZ4(3,2,3,1,11)
     $  + s12 * CsHYZ4(3,2,3,1,12)
     $  + s13 * CsHYZ4(3,2,3,1,13)
     $  + s14 * CsHYZ4(3,2,3,1,14)
     $  + s15 * CsHYZ4(3,2,3,1,15)
     $  + s16 * CsHYZ4(3,2,3,1,16)
     $  + s17 * CsHYZ4(3,2,3,1,17)
     $  + s18 * CsHYZ4(3,2,3,1,18)
     $  - HZ1( -1)*HYZ3(3,2,3)
     $  + HZ3( -1,0,-1)*HYZ1(3)
     $  - HZ3( -1,0,1)*HYZ1(3)
     $  + HZ3(0, -1,1)*HYZ1(3)
     $  + HZ3(0,1, -1)*HYZ1(3)
      HYZ4(3,3,2,1) =
     $  + s01 * CsHYZ4(3,3,2,1,01)
     $  + s02 * CsHYZ4(3,3,2,1,02)
     $  + s03 * CsHYZ4(3,3,2,1,03)
     $  + s04 * CsHYZ4(3,3,2,1,04)
     $  + s05 * CsHYZ4(3,3,2,1,05)
     $  + s06 * CsHYZ4(3,3,2,1,06)
     $  + s07 * CsHYZ4(3,3,2,1,07)
     $  + s08 * CsHYZ4(3,3,2,1,08)
     $  + s09 * CsHYZ4(3,3,2,1,09)
     $  + s10 * CsHYZ4(3,3,2,1,10)
     $  + s11 * CsHYZ4(3,3,2,1,11)
     $  + s12 * CsHYZ4(3,3,2,1,12)
     $  + s13 * CsHYZ4(3,3,2,1,13)
     $  + s14 * CsHYZ4(3,3,2,1,14)
     $  + s15 * CsHYZ4(3,3,2,1,15)
     $  + s16 * CsHYZ4(3,3,2,1,16)
     $  + s17 * CsHYZ4(3,3,2,1,17)
     $  + s18 * CsHYZ4(3,3,2,1,18)
     $  - HZ2(0, -1)*HYZ2(3,3)
     $  + HZ2(0,1) *HYZ2(3,3)
     $  + HZ3(0, -1,-1)*HYZ1(3)
     $  - 2.000000000000000d+00*HZ3(0,0,-1)*HYZ1(3)
     $  + 2.000000000000000d+00*HZ3(0,0,1)*HYZ1(3)
      endif

      return
      end

      subroutine fillirr2dhpl321e(iflag,n)
*********************************************************************
***** fillirr2dhpl321e(iflag,n) fills the irreducible 2dHPL     *****
***** with indices (3,2,1)                                      *****
***** up to weight n with EXPANDED z-dependent                  *****
***** expansion coefficients for all terms                      *****
***** which are analytic in z=1                                 *****
***** and EXACT z-dependent                                     *****
***** expansion coefficients for all terms                      *****
***** terms non-analytic in z=1,                                *****
*****       (z>0.5: call fillcoeff2dhpl321e)                    *****
*****   applicable for z>0.5                                    *****
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp)::
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22
      dimension CbHYZ1(3:3,22),CbHYZ2(0:3,0:3,22),
     $          CbHYZ3(0:3,0:3,0:3,22),CbHYZ4(0:3,0:3,0:3,0:3,22)
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
      common/aux/CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4
      common /b/
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22
!$omp threadprivate(/HPL2/,/b/,/aux/)

      iflag = -iflag
      call fillirr2dhpl321(iflag,n)
      iflag = -iflag
      if (iflag.eq.1) then
      call fillcoeff2dhplaux(iflag,n,CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4)
      call fillcoeff2dhpl321e(iflag,n,CbHYZ1,CbHYZ2,CbHYZ3,CbHYZ4)
      endif

* 2001-04-20:17:05:41.hpl
* <- HPLexpand.321enew
* produced by form-to-fortr for gehrt@pcth78

      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(3,2,1) = HYZ3(3,2,1)
     $  + b01 * CbHYZ3(3,2,1,01)
     $  + b02 * CbHYZ3(3,2,1,02)
     $  + b03 * CbHYZ3(3,2,1,03)
     $  + b04 * CbHYZ3(3,2,1,04)
     $  + b05 * CbHYZ3(3,2,1,05)
     $  + b06 * CbHYZ3(3,2,1,06)
     $  + b07 * CbHYZ3(3,2,1,07)
     $  + b08 * CbHYZ3(3,2,1,08)
     $  + b09 * CbHYZ3(3,2,1,09)
     $  + b10 * CbHYZ3(3,2,1,10)
     $  + b11 * CbHYZ3(3,2,1,11)
     $  + b12 * CbHYZ3(3,2,1,12)
     $  + b13 * CbHYZ3(3,2,1,13)
     $  + b14 * CbHYZ3(3,2,1,14)
     $  + b15 * CbHYZ3(3,2,1,15)
     $  + b16 * CbHYZ3(3,2,1,16)
     $  + b17 * CbHYZ3(3,2,1,17)
     $  + b18 * CbHYZ3(3,2,1,18)
     $  + b19 * CbHYZ3(3,2,1,19)
     $  + b20 * CbHYZ3(3,2,1,20)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,3,2,1) = HYZ4(0,3,2,1)
     $  + b01 * CbHYZ4(0,3,2,1,01)
     $  + b02 * CbHYZ4(0,3,2,1,02)
     $  + b03 * CbHYZ4(0,3,2,1,03)
     $  + b04 * CbHYZ4(0,3,2,1,04)
     $  + b05 * CbHYZ4(0,3,2,1,05)
     $  + b06 * CbHYZ4(0,3,2,1,06)
     $  + b07 * CbHYZ4(0,3,2,1,07)
     $  + b08 * CbHYZ4(0,3,2,1,08)
     $  + b09 * CbHYZ4(0,3,2,1,09)
     $  + b10 * CbHYZ4(0,3,2,1,10)
     $  + b11 * CbHYZ4(0,3,2,1,11)
     $  + b12 * CbHYZ4(0,3,2,1,12)
     $  + b13 * CbHYZ4(0,3,2,1,13)
     $  + b14 * CbHYZ4(0,3,2,1,14)
     $  + b15 * CbHYZ4(0,3,2,1,15)
     $  + b16 * CbHYZ4(0,3,2,1,16)
     $  + b17 * CbHYZ4(0,3,2,1,17)
     $  + b18 * CbHYZ4(0,3,2,1,18)
     $  + b19 * CbHYZ4(0,3,2,1,19)
      HYZ4(2,3,2,1) = HYZ4(2,3,2,1)
     $  + b01 * CbHYZ4(2,3,2,1,01)
     $  + b02 * CbHYZ4(2,3,2,1,02)
     $  + b03 * CbHYZ4(2,3,2,1,03)
     $  + b04 * CbHYZ4(2,3,2,1,04)
     $  + b05 * CbHYZ4(2,3,2,1,05)
     $  + b06 * CbHYZ4(2,3,2,1,06)
     $  + b07 * CbHYZ4(2,3,2,1,07)
     $  + b08 * CbHYZ4(2,3,2,1,08)
     $  + b09 * CbHYZ4(2,3,2,1,09)
     $  + b10 * CbHYZ4(2,3,2,1,10)
     $  + b11 * CbHYZ4(2,3,2,1,11)
     $  + b12 * CbHYZ4(2,3,2,1,12)
     $  + b13 * CbHYZ4(2,3,2,1,13)
     $  + b14 * CbHYZ4(2,3,2,1,14)
     $  + b15 * CbHYZ4(2,3,2,1,15)
     $  + b16 * CbHYZ4(2,3,2,1,16)
     $  + b17 * CbHYZ4(2,3,2,1,17)
     $  + b18 * CbHYZ4(2,3,2,1,18)
     $  + b19 * CbHYZ4(2,3,2,1,19)
      HYZ4(3,0,2,1) = HYZ4(3,0,2,1)
     $  + b01 * CbHYZ4(3,0,2,1,01)
     $  + b02 * CbHYZ4(3,0,2,1,02)
     $  + b03 * CbHYZ4(3,0,2,1,03)
     $  + b04 * CbHYZ4(3,0,2,1,04)
     $  + b05 * CbHYZ4(3,0,2,1,05)
     $  + b06 * CbHYZ4(3,0,2,1,06)
     $  + b07 * CbHYZ4(3,0,2,1,07)
     $  + b08 * CbHYZ4(3,0,2,1,08)
     $  + b09 * CbHYZ4(3,0,2,1,09)
     $  + b10 * CbHYZ4(3,0,2,1,10)
     $  + b11 * CbHYZ4(3,0,2,1,11)
     $  + b12 * CbHYZ4(3,0,2,1,12)
     $  + b13 * CbHYZ4(3,0,2,1,13)
     $  + b14 * CbHYZ4(3,0,2,1,14)
     $  + b15 * CbHYZ4(3,0,2,1,15)
     $  + b16 * CbHYZ4(3,0,2,1,16)
     $  + b17 * CbHYZ4(3,0,2,1,17)
     $  + b18 * CbHYZ4(3,0,2,1,18)
     $  + b19 * CbHYZ4(3,0,2,1,19)
      HYZ4(3,1,2,1) = HYZ4(3,1,2,1)
     $  + b01 * CbHYZ4(3,1,2,1,01)
     $  + b02 * CbHYZ4(3,1,2,1,02)
     $  + b03 * CbHYZ4(3,1,2,1,03)
     $  + b04 * CbHYZ4(3,1,2,1,04)
     $  + b05 * CbHYZ4(3,1,2,1,05)
     $  + b06 * CbHYZ4(3,1,2,1,06)
     $  + b07 * CbHYZ4(3,1,2,1,07)
     $  + b08 * CbHYZ4(3,1,2,1,08)
     $  + b09 * CbHYZ4(3,1,2,1,09)
     $  + b10 * CbHYZ4(3,1,2,1,10)
     $  + b11 * CbHYZ4(3,1,2,1,11)
     $  + b12 * CbHYZ4(3,1,2,1,12)
     $  + b13 * CbHYZ4(3,1,2,1,13)
     $  + b14 * CbHYZ4(3,1,2,1,14)
     $  + b15 * CbHYZ4(3,1,2,1,15)
     $  + b16 * CbHYZ4(3,1,2,1,16)
     $  + b17 * CbHYZ4(3,1,2,1,17)
     $  + b18 * CbHYZ4(3,1,2,1,18)
     $  + b19 * CbHYZ4(3,1,2,1,19)
     $  + b20 * CbHYZ4(3,1,2,1,20)
     $  + b21 * CbHYZ4(3,1,2,1,21)
     $  + b22 * CbHYZ4(3,1,2,1,22)
      HYZ4(3,2,0,1) = HYZ4(3,2,0,1)
     $  + b01 * CbHYZ4(3,2,0,1,01)
     $  + b02 * CbHYZ4(3,2,0,1,02)
     $  + b03 * CbHYZ4(3,2,0,1,03)
     $  + b04 * CbHYZ4(3,2,0,1,04)
     $  + b05 * CbHYZ4(3,2,0,1,05)
     $  + b06 * CbHYZ4(3,2,0,1,06)
     $  + b07 * CbHYZ4(3,2,0,1,07)
     $  + b08 * CbHYZ4(3,2,0,1,08)
     $  + b09 * CbHYZ4(3,2,0,1,09)
     $  + b10 * CbHYZ4(3,2,0,1,10)
     $  + b11 * CbHYZ4(3,2,0,1,11)
     $  + b12 * CbHYZ4(3,2,0,1,12)
     $  + b13 * CbHYZ4(3,2,0,1,13)
     $  + b14 * CbHYZ4(3,2,0,1,14)
     $  + b15 * CbHYZ4(3,2,0,1,15)
     $  + b16 * CbHYZ4(3,2,0,1,16)
     $  + b17 * CbHYZ4(3,2,0,1,17)
     $  + b18 * CbHYZ4(3,2,0,1,18)
     $  + b19 * CbHYZ4(3,2,0,1,19)
      HYZ4(3,2,1,1) = HYZ4(3,2,1,1)
     $  + b01 * CbHYZ4(3,2,1,1,01)
     $  + b02 * CbHYZ4(3,2,1,1,02)
     $  + b03 * CbHYZ4(3,2,1,1,03)
     $  + b04 * CbHYZ4(3,2,1,1,04)
     $  + b05 * CbHYZ4(3,2,1,1,05)
     $  + b06 * CbHYZ4(3,2,1,1,06)
     $  + b07 * CbHYZ4(3,2,1,1,07)
     $  + b08 * CbHYZ4(3,2,1,1,08)
     $  + b09 * CbHYZ4(3,2,1,1,09)
     $  + b10 * CbHYZ4(3,2,1,1,10)
     $  + b11 * CbHYZ4(3,2,1,1,11)
     $  + b12 * CbHYZ4(3,2,1,1,12)
     $  + b13 * CbHYZ4(3,2,1,1,13)
     $  + b14 * CbHYZ4(3,2,1,1,14)
     $  + b15 * CbHYZ4(3,2,1,1,15)
     $  + b16 * CbHYZ4(3,2,1,1,16)
     $  + b17 * CbHYZ4(3,2,1,1,17)
     $  + b18 * CbHYZ4(3,2,1,1,18)
     $  + b19 * CbHYZ4(3,2,1,1,19)
     $  + b20 * CbHYZ4(3,2,1,1,20)
     $  + b21 * CbHYZ4(3,2,1,1,21)
      HYZ4(3,2,2,1) = HYZ4(3,2,2,1)
     $  + b01 * CbHYZ4(3,2,2,1,01)
     $  + b02 * CbHYZ4(3,2,2,1,02)
     $  + b03 * CbHYZ4(3,2,2,1,03)
     $  + b04 * CbHYZ4(3,2,2,1,04)
     $  + b05 * CbHYZ4(3,2,2,1,05)
     $  + b06 * CbHYZ4(3,2,2,1,06)
     $  + b07 * CbHYZ4(3,2,2,1,07)
     $  + b08 * CbHYZ4(3,2,2,1,08)
     $  + b09 * CbHYZ4(3,2,2,1,09)
     $  + b10 * CbHYZ4(3,2,2,1,10)
     $  + b11 * CbHYZ4(3,2,2,1,11)
     $  + b12 * CbHYZ4(3,2,2,1,12)
     $  + b13 * CbHYZ4(3,2,2,1,13)
     $  + b14 * CbHYZ4(3,2,2,1,14)
     $  + b15 * CbHYZ4(3,2,2,1,15)
     $  + b16 * CbHYZ4(3,2,2,1,16)
     $  + b17 * CbHYZ4(3,2,2,1,17)
     $  + b18 * CbHYZ4(3,2,2,1,18)
     $  + b19 * CbHYZ4(3,2,2,1,19)
      HYZ4(3,2,3,1) = HYZ4(3,2,3,1)
     $  + b01 * CbHYZ4(3,2,3,1,01)
     $  + b02 * CbHYZ4(3,2,3,1,02)
     $  + b03 * CbHYZ4(3,2,3,1,03)
     $  + b04 * CbHYZ4(3,2,3,1,04)
     $  + b05 * CbHYZ4(3,2,3,1,05)
     $  + b06 * CbHYZ4(3,2,3,1,06)
     $  + b07 * CbHYZ4(3,2,3,1,07)
     $  + b08 * CbHYZ4(3,2,3,1,08)
     $  + b09 * CbHYZ4(3,2,3,1,09)
     $  + b10 * CbHYZ4(3,2,3,1,10)
     $  + b11 * CbHYZ4(3,2,3,1,11)
     $  + b12 * CbHYZ4(3,2,3,1,12)
     $  + b13 * CbHYZ4(3,2,3,1,13)
     $  + b14 * CbHYZ4(3,2,3,1,14)
     $  + b15 * CbHYZ4(3,2,3,1,15)
     $  + b16 * CbHYZ4(3,2,3,1,16)
     $  + b17 * CbHYZ4(3,2,3,1,17)
     $  + b18 * CbHYZ4(3,2,3,1,18)
     $  + b19 * CbHYZ4(3,2,3,1,19)
      HYZ4(3,3,2,1) = HYZ4(3,3,2,1)
     $  + b01 * CbHYZ4(3,3,2,1,01)
     $  + b02 * CbHYZ4(3,3,2,1,02)
     $  + b03 * CbHYZ4(3,3,2,1,03)
     $  + b04 * CbHYZ4(3,3,2,1,04)
     $  + b05 * CbHYZ4(3,3,2,1,05)
     $  + b06 * CbHYZ4(3,3,2,1,06)
     $  + b07 * CbHYZ4(3,3,2,1,07)
     $  + b08 * CbHYZ4(3,3,2,1,08)
     $  + b09 * CbHYZ4(3,3,2,1,09)
     $  + b10 * CbHYZ4(3,3,2,1,10)
     $  + b11 * CbHYZ4(3,3,2,1,11)
     $  + b12 * CbHYZ4(3,3,2,1,12)
     $  + b13 * CbHYZ4(3,3,2,1,13)
     $  + b14 * CbHYZ4(3,3,2,1,14)
     $  + b15 * CbHYZ4(3,3,2,1,15)
     $  + b16 * CbHYZ4(3,3,2,1,16)
     $  + b17 * CbHYZ4(3,3,2,1,17)
     $  + b18 * CbHYZ4(3,3,2,1,18)
     $  + b19 * CbHYZ4(3,3,2,1,19)
      endif

      return
      end


      subroutine fillred2dhpl(iflag,n)
*********************************************************************
*****  fillred2dhpl evaluates the reducible 2dhpl              ******
*****  irreducible 2dhpl need to be present in HPL2 already    ******
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,n
      real(dp):: HY1,HY2,HY3,HY4,HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      common /HPL2/
     $ HY1(-1:1),HY2(-1:1,-1:1),HY3(-1:1,-1:1,-1:1),
     $           HY4(-1:1,-1:1,-1:1,-1:1),
     $ HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $           HZ4(-1:1,-1:1,-1:1,-1:1),
     $ HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $           HYZ4(0:3,0:3,0:3,0:3)
!$omp threadprivate(/HPL2/)

* 2001-04-23:14:13:53.hpl
* <- red.out
* produced by form-to-fortr for gehrt@pcth62

      if (n.eq.2) then
      HYZ2(0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)
      HYZ2(1,0) =
     $  + HYZ1(0) *HYZ1(1)
     $  - HYZ2(0,1)
      HYZ2(1,1) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)
      HYZ2(1,2) =
     $  + HYZ1(1) *HYZ1(2)
     $  - HYZ2(2,1)
      HYZ2(1,3) =
     $  + HYZ1(1) *HYZ1(3)
     $  - HYZ2(3,1)
      HYZ2(2,0) =
     $  + HYZ1(0) *HYZ1(2)
     $  - HYZ2(0,2)
      HYZ2(2,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)
      HYZ2(2,3) =
     $  + HYZ1(2) *HYZ1(3)
     $  - HYZ2(3,2)
      HYZ2(3,0) =
     $  + HYZ1(0) *HYZ1(3)
     $  - HYZ2(0,3)
      HYZ2(3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)
      endif
      if (n.eq.2) return

      if (n.eq.3) then
      HYZ3(0,0,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(0)*HYZ1(0)
      HYZ3(0,1,0) =
     $  + HYZ1(0) *HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ3(0,0,1)
      HYZ3(0,1,2) =
     $  + HYZ1(2) *HYZ2(0,1)
     $  - HYZ3(0,2,1)
     $  - HYZ3(2,0,1)
      HYZ3(0,1,3) =
     $  + HYZ1(3) *HYZ2(0,1)
     $  - HYZ3(0,3,1)
     $  - HYZ3(3,0,1)
      HYZ3(0,2,0) =
     $  + HYZ1(0) *HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ3(0,0,2)
      HYZ3(0,2,3) =
     $  + HYZ1(3) *HYZ2(0,2)
     $  - HYZ3(0,3,2)
     $  - HYZ3(3,0,2)
      HYZ3(0,3,0) =
     $  + HYZ1(0) *HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ3(0,0,3)
      HYZ3(1,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(1)
     $  - HYZ1(0) *HYZ2(0,1)
     $  + HYZ3(0,0,1)
      HYZ3(1,0,1) =
     $  + HYZ1(1) *HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ3(0,1,1)
      HYZ3(1,0,2) =
     $  + HYZ1(1) *HYZ2(0,2)
     $  - HYZ1(2) *HYZ2(0,1)
     $  + HYZ3(2,0,1)
      HYZ3(1,0,3) =
     $  + HYZ1(1) *HYZ2(0,3)
     $  - HYZ1(3) *HYZ2(0,1)
     $  + HYZ3(3,0,1)
      HYZ3(1,1,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(1)*HYZ1(1)
     $  - HYZ1(1) *HYZ2(0,1)
     $  + HYZ3(0,1,1)
      HYZ3(1,1,1) =
     $  + 1.666666666666666d-01*HYZ1(1)*HYZ1(1)*HYZ1(1)
      HYZ3(1,1,2) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ1(2)
     $  - HYZ1(1) *HYZ2(2,1)
     $  + HYZ3(2,1,1)
      HYZ3(1,1,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ1(3)
     $  - HYZ1(1) *HYZ2(3,1)
     $  + HYZ3(3,1,1)
      HYZ3(1,2,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ1(2)
     $  - HYZ1(0) *HYZ2(2,1)
     $  - HYZ1(1) *HYZ2(0,2)
     $  + HYZ3(0,2,1)
      HYZ3(1,2,1) =
     $  + HYZ1(1) *HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ3(2,1,1)
      HYZ3(1,2,2) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(2)*HYZ1(2)
     $  - HYZ1(2) *HYZ2(2,1)
     $  + HYZ3(2,2,1)
      HYZ3(1,2,3) =
     $  + HYZ1(1) *HYZ1(2)*HYZ1(3)
     $  - HYZ1(1) *HYZ2(3,2)
     $  - HYZ1(3) *HYZ2(2,1)
     $  + HYZ3(3,2,1)
      HYZ3(1,3,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ1(3)
     $  - HYZ1(0) *HYZ2(3,1)
     $  - HYZ1(1) *HYZ2(0,3)
     $  + HYZ3(0,3,1)
      HYZ3(1,3,1) =
     $  + HYZ1(1) *HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ3(3,1,1)
      HYZ3(1,3,2) =
     $  + HYZ1(1) *HYZ2(3,2)
     $  - HYZ1(2) *HYZ2(3,1)
     $  + HYZ3(2,3,1)
      HYZ3(1,3,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(3) *HYZ2(3,1)
     $  + HYZ3(3,3,1)
      HYZ3(2,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(2)
     $  - HYZ1(0) *HYZ2(0,2)
     $  + HYZ3(0,0,2)
      HYZ3(2,0,2) =
     $  + HYZ1(2) *HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ3(0,2,2)
      HYZ3(2,0,3) =
     $  + HYZ1(2) *HYZ2(0,3)
     $  - HYZ1(3) *HYZ2(0,2)
     $  + HYZ3(3,0,2)
      HYZ3(2,1,0) =
     $  + HYZ1(0) *HYZ2(2,1)
     $  - HYZ3(0,2,1)
     $  - HYZ3(2,0,1)
      HYZ3(2,1,2) =
     $  + HYZ1(2) *HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ3(2,2,1)
      HYZ3(2,1,3) =
     $  + HYZ1(3) *HYZ2(2,1)
     $  - HYZ3(2,3,1)
     $  - HYZ3(3,2,1)
      HYZ3(2,2,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(2)*HYZ1(2)
     $  - HYZ1(2) *HYZ2(0,2)
     $  + HYZ3(0,2,2)
      HYZ3(2,2,2) =
     $  + 1.666666666666666d-01*HYZ1(2)*HYZ1(2)*HYZ1(2)
      HYZ3(2,2,3) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ1(3)
     $  - HYZ1(2) *HYZ2(3,2)
     $  + HYZ3(3,2,2)
      HYZ3(2,3,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ1(3)
     $  - HYZ1(0) *HYZ2(3,2)
     $  - HYZ1(2) *HYZ2(0,3)
     $  + HYZ3(0,3,2)
      HYZ3(2,3,2) =
     $  + HYZ1(2) *HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ3(3,2,2)
      HYZ3(2,3,3) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(3) *HYZ2(3,2)
     $  + HYZ3(3,3,2)
      HYZ3(3,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(3)
     $  - HYZ1(0) *HYZ2(0,3)
     $  + HYZ3(0,0,3)
      HYZ3(3,0,3) =
     $  + HYZ1(3) *HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ3(0,3,3)
      HYZ3(3,1,0) =
     $  + HYZ1(0) *HYZ2(3,1)
     $  - HYZ3(0,3,1)
     $  - HYZ3(3,0,1)
      HYZ3(3,1,2) =
     $  + HYZ1(2) *HYZ2(3,1)
     $  - HYZ3(2,3,1)
     $  - HYZ3(3,2,1)
      HYZ3(3,1,3) =
     $  + HYZ1(3) *HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ3(3,3,1)
      HYZ3(3,2,0) =
     $  + HYZ1(0) *HYZ2(3,2)
     $  - HYZ3(0,3,2)
     $  - HYZ3(3,0,2)
      HYZ3(3,2,3) =
     $  + HYZ1(3) *HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ3(3,3,2)
      HYZ3(3,3,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(3) *HYZ2(0,3)
     $  + HYZ3(0,3,3)
      HYZ3(3,3,3) =
     $  + 1.666666666666666d-01*HYZ1(3)*HYZ1(3)*HYZ1(3)
      endif
      if (n.eq.3) return

      if (n.eq.4) then
      HYZ4(0,0,0,0) =
     $  + 4.166666666666666d-02*HYZ1(0)*HYZ1(0)*HYZ1(0)*HYZ1(0)
      HYZ4(0,0,1,0) =
     $  + HYZ1(0) *HYZ3(0,0,1)
     $  - 3.000000000000000d+00*HYZ4(0,0,0,1)
      HYZ4(0,0,1,2) =
     $  + HYZ1(2) *HYZ3(0,0,1)
     $  - HYZ4(0,0,2,1)
     $  - HYZ4(0,2,0,1)
     $  - HYZ4(2,0,0,1)
      HYZ4(0,0,1,3) =
     $  + HYZ1(3) *HYZ3(0,0,1)
     $  - HYZ4(0,0,3,1)
     $  - HYZ4(0,3,0,1)
     $  - HYZ4(3,0,0,1)
      HYZ4(0,0,2,0) =
     $  + HYZ1(0) *HYZ3(0,0,2)
     $  - 3.000000000000000d+00*HYZ4(0,0,0,2)
      HYZ4(0,0,2,3) =
     $  + HYZ1(3) *HYZ3(0,0,2)
     $  - HYZ4(0,0,3,2)
     $  - HYZ4(0,3,0,2)
     $  - HYZ4(3,0,0,2)
      HYZ4(0,0,3,0) =
     $  + HYZ1(0) *HYZ3(0,0,3)
     $  - 3.000000000000000d+00*HYZ4(0,0,0,3)
      HYZ4(0,1,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,0,1)
     $  + 3.000000000000000d+00*HYZ4(0,0,0,1)
      HYZ4(0,1,0,1) =
     $  + 5.000000000000000d-01*HYZ2(0,1)*HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ4(0,0,1,1)
      HYZ4(0,1,0,2) =
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,0,1)
     $  + HYZ2(0,1) *HYZ2(0,2)
     $  + HYZ4(0,2,0,1)
     $  + 2.000000000000000d+00*HYZ4(2,0,0,1)
      HYZ4(0,1,0,3) =
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,0,1)
     $  + HYZ2(0,1) *HYZ2(0,3)
     $  + HYZ4(0,3,0,1)
     $  + 2.000000000000000d+00*HYZ4(3,0,0,1)
      HYZ4(0,1,1,0) =
     $  + HYZ1(0) *HYZ3(0,1,1)
     $  - 5.000000000000000d-01*HYZ2(0,1)*HYZ2(0,1)
      HYZ4(0,1,1,2) =
     $  + HYZ1(2) *HYZ3(0,1,1)
     $  - HYZ4(0,1,2,1)
     $  - HYZ4(0,2,1,1)
     $  - HYZ4(2,0,1,1)
      HYZ4(0,1,1,3) =
     $  + HYZ1(3) *HYZ3(0,1,1)
     $  - HYZ4(0,1,3,1)
     $  - HYZ4(0,3,1,1)
     $  - HYZ4(3,0,1,1)
      HYZ4(0,1,2,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(0,1)
     $  - HYZ1(0) *HYZ3(0,2,1)
     $  - HYZ1(0) *HYZ3(2,0,1)
     $  - HYZ2(0,1) *HYZ2(0,2)
     $  + 2.000000000000000d+00*HYZ4(0,0,2,1)
     $  + HYZ4(0,2,0,1)
      HYZ4(0,1,2,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,1)
     $  - HYZ1(2) *HYZ3(0,2,1)
     $  - HYZ1(2) *HYZ3(2,0,1)
     $  + HYZ4(0,2,2,1)
     $  + HYZ4(2,0,2,1)
     $  + HYZ4(2,2,0,1)
      HYZ4(0,1,2,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(0,1)
     $  - HYZ1(3) *HYZ3(0,2,1)
     $  - HYZ1(3) *HYZ3(2,0,1)
     $  - HYZ2(0,1) *HYZ2(3,2)
     $  + HYZ4(0,3,2,1)
     $  + HYZ4(3,0,2,1)
     $  + HYZ4(3,2,0,1)
      HYZ4(0,1,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(0,1)
     $  - HYZ1(0) *HYZ3(0,3,1)
     $  - HYZ1(0) *HYZ3(3,0,1)
     $  - HYZ2(0,1) *HYZ2(0,3)
     $  + 2.000000000000000d+00*HYZ4(0,0,3,1)
     $  + HYZ4(0,3,0,1)
      HYZ4(0,1,3,2) =
     $  - HYZ1(2) *HYZ3(0,3,1)
     $  - HYZ1(2) *HYZ3(3,0,1)
     $  + HYZ2(0,1) *HYZ2(3,2)
     $  + HYZ4(0,2,3,1)
     $  + HYZ4(2,0,3,1)
     $  + HYZ4(2,3,0,1)
      HYZ4(0,1,3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,1)
     $  - HYZ1(3) *HYZ3(0,3,1)
     $  - HYZ1(3) *HYZ3(3,0,1)
     $  + HYZ4(0,3,3,1)
     $  + HYZ4(3,0,3,1)
     $  + HYZ4(3,3,0,1)
      HYZ4(0,2,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,0,2)
     $  + 3.000000000000000d+00*HYZ4(0,0,0,2)
      HYZ4(0,2,0,2) =
     $  + 5.000000000000000d-01*HYZ2(0,2)*HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ4(0,0,2,2)
      HYZ4(0,2,0,3) =
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,0,2)
     $  + HYZ2(0,2) *HYZ2(0,3)
     $  + HYZ4(0,3,0,2)
     $  + 2.000000000000000d+00*HYZ4(3,0,0,2)
      HYZ4(0,2,1,0) =
     $  + HYZ1(0) *HYZ3(0,2,1)
     $  - 2.000000000000000d+00*HYZ4(0,0,2,1)
     $  - HYZ4(0,2,0,1)
      HYZ4(0,2,1,2) =
     $  + HYZ1(2) *HYZ3(0,2,1)
     $  - 2.000000000000000d+00*HYZ4(0,2,2,1)
     $  - HYZ4(2,0,2,1)
      HYZ4(0,2,1,3) =
     $  + HYZ1(3) *HYZ3(0,2,1)
     $  - HYZ4(0,2,3,1)
     $  - HYZ4(0,3,2,1)
     $  - HYZ4(3,0,2,1)
      HYZ4(0,2,2,0) =
     $  + HYZ1(0) *HYZ3(0,2,2)
     $  - 5.000000000000000d-01*HYZ2(0,2)*HYZ2(0,2)
      HYZ4(0,2,2,3) =
     $  + HYZ1(3) *HYZ3(0,2,2)
     $  - HYZ4(0,2,3,2)
     $  - HYZ4(0,3,2,2)
     $  - HYZ4(3,0,2,2)
      HYZ4(0,2,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(0,2)
     $  - HYZ1(0) *HYZ3(0,3,2)
     $  - HYZ1(0) *HYZ3(3,0,2)
     $  - HYZ2(0,2) *HYZ2(0,3)
     $  + 2.000000000000000d+00*HYZ4(0,0,3,2)
     $  + HYZ4(0,3,0,2)
      HYZ4(0,2,3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,2)
     $  - HYZ1(3) *HYZ3(0,3,2)
     $  - HYZ1(3) *HYZ3(3,0,2)
     $  + HYZ4(0,3,3,2)
     $  + HYZ4(3,0,3,2)
     $  + HYZ4(3,3,0,2)
      HYZ4(0,3,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,0,3)
     $  + 3.000000000000000d+00*HYZ4(0,0,0,3)
      HYZ4(0,3,0,3) =
     $  + 5.000000000000000d-01*HYZ2(0,3)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ4(0,0,3,3)
      HYZ4(0,3,1,0) =
     $  + HYZ1(0) *HYZ3(0,3,1)
     $  - 2.000000000000000d+00*HYZ4(0,0,3,1)
     $  - HYZ4(0,3,0,1)
      HYZ4(0,3,1,2) =
     $  + HYZ1(2) *HYZ3(0,3,1)
     $  - HYZ4(0,2,3,1)
     $  - HYZ4(0,3,2,1)
     $  - HYZ4(2,0,3,1)
      HYZ4(0,3,1,3) =
     $  + HYZ1(3) *HYZ3(0,3,1)
     $  - 2.000000000000000d+00*HYZ4(0,3,3,1)
     $  - HYZ4(3,0,3,1)
      HYZ4(0,3,2,0) =
     $  + HYZ1(0) *HYZ3(0,3,2)
     $  - 2.000000000000000d+00*HYZ4(0,0,3,2)
     $  - HYZ4(0,3,0,2)
      HYZ4(0,3,2,3) =
     $  + HYZ1(3) *HYZ3(0,3,2)
     $  - 2.000000000000000d+00*HYZ4(0,3,3,2)
     $  - HYZ4(3,0,3,2)
      HYZ4(0,3,3,0) =
     $  + HYZ1(0) *HYZ3(0,3,3)
     $  - 5.000000000000000d-01*HYZ2(0,3)*HYZ2(0,3)
      HYZ4(1,0,0,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(0)*HYZ1(0)*HYZ1(1)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,1)
     $  + HYZ1(0) *HYZ3(0,0,1)
     $  - HYZ4(0,0,0,1)
      HYZ4(1,0,0,1) =
     $  + HYZ1(1) *HYZ3(0,0,1)
     $  - 5.000000000000000d-01*HYZ2(0,1)*HYZ2(0,1)
      HYZ4(1,0,0,2) =
     $  + HYZ1(1) *HYZ3(0,0,2)
     $  + HYZ1(2) *HYZ3(0,0,1)
     $  - HYZ2(0,1) *HYZ2(0,2)
     $  - HYZ4(2,0,0,1)
      HYZ4(1,0,0,3) =
     $  + HYZ1(1) *HYZ3(0,0,3)
     $  + HYZ1(3) *HYZ3(0,0,1)
     $  - HYZ2(0,1) *HYZ2(0,3)
     $  - HYZ4(3,0,0,1)
      HYZ4(1,0,1,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,1,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,0,1)
     $  + 5.000000000000000d-01*HYZ2(0,1)*HYZ2(0,1)
     $  + 2.000000000000000d+00*HYZ4(0,0,1,1)
      HYZ4(1,0,1,1) =
     $  + HYZ1(1) *HYZ3(0,1,1)
     $  - 3.000000000000000d+00*HYZ4(0,1,1,1)
      HYZ4(1,0,1,2) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(0,1)
     $  - HYZ1(1) *HYZ3(0,2,1)
     $  - HYZ1(1) *HYZ3(2,0,1)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,1,1)
     $  + HYZ4(0,1,2,1)
     $  + 2.000000000000000d+00*HYZ4(0,2,1,1)
     $  + 2.000000000000000d+00*HYZ4(2,0,1,1)
      HYZ4(1,0,1,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(0,1)
     $  - HYZ1(1) *HYZ3(0,3,1)
     $  - HYZ1(1) *HYZ3(3,0,1)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,1,1)
     $  + HYZ4(0,1,3,1)
     $  + 2.000000000000000d+00*HYZ4(0,3,1,1)
     $  + 2.000000000000000d+00*HYZ4(3,0,1,1)
      HYZ4(1,0,2,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(0,2)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(0,1)
     $  + HYZ1(0) *HYZ3(2,0,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,0,2)
     $  + HYZ2(0,1) *HYZ2(0,2)
     $  - HYZ4(0,2,0,1)
      HYZ4(1,0,2,1) =
     $  + HYZ1(1) *HYZ3(0,2,1)
     $  - HYZ4(0,1,2,1)
     $  - 2.000000000000000d+00*HYZ4(0,2,1,1)
      HYZ4(1,0,2,2) =
     $  + HYZ1(1) *HYZ3(0,2,2)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,1)
     $  + HYZ1(2) *HYZ3(2,0,1)
     $  - HYZ4(2,2,0,1)
      HYZ4(1,0,2,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(0,2)
     $  - HYZ1(1) *HYZ3(0,3,2)
     $  - HYZ1(1) *HYZ3(3,0,2)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(0,1)
     $  + HYZ1(3) *HYZ3(2,0,1)
     $  + HYZ2(0,1) *HYZ2(3,2)
     $  - HYZ4(3,2,0,1)
      HYZ4(1,0,3,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(0,3)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(0,1)
     $  + HYZ1(0) *HYZ3(3,0,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,0,3)
     $  + HYZ2(0,1) *HYZ2(0,3)
     $  - HYZ4(0,3,0,1)
      HYZ4(1,0,3,1) =
     $  + HYZ1(1) *HYZ3(0,3,1)
     $  - HYZ4(0,1,3,1)
     $  - 2.000000000000000d+00*HYZ4(0,3,1,1)
      HYZ4(1,0,3,2) =
     $  + HYZ1(1) *HYZ3(0,3,2)
     $  + HYZ1(2) *HYZ3(3,0,1)
     $  - HYZ2(0,1) *HYZ2(3,2)
     $  - HYZ4(2,3,0,1)
      HYZ4(1,0,3,3) =
     $  + HYZ1(1) *HYZ3(0,3,3)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,1)
     $  + HYZ1(3) *HYZ3(3,0,1)
     $  - HYZ4(3,3,0,1)
      HYZ4(1,1,0,0) =
     $  + 2.500000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(1)*HYZ1(1)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(0,1)
     $  + HYZ1(0) *HYZ3(0,1,1)
     $  + HYZ1(1) *HYZ3(0,0,1)
     $  - HYZ4(0,0,1,1)
      HYZ4(1,1,0,1) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,1,1)
     $  + 3.000000000000000d+00*HYZ4(0,1,1,1)
      HYZ4(1,1,0,2) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,2)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(0,1)
     $  + HYZ1(1) *HYZ3(2,0,1)
     $  + HYZ1(2) *HYZ3(0,1,1)
     $  - HYZ4(2,0,1,1)
      HYZ4(1,1,0,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,3)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(0,1)
     $  + HYZ1(1) *HYZ3(3,0,1)
     $  + HYZ1(3) *HYZ3(0,1,1)
     $  - HYZ4(3,0,1,1)
      HYZ4(1,1,1,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(1)*HYZ1(1)*HYZ1(1)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,1)
     $  + HYZ1(1) *HYZ3(0,1,1)
     $  - HYZ4(0,1,1,1)
      HYZ4(1,1,1,1) =
     $  + 4.166666666666666d-02*HYZ1(1)*HYZ1(1)*HYZ1(1)*HYZ1(1)
      HYZ4(1,1,1,2) =
     $  + 1.666666666666666d-01*HYZ1(1)*HYZ1(1)*HYZ1(1)*HYZ1(2)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(2,1)
     $  + HYZ1(1) *HYZ3(2,1,1)
     $  - HYZ4(2,1,1,1)
      HYZ4(1,1,1,3) =
     $  + 1.666666666666666d-01*HYZ1(1)*HYZ1(1)*HYZ1(1)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(3,1)
     $  + HYZ1(1) *HYZ3(3,1,1)
     $  - HYZ4(3,1,1,1)
      HYZ4(1,1,2,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(1)*HYZ1(1)*HYZ1(2)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(2,1)
     $  + HYZ1(0) *HYZ3(2,1,1)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,2)
     $  + HYZ1(1) *HYZ3(0,2,1)
     $  - HYZ4(0,2,1,1)
      HYZ4(1,1,2,1) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(2,1,1)
     $  + 3.000000000000000d+00*HYZ4(2,1,1,1)
      HYZ4(1,1,2,2) =
     $  + 2.500000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ1(2)*HYZ1(2)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(2,1)
     $  + HYZ1(1) *HYZ3(2,2,1)
     $  + HYZ1(2) *HYZ3(2,1,1)
     $  - HYZ4(2,2,1,1)
      HYZ4(1,1,2,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ1(2)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(3,2)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(2,1)
     $  + HYZ1(1) *HYZ3(3,2,1)
     $  + HYZ1(3) *HYZ3(2,1,1)
     $  - HYZ4(3,2,1,1)
      HYZ4(1,1,3,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(1)*HYZ1(1)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(3,1)
     $  + HYZ1(0) *HYZ3(3,1,1)
     $  - 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(0,3)
     $  + HYZ1(1) *HYZ3(0,3,1)
     $  - HYZ4(0,3,1,1)
      HYZ4(1,1,3,1) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(3,1,1)
     $  + 3.000000000000000d+00*HYZ4(3,1,1,1)
      HYZ4(1,1,3,2) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ2(3,2)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(3,1)
     $  + HYZ1(1) *HYZ3(2,3,1)
     $  + HYZ1(2) *HYZ3(3,1,1)
     $  - HYZ4(2,3,1,1)
      HYZ4(1,1,3,3) =
     $  + 2.500000000000000d-01*HYZ1(1)*HYZ1(1)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(3,1)
     $  + HYZ1(1) *HYZ3(3,3,1)
     $  + HYZ1(3) *HYZ3(3,1,1)
     $  - HYZ4(3,3,1,1)
      HYZ4(1,2,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(1)*HYZ1(2)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(2,1)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(0,2)
     $  + HYZ1(0) *HYZ3(0,2,1)
     $  + HYZ1(1) *HYZ3(0,0,2)
     $  - HYZ4(0,0,2,1)
      HYZ4(1,2,0,1) =
     $  + HYZ1(1) *HYZ3(2,0,1)
     $  - HYZ2(0,1) *HYZ2(2,1)
     $  + HYZ4(0,1,2,1)
     $  + 2.000000000000000d+00*HYZ4(0,2,1,1)
      HYZ4(1,2,0,2) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,2,2)
     $  + HYZ1(2) *HYZ3(0,2,1)
     $  - HYZ2(0,2) *HYZ2(2,1)
     $  - HYZ4(2,0,2,1)
      HYZ4(1,2,0,3) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(0,3)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(0,2)
     $  + HYZ1(1) *HYZ3(3,0,2)
     $  + HYZ1(3) *HYZ3(0,2,1)
     $  - HYZ2(0,3) *HYZ2(2,1)
     $  - HYZ4(3,0,2,1)
      HYZ4(1,2,1,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(2,1,1)
     $  - HYZ1(1) *HYZ3(0,2,1)
     $  - HYZ1(1) *HYZ3(2,0,1)
     $  + HYZ2(0,1) *HYZ2(2,1)
     $  - HYZ4(0,1,2,1)
      HYZ4(1,2,1,1) =
     $  + HYZ1(1) *HYZ3(2,1,1)
     $  - 3.000000000000000d+00*HYZ4(2,1,1,1)
      HYZ4(1,2,1,2) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(2,2,1)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(2,1,1)
     $  + 5.000000000000000d-01*HYZ2(2,1)*HYZ2(2,1)
     $  + 2.000000000000000d+00*HYZ4(2,2,1,1)
      HYZ4(1,2,1,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(2,1)
     $  - HYZ1(1) *HYZ3(2,3,1)
     $  - HYZ1(1) *HYZ3(3,2,1)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(2,1,1)
     $  + HYZ2(2,1) *HYZ2(3,1)
     $  - HYZ4(3,1,2,1)
      HYZ4(1,2,2,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(1)*HYZ1(2)*HYZ1(2)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(2,1)
     $  + HYZ1(0) *HYZ3(2,2,1)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(0,2)
     $  + HYZ1(1) *HYZ3(0,2,2)
     $  + HYZ2(0,2) *HYZ2(2,1)
     $  - HYZ4(0,2,2,1)
      HYZ4(1,2,2,1) =
     $  + HYZ1(1) *HYZ3(2,2,1)
     $  - 5.000000000000000d-01*HYZ2(2,1)*HYZ2(2,1)
      HYZ4(1,2,2,2) =
     $  + 1.666666666666666d-01*HYZ1(1)*HYZ1(2)*HYZ1(2)*HYZ1(2)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(2,1)
     $  + HYZ1(2) *HYZ3(2,2,1)
     $  - HYZ4(2,2,2,1)
      HYZ4(1,2,2,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(2)*HYZ1(2)*HYZ1(3)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(3,2)
     $  + HYZ1(1) *HYZ3(3,2,2)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(2,1)
     $  + HYZ1(3) *HYZ3(2,2,1)
     $  + HYZ2(2,1) *HYZ2(3,2)
     $  - HYZ4(3,2,2,1)
      HYZ4(1,2,3,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ1(2)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(3,2)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(2,1)
     $  + HYZ1(0) *HYZ3(3,2,1)
     $  - HYZ1(1) *HYZ1(2)*HYZ2(0,3)
     $  + HYZ1(1) *HYZ3(0,3,2)
     $  + HYZ2(0,3) *HYZ2(2,1)
     $  - HYZ4(0,3,2,1)
      HYZ4(1,2,3,1) =
     $  + HYZ1(1) *HYZ3(2,3,1)
     $  - HYZ2(2,1) *HYZ2(3,1)
     $  + HYZ4(3,1,2,1)
     $  + 2.000000000000000d+00*HYZ4(3,2,1,1)
      HYZ4(1,2,3,2) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(3,2,2)
     $  + HYZ1(2) *HYZ3(3,2,1)
     $  - HYZ2(2,1) *HYZ2(3,2)
     $  - HYZ4(2,3,2,1)
      HYZ4(1,2,3,3) =
     $  + 5.000000000000000d-01*HYZ1(1)*HYZ1(2)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(3,2)
     $  + HYZ1(1) *HYZ3(3,3,2)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(2,1)
     $  + HYZ1(3) *HYZ3(3,2,1)
     $  - HYZ4(3,3,2,1)
      HYZ4(1,3,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(1)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(3,1)
     $  - HYZ1(0) *HYZ1(1)*HYZ2(0,3)
     $  + HYZ1(0) *HYZ3(0,3,1)
     $  + HYZ1(1) *HYZ3(0,0,3)
     $  - HYZ4(0,0,3,1)
      HYZ4(1,3,0,1) =
     $  + HYZ1(1) *HYZ3(3,0,1)
     $  - HYZ2(0,1) *HYZ2(3,1)
     $  + HYZ4(0,1,3,1)
     $  + 2.000000000000000d+00*HYZ4(0,3,1,1)
      HYZ4(1,3,0,2) =
     $  + HYZ1(1) *HYZ3(3,0,2)
     $  + HYZ1(2) *HYZ3(0,3,1)
     $  - HYZ2(0,2) *HYZ2(3,1)
     $  - HYZ4(2,0,3,1)
      HYZ4(1,3,0,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(0,3,3)
     $  + HYZ1(3) *HYZ3(0,3,1)
     $  - HYZ2(0,3) *HYZ2(3,1)
     $  - HYZ4(3,0,3,1)
      HYZ4(1,3,1,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(3,1,1)
     $  - HYZ1(1) *HYZ3(0,3,1)
     $  - HYZ1(1) *HYZ3(3,0,1)
     $  + HYZ2(0,1) *HYZ2(3,1)
     $  - HYZ4(0,1,3,1)
      HYZ4(1,3,1,1) =
     $  + HYZ1(1) *HYZ3(3,1,1)
     $  - 3.000000000000000d+00*HYZ4(3,1,1,1)
      HYZ4(1,3,1,2) =
     $  + HYZ1(1) *HYZ1(2)*HYZ2(3,1)
     $  - HYZ1(1) *HYZ3(2,3,1)
     $  - HYZ1(1) *HYZ3(3,2,1)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(3,1,1)
     $  + 2.000000000000000d+00*HYZ4(2,3,1,1)
     $  + HYZ4(3,1,2,1)
     $  + 2.000000000000000d+00*HYZ4(3,2,1,1)
      HYZ4(1,3,1,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(3,3,1)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(3,1,1)
     $  + 5.000000000000000d-01*HYZ2(3,1)*HYZ2(3,1)
     $  + 2.000000000000000d+00*HYZ4(3,3,1,1)
      HYZ4(1,3,2,0) =
     $  + HYZ1(0) *HYZ1(1)*HYZ2(3,2)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(3,1)
     $  + HYZ1(0) *HYZ3(2,3,1)
     $  - HYZ1(1) *HYZ3(0,3,2)
     $  - HYZ1(1) *HYZ3(3,0,2)
     $  + HYZ2(0,2) *HYZ2(3,1)
     $  - HYZ4(0,2,3,1)
      HYZ4(1,3,2,1) =
     $  + HYZ1(1) *HYZ3(3,2,1)
     $  - HYZ4(3,1,2,1)
     $  - 2.000000000000000d+00*HYZ4(3,2,1,1)
      HYZ4(1,3,2,2) =
     $  + HYZ1(1) *HYZ3(3,2,2)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(3,1)
     $  + HYZ1(2) *HYZ3(2,3,1)
     $  - HYZ4(2,2,3,1)
      HYZ4(1,3,2,3) =
     $  + HYZ1(1) *HYZ1(3)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(1)*HYZ3(3,3,2)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(3,1)
     $  + HYZ1(3) *HYZ3(2,3,1)
     $  + HYZ2(3,1) *HYZ2(3,2)
     $  - HYZ4(3,2,3,1)
      HYZ4(1,3,3,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(1)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(3,1)
     $  + HYZ1(0) *HYZ3(3,3,1)
     $  - HYZ1(1) *HYZ1(3)*HYZ2(0,3)
     $  + HYZ1(1) *HYZ3(0,3,3)
     $  + HYZ2(0,3) *HYZ2(3,1)
     $  - HYZ4(0,3,3,1)
      HYZ4(1,3,3,1) =
     $  + HYZ1(1) *HYZ3(3,3,1)
     $  - 5.000000000000000d-01*HYZ2(3,1)*HYZ2(3,1)
      HYZ4(1,3,3,2) =
     $  + HYZ1(1) *HYZ3(3,3,2)
     $  + HYZ1(2) *HYZ3(3,3,1)
     $  - HYZ2(3,1) *HYZ2(3,2)
     $  - HYZ4(2,3,3,1)
      HYZ4(1,3,3,3) =
     $  + 1.666666666666666d-01*HYZ1(1)*HYZ1(3)*HYZ1(3)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(3,1)
     $  + HYZ1(3) *HYZ3(3,3,1)
     $  - HYZ4(3,3,3,1)
      HYZ4(2,0,0,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(0)*HYZ1(0)*HYZ1(2)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,2)
     $  + HYZ1(0) *HYZ3(0,0,2)
     $  - HYZ4(0,0,0,2)
      HYZ4(2,0,0,2) =
     $  + HYZ1(2) *HYZ3(0,0,2)
     $  - 5.000000000000000d-01*HYZ2(0,2)*HYZ2(0,2)
      HYZ4(2,0,0,3) =
     $  + HYZ1(2) *HYZ3(0,0,3)
     $  + HYZ1(3) *HYZ3(0,0,2)
     $  - HYZ2(0,2) *HYZ2(0,3)
     $  - HYZ4(3,0,0,2)
      HYZ4(2,0,1,0) =
     $  + HYZ1(0) *HYZ3(2,0,1)
     $  - HYZ4(0,2,0,1)
     $  - 2.000000000000000d+00*HYZ4(2,0,0,1)
      HYZ4(2,0,1,2) =
     $  + HYZ1(2) *HYZ3(2,0,1)
     $  - HYZ4(2,0,2,1)
     $  - 2.000000000000000d+00*HYZ4(2,2,0,1)
      HYZ4(2,0,1,3) =
     $  + HYZ1(3) *HYZ3(2,0,1)
     $  - HYZ4(2,0,3,1)
     $  - HYZ4(2,3,0,1)
     $  - HYZ4(3,2,0,1)
      HYZ4(2,0,2,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,2,2)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,0,2)
     $  + 5.000000000000000d-01*HYZ2(0,2)*HYZ2(0,2)
     $  + 2.000000000000000d+00*HYZ4(0,0,2,2)
      HYZ4(2,0,2,2) =
     $  + HYZ1(2) *HYZ3(0,2,2)
     $  - 3.000000000000000d+00*HYZ4(0,2,2,2)
      HYZ4(2,0,2,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(0,2)
     $  - HYZ1(2) *HYZ3(0,3,2)
     $  - HYZ1(2) *HYZ3(3,0,2)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,2,2)
     $  + HYZ4(0,2,3,2)
     $  + 2.000000000000000d+00*HYZ4(0,3,2,2)
     $  + 2.000000000000000d+00*HYZ4(3,0,2,2)
      HYZ4(2,0,3,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(0,3)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(0,2)
     $  + HYZ1(0) *HYZ3(3,0,2)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,0,3)
     $  + HYZ2(0,2) *HYZ2(0,3)
     $  - HYZ4(0,3,0,2)
      HYZ4(2,0,3,2) =
     $  + HYZ1(2) *HYZ3(0,3,2)
     $  - HYZ4(0,2,3,2)
     $  - 2.000000000000000d+00*HYZ4(0,3,2,2)
      HYZ4(2,0,3,3) =
     $  + HYZ1(2) *HYZ3(0,3,3)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,2)
     $  + HYZ1(3) *HYZ3(3,0,2)
     $  - HYZ4(3,3,0,2)
      HYZ4(2,1,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(2,1)
     $  - HYZ1(0) *HYZ3(0,2,1)
     $  - HYZ1(0) *HYZ3(2,0,1)
     $  + HYZ4(0,0,2,1)
     $  + HYZ4(0,2,0,1)
     $  + HYZ4(2,0,0,1)
      HYZ4(2,1,0,1) =
     $  + HYZ2(0,1) *HYZ2(2,1)
     $  - HYZ4(0,1,2,1)
     $  - 2.000000000000000d+00*HYZ4(0,2,1,1)
     $  - 2.000000000000000d+00*HYZ4(2,0,1,1)
      HYZ4(2,1,0,2) =
     $  - HYZ1(2) *HYZ3(0,2,1)
     $  - HYZ1(2) *HYZ3(2,0,1)
     $  + HYZ2(0,2) *HYZ2(2,1)
     $  + HYZ4(2,0,2,1)
     $  + 2.000000000000000d+00*HYZ4(2,2,0,1)
      HYZ4(2,1,0,3) =
     $  - HYZ1(3) *HYZ3(0,2,1)
     $  - HYZ1(3) *HYZ3(2,0,1)
     $  + HYZ2(0,3) *HYZ2(2,1)
     $  + HYZ4(2,3,0,1)
     $  + HYZ4(3,0,2,1)
     $  + HYZ4(3,2,0,1)
      HYZ4(2,1,1,0) =
     $  + HYZ1(0) *HYZ3(2,1,1)
     $  - HYZ2(0,1) *HYZ2(2,1)
     $  + HYZ4(0,1,2,1)
     $  + HYZ4(0,2,1,1)
     $  + HYZ4(2,0,1,1)
      HYZ4(2,1,1,2) =
     $  + HYZ1(2) *HYZ3(2,1,1)
     $  - 5.000000000000000d-01*HYZ2(2,1)*HYZ2(2,1)
      HYZ4(2,1,1,3) =
     $  + HYZ1(3) *HYZ3(2,1,1)
     $  - HYZ2(2,1) *HYZ2(3,1)
     $  + HYZ4(2,3,1,1)
     $  + HYZ4(3,1,2,1)
     $  + HYZ4(3,2,1,1)
      HYZ4(2,1,2,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(2,2,1)
     $  - HYZ2(0,2) *HYZ2(2,1)
     $  + 2.000000000000000d+00*HYZ4(0,2,2,1)
     $  + HYZ4(2,0,2,1)
      HYZ4(2,1,2,1) =
     $  + 5.000000000000000d-01*HYZ2(2,1)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ4(2,2,1,1)
      HYZ4(2,1,2,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(2,2,1)
     $  + 3.000000000000000d+00*HYZ4(2,2,2,1)
      HYZ4(2,1,2,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(2,1)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(2,2,1)
     $  - HYZ2(2,1) *HYZ2(3,2)
     $  + HYZ4(2,3,2,1)
     $  + 2.000000000000000d+00*HYZ4(3,2,2,1)
      HYZ4(2,1,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(2,1)
     $  - HYZ1(0) *HYZ3(2,3,1)
     $  - HYZ1(0) *HYZ3(3,2,1)
     $  - HYZ2(0,3) *HYZ2(2,1)
     $  + HYZ4(0,2,3,1)
     $  + HYZ4(0,3,2,1)
     $  + HYZ4(2,0,3,1)
      HYZ4(2,1,3,1) =
     $  + HYZ2(2,1) *HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ4(2,3,1,1)
     $  - HYZ4(3,1,2,1)
     $  - 2.000000000000000d+00*HYZ4(3,2,1,1)
      HYZ4(2,1,3,2) =
     $  - HYZ1(2) *HYZ3(2,3,1)
     $  - HYZ1(2) *HYZ3(3,2,1)
     $  + HYZ2(2,1) *HYZ2(3,2)
     $  + 2.000000000000000d+00*HYZ4(2,2,3,1)
     $  + HYZ4(2,3,2,1)
      HYZ4(2,1,3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(2,1)
     $  - HYZ1(3) *HYZ3(2,3,1)
     $  - HYZ1(3) *HYZ3(3,2,1)
     $  + HYZ4(2,3,3,1)
     $  + HYZ4(3,2,3,1)
     $  + HYZ4(3,3,2,1)
      HYZ4(2,2,0,0) =
     $  + 2.500000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(2)*HYZ1(2)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(0,2)
     $  + HYZ1(0) *HYZ3(0,2,2)
     $  + HYZ1(2) *HYZ3(0,0,2)
     $  - HYZ4(0,0,2,2)
      HYZ4(2,2,0,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,2)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,2,2)
     $  + 3.000000000000000d+00*HYZ4(0,2,2,2)
      HYZ4(2,2,0,3) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,3)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(0,2)
     $  + HYZ1(2) *HYZ3(3,0,2)
     $  + HYZ1(3) *HYZ3(0,2,2)
     $  - HYZ4(3,0,2,2)
      HYZ4(2,2,1,0) =
     $  + HYZ1(0) *HYZ3(2,2,1)
     $  - HYZ4(0,2,2,1)
     $  - HYZ4(2,0,2,1)
     $  - HYZ4(2,2,0,1)
      HYZ4(2,2,1,2) =
     $  + HYZ1(2) *HYZ3(2,2,1)
     $  - 3.000000000000000d+00*HYZ4(2,2,2,1)
      HYZ4(2,2,1,3) =
     $  + HYZ1(3) *HYZ3(2,2,1)
     $  - HYZ4(2,2,3,1)
     $  - HYZ4(2,3,2,1)
     $  - HYZ4(3,2,2,1)
      HYZ4(2,2,2,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(2)*HYZ1(2)*HYZ1(2)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,2)
     $  + HYZ1(2) *HYZ3(0,2,2)
     $  - HYZ4(0,2,2,2)
      HYZ4(2,2,2,2) =
     $  + 4.166666666666666d-02*HYZ1(2)*HYZ1(2)*HYZ1(2)*HYZ1(2)
      HYZ4(2,2,2,3) =
     $  + 1.666666666666666d-01*HYZ1(2)*HYZ1(2)*HYZ1(2)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(3,2)
     $  + HYZ1(2) *HYZ3(3,2,2)
     $  - HYZ4(3,2,2,2)
      HYZ4(2,2,3,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(2)*HYZ1(2)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(3,2)
     $  + HYZ1(0) *HYZ3(3,2,2)
     $  - 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(0,3)
     $  + HYZ1(2) *HYZ3(0,3,2)
     $  - HYZ4(0,3,2,2)
      HYZ4(2,2,3,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(3,2,2)
     $  + 3.000000000000000d+00*HYZ4(3,2,2,2)
      HYZ4(2,2,3,3) =
     $  + 2.500000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(3,2)
     $  + HYZ1(2) *HYZ3(3,3,2)
     $  + HYZ1(3) *HYZ3(3,2,2)
     $  - HYZ4(3,3,2,2)
      HYZ4(2,3,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(2)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(3,2)
     $  - HYZ1(0) *HYZ1(2)*HYZ2(0,3)
     $  + HYZ1(0) *HYZ3(0,3,2)
     $  + HYZ1(2) *HYZ3(0,0,3)
     $  - HYZ4(0,0,3,2)
      HYZ4(2,3,0,2) =
     $  + HYZ1(2) *HYZ3(3,0,2)
     $  - HYZ2(0,2) *HYZ2(3,2)
     $  + HYZ4(0,2,3,2)
     $  + 2.000000000000000d+00*HYZ4(0,3,2,2)
      HYZ4(2,3,0,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(0,3,3)
     $  + HYZ1(3) *HYZ3(0,3,2)
     $  - HYZ2(0,3) *HYZ2(3,2)
     $  - HYZ4(3,0,3,2)
      HYZ4(2,3,1,0) =
     $  + HYZ1(0) *HYZ3(2,3,1)
     $  - HYZ4(0,2,3,1)
     $  - HYZ4(2,0,3,1)
     $  - HYZ4(2,3,0,1)
      HYZ4(2,3,1,2) =
     $  + HYZ1(2) *HYZ3(2,3,1)
     $  - 2.000000000000000d+00*HYZ4(2,2,3,1)
     $  - HYZ4(2,3,2,1)
      HYZ4(2,3,1,3) =
     $  + HYZ1(3) *HYZ3(2,3,1)
     $  - 2.000000000000000d+00*HYZ4(2,3,3,1)
     $  - HYZ4(3,2,3,1)
      HYZ4(2,3,2,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(3,2,2)
     $  - HYZ1(2) *HYZ3(0,3,2)
     $  - HYZ1(2) *HYZ3(3,0,2)
     $  + HYZ2(0,2) *HYZ2(3,2)
     $  - HYZ4(0,2,3,2)
      HYZ4(2,3,2,2) =
     $  + HYZ1(2) *HYZ3(3,2,2)
     $  - 3.000000000000000d+00*HYZ4(3,2,2,2)
      HYZ4(2,3,2,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(3,3,2)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(3,2,2)
     $  + 5.000000000000000d-01*HYZ2(3,2)*HYZ2(3,2)
     $  + 2.000000000000000d+00*HYZ4(3,3,2,2)
      HYZ4(2,3,3,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(2)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(3,2)
     $  + HYZ1(0) *HYZ3(3,3,2)
     $  - HYZ1(2) *HYZ1(3)*HYZ2(0,3)
     $  + HYZ1(2) *HYZ3(0,3,3)
     $  + HYZ2(0,3) *HYZ2(3,2)
     $  - HYZ4(0,3,3,2)
      HYZ4(2,3,3,2) =
     $  + HYZ1(2) *HYZ3(3,3,2)
     $  - 5.000000000000000d-01*HYZ2(3,2)*HYZ2(3,2)
      HYZ4(2,3,3,3) =
     $  + 1.666666666666666d-01*HYZ1(2)*HYZ1(3)*HYZ1(3)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(3,2)
     $  + HYZ1(3) *HYZ3(3,3,2)
     $  - HYZ4(3,3,3,2)
      HYZ4(3,0,0,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(0)*HYZ1(0)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(0,3)
     $  + HYZ1(0) *HYZ3(0,0,3)
     $  - HYZ4(0,0,0,3)
      HYZ4(3,0,0,3) =
     $  + HYZ1(3) *HYZ3(0,0,3)
     $  - 5.000000000000000d-01*HYZ2(0,3)*HYZ2(0,3)
      HYZ4(3,0,1,0) =
     $  + HYZ1(0) *HYZ3(3,0,1)
     $  - HYZ4(0,3,0,1)
     $  - 2.000000000000000d+00*HYZ4(3,0,0,1)
      HYZ4(3,0,1,2) =
     $  + HYZ1(2) *HYZ3(3,0,1)
     $  - HYZ4(2,3,0,1)
     $  - HYZ4(3,0,2,1)
     $  - HYZ4(3,2,0,1)
      HYZ4(3,0,1,3) =
     $  + HYZ1(3) *HYZ3(3,0,1)
     $  - HYZ4(3,0,3,1)
     $  - 2.000000000000000d+00*HYZ4(3,3,0,1)
      HYZ4(3,0,2,0) =
     $  + HYZ1(0) *HYZ3(3,0,2)
     $  - HYZ4(0,3,0,2)
     $  - 2.000000000000000d+00*HYZ4(3,0,0,2)
      HYZ4(3,0,2,3) =
     $  + HYZ1(3) *HYZ3(3,0,2)
     $  - HYZ4(3,0,3,2)
     $  - 2.000000000000000d+00*HYZ4(3,3,0,2)
      HYZ4(3,0,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(0,3,3)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,0,3)
     $  + 5.000000000000000d-01*HYZ2(0,3)*HYZ2(0,3)
     $  + 2.000000000000000d+00*HYZ4(0,0,3,3)
      HYZ4(3,0,3,3) =
     $  + HYZ1(3) *HYZ3(0,3,3)
     $  - 3.000000000000000d+00*HYZ4(0,3,3,3)
      HYZ4(3,1,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(3,1)
     $  - HYZ1(0) *HYZ3(0,3,1)
     $  - HYZ1(0) *HYZ3(3,0,1)
     $  + HYZ4(0,0,3,1)
     $  + HYZ4(0,3,0,1)
     $  + HYZ4(3,0,0,1)
      HYZ4(3,1,0,1) =
     $  + HYZ2(0,1) *HYZ2(3,1)
     $  - HYZ4(0,1,3,1)
     $  - 2.000000000000000d+00*HYZ4(0,3,1,1)
     $  - 2.000000000000000d+00*HYZ4(3,0,1,1)
      HYZ4(3,1,0,2) =
     $  - HYZ1(2) *HYZ3(0,3,1)
     $  - HYZ1(2) *HYZ3(3,0,1)
     $  + HYZ2(0,2) *HYZ2(3,1)
     $  + HYZ4(2,0,3,1)
     $  + HYZ4(2,3,0,1)
     $  + HYZ4(3,2,0,1)
      HYZ4(3,1,0,3) =
     $  - HYZ1(3) *HYZ3(0,3,1)
     $  - HYZ1(3) *HYZ3(3,0,1)
     $  + HYZ2(0,3) *HYZ2(3,1)
     $  + HYZ4(3,0,3,1)
     $  + 2.000000000000000d+00*HYZ4(3,3,0,1)
      HYZ4(3,1,1,0) =
     $  + HYZ1(0) *HYZ3(3,1,1)
     $  - HYZ2(0,1) *HYZ2(3,1)
     $  + HYZ4(0,1,3,1)
     $  + HYZ4(0,3,1,1)
     $  + HYZ4(3,0,1,1)
      HYZ4(3,1,1,2) =
     $  + HYZ1(2) *HYZ3(3,1,1)
     $  - HYZ4(2,3,1,1)
     $  - HYZ4(3,1,2,1)
     $  - HYZ4(3,2,1,1)
      HYZ4(3,1,1,3) =
     $  + HYZ1(3) *HYZ3(3,1,1)
     $  - 5.000000000000000d-01*HYZ2(3,1)*HYZ2(3,1)
      HYZ4(3,1,2,0) =
     $  + HYZ1(0) *HYZ1(2)*HYZ2(3,1)
     $  - HYZ1(0) *HYZ3(2,3,1)
     $  - HYZ1(0) *HYZ3(3,2,1)
     $  - HYZ2(0,2) *HYZ2(3,1)
     $  + HYZ4(0,2,3,1)
     $  + HYZ4(0,3,2,1)
     $  + HYZ4(3,0,2,1)
      HYZ4(3,1,2,2) =
     $  + 5.000000000000000d-01*HYZ1(2)*HYZ1(2)*HYZ2(3,1)
     $  - HYZ1(2) *HYZ3(2,3,1)
     $  - HYZ1(2) *HYZ3(3,2,1)
     $  + HYZ4(2,2,3,1)
     $  + HYZ4(2,3,2,1)
     $  + HYZ4(3,2,2,1)
      HYZ4(3,1,2,3) =
     $  + HYZ1(2) *HYZ1(3)*HYZ2(3,1)
     $  - HYZ1(3) *HYZ3(2,3,1)
     $  - HYZ1(3) *HYZ3(3,2,1)
     $  - HYZ2(3,1) *HYZ2(3,2)
     $  + HYZ4(3,2,3,1)
     $  + 2.000000000000000d+00*HYZ4(3,3,2,1)
      HYZ4(3,1,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(3,3,1)
     $  - HYZ2(0,3) *HYZ2(3,1)
     $  + 2.000000000000000d+00*HYZ4(0,3,3,1)
     $  + HYZ4(3,0,3,1)
      HYZ4(3,1,3,1) =
     $  + 5.000000000000000d-01*HYZ2(3,1)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ4(3,3,1,1)
      HYZ4(3,1,3,2) =
     $  - 2.000000000000000d+00*HYZ1(2)*HYZ3(3,3,1)
     $  + HYZ2(3,1) *HYZ2(3,2)
     $  + 2.000000000000000d+00*HYZ4(2,3,3,1)
     $  + HYZ4(3,2,3,1)
      HYZ4(3,1,3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(3,1)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(3,3,1)
     $  + 3.000000000000000d+00*HYZ4(3,3,3,1)
      HYZ4(3,2,0,0) =
     $  + 5.000000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ2(3,2)
     $  - HYZ1(0) *HYZ3(0,3,2)
     $  - HYZ1(0) *HYZ3(3,0,2)
     $  + HYZ4(0,0,3,2)
     $  + HYZ4(0,3,0,2)
     $  + HYZ4(3,0,0,2)
      HYZ4(3,2,0,2) =
     $  + HYZ2(0,2) *HYZ2(3,2)
     $  - HYZ4(0,2,3,2)
     $  - 2.000000000000000d+00*HYZ4(0,3,2,2)
     $  - 2.000000000000000d+00*HYZ4(3,0,2,2)
      HYZ4(3,2,0,3) =
     $  - HYZ1(3) *HYZ3(0,3,2)
     $  - HYZ1(3) *HYZ3(3,0,2)
     $  + HYZ2(0,3) *HYZ2(3,2)
     $  + HYZ4(3,0,3,2)
     $  + 2.000000000000000d+00*HYZ4(3,3,0,2)
      HYZ4(3,2,1,0) =
     $  + HYZ1(0) *HYZ3(3,2,1)
     $  - HYZ4(0,3,2,1)
     $  - HYZ4(3,0,2,1)
     $  - HYZ4(3,2,0,1)
      HYZ4(3,2,1,2) =
     $  + HYZ1(2) *HYZ3(3,2,1)
     $  - HYZ4(2,3,2,1)
     $  - 2.000000000000000d+00*HYZ4(3,2,2,1)
      HYZ4(3,2,1,3) =
     $  + HYZ1(3) *HYZ3(3,2,1)
     $  - HYZ4(3,2,3,1)
     $  - 2.000000000000000d+00*HYZ4(3,3,2,1)
      HYZ4(3,2,2,0) =
     $  + HYZ1(0) *HYZ3(3,2,2)
     $  - HYZ2(0,2) *HYZ2(3,2)
     $  + HYZ4(0,2,3,2)
     $  + HYZ4(0,3,2,2)
     $  + HYZ4(3,0,2,2)
      HYZ4(3,2,2,3) =
     $  + HYZ1(3) *HYZ3(3,2,2)
     $  - 5.000000000000000d-01*HYZ2(3,2)*HYZ2(3,2)
      HYZ4(3,2,3,0) =
     $  + HYZ1(0) *HYZ1(3)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(0)*HYZ3(3,3,2)
     $  - HYZ2(0,3) *HYZ2(3,2)
     $  + 2.000000000000000d+00*HYZ4(0,3,3,2)
     $  + HYZ4(3,0,3,2)
      HYZ4(3,2,3,2) =
     $  + 5.000000000000000d-01*HYZ2(3,2)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ4(3,3,2,2)
      HYZ4(3,2,3,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(3,2)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(3,3,2)
     $  + 3.000000000000000d+00*HYZ4(3,3,3,2)
      HYZ4(3,3,0,0) =
     $  + 2.500000000000000d-01*HYZ1(0)*HYZ1(0)*HYZ1(3)*HYZ1(3)
     $  - HYZ1(0) *HYZ1(3)*HYZ2(0,3)
     $  + HYZ1(0) *HYZ3(0,3,3)
     $  + HYZ1(3) *HYZ3(0,0,3)
     $  - HYZ4(0,0,3,3)
      HYZ4(3,3,0,3) =
     $  + 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,3)
     $  - 2.000000000000000d+00*HYZ1(3)*HYZ3(0,3,3)
     $  + 3.000000000000000d+00*HYZ4(0,3,3,3)
      HYZ4(3,3,1,0) =
     $  + HYZ1(0) *HYZ3(3,3,1)
     $  - HYZ4(0,3,3,1)
     $  - HYZ4(3,0,3,1)
     $  - HYZ4(3,3,0,1)
      HYZ4(3,3,1,2) =
     $  + HYZ1(2) *HYZ3(3,3,1)
     $  - HYZ4(2,3,3,1)
     $  - HYZ4(3,2,3,1)
     $  - HYZ4(3,3,2,1)
      HYZ4(3,3,1,3) =
     $  + HYZ1(3) *HYZ3(3,3,1)
     $  - 3.000000000000000d+00*HYZ4(3,3,3,1)
      HYZ4(3,3,2,0) =
     $  + HYZ1(0) *HYZ3(3,3,2)
     $  - HYZ4(0,3,3,2)
     $  - HYZ4(3,0,3,2)
     $  - HYZ4(3,3,0,2)
      HYZ4(3,3,2,3) =
     $  + HYZ1(3) *HYZ3(3,3,2)
     $  - 3.000000000000000d+00*HYZ4(3,3,3,2)
      HYZ4(3,3,3,0) =
     $  + 1.666666666666666d-01*HYZ1(0)*HYZ1(3)*HYZ1(3)*HYZ1(3)
     $  - 5.000000000000000d-01*HYZ1(3)*HYZ1(3)*HYZ2(0,3)
     $  + HYZ1(3) *HYZ3(0,3,3)
     $  - HYZ4(0,3,3,3)
      HYZ4(3,3,3,3) =
     $  + 4.166666666666666d-02*HYZ1(3)*HYZ1(3)*HYZ1(3)*HYZ1(3)
      endif

      return
      end

      subroutine swap2dhplxy(iflag,nmax,HXZ1,HXZ2,HXZ3,HXZ4,
     $                       HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4)
*********************************************************************
*****  swap2dhpl evaluates H(..;y) in terms of H(..;x=1-y-z)   ******
*****            and H(..;z)                                   ******
***** applied for y>(1-z)/2                                    ******
***** Input: HXZn,HZn; Output: HYZn                            ******
*********************************************************************
      implicit none
      include 'types.f'
      integer iflag,nmax
      real(dp):: Zeta2,Zeta3,Zeta4
      real(dp):: HZ1,HZ2,HZ3,HZ4,HYZ1,HYZ2,HYZ3,HYZ4
      real(dp):: HXZ1,HXZ2,HXZ3,HXZ4
      dimension HYZ1(0:3),HYZ2(0:3,0:3),HYZ3(0:3,0:3,0:3),
     $          HYZ4(0:3,0:3,0:3,0:3)
      dimension HXZ1(0:3),HXZ2(0:3,0:3),HXZ3(0:3,0:3,0:3),
     $          HXZ4(0:3,0:3,0:3,0:3)
      dimension HZ1(-1:1),HZ2(-1:1,-1:1),HZ3(-1:1,-1:1,-1:1),
     $          HZ4(-1:1,-1:1,-1:1,-1:1)
      parameter (Zeta2 = 1.6449340668482264365d0)
      parameter (Zeta3 = 1.2020569031595942854d0)
      parameter (Zeta4 = 1.0823232337111381915d0)

* 2001-04-23:17:22:16.hpl
* <- swap.out
* produced by form-to-fortr for gehrt@pcth62

      HYZ1(0) = -HXZ1(2)-HZ1(1)
      HYZ1(1) = -HXZ1(3)-HZ1(0)
      HYZ1(2) = -HXZ1(0)-HZ1(1)
      HYZ1(3) = -HXZ1(1)-HZ1(0)

      if (nmax.eq.1) return
      HYZ2(0,0) =
     $  + HXZ1(2) *HZ1(1)
     $  + HXZ2(2,2)
     $  + HZ2(1,1)
      HYZ2(0,1) =
     $  + Zeta2
     $  + HXZ1(2) *HZ1(0)
     $  + HXZ2(2,3)
     $  + HZ2(1,0)
      HYZ2(0,2) =
     $  + Zeta2
     $  + HXZ1(2) *HZ1(1)
     $  + HXZ2(2,0)
      HYZ2(0,3) =
     $  + Zeta2
     $  + HXZ1(2) *HZ1(0)
     $  + HXZ2(2,1)
     $  + HZ2(0,0)
     $  + HZ2(1,0)
      HYZ2(1,0) =
     $  - Zeta2
     $  + HXZ1(3) *HZ1(1)
     $  + HXZ2(3,2)
     $  + HZ2(0,1)
      HYZ2(1,1) =
     $  + HXZ1(3) *HZ1(0)
     $  + HXZ2(3,3)
     $  + HZ2(0,0)
      HYZ2(1,2) =
     $  + Zeta2
     $  + HXZ1(3) *HZ1(1)
     $  + HXZ2(3,0)
     $  + HZ2(0,0)
     $  + HZ2(1,0)
      HYZ2(1,3) =
     $  - Zeta2
     $  + HXZ1(3) *HZ1(0)
     $  + HXZ2(3,1)
     $  - 2.000000000000000d+00*HZ2(-1,0)
     $  + 2.000000000000000d+00*HZ2(0,0)
      HYZ2(2,0) =
     $  - Zeta2
     $  + HXZ1(0) *HZ1(1)
     $  + HXZ2(0,2)
     $  + 2.000000000000000d+00*HZ2(1,1)
      HYZ2(2,1) =
     $  - Zeta2
     $  + HXZ1(0) *HZ1(0)
     $  + HXZ2(0,3)
     $  - HZ2(0,0)
     $  + HZ2(0,1)
      HYZ2(2,2) =
     $  + HXZ1(0) *HZ1(1)
     $  + HXZ2(0,0)
     $  + HZ2(1,1)
      HYZ2(2,3) =
     $  - Zeta2
     $  + HXZ1(0) *HZ1(0)
     $  + HXZ2(0,1)
     $  + HZ2(0,1)
      HYZ2(3,0) =
     $  - Zeta2
     $  + HXZ1(1) *HZ1(1)
     $  + HXZ2(1,2)
     $  - HZ2(0,0)
     $  + HZ2(0,1)
      HYZ2(3,1) =
     $  + Zeta2
     $  + HXZ1(1) *HZ1(0)
     $  + HXZ2(1,3)
     $  + 2.000000000000000d+00*HZ2(-1,0)
      HYZ2(3,2) =
     $  + Zeta2
     $  + HXZ1(1) *HZ1(1)
     $  + HXZ2(1,0)
     $  + HZ2(1,0)
      HYZ2(3,3) =
     $  + HXZ1(1) *HZ1(0)
     $  + HXZ2(1,1)
     $  + HZ2(0,0)
      if (nmax.eq.2) return
      HYZ3(0,0,0) =
     $  - HXZ1(2) *HZ2(1,1)
     $  - HXZ2(2,2) *HZ1(1)
     $  - HXZ3(2,2,2)
     $  - HZ3(1,1,1)
      HYZ3(0,0,1) =
     $  + Zeta3
     $  - HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(1,0)
     $  - HXZ2(2,2) *HZ1(0)
     $  - HXZ3(2,2,3)
     $  - HZ1(1) *Zeta2
     $  - HZ3(1,1,0)
      HYZ3(0,0,2) =
     $  + Zeta3
     $  - HXZ1(2) *Zeta2
     $  - HXZ2(2,2) *HZ1(1)
     $  - HXZ3(2,2,0)
      HYZ3(0,0,3) =
     $  - HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(0,0)
     $  - HXZ1(2) *HZ2(1,0)
     $  - HXZ2(2,2) *HZ1(0)
     $  - HXZ3(2,2,1)
     $  - HZ1(0) *Zeta2
     $  - HZ1(1) *Zeta2
     $  - HZ3(0,0,0)
     $  - HZ3(0,1,0)
     $  - HZ3(1,0,0)
     $  - HZ3(1,1,0)
      HYZ3(0,1,0) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(0,1)
     $  - HXZ2(2,3) *HZ1(1)
     $  - HXZ3(2,3,2)
     $  + HZ1(1) *Zeta2
     $  - HZ3(1,0,1)
      HYZ3(0,1,1) =
     $  + Zeta3
     $  - HXZ1(2) *HZ2(0,0)
     $  - HXZ2(2,3) *HZ1(0)
     $  - HXZ3(2,3,3)
     $  - HZ3(1,0,0)
      HYZ3(0,1,2) =
     $  + Zeta3
     $  - HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(0,0)
     $  - HXZ1(2) *HZ2(1,0)
     $  - HXZ2(2,3) *HZ1(1)
     $  - HXZ3(2,3,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta2
     $  - HZ3(0,1,0)
     $  - HZ3(1,0,0)
     $  - 2.000000000000000d+00*HZ3(1,1,0)
      HYZ3(0,1,3) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(2) *Zeta2
     $  + 2.000000000000000d+00*HXZ1(2)*HZ2(-1,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ2(0,0)
     $  - HXZ2(2,3) *HZ1(0)
     $  - HXZ3(2,3,1)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(0) *Zeta2
     $  + HZ1(1) *Zeta2
     $  + HZ3( -1,0,0)
     $  + 2.000000000000000d+00*HZ3(-1,1,0)
     $  - HZ3(0,1,0)
     $  + 2.000000000000000d+00*HZ3(1,-1,0)
     $  - 2.000000000000000d+00*HZ3(1,0,0)
      HYZ3(0,2,0) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(2) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ2(1,1)
     $  - HXZ2(2,0) *HZ1(1)
     $  - HXZ3(2,0,2)
     $  - HZ1(1) *Zeta2
      HYZ3(0,2,1) =
     $  + Zeta3
     $  + HXZ1(2) *Zeta2
     $  + HXZ1(2) *HZ2(0,0)
     $  - HXZ1(2) *HZ2(0,1)
     $  - HXZ2(2,0) *HZ1(0)
     $  - HXZ3(2,0,3)
     $  + HZ1(1) *Zeta2
     $  + HZ3(0,1,0)
     $  + HZ3(1,1,0)
      HYZ3(0,2,2) =
     $  + Zeta3
     $  - HXZ1(2) *HZ2(1,1)
     $  - HXZ2(2,0) *HZ1(1)
     $  - HXZ3(2,0,0)
      HYZ3(0,2,3) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(0,1)
     $  - HXZ2(2,0) *HZ1(0)
     $  - HXZ3(2,0,1)
     $  - HZ1(0) *Zeta2
     $  + HZ1(1) *Zeta2
     $  + HZ3(1,0,0)
     $  + HZ3(1,1,0)
      HYZ3(0,3,0) =
     $  + HXZ1(2) *Zeta2
     $  + HXZ1(2) *HZ2(0,0)
     $  - HXZ1(2) *HZ2(0,1)
     $  - HXZ2(2,1) *HZ1(1)
     $  - HXZ3(2,1,2)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta2
     $  + HZ1(1) *Zeta2
     $  + 2.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,0,1)
     $  + HZ3(0,1,0)
     $  + HZ3(1,0,0)
     $  - HZ3(1,0,1)
      HYZ3(0,3,1) =
     $  + Zeta3
     $  - HXZ1(2) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ2(-1,0)
     $  - HXZ2(2,1) *HZ1(0)
     $  - HXZ3(2,1,3)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(1) *Zeta2
     $  - HZ3( -1,0,0)
     $  - 2.000000000000000d+00*HZ3(-1,1,0)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
     $  + HZ3(0,1,0)
     $  - 2.000000000000000d+00*HZ3(1,-1,0)
      HYZ3(0,3,2) =
     $  + Zeta3
     $  - HXZ1(2) *Zeta2
     $  - HXZ1(2) *HZ2(1,0)
     $  - HXZ2(2,1) *HZ1(1)
     $  - HXZ3(2,1,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta2
     $  - HZ3(0,1,0)
     $  - HZ3(1,0,0)
     $  - 2.000000000000000d+00*HZ3(1,1,0)
      HYZ3(0,3,3) =
     $  + Zeta3
     $  - HXZ1(2) *HZ2(0,0)
     $  - HXZ2(2,1) *HZ1(0)
     $  - HXZ3(2,1,1)
     $  - HZ3(0,0,0)
     $  - HZ3(1,0,0)
      HYZ3(1,0,0) =
     $  + Zeta3
     $  - HXZ1(3) *HZ2(1,1)
     $  - HXZ2(3,2) *HZ1(1)
     $  - HXZ3(3,2,2)
     $  - HZ3(0,1,1)
      HYZ3(1,0,1) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(1,0)
     $  - HXZ2(3,2) *HZ1(0)
     $  - HXZ3(3,2,3)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,1,0)
      HYZ3(1,0,2) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(3) *Zeta2
     $  - HXZ2(3,2) *HZ1(1)
     $  - HXZ3(3,2,0)
     $  - HZ1(0) *Zeta2
     $  + HZ1(1) *Zeta2
     $  + HZ3(1,0,0)
     $  + HZ3(1,1,0)
      HYZ3(1,0,3) =
     $  + Zeta3
     $  - HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(0,0)
     $  - HXZ1(3) *HZ2(1,0)
     $  - HXZ2(3,2) *HZ1(0)
     $  - HXZ3(3,2,1)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,1,0)
      HYZ3(1,1,0) =
     $  + Zeta3
     $  + HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(0,1)
     $  - HXZ2(3,3) *HZ1(1)
     $  - HXZ3(3,3,2)
     $  + HZ1(0) *Zeta2
     $  - HZ3(0,0,1)
      HYZ3(1,1,1) =
     $  - HXZ1(3) *HZ2(0,0)
     $  - HXZ2(3,3) *HZ1(0)
     $  - HXZ3(3,3,3)
     $  - HZ3(0,0,0)
      HYZ3(1,1,2) =
     $  + Zeta3
     $  - HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(0,0)
     $  - HXZ1(3) *HZ2(1,0)
     $  - HXZ2(3,3) *HZ1(1)
     $  - HXZ3(3,3,0)
     $  - HZ3(0,0,0)
     $  - HZ3(1,0,0)
      HYZ3(1,1,3) =
     $  + Zeta3
     $  + HXZ1(3) *Zeta2
     $  + 2.000000000000000d+00*HXZ1(3)*HZ2(-1,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ2(0,0)
     $  - HXZ2(3,3) *HZ1(0)
     $  - HXZ3(3,3,1)
     $  - HZ1( -1)*Zeta2
     $  + HZ1(0) *Zeta2
     $  - 2.000000000000000d+00*HZ3(-1,-1,0)
     $  + 3.000000000000000d+00*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HZ3(0,0,0)
      HYZ3(1,2,0) =
     $  + Zeta3
     $  + HXZ1(3) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ2(1,1)
     $  - HXZ2(3,0) *HZ1(1)
     $  - HXZ3(3,0,2)
     $  + HZ1(0) *Zeta2
     $  - HZ3(0,0,1)
     $  - HZ3(1,0,0)
     $  - HZ3(1,0,1)
     $  - HZ3(1,1,0)
      HYZ3(1,2,1) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(3) *Zeta2
     $  + HXZ1(3) *HZ2(0,0)
     $  - HXZ1(3) *HZ2(0,1)
     $  - HXZ2(3,0) *HZ1(0)
     $  - HXZ3(3,0,3)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,0,0)
     $  - HZ3(0,1,0)
      HYZ3(1,2,2) =
     $  - HXZ1(3) *HZ2(1,1)
     $  - HXZ2(3,0) *HZ1(1)
     $  - HXZ3(3,0,0)
     $  - HZ1(0) *Zeta2
     $  - HZ1(1) *Zeta2
     $  - HZ3(0,0,0)
     $  - HZ3(0,1,0)
     $  - HZ3(1,0,0)
     $  - HZ3(1,1,0)
      HYZ3(1,2,3) =
     $  + Zeta3
     $  + HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(0,1)
     $  - HXZ2(3,0) *HZ1(0)
     $  - HXZ3(3,0,1)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,1,0)
      HYZ3(1,3,0) =
     $  + Zeta3
     $  + HXZ1(3) *Zeta2
     $  + HXZ1(3) *HZ2(0,0)
     $  - HXZ1(3) *HZ2(0,1)
     $  - HXZ2(3,1) *HZ1(1)
     $  - HXZ3(3,1,2)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  + HZ1(0) *Zeta2
     $  - HZ3( -1,0,0)
     $  + 2.000000000000000d+00*HZ3(-1,0,1)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HZ3(0,0,0)
     $  - 2.000000000000000d+00*HZ3(0,0,1)
      HYZ3(1,3,1) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(3) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ2(-1,0)
     $  - HXZ2(3,1) *HZ1(0)
     $  - HXZ3(3,1,3)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(0) *Zeta2
     $  + 4.000000000000000d+00*HZ3(-1,-1,0)
     $  - 2.000000000000000d+00*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
      HYZ3(1,3,2) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(3) *Zeta2
     $  - HXZ1(3) *HZ2(1,0)
     $  - HXZ2(3,1) *HZ1(1)
     $  - HXZ3(3,1,0)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(0) *Zeta2
     $  + HZ1(1) *Zeta2
     $  + HZ3( -1,0,0)
     $  + 2.000000000000000d+00*HZ3(-1,1,0)
     $  - HZ3(0,1,0)
     $  + 2.000000000000000d+00*HZ3(1,-1,0)
     $  - 2.000000000000000d+00*HZ3(1,0,0)
      HYZ3(1,3,3) =
     $  + Zeta3
     $  - HXZ1(3) *HZ2(0,0)
     $  - HXZ2(3,1) *HZ1(0)
     $  - HXZ3(3,1,1)
     $  - HZ1( -1)*Zeta2
     $  + HZ1(0) *Zeta2
     $  - 2.000000000000000d+00*HZ3(-1,-1,0)
     $  + 3.000000000000000d+00*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HZ3(0,0,0)
      HYZ3(2,0,0) =
     $  + Zeta3
     $  - HXZ1(0) *HZ2(1,1)
     $  - HXZ2(0,2) *HZ1(1)
     $  - HXZ3(0,2,2)
     $  + HZ1(1) *Zeta2
     $  - 3.000000000000000d+00*HZ3(1,1,1)
      HYZ3(2,0,1) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(1,0)
     $  - HXZ2(0,2) *HZ1(0)
     $  - HXZ3(0,2,3)
     $  + HZ3(1,0,0)
     $  - HZ3(1,0,1)
     $  - HZ3(1,1,0)
      HYZ3(2,0,2) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(0) *Zeta2
     $  - HXZ2(0,2) *HZ1(1)
     $  - HXZ3(0,2,0)
     $  - HZ1(1) *Zeta2
      HYZ3(2,0,3) =
     $  + Zeta3
     $  - HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(0,0)
     $  - HXZ1(0) *HZ2(1,0)
     $  - HXZ2(0,2) *HZ1(0)
     $  - HXZ3(0,2,1)
     $  + HZ1(0) *Zeta2
     $  - HZ3(0,0,1)
     $  - HZ3(1,0,0)
     $  - HZ3(1,0,1)
     $  - HZ3(1,1,0)
      HYZ3(2,1,0) =
     $  + Zeta3
     $  + HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(0,1)
     $  - HXZ2(0,3) *HZ1(1)
     $  - HXZ3(0,3,2)
     $  + HZ3(0,0,1)
     $  - 2.000000000000000d+00*HZ3(0,1,1)
      HYZ3(2,1,1) =
     $  + Zeta3
     $  - HXZ1(0) *HZ2(0,0)
     $  - HXZ2(0,3) *HZ1(0)
     $  - HXZ3(0,3,3)
     $  + HZ1(0) *Zeta2
     $  + 2.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,0,1)
      HYZ3(2,1,2) =
     $  - HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(0,0)
     $  - HXZ1(0) *HZ2(1,0)
     $  - HXZ2(0,3) *HZ1(1)
     $  - HXZ3(0,3,0)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta2
     $  + HZ1(1) *Zeta2
     $  + 2.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,0,1)
     $  + HZ3(0,1,0)
     $  + HZ3(1,0,0)
     $  - HZ3(1,0,1)
      HYZ3(2,1,3) =
     $  + Zeta3
     $  + HXZ1(0) *Zeta2
     $  + 2.000000000000000d+00*HXZ1(0)*HZ2(-1,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ2(0,0)
     $  - HXZ2(0,3) *HZ1(0)
     $  - HXZ3(0,3,1)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  + HZ1(0) *Zeta2
     $  - HZ3( -1,0,0)
     $  + 2.000000000000000d+00*HZ3(-1,0,1)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HZ3(0,0,0)
     $  - 2.000000000000000d+00*HZ3(0,0,1)
      HYZ3(2,2,0) =
     $  + Zeta3
     $  + HXZ1(0) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ2(1,1)
     $  - HXZ2(0,0) *HZ1(1)
     $  - HXZ3(0,0,2)
     $  + HZ1(1) *Zeta2
     $  - 3.000000000000000d+00*HZ3(1,1,1)
      HYZ3(2,2,1) =
     $  + HXZ1(0) *Zeta2
     $  + HXZ1(0) *HZ2(0,0)
     $  - HXZ1(0) *HZ2(0,1)
     $  - HXZ2(0,0) *HZ1(0)
     $  - HXZ3(0,0,3)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,0,0)
     $  + HZ3(0,0,1)
     $  - HZ3(0,1,1)
      HYZ3(2,2,2) =
     $  - HXZ1(0) *HZ2(1,1)
     $  - HXZ2(0,0) *HZ1(1)
     $  - HXZ3(0,0,0)
     $  - HZ3(1,1,1)
      HYZ3(2,2,3) =
     $  + Zeta3
     $  + HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(0,1)
     $  - HXZ2(0,0) *HZ1(0)
     $  - HXZ3(0,0,1)
     $  - HZ3(0,1,1)
      HYZ3(2,3,0) =
     $  + Zeta3
     $  + HXZ1(0) *Zeta2
     $  + HXZ1(0) *HZ2(0,0)
     $  - HXZ1(0) *HZ2(0,1)
     $  - HXZ2(0,1) *HZ1(1)
     $  - HXZ3(0,1,2)
     $  + HZ3(0,0,1)
     $  - 2.000000000000000d+00*HZ3(0,1,1)
      HYZ3(2,3,1) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(0) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ2(-1,0)
     $  - HXZ2(0,1) *HZ1(0)
     $  - HXZ3(0,1,3)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  + HZ3( -1,0,0)
     $  - 2.000000000000000d+00*HZ3(-1,0,1)
      HYZ3(2,3,2) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(0) *Zeta2
     $  - HXZ1(0) *HZ2(1,0)
     $  - HXZ2(0,1) *HZ1(1)
     $  - HXZ3(0,1,0)
     $  + HZ1(1) *Zeta2
     $  - HZ3(1,0,1)
      HYZ3(2,3,3) =
     $  + Zeta3
     $  - HXZ1(0) *HZ2(0,0)
     $  - HXZ2(0,1) *HZ1(0)
     $  - HXZ3(0,1,1)
     $  + HZ1(0) *Zeta2
     $  - HZ3(0,0,1)
      HYZ3(3,0,0) =
     $  - HXZ1(1) *HZ2(1,1)
     $  - HXZ2(1,2) *HZ1(1)
     $  - HXZ3(1,2,2)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,0,0)
     $  + HZ3(0,0,1)
     $  - HZ3(0,1,1)
      HYZ3(3,0,1) =
     $  + Zeta3
     $  - HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(1,0)
     $  - HXZ2(1,2) *HZ1(0)
     $  - HXZ3(1,2,3)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - HZ3(0,1,0)
      HYZ3(3,0,2) =
     $  + Zeta3
     $  - HXZ1(1) *Zeta2
     $  - HXZ2(1,2) *HZ1(1)
     $  - HXZ3(1,2,0)
     $  + HZ1(1) *Zeta2
     $  + HZ3(0,1,0)
     $  + HZ3(1,1,0)
      HYZ3(3,0,3) =
     $  - 2.000000000000000d+00*Zeta3
     $  - HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(0,0)
     $  - HXZ1(1) *HZ2(1,0)
     $  - HXZ2(1,2) *HZ1(0)
     $  - HXZ3(1,2,1)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,0,0)
     $  - HZ3(0,1,0)
      HYZ3(3,1,0) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(0,1)
     $  - HXZ2(1,3) *HZ1(1)
     $  - HXZ3(1,3,2)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  + HZ3( -1,0,0)
     $  - 2.000000000000000d+00*HZ3(-1,0,1)
      HYZ3(3,1,1) =
     $  + Zeta3
     $  - HXZ1(1) *HZ2(0,0)
     $  - HXZ2(1,3) *HZ1(0)
     $  - HXZ3(1,3,3)
     $  - HZ1( -1)*Zeta2
     $  - 2.000000000000000d+00*HZ3(-1,-1,0)
     $  - HZ3( -1,0,0)
      HYZ3(3,1,2) =
     $  + Zeta3
     $  - HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(0,0)
     $  - HXZ1(1) *HZ2(1,0)
     $  - HXZ2(1,3) *HZ1(1)
     $  - HXZ3(1,3,0)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(1) *Zeta2
     $  - HZ3( -1,0,0)
     $  - 2.000000000000000d+00*HZ3(-1,1,0)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
     $  + HZ3(0,1,0)
     $  - 2.000000000000000d+00*HZ3(1,-1,0)
      HYZ3(3,1,3) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(1) *Zeta2
     $  + 2.000000000000000d+00*HXZ1(1)*HZ2(-1,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ2(0,0)
     $  - HXZ2(1,3) *HZ1(0)
     $  - HXZ3(1,3,1)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta2
     $  - HZ1(0) *Zeta2
     $  + 4.000000000000000d+00*HZ3(-1,-1,0)
     $  - 2.000000000000000d+00*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HZ3(0,-1,0)
      HYZ3(3,2,0) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(1) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ2(1,1)
     $  - HXZ2(1,0) *HZ1(1)
     $  - HXZ3(1,0,2)
     $  + HZ3(1,0,0)
     $  - HZ3(1,0,1)
     $  - HZ3(1,1,0)
      HYZ3(3,2,1) =
     $  + Zeta3
     $  + HXZ1(1) *Zeta2
     $  + HXZ1(1) *HZ2(0,0)
     $  - HXZ1(1) *HZ2(0,1)
     $  - HXZ2(1,0) *HZ1(0)
     $  - HXZ3(1,0,3)
     $  + 2.000000000000000d+00*HZ3(0,-1,0)
     $  - HZ3(0,1,0)
      HYZ3(3,2,2) =
     $  + Zeta3
     $  - HXZ1(1) *HZ2(1,1)
     $  - HXZ2(1,0) *HZ1(1)
     $  - HXZ3(1,0,0)
     $  - HZ1(1) *Zeta2
     $  - HZ3(1,1,0)
      HYZ3(3,2,3) =
     $  - 2.000000000000000d+00*Zeta3
     $  + HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(0,1)
     $  - HXZ2(1,0) *HZ1(0)
     $  - HXZ3(1,0,1)
     $  - HZ1(0) *Zeta2
     $  - HZ3(0,1,0)
      HYZ3(3,3,0) =
     $  + Zeta3
     $  + HXZ1(1) *Zeta2
     $  + HXZ1(1) *HZ2(0,0)
     $  - HXZ1(1) *HZ2(0,1)
     $  - HXZ2(1,1) *HZ1(1)
     $  - HXZ3(1,1,2)
     $  + HZ1(0) *Zeta2
     $  + 2.000000000000000d+00*HZ3(0,0,0)
     $  - HZ3(0,0,1)
      HYZ3(3,3,1) =
     $  + Zeta3
     $  - HXZ1(1) *Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ2(-1,0)
     $  - HXZ2(1,1) *HZ1(0)
     $  - HXZ3(1,1,3)
     $  - HZ1( -1)*Zeta2
     $  - 2.000000000000000d+00*HZ3(-1,-1,0)
     $  - HZ3( -1,0,0)
      HYZ3(3,3,2) =
     $  + Zeta3
     $  - HXZ1(1) *Zeta2
     $  - HXZ1(1) *HZ2(1,0)
     $  - HXZ2(1,1) *HZ1(1)
     $  - HXZ3(1,1,0)
     $  - HZ3(1,0,0)
      HYZ3(3,3,3) =
     $  - HXZ1(1) *HZ2(0,0)
     $  - HXZ2(1,1) *HZ1(0)
     $  - HXZ3(1,1,1)
     $  - HZ3(0,0,0)
      if (nmax.eq.3) return
      HYZ4(0,0,0,0) =
     $  + HXZ1(2) *HZ3(1,1,1)
     $  + HXZ2(2,2) *HZ2(1,1)
     $  + HXZ3(2,2,2) *HZ1(1)
     $  + HXZ4(2,2,2,2)
     $  + HZ4(1,1,1,1)
      HYZ4(0,0,0,1) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(1,0)
     $  + HXZ3(2,2,2) *HZ1(0)
     $  + HXZ4(2,2,2,3)
     $  - HZ1(1) *Zeta3
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(1,1,1,0)
      HYZ4(0,0,0,2) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ2(2,2) *Zeta2
     $  + HXZ3(2,2,2) *HZ1(1)
     $  + HXZ4(2,2,2,0)
      HYZ4(0,0,0,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(0,0)
     $  + HXZ2(2,2) *HZ2(1,0)
     $  + HXZ3(2,2,2) *HZ1(0)
     $  + HXZ4(2,2,2,1)
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,1,0)
      HYZ4(0,0,1,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(1,0,1)
     $  - HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(0,1)
     $  + HXZ3(2,2,3) *HZ1(1)
     $  + HXZ4(2,2,3,2)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,1) *Zeta2
     $  + HZ4(1,1,0,1)
      HYZ4(0,0,1,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ2(2,2) *HZ2(0,0)
     $  + HXZ3(2,2,3) *HZ1(0)
     $  + HXZ4(2,2,3,3)
     $  - HZ1(1) *Zeta3
     $  + HZ4(1,1,0,0)
      HYZ4(0,0,1,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,1,0)
     $  + HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(0,0)
     $  + HXZ2(2,2) *HZ2(1,0)
     $  + HXZ3(2,2,3) *HZ1(1)
     $  + HXZ4(2,2,3,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,0,1,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(-1,1,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,0,0)
     $  - HXZ2(2,2) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(2,2)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(2,2)*HZ2(0,0)
     $  + HXZ3(2,2,3) *HZ1(0)
     $  + HXZ4(2,2,3,1)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0,1,1,0)
     $  - HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(0,0,2,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ2(2,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,2)*HZ2(1,1)
     $  + HXZ3(2,2,0) *HZ1(1)
     $  + HXZ4(2,2,0,2)
     $  - HZ1(1) *Zeta3
      HYZ4(0,0,2,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(0,1,0)
     $  - HXZ1(2) *HZ3(1,1,0)
     $  - HXZ2(2,2) *Zeta2
     $  - HXZ2(2,2) *HZ2(0,0)
     $  + HXZ2(2,2) *HZ2(0,1)
     $  + HXZ3(2,2,0) *HZ1(0)
     $  + HXZ4(2,2,0,3)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,1,0)
      HYZ4(0,0,2,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ2(2,2) *HZ2(1,1)
     $  + HXZ3(2,2,0) *HZ1(1)
     $  + HXZ4(2,2,0,0)
      HYZ4(0,0,2,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(1,0,0)
     $  - HXZ1(2) *HZ3(1,1,0)
     $  - HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(0,1)
     $  + HXZ3(2,2,0) *HZ1(0)
     $  + HXZ4(2,2,0,1)
     $  - HZ1(0) *Zeta3
     $  - HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,0)
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  - HZ4(1,1,1,0)
      HYZ4(0,0,3,0) =
     $  - 5.250000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,0,1)
     $  - HXZ1(2) *HZ3(0,1,0)
     $  - HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  - HXZ2(2,2) *Zeta2
     $  - HXZ2(2,2) *HZ2(0,0)
     $  + HXZ2(2,2) *HZ2(0,1)
     $  + HXZ3(2,2,1) *HZ1(1)
     $  + HXZ4(2,2,1,2)
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
      HYZ4(0,0,3,1) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  - HXZ1(2) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,-1,0)
     $  + HXZ2(2,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,2)*HZ2(-1,0)
     $  + HXZ3(2,2,1) *HZ1(0)
     $  + HXZ4(2,2,1,3)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,0,0)
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  - HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - HZ4(0,1,1,0)
     $  + HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(0,0,3,2) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,1,0)
     $  + HXZ2(2,2) *Zeta2
     $  + HXZ2(2,2) *HZ2(1,0)
     $  + HXZ3(2,2,1) *HZ1(1)
     $  + HXZ4(2,2,1,0)
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,0,3,3) =
     $  - 1.750000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ2(2,2) *HZ2(0,0)
     $  + HXZ3(2,2,1) *HZ1(0)
     $  + HXZ4(2,2,1,1)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,1,0,0)
      HYZ4(0,1,0,0) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(0,1,1)
     $  + HXZ2(2,3) *HZ2(1,1)
     $  + HXZ3(2,3,2) *HZ1(1)
     $  + HXZ4(2,3,2,2)
     $  - HZ1(1) *Zeta3
     $  + HZ4(1,0,1,1)
      HYZ4(0,1,0,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(1,0)
     $  + HXZ3(2,3,2) *HZ1(0)
     $  + HXZ4(2,3,2,3)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(1,0,1,0)
      HYZ4(0,1,0,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(1,0,0)
     $  - HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,3) *Zeta2
     $  + HXZ3(2,3,2) *HZ1(1)
     $  + HXZ4(2,3,2,0)
     $  + 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  - 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  - HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
     $  - 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,1,0,3) =
     $  + 5.500000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(0,0)
     $  + HXZ2(2,3) *HZ2(1,0)
     $  + HXZ3(2,3,2) *HZ1(0)
     $  + HXZ4(2,3,2,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,1,1,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,1)
     $  - HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(0,1)
     $  + HXZ3(2,3,3) *HZ1(1)
     $  + HXZ4(2,3,3,2)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(1,0,0,1)
      HYZ4(0,1,1,1) =
     $  + Zeta4
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ2(2,3) *HZ2(0,0)
     $  + HXZ3(2,3,3) *HZ1(0)
     $  + HXZ4(2,3,3,3)
     $  + HZ4(1,0,0,0)
      HYZ4(0,1,1,2) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(0,0)
     $  + HXZ2(2,3) *HZ2(1,0)
     $  + HXZ3(2,3,3) *HZ1(1)
     $  + HXZ4(2,3,3,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(0,1,1,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(2)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  - HXZ2(2,3) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(2,3)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(2,3)*HZ2(0,0)
     $  + HXZ3(2,3,3) *HZ1(0)
     $  + HXZ4(2,3,3,1)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,1,0,0)
     $  + HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
      HYZ4(0,1,2,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,1)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  - HXZ2(2,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,3)*HZ2(1,1)
     $  + HXZ3(2,3,0) *HZ1(1)
     $  + HXZ4(2,3,0,2)
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,1)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,1)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,1,2,1) =
     $  + Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - HXZ2(2,3) *Zeta2
     $  - HXZ2(2,3) *HZ2(0,0)
     $  + HXZ2(2,3) *HZ2(0,1)
     $  + HXZ3(2,3,0) *HZ1(0)
     $  + HXZ4(2,3,0,3)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,1,2,2) =
     $  + Zeta4
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,3) *HZ2(1,1)
     $  + HXZ3(2,3,0) *HZ1(1)
     $  + HXZ4(2,3,0,0)
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,1,2,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(0,1)
     $  + HXZ3(2,3,0) *HZ1(0)
     $  + HXZ4(2,3,0,1)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,1,3,0) =
     $  + 5.000000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,1)
     $  - HXZ2(2,3) *Zeta2
     $  - HXZ2(2,3) *HZ2(0,0)
     $  + HXZ2(2,3) *HZ2(0,1)
     $  + HXZ3(2,3,1) *HZ1(1)
     $  + HXZ4(2,3,1,2)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  - HZ4( -1,0,0,1)
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,1)
      HYZ4(0,1,3,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + HXZ2(2,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,3)*HZ2(-1,0)
     $  + HXZ3(2,3,1) *HZ1(0)
     $  + HXZ4(2,3,1,3)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  - 4.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
      HYZ4(0,1,3,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(-1,1,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,0,0)
     $  + HXZ2(2,3) *Zeta2
     $  + HXZ2(2,3) *HZ2(1,0)
     $  + HXZ3(2,3,1) *HZ1(1)
     $  + HXZ4(2,3,1,0)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 4.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 4.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(1,1,-1,0)
     $  + 4.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(0,1,3,3) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(2)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ2(2,3) *HZ2(0,0)
     $  + HXZ3(2,3,1) *HZ1(0)
     $  + HXZ4(2,3,1,1)
     $  - HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,1,0,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
      HYZ4(0,2,0,0) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(1,1,1)
     $  + HXZ2(2,0) *HZ2(1,1)
     $  + HXZ3(2,0,2) *HZ1(1)
     $  + HXZ4(2,0,2,2)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,1) *Zeta2
      HYZ4(0,2,0,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(1,0)
     $  + HXZ3(2,0,2) *HZ1(0)
     $  + HXZ4(2,0,2,3)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,1,0)
      HYZ4(0,2,0,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ2(2,0) *Zeta2
     $  + HXZ3(2,0,2) *HZ1(1)
     $  + HXZ4(2,0,2,0)
      HYZ4(0,2,0,3) =
     $  + 5.500000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,1)
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(0,0)
     $  + HXZ2(2,0) *HZ2(1,0)
     $  + HXZ3(2,0,2) *HZ1(0)
     $  + HXZ4(2,0,2,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,0,0)
     $  - HZ4(1,1,1,0)
      HYZ4(0,2,1,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,1,1)
     $  - HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(0,1)
     $  + HXZ3(2,0,3) *HZ1(1)
     $  + HXZ4(2,0,3,2)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,1)
     $  - HZ4(1,1,0,1)
      HYZ4(0,2,1,1) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,0,1)
     $  + HXZ2(2,0) *HZ2(0,0)
     $  + HXZ3(2,0,3) *HZ1(0)
     $  + HXZ4(2,0,3,3)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,0,0)
      HYZ4(0,2,1,2) =
     $  + Zeta4
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(0)*Zeta2
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,0,1)
     $  - HXZ1(2) *HZ3(0,1,0)
     $  - HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(0,0)
     $  + HXZ2(2,0) *HZ2(1,0)
     $  + HXZ3(2,0,3) *HZ1(1)
     $  + HXZ4(2,0,3,0)
     $  - 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  - 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,2,1,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,1)
     $  - HXZ2(2,0) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(2,0)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(2,0)*HZ2(0,0)
     $  + HXZ3(2,0,3) *HZ1(0)
     $  + HXZ4(2,0,3,1)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(0,2,2,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(2)*HZ3(1,1,1)
     $  - HXZ2(2,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,0)*HZ2(1,1)
     $  + HXZ3(2,0,0) *HZ1(1)
     $  + HXZ4(2,0,0,2)
     $  - HZ1(1) *Zeta3
      HYZ4(0,2,2,1) =
     $  + Zeta4
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  - HXZ1(2) *HZ3(0,0,1)
     $  + HXZ1(2) *HZ3(0,1,1)
     $  - HXZ2(2,0) *Zeta2
     $  - HXZ2(2,0) *HZ2(0,0)
     $  + HXZ2(2,0) *HZ2(0,1)
     $  + HXZ3(2,0,0) *HZ1(0)
     $  + HXZ4(2,0,0,3)
     $  + HZ1(1) *Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,1,0)
      HYZ4(0,2,2,2) =
     $  + Zeta4
     $  + HXZ1(2) *HZ3(1,1,1)
     $  + HXZ2(2,0) *HZ2(1,1)
     $  + HXZ3(2,0,0) *HZ1(1)
     $  + HXZ4(2,0,0,0)
      HYZ4(0,2,2,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(0,1,1)
     $  - HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(0,1)
     $  + HXZ3(2,0,0) *HZ1(0)
     $  + HXZ4(2,0,0,1)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,1,0)
      HYZ4(0,2,3,0) =
     $  + 5.000000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,1,1)
     $  - HXZ2(2,0) *Zeta2
     $  - HXZ2(2,0) *HZ2(0,0)
     $  + HXZ2(2,0) *HZ2(0,1)
     $  + HXZ3(2,0,1) *HZ1(1)
     $  + HXZ4(2,0,1,2)
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
      HYZ4(0,2,3,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,1)
     $  + HXZ2(2,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,0)*HZ2(-1,0)
     $  + HXZ3(2,0,1) *HZ1(0)
     $  + HXZ4(2,0,1,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(0,2,3,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ2(2,0) *Zeta2
     $  + HXZ2(2,0) *HZ2(1,0)
     $  + HXZ3(2,0,1) *HZ1(1)
     $  + HXZ4(2,0,1,0)
     $  + 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  - 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  - HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
     $  - 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,2,3,3) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,1)
     $  + HXZ2(2,0) *HZ2(0,0)
     $  + HXZ3(2,0,1) *HZ1(0)
     $  + HXZ4(2,0,1,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,0,0,0)
     $  - HZ4(1,1,0,0)
      HYZ4(0,3,0,0) =
     $  + 5.250000000000000d+00*Zeta4
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  - HXZ1(2) *HZ3(0,0,1)
     $  + HXZ1(2) *HZ3(0,1,1)
     $  + HXZ2(2,1) *HZ2(1,1)
     $  + HXZ3(2,1,2) *HZ1(1)
     $  + HXZ4(2,1,2,2)
     $  + 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 3.000000000000000d+00*HZ4(0,0,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,0,1,1)
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  + HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
      HYZ4(0,3,0,1) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(1,0)
     $  + HXZ3(2,1,2) *HZ1(0)
     $  + HXZ4(2,1,2,3)
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,3,0,2) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(1)*Zeta2
     $  - HXZ1(2) *HZ3(0,1,0)
     $  - HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,1) *Zeta2
     $  + HXZ3(2,1,2) *HZ1(1)
     $  + HXZ4(2,1,2,0)
     $  - 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  - 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,3,0,3) =
     $  + 4.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  + HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(0,0)
     $  + HXZ2(2,1) *HZ2(1,0)
     $  + HXZ3(2,1,2) *HZ1(0)
     $  + HXZ4(2,1,2,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,3,1,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  - HXZ1(2) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,1)
     $  - HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(0,1)
     $  + HXZ3(2,1,3) *HZ1(1)
     $  + HXZ4(2,1,3,2)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  + HZ4( -1,0,0,1)
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - HZ4(0,1,0,1)
     $  - HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,1)
      HYZ4(0,3,1,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  + HXZ2(2,1) *HZ2(0,0)
     $  + HXZ3(2,1,3) *HZ1(0)
     $  + HXZ4(2,1,3,3)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
     $  - HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + HZ4(1, -1,0,0)
      HYZ4(0,3,1,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  - HXZ1(2) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(1,-1,0)
     $  + HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(0,0)
     $  + HXZ2(2,1) *HZ2(1,0)
     $  + HXZ3(2,1,3) *HZ1(1)
     $  + HXZ4(2,1,3,0)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 4.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + 4.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  + 4.000000000000000d+00*HZ4(-1,1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
     $  + 4.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(0,3,1,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ1(-1)*Zeta2
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  - HXZ2(2,1) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(2,1)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(2,1)*HZ2(0,0)
     $  + HXZ3(2,1,3) *HZ1(0)
     $  + HXZ4(2,1,3,1)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 4.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
      HYZ4(0,3,2,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  - HXZ1(2) *HZ3(1,0,0)
     $  + HXZ1(2) *HZ3(1,0,1)
     $  + HXZ1(2) *HZ3(1,1,0)
     $  - HXZ2(2,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,1)*HZ2(1,1)
     $  + HXZ3(2,1,0) *HZ1(1)
     $  + HXZ4(2,1,0,2)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  - HZ4(1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,1)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,3,2,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,-1,0)
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - HXZ2(2,1) *Zeta2
     $  - HXZ2(2,1) *HZ2(0,0)
     $  + HXZ2(2,1) *HZ2(0,1)
     $  + HXZ3(2,1,0) *HZ1(0)
     $  + HXZ4(2,1,0,3)
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,3,2,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(1)*Zeta2
     $  + HXZ1(2) *HZ3(1,1,0)
     $  + HXZ2(2,1) *HZ2(1,1)
     $  + HXZ3(2,1,0) *HZ1(1)
     $  + HXZ4(2,1,0,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + 3.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(0,3,2,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(2)*Zeta3
     $  + HXZ1(2) *HZ1(0)*Zeta2
     $  + HXZ1(2) *HZ3(0,1,0)
     $  - HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(0,1)
     $  + HXZ3(2,1,0) *HZ1(0)
     $  + HXZ4(2,1,0,1)
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(0,3,3,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(2) *Zeta3
     $  - HXZ1(2) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(2)*HZ3(0,0,0)
     $  + HXZ1(2) *HZ3(0,0,1)
     $  - HXZ2(2,1) *Zeta2
     $  - HXZ2(2,1) *HZ2(0,0)
     $  + HXZ2(2,1) *HZ2(0,1)
     $  + HXZ3(2,1,1) *HZ1(1)
     $  + HXZ4(2,1,1,2)
     $  - HZ1(1) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - HZ2(1,0) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
      HYZ4(0,3,3,1) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(2)*HZ3(-1,-1,0)
     $  + HXZ1(2) *HZ3(-1,0,0)
     $  + HXZ2(2,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(2,1)*HZ2(-1,0)
     $  + HXZ3(2,1,1) *HZ1(0)
     $  + HXZ4(2,1,1,3)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  - HZ2(0, -1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + HZ4(1, -1,0,0)
      HYZ4(0,3,3,2) =
     $  + Zeta4
     $  - HXZ1(2) *Zeta3
     $  + HXZ1(2) *HZ3(1,0,0)
     $  + HXZ2(2,1) *Zeta2
     $  + HXZ2(2,1) *HZ2(1,0)
     $  + HXZ3(2,1,1) *HZ1(1)
     $  + HXZ4(2,1,1,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(0,3,3,3) =
     $  + Zeta4
     $  + HXZ1(2) *HZ3(0,0,0)
     $  + HXZ2(2,1) *HZ2(0,0)
     $  + HXZ3(2,1,1) *HZ1(0)
     $  + HXZ4(2,1,1,1)
     $  + HZ4(0,0,0,0)
     $  + HZ4(1,0,0,0)
      HYZ4(1,0,0,0) =
     $  - Zeta4
     $  + HXZ1(3) *HZ3(1,1,1)
     $  + HXZ2(3,2) *HZ2(1,1)
     $  + HXZ3(3,2,2) *HZ1(1)
     $  + HXZ4(3,2,2,2)
     $  + HZ4(0,1,1,1)
      HYZ4(1,0,0,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(1,0)
     $  + HXZ3(3,2,2) *HZ1(0)
     $  + HXZ4(3,2,2,3)
     $  - HZ1(0) *Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + HZ4(0,1,1,0)
      HYZ4(1,0,0,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ2(3,2) *Zeta2
     $  + HXZ3(3,2,2) *HZ1(1)
     $  + HXZ4(3,2,2,0)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,1,0)
      HYZ4(1,0,0,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(0,0)
     $  + HXZ2(3,2) *HZ2(1,0)
     $  + HXZ3(3,2,2) *HZ1(0)
     $  + HXZ4(3,2,2,1)
     $  - HZ1(0) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
      HYZ4(1,0,1,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(1,0,1)
     $  - HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(0,1)
     $  + HXZ3(3,2,3) *HZ1(1)
     $  + HXZ4(3,2,3,2)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0,1,0,1)
      HYZ4(1,0,1,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ2(3,2) *HZ2(0,0)
     $  + HXZ3(3,2,3) *HZ1(0)
     $  + HXZ4(3,2,3,3)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,1,0,0)
      HYZ4(1,0,1,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,1,0)
     $  + HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(0,0)
     $  + HXZ2(3,2) *HZ2(1,0)
     $  + HXZ3(3,2,3) *HZ1(1)
     $  + HXZ4(3,2,3,0)
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(1,0,1,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(-1,1,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,0,0)
     $  - HXZ2(3,2) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(3,2)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(3,2)*HZ2(0,0)
     $  + HXZ3(3,2,3) *HZ1(0)
     $  + HXZ4(3,2,3,1)
     $  - HZ1( -1)*Zeta3
     $  + 4.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
      HYZ4(1,0,2,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ2(3,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,2)*HZ2(1,1)
     $  + HXZ3(3,2,0) *HZ1(1)
     $  + HXZ4(3,2,0,2)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(1,0,2,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(0,1,0)
     $  - HXZ1(3) *HZ3(1,1,0)
     $  - HXZ2(3,2) *Zeta2
     $  - HXZ2(3,2) *HZ2(0,0)
     $  + HXZ2(3,2) *HZ2(0,1)
     $  + HXZ3(3,2,0) *HZ1(0)
     $  + HXZ4(3,2,0,3)
     $  - HZ1(0) *Zeta3
     $  - 4.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  - HZ4(1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(1,0,2,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ2(3,2) *HZ2(1,1)
     $  + HXZ3(3,2,0) *HZ1(1)
     $  + HXZ4(3,2,0,0)
     $  - HZ1(0) *Zeta3
     $  - HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,0)
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  - HZ4(1,1,1,0)
      HYZ4(1,0,2,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(1,0,0)
     $  - HXZ1(3) *HZ3(1,1,0)
     $  - HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(0,1)
     $  + HXZ3(3,2,0) *HZ1(0)
     $  + HXZ4(3,2,0,1)
     $  + 4.000000000000000d+00*HZ1(0)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + 4.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(1,0,3,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,0,1)
     $  - HXZ1(3) *HZ3(0,1,0)
     $  - HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  - HXZ2(3,2) *Zeta2
     $  - HXZ2(3,2) *HZ2(0,0)
     $  + HXZ2(3,2) *HZ2(0,1)
     $  + HXZ3(3,2,1) *HZ1(1)
     $  + HXZ4(3,2,1,2)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 8.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
      HYZ4(1,0,3,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  - HXZ1(3) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,-1,0)
     $  + HXZ2(3,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,2)*HZ2(-1,0)
     $  + HXZ3(3,2,1) *HZ1(0)
     $  + HXZ4(3,2,1,3)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
      HYZ4(1,0,3,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,1,0)
     $  + HXZ2(3,2) *Zeta2
     $  + HXZ2(3,2) *HZ2(1,0)
     $  + HXZ3(3,2,1) *HZ1(1)
     $  + HXZ4(3,2,1,0)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(1,0,3,3) =
     $  - Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ2(3,2) *HZ2(0,0)
     $  + HXZ3(3,2,1) *HZ1(0)
     $  + HXZ4(3,2,1,1)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(1,1,0,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(0,1,1)
     $  + HXZ2(3,3) *HZ2(1,1)
     $  + HXZ3(3,3,2) *HZ1(1)
     $  + HXZ4(3,3,2,2)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,0,1,1)
      HYZ4(1,1,0,1) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(1,0)
     $  + HXZ3(3,3,2) *HZ1(0)
     $  + HXZ4(3,3,2,3)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ4(0,0,1,0)
      HYZ4(1,1,0,2) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(1,0,0)
     $  - HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,3) *Zeta2
     $  + HXZ3(3,3,2) *HZ1(1)
     $  + HXZ4(3,3,2,0)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,0,0,0)
     $  - HZ4(1,1,0,0)
      HYZ4(1,1,0,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(0,0)
     $  + HXZ2(3,3) *HZ2(1,0)
     $  + HXZ3(3,3,2) *HZ1(0)
     $  + HXZ4(3,3,2,1)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
      HYZ4(1,1,1,0) =
     $  - Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,1)
     $  - HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(0,1)
     $  + HXZ3(3,3,3) *HZ1(1)
     $  + HXZ4(3,3,3,2)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,1)
      HYZ4(1,1,1,1) =
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ2(3,3) *HZ2(0,0)
     $  + HXZ3(3,3,3) *HZ1(0)
     $  + HXZ4(3,3,3,3)
     $  + HZ4(0,0,0,0)
      HYZ4(1,1,1,2) =
     $  + Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(0,0)
     $  + HXZ2(3,3) *HZ2(1,0)
     $  + HXZ3(3,3,3) *HZ1(1)
     $  + HXZ4(3,3,3,0)
     $  + HZ4(0,0,0,0)
     $  + HZ4(1,0,0,0)
      HYZ4(1,1,1,3) =
     $  - Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(3)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  - HXZ2(3,3) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(3,3)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(3,3)*HZ2(0,0)
     $  + HXZ3(3,3,3) *HZ1(0)
     $  + HXZ4(3,3,3,1)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
      HYZ4(1,1,2,0) =
     $  - Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,1)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  - HXZ2(3,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,3)*HZ2(1,1)
     $  + HXZ3(3,3,0) *HZ1(1)
     $  + HXZ4(3,3,0,2)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,1,0,0)
      HYZ4(1,1,2,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - HXZ2(3,3) *Zeta2
     $  - HXZ2(3,3) *HZ2(0,0)
     $  + HXZ2(3,3) *HZ2(0,1)
     $  + HXZ3(3,3,0) *HZ1(0)
     $  + HXZ4(3,3,0,3)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(1,1,2,2) =
     $  - 1.750000000000000d+00*Zeta4
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,3) *HZ2(1,1)
     $  + HXZ3(3,3,0) *HZ1(1)
     $  + HXZ4(3,3,0,0)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,1,0,0)
      HYZ4(1,1,2,3) =
     $  - Zeta4
     $  - HXZ1(3) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(0,1)
     $  + HXZ3(3,3,0) *HZ1(0)
     $  + HXZ4(3,3,0,1)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(1,1,3,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,1)
     $  - HXZ2(3,3) *Zeta2
     $  - HXZ2(3,3) *HZ2(0,0)
     $  + HXZ2(3,3) *HZ2(0,1)
     $  + HXZ3(3,3,1) *HZ1(1)
     $  + HXZ4(3,3,1,2)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
      HYZ4(1,1,3,1) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + HXZ2(3,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,3)*HZ2(-1,0)
     $  + HXZ3(3,3,1) *HZ1(0)
     $  + HXZ4(3,3,1,3)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 3.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  + 6.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  - 5.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(1,1,3,2) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(-1,1,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,0,0)
     $  + HXZ2(3,3) *Zeta2
     $  + HXZ2(3,3) *HZ2(1,0)
     $  + HXZ3(3,3,1) *HZ1(1)
     $  + HXZ4(3,3,1,0)
     $  - HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,1,0,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
      HYZ4(1,1,3,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(3)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ2(3,3) *HZ2(0,0)
     $  + HXZ3(3,3,1) *HZ1(0)
     $  + HXZ4(3,3,1,1)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - 4.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + 6.000000000000000d+00*HZ4(-1,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 6.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(0,0,0,0)
      HYZ4(1,2,0,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(1,1,1)
     $  + HXZ2(3,0) *HZ2(1,1)
     $  + HXZ3(3,0,2) *HZ1(1)
     $  + HXZ4(3,0,2,2)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(0,0,1,1)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + HZ4(1,1,1,0)
      HYZ4(1,2,0,1) =
     $  + 5.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(1,0)
     $  + HXZ3(3,0,2) *HZ1(0)
     $  + HXZ4(3,0,2,3)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 4.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(1,2,0,2) =
     $  + 5.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ2(3,0) *Zeta2
     $  + HXZ3(3,0,2) *HZ1(1)
     $  + HXZ4(3,0,2,0)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,0,0)
     $  - HZ4(1,1,1,0)
      HYZ4(1,2,0,3) =
     $  + 2.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,1)
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(0,0)
     $  + HXZ2(3,0) *HZ2(1,0)
     $  + HXZ3(3,0,2) *HZ1(0)
     $  + HXZ4(3,0,2,1)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 4.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(1,2,1,0) =
     $  - 3.500000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,1,1)
     $  - HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(0,1)
     $  + HXZ3(3,0,3) *HZ1(1)
     $  + HXZ4(3,0,3,2)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(1,2,1,1) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,0,1)
     $  + HXZ2(3,0) *HZ2(0,0)
     $  + HXZ3(3,0,3) *HZ1(0)
     $  + HXZ4(3,0,3,3)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
      HYZ4(1,2,1,2) =
     $  + 4.750000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(0)*Zeta2
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,0,1)
     $  - HXZ1(3) *HZ3(0,1,0)
     $  - HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(0,0)
     $  + HXZ2(3,0) *HZ2(1,0)
     $  + HXZ3(3,0,3) *HZ1(1)
     $  + HXZ4(3,0,3,0)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(1,2,1,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,1)
     $  - HXZ2(3,0) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(3,0)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(3,0)*HZ2(0,0)
     $  + HXZ3(3,0,3) *HZ1(0)
     $  + HXZ4(3,0,3,1)
     $  - HZ1( -1)*Zeta3
     $  + HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
      HYZ4(1,2,2,0) =
     $  - 3.500000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(3)*HZ3(1,1,1)
     $  - HXZ2(3,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,0)*HZ2(1,1)
     $  + HXZ3(3,0,0) *HZ1(1)
     $  + HXZ4(3,0,0,2)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(1,2,2,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  - HXZ1(3) *HZ3(0,0,1)
     $  + HXZ1(3) *HZ3(0,1,1)
     $  - HXZ2(3,0) *Zeta2
     $  - HXZ2(3,0) *HZ2(0,0)
     $  + HXZ2(3,0) *HZ2(0,1)
     $  + HXZ3(3,0,0) *HZ1(0)
     $  + HXZ4(3,0,0,3)
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
      HYZ4(1,2,2,2) =
     $  + 1.750000000000000d+00*Zeta4
     $  + HXZ1(3) *HZ3(1,1,1)
     $  + HXZ2(3,0) *HZ2(1,1)
     $  + HXZ3(3,0,0) *HZ1(1)
     $  + HXZ4(3,0,0,0)
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,1,0)
      HYZ4(1,2,2,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(0,1,1)
     $  - HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(0,1)
     $  + HXZ3(3,0,0) *HZ1(0)
     $  + HXZ4(3,0,0,1)
     $  - HZ1(0) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
      HYZ4(1,2,3,0) =
     $  - 2.750000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,1,1)
     $  - HXZ2(3,0) *Zeta2
     $  - HXZ2(3,0) *HZ2(0,0)
     $  + HXZ2(3,0) *HZ2(0,1)
     $  + HXZ3(3,0,1) *HZ1(1)
     $  + HXZ4(3,0,1,2)
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
     $  - HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(1,2,3,1) =
     $  + 5.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,1)
     $  + HXZ2(3,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,0)*HZ2(-1,0)
     $  + HXZ3(3,0,1) *HZ1(0)
     $  + HXZ4(3,0,1,3)
     $  + HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(1,2,3,2) =
     $  + 5.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ2(3,0) *Zeta2
     $  + HXZ2(3,0) *HZ2(1,0)
     $  + HXZ3(3,0,1) *HZ1(1)
     $  + HXZ4(3,0,1,0)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(1,2,3,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,1)
     $  + HXZ2(3,0) *HZ2(0,0)
     $  + HXZ3(3,0,1) *HZ1(0)
     $  + HXZ4(3,0,1,1)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
      HYZ4(1,3,0,0) =
     $  - Zeta4
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  - HXZ1(3) *HZ3(0,0,1)
     $  + HXZ1(3) *HZ3(0,1,1)
     $  + HXZ2(3,1) *HZ2(1,1)
     $  + HXZ3(3,1,2) *HZ1(1)
     $  + HXZ4(3,1,2,2)
     $  + HZ1( -1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ4( -1,0,0,0)
     $  + HZ4( -1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  - 3.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
      HYZ4(1,3,0,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(1,0)
     $  + HXZ3(3,1,2) *HZ1(0)
     $  + HXZ4(3,1,2,3)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
      HYZ4(1,3,0,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(1)*Zeta2
     $  - HXZ1(3) *HZ3(0,1,0)
     $  - HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,1) *Zeta2
     $  + HXZ3(3,1,2) *HZ1(1)
     $  + HXZ4(3,1,2,0)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(1,3,0,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  + HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(0,0)
     $  + HXZ2(3,1) *HZ2(1,0)
     $  + HXZ3(3,1,2) *HZ1(0)
     $  + HXZ4(3,1,2,1)
     $  - HZ1( -1)*Zeta3
     $  + HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
      HYZ4(1,3,1,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  - HXZ1(3) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,1)
     $  - HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(0,1)
     $  + HXZ3(3,1,3) *HZ1(1)
     $  + HXZ4(3,1,3,2)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
      HYZ4(1,3,1,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  + HXZ2(3,1) *HZ2(0,0)
     $  + HXZ3(3,1,3) *HZ1(0)
     $  + HXZ4(3,1,3,3)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 3.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - 6.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
      HYZ4(1,3,1,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  - HXZ1(3) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(1,-1,0)
     $  + HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(0,0)
     $  + HXZ2(3,1) *HZ2(1,0)
     $  + HXZ3(3,1,3) *HZ1(1)
     $  + HXZ4(3,1,3,0)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  - 4.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
      HYZ4(1,3,1,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ1(-1)*Zeta2
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  - HXZ2(3,1) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(3,1)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(3,1)*HZ2(0,0)
     $  + HXZ3(3,1,3) *HZ1(0)
     $  + HXZ4(3,1,3,1)
     $  - 4.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 4.000000000000000d+00*HZ1(0)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 4.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + 8.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 8.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(1,3,2,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  - HXZ1(3) *HZ3(1,0,0)
     $  + HXZ1(3) *HZ3(1,0,1)
     $  + HXZ1(3) *HZ3(1,1,0)
     $  - HXZ2(3,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,1)*HZ2(1,1)
     $  + HXZ3(3,1,0) *HZ1(1)
     $  + HXZ4(3,1,0,2)
     $  - HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(1,3,2,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,-1,0)
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - HXZ2(3,1) *Zeta2
     $  - HXZ2(3,1) *HZ2(0,0)
     $  + HXZ2(3,1) *HZ2(0,1)
     $  + HXZ3(3,1,0) *HZ1(0)
     $  + HXZ4(3,1,0,3)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - HZ4( -1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
      HYZ4(1,3,2,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(1)*Zeta2
     $  + HXZ1(3) *HZ3(1,1,0)
     $  + HXZ2(3,1) *HZ2(1,1)
     $  + HXZ3(3,1,0) *HZ1(1)
     $  + HXZ4(3,1,0,0)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0,1,1,0)
     $  - HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(1,3,2,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(3)*Zeta3
     $  + HXZ1(3) *HZ1(0)*Zeta2
     $  + HXZ1(3) *HZ3(0,1,0)
     $  - HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(0,1)
     $  + HXZ3(3,1,0) *HZ1(0)
     $  + HXZ4(3,1,0,1)
     $  - HZ1( -1)*Zeta3
     $  + 4.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
      HYZ4(1,3,3,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  - HXZ1(3) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(3)*HZ3(0,0,0)
     $  + HXZ1(3) *HZ3(0,0,1)
     $  - HXZ2(3,1) *Zeta2
     $  - HXZ2(3,1) *HZ2(0,0)
     $  + HXZ2(3,1) *HZ2(0,1)
     $  + HXZ3(3,1,1) *HZ1(1)
     $  + HXZ4(3,1,1,2)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 5.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 8.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
      HYZ4(1,3,3,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(3)*HZ3(-1,-1,0)
     $  + HXZ1(3) *HZ3(-1,0,0)
     $  + HXZ2(3,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(3,1)*HZ2(-1,0)
     $  + HXZ3(3,1,1) *HZ1(0)
     $  + HXZ4(3,1,1,3)
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
      HYZ4(1,3,3,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(3) *Zeta3
     $  + HXZ1(3) *HZ3(1,0,0)
     $  + HXZ2(3,1) *Zeta2
     $  + HXZ2(3,1) *HZ2(1,0)
     $  + HXZ3(3,1,1) *HZ1(1)
     $  + HXZ4(3,1,1,0)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,1,0,0)
     $  + HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(1,0,0,0)
      HYZ4(1,3,3,3) =
     $  - Zeta4
     $  + HXZ1(3) *HZ3(0,0,0)
     $  + HXZ2(3,1) *HZ2(0,0)
     $  + HXZ3(3,1,1) *HZ1(0)
     $  + HXZ4(3,1,1,1)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
      HYZ4(2,0,0,0) =
     $  - Zeta4
     $  + HXZ1(0) *HZ3(1,1,1)
     $  + HXZ2(0,2) *HZ2(1,1)
     $  + HXZ3(0,2,2) *HZ1(1)
     $  + HXZ4(0,2,2,2)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  + 4.000000000000000d+00*HZ4(1,1,1,1)
      HYZ4(2,0,0,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(1,0)
     $  + HXZ3(0,2,2) *HZ1(0)
     $  + HXZ4(0,2,2,3)
     $  + HZ1(1) *Zeta3
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,0,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ2(0,2) *Zeta2
     $  + HXZ3(0,2,2) *HZ1(1)
     $  + HXZ4(0,2,2,0)
     $  - HZ1(1) *Zeta3
      HYZ4(2,0,0,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(0,0)
     $  + HXZ2(0,2) *HZ2(1,0)
     $  + HXZ3(0,2,2) *HZ1(0)
     $  + HXZ4(0,2,2,1)
     $  - HZ1(0) *Zeta3
     $  - HZ1(1) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,1,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(1,0,1)
     $  - HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(0,1)
     $  + HXZ3(0,2,3) *HZ1(1)
     $  + HXZ4(0,2,3,2)
     $  + HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(1,0,1,1)
     $  + HZ4(1,1,0,1)
      HYZ4(2,0,1,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ2(0,2) *HZ2(0,0)
     $  + HXZ3(0,2,3) *HZ1(0)
     $  + HXZ4(0,2,3,3)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,1,0,0)
      HYZ4(2,0,1,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,1,0)
     $  + HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(0,0)
     $  + HXZ2(0,2) *HZ2(1,0)
     $  + HXZ3(0,2,3) *HZ1(1)
     $  + HXZ4(0,2,3,0)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  - HZ4(1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,1)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,1,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(-1,1,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,0,0)
     $  - HXZ2(0,2) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(0,2)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(0,2)*HZ2(0,0)
     $  + HXZ3(0,2,3) *HZ1(0)
     $  + HXZ4(0,2,3,1)
     $  - HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
      HYZ4(2,0,2,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ2(0,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,2)*HZ2(1,1)
     $  + HXZ3(0,2,0) *HZ1(1)
     $  + HXZ4(0,2,0,2)
     $  + 4.000000000000000d+00*HZ1(1)*Zeta3
     $  + 2.000000000000000d+00*HZ2(1,1)*Zeta2
      HYZ4(2,0,2,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(0,1,0)
     $  - HXZ1(0) *HZ3(1,1,0)
     $  - HXZ2(0,2) *Zeta2
     $  - HXZ2(0,2) *HZ2(0,0)
     $  + HXZ2(0,2) *HZ2(0,1)
     $  + HXZ3(0,2,0) *HZ1(0)
     $  + HXZ4(0,2,0,3)
     $  - 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,2,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ2(0,2) *HZ2(1,1)
     $  + HXZ3(0,2,0) *HZ1(1)
     $  + HXZ4(0,2,0,0)
     $  - HZ1(1) *Zeta3
      HYZ4(2,0,2,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(1,0,0)
     $  - HXZ1(0) *HZ3(1,1,0)
     $  - HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(0,1)
     $  + HXZ3(0,2,0) *HZ1(0)
     $  + HXZ4(0,2,0,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,3,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,0,1)
     $  - HXZ1(0) *HZ3(0,1,0)
     $  - HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  - HXZ2(0,2) *Zeta2
     $  - HXZ2(0,2) *HZ2(0,0)
     $  + HXZ2(0,2) *HZ2(0,1)
     $  + HXZ3(0,2,1) *HZ1(1)
     $  + HXZ4(0,2,1,2)
     $  - HZ1(1) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,1)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
      HYZ4(2,0,3,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  - HXZ1(0) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,-1,0)
     $  + HXZ2(0,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,2)*HZ2(-1,0)
     $  + HXZ3(0,2,1) *HZ1(0)
     $  + HXZ4(0,2,1,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ1(1) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(2,0,3,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,1,0)
     $  + HXZ2(0,2) *Zeta2
     $  + HXZ2(0,2) *HZ2(1,0)
     $  + HXZ3(0,2,1) *HZ1(1)
     $  + HXZ4(0,2,1,0)
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,1)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,1,0,1)
     $  + 3.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(2,0,3,3) =
     $  - Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ2(0,2) *HZ2(0,0)
     $  + HXZ3(0,2,1) *HZ1(0)
     $  + HXZ4(0,2,1,1)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,1,0,0)
      HYZ4(2,1,0,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(0,1,1)
     $  + HXZ2(0,3) *HZ2(1,1)
     $  + HXZ3(0,3,2) *HZ1(1)
     $  + HXZ4(0,3,2,2)
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,0,1,1)
     $  + 3.000000000000000d+00*HZ4(0,1,1,1)
      HYZ4(2,1,0,1) =
     $  + 5.000000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(1,0)
     $  + HXZ3(0,3,2) *HZ1(0)
     $  + HXZ4(0,3,2,3)
     $  - HZ2(0,0) *Zeta2
     $  - HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(2,1,0,2) =
     $  + 5.000000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(1,0,0)
     $  - HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,3) *Zeta2
     $  + HXZ3(0,3,2) *HZ1(1)
     $  + HXZ4(0,3,2,0)
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
      HYZ4(2,1,0,3) =
     $  - 2.750000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(0,0)
     $  + HXZ2(0,3) *HZ2(1,0)
     $  + HXZ3(0,3,2) *HZ1(0)
     $  + HXZ4(0,3,2,1)
     $  - HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
     $  - HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(2,1,1,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,1)
     $  - HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(0,1)
     $  + HXZ3(0,3,3) *HZ1(1)
     $  + HXZ4(0,3,3,2)
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
      HYZ4(2,1,1,1) =
     $  - Zeta4
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ2(0,3) *HZ2(0,0)
     $  + HXZ3(0,3,3) *HZ1(0)
     $  + HXZ4(0,3,3,3)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
      HYZ4(2,1,1,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(0,0)
     $  + HXZ2(0,3) *HZ2(1,0)
     $  + HXZ3(0,3,3) *HZ1(1)
     $  + HXZ4(0,3,3,0)
     $  - HZ1(1) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - HZ2(1,0) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
      HYZ4(2,1,1,3) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(0)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  - HXZ2(0,3) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(0,3)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(0,3)*HZ2(0,0)
     $  + HXZ3(0,3,3) *HZ1(0)
     $  + HXZ4(0,3,3,1)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 5.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 8.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
      HYZ4(2,1,2,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,1)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  - HXZ2(0,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,3)*HZ2(1,1)
     $  + HXZ3(0,3,0) *HZ1(1)
     $  + HXZ4(0,3,0,2)
     $  - HZ1(1) *Zeta3
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,1)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
      HYZ4(2,1,2,1) =
     $  - 2.250000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - HXZ2(0,3) *Zeta2
     $  - HXZ2(0,3) *HZ2(0,0)
     $  + HXZ2(0,3) *HZ2(0,1)
     $  + HXZ3(0,3,0) *HZ1(0)
     $  + HXZ4(0,3,0,3)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
      HYZ4(2,1,2,2) =
     $  - 5.250000000000000d+00*Zeta4
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,3) *HZ2(1,1)
     $  + HXZ3(0,3,0) *HZ1(1)
     $  + HXZ4(0,3,0,0)
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
      HYZ4(2,1,2,3) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(0,1)
     $  + HXZ3(0,3,0) *HZ1(0)
     $  + HXZ4(0,3,0,1)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 8.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
      HYZ4(2,1,3,0) =
     $  + 2.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,1)
     $  - HXZ2(0,3) *Zeta2
     $  - HXZ2(0,3) *HZ2(0,0)
     $  + HXZ2(0,3) *HZ2(0,1)
     $  + HXZ3(0,3,1) *HZ1(1)
     $  + HXZ4(0,3,1,2)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  - 4.000000000000000d+00*HZ4(-1,0,1,1)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,1)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 6.000000000000000d+00*HZ4(0,0,0,0)
     $  - 6.000000000000000d+00*HZ4(0,0,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,1,1)
      HYZ4(2,1,3,1) =
     $  + 5.000000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + HXZ2(0,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,3)*HZ2(-1,0)
     $  + HXZ3(0,3,1) *HZ1(0)
     $  + HXZ4(0,3,1,3)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(2,1,3,2) =
     $  + 5.000000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(-1,1,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,0,0)
     $  + HXZ2(0,3) *Zeta2
     $  + HXZ2(0,3) *HZ2(1,0)
     $  + HXZ3(0,3,1) *HZ1(1)
     $  + HXZ4(0,3,1,0)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  - HZ4( -1,0,0,1)
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,0,1)
      HYZ4(2,1,3,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(0)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ2(0,3) *HZ2(0,0)
     $  + HXZ3(0,3,1) *HZ1(0)
     $  + HXZ4(0,3,1,1)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 6.000000000000000d+00*HZ4(0,0,0,0)
     $  + 3.000000000000000d+00*HZ4(0,0,0,1)
      HYZ4(2,2,0,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(1,1,1)
     $  + HXZ2(0,0) *HZ2(1,1)
     $  + HXZ3(0,0,2) *HZ1(1)
     $  + HXZ4(0,0,2,2)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(1,1)*Zeta2
     $  + 6.000000000000000d+00*HZ4(1,1,1,1)
      HYZ4(2,2,0,1) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(1,0)
     $  + HXZ3(0,0,2) *HZ1(0)
     $  + HXZ4(0,0,2,3)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + HZ4(1,1,1,0)
      HYZ4(2,2,0,2) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ2(0,0) *Zeta2
     $  + HXZ3(0,0,2) *HZ1(1)
     $  + HXZ4(0,0,2,0)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,1) *Zeta2
      HYZ4(2,2,0,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,1)
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(0,0)
     $  + HXZ2(0,0) *HZ2(1,0)
     $  + HXZ3(0,0,2) *HZ1(0)
     $  + HXZ4(0,0,2,1)
     $  - HZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(0,0,1,1)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
     $  + HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + HZ4(1,1,1,0)
      HYZ4(2,2,1,0) =
     $  - Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,1,1)
     $  - HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(0,1)
     $  + HXZ3(0,0,3) *HZ1(1)
     $  + HXZ4(0,0,3,2)
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,1)
     $  + 3.000000000000000d+00*HZ4(0,1,1,1)
      HYZ4(2,2,1,1) =
     $  + 1.750000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,0,1)
     $  + HXZ2(0,0) *HZ2(0,0)
     $  + HXZ3(0,0,3) *HZ1(0)
     $  + HXZ4(0,0,3,3)
     $  + HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + 3.000000000000000d+00*HZ4(0,0,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + HZ4(0,0,1,1)
      HYZ4(2,2,1,2) =
     $  + 5.250000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(0)*Zeta2
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,0,1)
     $  - HXZ1(0) *HZ3(0,1,0)
     $  - HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(0,0)
     $  + HXZ2(0,0) *HZ2(1,0)
     $  + HXZ3(0,0,3) *HZ1(1)
     $  + HXZ4(0,0,3,0)
     $  + 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 3.000000000000000d+00*HZ4(0,0,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,0,1,1)
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  + HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
      HYZ4(2,2,1,3) =
     $  - Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,1)
     $  - HXZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(0,0)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(0,0)*HZ2(0,0)
     $  + HXZ3(0,0,3) *HZ1(0)
     $  + HXZ4(0,0,3,1)
     $  + HZ1( -1)*Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ4( -1,0,0,0)
     $  + HZ4( -1,0,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,0,0)
     $  - 3.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
      HYZ4(2,2,2,0) =
     $  - Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(0)*HZ3(1,1,1)
     $  - HXZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,0)*HZ2(1,1)
     $  + HXZ3(0,0,0) *HZ1(1)
     $  + HXZ4(0,0,0,2)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  + 4.000000000000000d+00*HZ4(1,1,1,1)
      HYZ4(2,2,2,1) =
     $  - 1.750000000000000d+00*Zeta4
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  - HXZ1(0) *HZ3(0,0,1)
     $  + HXZ1(0) *HZ3(0,1,1)
     $  - HXZ2(0,0) *Zeta2
     $  - HXZ2(0,0) *HZ2(0,0)
     $  + HXZ2(0,0) *HZ2(0,1)
     $  + HXZ3(0,0,0) *HZ1(0)
     $  + HXZ4(0,0,0,3)
     $  - HZ2(0,0) *Zeta2
     $  - HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - HZ4(0,0,1,1)
     $  + HZ4(0,1,1,1)
      HYZ4(2,2,2,2) =
     $  + HXZ1(0) *HZ3(1,1,1)
     $  + HXZ2(0,0) *HZ2(1,1)
     $  + HXZ3(0,0,0) *HZ1(1)
     $  + HXZ4(0,0,0,0)
     $  + HZ4(1,1,1,1)
      HYZ4(2,2,2,3) =
     $  - Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(0,1,1)
     $  - HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(0,1)
     $  + HXZ3(0,0,0) *HZ1(0)
     $  + HXZ4(0,0,0,1)
     $  + HZ4(0,1,1,1)
      HYZ4(2,2,3,0) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,1,1)
     $  - HXZ2(0,0) *Zeta2
     $  - HXZ2(0,0) *HZ2(0,0)
     $  + HXZ2(0,0) *HZ2(0,1)
     $  + HXZ3(0,0,1) *HZ1(1)
     $  + HXZ4(0,0,1,2)
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,0,1,1)
     $  + 3.000000000000000d+00*HZ4(0,1,1,1)
      HYZ4(2,2,3,1) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,1)
     $  + HXZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,0)*HZ2(-1,0)
     $  + HXZ3(0,0,1) *HZ1(0)
     $  + HXZ4(0,0,1,3)
     $  - HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + HZ4( -1,0,0,0)
     $  - HZ4( -1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,1)
      HYZ4(2,2,3,2) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ2(0,0) *Zeta2
     $  + HXZ2(0,0) *HZ2(1,0)
     $  + HXZ3(0,0,1) *HZ1(1)
     $  + HXZ4(0,0,1,0)
     $  - HZ1(1) *Zeta3
     $  + HZ4(1,0,1,1)
      HYZ4(2,2,3,3) =
     $  - 2.500000000000000d-01*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,1)
     $  + HXZ2(0,0) *HZ2(0,0)
     $  + HXZ3(0,0,1) *HZ1(0)
     $  + HXZ4(0,0,1,1)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,0,1,1)
      HYZ4(2,3,0,0) =
     $  - Zeta4
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  - HXZ1(0) *HZ3(0,0,1)
     $  + HXZ1(0) *HZ3(0,1,1)
     $  + HXZ2(0,1) *HZ2(1,1)
     $  + HXZ3(0,1,2) *HZ1(1)
     $  + HXZ4(0,1,2,2)
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,1)
     $  + 3.000000000000000d+00*HZ4(0,1,1,1)
      HYZ4(2,3,0,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(1,0)
     $  + HXZ3(0,1,2) *HZ1(0)
     $  + HXZ4(0,1,2,3)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(2,3,0,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(1)*Zeta2
     $  - HXZ1(0) *HZ3(0,1,0)
     $  - HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,1) *Zeta2
     $  + HXZ3(0,1,2) *HZ1(1)
     $  + HXZ4(0,1,2,0)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(0,1,0,1)
     $  - HZ4(1,1,0,1)
      HYZ4(2,3,0,3) =
     $  - 3.500000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  + HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(0,0)
     $  + HXZ2(0,1) *HZ2(1,0)
     $  + HXZ3(0,1,2) *HZ1(0)
     $  + HXZ4(0,1,2,1)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,1)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(2,3,1,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  - HXZ1(0) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,1)
     $  - HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(0,1)
     $  + HXZ3(0,1,3) *HZ1(1)
     $  + HXZ4(0,1,3,2)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  + 4.000000000000000d+00*HZ4(-1,0,1,1)
      HYZ4(2,3,1,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  + HXZ2(0,1) *HZ2(0,0)
     $  + HXZ3(0,1,3) *HZ1(0)
     $  + HXZ4(0,1,3,3)
     $  + HZ1( -1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  + HZ4( -1,0,0,1)
      HYZ4(2,3,1,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  - HXZ1(0) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(1,-1,0)
     $  + HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(0,0)
     $  + HXZ2(0,1) *HZ2(1,0)
     $  + HXZ3(0,1,3) *HZ1(1)
     $  + HXZ4(0,1,3,0)
     $  + 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  + HZ4( -1,0,0,1)
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - HZ4(0,1,0,1)
     $  - HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,1)
      HYZ4(2,3,1,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ1(-1)*Zeta2
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  - HXZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(0,1)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(0,1)*HZ2(0,0)
     $  + HXZ3(0,1,3) *HZ1(0)
     $  + HXZ4(0,1,3,1)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
      HYZ4(2,3,2,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  - HXZ1(0) *HZ3(1,0,0)
     $  + HXZ1(0) *HZ3(1,0,1)
     $  + HXZ1(0) *HZ3(1,1,0)
     $  - HXZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,1)*HZ2(1,1)
     $  + HXZ3(0,1,0) *HZ1(1)
     $  + HXZ4(0,1,0,2)
     $  + HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(1,0,1,1)
     $  + HZ4(1,1,0,1)
      HYZ4(2,3,2,1) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,-1,0)
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - HXZ2(0,1) *Zeta2
     $  - HXZ2(0,1) *HZ2(0,0)
     $  + HXZ2(0,1) *HZ2(0,1)
     $  + HXZ3(0,1,0) *HZ1(0)
     $  + HXZ4(0,1,0,3)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + HZ4(0,1,0,1)
      HYZ4(2,3,2,2) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(1)*Zeta2
     $  + HXZ1(0) *HZ3(1,1,0)
     $  + HXZ2(0,1) *HZ2(1,1)
     $  + HXZ3(0,1,0) *HZ1(1)
     $  + HXZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,1) *Zeta2
     $  + HZ4(1,1,0,1)
      HYZ4(2,3,2,3) =
     $  + 1.750000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(0)*Zeta3
     $  + HXZ1(0) *HZ1(0)*Zeta2
     $  + HXZ1(0) *HZ3(0,1,0)
     $  - HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(0,1)
     $  + HXZ3(0,1,0) *HZ1(0)
     $  + HXZ4(0,1,0,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0,1,0,1)
      HYZ4(2,3,3,0) =
     $  + 1.500000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  - HXZ1(0) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(0)*HZ3(0,0,0)
     $  + HXZ1(0) *HZ3(0,0,1)
     $  - HXZ2(0,1) *Zeta2
     $  - HXZ2(0,1) *HZ2(0,0)
     $  + HXZ2(0,1) *HZ2(0,1)
     $  + HXZ3(0,1,1) *HZ1(1)
     $  + HXZ4(0,1,1,2)
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + 2.000000000000000d+00*HZ4(0,0,1,1)
      HYZ4(2,3,3,1) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(0)*HZ3(-1,-1,0)
     $  + HXZ1(0) *HZ3(-1,0,0)
     $  + HXZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(0,1)*HZ2(-1,0)
     $  + HXZ3(0,1,1) *HZ1(0)
     $  + HXZ4(0,1,1,3)
     $  + HZ1( -1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + HZ4( -1,0,0,1)
      HYZ4(2,3,3,2) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(0) *Zeta3
     $  + HXZ1(0) *HZ3(1,0,0)
     $  + HXZ2(0,1) *Zeta2
     $  + HXZ2(0,1) *HZ2(1,0)
     $  + HXZ3(0,1,1) *HZ1(1)
     $  + HXZ4(0,1,1,0)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,0) *Zeta2
     $  + HZ4(1,0,0,1)
      HYZ4(2,3,3,3) =
     $  - Zeta4
     $  + HXZ1(0) *HZ3(0,0,0)
     $  + HXZ2(0,1) *HZ2(0,0)
     $  + HXZ3(0,1,1) *HZ1(0)
     $  + HXZ4(0,1,1,1)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,1)
      HYZ4(3,0,0,0) =
     $  - 1.750000000000000d+00*Zeta4
     $  + HXZ1(1) *HZ3(1,1,1)
     $  + HXZ2(1,2) *HZ2(1,1)
     $  + HXZ3(1,2,2) *HZ1(1)
     $  + HXZ4(1,2,2,2)
     $  - HZ2(0,0) *Zeta2
     $  - HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - HZ4(0,0,1,1)
     $  + HZ4(0,1,1,1)
      HYZ4(3,0,0,1) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(1,0)
     $  + HXZ3(1,2,2) *HZ1(0)
     $  + HXZ4(1,2,2,3)
     $  + HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  - HZ4(0,0,1,0)
     $  + HZ4(0,1,1,0)
      HYZ4(3,0,0,2) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ2(1,2) *Zeta2
     $  + HXZ3(1,2,2) *HZ1(1)
     $  + HXZ4(1,2,2,0)
     $  + HZ1(1) *Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,1,0)
     $  + HZ4(1,0,1,0)
     $  + HZ4(1,1,1,0)
      HYZ4(3,0,0,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(0,0)
     $  + HXZ2(1,2) *HZ2(1,0)
     $  + HXZ3(1,2,2) *HZ1(0)
     $  + HXZ4(1,2,2,1)
     $  + HZ2(0,0) *Zeta2
     $  + HZ2(0,1) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + HZ4(0,1,1,0)
      HYZ4(3,0,1,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(1,0,1)
     $  - HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(0,1)
     $  + HXZ3(1,2,3) *HZ1(1)
     $  + HXZ4(1,2,3,2)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + HZ4(0,1,0,1)
      HYZ4(3,0,1,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ2(1,2) *HZ2(0,0)
     $  + HXZ3(1,2,3) *HZ1(0)
     $  + HXZ4(1,2,3,3)
     $  - HZ2(0, -1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - HZ4(0, -1,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(3,0,1,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,1,0)
     $  + HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(0,0)
     $  + HXZ2(1,2) *HZ2(1,0)
     $  + HXZ3(1,2,3) *HZ1(1)
     $  + HXZ4(1,2,3,0)
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(3,0,1,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(-1,1,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,0,0)
     $  - HXZ2(1,2) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(1,2)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(1,2)*HZ2(0,0)
     $  + HXZ3(1,2,3) *HZ1(0)
     $  + HXZ4(1,2,3,1)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - HZ4( -1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,0,0)
      HYZ4(3,0,2,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ2(1,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,2)*HZ2(1,1)
     $  + HXZ3(1,2,0) *HZ1(1)
     $  + HXZ4(1,2,0,2)
     $  - 3.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,1) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,0,1,0)
     $  + HZ4(1,1,0,0)
     $  - HZ4(1,1,0,1)
     $  - 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(3,0,2,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(0,1,0)
     $  - HXZ1(1) *HZ3(1,1,0)
     $  - HXZ2(1,2) *Zeta2
     $  - HXZ2(1,2) *HZ2(0,0)
     $  + HXZ2(1,2) *HZ2(0,1)
     $  + HXZ3(1,2,0) *HZ1(0)
     $  + HXZ4(1,2,0,3)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  + 4.000000000000000d+00*HZ4(1,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(3,0,2,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ2(1,2) *HZ2(1,1)
     $  + HXZ3(1,2,0) *HZ1(1)
     $  + HXZ4(1,2,0,0)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,1,0)
      HYZ4(3,0,2,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(1,0,0)
     $  - HXZ1(1) *HZ3(1,1,0)
     $  - HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(0,1)
     $  + HXZ3(1,2,0) *HZ1(0)
     $  + HXZ4(1,2,0,1)
     $  - HZ1(0) *Zeta3
     $  - 4.000000000000000d+00*HZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(0,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,1,1,0)
     $  - HZ4(1,0,0,0)
     $  - 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(3,0,3,0) =
     $  - 2.250000000000000d+00*Zeta4
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,0,1)
     $  - HXZ1(1) *HZ3(0,1,0)
     $  - HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  - HXZ2(1,2) *Zeta2
     $  - HXZ2(1,2) *HZ2(0,0)
     $  + HXZ2(1,2) *HZ2(0,1)
     $  + HXZ3(1,2,1) *HZ1(1)
     $  + HXZ4(1,2,1,2)
     $  - 2.000000000000000d+00*HZ1(0)*Zeta3
     $  - 3.000000000000000d+00*HZ2(0,0)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
      HYZ4(3,0,3,1) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  - HXZ1(1) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,-1,0)
     $  + HXZ2(1,2) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,2)*HZ2(-1,0)
     $  + HXZ3(1,2,1) *HZ1(0)
     $  + HXZ4(1,2,1,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + HZ4( -1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
      HYZ4(3,0,3,2) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,1,0)
     $  + HXZ2(1,2) *Zeta2
     $  + HXZ2(1,2) *HZ2(1,0)
     $  + HXZ3(1,2,1) *HZ1(1)
     $  + HXZ4(1,2,1,0)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + 3.000000000000000d+00*HZ2(0,1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + HZ4(0,1,0,0)
     $  + 3.000000000000000d+00*HZ4(0,1,1,0)
     $  + HZ4(1,0,0,0)
     $  + HZ4(1,0,1,0)
      HYZ4(3,0,3,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ2(1,2) *HZ2(0,0)
     $  + HXZ3(1,2,1) *HZ1(0)
     $  + HXZ4(1,2,1,1)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(3,1,0,0) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(0,1,1)
     $  + HXZ2(1,3) *HZ2(1,1)
     $  + HXZ3(1,3,2) *HZ1(1)
     $  + HXZ4(1,3,2,2)
     $  - HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + HZ4( -1,0,0,0)
     $  - HZ4( -1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,1)
      HYZ4(3,1,0,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(1,0)
     $  + HXZ3(1,3,2) *HZ1(0)
     $  + HXZ4(1,3,2,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
      HYZ4(3,1,0,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(1,0,0)
     $  - HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,3) *Zeta2
     $  + HXZ3(1,3,2) *HZ1(1)
     $  + HXZ4(1,3,2,0)
     $  + HZ1( -1)*Zeta3
     $  + HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  - HZ2(1,1) *Zeta2
     $  - HZ4( -1,0,1,0)
     $  - HZ4( -1,1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(1, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
     $  - 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(3,1,0,3) =
     $  + 5.500000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(0,0)
     $  + HXZ2(1,3) *HZ2(1,0)
     $  + HXZ3(1,3,2) *HZ1(0)
     $  + HXZ4(1,3,2,1)
     $  + HZ1( -1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 4.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(3,1,1,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,1)
     $  - HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(0,1)
     $  + HXZ3(1,3,3) *HZ1(1)
     $  + HXZ4(1,3,3,2)
     $  + HZ1( -1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + HZ4( -1,0,0,1)
      HYZ4(3,1,1,1) =
     $  + Zeta4
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ2(1,3) *HZ2(0,0)
     $  + HXZ3(1,3,3) *HZ1(0)
     $  + HXZ4(1,3,3,3)
     $  - HZ1( -1)*Zeta3
     $  + HZ2( -1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + HZ4( -1,-1,0,0)
     $  + HZ4( -1,0,0,0)
      HYZ4(3,1,1,2) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(0,0)
     $  + HXZ2(1,3) *HZ2(1,0)
     $  + HXZ3(1,3,3) *HZ1(1)
     $  + HXZ4(1,3,3,0)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  - HZ2(0, -1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + HZ4( -1,0,0,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + HZ4(1, -1,0,0)
      HYZ4(3,1,1,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(1)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  - HXZ2(1,3) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(1,3)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(1,3)*HZ2(0,0)
     $  + HXZ3(1,3,3) *HZ1(0)
     $  + HXZ4(1,3,3,1)
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
      HYZ4(3,1,2,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,1)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  - HXZ2(1,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,3)*HZ2(1,1)
     $  + HXZ3(1,3,0) *HZ1(1)
     $  + HXZ4(1,3,0,2)
     $  + HZ1( -1)*Zeta3
     $  + HZ1(1) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  - HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  + HZ4(0,1,0,0)
     $  - HZ4(0,1,0,1)
     $  - HZ4(0,1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(3,1,2,1) =
     $  + Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - HXZ2(1,3) *Zeta2
     $  - HXZ2(1,3) *HZ2(0,0)
     $  + HXZ2(1,3) *HZ2(0,1)
     $  + HXZ3(1,3,0) *HZ1(0)
     $  + HXZ4(1,3,0,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + HZ4( -1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
      HYZ4(3,1,2,2) =
     $  + Zeta4
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,3) *HZ2(1,1)
     $  + HXZ3(1,3,0) *HZ1(1)
     $  + HXZ4(1,3,0,0)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,1) *Zeta2
     $  + HZ4( -1,0,0,0)
     $  + HZ4( -1,0,1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,1,0)
     $  + HZ4(0, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  - HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - HZ4(0,1,1,0)
     $  + HZ4(1, -1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  - HZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(1,1,-1,0)
      HYZ4(3,1,2,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(0,1)
     $  + HXZ3(1,3,0) *HZ1(0)
     $  + HXZ4(1,3,0,1)
     $  + HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0,1) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  + 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,0,1,0)
     $  + 2.000000000000000d+00*HZ4(0,1,-1,0)
      HYZ4(3,1,3,0) =
     $  + 5.000000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,1)
     $  - HXZ2(1,3) *Zeta2
     $  - HXZ2(1,3) *HZ2(0,0)
     $  + HXZ2(1,3) *HZ2(0,1)
     $  + HXZ3(1,3,1) *HZ1(1)
     $  + HXZ4(1,3,1,2)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,0,1)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,0,1)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 3.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(3,1,3,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + HXZ2(1,3) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,3)*HZ2(-1,0)
     $  + HXZ3(1,3,1) *HZ1(0)
     $  + HXZ4(1,3,1,3)
     $  + 4.000000000000000d+00*HZ1(-1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 8.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + 4.000000000000000d+00*HZ4(-1,-1,0,0)
     $  + 4.000000000000000d+00*HZ4(-1,0,-1,0)
      HYZ4(3,1,3,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(-1,1,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,0,0)
     $  + HXZ2(1,3) *Zeta2
     $  + HXZ2(1,3) *HZ2(1,0)
     $  + HXZ3(1,3,1) *HZ1(1)
     $  + HXZ4(1,3,1,0)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(1,-1)*Zeta2
     $  + HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,-1,1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  - 4.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,0,0)
     $  - 4.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,-1,0)
      HYZ4(3,1,3,3) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  - 3.000000000000000d+00*HXZ1(1)*HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ2(1,3) *HZ2(0,0)
     $  + HXZ3(1,3,1) *HZ1(0)
     $  + HXZ4(1,3,1,1)
     $  - 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 3.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - 2.000000000000000d+00*HZ2(-1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ2(0,0) *Zeta2
     $  + 6.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  - 5.000000000000000d+00*HZ4(-1,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 3.000000000000000d+00*HZ4(-1,0,0,0)
     $  - 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
      HYZ4(3,2,0,0) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(1,1,1)
     $  + HXZ2(1,0) *HZ2(1,1)
     $  + HXZ3(1,0,2) *HZ1(1)
     $  + HXZ4(1,0,2,2)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(1,0,0,0)
     $  - HZ4(1,0,0,1)
     $  + HZ4(1,0,1,1)
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + HZ4(1,1,1,0)
      HYZ4(3,2,0,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(1,0)
     $  + HXZ3(1,0,2) *HZ1(0)
     $  + HXZ4(1,0,2,3)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - 4.000000000000000d+00*HZ4(1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(3,2,0,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ2(1,0) *Zeta2
     $  + HXZ3(1,0,2) *HZ1(1)
     $  + HXZ4(1,0,2,0)
     $  - HZ1(1) *Zeta3
     $  - HZ2(1,1) *Zeta2
     $  - HZ4(1,0,1,0)
     $  - HZ4(1,1,1,0)
      HYZ4(3,2,0,3) =
     $  + 5.500000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,1)
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(0,0)
     $  + HXZ2(1,0) *HZ2(1,0)
     $  + HXZ3(1,0,2) *HZ1(0)
     $  + HXZ4(1,0,2,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + 4.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + 2.000000000000000d+00*HZ2(1,0)*Zeta2
     $  + HZ4(0,0,1,0)
     $  + HZ4(1,0,0,0)
     $  + 2.000000000000000d+00*HZ4(1,0,1,0)
      HYZ4(3,2,1,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,1,1)
     $  - HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(0,1)
     $  + HXZ3(1,0,3) *HZ1(1)
     $  + HXZ4(1,0,3,2)
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,1)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(3,2,1,1) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,0,1)
     $  + HXZ2(1,0) *HZ2(0,0)
     $  + HXZ3(1,0,3) *HZ1(0)
     $  + HXZ4(1,0,3,3)
     $  - HZ2(0, -1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + HZ4(0,0,1,0)
      HYZ4(3,2,1,2) =
     $  + Zeta4
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(0)*Zeta2
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,0,1)
     $  - HXZ1(1) *HZ3(0,1,0)
     $  - HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(0,0)
     $  + HXZ2(1,0) *HZ2(1,0)
     $  + HXZ3(1,0,3) *HZ1(1)
     $  + HXZ4(1,0,3,0)
     $  - HZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,1,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
     $  - 2.000000000000000d+00*HZ4(0,1,-1,0)
     $  - 2.000000000000000d+00*HZ4(1,0,-1,0)
     $  + HZ4(1,0,1,0)
      HYZ4(3,2,1,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,1)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  - 3.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,1)
     $  - HXZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(1,0)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(1,0)*HZ2(0,0)
     $  + HXZ3(1,0,3) *HZ1(0)
     $  + HXZ4(1,0,3,1)
     $  - HZ1( -1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - HZ2( -1,0)*Zeta2
     $  + 2.000000000000000d+00*HZ2(0,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - 2.000000000000000d+00*HZ4(-1,0,1,0)
     $  + 4.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - 2.000000000000000d+00*HZ4(0,-1,0,0)
     $  - 4.000000000000000d+00*HZ4(0,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,0,1,0)
      HYZ4(3,2,2,0) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  + 3.000000000000000d+00*HXZ1(1)*HZ3(1,1,1)
     $  - HXZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,0)*HZ2(1,1)
     $  + HXZ3(1,0,0) *HZ1(1)
     $  + HXZ4(1,0,0,2)
     $  + HZ1(1) *Zeta3
     $  + HZ2(1,1) *Zeta2
     $  - HZ4(1,1,0,0)
     $  + HZ4(1,1,0,1)
     $  + 2.000000000000000d+00*HZ4(1,1,1,0)
      HYZ4(3,2,2,1) =
     $  + Zeta4
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  - HXZ1(1) *HZ3(0,0,1)
     $  + HXZ1(1) *HZ3(0,1,1)
     $  - HXZ2(1,0) *Zeta2
     $  - HXZ2(1,0) *HZ2(0,0)
     $  + HXZ2(1,0) *HZ2(0,1)
     $  + HXZ3(1,0,0) *HZ1(0)
     $  + HXZ4(1,0,0,3)
     $  + HZ2(0,1) *Zeta2
     $  + 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  - HZ4(0,0,1,0)
     $  + HZ4(0,1,1,0)
      HYZ4(3,2,2,2) =
     $  + Zeta4
     $  + HXZ1(1) *HZ3(1,1,1)
     $  + HXZ2(1,0) *HZ2(1,1)
     $  + HXZ3(1,0,0) *HZ1(1)
     $  + HXZ4(1,0,0,0)
     $  - HZ1(1) *Zeta3
     $  + HZ2(1,1) *Zeta2
     $  + HZ4(1,1,1,0)
      HYZ4(3,2,2,3) =
     $  - 1.250000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(0,1,1)
     $  - HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(0,1)
     $  + HXZ3(1,0,0) *HZ1(0)
     $  + HXZ4(1,0,0,1)
     $  - HZ1(0) *Zeta3
     $  + HZ2(0,1) *Zeta2
     $  + HZ4(0,1,1,0)
      HYZ4(3,2,3,0) =
     $  + 5.000000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ3(0,0,1)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,1,1)
     $  - HXZ2(1,0) *Zeta2
     $  - HXZ2(1,0) *HZ2(0,0)
     $  + HXZ2(1,0) *HZ2(0,1)
     $  + HXZ3(1,0,1) *HZ1(1)
     $  + HXZ4(1,0,1,2)
     $  - HZ2(0,0) *Zeta2
     $  - HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  + HZ4(0,1,0,1)
     $  + HZ4(0,1,1,0)
      HYZ4(3,2,3,1) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,1)
     $  + HXZ2(1,0) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,0)*HZ2(-1,0)
     $  + HXZ3(1,0,1) *HZ1(0)
     $  + HXZ4(1,0,1,3)
     $  + HZ1( -1)*Zeta3
     $  + HZ2( -1,0)*Zeta2
     $  - 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,1,0)
      HYZ4(3,2,3,2) =
     $  + 7.500000000000000d-01*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ2(1,0) *Zeta2
     $  + HXZ2(1,0) *HZ2(1,0)
     $  + HXZ3(1,0,1) *HZ1(1)
     $  + HXZ4(1,0,1,0)
     $  + 2.000000000000000d+00*HZ1(1)*Zeta3
     $  + HZ2(1,0) *Zeta2
     $  + HZ4(1,0,1,0)
      HYZ4(3,2,3,3) =
     $  + 3.000000000000000d+00*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,1)
     $  + HXZ2(1,0) *HZ2(0,0)
     $  + HXZ3(1,0,1) *HZ1(0)
     $  + HXZ4(1,0,1,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ4(0,0,1,0)
      HYZ4(3,3,0,0) =
     $  + 1.750000000000000d+00*Zeta4
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  - HXZ1(1) *HZ3(0,0,1)
     $  + HXZ1(1) *HZ3(0,1,1)
     $  + HXZ2(1,1) *HZ2(1,1)
     $  + HXZ3(1,1,2) *HZ1(1)
     $  + HXZ4(1,1,2,2)
     $  + HZ1(0) *Zeta3
     $  + 2.000000000000000d+00*HZ2(0,0)*Zeta2
     $  + 3.000000000000000d+00*HZ4(0,0,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,0,1)
     $  + HZ4(0,0,1,1)
      HYZ4(3,3,0,1) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(1,0)
     $  + HXZ3(1,1,2) *HZ1(0)
     $  + HXZ4(1,1,2,3)
     $  - HZ2(0, -1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - HZ4(0, -1,0,0)
     $  - 2.000000000000000d+00*HZ4(0,0,-1,0)
     $  + HZ4(0,0,1,0)
      HYZ4(3,3,0,2) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(1)*Zeta2
     $  - HXZ1(1) *HZ3(0,1,0)
     $  - HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,1) *Zeta2
     $  + HXZ3(1,1,2) *HZ1(1)
     $  + HXZ4(1,1,2,0)
     $  + HZ1(1) *Zeta3
     $  - HZ2(0,1) *Zeta2
     $  - HZ4(0,0,1,0)
     $  - HZ4(0,1,0,0)
     $  - HZ4(0,1,1,0)
     $  - HZ4(1,1,0,0)
      HYZ4(3,3,0,3) =
     $  + 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  + HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(0,0)
     $  + HXZ2(1,1) *HZ2(1,0)
     $  + HXZ3(1,1,2) *HZ1(0)
     $  + HXZ4(1,1,2,1)
     $  + 2.000000000000000d+00*HZ1(0)*Zeta3
     $  + HZ2(0,0) *Zeta2
     $  + HZ4(0,0,0,0)
     $  + HZ4(0,0,1,0)
      HYZ4(3,3,1,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  - HXZ1(1) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,1)
     $  - HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(0,1)
     $  + HXZ3(1,1,3) *HZ1(1)
     $  + HXZ4(1,1,3,2)
     $  + HZ1( -1)*Zeta3
     $  - 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  - HZ2( -1,0)*Zeta2
     $  - HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,1)
     $  - 2.000000000000000d+00*HZ4(-1,0,0,0)
     $  + HZ4( -1,0,0,1)
      HYZ4(3,3,1,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  + HXZ2(1,1) *HZ2(0,0)
     $  + HXZ3(1,1,3) *HZ1(0)
     $  + HXZ4(1,1,3,3)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + 4.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,0,0)
      HYZ4(3,3,1,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  - HXZ1(1) *HZ3(0,1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(1,-1,0)
     $  + HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(0,0)
     $  + HXZ2(1,1) *HZ2(1,0)
     $  + HXZ3(1,1,3) *HZ1(1)
     $  + HXZ4(1,1,3,0)
     $  - 2.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(1) *Zeta3
     $  + 2.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,1)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  + HZ2(1, -1)*Zeta2
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,-1,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  - HZ4( -1,0,1,0)
     $  + 2.000000000000000d+00*HZ4(-1,1,-1,0)
     $  + HZ4( -1,1,0,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
     $  - HZ4(0,1,0,0)
     $  + 2.000000000000000d+00*HZ4(1,-1,-1,0)
     $  + HZ4(1, -1,0,0)
      HYZ4(3,3,1,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ1(-1)*Zeta2
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  - 4.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,0,0)
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  - HXZ2(1,1) *Zeta2
     $  - 2.000000000000000d+00*HXZ2(1,1)*HZ2(-1,0)
     $  + 2.000000000000000d+00*HXZ2(1,1)*HZ2(0,0)
     $  + HXZ3(1,1,3) *HZ1(0)
     $  + HXZ4(1,1,3,1)
     $  + 3.000000000000000d+00*HZ1(-1)*Zeta3
     $  - HZ1(0) *Zeta3
     $  - 3.000000000000000d+00*HZ2(-1,-1)*Zeta2
     $  + HZ2( -1,0)*Zeta2
     $  + HZ2(0, -1)*Zeta2
     $  - 6.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + HZ4( -1,-1,0,0)
     $  + 2.000000000000000d+00*HZ4(-1,0,-1,0)
     $  + 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  + HZ4(0, -1,0,0)
      HYZ4(3,3,2,0) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  - HXZ1(1) *HZ3(1,0,0)
     $  + HXZ1(1) *HZ3(1,0,1)
     $  + HXZ1(1) *HZ3(1,1,0)
     $  - HXZ2(1,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,1)*HZ2(1,1)
     $  + HXZ3(1,1,0) *HZ1(1)
     $  + HXZ4(1,1,0,2)
     $  - 2.000000000000000d+00*HZ1(1)*Zeta3
     $  - HZ2(1,0) *Zeta2
     $  - 2.000000000000000d+00*HZ4(1,0,0,0)
     $  + HZ4(1,0,0,1)
     $  + HZ4(1,1,0,0)
      HYZ4(3,3,2,1) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,-1,0)
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - HXZ2(1,1) *Zeta2
     $  - HXZ2(1,1) *HZ2(0,0)
     $  + HXZ2(1,1) *HZ2(0,1)
     $  + HXZ3(1,1,0) *HZ1(0)
     $  + HXZ4(1,1,0,3)
     $  - HZ2(0, -1)*Zeta2
     $  - 2.000000000000000d+00*HZ4(0,-1,-1,0)
     $  - HZ4(0, -1,0,0)
     $  + HZ4(0,1,0,0)
      HYZ4(3,3,2,2) =
     $  + 2.500000000000000d-01*Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(1)*Zeta2
     $  + HXZ1(1) *HZ3(1,1,0)
     $  + HXZ2(1,1) *HZ2(1,1)
     $  + HXZ3(1,1,0) *HZ1(1)
     $  + HXZ4(1,1,0,0)
     $  - HZ1(1) *Zeta3
     $  + HZ4(1,1,0,0)
      HYZ4(3,3,2,3) =
     $  - 3.000000000000000d+00*Zeta4
     $  + 2.000000000000000d+00*HXZ1(1)*Zeta3
     $  + HXZ1(1) *HZ1(0)*Zeta2
     $  + HXZ1(1) *HZ3(0,1,0)
     $  - HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(0,1)
     $  + HXZ3(1,1,0) *HZ1(0)
     $  + HXZ4(1,1,0,1)
     $  - HZ1(0) *Zeta3
     $  + HZ4(0,1,0,0)
      HYZ4(3,3,3,0) =
     $  - Zeta4
     $  - HXZ1(1) *Zeta3
     $  - HXZ1(1) *HZ1(0)*Zeta2
     $  - 2.000000000000000d+00*HXZ1(1)*HZ3(0,0,0)
     $  + HXZ1(1) *HZ3(0,0,1)
     $  - HXZ2(1,1) *Zeta2
     $  - HXZ2(1,1) *HZ2(0,0)
     $  + HXZ2(1,1) *HZ2(0,1)
     $  + HXZ3(1,1,1) *HZ1(1)
     $  + HXZ4(1,1,1,2)
     $  - HZ1(0) *Zeta3
     $  - HZ2(0,0) *Zeta2
     $  - 3.000000000000000d+00*HZ4(0,0,0,0)
     $  + HZ4(0,0,0,1)
      HYZ4(3,3,3,1) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ1(-1)*Zeta2
     $  + 2.000000000000000d+00*HXZ1(1)*HZ3(-1,-1,0)
     $  + HXZ1(1) *HZ3(-1,0,0)
     $  + HXZ2(1,1) *Zeta2
     $  + 2.000000000000000d+00*HXZ2(1,1)*HZ2(-1,0)
     $  + HXZ3(1,1,1) *HZ1(0)
     $  + HXZ4(1,1,1,3)
     $  - HZ1( -1)*Zeta3
     $  + HZ2( -1,-1)*Zeta2
     $  + 2.000000000000000d+00*HZ4(-1,-1,-1,0)
     $  + HZ4( -1,-1,0,0)
     $  + HZ4( -1,0,0,0)
      HYZ4(3,3,3,2) =
     $  + Zeta4
     $  - HXZ1(1) *Zeta3
     $  + HXZ1(1) *HZ3(1,0,0)
     $  + HXZ2(1,1) *Zeta2
     $  + HXZ2(1,1) *HZ2(1,0)
     $  + HXZ3(1,1,1) *HZ1(1)
     $  + HXZ4(1,1,1,0)
     $  + HZ4(1,0,0,0)
      HYZ4(3,3,3,3) =
     $  + HXZ1(1) *HZ3(0,0,0)
     $  + HXZ2(1,1) *HZ2(0,0)
     $  + HXZ3(1,1,1) *HZ1(0)
     $  + HXZ4(1,1,1,1)
     $  + HZ4(0,0,0,0)
      return
      end


      function dlogone(x)
*********************************************************************
***** evaluates log(1+x) using the expansion                    *****
***** of the logarithm around 1                                 *****
***** avoids rounding errors from dlog(1d0+x)                   *****
*********************************************************************
** 1+x = (1+ep)/(1-ep), ep = x/(2+x)
** log(1+x) = log((1+x)/(1-x)) = 2*ep*(1+ep^2/3+ep^4/5+.....)
** at x= -1/2, ep = -1/3
** ep2 = 1/9, ep2^16 = 5.4 x 10^(-16)
      implicit none
      include 'types.f'
      real(dp):: dlogone
      real(dp):: x,ep,e2
      ep = x/(2.d0+x)
      e2 = ep*ep
      dlogone =2*ep*(1+e2*(1.d0/3+e2*(1.d0/5+e2*(1.d0/ 7+e2*(1.d0/ 9
     $              +e2*(1.d0/11+e2*(1.d0/13+e2*(1.d0/15+e2*(1.d0/17
     $              +e2*(1.d0/19+e2*(1.d0/21+e2*(1.d0/23+e2*(1.d0/25
     $              +e2*(1.d0/27+e2*(1.d0/29+e2*(1.d0/31+e2*(1.d0/33
     $              )))))))))))))))))
      return
      end

      subroutine makeexpar(y,z)
*********************************************************************
*****  makeexpar(y,z) creates the expansion parameters          *****
*****  for the evaluation of the coefficients                   *****
*****  tz#n:   n-th Chebyshev polynomial in 2*z-1               *****
*****      u =  4/3*H(1;z)                                      *****
*****      v =  4/3*H(-1;z)                                     *****
*****  tu#n:   n-th Chebyshev polynomial in 2*u-1               *****
*****  tv#n:   n-th Chebyshev polynomial in 2*v-1               *****
*****  zm#n:   (1-z)^n                                          *****
*****      r =  4/3*H(1;y)                                      *****
*****      s =  4/3*H(2;y)                                      *****
*****  s#n:    s^n                                              *****
*****  r#n:    s^n                                              *****
*****  b#n:    (s/z)^n                                          *****
*********************************************************************
      implicit none
      include 'types.f'
      real(dp):: dlogone
      real(dp):: y,z
      real(dp):: u,v
      real(dp)::
     $     tz01,tz02,tz03,tz04,tz05,tz06,tz07,tz08,tz09,tz10,
     $     tz11,tz12,tz13,tz14,tz15,tz16,tz17,tz18,tz19,tz20,tz21
      real(dp)::
     $     tu01,tu02,tu03,tu04,tu05,tu06,tu07,tu08,tu09,tu10,
     $     tu11,tu12,tu13,tu14,tu15,tu16,tu17,tu18,tu19,tu20,tu21
      real(dp)::
     $     tv01,tv02,tv03,tv04,tv05,tv06,tv07,tv08,tv09,tv10,
     $     tv11,tv12,tv13,tv14,tv15,tv16,tv17,tv18,tv19,tv20,tv21
      real(dp)::
     $     zm01,zm02,zm03,zm04,zm05,zm06,zm07,zm08,zm09,zm10,
     $     zm11,zm12,zm13,zm14,zm15,zm16,zm17,zm18,zm19,zm20,zm21
      real(dp)::
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      real(dp)::
     $     r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,
     $     r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22
      real(dp)::
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22

      common /tpz/
     $     tz01,tz02,tz03,tz04,tz05,tz06,tz07,tz08,tz09,tz10,
     $     tz11,tz12,tz13,tz14,tz15,tz16,tz17,tz18,tz19,tz20,tz21
      common /tpv/
     $     tv01,tv02,tv03,tv04,tv05,tv06,tv07,tv08,tv09,tv10,
     $     tv11,tv12,tv13,tv14,tv15,tv16,tv17,tv18,tv19,tv20,tv21
      common /tpu/
     $     tu01,tu02,tu03,tu04,tu05,tu06,tu07,tu08,tu09,tu10,
     $     tu11,tu12,tu13,tu14,tu15,tu16,tu17,tu18,tu19,tu20,tu21
      common /zm/
     $     zm01,zm02,zm03,zm04,zm05,zm06,zm07,zm08,zm09,zm10,
     $     zm11,zm12,zm13,zm14,zm15,zm16,zm17,zm18,zm19,zm20,zm21
      common /s/
     $     s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,
     $     s11,s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22
      common /rtdhpl/
     $     r01,r02,r03,r04,r05,r06,r07,r08,r09,r10,
     $     r11,r12,r13,r14,r15,r16,r17,r18,r19,r20,r21,r22
      common /b/
     $     b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,
     $     b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22
!$omp threadprivate(/s/,/rtdhpl/,/b/,/tpz/,/tpv/,/tpu/,/zm/)

      zm01 = 1d0-z
      zm02 = zm01*zm01
      zm03 = zm01*zm02
      zm04 = zm01*zm03
      zm05 = zm01*zm04
      zm06 = zm01*zm05
      zm07 = zm01*zm06
      zm08 = zm01*zm07
      zm09 = zm01*zm08
      zm10 = zm01*zm09
      zm11 = zm01*zm10
      zm12 = zm01*zm11
      zm13 = zm01*zm12
      zm14 = zm01*zm13
      zm15 = zm01*zm14
      zm16 = zm01*zm15
      zm17 = zm01*zm16
      zm18 = zm01*zm17
      zm19 = zm01*zm18
      zm20 = zm01*zm19
      zm21 = zm01*zm20

      tz01 = 2d0*z-1d0
      tz02 = 2d0*tz01*tz01 - 1d0
      tz03 = 2d0*tz01*tz02 - tz01
      tz04 = 2d0*tz01*tz03 - tz02
      tz05 = 2d0*tz01*tz04 - tz03
      tz06 = 2d0*tz01*tz05 - tz04
      tz07 = 2d0*tz01*tz06 - tz05
      tz08 = 2d0*tz01*tz07 - tz06
      tz09 = 2d0*tz01*tz08 - tz07
      tz10 = 2d0*tz01*tz09 - tz08
      tz11 = 2d0*tz01*tz10 - tz09
      tz12 = 2d0*tz01*tz11 - tz10
      tz13 = 2d0*tz01*tz12 - tz11
      tz14 = 2d0*tz01*tz13 - tz12
      tz15 = 2d0*tz01*tz14 - tz13
      tz16 = 2d0*tz01*tz15 - tz14
      tz17 = 2d0*tz01*tz16 - tz15
      tz18 = 2d0*tz01*tz17 - tz16
      tz19 = 2d0*tz01*tz18 - tz17
      tz20 = 2d0*tz01*tz19 - tz18
      tz21 = 2d0*tz01*tz20 - tz19

      u = -4d0/3d0*dlogone(-z)
      v =  4d0/3d0*dlogone(z)

      tu01 = 2d0*u-1d0
      tu02 = 2d0*tu01*tu01 - 1d0
      tu03 = 2d0*tu01*tu02 - tu01
      tu04 = 2d0*tu01*tu03 - tu02
      tu05 = 2d0*tu01*tu04 - tu03
      tu06 = 2d0*tu01*tu05 - tu04
      tu07 = 2d0*tu01*tu06 - tu05
      tu08 = 2d0*tu01*tu07 - tu06
      tu09 = 2d0*tu01*tu08 - tu07
      tu10 = 2d0*tu01*tu09 - tu08
      tu11 = 2d0*tu01*tu10 - tu09
      tu12 = 2d0*tu01*tu11 - tu10
      tu13 = 2d0*tu01*tu12 - tu11
      tu14 = 2d0*tu01*tu13 - tu12
      tu15 = 2d0*tu01*tu14 - tu13
      tu16 = 2d0*tu01*tu15 - tu14
      tu17 = 2d0*tu01*tu16 - tu15
      tu18 = 2d0*tu01*tu17 - tu16
      tu19 = 2d0*tu01*tu18 - tu17
      tu20 = 2d0*tu01*tu19 - tu18
      tu21 = 2d0*tu01*tu20 - tu19

      tv01 = 2d0*v-1d0
      tv02 = 2d0*tv01*tv01 - 1d0
      tv03 = 2d0*tv01*tv02 - tv01
      tv04 = 2d0*tv01*tv03 - tv02
      tv05 = 2d0*tv01*tv04 - tv03
      tv06 = 2d0*tv01*tv05 - tv04
      tv07 = 2d0*tv01*tv06 - tv05
      tv08 = 2d0*tv01*tv07 - tv06
      tv09 = 2d0*tv01*tv08 - tv07
      tv10 = 2d0*tv01*tv09 - tv08
      tv11 = 2d0*tv01*tv10 - tv09
      tv12 = 2d0*tv01*tv11 - tv10
      tv13 = 2d0*tv01*tv12 - tv11
      tv14 = 2d0*tv01*tv13 - tv12
      tv15 = 2d0*tv01*tv14 - tv13
      tv16 = 2d0*tv01*tv15 - tv14
      tv17 = 2d0*tv01*tv16 - tv15
      tv18 = 2d0*tv01*tv17 - tv16
      tv19 = 2d0*tv01*tv18 - tv17
      tv20 = 2d0*tv01*tv19 - tv18
      tv21 = 2d0*tv01*tv20 - tv19

      s01 = -4d0/3d0*dlogone(-y/(1d0-z))
      b01 = s01/z
      r01 = -4d0/3d0*dlogone(-y)

      s02 = s01*s01
      s03 = s01*s02
      s04 = s01*s03
      s05 = s01*s04
      s06 = s01*s05
      s07 = s01*s06
      s08 = s01*s07
      s09 = s01*s08
      s10 = s01*s09
      s11 = s01*s10
      s12 = s01*s11
      s13 = s01*s12
      s14 = s01*s13
      s15 = s01*s14
      s16 = s01*s15
      s17 = s01*s16
      s18 = s01*s17
      s19 = s01*s18
      s20 = s01*s19
      s21 = s01*s20
      s22 = s01*s21

      r02 = r01*r01
      r03 = r01*r02
      r04 = r01*r03
      r05 = r01*r04
      r06 = r01*r05
      r07 = r01*r06
      r08 = r01*r07
      r09 = r01*r08
      r10 = r01*r09
      r11 = r01*r10
      r12 = r01*r11
      r13 = r01*r12
      r14 = r01*r13
      r15 = r01*r14
      r16 = r01*r15
      r17 = r01*r16
      r18 = r01*r17
      r19 = r01*r18
      r20 = r01*r19
      r21 = r01*r20
      r22 = r01*r21

      b02 = b01*b01
      b03 = b01*b02
      b04 = b01*b03
      b05 = b01*b04
      b06 = b01*b05
      b07 = b01*b06
      b08 = b01*b07
      b09 = b01*b08
      b10 = b01*b09
      b11 = b01*b10
      b12 = b01*b11
      b13 = b01*b12
      b14 = b01*b13
      b15 = b01*b14
      b16 = b01*b15
      b17 = b01*b16
      b18 = b01*b17
      b19 = b01*b18
      b20 = b01*b19
      b21 = b01*b20
      b22 = b01*b21

      return
      end
