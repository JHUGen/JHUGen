******************************************************************************
**  hplog: a subroutine for the evaluation of harmonic polylogarithms
**         Version 1.0         12/07/2001
**  described in: 
**  T.Gehrmann and E.Remiddi: Numerical Evaluation of the Harmonic 
**                            Polylogarithms up to Weight 4
**                            (hep-ph/0107173; CERN-TH/2001/188)
**  the harmonic polylogarithms are defined in: 
**  E.Remiddi and J.Vermaseren: Harmonic Polylogarithms
**                            (hep-ph/9905237; Int.J.Mod.Phys. A15 (2000) 725)
**  email:
**  Thomas.Gehrmann@cern.ch and Ettore.Remiddi@bo.infn.it
**
******************************************************************************
      subroutine hplog(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                      Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
****** 
** x is the argument of the 1dHPL's (1 dimensional Harmonic PolyLogarithms) 
**   to be evaluated; 
** nw is the maximum weight of the required 1dHPL's; 
**    the maximum allowed value of nw of this implementation is 4; 
** Hc1,Hc2,Hc3,Hc4 are the complex*16 values of the 1dHPL; 
**    they must all be supplied in the arguments even if some of them 
**    are not to be evaluated; 
** Hr1,Hr2,Hr3,Hr4 are the double precision real parts of 
**    Hc1,Hc2,Hc3,Hc4; 
** Hi1,Hi2,Hi3,Hi4 are the double precision immaginary parts of 
**    Hc1,Hc2,Hc3,Hc4 divided by pi=3.114159.... 
** n1,n2 is the required range of indices, the allowed ranges are 
**    (0,1), (-1,0), (-1,1) ; 
****** 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: x,r2m1,r2p1
      complex(dp):: Hc1(n1:n2),Hc2(n1:n2,n1:n2),Hc3(n1:n2,n1:n2,n1:n2), 
     $          Hc4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hr1(n1:n2),Hr2(n1:n2,n1:n2),Hr3(n1:n2,n1:n2,n1:n2), 
     $          Hr4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2)
      integer::infilldim,infill(3)
      common /fillred/infilldim,infill 
!$omp threadprivate(/fillred/)
      real(dp),parameter::r2   = 1.4142135623730950488_dp
** check on the weight nw 
      if ( (nw.lt.1).or.(nw.gt.4) ) then 
        print*, ' illegal call of eval1dhpl with second argument', 
     $          ' (the weight) = ',nw 
        print*, ' the allowed values of the weight are 1,2,3,4 ' 
        stop
      endif
** check on the range n1:n2 
      if ( (n1.eq.-1).and.(n2.eq.0) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) = -1  
      elseif ( (n1.eq.0).and.(n2.eq.1) ) then 
        infilldim =  2 
        infill(1) =  0 
        infill(2) =  1  
      elseif ( (n1.eq.-1).and.(n2.eq.1) ) then 
        infilldim =  3 
        infill(1) =  0 
        infill(2) = -1  
        infill(3) =  1  
      else 
        print*, ' illegal call of eval1dhpl with the two last ', 
     $          'arguments = (',n1,',',n2,')' 
        print*, ' the allowed values are (-1,0), (0,1), (-1,1) ' 
        stop 
      endif 
** setting the immaginary parts equal to zero 
      call setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** looking at the range of the argument 
*      r2 = sqrt(2._dp) 
      r2m1 = r2 - 1 
      r2p1 = r2 + 1 
      if ( ( x.gt.-r2m1 ).and.( x.le.r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat0 ' 
        call eval1dhplat0(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.1_dp ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplin1 ' 
        call eval1dhplin1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2m1 ).and.( x.le.r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplat1 ' 
        call eval1dhplat1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.r2p1 ) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatinf ' 
        call eval1dhplatinf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.le.-r2p1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatminf ' 
        call eval1dhplatminf(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( x.eq.-1_dp ) then
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplinm1 ' 
        call eval1dhplinm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      elseif ( ( x.gt.-r2p1 ).and.( x.le.-r2m1) ) then 
*        print*, ' eval1dhpl:      x = ',x,', call eval1dhplatm1 ' 
        call eval1dhplatm1(x,nw,Hc1,Hc2,Hc3,Hc4, 
     $                          Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
        return 
      endif 
** 
      end 
************************************************************************ 
      subroutine eval1dhplat0(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 0-range  -(r2-1) < y <= (r2-1) 
** by direct series expansion (Bernoulli-accelerated) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp)::y
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplin1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=1 (explicit values are tabulated)
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: y
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp),parameter:: pi   = 3.14159265358979324_dp
** evaluate the irreducible 1dHPL's first 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
      call fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      if (n2.eq.0) return
** correct the ill-defined entries
      HY2(1,0) = - HY2(0,1) 
      Hi2(1,0) = 0._dp 
      H2(1,0) = cmplx(HY2(1,0),Hi2(1,0)*pi,dp)
      if ( nw.eq.2 ) return
      HY3(1,0,0) = HY3(0,0,1) 
      Hi3(1,0,0) = 0._dp 
      H3(1,0,0) = cmplx(HY3(1,0,0),Hi3(1,0,0)*pi,dp)
      if ( nw.eq.3 ) return
      HY4(1,0,0,0) = -HY4(0,0,0,1) 
      Hi4(1,0,0,0) = 0._dp 
      H4(1,0,0,0) = cmplx(HY4(1,0,0,0),Hi4(1,0,0,0)*pi,dp)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplat1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the 1-range  (r2-1) < y <= (r2+1) 
** evaluating first the H(..,r=(1-y)/(1+y)) by calling eval1dhplat0(r)  
** and then expressing H(..,y=(1-r)/(1+r)) in terms of H(..,r) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: y,r
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      real(dp):: HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      r = (1._dp-y)/(1._dp+y) 
*      print*,' eval1dhplat1: y = ',y,', r = ',r 
** the whole (-1,1) range is in general needed for any pair (n1,n2)
      call fillirr1dhplat0(r,nw,HR1,HR2,HR3,HR4,-1,1) 
** fillirr1dhplat1 takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                          HY1,HY2,HY3,HY4, 
     $                          Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatinf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the inf-range  (r2+1) < abs(y) 
** evaluating first the H(..,x=1/y) by calling eval1dhplat0(x)  
** and then expressing H(..,y=1/x) in terms of H(..,x) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: y,x
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      real(dp):: HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp),parameter:: pi   = 3.14159265358979324_dp
** the nw = 1 case 
      call fillh1(y,H1,HY1,Hi1,n1,n2) 
      if ( nw.eq.1 ) return 
** the nw > 1 case 
      x = 1._dp/y 
*      print*,' eval1dhplatinf: y = ',y,', x = ',x 
      call fillirr1dhplat0(x,nw,HX1,HX2,HX3,HX4,n1,n2) 
** fillirr1dhplatinf takes care automatically of all the immaginary 
** parts as well as of the jump across y=1 
      call fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                            HY1,HY2,HY3,HY4, 
     $                            Hi1,Hi2,Hi3,Hi4,n1,n2) 
** then the reducible 1dHPL's 
      call fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
      return 
      end 
************************************************************************ 
      subroutine eval1dhplinm1(y,nw,H1,H2,H3,H4, 
     $                           HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's for y=-1 (explicit values are tabulated)
      implicit none 
      include 'types.f'
      integer::n1,n2,nw,i,k1,k2,k3,k4,nph1,nph2,nph3,nph4
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: y
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      complex(dp):: G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      integer::infilldim,infill(3)
      common /fillred/infilldim,infill 
!$omp threadprivate(/fillred/)
      integer::istorfill(3)
      integer,parameter::nphase(-1:1)=(/-1,1,-1/) 
      real(dp),parameter:: pi   = 3.14159265358979324_dp
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplin1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = cmplx(HY1(k1),Hi1(k1)*pi,dp) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = cmplx(HY2(k1,k2),Hi2(k1,k2)*pi,dp) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = cmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi,dp) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    cmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi,dp) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      if (n1.eq.0) return
** correct the ill-defined entries
      HY2(-1,0) = - HY2(0,-1) 
      Hi2(-1,0) = Hi1(0)*HY1(-1)
      H2(-1,0) = cmplx(HY2(-1,0),Hi2(-1,0)*pi,dp)
      if ( nw.eq.2 ) return
      HY3(-1,0,0) = HY1(-1)*HY2(0,0)+HY3(0,0,-1) 
      Hi3(-1,0,0) = HY1(-1)*Hi2(0,0)-HY2(0,-1)*Hi1(0)
      H3(-1,0,0) = cmplx(HY3(-1,0,0),Hi3(-1,0,0)*pi,dp)
      if ( nw.eq.3 ) return
      HY4(-1,0,0,0) = -HY2(0,-1)*HY2(0,0)-HY4(0,0,0,-1) 
      Hi4(-1,0,0,0) = HY1(-1)*Hi3(0,0,0)+HY3(0,0,-1)*Hi1(0)
      H4(-1,0,0,0) = cmplx(HY4(-1,0,0,0),Hi4(-1,0,0,0)*pi,dp)
      return 
      end 
************************************************************************ 
      subroutine eval1dhplatm1(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range  -(r2+1) < y <= -(r2-1) 
** evaluating first the H(..,-y) by calling eval1dhplat1(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw,nph1,nph2,nph3,nph4
      real(dp):: y
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      complex(dp):: G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      integer::i,k1,k2,k3,k4,istorfill(3)
      integer::infilldim,infill(3)
      common /fillred/infilldim,infill 
!$omp threadprivate(/fillred/)
      integer,parameter::nphase(-1:1)=(/-1,1,-1/) 

      real(dp),parameter:: pi   = 3.14159265358979324_dp
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplat1(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = cmplx(HY1(k1),Hi1(k1)*pi,dp) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = cmplx(HY2(k1,k2),Hi2(k1,k2)*pi,dp) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = cmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi,dp) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    cmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi,dp) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 


      subroutine eval1dhplatminf(y,nw,H1,H2,H3,H4, 
     $                          HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates 1dhpl's in the (-1)-range y  <= -(r2+1) 
** evaluating first the H(..,-y) by calling eval1dhplatinf(-y), 
** and then expressing H(..,y) in terms of H(..,-y) 
      implicit none 
      include 'types.f'
      real(dp):: y
      integer::n1,n2,nw,k1,k2,k3,k4,i
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** additional arrays required within this routine 
      complex(dp):: G1(-n2:-n1),G2(-n2:-n1,-n2:-n1),
     $          G3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          G4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: GY1(-n2:-n1),GY2(-n2:-n1,-n2:-n1),
     $          GY3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          GY4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
      real(dp):: Gi1(-n2:-n1),Gi2(-n2:-n1,-n2:-n1),
     $          Gi3(-n2:-n1,-n2:-n1,-n2:-n1), 
     $          Gi4(-n2:-n1,-n2:-n1,-n2:-n1,-n2:-n1) 
** 
      integer::nph1,nph2,nph3,nph4
      integer::infilldim,infill(3)
      common /fillred/infilldim,infill 
!$omp threadprivate(/fillred/)
      integer::istorfill(3)
      integer,parameter::nphase(-1:1)=(/-1,1,-1/) 

      real(dp),parameter:: pi   = 3.14159265358979324_dp
*      print*,' eval1dhplatm1: y = ',y 
      if (infilldim.eq.2) then
         do i=1,2
            istorfill(i) = infill(i)
            infill(i) = -istorfill(i)
         enddo
      endif
** evaluate H(...,-y) 
      call setzero(nw,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      Gi1(0) = -1
      call eval1dhplatinf(-y,nw,G1,G2,G3,G4, 
     $                        GY1,GY2,GY3,GY4,Gi1,Gi2,Gi3,Gi4,-n2,-n1) 
      if (infilldim.eq.2) then
         do i=1,2
            infill(i) = istorfill(i)
         enddo
      endif
** fill the arrays H's 
      do k1=n1,n2 
        nph1 = nphase(k1) 
        HY1(k1) =   nph1*GY1(-k1) 
        Hi1(k1) = - nph1*Gi1(-k1) 
        H1(k1)  = cmplx(HY1(k1),Hi1(k1)*pi,dp) 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            nph2 = nph1*nphase(k2) 
            HY2(k1,k2) =   nph2*GY2(-k1,-k2) 
            Hi2(k1,k2) = - nph2*Gi2(-k1,-k2) 
            H2(k1,k2)  = cmplx(HY2(k1,k2),Hi2(k1,k2)*pi,dp) 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                nph3 = nph2*nphase(k3) 
                HY3(k1,k2,k3) =   nph3*GY3(-k1,-k2,-k3) 
                Hi3(k1,k2,k3) = - nph3*Gi3(-k1,-k2,-k3) 
                H3(k1,k2,k3)  = cmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi,dp) 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    nph4 = nph3*nphase(k4) 
                    HY4(k1,k2,k3,k4) =   nph4*GY4(-k1,-k2,-k3,-k4) 
                    Hi4(k1,k2,k3,k4) = - nph4*Gi4(-k1,-k2,-k3,-k4) 
                    H4(k1,k2,k3,k4)  = 
     $                    cmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi,dp) 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine setzero(nw,Hi1,Hi2,Hi3,Hi4,n1,n2) 
** initializes with 0 the elements of the arrays 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw,k1,k2,k3,k4
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      do k1=n1,n2 
        Hi1(k1) = 0._dp 
        if ( nw.gt.1 ) then 
          do k2=n1,n2 
            Hi2(k1,k2) = 0._dp 
            if ( nw.gt.2 ) then 
              do k3=n1,n2 
                Hi3(k1,k2,k3) = 0._dp 
                if ( nw.gt.3 ) then 
                  do k4=n1,n2 
                    Hi4(k1,k2,k3,k4) = 0._dp 
                  enddo 
                endif 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine fillred1dhpl(nw,H1,H2,H3,H4, 
     $                     HY1,HY2,HY3,HY4,Hi1,Hi2,Hi3,Hi4,n1,n2) 
* fills the reducible 1dhpl from the irreducible set
      implicit none 
      include 'types.f'
      integer::nw,n1,n2,k1,k2,k3,k4,ia,ib,ic,id,iflag
      complex(dp):: H1(n1:n2),H2(n1:n2,n1:n2),H3(n1:n2,n1:n2,n1:n2), 
     $          H4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
      integer::infilldim,infill(3)
      common /fillred/infilldim,infill 
!$omp threadprivate(/fillred/)
      real(dp),parameter::pinv = 0.318309886183790672_dp
      real(dp),parameter:: pi   = 3.14159265358979324_dp
** combining real and immaginary into the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        H2(k1,k2) = cmplx(HY2(k1,k2),Hi2(k1,k2)*pi,dp) 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            H3(k1,k2,k3) = cmplx(HY3(k1,k2,k3),Hi3(k1,k2,k3)*pi,dp) 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                H4(k1,k2,k3,k4) = 
     $               cmplx(HY4(k1,k2,k3,k4),Hi4(k1,k2,k3,k4)*pi,dp) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** evaluating the reduced HPL's 
** iflag = 0 to suppress auxiliary printings of FILLREDHPLx 
      iflag = 0 
      do ia =  1,infilldim 
      do ib = ia,infilldim 
        call FILLREDHPL2(iflag,H1,H2,n1,n2,infill(ia),infill(ib)) 
        if ( nw.gt.2 ) then 
          do ic = ib,infilldim 
            call FILLREDHPL3(iflag,H1,H2,H3,n1,n2, 
     $                          infill(ia),infill(ib),infill(ic)) 
            if ( nw.gt.3 ) then 
              do id = ic,infilldim 
                call FILLREDHPL4(iflag,H1,H2,H3,H4,n1,n2, 
     $               infill(ia),infill(ib),infill(ic),infill(id)) 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
** extractin real and immaginary parts from the complex value 
      do k1=n1,n2 
      do k2=n1,n2 
        HY2(k1,k2) =  real(H2(k1,k2),dp) 
        Hi2(k1,k2) = aimag(H2(k1,k2))*pinv 
        if ( nw.gt.2 ) then 
          do k3=n1,n2 
            HY3(k1,k2,k3) =  real(H3(k1,k2,k3),dp) 
            Hi3(k1,k2,k3) = aimag(H3(k1,k2,k3))*pinv 
            if ( nw.gt.3 ) then 
              do k4=n1,n2 
                HY4(k1,k2,k3,k4) =  real(H4(k1,k2,k3,k4),dp) 
                Hi4(k1,k2,k3,k4) = aimag(H4(k1,k2,k3,k4))*pinv 
              enddo 
            endif 
          enddo 
        endif 
      enddo 
      enddo 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL2(iflag,H1,H2,i1,i2,na,nb) 
      implicit none 
      include 'types.f'
      integer::iflag,i1,i2,na,nb
      complex(dp):: H1(i1:i2),H2(i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with ordered indices na <= nb 
*      print*,' FILLREDHPL2, iflag =',iflag 
      if ( na.eq.nb ) then 
        H2(na,na) = 1._dp/2*( H1(na) )**2 
      else 
        H2(nb,na) = + H1(na)*H1(nb) - H2(na,nb) 
        if ( iflag.eq.1 ) then 
          call printer2(na,nb) 
        endif 
      endif 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL3(iflag,H1,H2,H3,i1,i2,ia,ib,ic) 
      implicit none 
      include 'types.f'
      integer::iflag,i1,i2,ia,ib,ic,na,nb,nc
      complex(dp):: H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2) 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na) 
        H3(na,na,na) = 1._dp/6*( H1(na) )**3 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL3, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ia.eq.ib ) then 
* case (na,na,nb) 
        nb = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,na,nb) 
        endif 
        H3(na,nb,na) = + H1(na)*H2(na,nb) - 2*H3(na,na,nb) 
        H3(nb,na,na) = + 1._dp/2*H1(na)*H1(na)*H1(nb) 
     $                 - H1(na)*H2(na,nb) + H3(na,na,nb) 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL3, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb) 
        nb = ib 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nb) 
        endif 
        H3(nb,na,nb) = + H1(nb)*H2(na,nb) - 2*H3(na,nb,nb) 
        H3(nb,nb,na) = + 1._dp/2*H1(na)*H1(nb)*H1(nb) 
     $                 - H1(nb)*H2(na,nb) + H3(na,nb,nb) 
* no need to protect against ic.eq.ib 
* when arriving here all indices are different 
      else 
* case (na,nb,nc)    all indices are different 
        nb = ib 
        nc = ic 
        if ( iflag.eq.1 ) then 
          call printer3(na,nb,nc) 
          call printer3(na,nc,nb) 
        endif 
        H3(nb,na,nc) = + H1(nb)*H2(na,nc) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nb,nc,na) = + H1(na)*H2(nb,nc) 
     $                 - H1(nb)*H2(na,nc) + H3(na,nc,nb) 
        H3(nc,na,nb) = + H1(nc)*H2(na,nb) 
     $                 - H3(na,nb,nc) - H3(na,nc,nb) 
        H3(nc,nb,na) = + H1(na)*H1(nb)*H1(nc) - H1(na)*H2(nb,nc) 
     $                 - H1(nc)*H2(na,nb) + H3(na,nb,nc) 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine FILLREDHPL4(iflag,H1,H2,H3,H4,i1,i2,ia,ib,ic,id) 
      implicit none 
      include 'types.f'
      integer::iflag,i1,i2,ia,ib,ic,id,na,nb,nc,nd
      complex(dp):: H1(i1:i2),H2(i1:i2,i1:i2),H3(i1:i2,i1:i2,i1:i2),
     & H4(i1:i2,i1:i2,i1:i2,i1:i2) 
*23456789012345678901234567890123456789012345678901234567890123456789012 
* must be called with "properly ordered" indices 
* note in particular the remapping of, say, (ia,ib,ic) into 
* (na,na,nb) of ReducerTest.out 
      na = ia 
      if ( (ia.eq.ib).and.(ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,na,na,na) 
        H4(na,na,na,na) = 1._dp/24*( H1(na) )**4 
* id cannot be anymore equal to ia 
      else if ( id.eq.ia ) then 
        print*,' FILLREDHPL4, error 1, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ib.eq.ic) ) then 
* case (na,na,na,nb) 
        nb = id 
        H4(na,na,nb,na) = + H1(na)*H3(na,na,nb) - 3*H4(na,na,na,nb) 
        H4(na,nb,na,na) = + 1._dp/2*H1(na)*H1(na)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,na,nb) + 3*H4(na,na,na,nb) 
        H4(nb,na,na,na) = + 1._dp/6*H1(na)*H1(na)*H1(na)*H1(nb) 
     $                    - 1._dp/2*H1(na)*H1(na)*H2(na,nb) 
     $                    + H1(na)*H3(na,na,nb) - H4(na,na,na,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,na,nb) 
        endif 
* ic cannot be anymore equal to ia 
      else if ( ic.eq.ia ) then 
        print*,' FILLREDHPL4, error 2, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ia.eq.ib).and.(ic.eq.id) ) then 
* case (na,na,nb,nb) 
        nb = ic 
        H4(na,nb,na,nb) = + 1._dp/2*H2(na,nb)*H2(na,nb) 
     $                    - 2*H4(na,na,nb,nb) 
        H4(na,nb,nb,na) = + H1(na)*H3(na,nb,nb) 
     $                    - 1._dp/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,na,nb) = + H1(nb)*H3(na,na,nb) 
     $                    - 1._dp/2*H2(na,nb)*H2(na,nb) 
        H4(nb,na,nb,na) = + H1(na)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(na)*H3(na,nb,nb) 
     $                    - 2*H1(nb)*H3(na,na,nb) 
     $                    + 1._dp/2*H2(na,nb)*H2(na,nb) 
     $                    + 2*H4(na,na,nb,nb) 
        H4(nb,nb,na,na) = + 1._dp/4*H1(na)*H1(na)*H1(nb)*H1(nb) 
     $                    - H1(na)*H1(nb)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nb) 
     $                    + H1(nb)*H3(na,na,nb) - H4(na,na,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nb) 
        endif 
      else if ( ia.eq.ib ) then 
* case (na,na,nb,nc) 
        nb = ic 
        nc = id 
        H4(na,nb,nc,na) = + H1(na)*H3(na,nb,nc) - 2*H4(na,na,nb,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(na,nc,na,nb) = + H2(na,nb)*H2(na,nc) - 2*H4(na,na,nb,nc) 
     $                    - 2*H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(na,nc,nb,na) = + H1(na)*H3(na,nc,nb) - H2(na,nb)*H2(na,nc) 
     $                    + 2*H4(na,na,nb,nc) + H4(na,nb,na,nc) 
        H4(nb,na,na,nc) = + H1(nb)*H3(na,na,nc) - H4(na,na,nb,nc) 
     $                    - H4(na,na,nc,nb) - H4(na,nb,na,nc) 
        H4(nb,na,nc,na) = + H1(na)*H1(nb)*H2(na,nc) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nb)*H3(na,na,nc) + 2*H4(na,na,nb,nc) 
     $                    + 2*H4(na,na,nc,nb) + H4(na,nb,na,nc) 
        H4(nb,nc,na,na) = + 1._dp/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nb)*H2(na,nc) 
     $                    + H1(na)*H3(na,nc,nb) + H1(nb)*H3(na,na,nc) 
     $                    - H4(na,na,nc,nb) 
        H4(nc,na,na,nb) = + H1(nc)*H3(na,na,nb) - H2(na,nb)*H2(na,nc) 
     $                    + H4(na,na,nb,nc) + H4(na,na,nc,nb) 
     $                    + H4(na,nb,na,nc) 
        H4(nc,na,nb,na) = + H1(na)*H1(nc)*H2(na,nb) 
     $                    - H1(na)*H3(na,nb,nc) - H1(na)*H3(na,nc,nb) 
     $                    - 2*H1(nc)*H3(na,na,nb) + H2(na,nb)*H2(na,nc) 
     $                    - H4(na,nb,na,nc) 
        H4(nc,nb,na,na) = + 1._dp/2*H1(na)*H1(na)*H1(nb)*H1(nc) 
     $                    - 1._dp/2*H1(na)*H1(na)*H2(nb,nc) 
     $                    - H1(na)*H1(nc)*H2(na,nb) 
     $                    + H1(na)*H3(na,nb,nc) + H1(nc)*H3(na,na,nb) 
     $                    - H4(na,na,nb,nc)  
        if ( iflag.eq.1 ) then 
          call printer4(na,na,nb,nc) 
          call printer4(na,na,nc,nb) 
          call printer4(na,nb,na,nc) 
        endif 
* ib cannot be anymore equal to ia 
      else if ( ib.eq.ia ) then 
        print*,' FILLREDHPL4, error 3, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( (ib.eq.ic).and.(ic.eq.id) ) then 
* case (na,nb,nb,nb) 
        nb = ib 
        H4(nb,na,nb,nb) = + H1(nb)*H3(na,nb,nb) - 3*H4(na,nb,nb,nb) 
        H4(nb,nb,na,nb) = + 1._dp/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    - 2*H1(nb)*H3(na,nb,nb) + 3*H4(na,nb,nb,nb) 
        H4(nb,nb,nb,na) = + 1._dp/6*H1(na)*H1(nb)*H1(nb)*H1(nb) 
     $                    - 1._dp/2*H1(nb)*H1(nb)*H2(na,nb) 
     $                    + H1(nb)*H3(na,nb,nb) - H4(na,nb,nb,nb) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nb) 
        endif 
* id cannot be anymore equal to ib 
      else if ( id.eq.ib ) then 
        print*,' FILLREDHPL4, error 4, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ib.eq.ic ) then 
* case (na,nb,nb,nc) 
        nb = ib 
        nc = id 
        H4(nb,na,nb,nc) = + H1(nb)*H3(na,nb,nc) 
     $                    - 2*H4(na,nb,nb,nc) - H4(na,nb,nc,nb) 
        H4(nb,na,nc,nb) = + H1(nb)*H3(na,nc,nb) - H4(na,nb,nc,nb) 
     $                    - 2*H4(na,nc,nb,nb) 
        H4(nb,nb,na,nc) = + 1._dp/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H4(na,nb,nb,nc) + H4(na,nb,nc,nb) 
     $                    + H4(na,nc,nb,nb) 
        H4(nb,nb,nc,na) = + H1(na)*H3(nb,nb,nc) 
     $                    - 1._dp/2*H1(nb)*H1(nb)*H2(na,nc) 
     $                    + H1(nb)*H3(na,nc,nb) - H4(na,nc,nb,nb) 
        H4(nb,nc,na,nb) = - H1(nb)*H3(na,nb,nc) - H1(nb)*H3(na,nc,nb) 
     $                    + H2(na,nb)*H2(nb,nc) + H4(na,nb,nc,nb) 
     $                    + 2*H4(na,nc,nb,nb) 
        H4(nb,nc,nb,na) = + H1(na)*H1(nb)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nb,nc) 
     $                    + H1(nb)*H3(na,nb,nc) 
     $                    - H2(na,nb)*H2(nb,nc) - H4(na,nb,nc,nb) 
        H4(nc,na,nb,nb) = + H1(nc)*H3(na,nb,nb) - H4(na,nb,nb,nc) 
     $                    - H4(na,nb,nc,nb) - H4(na,nc,nb,nb) 
        H4(nc,nb,na,nb) = + H1(nb)*H1(nc)*H2(na,nb) 
     $                    - 2*H1(nc)*H3(na,nb,nb) 
     $                    - H2(na,nb)*H2(nb,nc) + 2*H4(na,nb,nb,nc) 
     $                    + H4(na,nb,nc,nb) 
        H4(nc,nb,nb,na) = + 1._dp/2*H1(na)*H1(nb)*H1(nb)*H1(nc) 
     $                    - H1(na)*H1(nb)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nb,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nb) + H2(na,nb)*H2(nb,nc) 
     $                    - H4(na,nb,nb,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nb,nc) 
          call printer4(na,nb,nc,nb) 
          call printer4(na,nc,nb,nb) 
        endif 
* ic cannot be anymore equal to ib 
      else if ( ic.eq.ib ) then 
        print*,' FILLREDHPL4, error 5, called with arguments ' 
        print*,'               ',ia,ib,ic,id 
        stop 
      else if ( ic.eq.id ) then 
* case (na,nb,nc,nc) 
        nb = ib 
        nc = ic 
        H4(nb,na,nc,nc) = + H1(nb)*H3(na,nc,nc) - H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) - H4(na,nc,nc,nb) 
        H4(nb,nc,na,nc) = - 2*H1(nb)*H3(na,nc,nc) + H2(na,nc)*H2(nb,nc) 
     $                    + H4(na,nc,nb,nc) + 2*H4(na,nc,nc,nb) 
        H4(nb,nc,nc,na) = + H1(na)*H3(nb,nc,nc) + H1(nb)*H3(na,nc,nc) 
     $                    - H2(na,nc)*H2(nb,nc) - H4(na,nc,nc,nb) 
        H4(nc,na,nb,nc) = + H1(nc)*H3(na,nb,nc) - 2*H4(na,nb,nc,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,na,nc,nb) = + H1(nc)*H3(na,nc,nb) - H4(na,nc,nb,nc) 
     $                    - 2*H4(na,nc,nc,nb) 
        H4(nc,nb,na,nc) = + H1(nb)*H1(nc)*H2(na,nc) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nc) + 2*H4(na,nb,nc,nc) 
     $                    + H4(na,nc,nb,nc) 
        H4(nc,nb,nc,na) = + H1(na)*H1(nc)*H2(nb,nc) 
     $                    - 2*H1(na)*H3(nb,nc,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nc) 
     $                    + H1(nc)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nc) 
     $                    - H4(na,nc,nb,nc) 
        H4(nc,nc,na,nb) = + 1._dp/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    - H1(nc)*H3(na,nb,nc) - H1(nc)*H3(na,nc,nb) 
     $                    + H4(na,nb,nc,nc) + H4(na,nc,nb,nc) 
     $                    + H4(na,nc,nc,nb) 
        H4(nc,nc,nb,na) = + 1._dp/2*H1(na)*H1(nb)*H1(nc)*H1(nc) 
     $                    - H1(na)*H1(nc)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nc) 
     $                    - 1._dp/2*H1(nc)*H1(nc)*H2(na,nb) 
     $                    + H1(nc)*H3(na,nb,nc) - H4(na,nb,nc,nc) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nc) 
          call printer4(na,nc,nb,nc) 
          call printer4(na,nc,nc,nb) 
        endif 
* no need to protect against id.eq.ic 
* when arriving here all indices are different 
      else 
* case (na,nb,nc,nd) all indices are different 
        nb = ib 
        nc = ic 
        nd = id 
        H4(nb,na,nc,nd) = + H1(nb)*H3(na,nc,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nc,nb,nd) - H4(na,nc,nd,nb) 
        H4(nb,na,nd,nc) = + H1(nb)*H3(na,nd,nc) - H4(na,nb,nd,nc) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nc,na,nd) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nd)*H2(nb,nc) + H4(na,nc,nb,nd) 
     $                    + H4(na,nc,nd,nb) + H4(na,nd,nc,nb) 
        H4(nb,nc,nd,na) = + H1(na)*H3(nb,nc,nd) + H1(nb)*H3(na,nd,nc) 
     $                    - H2(na,nd)*H2(nb,nc) - H4(na,nd,nc,nb) 
        H4(nb,nd,na,nc) = - H1(nb)*H3(na,nc,nd) - H1(nb)*H3(na,nd,nc) 
     $                    + H2(na,nc)*H2(nb,nd) + H4(na,nc,nd,nb) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nb,nd,nc,na) = + H1(na)*H3(nb,nd,nc) + H1(nb)*H3(na,nc,nd) 
     $                    - H2(na,nc)*H2(nb,nd) - H4(na,nc,nd,nb) 
        H4(nc,na,nb,nd) = + H1(nc)*H3(na,nb,nd) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nc,nb,nd) 
        H4(nc,na,nd,nb) = + H1(nc)*H3(na,nd,nb) - H4(na,nc,nd,nb) 
     $                    - H4(na,nd,nb,nc) - H4(na,nd,nc,nb) 
        H4(nc,nb,na,nd) = + H1(nb)*H1(nc)*H2(na,nd) 
     $                    - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    - H2(na,nd)*H2(nb,nc) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nd,nb,nc) 
        H4(nc,nb,nd,na) = + H1(na)*H1(nc)*H2(nb,nd) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nc)*H2(na,nd) 
     $                    + H1(nc)*H3(na,nd,nb) + H2(na,nd)*H2(nb,nc) 
     $                    - H4(na,nd,nb,nc) 
        H4(nc,nd,na,nb) = - H1(nc)*H3(na,nb,nd) - H1(nc)*H3(na,nd,nb) 
     $                    + H2(na,nb)*H2(nc,nd) + H4(na,nb,nd,nc) 
     $                    + H4(na,nd,nb,nc) + H4(na,nd,nc,nb) 
        H4(nc,nd,nb,na) = + H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nc)*H2(nb,nd) 
     $                    + H1(na)*H3(nb,nd,nc) + H1(nc)*H3(na,nb,nd) 
     $                    - H2(na,nb)*H2(nc,nd) - H4(na,nb,nd,nc) 
        H4(nd,na,nb,nc) = + H1(nd)*H3(na,nb,nc) - H4(na,nb,nc,nd) 
     $                    - H4(na,nb,nd,nc) - H4(na,nd,nb,nc) 
        H4(nd,na,nc,nb) = + H1(nd)*H3(na,nc,nb) - H4(na,nc,nb,nd) 
     $                    - H4(na,nc,nd,nb) - H4(na,nd,nc,nb) 
        H4(nd,nb,na,nc) = + H1(nb)*H1(nd)*H2(na,nc) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nc)*H2(nb,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nb,nd,nc) + H4(na,nc,nb,nd) 
        H4(nd,nb,nc,na) = + H1(na)*H1(nd)*H2(nb,nc) 
     $                    - H1(na)*H3(nb,nc,nd) - H1(na)*H3(nb,nd,nc) 
     $                    - H1(nb)*H1(nd)*H2(na,nc) 
     $                    + H1(nd)*H3(na,nc,nb) + H2(na,nc)*H2(nb,nd) 
     $                    - H4(na,nc,nb,nd) 
        H4(nd,nc,na,nb) = + H1(nc)*H1(nd)*H2(na,nb) 
     $                    - H1(nd)*H3(na,nb,nc) - H1(nd)*H3(na,nc,nb) 
     $                    - H2(na,nb)*H2(nc,nd) + H4(na,nb,nc,nd) 
     $                    + H4(na,nc,nb,nd) + H4(na,nc,nd,nb) 
        H4(nd,nc,nb,na) = + H1(na)*H1(nb)*H1(nc)*H1(nd) 
     $                    - H1(na)*H1(nb)*H2(nc,nd) 
     $                    - H1(na)*H1(nd)*H2(nb,nc) 
     $                    + H1(na)*H3(nb,nc,nd) 
     $                    - H1(nc)*H1(nd)*H2(na,nb) 
     $                    + H1(nd)*H3(na,nb,nc) 
     $                    + H2(na,nb)*H2(nc,nd) - H4(na,nb,nc,nd) 
        if ( iflag.eq.1 ) then 
          call printer4(na,nb,nc,nd) 
          call printer4(na,nb,nd,nc) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nc,nb,nd) 
          call printer4(na,nd,nb,nc) 
          call printer4(na,nd,nc,nb) 
        endif 
      endif 
*23456789012345678901234567890123456789012345678901234567890123456789012 
      return 
      end 
************************************************************************ 
      subroutine printer2(na,nb) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer3(na,nb,nc) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine printer4(na,nb,nc,nd) 

      write(11,'(''g [H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y)] = H('',$)') 
      call subprint(11,na) 
      write(11,'('','',$)') 
      call subprint(11,nb) 
      write(11,'('','',$)') 
      call subprint(11,nc) 
      write(11,'('','',$)') 
      call subprint(11,nd) 
      write(11,'('',y) ; '')') 

      write(12,'(''id H('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y) = H[('',$)') 
      call subprint(12,na) 
      write(12,'('','',$)') 
      call subprint(12,nb) 
      write(12,'('','',$)') 
      call subprint(12,nc) 
      write(12,'('','',$)') 
      call subprint(12,nd) 
      write(12,'('',y)] ; '')') 

      return 
      end 
*** 
      subroutine subprint(n,na) 
      if ( na.lt.0 ) then 
        write (n,102) na 
      else 
        write (n,101) na 
      endif 
      return 
  101 format(i1,$) 
  102 format(i2,$) 
      end 

************************************************************************
** the following routines contain th set of routines evaluating 
** irreducible 1dhpl's for various values of the arguments 
************************************************************************ 
      subroutine fillh1(y,H1,HY1,Hi1,n1,n2) 
** fillh1 evaluates the 1dhpl's of weight 1 
      implicit none 
      include 'types.f'
      integer::n1,n2
      real(dp):: y
      complex(dp):: H1(n1:n2) 
      real(dp):: HY1(n1:n2) 
      real(dp):: Hi1(n1:n2) 
      real(dp),parameter:: pi   = 3.14159265358979324_dp
      if ( n1.eq.-1) then 
        if ( y.ge.-1._dp ) then 
          HY1(-1) = log(1._dp+y) 
          Hi1(-1) = 0._dp 
        elseif ( y.lt.-1._dp ) then 
          HY1(-1) = log(-1._dp-y) 
          Hi1(-1) = 1._dp 
        endif 
        H1(-1) = cmplx(HY1(-1),pi*Hi1(-1),dp) 
      endif 
      if ( y.ge.0._dp ) then 
        HY1(0) = log(y) 
*        Hi1(0) = 0._dp 
      elseif ( y.lt.0._dp ) then 
        HY1(0) = log(-y) 
        Hi1(0) = 1._dp 
      endif 
      H1(0) = cmplx(HY1(0),pi*Hi1(0),dp) 
      if ( n2.eq.1 ) then 
        if ( y.ge.1._dp ) then 
          HY1(1) = - log(-1._dp+y) 
          Hi1(1) = 1._dp 
        elseif ( y.lt.1._dp ) then 
          HY1(1) = - log(1._dp-y) 
          Hi1(1) = 0._dp 
        endif 
        H1(1) = cmplx(HY1(1),pi*Hi1(1),dp) 
      endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat0(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluate the HPL from their power series expansions
** fillirr1dhplat0 is called by eval1dhplat0; 
** it is guaranteed that nw is in the range 1:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
** 
** for y < 0 DOES NOT evaluates the immaginary part of H(0,y) = log(y) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: y,ep,e2
      real(dp)::u,tu01,tu02,tu03,tu04,tu05,tu06,
     &     tu07,tu08,tu09,tu10,tu11,tu12
      real(dp)::v,tv01,tv02,tv03,tv04,tv05,tv06,
     & tv07,tv08,tv09,tv10,tv11,tv12
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** evaluating the required 1dHPL of weight 1 
      if ( n1.eq.-1) then 
** 1+y = (1+ep)/(1-ep), ep = y/(2+y) 
** log(1+y) = log((1+y)/(1-y)) = 2*ep*(1+ep^2/3+ep^4/5+.....) 
** at y= -(r2-1) = - 0.4142135624, ep = - 0.26120387496 
** ep2 = 0.068227464296, ep2^13 = 6.9 x 10^(-16) 
         ep = y/(2._dp+y) 
         e2 = ep*ep 
*         v = log(1._dp+y) 
         v = 2*ep*(1+e2*(1._dp/ 3+e2*(1._dp/ 5+e2*(1._dp/ 7+e2*(1._dp/ 9 
     $              +e2*(1._dp/11+e2*(1._dp/13+e2*(1._dp/15+e2*(1._dp/17    
     $              +e2*(1._dp/19+e2*(1._dp/21+e2*(1._dp/23+e2*(1._dp/25    
     $              )))))))))))))
         HY1(-1) = v 
      endif 
      if (y.ge.0._dp) then 
         HY1(0) = log(y) 
      else 
         HY1(0) = log(-y) 
** the immaginary part is evaluated in the calling routine eval1dhplat0 
**       Hi1(0) = 1_dp 
      endif 
      if ( n2.eq.1) then 
** 1-y = (1-ep)/(1+ep), ep = y/(2-y) 
         ep = y/(2._dp-y) 
         e2 = ep*ep 
*         u = - log(1._dp-y) 
         u = 2*ep*(1+e2*(1._dp/ 3+e2*(1._dp/ 5+e2*(1._dp/ 7+e2*(1._dp/ 9 
     $              +e2*(1._dp/11+e2*(1._dp/13+e2*(1._dp/15+e2*(1._dp/17    
     $              +e2*(1._dp/19+e2*(1._dp/21+e2*(1._dp/23+e2*(1._dp/25    
     $              )))))))))))))
         HY1(1) = u 
      endif 
      if ( nw.eq.1 ) return 
** from now on nw > 1 
** evaluating the Cebyshev polynomials for the expansions 
      ep = y 
      if ( n2.eq.1) then 
        tu01 = 20._dp/11._dp*u 
        tu02 = 2._dp*tu01*tu01 - 1._dp 
        tu03 = 2._dp*tu01*tu02 - tu01 
        tu04 = 2._dp*tu01*tu03 - tu02 
        tu05 = 2._dp*tu01*tu04 - tu03 
        tu06 = 2._dp*tu01*tu05 - tu04 
        tu07 = 2._dp*tu01*tu06 - tu05 
        tu08 = 2._dp*tu01*tu07 - tu06 
        tu09 = 2._dp*tu01*tu08 - tu07 
        tu10 = 2._dp*tu01*tu09 - tu08 
        tu11 = 2._dp*tu01*tu10 - tu09 
        tu12 = 2._dp*tu01*tu11 - tu10 
      endif 
      if ( n1.eq.-1 ) then 
        tv01 = 20._dp/11._dp*v 
        tv02 = 2._dp*tv01*tv01 - 1._dp 
        tv03 = 2._dp*tv01*tv02 - tv01 
        tv04 = 2._dp*tv01*tv03 - tv02 
        tv05 = 2._dp*tv01*tv04 - tv03 
        tv06 = 2._dp*tv01*tv05 - tv04 
        tv07 = 2._dp*tv01*tv06 - tv05 
        tv08 = 2._dp*tv01*tv07 - tv06 
        tv09 = 2._dp*tv01*tv08 - tv07 
        tv10 = 2._dp*tv01*tv09 - tv08 
        tv11 = 2._dp*tv01*tv10 - tv09 
        tv12 = 2._dp*tv01*tv11 - tv10 
      endif 
** evaluating the expansions 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  - 3.781250000000000e-02_dp 
     $  + 5.534574473824441e-01_dp*tu01 
     $  - 3.781250000000000e-02_dp*tu02 
     $  + 1.151036617760703e-03_dp*tu03 
     $  - 8.659502433858922e-07_dp*tu05 
     $  + 1.109042494804544e-09_dp*tu07 
     $  - 1.624415058184216e-12_dp*tu09 
     $  + 2.528376460336939e-15_dp*tu11 
**    it would be wrong to write 
**    if ( nw.eq.2 ) return 
**    because the (n1.eq.-1).and.(n2.eq.1) case is not yet complete 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 5.701592410758114e-02_dp 
     $  + 5.598247957892565e-01_dp*tu01 
     $  - 5.711486614505007e-02_dp*tu02 
     $  + 3.275603992203700e-03_dp*tu03 
     $  - 9.887255877938583e-05_dp*tu04 
     $  + 4.021153684652295e-07_dp*tu05 
     $  + 6.939288687864526e-08_dp*tu06 
     $  - 7.995347631322020e-10_dp*tu07 
     $  - 8.567978673919505e-11_dp*tu08 
     $  + 1.526387027481200e-12_dp*tu09 
     $  + 1.226899454816980e-13_dp*tu10 
     $  - 2.848614761014972e-15_dp*tu11 
     $  - 1.880542777479446e-16_dp*tu12 
      HY3(0,1,1) = 
     $  + 3.816894981500984e-02_dp 
     $  - 1.039843750000000e-02_dp*tu01 
     $  + 3.828760080995617e-02_dp*tu02 
     $  - 3.466145833333333e-03_dp*tu03 
     $  + 1.185518160084905e-04_dp*tu04 
     $  - 9.904555648775859e-08_dp*tu06 
     $  + 1.331803984518588e-10_dp*tu08 
     $  - 2.006389465106708e-13_dp*tu10 
     $  + 3.180731062055677e-16_dp*tu12 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  - 6.685228257646101e-02_dp 
     $  + 5.645990701998083e-01_dp*tu01 
     $  - 6.707912936340146e-02_dp*tu02 
     $  + 4.876429488624746e-03_dp*tu03 
     $  - 2.268732672568699e-04_dp*tu04 
     $  + 6.038494106229146e-06_dp*tu05 
     $  - 2.642577015932576e-08_dp*tu06 
     $  - 3.679843316593900e-09_dp*tu07 
     $  + 5.444046563879984e-11_dp*tu08 
     $  + 4.063821221202881e-12_dp*tu09 
     $  - 1.055985864474070e-13_dp*tu10 
     $  - 5.190408125225683e-15_dp*tu11 
     $  + 1.985464489219049e-16_dp*tu12 
      HY4(0,0,1,1) = 
     $  + 1.953236111099851e-02_dp 
     $  - 8.741612828671381e-03_dp*tu01 
     $  + 1.974116110893196e-02_dp*tu02 
     $  - 2.926558492394004e-03_dp*tu03 
     $  + 2.088576190269387e-04_dp*tu04 
     $  - 7.604351107741397e-06_dp*tu05 
     $  + 5.751031394942524e-08_dp*tu06 
     $  + 5.832253077603139e-09_dp*tu07 
     $  - 1.105713721511985e-10_dp*tu08 
     $  - 7.453416210082473e-12_dp*tu09 
     $  + 2.077758906032370e-13_dp*tu10 
     $  + 1.085601519719514e-14_dp*tu11 
     $  - 3.848312092795918e-16_dp*tu12 
      HY4(0,1,1,1) = 
     $  - 7.148925781250000e-04_dp 
     $  + 7.019393481825299e-03_dp*tu01 
     $  - 9.531901041666666e-04_dp*tu02 
     $  + 2.354287493676137e-03_dp*tu03 
     $  - 2.382975260416666e-04_dp*tu04 
     $  + 8.682904829408987e-06_dp*tu05 
     $  - 7.768198634676578e-09_dp*tu07 
     $  + 1.083130072188330e-11_dp*tu09 
     $  - 1.668810490326842e-14_dp*tu11 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 3.781250000000000e-02_dp 
     $  + 5.534574473824441e-01_dp*tv01 
     $  + 3.781250000000000e-02_dp*tv02 
     $  + 1.151036617760703e-03_dp*tv03 
     $  - 8.659502433858922e-07_dp*tv05 
     $  + 1.109042494804544e-09_dp*tv07 
     $  - 1.624415058184216e-12_dp*tv09 
     $  + 2.528376460336939e-15_dp*tv11 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 5.701592410758114e-02_dp 
     $  + 5.598247957892565e-01_dp*tv01 
     $  + 5.711486614505007e-02_dp*tv02 
     $  + 3.275603992203700e-03_dp*tv03 
     $  + 9.887255877938583e-05_dp*tv04 
     $  + 4.021153684652295e-07_dp*tv05 
     $  - 6.939288687864526e-08_dp*tv06 
     $  - 7.995347631322020e-10_dp*tv07 
     $  + 8.567978673919505e-11_dp*tv08 
     $  + 1.526387027481200e-12_dp*tv09 
     $  - 1.226899454816980e-13_dp*tv10 
     $  - 2.848614761014972e-15_dp*tv11 
     $  + 1.880542777479446e-16_dp*tv12 
      HY3(0,-1,-1) = 
     $  + 3.816894981500984e-02_dp 
     $  + 1.039843750000000e-02_dp*tv01 
     $  + 3.828760080995617e-02_dp*tv02 
     $  + 3.466145833333333e-03_dp*tv03 
     $  + 1.185518160084905e-04_dp*tv04 
     $  - 9.904555648775859e-08_dp*tv06 
     $  + 1.331803984518588e-10_dp*tv08 
     $  - 2.006389465106708e-13_dp*tv10 
     $  + 3.180731062055677e-16_dp*tv12 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 6.685228257646101e-02_dp 
     $  + 5.645990701998083e-01_dp*tv01 
     $  + 6.707912936340146e-02_dp*tv02 
     $  + 4.876429488624746e-03_dp*tv03 
     $  + 2.268732672568699e-04_dp*tv04 
     $  + 6.038494106229146e-06_dp*tv05 
     $  + 2.642577015932576e-08_dp*tv06 
     $  - 3.679843316593900e-09_dp*tv07 
     $  - 5.444046563879984e-11_dp*tv08 
     $  + 4.063821221202881e-12_dp*tv09 
     $  + 1.055985864474070e-13_dp*tv10 
     $  - 5.190408125225683e-15_dp*tv11 
     $  - 1.985464489219049e-16_dp*tv12 
      HY4(0,0,-1,-1) = 
     $  + 1.953236111099851e-02_dp 
     $  + 8.741612828671381e-03_dp*tv01 
     $  + 1.974116110893196e-02_dp*tv02 
     $  + 2.926558492394004e-03_dp*tv03 
     $  + 2.088576190269387e-04_dp*tv04 
     $  + 7.604351107741397e-06_dp*tv05 
     $  + 5.751031394942524e-08_dp*tv06 
     $  - 5.832253077603139e-09_dp*tv07 
     $  - 1.105713721511985e-10_dp*tv08 
     $  + 7.453416210082473e-12_dp*tv09 
     $  + 2.077758906032370e-13_dp*tv10 
     $  - 1.085601519719514e-14_dp*tv11 
     $  - 3.848312092795918e-16_dp*tv12 
      HY4(0,-1,-1,-1) = 
     $  + 7.148925781250000e-04_dp 
     $  + 7.019393481825299e-03_dp*tv01 
     $  + 9.531901041666666e-04_dp*tv02 
     $  + 2.354287493676137e-03_dp*tv03 
     $  + 2.382975260416666e-04_dp*tv04 
     $  + 8.682904829408987e-06_dp*tv05 
     $  - 7.768198634676578e-09_dp*tv07 
     $  + 1.083130072188330e-11_dp*tv09 
     $  - 1.668810490326842e-14_dp*tv11 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  - 2.924454241163343e-02_dp 
     $  + 3.845279287117326e-01_dp*tu01 
     $  - 2.925485694830038e-02_dp*tu02 
     $  + 1.097780471057338e-03_dp*tu03 
     $  - 1.029703135442673e-05_dp*tu04 
     $  - 7.265175511511970e-07_dp*tu05 
     $  + 1.747461299829753e-08_dp*tu06 
     $  + 7.707353556013722e-10_dp*tu07 
     $  - 3.064611747990741e-11_dp*tu08 
     $  - 8.531228176305706e-13_dp*tu09 
     $  + 5.331187822989144e-14_dp*tu10 
     $  + 8.500141365188675e-16_dp*tu11 
     $  - 6.931471805599453e-01_dp*HY1(-1) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 4.107537580582269e-02_dp 
     $  + 3.887609555197323e-01_dp*tu01 
     $  - 4.116162793629221e-02_dp*tu02 
     $  + 2.511526558054413e-03_dp*tu03 
     $  - 8.620496933228561e-05_dp*tu04 
     $  + 9.128023201466990e-07_dp*tu05 
     $  + 4.711634663963971e-08_dp*tu06 
     $  - 1.347359673414334e-09_dp*tu07 
     $  - 4.474345520888852e-11_dp*tu08 
     $  + 2.138249646727980e-12_dp*tu09 
     $  + 4.709915818801180e-14_dp*tu10 
     $  - 3.454431385666621e-15_dp*tu11 
     $  - 6.931471805599453e-01_dp*HY2(0,-1) 
      HY3(0,1,-1) = 
     $  - 4.107537580582269e-02_dp 
     $  - 3.887609555197323e-01_dp*tv01 
     $  - 4.116162793629221e-02_dp*tv02 
     $  - 2.511526558054413e-03_dp*tv03 
     $  - 8.620496933228561e-05_dp*tv04 
     $  - 9.128023201466990e-07_dp*tv05 
     $  + 4.711634663963971e-08_dp*tv06 
     $  + 1.347359673414334e-09_dp*tv07 
     $  - 4.474345520888852e-11_dp*tv08 
     $  - 2.138249646727980e-12_dp*tv09 
     $  + 4.709915818801180e-14_dp*tv10 
     $  + 3.454431385666621e-15_dp*tv11 
     $  + 6.931471805599453e-01_dp*HY2(0,1) 
      HY3(-1,-1,1) = 
     $  - 3.590863871372201e-02_dp 
     $  + 3.272029419300922e-01_dp*tu01 
     $  - 3.599657175069328e-02_dp*tu02 
     $  + 2.325685169395631e-03_dp*tu03 
     $  - 8.788997314012583e-05_dp*tu04 
     $  + 1.277831858501559e-06_dp*tu05 
     $  + 4.303730428865162e-08_dp*tu06 
     $  - 1.992295216809703e-09_dp*tu07 
     $  - 2.652932076676834e-11_dp*tu08 
     $  + 3.159865930142703e-12_dp*tu09 
     $  - 2.395589527593406e-15_dp*tu10 
     $  - 4.870947810519399e-15_dp*tu11 
     $  - 5.822405264650125e-01_dp*HY1(-1) 
     $  - 3.465735902799726e-01_dp*HY1(-1)*HY1(-1) 
      HY3(-1,1,1) = 
     $  + 3.668493142404161e-02_dp 
     $  - 1.413123104773291e-01_dp*tu01 
     $  + 3.680167312678666e-02_dp*tu02 
     $  - 3.064044728536094e-03_dp*tu03 
     $  + 1.166524199994130e-04_dp*tu04 
     $  - 8.779983417383380e-07_dp*tu05 
     $  - 8.917940330502000e-08_dp*tu06 
     $  + 1.787575622706040e-09_dp*tu07 
     $  + 1.032182649980912e-10_dp*tu08 
     $  - 3.441821872732193e-12_dp*tu09 
     $  - 1.239218730863368e-13_dp*tu10 
     $  + 6.355731482672869e-15_dp*tu11 
     $  + 1.386175839607904e-16_dp*tu12 
     $  + 2.402265069591007e-01_dp*HY1(-1)       
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  - 4.713463351559199e-02_dp 
     $  + 3.918037828258655e-01_dp*tu01 
     $  - 4.730698763577787e-02_dp*tu02 
     $  + 3.532784273601097e-03_dp*tu03 
     $  - 1.724036773635937e-04_dp*tu04 
     $  + 5.100573466380115e-06_dp*tu05 
     $  - 4.948996960052575e-08_dp*tu06 
     $  - 2.345390965359666e-09_dp*tu07 
     $  + 6.710522628543514e-11_dp*tu08 
     $  + 1.979867116023822e-12_dp*tu09 
     $  - 1.027163441987459e-13_dp*tu10 
     $  - 1.836436639605094e-15_dp*tu11 
     $  + 1.633620651699784e-16_dp*tu12 
     $  - 6.931471805599453e-01_dp*HY3(0,0,-1) 
      HY4(0,0,1,-1) = 
     $  - 4.713463351559199e-02_dp 
     $  - 3.918037828258655e-01_dp*tv01 
     $  - 4.730698763577787e-02_dp*tv02 
     $  - 3.532784273601097e-03_dp*tv03 
     $  - 1.724036773635937e-04_dp*tv04 
     $  - 5.100573466380115e-06_dp*tv05 
     $  - 4.948996960052575e-08_dp*tv06 
     $  + 2.345390965359666e-09_dp*tv07 
     $  + 6.710522628543514e-11_dp*tv08 
     $  - 1.979867116023822e-12_dp*tv09 
     $  - 1.027163441987459e-13_dp*tv10  
     $  + 1.836436639605094e-15_dp*tv11 
     $  + 1.633620651699784e-16_dp*tv12 
     $  + 6.931471805599453e-01_dp*HY3(0,0,1) 
      HY4(0,-1,0,1) = 
     $  - 5.610575179941452e-02_dp 
     $  + 4.649892609082033e-01_dp*tu01 
     $  - 5.631239161843284e-02_dp*tu02 
     $  + 4.220972769653239e-03_dp*tu03 
     $  - 2.066940413626322e-04_dp*tu04 
     $  + 6.100628682175971e-06_dp*tu05 
     $  - 5.412969106099992e-08_dp*tu06 
     $  - 3.230915912784154e-09_dp*tu07 
     $  + 9.249866333323043e-11_dp*tu08 
     $  + 2.685990764581699e-12_dp*tu09 
     $  - 1.543312114608473e-13_dp*tu10 
     $  - 2.036971731594398e-15_dp*tu11 
     $  + 2.517450307574790e-16_dp*tu12 
     $  - 8.224670334241132e-01_dp*HY2(0,-1) 
      HY4(0,-1,-1,1) = 
     $  - 4.031271939759038e-02_dp 
     $  + 3.295217254379970e-01_dp*tu01 
     $  - 4.047097737450547e-02_dp*tu02 
     $  + 3.104955391145708e-03_dp*tu03 
     $  - 1.583251510732719e-04_dp*tu04 
     $  + 5.083334568184305e-06_dp*tu05 
     $  - 6.708598619683341e-08_dp*tu06 
     $  - 1.944278941559733e-09_dp*tu07 
     $  + 8.804863765356287e-11_dp*tu08 
     $  + 9.341312729419985e-13_dp*tu09 
     $  - 1.231746977889946e-13_dp*tu10 
     $  + 3.370647349658755e-16_dp*tu11 
     $  + 1.718647072955689e-16_dp*tu12 
     $  - 5.822405264650125e-01_dp*HY2(0,-1) 
     $  - 6.931471805599453e-01_dp*HY3(0,-1,-1) 
      HY4(0,-1,1,-1) = 
     $  - 4.495764739674318e-02_dp 
     $  - 2.758514579198452e-01_dp*tv01 
     $  - 4.515130668959398e-02_dp*tv02 
     $  - 3.875995092451054e-03_dp*tv03 
     $  - 1.936768370518385e-04_dp*tv04 
     $  - 5.133195476137788e-06_dp*tv05 
     $  - 1.752786900562004e-08_dp*tv06 
     $  + 2.715518363893619e-09_dp*tv07 
     $  + 1.631155670579918e-11_dp*tv08 
     $  - 2.940721244025822e-12_dp*tv09 
     $  - 2.045219059123054e-14_dp*tv10 
     $  + 3.895696592051861e-15_dp*tv11 
     $  + 4.804530139182014e-01_dp*HY2(0,-1) 
     $  + 6.931471805599453e-01_dp*HY3(0,-1,1) 
      HY4(0,1,-1,-1) = 
     $  - 2.782664607935622e-02_dp 
     $  - 1.410831481728889e-01_dp*tv01 
     $  - 2.801876266982354e-02_dp*tv02 
     $  - 2.997894208020603e-03_dp*tv03 
     $  - 1.921960113936824e-04_dp*tv04 
     $  - 7.016503666427137e-06_dp*tv05 
     $  - 7.928257765061337e-08_dp*tv06 
     $  + 4.388745575295455e-09_dp*tv07 
     $  + 1.381107719492586e-10_dp*tv08 
     $  - 4.341921500497716e-12_dp*tv09 
     $  - 2.375364913875066e-13_dp*tv10 
     $  + 4.522044546598701e-15_dp*tv11 
     $  + 4.033357472727688e-16_dp*tv12 
     $  + 2.402265069591007e-01_dp*HY2(0,1) 
      HY4(0,-1,1,1) = 
     $  + 2.782664607935622e-02_dp 
     $  - 1.410831481728889e-01_dp*tu01 
     $  + 2.801876266982354e-02_dp*tu02 
     $  - 2.997894208020603e-03_dp*tu03 
     $  + 1.921960113936824e-04_dp*tu04 
     $  - 7.016503666427137e-06_dp*tu05 
     $  + 7.928257765061337e-08_dp*tu06 
     $  + 4.388745575295455e-09_dp*tu07 
     $  - 1.381107719492586e-10_dp*tu08 
     $  - 4.341921500497716e-12_dp*tu09 
     $  + 2.375364913875066e-13_dp*tu10 
     $  + 4.522044546598701e-15_dp*tu11 
     $  - 4.033357472727688e-16_dp*tu12 
     $  + 2.402265069591007e-01_dp*HY2(0,-1) 
      HY4(0,1,-1,1) = 
     $  + 4.495764739674318e-02_dp 
     $  - 2.758514579198452e-01_dp*tu01 
     $  + 4.515130668959398e-02_dp*tu02 
     $  - 3.875995092451054e-03_dp*tu03 
     $  + 1.936768370518385e-04_dp*tu04 
     $  - 5.133195476137788e-06_dp*tu05 
     $  + 1.752786900562004e-08_dp*tu06 
     $  + 2.715518363893619e-09_dp*tu07 
     $  - 1.631155670579918e-11_dp*tu08 
     $  - 2.940721244025822e-12_dp*tu09 
     $  + 2.045219059123054e-14_dp*tu10 
     $  + 3.895696592051861e-15_dp*tu11 
     $  + 4.804530139182014e-01_dp*HY2(0,1) 
     $  - 6.931471805599453e-01_dp*HY3(0,1,-1) 
      HY4(0,1,1,-1) = 
     $  + 4.031271939759038e-02_dp 
     $  + 3.295217254379970e-01_dp*tv01 
     $  + 4.047097737450547e-02_dp*tv02 
     $  + 3.104955391145708e-03_dp*tv03 
     $  + 1.583251510732719e-04_dp*tv04 
     $  + 5.083334568184305e-06_dp*tv05 
     $  + 6.708598619683341e-08_dp*tv06 
     $  - 1.944278941559733e-09_dp*tv07 
     $  - 8.804863765356287e-11_dp*tv08 
     $  + 9.341312729419985e-13_dp*tv09 
     $  + 1.231746977889946e-13_dp*tv10 
     $  + 3.370647349658755e-16_dp*tv11 
     $  - 1.718647072955689e-16_dp*tv12 
     $  - 5.822405264650125e-01_dp*HY2(0,1) 
     $  + 6.931471805599453e-01_dp*HY3(0,1,1) 
      HY4(-1,-1,-1,1) = 
     $  - 3.768651335815766e-02_dp 
     $  + 3.043162147119780e-01_dp*tu01 
     $  - 3.784162844891144e-02_dp*tu02 
     $  + 2.958351024362477e-03_dp*tu03 
     $  - 1.551924666783514e-04_dp*tu04 
     $  + 5.216293832777793e-06_dp*tu05 
     $  - 7.726843592398867e-08_dp*tu06 
     $  - 1.910379383726989e-09_dp*tu07 
     $  + 1.073377838077624e-10_dp*tu08 
     $  + 4.147979000313175e-13_dp*tu09 
     $  - 1.506593045440627e-13_dp*tu10 
     $  + 1.921276747438603e-15_dp*tu11 
     $  + 1.977332880766160e-16_dp*tu12 
     $  - 5.372131936080402e-01_dp*HY1(-1) 
     $  - 2.911202632325062e-01_dp*HY1(-1)*HY1(-1) 
     $  - 1.155245300933242e-01_dp*HY1(-1)*HY1(-1)*HY1(-1) 
      HY4(-1,-1,1,1) = 
     $  + 2.908893189635991e-02_dp 
     $  - 1.784837106345115e-01_dp*tu01 
     $  + 2.927117884632272e-02_dp*tu02 
     $  - 2.888221776586007e-03_dp*tu03 
     $  + 1.823501630828519e-04_dp*tu04 
     $  - 6.976883920991888e-06_dp*tu05 
     $  + 1.030302948541690e-07_dp*tu06 
     $  + 3.794029548474434e-09_dp*tu07 
     $  - 1.825184393299693e-10_dp*tu08 
     $  - 2.300206200729610e-12_dp*tu09 
     $  + 3.062629564489397e-13_dp*tu10 
     $  - 7.629393984387632e-16_dp*tu11 
     $  - 4.860728618463296e-16_dp*tu12 
     $  + 3.088253750968339e-01_dp*HY1(-1) 
     $  + 1.201132534795503e-01_dp*HY1(-1)*HY1(-1) 
      HY4(-1,1,1,1) = 
     $  - 9.029205146496301e-03_dp 
     $  + 3.753824045412342e-02_dp*tu01 
     $  - 9.240717745810759e-03_dp*tu02 
     $  + 2.351153976182453e-03_dp*tu03 
     $  - 2.115782190216214e-04_dp*tu04 
     $  + 8.486524807740892e-06_dp*tu05 
     $  - 6.547885807612483e-08_dp*tu06 
     $  - 6.934422754020238e-09_dp*tu07 
     $  + 1.405695202725693e-10_dp*tu08 
     $  + 8.329441237576153e-12_dp*tu09 
     $  - 2.790404594803712e-13_dp*tu10 
     $  - 1.024489568815216e-14_dp*tu11 
     $  + 5.256388245544115e-16_dp*tu12 
     $  - 5.550410866482157e-02_dp*HY1(-1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplat1(r,nw,HR1,HR2,HR3,HR4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for r2m1 < y < r2p1
** fillirr1dhplat1 is called by eval1dhplat1 after calling 
** fillirr1dhplat0 with argument r=(1-y)/(1+y) 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: r
      real(dp):: HR1(-1:1),HR2(-1:1,-1:1),HR3(-1:1,-1:1,-1:1), 
     $          HR4(-1:1,-1:1,-1:1,-1:1) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264e+00_dp
     $  + 6.9314718055994530e-01_dp*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      if (r.lt.0._dp) then 
      Hi2(0,1) = 
     $  - HR1( -1) 
     $  - HR1(1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942e+00_dp
     $  - 1.6449340668482264e+00_dp*HR1(-1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1( -1)*HR2(0,1)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(1)*HR1(1)
     $  - 1.6449340668482264e+00_dp*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(1)*HR1(1)
     $  - HR1(1) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + HR1(1) *HR2(0,1)
     $  - HR3( -1,-1,1)
     $  + HR3( -1,1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
     $  - HR3(0,1, -1)
     $  - HR3(0,1,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942e+00_dp
     $  - 2.4022650695910071e-01_dp*HR1(-1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR1(1)
     $  - HR1(0) *HR2(-1,1)
     $  + HR1(0) *HR2(0,-1)
     $  + HR1(0) *HR2(0,1)
     $  - 2.4022650695910071e-01_dp*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR2(-1,1)
     $  + 6.9314718055994530e-01_dp*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(0,1)
     $  - HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
     $  - HR3(0,0,1) 
     $  - HR3(0,1, -1)
      if (r.lt.0._dp) then 
      HY3(0,1,1) = HY3(0,1,1) 
     $  + 4.9348022005446793e+00_dp*HR1(-1)
     $  + 4.9348022005446793e+00_dp*HR1(1)
      Hi3(0,0,1) = 
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  + 5.0000000000000000e-01_dp*HR1(1)*HR1(1) 
      Hi3(0,1,1) = 
     $  + 6.9314718055994530e-01_dp*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR1( -1)*HR1(1)
     $  + HR1(0) *HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)
     $  + HR2( -1,1)
     $  - HR2(0, -1)
     $  - HR2(0,1) 
      endif
      endif
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381e+00_dp
     $  - 1.2020569031595942e+00_dp*HR1(-1)
     $  + 8.2246703342411321e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 1.6449340668482264e+00_dp*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1( -1)*HR3(0,1,1)
     $  + 1.6666666666666666e-01_dp*HR1(0)*HR1(1)*HR1(1)*HR1(1)
     $  - 1.2020569031595942e+00_dp*HR1(1)
     $  + 8.2246703342411321e-01_dp*HR1(1)*HR1(1)
     $  + 1.1552453009332421e-01_dp*HR1(1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  - HR1(1) *HR3(-1,1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + HR1(1) *HR3(0,1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  + HR4( -1,1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
     $  - HR4(0,1,1,1) 
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454e-01_dp
     $  - 1.2020569031595942e+00_dp*HR1(-1)
     $  + 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1( -1)*HR3(0,0,1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + 2.5000000000000000e-01_dp*HR1(0)*HR1(0)*HR1(1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(0)*HR1(1)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR1(1)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - HR1(0) *HR3(-1,1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + HR1(0) *HR3(0,1,1)
     $  - 1.2020569031595942e+00_dp*HR1(1)
     $  + 1.2011325347955035e-01_dp*HR1(1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(0,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  + HR1(1) *HR3(0,0,1)
     $  + HR1(1) *HR3(0,1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,1,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,1,1)
     $  + 2.0000000000000000e+00_dp*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
     $  - 2.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
     $  - HR4(0,0,1, -1)
     $  - HR4(0,0,1,1) 
     $  - 2.0000000000000000e+00_dp*HR4(0,1,-1,-1)
     $  - HR4(0,1, -1,1)
     $  - HR4(0,1,1, -1)
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381e+00_dp
     $  + 5.5504108664821579e-02_dp*HR1(-1)
     $  - 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(0)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  - 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + 1.6666666666666666e-01_dp*HR1(0)*HR1(0)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(0)*HR1(0)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR2(0,1)
     $  + 2.4022650695910071e-01_dp*HR1(0)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR2(0,1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  + HR1(0) *HR3(0,0,1)
     $  + HR1(0) *HR3(0,1,-1)
     $  + 5.5504108664821579e-02_dp*HR1(1)
     $  + 2.4022650695910071e-01_dp*HR2(-1,1)
     $  - 2.4022650695910071e-01_dp*HR2(0,-1)
     $  - 2.4022650695910071e-01_dp*HR2(0,1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,0,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,0,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
     $  - HR4(0,0,0,1) 
     $  - HR4(0,0,1, -1)
     $  - HR4(0,1, -1,-1)
      if (r.lt.0._dp) then 
      HY4(0,0,1,1) = HY4(0,0,1,1) 
     $  - 2.4674011002723396e+00_dp*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793e+00_dp*HR1(-1)*HR1(1)
     $  - 2.4674011002723396e+00_dp*HR1(1)*HR1(1)
      HY4(0,1,1,1) = HY4(0,1,1,1) 
     $  - 3.4205442319285582e+00_dp*HR1(-1)
     $  + 2.4674011002723396e+00_dp*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793e+00_dp*HR1(-1)*HR1(0)
     $  + 4.9348022005446793e+00_dp*HR1(-1)*HR1(1)
     $  - 4.9348022005446793e+00_dp*HR1(0)*HR1(1)
     $  - 3.4205442319285582e+00_dp*HR1(1)
     $  - 4.9348022005446793e+00_dp*HR2(-1,1)
     $  + 4.9348022005446793e+00_dp*HR2(0,-1)
     $  + 4.9348022005446793e+00_dp*HR2(0,1)
      Hi4(0,0,0,1) = 
     $  - 1.666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.000000000000000e-01_dp*HR1(-1)*HR1(1)*HR1(1) 
     $  - 1.666666666666666e-01_dp*HR1(1)*HR1(1)*HR1(1) 
      Hi4(0,0,1,1) = 
     $  - 3.465735902799726e-01_dp*HR1(-1)*HR1(-1) 
     $  + 1.666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0) 
     $  + 5.000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  - HR1( -1)*HR1(0)*HR1(1) 
     $  - 6.931471805599453e-01_dp*HR1(-1)*HR1(1) 
     $  + 5.000000000000000e-01_dp*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1( -1)*HR2(0,-1) 
     $  + HR1( -1)*HR2(0,1) 
     $  - 5.000000000000000e-01_dp*HR1(0)*HR1(1)*HR1(1) 
     $  - 3.465735902799726e-01_dp*HR1(1)*HR1(1) 
     $  - HR1(1) *HR2(-1,1) 
     $  + HR1(1) *HR2(0,-1) 
     $  + HR1(1) *HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  + HR3( -1,1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0, -1,1) 
     $  - HR3(0,1, -1) 
     $  - HR3(0,1,1) 
      Hi4(0,1,1,1) = 
     $  + 1.404707559889125e+00_dp*HR1(-1) 
     $  + 3.465735902799726e-01_dp*HR1(-1)*HR1(-1) 
     $  - 1.666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  + 5.000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0) 
     $  - 5.000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 6.931471805599453e-01_dp*HR1(-1)*HR1(0) 
     $  - 5.000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0) 
     $  + HR1( -1)*HR1(0)*HR1(1) 
     $  + 6.931471805599453e-01_dp*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - 5.000000000000000e-01_dp*HR1(0)*HR1(0)*HR1(1) 
     $  - 6.931471805599453e-01_dp*HR1(0)*HR1(1) 
     $  - HR1(0) *HR2(-1,1) 
     $  + HR1(0) *HR2(0,-1) 
     $  + HR1(0) *HR2(0,1) 
     $  + 1.404707559889125e+00_dp*HR1(1) 
     $  - 6.931471805599453e-01_dp*HR2(-1,1) 
     $  + 6.931471805599453e-01_dp*HR2(0,-1) 
     $  + 6.931471805599453e-01_dp*HR2(0,1) 
     $  - HR3( -1,-1,1) 
     $  - HR3(0, -1,-1) 
     $  - HR3(0,0, -1) 
     $  - HR3(0,0,1)  
     $  - HR3(0,1, -1) 
      endif 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
       HY2(0,-1) = 
     $  + 8.2246703342411321e-01_dp
     $  - 6.9314718055994530e-01_dp*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(1)
     $  - HR2( -1,1)
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  + 9.0154267736969571e-01_dp
     $  - 8.2246703342411321e-01_dp*HR1(-1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  - 8.2246703342411321e-01_dp*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(1)*HR1(1)
     $  + HR1(1) *HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3( -1,1,1)
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928e-01_dp
     $  - 2.4022650695910071e-01_dp*HR1(-1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR2(-1,1)
     $  - 2.4022650695910071e-01_dp*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR2(-1,1)
     $  - HR3( -1,-1,1)
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591e-01_dp
     $  - 9.0154267736969571e-01_dp*HR1(-1)
     $  + 4.1123351671205660e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 8.2246703342411321e-01_dp*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(1)*HR1(1)*HR1(1)
     $  - 9.0154267736969571e-01_dp*HR1(1)
     $  + 4.1123351671205660e-01_dp*HR1(1)*HR1(1)
     $  - 1.1552453009332421e-01_dp*HR1(1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(-1,1,1)
     $  - HR4( -1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4( -1,1,1,1)
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302e-02_dp
     $  - 1.5025711289494928e-01_dp*HR1(-1)
     $  + 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(-1,1,1)
     $  - 1.5025711289494928e-01_dp*HR1(1)
     $  + 1.2011325347955035e-01_dp*HR1(1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(-1,1)
     $  + HR1(1) *HR3(-1,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,1,1)
     $  + 2.0000000000000000e+00_dp*HR4(-1,-1,-1,1)
     $  - HR4( -1,-1,1,1)
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485e-02_dp
     $  - 5.5504108664821579e-02_dp*HR1(-1)
     $  + 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 5.5504108664821579e-02_dp*HR1(1)
     $  - 2.4022650695910071e-01_dp*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - HR4( -1,-1,-1,1)
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250e-01_dp
     $  + 6.9314718055994530e-01_dp*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)
     $  + HR1( -1)*HR1(0)
     $  - HR2(0, -1)
      if (r.lt.0._dp) then 
      Hi2(-1,1) = 
     $  - HR1( -1) 
      endif 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  + 2.4307035167006157e-01_dp
     $  - 5.8224052646501250e-01_dp*HR1(-1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  - 5.8224052646501250e-01_dp*HR1(1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      HY3(0,1,-1) = 
     $  + 5.0821521280468485e-01_dp
     $  + 1.0626935403832139e+00_dp*HR1(-1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR1(1)
     $  + 1.0626935403832139e+00_dp*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR2(0,1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000e+00_dp*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705e-02_dp
     $  - 5.8224052646501250e-01_dp*HR1(-1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      HY3(-1,1,1) = 
     $  + 5.3721319360804020e-01_dp
     $  - 2.4022650695910071e-01_dp*HR1(-1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      if (r.lt.0._dp) then 
      HY3(-1,1,1) = HY3(-1,1,1) 
     $  + 4.9348022005446793e+00_dp*HR1(-1)
      Hi3(0,-1,1) = 
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(1) 
     $  - HR2( -1,1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530e-01_dp*HR1(-1) 
     $  - 6.9314718055994530e-01_dp*HR1(1) 
      Hi3(-1,-1,1) = 
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1) 
      Hi3(-1,1,1) = 
     $  + 6.9314718055994530e-01_dp*HR1(-1) 
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1) 
     $  + HR1( -1)*HR1(0) 
     $  - HR2(0, -1) 
      endif 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932e-01_dp
     $  - 2.4307035167006157e-01_dp*HR1(-1)
     $  + 2.9112026323250625e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(1)*HR1(1)
     $  + 5.8224052646501250e-01_dp*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,-1,1)
     $  - HR1(0) *HR1(1)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(-1,1,1)
     $  - 2.4307035167006157e-01_dp*HR1(1)
     $  + 2.9112026323250625e-01_dp*HR1(1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(-1,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,1,1)
     $  - 2.0000000000000000e+00_dp*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0, -1,1,1)
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438e-01_dp
     $  - 5.0821521280468485e-01_dp*HR1(-1)
     $  - 5.3134677019160696e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139e+00_dp*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(-1,1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - HR1( -1)*HR3(-1,1,1)
     $  - 2.0000000000000000e+00_dp*HR1(-1)*HR3(0,-1,-1)
     $  - HR1( -1)*HR3(0,-1,1)
     $  - HR1( -1)*HR3(0,1,-1)
     $  - 3.4657359027997265e-01_dp*HR1(0)*HR1(1)*HR1(1)
     $  - 5.0821521280468485e-01_dp*HR1(1)
     $  - 5.3134677019160696e-01_dp*HR1(1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(-1,1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(0,1)
     $  - HR1(1) *HR3(-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(1)*HR3(0,-1,-1)
     $  - HR1(1) *HR3(0,-1,1)
     $  - HR1(1) *HR3(0,1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,1,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,1,1)
     $  - 2.0000000000000000e+00_dp*HR4(-1,-1,-1,1)
     $  + HR4( -1,-1,1,1)
     $  + 3.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000e+00_dp*HR4(0,-1,-1,1)
     $  + 2.0000000000000000e+00_dp*HR4(0,-1,1,-1)
     $  + HR4(0, -1,1,1)
     $  + 2.0000000000000000e+00_dp*HR4(0,1,-1,-1)
     $  + HR4(0,1, -1,1)
     $  + HR4(0,1,1, -1)
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841e-01_dp
     $  - 3.8889584616810632e-01_dp*HR1(-1)
     $  + 8.2246703342411321e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 1.6449340668482264e+00_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(-1,1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - HR1( -1)*HR1(1)*HR2(0,1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(-1)*HR3(-1,-1,1)
     $  + 2.0000000000000000e+00_dp*HR1(-1)*HR3(-1,1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,1,-1)
     $  + HR1(0) *HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000e+00_dp*HR1(0)*HR3(-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(0)*HR3(-1,1,1)
     $  - 3.8889584616810632e-01_dp*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(-1,1)
     $  + 2.0000000000000000e+00_dp*HR1(1)*HR3(-1,-1,1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,1,-1)
     $  - 1.6449340668482264e+00_dp*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR2(-1,1)*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + HR2( -1,1)*HR2(0,1)
     $  + 1.3862943611198906e+00_dp*HR3(-1,-1,1)
     $  - 1.3862943611198906e+00_dp*HR3(-1,1,1)
     $  + 4.0000000000000000e+00_dp*HR4(-1,-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR4(-1,-1,1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0,1, -1,-1)
     $  - HR4(0,1, -1,1)
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913e-02_dp
     $  - 9.4753004230127705e-02_dp*HR1(-1)
     $  + 2.9112026323250625e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR1(0)*HR2(-1,1)
     $  + 5.8224052646501250e-01_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1(0) *HR3(-1,-1,1)
     $  - 9.4753004230127705e-02_dp*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  - 5.8224052646501250e-01_dp*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0, -1,-1,1)
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652e-02_dp
     $  - 2.1407237086670622e-01_dp*HR1(-1)
     $  - 5.3134677019160696e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  - 1.0626935403832139e+00_dp*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(-1)*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR2(-1,1)
     $  - 2.1407237086670622e-01_dp*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000e+00_dp*HR1(1)*HR3(0,-1,-1)
     $  + 1.0626935403832139e+00_dp*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  + HR4( -1,-1,-1,1)
     $  + 3.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000e+00_dp*HR4(0,-1,-1,1)
     $  + HR4(0, -1,1,-1)
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084e-01_dp
     $  + 4.7533770109129867e-01_dp*HR1(-1)
     $  - 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(0)
     $  - 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  - HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + 2.4022650695910071e-01_dp*HR1(0)*HR1(1)
     $  + 4.7533770109129867e-01_dp*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + 2.4022650695910071e-01_dp*HR2(-1,1)
     $  - 2.4022650695910071e-01_dp*HR2(0,-1)
     $  - 2.4022650695910071e-01_dp*HR2(0,1)
     $  + 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + 1.3862943611198906e+00_dp*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  + HR4( -1,-1,-1,1)
     $  - 3.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577e-02_dp
     $  - 5.3721319360804020e-01_dp*HR1(-1)
     $  + 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR1(1)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR2(-1,1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,-1,1)
     $  - 5.3721319360804020e-01_dp*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 2.4022650695910071e-01_dp*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  - HR4( -1,-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  - HR4(0, -1,-1,1)
     $  - HR4(0, -1,1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0, -1,1)
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524e-01_dp
     $  + 1.4780047665430420e+00_dp*HR1(-1)
     $  - 2.9112026323250625e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + 5.8224052646501250e-01_dp*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR1(0)*HR2(-1,1)
     $  + HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 5.8224052646501250e-01_dp*HR1(-1)*HR1(1)
     $  + HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(-1)*HR3(0,-1,-1)
     $  - 2.0000000000000000e+00_dp*HR1(-1)*HR3(0,0,-1)
     $  + 5.8224052646501250e-01_dp*HR1(0)*HR1(1)
     $  + HR1(0) *HR1(1)*HR2(0,-1)
     $  - HR1(0) *HR3(-1,-1,1)
     $  - 2.0000000000000000e+00_dp*HR1(0)*HR3(0,-1,-1)
     $  - HR1(0) *HR3(0,-1,1)
     $  - HR1(0) *HR3(0,1,-1)
     $  + 1.4780047665430420e+00_dp*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(1)*HR2(0,-1)
     $  - 2.0000000000000000e+00_dp*HR1(1)*HR3(0,-1,-1)
     $  - 2.0000000000000000e+00_dp*HR1(1)*HR3(0,0,-1)
     $  + 5.8224052646501250e-01_dp*HR2(-1,1)
     $  - HR2( -1,1)*HR2(0,-1)
     $  - 5.8224052646501250e-01_dp*HR2(0,-1)
     $  + 5.0000000000000000e-01_dp*HR2(0,-1)*HR2(0,-1)
     $  + HR2(0, -1)*HR2(0,1)
     $  - 5.8224052646501250e-01_dp*HR2(0,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - 1.3862943611198906e+00_dp*HR3(0,-1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  + 4.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  + 2.0000000000000000e+00_dp*HR4(0,-1,-1,1)
     $  - HR4(0, -1,0,1)
     $  + HR4(0, -1,1,-1)
     $  + 2.0000000000000000e+00_dp*HR4(0,0,-1,-1)
     $  + HR4(0,1, -1,-1)
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519e-01_dp
     $  - 1.1073038989294665e+00_dp*HR1(-1)
     $  + 5.3134677019160696e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(-1,1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  - 1.0626935403832139e+00_dp*HR1(-1)*HR1(0)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)*HR1(1)
     $  + 1.0626935403832139e+00_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(1)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR2(-1,1)
     $  + HR1( -1)*HR3(-1,-1,1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  - 3.4657359027997265e-01_dp*HR1(0)*HR1(0)*HR1(1)
     $  - 1.0626935403832139e+00_dp*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR2(-1,1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR2(0,1)
     $  - 1.1073038989294665e+00_dp*HR1(1)
     $  + HR1(1) *HR3(0,-1,-1)
     $  + HR1(1) *HR3(0,0,-1)
     $  - 1.0626935403832139e+00_dp*HR2(-1,1)
     $  + HR2( -1,1)*HR2(0,-1)
     $  + 1.0626935403832139e+00_dp*HR2(0,-1)
     $  - 5.0000000000000000e-01_dp*HR2(0,-1)*HR2(0,-1)
     $  - HR2(0, -1)*HR2(0,1)
     $  + 1.0626935403832139e+00_dp*HR2(0,1)
     $  - 6.9314718055994530e-01_dp*HR3(-1,-1,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,0,-1)
     $  - 6.9314718055994530e-01_dp*HR3(0,0,1)
     $  - 6.9314718055994530e-01_dp*HR3(0,1,-1)
     $  - HR4( -1,-1,-1,1)
     $  - HR4(0, -1,-1,1)
     $  + HR4(0, -1,0,1)
     $  + HR4(0,0, -1,1)
     $  + HR4(0,0,1, -1)
     $  + HR4(0,1, -1,-1)
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008e-02_dp
     $  - 9.4753004230127705e-02_dp*HR1(-1)
     $  + 2.9112026323250625e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  - HR4(0, -1,-1,-1)
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251e-02_dp
     $  - 5.3721319360804020e-01_dp*HR1(-1)
     $  + 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  - HR1( -1)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR2(0,-1)
     $  + HR1( -1)*HR3(0,-1,-1)
     $  + HR1( -1)*HR3(0,0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  - 2.0000000000000000e+00_dp*HR4(0,-1,-1,-1)
     $  - HR4(0,0, -1,-1)
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938e-01_dp
     $  + 5.5504108664821579e-02_dp*HR1(-1)
     $  - 1.2011325347955035e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.1552453009332421e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 4.1666666666666666e-02_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)*HR1(0)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 2.5000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)*HR1(0)
     $  + 2.4022650695910071e-01_dp*HR1(-1)*HR1(0)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(0)*HR1(0)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(0)*HR1(0)*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR1(0)*HR2(0,-1)
     $  + HR1(0) *HR3(0,-1,-1)
     $  + HR1(0) *HR3(0,0,-1)
     $  - 2.4022650695910071e-01_dp*HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,-1,-1)
     $  + 6.9314718055994530e-01_dp*HR3(0,0,-1)
     $  - HR4(0, -1,-1,-1)
     $  - HR4(0,0, -1,-1)
     $  - HR4(0,0,0, -1)
      if (r.lt.0._dp) then 
      HY4(0,-1,1,1) = HY4(0,-1,1,1) 
     $  - 2.4674011002723396e+00_dp*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793e+00_dp*HR1(-1)*HR1(1)
     $  + 4.9348022005446793e+00_dp*HR2(-1,1)
      HY4(0,1,1,-1) = HY4(0,1,1,-1) 
     $  + 3.4205442319285582e+00_dp*HR1(-1)
     $  + 3.4205442319285582e+00_dp*HR1(1)
      HY4(-1,-1,1,1) = HY4(-1,-1,1,1) 
     $  - 2.4674011002723396e+00_dp*HR1(-1)*HR1(-1) 
      HY4(-1,1,1,1) = HY4(-1,1,1,1)
     $  - 3.4205442319285582e+00_dp*HR1(-1)
     $  + 2.4674011002723396e+00_dp*HR1(-1)*HR1(-1)
     $  - 4.9348022005446793e+00_dp*HR1(-1)*HR1(0)
     $  + 4.9348022005446793e+00_dp*HR2(0,-1)
      Hi4(0,0,-1,1) = 
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(1)*HR1(1) 
     $  + HR1(1) *HR2(-1,1) 
     $  + HR3( -1,-1,1) 
     $  - HR3( -1,1,1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(1) 
     $  + 3.4657359027997265e-01_dp*HR1(1)*HR1(1) 
      Hi4(0,-1,0,1) = 
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR1(1) *HR2(-1,1) 
     $  - 2.0000000000000000e+00_dp*HR3(-1,-1,1) 
     $  + 2.0000000000000000e+00_dp*HR3(-1,1,1) 
      Hi4(0,-1,-1,1) = 
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1) 
     $  + HR1( -1)*HR2(-1,1) 
     $  - HR3( -1,-1,1) 
      Hi4(0,-1,1,-1) = 
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1) 
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(1) 
     $  - 6.9314718055994530e-01_dp*HR2(-1,1) 
      Hi4(0,1,-1,-1) =  
     $  - 2.4022650695910071e-01_dp*HR1(-1) 
     $  - 2.4022650695910071e-01_dp*HR1(1) 
      Hi4(0,-1,1,1) = 
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR1(0)*HR1(1)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  + HR1( -1)*HR2(0,-1)
     $  + HR1(0) *HR2(-1,1)
     $  + HR1(1) *HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(-1,1)
     $  + HR3( -1,-1,1)
     $  - HR3(0, -1,-1)
     $  - HR3(0, -1,1)
      Hi4(0,1,-1,1) = 
     $  - 5.8224052646501250e-01_dp*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(1)
     $  - HR1( -1)*HR2(-1,1)
     $  - HR1( -1)*HR2(0,-1)
     $  - 5.8224052646501250e-01_dp*HR1(1)
     $  - HR1(1) *HR2(0,-1)
     $  + HR3( -1,-1,1)
     $  + 2.0000000000000000e+00_dp*HR3(0,-1,-1)
     $  + HR3(0, -1,1)
     $  + HR3(0,1, -1)
      Hi4(0,1,1,-1) = 
     $  + 1.0626935403832139e+00_dp*HR1(-1)
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR1(0)*HR1(1)
     $  + 1.0626935403832139e+00_dp*HR1(1)
     $  + 6.9314718055994530e-01_dp*HR2(-1,1)
     $  - 6.9314718055994530e-01_dp*HR2(0,-1)
     $  - 6.9314718055994530e-01_dp*HR2(0,1)
      Hi4(-1,-1,-1,1) = 
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1) 
      Hi4(-1,-1,1,1) = 
     $  - 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  + 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  + HR1( -1)*HR2(0,-1)
     $  - HR3(0, -1,-1)
      Hi4(-1,1,1,1) = 
     $  + 1.4047075598891257e+00_dp*HR1(-1)
     $  + 3.4657359027997265e-01_dp*HR1(-1)*HR1(-1)
     $  - 1.6666666666666666e-01_dp*HR1(-1)*HR1(-1)*HR1(-1)
     $  + 5.0000000000000000e-01_dp*HR1(-1)*HR1(-1)*HR1(0)
     $  - 6.9314718055994530e-01_dp*HR1(-1)*HR1(0)
     $  - 5.0000000000000000e-01_dp*HR1(-1)*HR1(0)*HR1(0)
     $  + HR1(0) *HR2(0,-1)
     $  + 6.9314718055994530e-01_dp*HR2(0,-1)
     $  - HR3(0, -1,-1)
     $  - HR3(0,0, -1)
      endif  
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplatinf(x,nw,HX1,HX2,HX3,HX4, 
     $                                HY1,HY2,HY3,HY4, 
     $                                Hi1,Hi2,Hi3,Hi4,n1,n2) 
** evaluates the HPL for y > r2p1
** fillirr1dhplatinf is called by eval1dhplatinf after calling 
** fillirr1dhplat0 with argument r=1/y 
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: x
      real(dp):: HX1(n1:n2),HX2(n1:n2,n1:n2),HX3(n1:n2,n1:n2,n1:n2), 
     $          HX4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
      real(dp):: Hi1(n1:n2),Hi2(n1:n2,n1:n2),Hi3(n1:n2,n1:n2,n1:n2), 
     $          Hi4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 3.2898681336964528e+00_dp 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      Hi2(0,1) = 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,0,1) = 
     $  - 3.2898681336964528e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      HY3(0,1,1) = 
     $  + 1.2020569031595942e+00_dp 
     $  + 4.9348022005446793e+00_dp*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,0,1) = 
     $  + 5.000000000000000e-01_dp*HX1(0)*HX1(0) 
      Hi3(0,1,1) = 
     $  + 1.6449340668482264e+00_dp 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,1) = 
     $  + 2.1646464674222763e+00_dp 
     $  + 1.6449340668482264e+00_dp*HX1(0)*HX1(0) 
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,1) = 
     $  + 2.1646464674222763e+00_dp 
     $  - 1.2020569031595942e+00_dp*HX1(0) 
     $  - 2.4674011002723396e+00_dp*HX1(0)*HX1(0) 
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - 2.0000000000000000e+00_dp*HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
      HY4(0,1,1,1) = 
     $  - 5.1410353601279064e+00_dp 
     $  + 2.4674011002723396e+00_dp*HX1(0)*HX1(0) 
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,1) 
     $  + HX1(0) *HX3(0,0,1) 
     $  - HX1(0) *HX3(0,1,1) 
     $  + 4.9348022005446793e+00_dp*HX2(0,1) 
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,0,1) = 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
      Hi4(0,0,1,1) = 
     $  - 1.2020569031595942e+00_dp 
     $  - 1.6449340668482264e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0,1) 
      Hi4(0,1,1,1) = 
     $  + 1.6449340668482264e+00_dp*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 1.6449340668482264e+00_dp 
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      if ( nw.gt.2 ) then 
      HY3(0,0,-1) = 
     $  - 1.6449340668482264e+00_dp*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      HY3(0,-1,-1) = 
     $  + 1.2020569031595942e+00_dp 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,0,-1) = 
     $  + 1.8940656589944918e+00_dp 
     $  + 8.2246703342411321e-01_dp*HX1(0)*HX1(0) 
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX4(0,0,0, -1) 
      HY4(0,0,-1,-1) = 
     $  - 1.8940656589944918e+00_dp 
     $  - 1.2020569031595942e+00_dp*HX1(0) 
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX3(0,0,-1) 
     $  + HX4(0,0, -1,-1) 
     $  + 2.0000000000000000e+00_dp*HX4(0,0,0,-1) 
      HY4(0,-1,-1,-1) = 
     $  + 1.0823232337111381e+00_dp 
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0) 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1) 
     $  + HX1(0) *HX3(0,-1,-1) 
     $  + HX1(0) *HX3(0,0,-1) 
     $  - HX4(0, -1,-1,-1) 
     $  - HX4(0,0, -1,-1) 
     $  - HX4(0,0,0, -1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 2.4674011002723396e+00_dp 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      Hi2(-1,1) = 
     $  - 6.9314718055994530e-01_dp 
     $  + HX1( -1) 
     $  - HX1(0) 
      if ( nw.gt.2 ) then 
      HY3(0,-1,1) = 
     $  - 2.5190015545588625e+00_dp 
     $  - 2.4674011002723396e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000e+00_dp*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      HY3(0,1,-1) = 
     $  + 4.3220869092982539e+00_dp 
     $  + 2.4674011002723396e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000e+00_dp*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      HY3(-1,-1,1) = 
     $  - 2.7620719062289241e+00_dp 
     $  + 2.4674011002723396e+00_dp*HX1(-1) 
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX1(0) 
     $  - 5.0000000000000000e-01_dp*HX1(-1)*HX1(0)*HX1(0) 
     $  - HX1( -1)*HX2(0,-1) 
     $  - HX1( -1)*HX2(0,1) 
     $  - 2.4674011002723396e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3( -1,-1,1) 
     $  + HX3(0, -1,-1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1) 
      HY3(-1,1,1) = 
     $  + 2.7620719062289241e+00_dp 
     $  - 4.9348022005446793e+00_dp*HX1(-1) 
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(0)*HX1(0) 
     $  + 4.9348022005446793e+00_dp*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(-1,1) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX1(0) *HX2(0,1) 
     $  + HX3( -1,1,1) 
     $  - HX3(0, -1,1) 
     $  + HX3(0,0, -1) 
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      Hi3(0,-1,1) = 
     $  + 8.2246703342411321e-01_dp 
     $  + 6.9314718055994530e-01_dp*HX1(0) 
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  - HX2(0, -1) 
      Hi3(0,1,-1) = 
     $  - 6.9314718055994530e-01_dp*HX1(0) 
      Hi3(-1,-1,1) = 
     $  + 2.4022650695910071e-01_dp 
     $  - 6.9314718055994530e-01_dp*HX1(-1) 
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1) 
     $  - HX1( -1)*HX1(0) 
     $  + 6.9314718055994530e-01_dp*HX1(0) 
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
      Hi3(-1,1,1) = 
     $  + 1.8851605738073271e+00_dp 
     $  + HX1( -1)*HX1(0) 
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0) 
     $  + HX2( -1,1) 
     $  - HX2(0, -1) 
     $  - HX2(0,1) 
      endif 
      if ( nw.gt.3 ) then 
      HY4(0,0,-1,1) = 
     $  + 3.9234217222028759e+00_dp
     $  + 2.5190015545588625e+00_dp*HX1(0)
     $  + 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX4(0,0, -1,1)
     $  - 3.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,0,1,-1) = 
     $  - 4.1940025306306604e+00_dp
     $  - 4.3220869092982539e+00_dp*HX1(0)
     $  - 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX4(0,0,0, -1)
     $  + 3.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,-1,0,1) = 
     $  + 9.4703282949724591e-01_dp
     $  + 1.8030853547393914e+00_dp*HX1(0)
     $  + 1.6449340668482264e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 2.0000000000000000e+00_dp*HX1(0)*HX3(0,0,-1)
     $  - 3.2898681336964528e+00_dp*HX2(0,-1)
     $  + HX4(0, -1,0,1)
     $  + 3.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,-1,1) = 
     $  + 2.5209599327464717e+00_dp
     $  + 2.7620719062289241e+00_dp*HX1(0)
     $  + 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,-1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.4674011002723396e+00_dp*HX2(0,-1)
     $  + 5.0000000000000000e-01_dp*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,-1,1)
     $  + HX4(0, -1,0,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0,1) 
      HY4(0,-1,1,-1) = 
     $  - 8.5266539820739622e+00_dp
     $  - 5.5241438124578482e+00_dp*HX1(0)
     $  - 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1)
     $  + HX1(0) *HX3(0,-1,1)
     $  - 2.0000000000000000e+00_dp*HX1(0)*HX3(0,0,-1)
     $  - HX1(0) *HX3(0,0,1)
     $  + 2.4674011002723396e+00_dp*HX2(0,-1)
     $  - 5.0000000000000000e-01_dp*HX2(0,-1)*HX2(0,-1)
     $  - HX4(0, -1,0,1)
     $  - HX4(0, -1,1,-1)
     $  + 2.0000000000000000e+00_dp*HX4(0,0,-1,-1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,-1,1)
     $  + 4.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  + 3.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  + HX4(0,0,1, -1)
      HY4(0,1,-1,-1) = 
     $  + 5.8027584430066521e+00_dp
     $  + 2.7620719062289241e+00_dp*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + 2.0000000000000000e+00_dp*HX1(0)*HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  - HX4(0,0, -1,-1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  - 3.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,1,-1)
     $  - HX4(0,1, -1,-1)
      HY4(0,-1,1,1) = 
     $  + 6.2689427375197987e-01_dp
     $  - 2.7620719062289241e+00_dp*HX1(0)
     $  - 2.4674011002723396e+00_dp*HX1(0)*HX1(0)
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + 2.0000000000000000e+00_dp*HX1(0)*HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + 4.9348022005446793e+00_dp*HX2(0,-1)
     $  - HX4(0, -1,1,1)
     $  + 2.0000000000000000e+00_dp*HX4(0,0,-1,1)
     $  - 3.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  + HX4(0,0,1,1) 
      HY4(0,1,-1,1) = 
     $  - 4.3326514514433017e+00_dp
     $  - 1.3169446513992682e+00_dp*HX1(0)
     $  - 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,-1)
     $  - 2.0000000000000000e+00_dp*HX1(0)*HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,-1)
     $  + HX2(0, -1)*HX2(0,1)
     $  - 2.4674011002723396e+00_dp*HX2(0,1)
     $  + 5.0000000000000000e-01_dp*HX2(0,1)*HX2(0,1)
     $  - HX4(0, -1,0,1)
     $  - 3.0000000000000000e+00_dp*HX4(0,0,-1,1)
     $  + 3.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  + 4.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,1,1)
     $  - HX4(0,1, -1,1)
      HY4(0,1,1,-1) = 
     $  - 1.5001934240460787e-01_dp
     $  + 4.0790165576281924e+00_dp*HX1(0)
     $  + 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,1)
     $  - HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,1)
     $  - HX2(0, -1)*HX2(0,1)
     $  + 2.4674011002723396e+00_dp*HX2(0,1)
     $  - 5.0000000000000000e-01_dp*HX2(0,1)*HX2(0,1)
     $  + HX4(0, -1,0,1)
     $  + 2.0000000000000000e+00_dp*HX4(0,0,-1,1)
     $  - HX4(0,0,0, -1)
     $  + HX4(0,0,1, -1)
     $  - HX4(0,1,1, -1)
      HY4(-1,-1,-1,1) = 
     $  + 2.4278628067547031e+00_dp
     $  - 2.7620719062289241e+00_dp*HX1(-1)
     $  + 1.2337005501361698e+00_dp*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666e-01_dp*HX1(-1)*HX1(-1)*HX1(-1)*HX1(0)
     $  - 2.5000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  - 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX2(0,-1)
     $  - 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX2(0,1)
     $  - 2.4674011002723396e+00_dp*HX1(-1)*HX1(0)
     $  + 1.6666666666666666e-01_dp*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1( -1)*HX3(0,-1,-1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  + HX1( -1)*HX3(0,1,-1)
     $  + 2.7620719062289241e+00_dp*HX1(0)
     $  + 1.2337005501361698e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX4( -1,-1,-1,1)
     $  - HX4(0, -1,-1,-1)
     $  - HX4(0,0, -1,-1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  - HX4(0,0,1, -1)
     $  - HX4(0,1, -1,-1)
      HY4(-1,-1,1,1) = 
     $  + 2.0293560632083841e+00_dp
     $  + 2.7620719062289241e+00_dp*HX1(-1)
     $  - 2.4674011002723396e+00_dp*HX1(-1)*HX1(-1)
     $  + 2.5000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX1(0)*HX1(0)
     $  + 4.9348022005446793e+00_dp*HX1(-1)*HX1(0)
     $  - 1.6666666666666666e-01_dp*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX1(0)*HX2(0,-1)
     $  - HX1( -1)*HX1(0)*HX2(0,1)
     $  - HX1( -1)*HX3(0,-1,1)
     $  + HX1( -1)*HX3(0,0,-1)
     $  + HX1( -1)*HX3(0,0,1)
     $  - HX1( -1)*HX3(0,1,1)
     $  - 2.7620719062289241e+00_dp*HX1(0)
     $  - 2.4674011002723396e+00_dp*HX1(0)*HX1(0)
     $  + 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX3(-1,-1,1)
     $  + HX1(0) *HX3(0,-1,-1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  + HX1(0) *HX3(0,1,-1)
     $  + HX4( -1,-1,1,1)
     $  + HX4(0, -1,-1,1)
     $  + HX4(0, -1,1,-1)
     $  - HX4(0,0, -1,-1)
     $  + HX4(0,0, -1,1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,0,-1)
     $  - 2.0000000000000000e+00_dp*HX4(0,0,0,1)
     $  - HX4(0,0,1, -1)
     $  + HX4(0,0,1,1) 
     $  + HX4(0,1, -1,1)
     $  + HX4(0,1,1, -1)
      HY4(-1,1,1,1) = 
     $  - 6.4865749331714713e+00_dp
     $  - 4.9348022005446793e+00_dp*HX1(-1)*HX1(0)
     $  + 1.6666666666666666e-01_dp*HX1(-1)*HX1(0)*HX1(0)*HX1(0)
     $  + 2.4674011002723396e+00_dp*HX1(0)*HX1(0)
     $  - 4.1666666666666666e-02_dp*HX1(0)*HX1(0)*HX1(0)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(-1,1)
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,-1)
     $  - 5.0000000000000000e-01_dp*HX1(0)*HX1(0)*HX2(0,1)
     $  + HX1(0) *HX3(-1,1,1)
     $  - HX1(0) *HX3(0,-1,1)
     $  + HX1(0) *HX3(0,0,-1)
     $  + HX1(0) *HX3(0,0,1)
     $  - HX1(0) *HX3(0,1,1)
     $  - 4.9348022005446793e+00_dp*HX2(-1,1)
     $  + 4.9348022005446793e+00_dp*HX2(0,-1)
     $  + 4.9348022005446793e+00_dp*HX2(0,1)
     $  + HX4( -1,1,1,1)
     $  - HX4(0, -1,1,1)
     $  + HX4(0,0, -1,1)
     $  - HX4(0,0,0, -1)
     $  - HX4(0,0,0,1) 
     $  + HX4(0,0,1,1) 
     $  - HX4(0,1,1,1) 
      Hi4(0,0,-1,1) = 
     $  - 9.0154267736969571e-01_dp 
     $  - 8.2246703342411321e-01_dp*HX1(0) 
     $  - 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX3(0,0, -1) 
      Hi4(0,0,1,-1) = 
     $  + 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
      Hi4(0,-1,0,1) = 
     $  + 1.8030853547393914e+00_dp 
     $  + 8.2246703342411321e-01_dp*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  - 2.0000000000000000e+00_dp*HX3(0,0,-1) 
      Hi4(0,-1,-1,1) = 
     $  + 4.8170908494321862e-01_dp 
     $  - 2.4022650695910071e-01_dp*HX1(0) 
     $  - 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,-1) 
     $  + 6.9314718055994530e-01_dp*HX2(0,-1) 
     $  - HX3(0, -1,-1) 
     $  - HX3(0,0, -1) 
      Hi4(0,-1,1,-1) = 
     $  + 5.7009070532142637e-01_dp 
     $  + 4.8045301391820142e-01_dp*HX1(0) 
     $  + 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
     $  - 6.9314718055994530e-01_dp*HX2(0,-1) 
      Hi4(0,1,-1,-1) = 
     $  - 2.4022650695910071e-01_dp*HX1(0) 
      Hi4(0,-1,1,1) = 
     $  - 2.7620719062289241e+00_dp 
     $  - 1.8851605738073271e+00_dp*HX1(0) 
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  - HX1(0) *HX2(0,-1) 
     $  - HX3(0, -1,1) 
     $  + 2.0000000000000000e+00_dp*HX3(0,0,-1) 
     $  + HX3(0,0,1) 
      Hi4(0,1,-1,1) = 
     $  + 2.6736902858507163e+00_dp 
     $  + 1.3029200473423146e+00_dp*HX1(0) 
     $  + 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
     $  + 1.6666666666666665e-01_dp*HX1(0)*HX1(0)*HX1(0) 
     $  + HX1(0) *HX2(0,1) 
     $  + 6.9314718055994530e-01_dp*HX2(0,1) 
     $  - HX3(0,0, -1) 
     $  - 2.0000000000000000e+00_dp*HX3(0,0,1) 
     $  - HX3(0,1, -1) 
      Hi4(0,1,1,-1) =  
     $  + 1.1401814106428527e+00_dp 
     $  + 5.8224052646501250e-01_dp*HX1(0) 
     $  - 3.4657359027997265e-01_dp*HX1(0)*HX1(0) 
     $  - 6.9314718055994530e-01_dp*HX2(0,1) 
      Hi4(-1,-1,-1,1) = 
     $  - 5.5504108664821579e-02_dp
     $  + 2.4022650695910071e-01_dp*HX1(-1)
     $  - 3.4657359027997265e-01_dp*HX1(-1)*HX1(-1)
     $  + 1.6666666666666666e-01_dp*HX1(-1)*HX1(-1)*HX1(-1)
     $  - 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX1(0)
     $  + 6.9314718055994530e-01_dp*HX1(-1)*HX1(0)
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(0)*HX1(0)
     $  - 2.4022650695910071e-01_dp*HX1(0)
     $  - 3.4657359027997265e-01_dp*HX1(0)*HX1(0)
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0)
      Hi4(-1,-1,1,1) = 
     $  - 2.4532465311320902e+00_dp
     $  + 1.8851605738073271e+00_dp*HX1(-1)
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(-1)*HX1(0)
     $  - 5.0000000000000000e-01_dp*HX1(-1)*HX1(0)*HX1(0)
     $  - HX1( -1)*HX2(0,-1)
     $  - HX1( -1)*HX2(0,1)
     $  - 1.8851605738073271e+00_dp*HX1(0)
     $  + 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0)
     $  + HX3( -1,-1,1)
     $  + HX3(0, -1,-1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  + HX3(0,1, -1)
      Hi4(-1,1,1,1) = 
     $  - 5.5504108664821579e-02_dp
     $  - 1.6449340668482264e+00_dp*HX1(-1)
     $  + 5.0000000000000000e-01_dp*HX1(-1)*HX1(0)*HX1(0)
     $  + 1.6449340668482264e+00_dp*HX1(0)
     $  - 1.6666666666666666e-01_dp*HX1(0)*HX1(0)*HX1(0)
     $  + HX1(0) *HX2(-1,1)
     $  - HX1(0) *HX2(0,-1)
     $  - HX1(0) *HX2(0,1)
     $  + HX3( -1,1,1)
     $  - HX3(0, -1,1)
     $  + HX3(0,0, -1)
     $  + HX3(0,0,1) 
     $  - HX3(0,1,1) 
      endif 
** nw > 3 endif 
      endif 
** (n1,n2) = (-1,1) -- completion endif 
      return 
      end 
************************************************************************ 
      subroutine fillirr1dhplin1(y,nw,HY1,HY2,HY3,HY4,n1,n2) 
** evaluates the irreducible HPL for y =1
** it is guaranteed that nw is in the range 2:4, and that (n1,n2) 
** take one of the pairs of values (0,1), (-1,0) or (-1,1) 
      implicit none 
      include 'types.f'
      integer::n1,n2,nw
      real(dp):: y
      real(dp):: HY1(n1:n2),HY2(n1:n2,n1:n2),HY3(n1:n2,n1:n2,n1:n2), 
     $          HY4(n1:n2,n1:n2,n1:n2,n1:n2) 
** (n1,n2) = (0,1) or (-1,1) 
      if (    ( (n1.eq.0).and.(n2.eq.1) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,1) = 
     $  + 1.6449340668482264e+00_dp
      if (nw.gt.2) then 
      HY3(0,0,1) = 
     $  + 1.2020569031595942e+00_dp
      HY3(0,1,1) = 
     $  + 1.2020569031595942e+00_dp
      endif
      if (nw.gt.3) then 
      HY4(0,0,0,1) = 
     $  + 1.0823232337111381e+00_dp
      HY4(0,0,1,1) = 
     $  + 2.7058080842778454e-01_dp
      HY4(0,1,1,1) = 
     $  + 1.0823232337111381e+00_dp
      endif
      endif
** (n1,n2) = (0,1) or (-1,1) endif 
************ 
** (n1,n2) = (-1,0) or (-1,1) 
      if (    ( (n1.eq.-1).and.(n2.eq.0) ) 
     $    .or.( (n1.eq.-1).and.(n2.eq.1) ) ) then 
      HY2(0,-1) = 
     $  + 8.2246703342411321e-01_dp
      if (nw.gt.2) then
      HY3(0,-1,-1) = 
     $  + 1.5025711289494928e-01_dp
      HY3(0,0,-1) = 
     $  + 9.0154267736969571e-01_dp
      endif
      if (nw.gt.3) then
      HY4(0,-1,-1,-1) = 
     $  + 2.3752366322618485e-02_dp
      HY4(0,0,-1,-1) = 
     $  + 8.7785671568655302e-02_dp
      HY4(0,0,0,-1) = 
     $  + 9.4703282949724591e-01_dp
      endif
      endif 
** (n1,n2) = (-1,0) or (-1,1) endif 
** (n1,n2) = (-1,1) -- completion 
      if ( (n1.eq.-1).and.(n2.eq.1) ) then 
      HY2(-1,1) = 
     $  + 5.8224052646501250e-01_dp
      if (nw.gt.2) then
      HY3(0,-1,1) = 
     $  + 2.4307035167006157e-01_dp
      HY3(0,1,-1) = 
     $  + 5.0821521280468485e-01_dp
      HY3(-1,-1,1) = 
     $  + 9.4753004230127705e-02_dp
      HY3(-1,1,1) = 
     $  + 5.3721319360804020e-01_dp
      endif
      if (nw.gt.3) then 
      HY4(0,0,-1,1) = 
     $  + 1.1787599965050932e-01_dp
      HY4(0,0,1,-1) = 
     $  + 1.7284527823898438e-01_dp
      HY4(0,-1,0,1) = 
     $  + 2.0293560632083841e-01_dp
      HY4(0,-1,-1,1) = 
     $  + 3.4159126166513913e-02_dp
      HY4(0,-1,1,-1) = 
     $  + 5.4653052738263652e-02_dp
      HY4(0,1,-1,-1) = 
     $  + 1.1412342741606084e-01_dp
      HY4(0,-1,1,1) = 
     $  + 9.3097125991768577e-02_dp
      HY4(0,1,-1,1) = 
     $  + 1.9355535381306524e-01_dp
      HY4(0,1,1,-1) = 
     $  + 4.3369237704895519e-01_dp
      HY4(-1,-1,-1,1) = 
     $  + 1.4134237214990008e-02_dp
      HY4(-1,-1,1,1) = 
     $  + 4.0758239159309251e-02_dp
      HY4(-1,1,1,1) = 
     $  + 5.1747906167389938e-01_dp
      endif
      endif
** (n1,n2) = (-1,1) -- completion endif 
      return
      end
