      function Ampvirtsq_AQaq_nonid(p1,p2,p3,p4)
      implicit none
      include 'types.f'
      real(dp):: Ampvirtsq_AQaq_nonid
c--- this routine is a wrapper to the 1-loop/Born interference
c--- of D.&S. 0906.0008; it performs the interference and sum over
c--- helicities and returns the ME in the same fashion as Hqarbvsq
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3,j4,h1,h2
      real(dp):: sum,ren,H4prenorm
      complex(dp):: amplo(2,2),ampv(2,2),
     & A0Hqarbmppm,A0Hqarbmpmp,A41Hqarbmppm,A41Hqarbmpmp

c--- the routine expects momenta to be labelled:
c---                                 a(p1)+q(p2) --> Q(p3)+A(p4)
c---   and in this routine we use:   0 --> A(j1)+Q(j2)+a(j3)+q(j4)
c--- 
c--- so we permute appropriately (NOTE: different from LO routine)
      j1=p4
      j2=p3
      j3=p2
      j4=p1

c--- Born terms
c--- coded amplitudes      
      amplo(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amplo(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
c      amplo(2,1)=A0Hqarbmppm(j2,j1,j4,j3,za,zb)
c      amplo(2,2)=A0Hqarbmpmp(j2,j1,j4,j3,za,zb)
c--- equivalent to complex conjugation for Born amplitudes    
      amplo(2,1)=conjg(amplo(1,2))
      amplo(2,2)=conjg(amplo(1,1))

c--- 1-loop terms
c--- coded amplitudes      
      ampv(1,2) =A41Hqarbmppm(j1,j2,j3,j4,za,zb)
      ampv(1,1) =A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
      ampv(2,1) =A41Hqarbmppm(j2,j1,j4,j3,za,zb)
      ampv(2,2) =A41Hqarbmpmp(j2,j1,j4,j3,za,zb)

c--- get renormalization factor
      ren=H4prenorm()

c--- square, with appropriate overall factor (in D.&S. notation, 
c---   overall 4*C^2*gs^4*ason2pi is applied in gg_hgg_v.f)
      sum=zip
      do h1=1,2
      do h2=1,2
      ampv(h1,h2)=ampv(h1,h2)+ren*amplo(h1,h2)

      sum=sum+real(ampv(h1,h2)*conjg(amplo(h1,h2)))
      enddo
      enddo
      
      Ampvirtsq_AQaq_nonid=sum*V*xn/4._dp
                  
      return
      end


      function Ampvirtsq_AQaq_ident(p1,p2,p3,p4)
      implicit none
      include 'types.f'
      real(dp):: Ampvirtsq_AQaq_ident
c--- this routine is a wrapper to the 1-loop/Born interference
c--- of D.&S. 0906.0008; it performs the interference and sum over
c--- helicities and returns the ME in the same fashion as Hqarbvsq
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3,j4,h1,h2
      real(dp):: sum,ren,H4prenorm
      complex(dp):: amplo_a(2,2),amplo_b(2,2),
     & amp42_a(2,2),amp42_b(2,2),
     & A0Hqarbmpmp,A42Hqarbmpmp

c--- the routine expects momenta to be labelled:
c---                                 a(p1)+q(p2) --> Q(p3)+A(p4)
c---   and in this routine we use:   0 --> A(j1)+Q(j2)+a(j3)+q(j4)
c--- 
c--- so we permute appropriately (NOTE: different from LO routine)
      j1=p4
      j2=p3
      j3=p2
      j4=p1

c--- Note that the A42 amplitude only contributes for identical
c--- quarks and, therefore, we only need the case of equal helicities

c--- First calculate with the normal ordering, (1,2,3,4)
c--- Born terms
c--- coded amplitudes      
c      amplo_a(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amplo_a(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
c      amplo_a(2,1)=A0Hqarbmppm(j2,j1,j4,j3,za,zb)
c      amplo_a(2,2)=A0Hqarbmpmp(j2,j1,j4,j3,za,zb)
c--- equivalent to complex conjugation for Born amplitudes    
      amplo_a(2,2)=conjg(amplo_a(1,1))

c--- 1-loop terms
c--- coded amplitudes      
c      amp41_a(1,2) =A41Hqarbmppm(j1,j2,j3,j4,za,zb)
c      amp41_a(1,1) =A41Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
c      amp41_a(2,1) =A41Hqarbmppm(j2,j1,j4,j3,za,zb)
c      amp41_a(2,2) =A41Hqarbmpmp(j2,j1,j4,j3,za,zb)

c--- coded amplitudes      
c      amp42_a(1,2) =A42Hqarbmppm(j1,j2,j3,j4,za,zb)
      amp42_a(1,1) =A42Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
c      amp42_a(2,1) =A42Hqarbmppm(j2,j1,j4,j3,za,zb)
      amp42_a(2,2) =A42Hqarbmpmp(j2,j1,j4,j3,za,zb)

c--- Now calculate with the exhanged order,    (1,4,3,2) 
c--- Born terms
c--- coded amplitudes      
c      amplo_b(1,2)=A0Hqarbmppm(j1,j4,j3,j2,za,zb)
      amplo_b(1,1)=A0Hqarbmpmp(j1,j4,j3,j2,za,zb)
c--- charge conjugation      
c      amplo_b(2,1)=A0Hqarbmppm(j4,j1,j2,j3,za,zb)
c      amplo_b(2,2)=A0Hqarbmpmp(j4,j1,j2,j3,za,zb)
c--- equivalent to complex conjugation for Born amplitudes    
      amplo_b(2,2)=conjg(amplo_b(1,1))

c--- 1-loop terms
c--- coded amplitudes      
c      amp41_b(1,2) =A41Hqarbmppm(j1,j4,j3,j2,za,zb)
c      amp41_b(1,1) =A41Hqarbmpmp(j1,j4,j3,j2,za,zb)
c--- charge conjugation      
c      amp41_b(2,1) =A41Hqarbmppm(j4,j1,j2,j3,za,zb)
c      amp41_b(2,2) =A41Hqarbmpmp(j4,j1,j2,j3,za,zb)

c--- coded amplitudes      
c      amp42_b(1,2) =A42Hqarbmppm(j1,j4,j3,j2,za,zb)
      amp42_b(1,1) =A42Hqarbmpmp(j1,j4,j3,j2,za,zb)
c--- charge conjugation      
c      amp42_b(2,1) =A42Hqarbmppm(j4,j1,j2,j3,za,zb)
      amp42_b(2,2) =A42Hqarbmpmp(j4,j1,j2,j3,za,zb)


c--- get renormalization factor
      ren=H4prenorm()

c--- square, with appropriate overall factor (in D.&S. notation, 
c---   overall 4*C^2*gs^4*ason2pi is applied in gg_hgg_v.f)
      sum=zip
      do h1=1,2
c      do h2=1,2
      h2=h1
      
c--- note: A42 is renormalized with the opposite sign
c      amp41_a(h1,h2)=amp41_a(h1,h2)+ren*amplo_a(h1,h2)
c      amp41_b(h1,h2)=amp41_b(h1,h2)+ren*amplo_b(h1,h2)
      amp42_a(h1,h2)=amp42_a(h1,h2)-ren*amplo_a(h1,h2)
      amp42_b(h1,h2)=amp42_b(h1,h2)-ren*amplo_b(h1,h2)

c--- We commment out these two contributions in order that we can
c--- take over the existing structure in gg_hgg_v.f. This routine
c--- then returns the same quantity as Hqaqavsq.f
c      sum=sum+real(amp41_a(h1,h2)*conjg(amplo_a(h1,h2)))
c     &       +real(amp41_b(h1,h2)*conjg(amplo_b(h1,h2)))
c      if (h1 == h2) then
      sum=sum-1._dp/xn*real(amp42_a(h1,h2)*conjg(amplo_b(h1,h2)))
     &       -1._dp/xn*real(amp42_b(h1,h2)*conjg(amplo_a(h1,h2)))      
c      endif
c      enddo
      enddo
      
      Ampvirtsq_AQaq_ident=sum*V*xn/4._dp
                  
      return
      end
