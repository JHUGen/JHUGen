      subroutine Ampsq_AQaq_nonid(p1,p2,p3,p4,ampsq)
      implicit none
      include 'types.f'
c--- this routine is a wrapper to the leading order amplitudes
c--- of D.&S. 0906.0008; it performs the squaring and sum over
c--- helicities and returns the ME in the same fashion as h4qn
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3,j4,h1,h2
      real(dp):: ampsq
      complex(dp):: amp(2,2),
     & A0Hqarbmppm,A0Hqarbmpmp

c--- the routine expects momenta to be labelled:
c---                                 q(p1)+Q(p2) --> q(p3)+Q(p4)
c---   and in this routine we use:   0 --> A(j1)+Q(j2)+a(j3)+q(j4)
c--- 
c--- so we permute appropriately
      j1=p2
      j2=p4
      j3=p1
      j4=p3

c--- coded amplitudes      
      amp(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amp(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
      amp(2,1)=A0Hqarbmppm(j2,j1,j4,j3,za,zb)
      amp(2,2)=A0Hqarbmpmp(j2,j1,j4,j3,za,zb)

c--- square, with appropriate overall factor 
c--- (in D.&S. notation, overall 4*C^2*gs^4 is applied in gg_hgg.f)

      ampsq=zip
      do h1=1,2
      do h2=1,2
      ampsq=ampsq+abs(amp(h1,h2))**2
      enddo
      enddo
      
      ampsq=ampsq*V/4._dp
      
      return
      end
      
      
      subroutine Ampsq_AQaq_ident(p1,p2,p3,p4,
     &                            ampsq,ampsq_a,ampsq_b,ampsq_i)
      implicit none
      include 'types.f'
c--- this routine is a wrapper to the leading order amplitudes
c--- of D.&S. 0906.0008; it performs the squaring and sum over
c--- helicities and returns the ME in the same fashion as h4qn
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4,j1,j2,j3,j4,h1,h2
      real(dp):: ampsq,ampsq_a,ampsq_b,ampsq_i
      complex(dp):: amp_a(2,2),amp_b(2,2),
     & A0Hqarbmppm,A0Hqarbmpmp

c--- the routine expects momenta to be labelled:
c---                                 q(p1)+Q(p2) --> q(p3)+Q(p4)
c---   and in this routine we use:   0 --> A(j1)+Q(j2)+a(j3)+q(j4)
c--- 
c--- so we permute appropriately
      j1=p2
      j2=p4
      j3=p1
      j4=p3

c--- First calculate with the normal ordering, (1,2,3,4)
c--- coded amplitudes      
      amp_a(1,2)=A0Hqarbmppm(j1,j2,j3,j4,za,zb)
      amp_a(1,1)=A0Hqarbmpmp(j1,j2,j3,j4,za,zb)
c--- charge conjugation      
      amp_a(2,1)=A0Hqarbmppm(j2,j1,j4,j3,za,zb)
      amp_a(2,2)=A0Hqarbmpmp(j2,j1,j4,j3,za,zb)

c--- Now calculate with the exhanged order,    (1,4,3,2) 
c--- coded amplitudes      
      amp_b(1,2)=A0Hqarbmppm(j1,j4,j3,j2,za,zb)
      amp_b(1,1)=A0Hqarbmpmp(j1,j4,j3,j2,za,zb)
c--- charge conjugation      
      amp_b(2,1)=A0Hqarbmppm(j4,j1,j2,j3,za,zb)
      amp_b(2,2)=A0Hqarbmpmp(j4,j1,j2,j3,za,zb)

c--- square, with appropriate overall factor 
c--- (in D.&S. notation, overall 4*C^2*gs^4 is applied in gg_hgg.f)

      ampsq_a=zip
      ampsq_b=zip
      ampsq_i=zip
      do h1=1,2
      do h2=1,2
      ampsq_a=ampsq_a+abs(amp_a(h1,h2))**2
      ampsq_b=ampsq_b+abs(amp_b(h1,h2))**2
      if (h1 == h2) then
        ampsq_i=ampsq_i+2._dp/xn*real(amp_a(h1,h2)*conjg(amp_b(h1,h2)),dp)
      endif
      enddo
      enddo
      
      ampsq_a=ampsq_a*V/4._dp
      ampsq_b=ampsq_b*V/4._dp
      ampsq_i=ampsq_i*V/4._dp
      
      ampsq=ampsq_a+ampsq_b+ampsq_i
      
      return
      end
