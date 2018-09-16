C----These are reduced coefficients at times with
C----dimensional factors removed. For exact factors
C----see writecoeffs.f
      complex(dp):: coeff(0:4,20),scints(4,20,-2:0)

      integer:: d2x3x4,d1x2x3,d1x23x4,d1x2x34,d12x3x4,
     & d23x1x4,d23x4x1
      parameter (d2x3x4=1,d1x2x3=2,d1x23x4=3,d1x2x34=4,d12x3x4=5,
     & d23x1x4=6,d23x4x1=7)

      integer:: c23x4,c1x23,c2x3,c2x34,c3x4,c12x3,c1x2,
     & c1x234,c123x4,c12x34,c14x23,c1x4,c2x3x2m,c14x23bis

      parameter (c23x4=1,c1x23=2,c2x3=3,c2x34=4,c3x4=5,c12x3=6,c1x2=7,
     & c1x234=8,c123x4=9,c12x34=10,c14x23=11,c1x4=12,c2x3x2m=13,
     & c14x23bis=14)

      integer:: b123,b234,b23x2m,b1234,b14,b23,b2x1m,b12,b34
      parameter (b123=1,b234=2,b23x2m=3,b1234=4,b14=5,b23=6,b2x1m=7,
     & b12=8,b34=9)

      integer:: a0m
      parameter (a0m=1)

      integer:: irat 
      parameter (irat=1)
