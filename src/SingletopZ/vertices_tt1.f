      subroutine vertices_tt1(mtsq,ep,
     & vert1,vert2,vert3,vert4,vert5)
      include 'constants.f'
      include 'poles.f'
      include 'scale.f'
      include 'masses.f'
      include 'decl_kininv.f'
      
      double precision mtsq
      integer ep
      double complex vert1,vert2,vert3,vert4,vert5 
      double precision dentt1
      double complex qlI2diffs345s16,qlI2,qlI3
       mtsq=mt**2
       dentt1=1d0/( -s345 + s16)
       qlI2diffs345s16=qlI2(s345,0d0,mtsq,musq,ep)-
     & qlI2(s16,0d0,mtsq,musq,ep)
      
      vert1= + qlI2(s16,zip,mtsq,musq,ep) * ( 4.D0*mt*dentt1 )
      vert1 = vert1 + qlI3(zip,s345,s16,zip,zip,mtsq,musq,ep) * ( 4.D0*
     &    mt*dentt1*s16 - 4.D0*mt**3*dentt1 )
      vert1 = vert1 + 8.D0*qlI2diffs345s16*mt*dentt1 - 8.D0*
     &    qlI2diffs345s16*mt*dentt1**2*s16

       vert2= + 4.D0*qlI2diffs345s16*mt*dentt1

       vert3= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 2.D0*mt**2*dentt1*
     &    s16**(-1) )
      vert3 = vert3 + qlI2(s16,zip,mtsq,musq,ep) * (  - 4.D0*dentt1 + 2.
     &    D0*mt**2*dentt1*s16**(-1) )
      vert3 = vert3 + qlI3(zip,s345,s16,zip,zip,mtsq,musq,ep) * (  - 4.D
     &    0*dentt1*s16 + 4.D0*mt**2*dentt1 )
      vert3 = vert3 - 2.D0*fp(ep)*dentt1 + 2.D0*fp(ep)*mt**2*dentt1*
     &    s16**(-1) - 6.D0*qlI2diffs345s16*dentt1 + 10.D0*
     &    qlI2diffs345s16*dentt1**2*s16 - 2.D0*qlI2diffs345s16*mt**2*
     &    dentt1**2

       vert4= + qlI2(mtsq,zip,mtsq,musq,ep) * (  - 2.D0*mt**2*dentt1*
     &    s16**(-1) + 2.D0*mt**2*dentt1*s345**(-1) )
      vert4 = vert4 + qlI2(s16,zip,mtsq,musq,ep) * ( 2.D0*mt**2*dentt1*
     &    s16**(-1) - 2.D0*mt**2*dentt1*s345**(-1) )
      vert4 = vert4 + 2.D0*fp(ep)*mt**2*dentt1*s16**(-1) - 2.D0*fp(ep)*
     &    mt**2*dentt1*s345**(-1) - 2.D0*qlI2diffs345s16*dentt1 - 2.D0*
     &    qlI2diffs345s16*mt**2*dentt1*s345**(-1)

       vert5= + qlI2(s16,zip,mtsq,musq,ep) * ( 1.D0 )
      vert5 = vert5 + qlI3(zip,s345,s16,zip,zip,mtsq,musq,ep) * ( 2.D0*
     &    s16 - 2.D0*mt**2 )
      vert5 = vert5 + fp(ep) + qlI2diffs345s16 - 3.D0*qlI2diffs345s16*
     &    dentt1*s16 + qlI2diffs345s16*mt**2*dentt1

      
      return
      end
      
