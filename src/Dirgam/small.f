      double precision function smalla(ss,tt,uu)
      implicit none
      include 'constants.f'
      double precision ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smalla=2d0*V*(ssq+usq)/tsq
      return
      end

      double precision function smallb(ss,tt,uu)
      implicit none
      include 'constants.f'
      double precision ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smallb=2d0*V*((ssq+usq)/tsq+(ssq+tsq)/usq-2d0/xn*ssq/(tt*uu))
      return
      end

      double precision function smallc(ss,tt,uu)
      implicit none
      include 'constants.f'
      double precision ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smallc=2d0*V/xn*(V/(uu*tt)-2*xn**2/ssq)*(tsq+usq)
      return
      end

      double precision function smalld(ss,tt,uu)
      implicit none
      include 'constants.f'
      double precision ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smalld=16d0*V*xn**2*(3d0-tt*uu/ssq-ss*uu/tsq-ss*tt/usq)
      return
      end
