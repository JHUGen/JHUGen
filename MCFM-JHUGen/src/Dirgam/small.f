      function smalla(ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: smalla
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smalla=2._dp*V*(ssq+usq)/tsq
      return
      end

      function smallb(ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: smallb
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smallb=2._dp*V*((ssq+usq)/tsq+(ssq+tsq)/usq-2._dp/xn*ssq/(tt*uu))
      return
      end

      function smallc(ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: smallc
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smallc=2._dp*V/xn*(V/(uu*tt)-2*xn**2/ssq)*(tsq+usq)
      return
      end

      function smalld(ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: smalld
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: ss,tt,uu,ssq,tsq,usq
      ssq=ss**2
      tsq=tt**2
      usq=uu**2
      smalld=16._dp*V*xn**2*(3._dp-tt*uu/ssq-ss*uu/tsq-ss*tt/usq)
      return
      end
