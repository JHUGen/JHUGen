      subroutine Ampvirt_AQgg(p1,p2,p3,p4,ab41,ba41,ab43,ba43)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      integer:: p1,p2,p3,p4
      complex(dp)::ab41(2,2,2),ba41(2,2,2),ab43(2,2,2),ba43(2,2,2),
     & A41HAQggmppp,A41HAQggmpmm,A41HAQggmpmp,A41HAQggmppm,
     & A43HAQggmppp,A43HAQggmpmm,A43HAQggmpmp,A43HAQggmppm

c--- calculate all the A41 amplitudes
      ab41(1,2,2)=A41HAQggmppp(p1,p2,p3,p4,za,zb)
      ab41(1,1,1)=A41HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab41(1,1,2)=A41HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab41(1,2,1)=A41HAQggmppm(p1,p2,p3,p4,za,zb)

      ba41(1,2,2)=A41HAQggmppp(p1,p2,p4,p3,za,zb)
      ba41(1,1,1)=A41HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba41(1,1,2)=A41HAQggmppm(p1,p2,p4,p3,za,zb)
      ba41(1,2,1)=A41HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- Obtain opposite helicities through parity relation, c.f. Eq.(4.6) of DS 0906.0008
      ab41(2,1,1)=-A41HAQggmppp(p1,p2,p3,p4,zb,za)
      ab41(2,2,2)=-A41HAQggmpmm(p1,p2,p3,p4,zb,za)
      ab41(2,2,1)=-A41HAQggmpmp(p1,p2,p3,p4,zb,za)
      ab41(2,1,2)=-A41HAQggmppm(p1,p2,p3,p4,zb,za)

      ba41(2,1,1)=-A41HAQggmppp(p1,p2,p4,p3,zb,za)
      ba41(2,2,2)=-A41HAQggmpmm(p1,p2,p4,p3,zb,za)
      ba41(2,2,1)=-A41HAQggmppm(p1,p2,p4,p3,zb,za)
      ba41(2,1,2)=-A41HAQggmpmp(p1,p2,p4,p3,zb,za)


c--- calculate all the A43 amplitudes
      ab43(1,2,2)=A43HAQggmppp(p1,p2,p3,p4,za,zb)
      ab43(1,1,1)=A43HAQggmpmm(p1,p2,p3,p4,za,zb)
      ab43(1,1,2)=A43HAQggmpmp(p1,p2,p3,p4,za,zb)
      ab43(1,2,1)=A43HAQggmppm(p1,p2,p3,p4,za,zb)

      ba43(1,2,2)=A43HAQggmppp(p1,p2,p4,p3,za,zb)
      ba43(1,1,1)=A43HAQggmpmm(p1,p2,p4,p3,za,zb)
      ba43(1,1,2)=A43HAQggmppm(p1,p2,p4,p3,za,zb)
      ba43(1,2,1)=A43HAQggmpmp(p1,p2,p4,p3,za,zb)

c--- Obtain opposite helicities through parity relation, c.f. Eq.(4.6) of DS 0906.0008
      ab43(2,1,1)=-A43HAQggmppp(p1,p2,p3,p4,zb,za)
      ab43(2,2,2)=-A43HAQggmpmm(p1,p2,p3,p4,zb,za)
      ab43(2,2,1)=-A43HAQggmpmp(p1,p2,p3,p4,zb,za)
      ab43(2,1,2)=-A43HAQggmppm(p1,p2,p3,p4,zb,za)

      ba43(2,1,1)=-A43HAQggmppp(p1,p2,p4,p3,zb,za)
      ba43(2,2,2)=-A43HAQggmpmm(p1,p2,p4,p3,zb,za)
      ba43(2,2,1)=-A43HAQggmppm(p1,p2,p4,p3,zb,za)
      ba43(2,1,2)=-A43HAQggmpmp(p1,p2,p4,p3,zb,za)


      return
      end
