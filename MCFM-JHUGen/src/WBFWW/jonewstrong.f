      subroutine jonewstrong(p7,p3,p4,p1,za,zb,zab,jqcd)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'cmplxmass.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: p1,p2,p3,p4,p7,ro
      complex(dp):: zab(mxpart,4,mxpart),jqcd(2,4),propw34
      real(dp):: t3,s34,s134,s347
C-----Begin statement functions
      t3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
C-----end statement functions

      s34=s(p3,p4)
      s134=t3(p1,p3,p4)
      s347=t3(p3,p4,p7)

      propw34=s34-cwmass2

      do ro=1,4
      jqcd(1,ro)= + propw34**(-1) * ( za(p7,p3)*zb(p7,p4)*zab(p7,ro,p1)
     &    *s347**(-1) + za(p7,p3)*zb(p3,p4)*zab(p3,ro,p1)*s347**(-1) - 
     &    za(p1,p3)*zb(p1,p4)*zab(p7,ro,p1)*s134**(-1) + za(p3,p4)*zb(
     &    p1,p4)*zab(p7,ro,p4)*s134**(-1) )
      jqcd(2,ro)=jqcd(1,ro)
      enddo
      jqcd(:,:)=jqcd(:,:)/cxw
      return
      end
