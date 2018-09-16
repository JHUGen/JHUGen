      function BSYALslggpppp(e1,p2,e4,p3,za,zb,zab,zba)
      implicit none
      include 'types.f'
      complex(dp):: BSYALslggpppp

C-----Authors: John Campbell and Keith Ellis, March 2012
C---- arXiv:1101.5947 [hep-ph], Eq. (96)
C---- Modification wrt published version
C---- (These are twiddle functions, c.f.arXiv:1101.5947[hep-ph],Eq.(91))
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'zabprods_decl.f'
      include 'qdef.f'
      include 'massiveintegrals.f'
      real(dp):: s23,mt3
      integer:: e1,p2,p3,e4,j
      j=p2-1
      mt3=mt**3
      s23=s(p2,p3)
      BSYALslggpppp =  + I3m12x3x4(j) * ( 1/(za(p3,p2))/(za(p3,p2))*za(
     &    e1,p2)*za(e4,p3)*zb(p3,p2)*mt3 + 1/(za(p3,p2))/(za(p3,p2))*
     &    za(e1,p3)*za(e4,p2)*zb(p3,p2)*mt3 - 1/(zab(p2,q1,p2))/(za(p3,
     &    p2))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*zb(p3,p2)*zab(p3,q1,p2)*
     &    mt3 - 1.D0/2.D0/(zab(p3,q1,p2))/(za(p3,p2))/(za(p3,p2))*za(e1
     &    ,p3)*za(e4,p3)*zb(p3,p2)*zab(p2,q1,p2)*mt3 + 1/(zab(p3,q1,p2)
     &    )/(za(p3,p2))/(za(p3,p2))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,
     &    p2)*mt3 )
      BSYALslggpppp = BSYALslggpppp + I312x3x4(j) * ( 1.D0/2.D0/(zab(p3
     &    ,q1,p2))/(za(p3,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zb(p3,p2
     &    )*zab(p2,q1,p2)*mt3 + 1/(zab(p3,q1,p2))/(za(p3,p2))/(za(p3,p2
     &    ))*za(e4,e1)*zab(p2,q1,p2)*zab(p3,q1,p2)*mt3 )
      BSYALslggpppp = BSYALslggpppp + I3m13x2x4(j) * (  - 1/(za(p3,p2))
     &    /(za(p3,p2))*za(e1,p2)*za(e4,p3)*zb(p3,p2)*mt3 - 1/(za(p3,p2)
     &    )/(za(p3,p2))*za(e1,p3)*za(e4,p2)*zb(p3,p2)*mt3 - 1.D0/2.D0/(
     &    zab(p3,q1,p2))/(za(p3,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*
     &    zb(p3,p2)*zab(p3,q1,p3)*mt3 + 1/(zab(p3,q1,p2))/(za(p3,p2))/(
     &    za(p3,p2))*za(e4,e1)*zab(p3,q1,p2)*zab(p3,q1,p3)*mt3 - 1/(
     &    zab(p3,q1,p3))/(za(p3,p2))/(za(p3,p2))*za(e1,p2)*za(e4,p2)*
     &    zb(p3,p2)*zab(p3,q1,p2)*mt3 )
      BSYALslggpppp = BSYALslggpppp + I313x2x4(j) * ( 1.D0/2.D0/(zab(p3
     &    ,q1,p2))/(za(p3,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zb(p3,p2
     &    )*zab(p3,q1,p3)*mt3 + 1/(zab(p3,q1,p2))/(za(p3,p2))/(za(p3,p2
     &    ))*za(e4,e1)*zab(p3,q1,p2)*zab(p3,q1,p3)*mt3 )
      BSYALslggpppp = BSYALslggpppp + I41x2x4x3 * (  - 1/(za(p3,p2))/(
     &    za(p3,p2))*za(e4,e1)*zab(p2,q1,p2)*s23*mt3 + 1/(za(p3,p2))/(
     &    za(p3,p2))*za(e4,e1)*zab(p3,q1,p3)**2*mt3 - 1.D0/2.D0/(zab(p3
     &    ,q1,p2))/(za(p3,p2))/(za(p3,p2))*za(e1,p3)*za(e4,p3)*zb(p3,p2
     &    )*zab(p2,q1,p2)*zab(p3,q1,p3)*mt3 )

c---- SB sign change: in comparison to published version, the above
c---- SB expression has an additional overall factor of -1
      BSYALslggpppp=-BSYALslggpppp
      return
      end
