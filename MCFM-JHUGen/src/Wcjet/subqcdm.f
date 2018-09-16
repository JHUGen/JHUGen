      subroutine subqcdm(i1,i2,i3,i4,i5,i6,p156,p256,za,zb,
     & invtwog1Dc,invtwog2Dc,mc,aamp,bamp)
      implicit none
      include 'types.f'
c*******************************************************************
c     the matrix elements of the
C     helicity amplitudes for the QCD process
c     s(-p1)+cbar(-p2) --> l(p3)+abar(p4)+g(p5)+g(p6)
c     multiplied by ((a+l)^2-M**2)/g^4/gwsq^2/2
c*******************************************************************

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_decl.f'
      integer:: i1,i2,i3,i4,i5,i6,h1,h2,hf
      real(dp):: p156,p256,s34,mc,invtwog1Dc,invtwog2Dc
C     p156=2*p1.p5+2*p1.p6+2*p5.p6
C     p256=2*p2.p5+2*p2.p6+2*p5.p6
      complex(dp):: aamp(2,2,2),bamp(2,2,2)
      s34=real(za(i3,i4)*zb(i4,i3))
      aamp(1,1,1) =  + four*p256**(-1) * (  -1./(zb(i5,i1))/(zb(i5,i6))
     &    *za(i4,i3)*za(i6,i2)*zb(i1,i4)**2 -1./(zb(i5,i6))/(zb(i6,i1))
     &    *za(i4,i3)*za(i5,i2)*zb(i1,i4)**2 )
      aamp(1,1,1) = aamp(1,1,1) + four*p256**(-1)*invtwog2Dc * (1./(zb(
     &    i5,i1))*za(i4,i3)*za(i5,i6)*za(i6,i2)*zb(i1,i4)**2 -1./(zb(i5
     &    ,i1))/(zb(i6,i1))*za(i4,i3)*za(i5,i2)*za(i6,i2)*zb(i1,i2)*zb(
     &    i1,i4)**2 )

      aamp(2,2,2) =  + mc*four*p256**(-1) * (1./(za(i5,i1))/(za(i5,i6))
     &    *za(i1,i3)*zb(i1,i4)*zb(i6,i2) +1./(za(i5,i6))/(za(i6,i1))*
     &    za(i1,i3)*zb(i1,i4)*zb(i5,i2) )
      aamp(2,2,2) = aamp(2,2,2) + mc*four*p256**(-1)*invtwog2Dc * (  -
     &   1./(za(i5,i1))*za(i1,i3)*zb(i1,i4)*zb(i5,i6)*zb(i6,i2) +1./(
     &    za(i5,i1))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i5,i2
     &    )*zb(i6,i2) )
      aamp(2,2,2) = aamp(2,2,2) + mc*four*p156**(-1) * (1./(za(i1,i2))
     &    /(za(i5,i1))*za(i1,i3)*zb(i5,i6)*zb(i6,i4) +1./(za(i1,i2))/(
     &    za(i5,i1))/(za(i5,i6))*za(i1,i3)*za(i6,i1)*zb(i6,i1)*zb(i6,i4
     &    ) +1./(za(i1,i2))/(za(i5,i6))*za(i1,i3)*zb(i5,i1)*zb(i6,i4)
     &     +1./(za(i1,i2))/(za(i5,i6))*za(i1,i3)*zb(i5,i4)*zb(i6,i1) +
     &   1./(za(i1,i2))/(za(i5,i6))/(za(i6,i1))*za(i1,i3)*za(i5,i1)*zb(
     &    i5,i1)*zb(i5,i4) +1./(za(i1,i2))/(za(i6,i1))*za(i1,i3)*zb(i5,
     &    i4)*zb(i5,i6) )
      aamp(2,2,2) = aamp(2,2,2) + mc*four*invtwog2Dc * (  -1./(za(i5,i1
     &    ))/(za(i6,i1))*za(i1,i3)*zb(i5,i4)*zb(i6,i2) )

      aamp(2,1,1) =  + four*p256**(-1)*invtwog2Dc * (1./(za(i5,i6))/(
     &    zb(i5,i6))*za(i1,i3)*za(i6,i2)**2*zb(i1,i4)*zb(i5,i1)*zb(i5,
     &    i2) +1./(za(i5,i6))/(zb(i5,i6))*za(i4,i3)*za(i6,i2)**2*zb(i1,
     &    i4)*zb(i5,i2)*zb(i5,i4) )
      aamp(2,1,1) = aamp(2,1,1) + four*p156**(-1) * (  -1./(za(i5,i1))
     &    /(za(i5,i6))/(zb(i5,i6))*za(i3,i2)*za(i6,i1)**2*zb(i1,i4)*zb(
     &    i5,i1) +1./(za(i5,i1))/(zb(i5,i6))*za(i3,i2)*za(i6,i1)*zb(i5,
     &    i1)*zb(i5,i4) )
      aamp(2,1,1) = aamp(2,1,1) + four*invtwog2Dc * (1./(za(i5,i1))*za(
     &    i6,i2)*za(i6,i3)*zb(i5,i4) -1./(za(i5,i1))/(za(i5,i6))*za(i6,
     &    i1)*za(i6,i2)*za(i6,i3)*zb(i1,i4) +1./(za(i5,i1))/(za(i5,i6))
     &    /(zb(i5,i6))*za(i3,i2)*za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i2
     &    ) -1./(za(i5,i1))/(zb(i5,i6))*za(i3,i2)*za(i6,i2)*zb(i5,i2)*
     &    zb(i5,i4) )
      aamp(2,1,1) = aamp(2,1,1) + mc**2*four*p256**(-1)*invtwog2Dc * (
     &     -1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3
     &    )*za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i1)**2 -1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,i3)*za(i6,i1)*za(
     &    i6,i2)*zb(i1,i4)*zb(i5,i1)*zb(i5,i4) -1./(za(i5,i6))/(zb(i1,
     &    i2))/(zb(i5,i6))*za(i6,i2)*za(i6,i3)*zb(i1,i4)*zb(i5,i1)*zb(
     &    i5,i2) )
      aamp(2,1,1) = aamp(2,1,1) + mc**2*four*invtwog2Dc * (1./(za(i1,i2
     &    ))/(za(i5,i1))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*
     &    za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i1) -1./(za(i1,i2))/(za(
     &    i5,i1))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*za(i6,i2)*zb(i5,i1)
     &    *zb(i5,i4) -1./(za(i5,i1))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6)
     &    )*za(i6,i1)*za(i6,i3)*zb(i1,i4)*zb(i5,i1) +1./(za(i5,i1))/(
     &    zb(i1,i2))/(zb(i5,i6))*za(i6,i3)*zb(i5,i1)*zb(i5,i4) )
      aamp(2,1,1) = aamp(2,1,1) + mc**4*four*p256**(-1)*invtwog2Dc * (
     &   1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i1,i2))/(zb(i5,i6))
     &    *za(i6,i1)*za(i6,i3)*zb(i1,i4)*zb(i5,i1)**2 )

      aamp(1,2,1) =  + four*p256**(-1)*invtwog2Dc * (1./(za(i5,i6))/(
     &    zb(i5,i6))*za(i1,i3)*za(i5,i2)**2*zb(i1,i4)*zb(i6,i1)*zb(i6,
     &    i2) +1./(za(i5,i6))/(zb(i5,i6))*za(i4,i3)*za(i5,i2)**2*zb(i1,
     &    i4)*zb(i6,i2)*zb(i6,i4) )
      aamp(1,2,1) = aamp(1,2,1) + four*p156**(-1) * (  -1./(za(i5,i6))
     &    /(zb(i5,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i1)*zb(i1,i4)*zb(i6,
     &    i1)**2 -1./(zb(i5,i1))/(zb(i5,i6))*za(i3,i2)*zb(i6,i1)**2*zb(
     &    i6,i4) )
      aamp(1,2,1) = aamp(1,2,1) + four*invtwog2Dc * (1./(za(i5,i6))/(
     &    zb(i5,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i2)*zb(i1,i4)*zb(i6,i1
     &    )*zb(i6,i2) )
      aamp(1,2,1) = aamp(1,2,1) + mc**2*four*p256**(-1)*invtwog2Dc * (
     &     -1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3
     &    )*za(i5,i1)*za(i5,i2)*zb(i1,i4)*zb(i6,i1)**2 -1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,i3)*za(i5,i1)*za(
     &    i5,i2)*zb(i1,i4)*zb(i6,i1)*zb(i6,i4) -1./(za(i5,i6))/(zb(i1,
     &    i2))/(zb(i5,i6))*za(i5,i2)*za(i5,i3)*zb(i1,i4)*zb(i6,i1)*zb(
     &    i6,i2) )
      aamp(1,2,1) = aamp(1,2,1) + mc**2*four*invtwog2Dc * (1./(za(i1,i2
     &    ))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i1))/(zb(i5,i6))*za(i1,i3)*
     &    za(i5,i2)*zb(i1,i4)*zb(i6,i1)**2 -1./(za(i5,i6))/(zb(i1,i2))
     &    /(zb(i5,i1))/(zb(i5,i6))*za(i5,i3)*zb(i1,i4)*zb(i6,i1)**2 )
      aamp(1,2,1) = aamp(1,2,1) + mc**4*four*p256**(-1)*invtwog2Dc * (
     &   1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i1,i2))/(zb(i5,i6))
     &    *za(i5,i1)*za(i5,i3)*zb(i1,i4)*zb(i6,i1)**2 )

      aamp(1,1,2) =  + mc*four*p256**(-1) * (  -1./(za(i1,i2))/(zb(i5,
     &    i1))/(zb(i5,i6))*za(i4,i3)*za(i6,i1)*zb(i1,i4)**2 -1./(za(i1,
     &    i2))/(zb(i5,i6))/(zb(i6,i1))*za(i4,i3)*za(i5,i1)*zb(i1,i4)**2
     &     +1./(zb(i5,i1))/(zb(i5,i6))*za(i6,i3)*zb(i1,i2)*zb(i1,i4) +
     &   1./(zb(i5,i6))/(zb(i6,i1))*za(i5,i3)*zb(i1,i2)*zb(i1,i4) )
      aamp(1,1,2) = aamp(1,1,2) + mc*four*p256**(-1)*invtwog2Dc * (1./(
     &    za(i1,i2))/(zb(i5,i1))*za(i4,i3)*za(i5,i6)*za(i6,i1)*zb(i1,i4
     &    )**2 -1./(za(i1,i2))/(zb(i5,i1))/(zb(i6,i1))*za(i4,i3)*za(i5,
     &    i2)*za(i6,i1)*zb(i1,i2)*zb(i1,i4)**2 -1./(zb(i5,i1))/(zb(i6,
     &    i1))*za(i4,i3)*za(i5,i6)*zb(i1,i2)*zb(i1,i4)**2 +1./(zb(i5,i1
     &    ))/(zb(i6,i1))*za(i5,i3)*za(i6,i2)*zb(i1,i2)**2*zb(i1,i4) )

      aamp(1,2,2) =  + mc*four*p256**(-1)*invtwog2Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))*za(i1,i3)*za(i5,i1)*za(i5,i2)*zb(i1,
     &    i4)*zb(i6,i1)*zb(i6,i2) +1./(za(i1,i2))/(za(i5,i6))/(zb(i5,i6
     &    ))*za(i4,i3)*za(i5,i1)*za(i5,i2)*zb(i1,i4)*zb(i6,i2)*zb(i6,i4
     &    ) -1./(za(i5,i6))/(zb(i5,i6))*za(i5,i2)*za(i5,i3)*zb(i1,i4)*
     &    zb(i6,i2)**2 )
      aamp(1,2,2) = aamp(1,2,2) + mc*four*p156**(-1) * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i1))/(zb(i5,i6))*za(i1,i3)*za(i5,i1)*zb(
     &    i1,i4)*zb(i6,i1)**2 +1./(za(i1,i2))/(zb(i5,i1))/(zb(i5,i6))*
     &    za(i1,i3)*zb(i6,i1)**2*zb(i6,i4) )
      aamp(1,2,2) = aamp(1,2,2) + mc*four*invtwog2Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i1)*zb(
     &    i1,i4)*zb(i6,i1)*zb(i6,i2) -1./(za(i5,i6))/(zb(i5,i1))/(zb(i5
     &    ,i6))*za(i5,i3)*zb(i1,i4)*zb(i6,i1)*zb(i6,i2) )
      aamp(1,2,2) = aamp(1,2,2) + mc**3*four*p256**(-1)*invtwog2Dc * (
     &     -1./(za(i1,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,
     &    i6))*za(i1,i3)*za(i5,i1)**2*zb(i1,i4)*zb(i6,i1)**2 -1./(za(i1
     &    ,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,
     &    i3)*za(i5,i1)**2*zb(i1,i4)*zb(i6,i1)*zb(i6,i4) +1./(za(i1,i2)
     &    )/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i5,i1)*za(i5,i3)*zb(
     &    i1,i4)*zb(i6,i1)*zb(i6,i2) )
      aamp(1,2,2) = aamp(1,2,2) + mc**3*four*invtwog2Dc * (1./(za(i1,i2
     &    ))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i1))/(zb(i5,i6)
     &    )*za(i1,i3)*za(i5,i1)*zb(i1,i4)*zb(i6,i1)**2 )

      aamp(2,1,2) =  + mc*four*p256**(-1)*invtwog2Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))*za(i1,i3)*za(i6,i1)*za(i6,i2)*zb(i1,
     &    i4)*zb(i5,i1)*zb(i5,i2) +1./(za(i1,i2))/(za(i5,i6))/(zb(i5,i6
     &    ))*za(i4,i3)*za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i2)*zb(i5,i4
     &    ) -1./(za(i5,i6))/(zb(i5,i6))*za(i6,i2)*za(i6,i3)*zb(i1,i4)*
     &    zb(i5,i2)**2 )
      aamp(2,1,2) = aamp(2,1,2) + mc*four*p156**(-1) * (1./(za(i1,i2))
     &    /(za(i5,i1))/(za(i5,i6))/(zb(i5,i6))*za(i1,i3)*za(i6,i1)**2*
     &    zb(i1,i4)*zb(i5,i1) -1./(za(i1,i2))/(za(i5,i1))/(zb(i5,i6))*
     &    za(i1,i3)*za(i6,i1)*zb(i5,i1)*zb(i5,i4) )
      aamp(2,1,2) = aamp(2,1,2) + mc*four*invtwog2Dc * (1./(za(i1,i2))
     &    /(za(i5,i1))*za(i6,i1)*za(i6,i3)*zb(i5,i4) -1./(za(i1,i2))/(
     &    za(i5,i1))/(za(i5,i6))*za(i6,i1)**2*za(i6,i3)*zb(i1,i4) +1./(
     &    za(i1,i2))/(za(i5,i1))/(za(i5,i6))/(zb(i5,i6))*za(i3,i2)*za(
     &    i6,i1)**2*zb(i1,i4)*zb(i5,i2) -1./(za(i1,i2))/(za(i5,i1))/(
     &    zb(i5,i6))*za(i3,i2)*za(i6,i1)*zb(i5,i2)*zb(i5,i4) -1./(za(i5
     &    ,i1))/(za(i5,i6))/(zb(i5,i6))*za(i6,i1)*za(i6,i3)*zb(i1,i4)*
     &    zb(i5,i2) +1./(za(i5,i1))/(zb(i5,i6))*za(i6,i3)*zb(i5,i2)*zb(
     &    i5,i4) )
      aamp(2,1,2) = aamp(2,1,2) + mc**3*four*p256**(-1)*invtwog2Dc * (
     &     -1./(za(i1,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,
     &    i6))*za(i1,i3)*za(i6,i1)**2*zb(i1,i4)*zb(i5,i1)**2 -1./(za(i1
     &    ,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,
     &    i3)*za(i6,i1)**2*zb(i1,i4)*zb(i5,i1)*zb(i5,i4) +1./(za(i1,i2)
     &    )/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i6,i1)*za(i6,i3)*zb(
     &    i1,i4)*zb(i5,i1)*zb(i5,i2) )
      aamp(2,1,2) = aamp(2,1,2) + mc**3*four*invtwog2Dc * (1./(za(i1,i2
     &    ))/(za(i1,i2))/(za(i5,i1))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6)
     &    )*za(i1,i3)*za(i6,i1)**2*zb(i1,i4)*zb(i5,i1) -1./(za(i1,i2))
     &    /(za(i1,i2))/(za(i5,i1))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*
     &    za(i6,i1)*zb(i5,i1)*zb(i5,i4) )

      aamp(2,2,1) =  + four*p256**(-1) * (  -1./(za(i5,i1))/(za(i5,i6))
     &    *za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i6,i1) -1./(za(i5,i1))/(za(
     &    i5,i6))*za(i1,i2)*za(i4,i3)*zb(i1,i4)*zb(i6,i4) -1./(za(i5,i6
     &    ))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i5,i1) -1./(
     &    za(i5,i6))/(za(i6,i1))*za(i1,i2)*za(i4,i3)*zb(i1,i4)*zb(i5,i4
     &    ) )
      aamp(2,2,1) = aamp(2,2,1) + four*p256**(-1)*invtwog2Dc * (  -1./(
     &    za(i5,i1))/(za(i6,i1))*za(i1,i2)**2*za(i1,i3)*zb(i1,i4)*zb(i5
     &    ,i1)*zb(i6,i2) -1./(za(i5,i1))/(za(i6,i1))*za(i1,i2)**2*za(i4
     &    ,i3)*zb(i1,i4)*zb(i5,i4)*zb(i6,i2) )
      aamp(2,2,1) = aamp(2,2,1) + four*p156**(-1) * (  -1./(za(i5,i1))*
     &    za(i3,i2)*zb(i5,i6)*zb(i6,i4) -1./(za(i5,i1))/(za(i5,i6))*za(
     &    i3,i2)*za(i6,i1)*zb(i6,i1)*zb(i6,i4) -1./(za(i5,i6))*za(i3,i2
     &    )*zb(i5,i1)*zb(i6,i4) -1./(za(i5,i6))*za(i3,i2)*zb(i5,i4)*zb(
     &    i6,i1) -1./(za(i5,i6))/(za(i6,i1))*za(i3,i2)*za(i5,i1)*zb(i5,
     &    i1)*zb(i5,i4) -1./(za(i6,i1))*za(i3,i2)*zb(i5,i4)*zb(i5,i6) )
      aamp(2,2,1) = aamp(2,2,1) + four*invtwog2Dc * (1./(za(i5,i1))/(
     &    za(i6,i1))*za(i1,i2)*za(i3,i2)*zb(i5,i4)*zb(i6,i2) )
      aamp(2,2,1) = aamp(2,2,1) + mc**2*four*p256**(-1) * (1./(za(i5,i1
     &    ))/(za(i5,i6))/(zb(i1,i2))*za(i1,i3)*zb(i1,i4)*zb(i6,i1) +
     &   1./(za(i5,i6))/(za(i6,i1))/(zb(i1,i2))*za(i1,i3)*zb(i1,i4)*zb(
     &    i5,i1) )
      aamp(2,2,1) = aamp(2,2,1) + mc**2*four*p256**(-1)*invtwog2Dc * (
     &   1./(za(i5,i1))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i5
     &    ,i6) +1./(za(i5,i1))/(za(i6,i1))/(zb(i1,i2))*za(i1,i2)*za(i1,
     &    i3)*zb(i1,i4)*zb(i5,i2)*zb(i6,i1) -1./(za(i5,i1))/(zb(i1,i2))
     &    *za(i1,i3)*zb(i1,i4)*zb(i5,i6)*zb(i6,i1) )

      bamp(1,1,1) =  + four*p256**(-1) * (1./(zb(i5,i1))/(zb(i5,i6))*
     &    za(i4,i3)*za(i6,i2)*zb(i1,i4)**2 +1./(zb(i5,i6))/(zb(i6,i1))*
     &    za(i4,i3)*za(i5,i2)*zb(i1,i4)**2 )
      bamp(1,1,1) = bamp(1,1,1) + four*p256**(-1)*invtwog1Dc * (  -1./(
     &    zb(i5,i1))/(zb(i6,i1))*za(i4,i3)*za(i5,i2)*za(i6,i2)*zb(i1,i2
     &    )*zb(i1,i4)**2 -1./(zb(i6,i1))*za(i4,i3)*za(i5,i2)*za(i5,i6)*
     &    zb(i1,i4)**2 )

      bamp(2,2,2) =  + mc*four*p256**(-1) * (  -1./(za(i5,i1))/(za(i5,
     &    i6))*za(i1,i3)*zb(i1,i4)*zb(i6,i2) -1./(za(i5,i6))/(za(i6,i1)
     &    )*za(i1,i3)*zb(i1,i4)*zb(i5,i2) )
      bamp(2,2,2) = bamp(2,2,2) + mc*four*p256**(-1)*invtwog1Dc * (1./(
     &    za(i5,i1))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i5,i2
     &    )*zb(i6,i2) +1./(za(i6,i1))*za(i1,i3)*zb(i1,i4)*zb(i5,i2)*zb(
     &    i5,i6) )
      bamp(2,2,2) = bamp(2,2,2) + mc*four*p156**(-1) * (  -1./(za(i1,i2
     &    ))/(za(i5,i1))*za(i1,i3)*zb(i5,i6)*zb(i6,i4) -1./(za(i1,i2))
     &    /(za(i5,i1))/(za(i5,i6))*za(i1,i3)*za(i6,i1)*zb(i6,i1)*zb(i6,
     &    i4) -1./(za(i1,i2))/(za(i5,i6))*za(i1,i3)*zb(i5,i1)*zb(i6,i4)
     &     -1./(za(i1,i2))/(za(i5,i6))*za(i1,i3)*zb(i5,i4)*zb(i6,i1) -
     &   1./(za(i1,i2))/(za(i5,i6))/(za(i6,i1))*za(i1,i3)*za(i5,i1)*zb(
     &    i5,i1)*zb(i5,i4) -1./(za(i1,i2))/(za(i6,i1))*za(i1,i3)*zb(i5,
     &    i4)*zb(i5,i6) )
      bamp(2,2,2) = bamp(2,2,2) + mc*four*invtwog1Dc * (  -1./(za(i5,i1
     &    ))/(za(i6,i1))*za(i1,i3)*zb(i5,i2)*zb(i6,i4) )

      bamp(2,1,1) =  + four*p256**(-1)*invtwog1Dc * (1./(za(i5,i6))/(
     &    zb(i5,i6))*za(i1,i3)*za(i6,i2)**2*zb(i1,i4)*zb(i5,i1)*zb(i5,
     &    i2) +1./(za(i5,i6))/(zb(i5,i6))*za(i4,i3)*za(i6,i2)**2*zb(i1,
     &    i4)*zb(i5,i2)*zb(i5,i4) )
      bamp(2,1,1) = bamp(2,1,1) + four*p156**(-1) * (  -1./(za(i5,i6))
     &    /(zb(i5,i6))/(zb(i6,i1))*za(i3,i2)*za(i6,i1)*zb(i1,i4)*zb(i5,
     &    i1)**2 +1./(zb(i5,i6))/(zb(i6,i1))*za(i3,i2)*zb(i5,i1)**2*zb(
     &    i5,i4) )
      bamp(2,1,1) = bamp(2,1,1) + four*invtwog1Dc * (1./(za(i5,i6))/(
     &    zb(i5,i6))/(zb(i6,i1))*za(i3,i2)*za(i6,i2)*zb(i1,i4)*zb(i5,i1
     &    )*zb(i5,i2) )
      bamp(2,1,1) = bamp(2,1,1) + mc**2*four*p256**(-1)*invtwog1Dc * (
     &     -1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3
     &    )*za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i1)**2 -1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,i3)*za(i6,i1)*za(
     &    i6,i2)*zb(i1,i4)*zb(i5,i1)*zb(i5,i4) -1./(za(i5,i6))/(zb(i1,
     &    i2))/(zb(i5,i6))*za(i6,i2)*za(i6,i3)*zb(i1,i4)*zb(i5,i1)*zb(
     &    i5,i2) )
      bamp(2,1,1) = bamp(2,1,1) + mc**2*four*invtwog1Dc * (1./(za(i1,i2
     &    ))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))/(zb(i6,i1))*za(i1,i3)*
     &    za(i6,i2)*zb(i1,i4)*zb(i5,i1)**2 -1./(za(i5,i6))/(zb(i1,i2))
     &    /(zb(i5,i6))/(zb(i6,i1))*za(i6,i3)*zb(i1,i4)*zb(i5,i1)**2 )
      bamp(2,1,1) = bamp(2,1,1) + mc**4*four*p256**(-1)*invtwog1Dc * (
     &   1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i1,i2))/(zb(i5,i6))
     &    *za(i6,i1)*za(i6,i3)*zb(i1,i4)*zb(i5,i1)**2 )

      bamp(1,2,1) =  + four*p256**(-1)*invtwog1Dc * (1./(za(i5,i6))/(
     &    zb(i5,i6))*za(i1,i3)*za(i5,i2)**2*zb(i1,i4)*zb(i6,i1)*zb(i6,
     &    i2) +1./(za(i5,i6))/(zb(i5,i6))*za(i4,i3)*za(i5,i2)**2*zb(i1,
     &    i4)*zb(i6,i2)*zb(i6,i4) )
      bamp(1,2,1) = bamp(1,2,1) + four*p156**(-1) * (  -1./(za(i5,i6))
     &    /(za(i6,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i1)**2*zb(i1,i4)*zb(
     &    i6,i1) -1./(za(i6,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i1)*zb(i6,
     &    i1)*zb(i6,i4) )
      bamp(1,2,1) = bamp(1,2,1) + four*invtwog1Dc * (1./(za(i5,i6))/(
     &    za(i6,i1))*za(i5,i1)*za(i5,i2)*za(i5,i3)*zb(i1,i4) +1./(za(i5
     &    ,i6))/(za(i6,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i1)*za(i5,i2)*
     &    zb(i1,i4)*zb(i6,i2) +1./(za(i6,i1))*za(i5,i2)*za(i5,i3)*zb(i6
     &    ,i4) +1./(za(i6,i1))/(zb(i5,i6))*za(i3,i2)*za(i5,i2)*zb(i6,i2
     &    )*zb(i6,i4) )
      bamp(1,2,1) = bamp(1,2,1) + mc**2*four*p256**(-1)*invtwog1Dc * (
     &     -1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3
     &    )*za(i5,i1)*za(i5,i2)*zb(i1,i4)*zb(i6,i1)**2 -1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,i3)*za(i5,i1)*za(
     &    i5,i2)*zb(i1,i4)*zb(i6,i1)*zb(i6,i4) -1./(za(i5,i6))/(zb(i1,
     &    i2))/(zb(i5,i6))*za(i5,i2)*za(i5,i3)*zb(i1,i4)*zb(i6,i1)*zb(
     &    i6,i2) )
      bamp(1,2,1) = bamp(1,2,1) + mc**2*four*invtwog1Dc * (1./(za(i1,i2
     &    ))/(za(i5,i6))/(za(i6,i1))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*
     &    za(i5,i1)*za(i5,i2)*zb(i1,i4)*zb(i6,i1) +1./(za(i1,i2))/(za(
     &    i6,i1))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*za(i5,i2)*zb(i6,i1)
     &    *zb(i6,i4) -1./(za(i5,i6))/(za(i6,i1))/(zb(i1,i2))/(zb(i5,i6)
     &    )*za(i5,i1)*za(i5,i3)*zb(i1,i4)*zb(i6,i1) -1./(za(i6,i1))/(
     &    zb(i1,i2))/(zb(i5,i6))*za(i5,i3)*zb(i6,i1)*zb(i6,i4) )
      bamp(1,2,1) = bamp(1,2,1) + mc**4*four*p256**(-1)*invtwog1Dc * (
     &   1./(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i1,i2))/(zb(i5,i6))
     &    *za(i5,i1)*za(i5,i3)*zb(i1,i4)*zb(i6,i1)**2 )

      bamp(1,1,2) =  + mc*four*p256**(-1) * (1./(za(i1,i2))/(zb(i5,i1))
     &    /(zb(i5,i6))*za(i4,i3)*za(i6,i1)*zb(i1,i4)**2 +1./(za(i1,i2))
     &    /(zb(i5,i6))/(zb(i6,i1))*za(i4,i3)*za(i5,i1)*zb(i1,i4)**2 -
     &   1./(zb(i5,i1))/(zb(i5,i6))*za(i6,i3)*zb(i1,i2)*zb(i1,i4) -1./(
     &    zb(i5,i6))/(zb(i6,i1))*za(i5,i3)*zb(i1,i2)*zb(i1,i4) )
      bamp(1,1,2) = bamp(1,1,2) + mc*four*p256**(-1)*invtwog1Dc * (  -
     &   1./(za(i1,i2))/(zb(i5,i1))/(zb(i6,i1))*za(i4,i3)*za(i5,i1)*za(
     &    i6,i2)*zb(i1,i2)*zb(i1,i4)**2 -1./(za(i1,i2))/(zb(i6,i1))*za(
     &    i4,i3)*za(i5,i1)*za(i5,i6)*zb(i1,i4)**2 +1./(zb(i5,i1))/(zb(
     &    i6,i1))*za(i4,i3)*za(i5,i6)*zb(i1,i2)*zb(i1,i4)**2 +1./(zb(i5
     &    ,i1))/(zb(i6,i1))*za(i5,i2)*za(i6,i3)*zb(i1,i2)**2*zb(i1,i4)
     &     )

      bamp(1,2,2) =  + mc*four*p256**(-1)*invtwog1Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))*za(i1,i3)*za(i5,i1)*za(i5,i2)*zb(i1,
     &    i4)*zb(i6,i1)*zb(i6,i2) +1./(za(i1,i2))/(za(i5,i6))/(zb(i5,i6
     &    ))*za(i4,i3)*za(i5,i1)*za(i5,i2)*zb(i1,i4)*zb(i6,i2)*zb(i6,i4
     &    ) -1./(za(i5,i6))/(zb(i5,i6))*za(i5,i2)*za(i5,i3)*zb(i1,i4)*
     &    zb(i6,i2)**2 )
      bamp(1,2,2) = bamp(1,2,2) + mc*four*p156**(-1) * (1./(za(i1,i2))
     &    /(za(i5,i6))/(za(i6,i1))/(zb(i5,i6))*za(i1,i3)*za(i5,i1)**2*
     &    zb(i1,i4)*zb(i6,i1) +1./(za(i1,i2))/(za(i6,i1))/(zb(i5,i6))*
     &    za(i1,i3)*za(i5,i1)*zb(i6,i1)*zb(i6,i4) )
      bamp(1,2,2) = bamp(1,2,2) + mc*four*invtwog1Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(za(i6,i1))*za(i5,i1)**2*za(i5,i3)*zb(i1,i4) +
     &   1./(za(i1,i2))/(za(i5,i6))/(za(i6,i1))/(zb(i5,i6))*za(i3,i2)*
     &    za(i5,i1)**2*zb(i1,i4)*zb(i6,i2) +1./(za(i1,i2))/(za(i6,i1))*
     &    za(i5,i1)*za(i5,i3)*zb(i6,i4) +1./(za(i1,i2))/(za(i6,i1))/(
     &    zb(i5,i6))*za(i3,i2)*za(i5,i1)*zb(i6,i2)*zb(i6,i4) -1./(za(i5
     &    ,i6))/(za(i6,i1))/(zb(i5,i6))*za(i5,i1)*za(i5,i3)*zb(i1,i4)*
     &    zb(i6,i2) -1./(za(i6,i1))/(zb(i5,i6))*za(i5,i3)*zb(i6,i2)*zb(
     &    i6,i4) )
      bamp(1,2,2) = bamp(1,2,2) + mc**3*four*p256**(-1)*invtwog1Dc * (
     &     -1./(za(i1,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,
     &    i6))*za(i1,i3)*za(i5,i1)**2*zb(i1,i4)*zb(i6,i1)**2 -1./(za(i1
     &    ,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,
     &    i3)*za(i5,i1)**2*zb(i1,i4)*zb(i6,i1)*zb(i6,i4) +1./(za(i1,i2)
     &    )/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i5,i1)*za(i5,i3)*zb(
     &    i1,i4)*zb(i6,i1)*zb(i6,i2) )
      bamp(1,2,2) = bamp(1,2,2) + mc**3*four*invtwog1Dc * (1./(za(i1,i2
     &    ))/(za(i1,i2))/(za(i5,i6))/(za(i6,i1))/(zb(i1,i2))/(zb(i5,i6)
     &    )*za(i1,i3)*za(i5,i1)**2*zb(i1,i4)*zb(i6,i1) +1./(za(i1,i2))
     &    /(za(i1,i2))/(za(i6,i1))/(zb(i1,i2))/(zb(i5,i6))*za(i1,i3)*
     &    za(i5,i1)*zb(i6,i1)*zb(i6,i4) )

      bamp(2,1,2) =  + mc*four*p256**(-1)*invtwog1Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))*za(i1,i3)*za(i6,i1)*za(i6,i2)*zb(i1,
     &    i4)*zb(i5,i1)*zb(i5,i2) +1./(za(i1,i2))/(za(i5,i6))/(zb(i5,i6
     &    ))*za(i4,i3)*za(i6,i1)*za(i6,i2)*zb(i1,i4)*zb(i5,i2)*zb(i5,i4
     &    ) -1./(za(i5,i6))/(zb(i5,i6))*za(i6,i2)*za(i6,i3)*zb(i1,i4)*
     &    zb(i5,i2)**2 )
      bamp(2,1,2) = bamp(2,1,2) + mc*four*p156**(-1) * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))/(zb(i6,i1))*za(i1,i3)*za(i6,i1)*zb(
     &    i1,i4)*zb(i5,i1)**2 -1./(za(i1,i2))/(zb(i5,i6))/(zb(i6,i1))*
     &    za(i1,i3)*zb(i5,i1)**2*zb(i5,i4) )
      bamp(2,1,2) = bamp(2,1,2) + mc*four*invtwog1Dc * (1./(za(i1,i2))
     &    /(za(i5,i6))/(zb(i5,i6))/(zb(i6,i1))*za(i3,i2)*za(i6,i1)*zb(
     &    i1,i4)*zb(i5,i1)*zb(i5,i2) -1./(za(i5,i6))/(zb(i5,i6))/(zb(i6
     &    ,i1))*za(i6,i3)*zb(i1,i4)*zb(i5,i1)*zb(i5,i2) )
      bamp(2,1,2) = bamp(2,1,2) + mc**3*four*p256**(-1)*invtwog1Dc * (
     &     -1./(za(i1,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,
     &    i6))*za(i1,i3)*za(i6,i1)**2*zb(i1,i4)*zb(i5,i1)**2 -1./(za(i1
     &    ,i2))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i4,
     &    i3)*za(i6,i1)**2*zb(i1,i4)*zb(i5,i1)*zb(i5,i4) +1./(za(i1,i2)
     &    )/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))*za(i6,i1)*za(i6,i3)*zb(
     &    i1,i4)*zb(i5,i1)*zb(i5,i2) )
      bamp(2,1,2) = bamp(2,1,2) + mc**3*four*invtwog1Dc * (1./(za(i1,i2
     &    ))/(za(i1,i2))/(za(i5,i6))/(zb(i1,i2))/(zb(i5,i6))/(zb(i6,i1)
     &    )*za(i1,i3)*za(i6,i1)*zb(i1,i4)*zb(i5,i1)**2 )

      bamp(2,2,1) =  + four*p256**(-1) * (1./(za(i5,i1))/(za(i5,i6))*
     &    za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i6,i1) +1./(za(i5,i1))/(za(
     &    i5,i6))*za(i1,i2)*za(i4,i3)*zb(i1,i4)*zb(i6,i4) +1./(za(i5,i6
     &    ))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*zb(i5,i1) +1./(
     &    za(i5,i6))/(za(i6,i1))*za(i1,i2)*za(i4,i3)*zb(i1,i4)*zb(i5,i4
     &    ) )
      bamp(2,2,1) = bamp(2,2,1) + four*p256**(-1)*invtwog1Dc * (  -1./(
     &    za(i5,i1))/(za(i6,i1))*za(i1,i2)**2*za(i1,i3)*zb(i1,i4)*zb(i5
     &    ,i2)*zb(i6,i1) -1./(za(i5,i1))/(za(i6,i1))*za(i1,i2)**2*za(i4
     &    ,i3)*zb(i1,i4)*zb(i5,i2)*zb(i6,i4) )
      bamp(2,2,1) = bamp(2,2,1) + four*p156**(-1) * (1./(za(i5,i1))*za(
     &    i3,i2)*zb(i5,i6)*zb(i6,i4) +1./(za(i5,i1))/(za(i5,i6))*za(i3,
     &    i2)*za(i6,i1)*zb(i6,i1)*zb(i6,i4) +1./(za(i5,i6))*za(i3,i2)*
     &    zb(i5,i1)*zb(i6,i4) +1./(za(i5,i6))*za(i3,i2)*zb(i5,i4)*zb(i6
     &    ,i1) +1./(za(i5,i6))/(za(i6,i1))*za(i3,i2)*za(i5,i1)*zb(i5,i1
     &    )*zb(i5,i4) +1./(za(i6,i1))*za(i3,i2)*zb(i5,i4)*zb(i5,i6) )
      bamp(2,2,1) = bamp(2,2,1) + four*invtwog1Dc * (1./(za(i5,i1))/(
     &    za(i6,i1))*za(i1,i2)*za(i3,i2)*zb(i5,i2)*zb(i6,i4) )
      bamp(2,2,1) = bamp(2,2,1) + mc**2*four*p256**(-1) * (  -1./(za(i5
     &    ,i1))/(za(i5,i6))/(zb(i1,i2))*za(i1,i3)*zb(i1,i4)*zb(i6,i1)
     &     -1./(za(i5,i6))/(za(i6,i1))/(zb(i1,i2))*za(i1,i3)*zb(i1,i4)*
     &    zb(i5,i1) )
      bamp(2,2,1) = bamp(2,2,1) + mc**2*four*p256**(-1)*invtwog1Dc * (
     &     -1./(za(i5,i1))/(za(i6,i1))*za(i1,i2)*za(i1,i3)*zb(i1,i4)*
     &    zb(i5,i6) +1./(za(i5,i1))/(za(i6,i1))/(zb(i1,i2))*za(i1,i2)*
     &    za(i1,i3)*zb(i1,i4)*zb(i5,i1)*zb(i6,i2) +1./(za(i6,i1))/(zb(
     &    i1,i2))*za(i1,i3)*zb(i1,i4)*zb(i5,i1)*zb(i5,i6) )


      do h1=1,2
      do h2=1,2
      do hf=1,2
      aamp(h1,h2,hf)=aamp(h1,h2,hf)/s34
      bamp(h1,h2,hf)=bamp(h1,h2,hf)/s34
      enddo
      enddo
      enddo

      return
      end
