      subroutine ggttww1n(in,s1t,s2t,s12,c6,c8,res)
      implicit none
      include 'types.f'

c----Matrix element for tt production
C----averaged over initial colours and spins
c    line in contracted with the vector n(mu)
C in is the label of the contracted line
c     g(-p1)+g(-p2)--> t(p3,p4,p5)+tb(p6,p7,p8)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      include 'masses.f'
      complex(dp):: zanb(mxpart,mxpart)
      common/zprodsgvec/zanb
!$omp threadprivate(/zprodsgvec/)
      integer:: in,p1,p2,q3,q4,q5,q6,q7,q8,icol,j
      real(dp):: s1t,s2t,s12,c6,c8,mtsq,res(0:2)
      complex(dp):: loab(2),loba(2),loqed(2)
      parameter(p1=1,p2=2,q3=3,q4=4,q5=5,q6=6,q7=7,q8=8)

      mtsq=mt**2

      if  (in == 1) then
        loab(1)= + s2t**(-1) * ( za(q3,p2)*za(p1,p2)*za(q5,q7)*zb(q3,q4
     &    )*zb(p1,p2)*zanb(p2,q5)*s12**(-1)*c8 + za(q3,p2)*za(p1,p2)*
     &    za(q5,q7)*zb(q3,q4)*zb(p1,q8)*zanb(q8,q5)*s12**(-1) )
      loab(1) = loab(1) + s2t**(-1)*mtsq * (  - za(q3,p2)*za(p1,p2)*zb(
     &    q3,q4)*zanb(q7,p1)*s12**(-1) + za(p1,p2)*za(p2,q8)*zb(q4,p1)*
     &    zanb(q7,q8)*s12**(-1) - za(p1,p2)*za(q5,q7)*zb(q4,p1)*zanb(p2
     &    ,q5)*s12**(-1) )
      loab(1) = loab(1) + mtsq * ( za(p1,p2)*za(p1,q7)*zb(q4,p1)*zanb(
     &    p2,p1)*s12**(-2) + za(p1,p2)*za(p2,q7)*zb(q4,p1)*zanb(p2,p2)*
     &    s12**(-2) )
      loab(1) = loab(1) + za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     & ,q5)*zanb(p2,p1)*s12**(-2) + za(q3,p2)*za(p1,p2)*za(q5,q7)*zb(q3
     &    ,q4)*zb(p1,q5)*zanb(p2,p2)*s12**(-2)

        loab(2)= + s2t**(-1) * (  - za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p2,q8)*zanb(q8,q5)*s12**(-1) )
      loab(2) = loab(2) + s2t**(-1)*mtsq * ( za(q3,p1)*zb(q3,q4)*zb(p1,
     &    p2)*zanb(q7,p2)*s12**(-1) - za(p1,p2)*zb(q4,p2)*zb(p1,p2)*
     &    zanb(q7,p2)*s12**(-1)*c8 - za(p1,q8)*zb(q4,p2)*zb(p1,p2)*
     &    zanb(q7,q8)*s12**(-1) + za(q5,q7)*zb(q4,p2)*zb(p1,p2)*zanb(p1
     &    ,q5)*s12**(-1) )
      loab(2) = loab(2) + mtsq * (  - za(p1,q7)*zb(q4,p1)*zb(p1,p2)*
     &    zanb(p1,p2)*s12**(-2) - za(p1,q7)*zb(q4,p2)*zb(p1,p2)*zanb(p2
     &    ,p2)*s12**(-2) )
      loab(2) = loab(2) - za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1,p2)*zb(p1
     & ,q5)*zanb(p1,p2)*s12**(-2) - za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p2,q5)*zanb(p2,p2)*s12**(-2)

        loba(1)= + s1t**(-1) * ( za(p1,p2)**2*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,q5)*zanb(q3,p1)*s12**(-1)*c6 - za(p1,p2)*za(p2,q6)*za(q5,q7)
     &    *zb(q3,q4)*zb(p1,q5)*zanb(q3,q6)*s12**(-1) )
      loba(1) = loba(1) + s1t**(-1)*mtsq * (  - za(p1,p2)*za(p2,q7)*zb(
     &    q3,q4)*zanb(q3,p1)*s12**(-1) - za(p1,p2)*za(p2,q7)*zb(p1,q6)*
     &    zanb(q6,q4)*s12**(-1) - za(p1,p2)*za(q5,q7)*zb(p1,q5)*zanb(p2
     &    ,q4)*s12**(-1) )
      loba(1) = loba(1) + mtsq * (  - za(p1,p2)*za(p1,q7)*zb(q4,p1)*
     &    zanb(p2,p1)*s12**(-2) - za(p1,p2)*za(p2,q7)*zb(q4,p1)*zanb(p2
     &    ,p2)*s12**(-2) )
      loba(1) = loba(1) - za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     & ,q5)*zanb(p2,p1)*s12**(-2) - za(q3,p2)*za(p1,p2)*za(q5,q7)*zb(q3
     &    ,q4)*zb(p1,q5)*zanb(p2,p2)*s12**(-2)

        loba(2)= + s1t**(-1) * ( za(p1,q6)*za(q5,q7)*zb(q3,q4)*zb(p1,p2
     &    )*zb(p2,q5)*zanb(q3,q6)*s12**(-1) )
      loba(2) = loba(2) + s1t**(-1)*mtsq * ( za(p1,q7)*zb(q3,q4)*zb(p1,
     &    p2)*zanb(q3,p2)*s12**(-1) - za(p1,q7)*zb(p1,p2)**2*zanb(p1,q4
     &    )*s12**(-1)*c6 + za(p1,q7)*zb(p1,p2)*zb(p2,q6)*zanb(q6,q4)*
     &    s12**(-1) + za(q5,q7)*zb(p1,p2)*zb(p2,q5)*zanb(p1,q4)*
     &    s12**(-1) )
      loba(2) = loba(2) + mtsq * ( za(p1,q7)*zb(q4,p1)*zb(p1,p2)*zanb(
     &    p1,p2)*s12**(-2) + za(p1,q7)*zb(q4,p2)*zb(p1,p2)*zanb(p2,p2)*
     &    s12**(-2) )
      loba(2) = loba(2) + za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1,p2)*zb(p1
     & ,q5)*zanb(p1,p2)*s12**(-2) + za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p2,q5)*zanb(p2,p2)*s12**(-2)

      elseif (in == 2) then
        loab(1)= + s2t**(-1) * ( za(p1,p2)**2*za(q5,q7)*zb(q3,q4)*zb(p2
     &    ,q5)*zanb(q3,p2)*s12**(-1)*c8 + za(p1,p2)*za(p1,q8)*za(q5,q7)
     &    *zb(q3,q4)*zb(p2,q5)*zanb(q3,q8)*s12**(-1) )
      loab(1) = loab(1) + s2t**(-1)*mtsq * ( za(p1,p2)*za(p1,q7)*zb(q3,
     &    q4)*zanb(q3,p2)*s12**(-1) + za(p1,p2)*za(p1,q7)*zb(p2,q8)*
     &    zanb(q8,q4)*s12**(-1) + za(p1,p2)*za(q5,q7)*zb(p2,q5)*zanb(p1
     &    ,q4)*s12**(-1) )
      loab(1) = loab(1) + mtsq * (  - za(p1,p2)*za(p1,q7)*zb(q4,p1)*
     &    zanb(p1,p2)*s12**(-2) + za(p1,p2)*za(p1,q7)*zb(q4,p2)*zanb(p1
     &    ,p1)*s12**(-2) )
      loab(1) = loab(1) - za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     & ,q5)*zanb(p1,p2)*s12**(-2) + za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3
     &    ,q4)*zb(p2,q5)*zanb(p1,p1)*s12**(-2)

        loab(2)= + s2t**(-1) * (  - za(p2,q8)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p1,q5)*zanb(q3,q8)*s12**(-1) )
      loab(2) = loab(2) + s2t**(-1)*mtsq * (  - za(p2,q7)*zb(q3,q4)*zb(
     &    p1,p2)*zanb(q3,p1)*s12**(-1) - za(p2,q7)*zb(p1,p2)**2*zanb(p2
     &    ,q4)*s12**(-1)*c8 - za(p2,q7)*zb(p1,p2)*zb(p1,q8)*zanb(q8,q4)
     &    *s12**(-1) - za(q5,q7)*zb(p1,p2)*zb(p1,q5)*zanb(p2,q4)*
     &    s12**(-1) )
      loab(2) = loab(2) + mtsq * ( za(p1,q7)*zb(q4,p1)*zb(p1,p2)*zanb(
     &    p2,p1)*s12**(-2) - za(p2,q7)*zb(q4,p1)*zb(p1,p2)*zanb(p1,p1)*
     &    s12**(-2) )
      loab(2) = loab(2) + za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1,p2)*zb(p1
     & ,q5)*zanb(p2,p1)*s12**(-2) - za(q3,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p1,q5)*zanb(p1,p1)*s12**(-2)

        loba(1)= + s1t**(-1) * ( za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3,q4
     &    )*zb(p1,p2)*zanb(p1,q5)*s12**(-1)*c6 - za(q3,p1)*za(p1,p2)*
     &    za(q5,q7)*zb(q3,q4)*zb(p2,q6)*zanb(q6,q5)*s12**(-1) )
      loba(1) = loba(1) + s1t**(-1)*mtsq * ( za(q3,p1)*za(p1,p2)*zb(q3,
     &    q4)*zanb(q7,p2)*s12**(-1) - za(p1,p2)*za(p1,q6)*zb(q4,p2)*
     &    zanb(q7,q6)*s12**(-1) + za(p1,p2)*za(q5,q7)*zb(q4,p2)*zanb(p1
     &    ,q5)*s12**(-1) )
      loba(1) = loba(1) + mtsq * ( za(p1,p2)*za(p1,q7)*zb(q4,p1)*zanb(
     &    p1,p2)*s12**(-2) - za(p1,p2)*za(p1,q7)*zb(q4,p2)*zanb(p1,p1)*
     &    s12**(-2) )
      loba(1) = loba(1) + za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     & ,q5)*zanb(p1,p2)*s12**(-2) - za(q3,p1)*za(p1,p2)*za(q5,q7)*zb(q3
     &    ,q4)*zb(p2,q5)*zanb(p1,p1)*s12**(-2)

        loba(2)= + s1t**(-1) * ( za(q3,p2)*za(q5,q7)*zb(q3,q4)*zb(p1,p2
     &    )*zb(p1,q6)*zanb(q6,q5)*s12**(-1) )
      loba(2) = loba(2) + s1t**(-1)*mtsq * (  - za(q3,p2)*zb(q3,q4)*zb(
     &    p1,p2)*zanb(q7,p1)*s12**(-1) - za(p1,p2)*zb(q4,p1)*zb(p1,p2)*
     &    zanb(q7,p1)*s12**(-1)*c6 + za(p2,q6)*zb(q4,p1)*zb(p1,p2)*
     &    zanb(q7,q6)*s12**(-1) - za(q5,q7)*zb(q4,p1)*zb(p1,p2)*zanb(p2
     &    ,q5)*s12**(-1) )
      loba(2) = loba(2) + mtsq * (  - za(p1,q7)*zb(q4,p1)*zb(p1,p2)*
     &    zanb(p2,p1)*s12**(-2) + za(p2,q7)*zb(q4,p1)*zb(p1,p2)*zanb(p1
     &    ,p1)*s12**(-2) )
      loba(2) = loba(2) - za(q3,p1)*za(q5,q7)*zb(q3,q4)*zb(p1,p2)*zb(p1
     & ,q5)*zanb(p2,p1)*s12**(-2) + za(q3,p2)*za(q5,q7)*zb(q3,q4)*zb(p1
     &    ,p2)*zb(p1,q5)*zanb(p1,p1)*s12**(-2)

      else
        write(6,*) 'ggttww1n: Unimplemented value of in',in
        stop
      endif

      do icol=0,2
      res(icol)=zip
      enddo

      do j=1,2
        loqed(j)=loab(j)+loba(j)

        res(1)=res(1)+abs(loab(j))**2
        res(2)=res(2)+abs(loba(j))**2
        res(0)=res(0)-1._dp/xnsq*abs(loqed(j))**2
      enddo

      return
      end
