      subroutine spassign(p,k1,k2,k3,q,l,l_)
      implicit none
      include 'types.f'
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'debr.f'

      HQ1H=za(q,k1)
      HQ2H=za(q,k2)
      HQ3H=za(q,k3)
      HQPH=za(q,p)


      HP1H=za(p,k1)
      HP2H=za(p,k2)
      HP3H=za(p,k3)
      HPQH=za(p,q)

      H1QH=za(k1,q)
      H2QH=za(k2,q)
      H3QH=za(k3,q)

      H1PH=za(k1,p)
      H2PH=za(k2,p)
      H3PH=za(k3,p)

      H12H=za(k1,k2)
      H13H=za(k1,k3)
      H23H=za(k2,k3)

      H21H=za(k2,k1)
      H31H=za(k3,k1)
      H32H=za(k3,k2)


      H1LH=za(k1,l)
      H2LH=za(k2,l)
      H3LH=za(k3,l)

      HL1H=za(l,k1)
      HL2H=za(l,k2)
      HL3H=za(l,k3)

      H1L_H=za(k1,l_)
      H2L_H=za(k2,l_)
      H3L_H=za(k3,l_)

      HL_1H=za(l_,k1)
      HL_2H=za(l_,k2)
      HL_3H=za(l_,k3)

      HL_QH=za(l_,q)
      HL_PH=za(l_,p)

      HQL_H=za(q,l_)
      HPL_H=za(p,l_)


      TQ1T=zb(q,k1)
      TQ2T=zb(q,k2)
      TQ3T=zb(q,k3)
      TQPT=zb(q,p)

      TP1T=zb(p,k1)
      TP2T=zb(p,k2)
      TP3T=zb(p,k3)
      TPQT=zb(p,q)

      T1QT=zb(k1,q)
      T2QT=zb(k2,q)
      T3QT=zb(k3,q)

      T1PT=zb(k1,p)
      T2PT=zb(k2,p)
      T3PT=zb(k3,p)


      T12T=zb(k1,k2)
      T13T=zb(k1,k3)
      T23T=zb(k2,k3)

      T21T=zb(k2,k1)
      T31T=zb(k3,k1)
      T32T=zb(k3,k2)

      T1L_T=zb(k1,l_)
      T2L_T=zb(k2,l_)
      T3L_T=zb(k3,l_)

      T1LT=zb(k1,l)
      T2LT=zb(k2,l)
      T3LT=zb(k3,l)

      TL_1T=zb(l_,k1)
      TL_2T=zb(l_,k2)
      TL_3T=zb(l_,k3)

      TL_QT=zb(l_,q)
      TL_PT=zb(l_,p)

      TLQT=zb(l,q)
      TLPT=zb(l,p)

      TQLT=zb(q,l)
      TPLT=zb(p,l)


      S12=s(k1,k2)
      S13=s(k1,k3)
      S23=s(k2,k3)

      SQ1=s(q,k1)
      SQ2=s(q,k2)
      SQ3=s(q,k3)

      SP1=s(p,k1)
      SP2=s(p,k2)
      SP3=s(p,k3)
      return
      end
