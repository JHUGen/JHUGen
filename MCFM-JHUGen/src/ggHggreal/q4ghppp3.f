      function q4ghppp3(p1,p2,p3,p4,p5,za,zb)

C     This is the reduced matrix element squared
C     for the process
c     q(p1)+qbar(p2) --> H((p5+p6)+Q(p3)+qbar(p4)+g(p5)
c     calculated by the program Hq4g.frm
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp)::q4ghppp3
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5
      real(dp)::s5h,t134,t234,t125
      integer::i1,i2,i3,i4
      complex(dp)::zab,zaba,zbab
      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zaba(i1,i2,i3,i4)=+za(i1,i2)*zb(i2,i3)*za(i3,i4)
      zbab(i1,i2,i3,i4)=+zb(i1,i2)*za(i2,i3)*zb(i3,i4)

      t134=s(p1,p3)+s(p3,p4)+s(p4,p1)
      t234=s(p2,p3)+s(p3,p4)+s(p4,p2)
      t125=s(p1,p2)+s(p2,p5)+s(p5,p1)
      s5h=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &            +s(p2,p3)+s(p2,p4)
     &                     +s(p3,p4)
C %\cite{DelDuca:2004wt}
C \bibitem{DelDuca:2004wt}
C V.~Del Duca, A.~Frizzo and F.~Maltoni,
C %``Higgs boson production in association with three jets,''
C JHEP {\bf 0405}, 064 (2004)
C [arXiv:hep-ph/0404013].
C %%CITATION = HEP-PH 0404013;%%
C Eq. B11

      q4ghppp3=
     & -one/(za(p1,p5)*s(p3,p4)*t125)
     & *(za(p2,p4)*zb(p4,p3)*(+zab(p4,p1,p5)+zab(p4,p2,p5))
     & + za(p4,p3)*zb(p3,p5)*(+zab(p2,p1,p3)+zab(p2,p5,p3)))
      q4ghppp3= q4ghppp3
     & +za(p1,p2)/(za(p1,p5)*za(p2,p5)*s(p3,p4)*t125)
     & *(zb(p1,p3)*(zaba(p2,p1,p3,p4)+zaba(p2,p5,p3,p4))
     &  -za(p2,p4)*(zbab(p3,p4,p2,p1)+zbab(p3,p4,p5,p1)))
      q4ghppp3= q4ghppp3
     & +zb(p1,p3)/(s(p3,p4)*t134*s5h)*(zab(p4,p1,p5)+zab(p4,p3,p5))
     & *(-zab(p2,p1,p5)-zab(p2,p3,p5)-zab(p2,p4,p5))
      q4ghppp3= q4ghppp3
     & -za(p2,p4)*zb(p1,p5)/(s(p3,p4)*s5h)
     & *(zb(p5,p1)*(zab(p1,p2,p3)+zab(p1,p4,p3))/t234+zb(p5,p3))


      return
      end
