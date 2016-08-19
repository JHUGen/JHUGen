      double complex function q4ghpmp1(p1,p2,p3,p4,p5,za,zb)
      implicit none
C     This is the reduced matrix element squared
C     for the process
c     q(p1+)+qbar(p2) --> H((p5+p6)+Q(p3-)+qbar(p4)+g(p5+)
      include 'constants.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer p1,p2,p3,p4,p5
      double precision s5h,t123,t234,t345,t125
      integer i1,i2,i3,i4
      double complex zab,zaba

      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zaba(i1,i2,i3,i4)=+za(i1,i2)*zb(i2,i3)*za(i3,i4)


      t123=s(p1,p2)+s(p2,p3)+s(p3,p1)
      t234=s(p2,p3)+s(p3,p4)+s(p4,p2)
      t125=s(p1,p2)+s(p2,p5)+s(p5,p1)
      t345=s(p3,p4)+s(p4,p5)+s(p5,p3)
      s5h=s(p1,p2)+s(p1,p3)+s(p1,p4)
     .            +s(p2,p3)+s(p2,p4)
     .                     +s(p3,p4)


C %\cite{DelDuca:2004wt}
C \bibitem{DelDuca:2004wt}
C V.~Del Duca, A.~Frizzo and F.~Maltoni,
C %``Higgs boson production in association with three jets,''
C JHEP {\bf 0405}, 064 (2004)
C [arXiv:hep-ph/0404013].
C %%CITATION = HEP-PH 0404013;%%
C---B13
      q4ghpmp1=
     . +1d0/(za(p1,p5)*s(p1,p2)*s(p3,p4)*t345)
     . *(za(p2,p3)*zb(p1,p2)*(-za(p2,p3)*zb(p4,p5)
     . *(zab(p1,p4,p3)+zab(p1,p5,p3))
     . +za(p1,p3)*zb(p3,p5)*(zab(p2,p3,p4)+zab(p2,p5,p4)))
     . +za(p1,p2)*(zab(p3,p4,p1)+zab(p3,p5,p1))
     . *(-za(p1,p3)*zb(p1,p4)*zb(p3,p5)+zb(p4,p5)*(s(p1,p4)+s(p1,p5)))
     . -s(p1,p2)*s(p4,p5)*za(p2,p3)*zb(p4,p5))
      q4ghpmp1=q4ghpmp1
     . -1d0/(s(p1,p2)*s(p3,p4)*s5h)
     . *(-za(p1,p2)*zb(p1,p5)**2*(zab(p3,p1,p4)+zab(p3,p2,p4))
     . -zb(p1,p5)*zb(p4,p5)*(2d0*zaba(p2,p1,p4,p3)
     . -za(p2,p3)*(s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4)))
     . -za(p3,p4)*zb(p4,p5)**2*(zab(p2,p3,p1)+zab(p2,p4,p1)))
      q4ghpmp1=q4ghpmp1
     . +1d0/(za(p1,p5)*s(p1,p2)*s(p3,p4)*t125)
     . *(-za(p1,p2)**2*za(p3,p4)*zb(p1,p4)**2*zb(p2,p5)
     . +za(p1,p2)*(-za(p2,p3)**2*zb(p1,p2)*zb(p2,p5)*zb(p3,p4)
     . +za(p3,p4)*zb(p1,p4)*zb(p4,p5)*(s(p2,p5)-s(p1,p5)))
     . +zb(p1,p5)*(-za(p1,p3)*za(p2,p3)*zb(p3,p4)*(s(p1,p5)+s(p2,p5))
     . +za(p1,p5)*za(p2,p5)*za(p3,p4)*zb(p4,p5)**2))
      q4ghpmp1=q4ghpmp1
     . +1d0/(za(p1,p5)*s(p1,p2)*s(p3,p4))
     . *(za(p2,p3)*zb(p1,p5)*(zab(p1,p3,p4)-zab(p1,p5,p4))
     . +za(p1,p2)*(za(p1,p3)*zb(p1,p4)*zb(p1,p5)-zb(p4,p5)
     . *(zab(p3,p2,p1)-zab(p3,p4,p1))))
       q4ghpmp1=q4ghpmp1
     . -1d0/(za(p1,p5)*za(p4,p5)*s(p1,p2)*t345)
     . *(-za(p2,p3)**2*zb(p1,p2)*(zab(p1,p4,p3)+zab(p1,p5,p3))
     . +za(p1,p2)*(s(p1,p4)+s(p1,p5))*(zab(p3,p4,p1)+zab(p3,p5,p1))
     . -s(p1,p2)*s(p4,p5)*za(p2,p3))
       q4ghpmp1=q4ghpmp1
     . +za(p2,p3)*zb(p4,p5)/(s(p1,p2)*t123*s5h)
     . *(-zb(p1,p5)*t123-zb(p4,p5)*(zab(p4,p2,p1)+zab(p4,p3,p1)))
     . +za(p2,p3)*zb(p1,p5)/(s(p3,p4)*t234*s5h)
     . *(-zb(p1,p5)*(zab(p1,p2,p4)+zab(p1,p3,p4))-zb(p4,p5)*t234)
       q4ghpmp1=q4ghpmp1
     . +1d0/(za(p1,p5)*s(p3,p4)*t125)
     . *(za(p2,p3)*zb(p3,p4)*(zab(p3,p1,p5)+zab(p3,p2,p5))
     . +za(p3,p4)*zb(p4,p5)*(zab(p2,p1,p4)+zab(p2,p5,p4)))


      return
      end
