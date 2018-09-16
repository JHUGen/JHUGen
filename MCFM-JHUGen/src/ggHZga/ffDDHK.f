      function ffDDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: ffDDHK

C-----Author: R.K.Ellis July 2012
C-----The top loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
      real(dp):: s12,s34,msq
      complex(dp):: C0DDHK,C2DDHK
      ffDDHK=C0DDHK(s12,s34,msq)+4._dp*C2DDHK(s12,s34,msq)
      return
      end

      function fWDDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: fWDDHK

C-----Author: R.K.Ellis July 2012
C-----The W-loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
c----- msq here should only be equal to wmass**2
      real(dp):: s12,s34,msq
      complex(dp):: C0DDHK,C2DDHK
      fWDDHK=2._dp*(s12/msq*(1._dp-2._dp*msq/s34)
     & +2._dp*(1._dp-6._dp*msq/s34))*C2DDHK(s12,s34,msq)
     & +4._dp*(1._dp-4._dp*msq/s34)*C0DDHK(s12,s34,msq)
      return
      end

      function f0DDHK(s12,s34)
      implicit none
      include 'types.f'
      complex(dp):: f0DDHK
C-----The F0-function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(2)
C-----suitably generalized to allow off-shell Z-line
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'first.f'
      include 'msbarmasses.f'
      integer:: top
      complex(dp):: ffDDHK,fWDDHK
      real(dp):: s12,s34,cotw,mtsq,mwsq,mt_eff,massfrun
      parameter(top=2)
      save mt_eff
!$omp threadprivate(mt_eff)

      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

      cotw=sqrt((1._dp-xw)/xw)
      mtsq=mt**2
      mwsq=wmass**2
      f0DDHK=s34*(cotw*fWDDHK(s12,s34,mwsq)
     &+2._dp*Q(top)*xn*mt_eff**2/s34*(L(top)+R(top))*ffDDHK(s12,s34,mtsq))
      return
      end


      function qlC2DDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      include 'cplx.h'
      complex(dp):: qlC2DDHK

      include 'scale.f'
      complex(dp):: qlI2,qlI3
      real(dp):: s12,s34,msq
      qlC2DDHK=-cplx1(0.5_dp/(s12-s34))
     & +0.5_dp*s34/(s12-s34)**2
     & *(qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0))
     & -msq/(s12-s34)*qlI3(s12,s34,0._dp,msq,msq,msq,musq,0)
      return
      end
c       - 1/2*[s12-s34]^-1
c          + 1/2*B0f(p1,mt,mt)*s34*[s12-s34]^-2
c          - 1/2*B0f(p12,mt,mt)*s34*[s12-s34]^-2
c          - C0DDHK(p1,p2,mt,mt,mt)*mt^2*[s12-s34]^-1


      function qlC0DDHK(s12,s34,msq)
      implicit none
      include 'types.f'
      complex(dp):: qlC0DDHK

      include 'scale.f'
      complex(dp):: qlI3
      real(dp):: s12,s34,msq

      qlC0DDHK=qlI3(s12,s34,0._dp,msq,msq,msq,musq,0)

      return
      end
