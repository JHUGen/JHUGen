      double complex function ffDDHK(s12,s34,msq)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The top loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
      double precision s12,s34,msq
      double complex C0DDHK,C2DDHK
      ffDDHK=C0DDHK(s12,s34,msq)+4d0*C2DDHK(s12,s34,msq)
      return
      end

      double complex function fWDDHK(s12,s34,msq)
      implicit none
C-----Author: R.K.Ellis July 2012
C-----The W-loop function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(4)
c----- msq here should only be equal to wmass**2
      double precision s12,s34,msq
      double complex C0DDHK,C2DDHK
      fWDDHK=2d0*(s12/msq*(1d0-2d0*msq/s34)
     & +2d0*(1d0-6d0*msq/s34))*C2DDHK(s12,s34,msq)
     & +4d0*(1d0-4d0*msq/s34)*C0DDHK(s12,s34,msq)
      return
      end

      double complex function f0DDHK(s12,s34)
      implicit none
C-----The F0-function taken from
C-----Djouadi,Driesen,Hollik,Kraft arXiv:hep-ph/9701342v1, Eq(2)
C-----suitably generalized to allow off-shell Z-line
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      include 'couple.f'
      include 'part.f'
      include 'first.f'
      include 'msbarmasses.f'
      integer top
      double complex ffDDHK,fWDDHK
      double precision s12,s34,cotw,mtsq,mwsq,mt_eff,massfrun
      parameter(top=2)
      save mt_eff
!$omp threadprivate(mt_eff)
                       
      if (first) then           
c--- run mt to appropriate scale
        if (part .eq. 'lord') then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif  
        first=.false.
      endif
      
      cotw=sqrt((1d0-xw)/xw)
      mtsq=mt**2 
      mwsq=wmass**2 
      f0DDHK=s34*(cotw*fWDDHK(s12,s34,mwsq)
     &+2d0*Q(top)*xn*mt_eff**2/s34*(L(top)+R(top))*ffDDHK(s12,s34,mtsq))
      return
      end
      

      double complex function qlC2DDHK(s12,s34,msq)
      implicit none
      include 'scale.f'
      double complex qlI2,qlI3
      double precision s12,s34,msq 
      qlC2DDHK=-dcmplx(0.5d0/(s12-s34))
     & +0.5d0*s34/(s12-s34)**2
     & *(qlI2(s34,msq,msq,musq,0)-qlI2(s12,msq,msq,musq,0))
     & -msq/(s12-s34)*qlI3(s12,s34,0d0,msq,msq,msq,musq,0)
      return
      end
c       - 1/2*[s12-s34]^-1
c          + 1/2*B0f(p1,mt,mt)*s34*[s12-s34]^-2
c          - 1/2*B0f(p12,mt,mt)*s34*[s12-s34]^-2
c          - C0DDHK(p1,p2,mt,mt,mt)*mt^2*[s12-s34]^-1


      double complex function qlC0DDHK(s12,s34,msq)
      implicit none
      include 'scale.f'
      double complex qlI3
      double precision s12,s34,msq 
      
      qlC0DDHK=qlI3(s12,s34,0d0,msq,msq,msq,musq,0)
      
      return
      end
