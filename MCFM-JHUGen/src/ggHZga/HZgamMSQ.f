      function HZgamMSQ(p3,p4,p5)
      implicit none
      include 'types.f'
      real(dp):: HZgamMSQ
      
C-----Author: R.K.Ellis July 2012
C-----Matrix element squared for the Higgs going to l^-(p3) l^+(p5) Gamma(p5)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'couple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      include 'msbarmasses.f'
      include 'kpart.f'
      include 'first.f'
      real(dp):: s12,s34,s35,s45
      real(dp):: v1(2),vfbar,Qf,mtsq,mwsq,cotw,fac,
     & mt_eff,massfrun
      complex(dp):: F(2),FfDDHK,fWDDHK,prop34
      integer:: j,top,p3,p4,p5
      parameter(top=2)
      save mt_eff
                        
      if (first) then            
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif  
        first=.false.
      endif
      
      fac=esq**4/xw/(sixteen*pisq*wmass)**2
      mwsq=wmass**2
      mtsq=mt**2
      v1(1)=l1
      v1(2)=r1
      Qf=Q(top)
      cotw=sqrt((one-xw)/xw)
      vfbar=half*(L(top)+R(top))
      s34=s(p3,p4)
      s35=s(p3,p5)
      s45=s(p4,p5)
      s12=s34+s35+s45
      prop34=s34/cplx2(s34-zmass**2,zmass*zwidth)
C  -- j labels the helicity of the lepton line
      do j=1,2
      F(j)=fWDDHK(s12,s34,mwsq)
     & *(cplx1(q1)+cotw*v1(j)*prop34)
     & +ffDDHK(s12,s34,mtsq)
     & *four*Qf*xn*(mt_eff**2/s34)
     & *(cplx1(Qf*q1)+vfbar*v1(j)*prop34)
      enddo

      HZgamMSQ=half*fac*s34*(s35**2+s45**2)
     & *(abs(F(1))**2+abs(F(2))**2)          

      return
      end 
