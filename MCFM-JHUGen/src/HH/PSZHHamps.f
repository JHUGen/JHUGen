      subroutine PSZHHamps(pa,pb,pc,pd,amp)
C--    Formula taken from Plehn, Spira and Zerwas                        
C--    Nucl. Phys. B479 (1996) 46
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'scale.f'
      include 'first.f'
      double precision pa(4),pb(4),pc(4),pd(4),pcsq,pdsq,
     & ss,tt,uu,S,T,U,T1,T2,U1,U2,rhoc,rhod,tauQ,mQsq,lambdaHHH
      double complex Ftriangle,Fbox,Gbox,
     & Dabc,Dbac,Dacb,Cab,Cbc,Cac,Ccd,Cbd,Cad,qlI4,qlI3,
     & Cdelta,Cbox,amp(2)

c--- initialize QCDLoop, if necessary
      if (first) then
        call qlinit
        first=.false.
      endif


      mQsq=mt**2
      ss=(pa(4)+pb(4))**2
     &  -(pa(1)+pb(1))**2-(pa(2)+pb(2))**2-(pa(3)+pb(3))**2
      tt=(pa(4)+pc(4))**2
     &  -(pa(1)+pc(1))**2-(pa(2)+pc(2))**2-(pa(3)+pc(3))**2
      uu=(pb(4)+pc(4))**2
     &  -(pb(1)+pc(1))**2-(pb(2)+pc(2))**2-(pb(3)+pc(3))**2
      S=ss/mQsq
      T=tt/mQsq
      U=uu/mQsq
      lambdaHHH=3d0*(hmass/zmass)**2
      rhoc=hmass**2/mQsq
      rhod=hmass**2/mQsq
      T1=T-rhoc
      U1=U-rhoc
      T2=T-rhod
      U2=U-rhod
      tauQ=4d0/S
C      ss=(pa+pb)**2
C      tt=(pa+pc)**2
C      uu=(pb+pc)**2
      pcsq=pc(4)**2-pc(1)**2-pc(2)**2-pc(3)**2
      pdsq=pd(4)**2-pd(1)**2-pd(2)**2-pd(3)**2
      Cab=qlI3(0d0,0d0,ss,mQsq,mQsq,mQsq,musq,0)
      Cbc=qlI3(0d0,pcsq,uu,mQsq,mQsq,mQsq,musq,0)
      Cac=qlI3(0d0,pcsq,tt,mQsq,mQsq,mQsq,musq,0)
      Ccd=qlI3(pcsq,pdsq,ss,mQsq,mQsq,mQsq,musq,0)
      Cbd=qlI3(0d0,pdsq,tt,mQsq,mQsq,mQsq,musq,0)
      Cad=qlI3(0d0,pdsq,uu,mQsq,mQsq,mQsq,musq,0)
      Dabc=qlI4(0d0,0d0,pcsq,pdsq,ss,uu,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dbac=qlI4(0d0,0d0,pcsq,pdsq,ss,tt,mQsq,mQsq,mQsq,mQsq,musq,0)
      Dacb=qlI4(0d0,pcsq,0d0,pdsq,tt,uu,mQsq,mQsq,mQsq,mQsq,musq,0)


      Ftriangle=2d0/S*(2d0+(4d0-S)*mQsq*Cab)

      Fbox=(4d0*S+8d0*S*mQsq*Cab-2d0*S*(S+rhoc+rhod-8d0)
     & *mQsq**2*(Dabc+Dbac+Dacb)
     & +(rhoc+rhod-8d0)*mQsq*(T1*Cac+U1*Cbc+U2*Cad+T2*Cbd
     & -(T*U-rhoc*rhod)*mQsq*Dacb))/S**2
      Gbox=((T**2+rhoc*rhod-8d0*T)*mQsq
     & *(S*Cab+T1*Cac+T2*Cbd-S*T*mQsq*Dbac)
     & +(U**2+rhoc*rhod-8d0*U)*mQsq
     & *(S*Cab+U1*Cbc+U2*Cad-S*U*mQsq*Dabc)
     & -(T**2+U**2-2d0*rhoc*rhod)*(T+U-8d0)*mQsq*Ccd
     & -2d0*(T+U-8d0)*(T*U-rhoc*rhod)*mQsq**2
     & *(Dabc+Dbac+Dacb))/(S*(T*U-rhoc*rhod))

      Cdelta=lambdaHHH*zmass**2/dcmplx(ss-hmass**2,hmass*hwidth) 
      Cbox=cone

      amp(1)=Cdelta*Ftriangle+Cbox*Fbox
      amp(2)=Cbox*Gbox
      return
      end
