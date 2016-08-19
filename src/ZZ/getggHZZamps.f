c--- 
c--- MODIFICATION OF THE ORIGINAL MCFM SUBROUTINE TO ALLOW FOR ANOMALOUS H-Z-Z COUPLINGS  (e.g. nproc=128)
c--- SAME CHOICE OF CONVENTIONS AS IN JHUGEN
c--- 
      subroutine getggHZZamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      4d0*esq*gsq/(16d0*pisq)*esq * delta(a,b)
c---
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      integer h1,h34,h56
      double precision p(mxpart,4),mb2,mt2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmt(2,2),ggHmb(2,2),qlI3,C0mt,C0mb,prop12,prop34,prop56,
     & H4l(2,2),sinthw,higgsprop
      double precision rescale 
      double complex anomhzzamp

!==== for width studies rescale by appropriate factor 
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
      else
         rescale=1d0
      endif

      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
      if(hmass.lt.zip) then
         return
      endif
      call spinoru(6,p,za,zb)

c--- squared masses and sin(thetaw)     
      mt2=mt**2
      mb2=mb**2
      sinthw=dsqrt(xw)
      
      
c--- propagator factors
      prop12=higgsprop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Amplitudes for production 
      C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)
      C0mb=qlI3(zip,zip,s(1,2),mb2,mb2,mb2,musq,0)
   
c------ top quark in the loop
      ggHmt(2,2)=mt2*(2d0-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
     & /(2d0*wmass*sinthw)
      ggHmt(1,1)=ggHmt(2,2)*za(1,2)/zb(1,2)
      ggHmt(2,2)=ggHmt(2,2)*zb(1,2)/za(1,2)

c------ bottom quark in the loop
      ggHmb(2,2)=mb2*(2d0-s(1,2)*C0mb*(1d0-4d0*mb2/s(1,2)))
     & /(2d0*wmass*sinthw)
      ggHmb(1,1)=ggHmb(2,2)*za(1,2)/zb(1,2)
      ggHmb(2,2)=ggHmb(2,2)*zb(1,2)/za(1,2)


      H4l(1,1)=anomhzzamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,1)=anomhzzamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(1,2)=anomhzzamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,2)=anomhzzamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56


c--- Assemble: insert factor of (im) here      
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=im*ggHmb(h1,h1)*H4l(h34,h56)*prop12
      Mloop_tquark(h1,h1,h34,h56)=im*ggHmt(h1,h1)*H4l(h34,h56)*prop12
      enddo
      enddo
      enddo


c--- Rescale for width study
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=rescale*Mloop_bquark(h1,h1,h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=rescale*Mloop_tquark(h1,h1,h34,h56)
      enddo
      enddo
      enddo
      

      return
      end
      

      subroutine getggH2ZZamps(p,Mloop_bquark,Mloop_tquark)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'anom_higgs.f' 
      include 'spinzerohiggs_anomcoupl.f'
      integer h1,h34,h56
      double precision p(mxpart,4),mb2,mt2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmt(2,2),ggHmb(2,2),qlI3,C0mt,C0mb,prop12,prop34,prop56,
     & H4l(2,2),sinthw,higgs2prop
      double precision rescale 
      double complex anomhzzamp

!==== for width studies rescale by appropriate factor 
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
      else
         rescale=1d0
      endif

      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
      if(h2mass.lt.zip) then
         return
      endif
     
      call spinoru(6,p,za,zb)

c--- squared masses and sin(thetaw)     
      mt2=mt**2
      mb2=mb**2
      sinthw=dsqrt(xw)
      
      
c--- propagator factors
      prop12=higgs2prop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Amplitudes for production 
      C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)
      C0mb=qlI3(zip,zip,s(1,2),mb2,mb2,mb2,musq,0)
   
c------ top quark in the loop
      ggHmt(2,2)=mt2*(2d0-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
     & /(2d0*wmass*sinthw)
      ggHmt(1,1)=ggHmt(2,2)*za(1,2)/zb(1,2)
      ggHmt(2,2)=ggHmt(2,2)*zb(1,2)/za(1,2)

c------ bottom quark in the loop
      ggHmb(2,2)=mb2*(2d0-s(1,2)*C0mb*(1d0-4d0*mb2/s(1,2)))
     & /(2d0*wmass*sinthw)
      ggHmb(1,1)=ggHmb(2,2)*za(1,2)/zb(1,2)
      ggHmb(2,2)=ggHmb(2,2)*zb(1,2)/za(1,2)


      H4l(1,1)=anomhzzamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,1)=anomhzzamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(1,2)=anomhzzamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,2)=anomhzzamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56


c--- Assemble: insert factor of (im) here      
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=im*ggHmb(h1,h1)*H4l(h34,h56)*prop12
      Mloop_tquark(h1,h1,h34,h56)=im*ggHmt(h1,h1)*H4l(h34,h56)*prop12
      enddo
      enddo
      enddo


c--- Rescale for width study
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=rescale*Mloop_bquark(h1,h1,h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=rescale*Mloop_tquark(h1,h1,h34,h56)
      enddo
      enddo
      enddo
      

      return
      end

