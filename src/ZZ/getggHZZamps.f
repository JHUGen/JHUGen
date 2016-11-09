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
      double precision p(mxpart,4),mb2,mt2,mtX2,mbX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),facHiggs,higgsprop
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
      ggHmq(:,:,:)=czip
      H4l(:,:)=czip

      call spinoru(6,p,za,zb)

c--- propagator factors
      prop12=higgsprop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Factor
      facHiggs=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,1,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)

c--- Amplitudes for decay
      H4l(1,1)=anomhzzamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2
      H4l(2,1)=anomhzzamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2
      H4l(1,2)=anomhzzamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2
      H4l(2,2)=anomhzzamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2
      H4l(:,:) = H4l(:,:)*facHiggs

c--- Assemble
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=ggHmq(1,h1,h1)*H4l(h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=ggHmq(2,h1,h1)*H4l(h34,h56)
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
      double precision p(mxpart,4),mb2,mt2,mbX2,mtX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),facHiggs,higgs2prop
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
      ggHmq(:,:,:)=czip
      H4l(:,:)=czip

      call spinoru(6,p,za,zb)

c--- propagator factors
      prop12=higgs2prop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Factor
      facHiggs=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,2,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)

c--- Amplitudes for decay
      H4l(1,1)=anomhzzamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2
      H4l(2,1)=anomhzzamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2
      H4l(1,2)=anomhzzamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2
      H4l(2,2)=anomhzzamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2
      H4l(:,:) = H4l(:,:)*facHiggs

c--- Assemble
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=ggHmq(1,h1,h1)*H4l(h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=ggHmq(2,h1,h1)*H4l(h34,h56)
      enddo
      enddo
      enddo

      return
      end

