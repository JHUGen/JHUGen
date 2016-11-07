      subroutine getggHWWamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      (2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2 * delta(a,b)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      include 'spinzerohiggs_anomcoupl.f'
      integer h1
      double precision p(mxpart,4)
      double precision rescale
      double complex Mloop_tquark(2,2),Mloop_bquark(2,2),
     & fachiggs,ggHmq(2,2,2),e3De4,higgsprop,props
      double complex anomhwwamp

!==== for width studies rescale by appropriate factor
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
      else
         rescale=1d0
      endif

      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(hmass.lt.zip) then
         return
      endif
      ggHmq(:,:,:)=czip

      call spinoru(6,p,za,zb)

c--- Production amplitude
      ! Overall factor=1
      call anomhggvtxamp(1,2,1,za,zb,ggHmq)

c--- Decay amplitude
      e3De4=-anomhwwamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)
      props=
     &  cone/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     & *cone/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
      fachiggs=higgsprop(s(1,2))*e3De4*im*props*rescale

      print *,"props=",props*(s(3,4)*s(5,6))
      print *,"e3De4=",e3De4/(s(3,4)*s(5,6))
      print *,"fachiggs=",higgsprop(s(1,2))
      print *,"s(1,2)=",s(1,2)

      print *,"amphiggs_bot=",ggHmq(1,1,1)*zb(1,2)/za(1,2)*im*e3De4
     & /(s(3,4)*s(5,6))
      print *,"amphiggs_top=",ggHmq(2,1,1)*zb(1,2)/za(1,2)*im*e3De4
     & /(s(3,4)*s(5,6))


      do h1=1,2
         Mloop_bquark(h1,h1)=fachiggs*ggHmq(1,h1,h1)
         Mloop_tquark(h1,h1)=fachiggs*ggHmq(2,h1,h1)
      enddo

      print *,"Mloop_bquark(1,1)=",Mloop_bquark(1,1)
      print *,"Mloop_bquark(2,2)=",Mloop_bquark(2,2)
      print *,"Mloop_tquark(1,1)=",Mloop_tquark(1,1)
      print *,"Mloop_tquark(2,2)=",Mloop_tquark(2,2)

      return
      end



      subroutine getggH2WWamps(p,Mloop_bquark,Mloop_tquark)
c--- Returns a series of arrays representing the dressed amp[itudes
c--- for the process gg->Higgs->ZZ; there are:
c---        Mloop_bquark(h1,h2,h34,h56)   top quark mass=mt
c---        Mloop_tquark(h1,h2,h34,h56)   bottom quark mass=mb
c---
c--- The overall factor on the amplitude is:
c---
c---      (2d0*gwsq*gsq/(16d0*pisq)*gwsq/2d0)**2 * delta(a,b)
c---
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'anom_higgs.f'
      include 'spinzerohiggs_anomcoupl.f'
      integer h1
      double precision p(mxpart,4)
      double precision rescale
      double complex Mloop_tquark(2,2),Mloop_bquark(2,2),
     & fachiggs,ggHmq(2,2,2),e3De4,higgs2prop,props
      double complex anomhwwamp

!==== for width studies rescale by appropriate factor
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
      else
         rescale=1d0
      endif

      Mloop_tquark(:,:)=czip
      Mloop_bquark(:,:)=czip
      if(h2mass.lt.zip) then
         return
      endif
      ggHmq(:,:,:)=czip

      call spinoru(6,p,za,zb)

c--- Production amplitude
      ! Overall factor=1
      call anomhggvtxamp(1,2,2,za,zb,ggHmq)

c--- Decay amplitude
      e3De4=-anomhwwamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)
      props=
     &  cone/dcmplx(s(3,4)-wmass**2,wmass*wwidth)
     & *cone/dcmplx(s(5,6)-wmass**2,wmass*wwidth)
      fachiggs=higgs2prop(s(1,2))*e3De4*im*props*rescale

      do h1=1,2
         Mloop_bquark(h1,h1)=fachiggs*ggHmq(1,h1,h1)
         Mloop_tquark(h1,h1)=fachiggs*ggHmq(2,h1,h1)
      enddo

      return
      end
