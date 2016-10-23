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
     & ggHmt(2,2),ggHmb(2,2),qlI3,C0mt,C0mb,prop12,prop34,prop56,
     & H4l(2,2),sinthw,higgsprop,B0a,B0b,qlI2
      double precision rescale
      double complex anomhzzamp
      double complex a1,a3,a1_4gen,a3_4gen,C0mtX,C0mbX

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
      mtX2=mt_4gen**2
      mbX2=mb_4gen**2
      sinthw=dsqrt(xw)


c--- propagator factors
      prop12=higgsprop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Amplitudes for production
      C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)
      C0mb=qlI3(zip,zip,s(1,2),mb2,mb2,mb2,musq,0)
      C0mtX=qlI3(zip,zip,s(1,2),mtX2,mtX2,mtX2,musq,0)
      C0mbX=qlI3(zip,zip,s(1,2),mbX2,mbX2,mbX2,musq,0)


c    Couplings for point-like interactions
      a1 = ghg2+ghg3*s(1,2)/4d0/LambdaBSM**2
      a3 = -2d0*ghg4
      a1_4gen= ghg2_4gen+ghg3_4gen*s(1,2)/4d0/LambdaBSM**2
      a3_4gen= -2d0*ghg4_4gen


c------ top-flavor quarks
c        SM top
c        kappa couplings (quark loop interaction)
         ggHmt(2,2)=mt2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mt*kappa_top*(1d0-4d0*mt2/s(1,2))
     &   +2d0*kappa_top-kappa_tilde_top )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)

c        4th generation top
c        kappa couplings (quark loop interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              mtX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mtX*kappa_4gen_top*(1d0-4d0*mtX2/s(1,2))
     &   +2d0*kappa_4gen_top-kappa_tilde_4gen_top )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)


c        same as above for other helicity
         ggHmt(1,1)=mt2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mt*kappa_top*(1d0-4d0*mt2/s(1,2))
     &    +2d0*kappa_top+kappa_tilde_top )*za(1,2)/zb(1,2)

         ggHmt(1,1)=ggHmt(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

         ggHmt(1,1)=ggHmt(1,1) +
     &              mtX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mtX*kappa_4gen_top*(1d0-4d0*mtX2/s(1,2))
     &    +2d0*kappa_4gen_top+kappa_tilde_4gen_top )*za(1,2)/zb(1,2)

         ggHmt(1,1)=ggHmt(1,1) +
     &             s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)






c------ bot-flavor quarks
c        SM bot
c        kappa couplings (quark loop interaction)
         ggHmb(2,2)=mb2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mb*kappa_bot*(1d0-4d0*mb2/s(1,2))
     &   +2d0*kappa_bot-kappa_tilde_bot )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)

c        4th generation bot
c        kappa couplings (quark loop interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              mbX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mbX*kappa_4gen_bot*(1d0-4d0*mbX2/s(1,2))
     &   +2d0*kappa_4gen_bot-kappa_tilde_4gen_bot )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)



c        same as above for other helicity
         ggHmb(1,1)=mb2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mb*kappa_bot*(1d0-4d0*mb2/s(1,2))
     &    +2d0*kappa_bot+kappa_tilde_bot )*za(1,2)/zb(1,2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              mbX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mbX*kappa_4gen_bot*(1d0-4d0*mbX2/s(1,2))
     &    +2d0*kappa_4gen_bot+kappa_tilde_4gen_bot )*za(1,2)/zb(1,2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)










! c------ top quark in the loop
!       ggHmt(2,2)=mt2*(2d0-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
!      & /(2d0*wmass*sinthw)

c------ MARKUS new code: the approximation for the above ggHmt in the limit mt-->infinity is:
!       ggHmt(2,2)=s(1,2)/3d0/(2d0*wmass*sinthw)


! c------ MARKUS comment: below this is the multiplication with -ep_glu1(mu)*ep_glu2(nu) * g_{mu,nu}
!       ggHmt(1,1)=ggHmt(2,2)*za(1,2)/zb(1,2)
!       ggHmt(2,2)=ggHmt(2,2)*zb(1,2)/za(1,2)





! c_________________________________________________________________
!       print *, ""
!       print *, dsqrt(dabs(s(1,2)))
!       print *, "orig ggHmt(2,2)",ggHmt(2,2)
!       print *, ""
!
!       mt2=(mt*1000d0)**2
! c--- Amplitudes for production
!       C0mt=qlI3(zip,zip,s(1,2),mt2,mt2,mt2,musq,0)
!       C0mb=qlI3(zip,zip,s(1,2),mb2,mb2,mb2,musq,0)
!
! c------ top quark in the loop
!       ggHmt(2,2)=mt2*(2d0-s(1,2)*C0mt*(1d0-4d0*mt2/s(1,2)))
!      & /(2d0*wmass*sinthw)
!
!       print *, "resc ggHmt(2,2)",ggHmt(2,2)
!       print *, "resc approx    ",s(1,2)/(6d0*wmass*sinthw)
!       print *, "ratio    ",ggHmt(2,2)
!      &          /(s(1,2)/(6d0*wmass*sinthw))
!          pause
! c_________________________________________________________________


! c------ bottom quark in the loop
!       ggHmb(2,2)=mb2*(2d0-s(1,2)*C0mb*(1d0-4d0*mb2/s(1,2)))
!      & /(2d0*wmass*sinthw)
!       ggHmb(1,1)=ggHmb(2,2)*za(1,2)/zb(1,2)
!       ggHmb(2,2)=ggHmb(2,2)*zb(1,2)/za(1,2)


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
      double precision p(mxpart,4),mb2,mt2,mbX2,mtX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmt(2,2),ggHmb(2,2),qlI3,C0mt,C0mb,prop12,prop34,prop56,
     & H4l(2,2),sinthw,higgs2prop
      double precision rescale
      double complex anomhzzamp
      double complex a1,a3,a1_4gen,a3_4gen,C0mtX,C0mbX

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
      mtX2=mt_4gen**2
      mbX2=mb_4gen**2
      sinthw=dsqrt(xw)


c--- propagator factors
      prop12=higgs2prop(s(1,2))
      prop34=cone/dcmplx(s(3,4)-zmass**2,zmass*zwidth)
      prop56=cone/dcmplx(s(5,6)-zmass**2,zmass*zwidth)

c--- Amplitudes for production
      C0mt =qlI3(zip,zip,s(1,2),mt2, mt2, mt2, musq,0)
      C0mb =qlI3(zip,zip,s(1,2),mb2, mb2, mb2, musq,0)
      C0mtX=qlI3(zip,zip,s(1,2),mtX2,mtX2,mtX2,musq,0)
      C0mbX=qlI3(zip,zip,s(1,2),mbX2,mbX2,mbX2,musq,0)


c    Couplings for point-like interactions
      a1 = gh2g2+gh2g3*s(1,2)/4d0/Lambda2BSM**2
      a3 = -2d0*gh2g4
      a1_4gen= gh2g2_4gen+gh2g3_4gen*s(1,2)/4d0/Lambda2BSM**2
      a3_4gen= -2d0*gh2g4_4gen


c------ top-flavor quarks
c        SM top
c        kappa2 couplings (quark loop interaction)
         ggHmt(2,2)=mt2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mt*kappa2_top*(1d0-4d0*mt2/s(1,2))
     &   +2d0*kappa2_top-kappa2_tilde_top )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)

c        4th generation top
c        kappa2 couplings (quark loop interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              mtX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mtX*kappa2_4gen_top*(1d0-4d0*mtX2/s(1,2))
     &   +2d0*kappa2_4gen_top-kappa2_tilde_4gen_top )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmt(2,2)=ggHmt(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)


c        same as above for other helicity
         ggHmt(1,1)=mt2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mt*kappa2_top*(1d0-4d0*mt2/s(1,2))
     &    +2d0*kappa2_top+kappa2_tilde_top )*za(1,2)/zb(1,2)

         ggHmt(1,1)=ggHmt(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

         ggHmt(1,1)=ggHmt(1,1) +
     &              mtX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mtX*kappa2_4gen_top*(1d0-4d0*mtX2/s(1,2))
     &    +2d0*kappa2_4gen_top+kappa2_tilde_4gen_top )*za(1,2)/zb(1,2)

         ggHmt(1,1)=ggHmt(1,1) +
     &             s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)






c------ bot-flavor quarks
c        SM bot
c        kappa2 couplings (quark loop interaction)
         ggHmb(2,2)=mb2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mb*kappa2_bot*(1d0-4d0*mb2/s(1,2))
     &   +2d0*kappa2_bot-kappa2_tilde_bot )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*zb(1,2)/za(1,2)+a3*0.5d0*zb(1,2)**2)

c        4th generation bot
c        kappa2 couplings (quark loop interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              mbX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mbX*kappa2_4gen_bot*(1d0-4d0*mbX2/s(1,2))
     &   +2d0*kappa2_4gen_bot-kappa2_tilde_4gen_bot )*zb(1,2)/za(1,2)

c        ghg couplings (point-like interaction)
         ggHmb(2,2)=ggHmb(2,2) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*zb(1,2)/za(1,2)+a3_4gen*0.5d0*zb(1,2)**2)



c        same as above for other helicity
         ggHmb(1,1)=mb2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mb*kappa2_bot*(1d0-4d0*mb2/s(1,2))
     &    +2d0*kappa2_bot+kappa2_tilde_bot )*za(1,2)/zb(1,2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1*za(1,2)/zb(1,2)-a3*0.5d0*za(1,2)**2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              mbX2/(2d0*wmass*sinthw)*
     &   (-s(1,2)*C0mbX*kappa2_4gen_bot*(1d0-4d0*mbX2/s(1,2))
     &    +2d0*kappa2_4gen_bot+kappa2_tilde_4gen_bot )*za(1,2)/zb(1,2)

         ggHmb(1,1)=ggHmb(1,1) +
     &              s(1,2)/3d0/(2d0*wmass*sinthw)*
     &   (a1_4gen*za(1,2)/zb(1,2)-a3_4gen*0.5d0*za(1,2)**2)






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

