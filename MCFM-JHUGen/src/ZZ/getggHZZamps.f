c---
c--- MODIFICATION OF THE ORIGINAL MCFM SUBROUTINE TO ALLOW FOR ANOMALOUS H-Z-Z COUPLINGS  (e.g. nproc=128)
c--- SAME CHOICE OF CONVENTIONS AS IN JHUGEN
c---
      subroutine getggHZZamps(p,Mloop_bquark,Mloop_tquark,
     &Mloop_c6_propagator,Mloop_c6_decay,
     &Mloop_c6_production,Mloop_c6_width)
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
      include 'plabel.f'
      include 'AnomZffCouplings.f'
      include 'qlfirst.f'
      integer h1,h34,h56
      integer AllowAnomTriLinear, wfn
      double precision wf(11)
      double precision p(mxpart,4),mb2,mt2,mtX2,mbX2
      double precision MH2
      double precision dB0h, dZh
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & Mloop_c6_production(2,2,2,2),Mloop_c6_width(2,2,2,2),
     & Mloop_c6_propagator(2,2,2,2),Mloop_c6_decay(2,2,2,2),
     & ggHmt(2,2),ggHmt_c6(2,2),ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),H4lSM(2,2),facHZZ,facHZA,facHAZ,facHAA,higgsprop,
     & prop12_c6,sigmah,hwidth_c6,width_c6,
     & H4l_c6_gmunu(2,2), H4l_c6_qmuqnu(2,2)
      double precision rescale
      double complex anomhzzamp,anomhzaamp,anomhaaamp

!==== for width studies rescale by appropriate factor
      if((keep_smhiggs_norm).and.(anom_higgs)) then
         rescale=chi_higgs**2
      else
         rescale=1d0
      endif

c---  end width corrections
      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
c---  Set c6 corrections to 0 by default
      AllowAnomTriLinear = zip
      prop12_c6=zip
      width_c6=zip
      Mloop_c6_production(:,:,:,:)=czip
      Mloop_c6_width(:,:,:,:)=czip
      Mloop_c6_propagator(:,:,:,:)=czip
      Mloop_c6_decay(:,:,:,:)=czip
      ggHmt_c6(:,:)=czip
      ggHmt(:,:)=czip
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
      facHZZ=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))
      facHZA=-im*rescale*prop12*prop34/s(5,6)/(2d0*xw*(1d0-xw))
      facHAZ=-im*rescale*prop12/s(3,4)*prop56/(2d0*xw*(1d0-xw))
      facHAA=im*rescale*prop12/s(3,4)/s(5,6)/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,1,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)


c--- Check if any wavefunctions are non zero
      wf = (/ t1_c6,t2_c6,t3_c6,t4_c6,t5_c6,t6_c6,
     & w1_c6,w2_c6,w3_c6,w4_c6,w5_c6 /)
      do wfn = 1, size(wf)
        if(wf(wfn) .ne. 0) then 
          AllowAnomTriLinear = 1
        endif
      end do
c--- Do self couplings if allowed  
      if(AllowAnomTriLinear .eq. 1)then
c--- c6 correction to propagator      
        prop12_c6=-higgsprop(s(1,2))*(sigmah(s(1,2),c6,w1_c6))
                  
c--- for width corrections due to c6 operator
        MH2 = hmass**2
        dB0h = (-9 + 2*Sqrt(3.)*Pi)/(9.*MH2)
        dZh = (-9*c6*(2.d0 + c6)*dB0h*MH2**2)/(32.d0*Pi**2*vevsq)
        hwidth_c6 = 0.0023*c6*hwidth
        width_c6 = im*hmass*(t5_c6*w4_c6*dZh*hwidth -
     &  t6_c6*(w5_c6*dZh/2.d0*hwidth + hwidth_c6))*prop12

c--- c6 production corrections
        call anomhggvtxamp_c6(za,zb,ggHmt_c6)
c--- c6 decay correction
        call anomhzzamp_c6(prop34,prop56,za,zb,
     & H4l_c6_gmunu,H4l_c6_qmuqnu)

c---  Load SM gghtloop and for c6 corrections
        call SMggHmtvertex(za,zb,ggHmt)

c---  Load SM H4l amplitude for c6 corrections       
        H4lSM(1,1)=za(3,5)*zb(4,6)*l1*l2*facHZZ
        H4lSM(2,1)=za(4,5)*zb(3,6)*r1*l2*facHZZ
        H4lSM(1,2)=za(3,6)*zb(4,5)*l1*r2*facHZZ
        H4lSM(2,2)=za(4,6)*zb(3,5)*r1*r2*facHZZ
      endif

c--- Setting Anomalous Zff Couplings 
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(3) .eq. 'el') .or. (plabel(3) .eq. 'ml')
     &.or. (plabel(3) .eq. 'tl')) then
          l1 = leZ
          r1 = reZ 
        elseif (plabel(3) .eq. 'nl') then
          l1 = lnZ*dsqrt(3d0)
          r1 = rnZ*dsqrt(3d0) 
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
          l1=lqdZ*dsqrt(3d0)
          r1=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
          l1=lquZ*dsqrt(3d0)
          r1=rquZ*dsqrt(3d0)
        endif 
      endif
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(5) .eq. 'el') .or. (plabel(5) .eq. 'ml')
     &.or. (plabel(5) .eq. 'tl')) then
          l2 = leZ
          r2 = reZ 
        elseif (plabel(5) .eq. 'nl') then
          l2 = lnZ*dsqrt(3d0)
          r2 = rnZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
            l2=lqdZ*dsqrt(3d0)
            r2=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
            l2=lquZ*dsqrt(3d0)
            r2=rquZ*dsqrt(3d0)  
        endif
      endif

c--- Amplitudes for decay
      H4l(1,1)=
     &  anomhzzamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2*facHZZ
     & +anomhzaamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(5,6,3,4,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(3,4,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,1)=
     &  anomhzzamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2*facHZZ
     & +anomhzaamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(5,6,4,3,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(4,3,5,6,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(1,2)=
     &  anomhzzamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2*facHZZ
     & +anomhzaamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(6,5,3,4,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(3,4,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,2)=
     &  anomhzzamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2*facHZZ
     & +anomhzaamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(6,5,4,3,1,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(4,3,6,5,1,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA

c--- Assemble
      do h1=1,2
      do h34=1,2
      do h56=1,2
      Mloop_bquark(h1,h1,h34,h56)=ggHmq(1,h1,h1)*H4l(h34,h56)
      Mloop_tquark(h1,h1,h34,h56)=ggHmq(2,h1,h1)*H4l(h34,h56)
c--- Assemble c6 corrections
c-------------------------------------------------------------
c--- Could replace ggHmt(h1,h1) with ggHmq(2,h1,h1) to include anomalous effects
c--- Could replace H4lSM with H4l to include anomalous effects 
c--- Propbably not correct to do but may be worth considering
c-------------------------------------------------------------
c--- propagator correction
      Mloop_c6_propagator(h1,h1,h34,h56)=t1_c6*ggHmt(h1,h1)*
     &     H4lSM(h34,h56)*prop12_c6/prop12
c--- production correction
      Mloop_c6_production(h1,h1,h34,h56) = 
     &     t4_c6*ggHmt_c6(h1,h1)* H4lSM(h34,h56)
     &     * (2*wmass*sqrt(xw))
c--- decay correction  
      Mloop_c6_decay(h1,h1,h34,h56)=im*ggHmt(h1,h1)*
     &     (t2_c6*H4l_c6_gmunu(h34,h56)+t3_c6*H4l_c6_qmuqnu(h34,h56))*
     &     prop12*rescale*(1d0/(2*wmass*sqrt(xw)))
c---  width correction
      Mloop_c6_width(h1,h1,h34,h56)=ggHmt(h1,h1)*
     &     H4lSM(h34,h56)*width_c6
      enddo
      enddo
      enddo

c---  Assemble c6 ggH vertex corrections

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
      include 'plabel.f'
      include 'AnomZffCouplings.f'
      integer h1,h34,h56
      double precision p(mxpart,4),mb2,mt2,mbX2,mtX2
      double complex Mloop_bquark(2,2,2,2),Mloop_tquark(2,2,2,2),
     & ggHmq(2,2,2),prop12,prop34,prop56,
     & H4l(2,2),facHZZ,facHZA,facHAZ,facHAA,higgs2prop
      double precision rescale
      double complex anomhzzamp,anomhzaamp,anomhaaamp

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
      facHZZ=im*rescale*prop12*prop34*prop56/(2d0*xw*(1d0-xw))
      facHZA=-im*rescale*prop12*prop34/s(5,6)/(2d0*xw*(1d0-xw))
      facHAZ=-im*rescale*prop12/s(3,4)*prop56/(2d0*xw*(1d0-xw))
      facHAA=im*rescale*prop12/s(3,4)/s(5,6)/(2d0*xw*(1d0-xw))

c--- Amplitudes for production
      call anomhggvtxamp(1,2,2,za,zb,ggHmq)
      ! Overall factor=1
      !ggHmq(:,:,:) = ggHmq(:,:,:)

c--- Setting Anomalous Zff Couplings 
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(3) .eq. 'el') .or. (plabel(3) .eq. 'ml')
     &.or. (plabel(3) .eq. 'tl')) then
          l1 = leZ
          r1 = reZ 
        elseif (plabel(3) .eq. 'nl') then
          l1 = lnZ*dsqrt(3d0)
          r1 = rnZ*dsqrt(3d0) 
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
          l1=lqdZ*dsqrt(3d0)
          r1=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
          l1=lquZ*dsqrt(3d0)
          r1=rquZ*dsqrt(3d0)
        endif 
      endif
      if (AllowAnomalousZffCouplings .eq. 1) then
        if ((plabel(5) .eq. 'el') .or. (plabel(5) .eq. 'ml')
     &.or. (plabel(5) .eq. 'tl')) then
          l2 = leZ
          r2 = reZ 
        elseif (plabel(5) .eq. 'nl') then
          l2 = lnZ*dsqrt(3d0)
          r2 = rnZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'bq') .or. (plabel(5) .eq. 'sq')
     &.or. (plabel(5) .eq. 'dq')) then
            l2=lqdZ*dsqrt(3d0)
            r2=rqdZ*dsqrt(3d0)
        elseif ((plabel(5) .eq. 'uq') .or. (plabel(5) .eq. 'cq')) then
            l2=lquZ*dsqrt(3d0)
            r2=rquZ*dsqrt(3d0)  
        endif
      endif

c--- Amplitudes for decay
      H4l(1,1)=
     &  anomhzzamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*l2*facHZZ
     & +anomhzaamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(5,6,3,4,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(3,4,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,1)=
     &  anomhzzamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*l2*facHZZ
     & +anomhzaamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(5,6,4,3,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*l2*facHAZ
     & +anomhaaamp(4,3,5,6,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(1,2)=
     &  anomhzzamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*r2*facHZZ
     & +anomhzaamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*l1*q2*facHZA
     & +anomhzaamp(6,5,3,4,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(3,4,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA
      H4l(2,2)=
     &  anomhzzamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*r2*facHZZ
     & +anomhzaamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*r1*q2*facHZA
     & +anomhzaamp(6,5,4,3,2,s(1,2),s(5,6),s(3,4),za,zb)*q1*r2*facHAZ
     & +anomhaaamp(4,3,6,5,2,s(1,2),s(3,4),s(5,6),za,zb)*q1*q2*facHAA

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

