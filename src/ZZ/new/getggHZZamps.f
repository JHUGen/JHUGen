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
      double complex FFa1(2,2), FFa2(2,2), FFa3(2,2)
      double complex aa1,aa2,aa3
      double precision shat,q3_q3,q4_q4
      double complex ghz1_dyn,ghz2_dyn,ghz3_dyn,ghz4_dyn

!==== for width studies rescale by appropriate factor 
      if((keep_smhiggs_norm).and.(anom_higgs)) then 
         rescale=chi_higgs**2 
      else
         rescale=1d0
      endif

      Mloop_bquark(:,:,:,:)=czip
      Mloop_tquark(:,:,:,:)=czip
     
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


      IF( AllowAnomalousCouplings ) THEN
c----- MARKUS: ANOMALOUS COUPLINGS (factor out mV**2 wrt. JHUGen norm)

      shat  = s(1,2)
      q3_q3 = s(3,4)
      q4_q4 = s(5,6)

c------ MARKUS: NEW FORM FACTORS FOR ANOMALOUS COUPLINGS
c L1L2
      FFa1(1,1) = za(3,5)*zb(4,6)*shat
      FFa2(1,1) =-0.5d0*za(3,5)**2*zb(3,6)*zb(4,5) 
     &           -0.5d0*za(3,5)*za(3,6)*zb(3,6)*zb(4,6) 
     &           -0.5d0*za(3,5)*za(4,5)*zb(4,5)*zb(4,6) 
     &           -0.5d0*za(3,6)*za(4,5)*zb(4,6)**2
      FFa3(1,1) = 0.5d0*za(3,4)*za(5,6)*zb(4,6)**2 
     &           -0.5d0*za(3,5)**2*zb(3,4)*zb(5,6)


c R1L2:
      FFa1(2,1) = za(4,5)*zb(3,6)*shat
      FFa2(2,1) =-0.5d0*za(3,5)*za(4,5)*zb(3,5)*zb(3,6) 
     &           -0.5d0*za(3,5)*za(4,6)*zb(3,6)**2 
     &           -0.5d0*za(4,5)**2*zb(3,5)*zb(4,6) 
     &           -0.5d0*za(4,5)*za(4,6)*zb(3,6)*zb(4,6)
      FFa3(2,1) =-0.5d0*za(3,4)*za(5,6)*zb(3,6)**2 
     &           +0.5d0*za(4,5)**2*zb(3,4)*zb(5,6)


c L1R2:
      FFa1(1,2) = za(3,6)*zb(4,5)*shat
      FFa2(1,2) =-0.5d0*za(3,5)*za(3,6)*zb(3,5)*zb(4,5) 
     &           -0.5d0*za(3,5)*za(4,6)*zb(4,5)**2 
     &           -0.5d0*za(3,6)**2*zb(3,5)*zb(4,6) 
     &           -0.5d0*za(3,6)*za(4,6)*zb(4,5)*zb(4,6)
      FFa3(1,2) =-0.5d0*za(3,4)*za(5,6)*zb(4,5)**2 
     &           +0.5d0*za(3,6)**2*zb(3,4)*zb(5,6)


c R1R2:
      FFa1(2,2) = za(4,6)*zb(3,5)*shat
      FFa2(2,2) =-0.5d0*za(3,6)*za(4,5)*zb(3,5)**2 
     &           -0.5d0*za(3,6)*za(4,6)*zb(3,5)*zb(3,6) 
     &           -0.5d0*za(4,5)*za(4,6)*zb(3,5)*zb(4,5)
     &           -0.5d0*za(4,6)**2*zb(3,6)*zb(4,5)
      FFa3(2,2) = 0.5d0*za(3,4)*za(5,6)*zb(3,5)**2
     &           -0.5d0*za(4,6)**2*zb(3,4)*zb(5,6)

      FFa3 = FFa3 * (0d0,-1d0)!  phase convention to match JHUGen

c--- MARKUS: define q^2 dependent couplings
      ghz1_dyn = ghz1 
     &         + ghz1_prime * Lambda_z1**4/( Lambda_z1**2 
     &         + abs(q3_q3) )/(Lambda_z1**2 + abs(q4_q4))
     &         + ghz1_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z1**2                                   
     &         + ghz1_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z1**2                                                                                   
     &         + ghz1_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz1_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z1**4
     &         + ghz1_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z1**4
     &         + ghz1_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z1**4



      ghz2_dyn = ghz2 
     &         + ghz2_prime * Lambda_z2**4/( Lambda_z2**2 
     &         + abs(q3_q3) )/(Lambda_z2**2 + abs(q4_q4))
     &         + ghz2_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z2**2                                   
     &         + ghz2_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z2**2                                                                                   
     &         + ghz2_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz2_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z2**4
     &         + ghz2_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z2**4
     &         + ghz2_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z2**4


      ghz3_dyn = ghz3 
     &         + ghz3_prime * Lambda_z3**4/( Lambda_z3**2 
     &         + abs(q3_q3) )/(Lambda_z3**2 + abs(q4_q4))
     &         + ghz3_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z3**2                                   
     &         + ghz3_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z3**2                                                                                   
     &         + ghz3_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz3_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z3**4
     &         + ghz3_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z3**4
     &         + ghz3_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z3**4


      ghz4_dyn = ghz4 
     &         + ghz4_prime * Lambda_z4**4/( Lambda_z4**2 
     &         + abs(q3_q3) )/(Lambda_z4**2 + abs(q4_q4))
     &         + ghz4_prime2*(abs(q3_q3)+abs(q4_q4))/Lambda_z4**2                                   
     &         + ghz4_prime3*(abs(q3_q3)-abs(q4_q4))/Lambda_z4**2                                                                                   
     &         + ghz4_prime4*(        shat         )/Lambda_Q**2                                                                                    
     &         + ghz4_prime5*(abs(q3_q3)**2+abs(q4_q4)**2)/Lambda_z4**4
     &         + ghz4_prime6*(abs(q3_q3)**2-abs(q4_q4)**2)/Lambda_z4**4
     &         + ghz4_prime7*(abs(q3_q3)*abs(q4_q4))/Lambda_z4**4


      aa1 =ghz1_dyn*zmass**2/shat
     &     + (s(1,2)-s(3,4)-s(5,6))/shat*
     &       (ghz2_dyn
     &       +ghz3_dyn*(s(1,2)-s(3,4)-s(5,6))/4d0/LambdaBSM**2)
      aa2 =-2d0*ghz2_dyn
     &     -ghz3_dyn*(s(1,2)-s(3,4)-s(5,6))/2d0/LambdaBSM**2
      aa3 =-2d0*ghz4_dyn


      aa1 = aa1 / zmass**2 *wmass/(sinthw*(1d0-xw))*prop34*prop56 
      aa2 = aa2 / zmass**2 *wmass/(sinthw*(1d0-xw))*prop34*prop56 
      aa3 = aa3 / zmass**2 *wmass/(sinthw*(1d0-xw))*prop34*prop56 


c--- MARKUS: NEW amplitudes with anomalous couplings for Higgs decay
      H4l(1,1)=( aa1*FFa1(1,1) + aa2*FFa2(1,1) + aa3*FFa3(1,1) )*l1*l2
      H4l(2,1)=( aa1*FFa1(2,1) + aa2*FFa2(2,1) + aa3*FFa3(2,1) )*r1*l2
      H4l(1,2)=( aa1*FFa1(1,2) + aa2*FFa2(1,2) + aa3*FFa3(1,2) )*l1*r2
      H4l(2,2)=( aa1*FFa1(2,2) + aa2*FFa2(2,2) + aa3*FFa3(2,2) )*r1*r2     
 
      ELSE

c--- original MCFM amplitudes for the Higgs decay
c--- gives identical values to the above implementation if ghz1=1
      H4l(1,1)=za(3,5)*zb(4,6)*l1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,1)=za(4,5)*zb(3,6)*r1*l2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(1,2)=za(3,6)*zb(4,5)*l1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56
      H4l(2,2)=za(4,6)*zb(3,5)*r1*r2
     &        *wmass/(sinthw*(1d0-xw))*prop34*prop56


      ENDIF


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
      
