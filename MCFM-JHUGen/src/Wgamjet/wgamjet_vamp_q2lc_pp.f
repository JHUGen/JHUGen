!==== this routine builds the Leading Color one-loop amplitude for 
!==== Wgam jet ++ (34) pieces


      double complex function wgamjet_vamp_q2lc_pp(i1,i2,i3,i4,i5,i6
     &,za,zb)
      include 'constants.f' 
      include 'zprods_decl.f' 
      include 'sprods_com.f'
      integer i1,i2,i3,i4,i5,i6
      double complex Atree_q2pp
      double complex d2me1,d1m1,d1m2 
      double complex Dcoeff1m_q2lcpp_s134,Dcoeff1m_q2lcpp_s126
      double complex Dcoeff2me_q2lcpp_s26s34s126s134
      double complex Bcoeff_q2lcpp_s34
      double complex Bcoeff_q2lcpp_s16
      double complex b34,b126,b16

      logical do_check 
      common/do_check_wgamj/do_check
      double complex zab2
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      Atree_q2pp= (za(i1,i3)**2*zab2(i2,i3,i4,i5))/
     -  (za(i1,i6)*za(i2,i5)*za(i3,i4)*
     -    za(i6,i2)*zab2(i5,i3,i4,i5))

!===== KC normalization 
!      Atree_q2pp=Atree_q2pp

!======== box coefficients       
      d2me1=Dcoeff2me_q2lcpp_s26s34s126s134(i1,i2,i3,i4,i5,i6,za,zb)
      d1m1=Dcoeff1m_q2lcpp_s126(i1,i2,i3,i4,i5,i6,za,zb)
      d1m2=Dcoeff1m_q2lcpp_s134(i1,i2,i3,i4,i5,i6,za,zb)


!======== bubble coefficients 
      b34=Bcoeff_q2lcpp_s34(i1,i2,i3,i4,i5,i6,za,zb)
      b16=Bcoeff_q2lcpp_s16(i1,i2,i3,i4,i5,i6,za,zb)
      b126=3d0/2d0*Atree_q2pp-b34-b16

      if(do_check) then
!==== print out information for KC 
         write(6,*) ' LO ',im*Atree_q2pp
         write(6,*) 'im * Box 2me /LO ',im*d2me1/Atree_q2pp
         write(6,*) 'im * Box 1m 1 /LO ',im*d1m1/Atree_q2pp
         write(6,*) 'im*Box 1m 2 /LO ',im*d1m2/Atree_q2pp
         write(6,*) 'im*Bub (34) /LO ',im*b34/Atree_q2pp
         write(6,*) 'im*Bub (16) /LO ',im*b16/Atree_q2pp
         write(6,*) 'im*Bub (126) /LO ',im*b126/Atree_q2pp
         
      endif

      return 
      end
      
