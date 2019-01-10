      subroutine HggggLO(i1,i2,i3,i4,
     .                   Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625)
c--- This routine is simply a wrapper to the four
c--- gluon Born routines. By changing "imode" below, one can
c--- switch between the calculation of Kauffman, Desai and Risal
c--- and one corresponding directly to the virtual amps (from H4pCode) 
      implicit none
      integer i1,i2,i3,i4,imode
      double precision Hgggg,Hgggg_1256,Hgggg_1265,Hgggg_1625,
     .                 Hgggg_1,Hgggg_1256_1,Hgggg_1265_1,Hgggg_1625_1,
     .                 Hgggg_2,Hgggg_1256_2,Hgggg_1265_2,Hgggg_1625_2
      
      imode=1
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square using results from KDR
c--- imode=3   ! compare the two calculations 
      
      if ((imode .eq. 1) .or. (imode .eq. 3)) then
        call h4gnew(i1,i2,i3,i4,
     .              Hgggg_1,Hgggg_1256_1,Hgggg_1265_1,Hgggg_1625_1)
      endif
      
      if ((imode .eq. 2) .or. (imode .eq. 3)) then
        call h4g(i1,i2,i3,i4,
     .           Hgggg_2,Hgggg_1256_2,Hgggg_1265_2,Hgggg_1625_2)
      endif

c--- Checked that h4g == h4gnew on 24/8/09
      if (imode .eq. 3) then
       write(6,99) 'i1,i2,i3,i4 Hgggg',i1,i2,i3,i4,
     .  1d0-Hgggg_1/Hgggg_2,1d0-Hgggg_1256_1/Hgggg_1256_2,
     .  1d0-Hgggg_1265_1/Hgggg_1265_2,1d0-Hgggg_1625_1/Hgggg_1625_2
      endif
      
      if (imode .eq. 2) then
        Hgggg=Hgggg_2
        Hgggg_1256=Hgggg_1256_2
        Hgggg_1265=Hgggg_1265_2
        Hgggg_1625=Hgggg_1625_2
      else
        Hgggg=Hgggg_1
        Hgggg_1256=Hgggg_1256_1
        Hgggg_1265=Hgggg_1265_1
        Hgggg_1625=Hgggg_1625_1
      endif
      
      return
      
   99 format(a18,4i3,4e12.4) 
      
      end
      

      subroutine HQAggLO(i1,i2,i3,i4,Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym)
c--- This routine is simply a wrapper to the four
c--- gluon Born routines. By changing "imode" below, one can
c--- switch between the calculation of Del Duca, Frizzo and Maltoni
c--- and one corresponding directly to the virtual amps (from H4pCode) 
      implicit none
      integer i1,i2,i3,i4,imode
      double precision Hqagg,Hqagg_ab,Hqagg_ba,Hqagg_sym,
     .                 Hqagg_1,Hqagg_ab_1,Hqagg_ba_1,Hqagg_sym_1,
     .                 Hqagg_2,Hqagg_ab_2,Hqagg_ba_2,Hqagg_sym_2
      
      imode=2
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square using results from DFM
c--- imode=3   ! compare the two calculations 
      
      if ((imode .eq. 1) .or. (imode .eq. 3)) then
        call HAQggnew(i1,i2,i3,i4,
     .                Hqagg_1,Hqagg_ab_1,Hqagg_ba_1,Hqagg_sym_1)
      endif
      
      if ((imode .eq. 2) .or. (imode .eq. 3)) then
        call hqqggdfm(i1,i2,i3,i4,
     .                Hqagg_2,Hqagg_ab_2,Hqagg_ba_2,Hqagg_sym_2)
      endif

c--- Checked that hqqggdfm == HAQggnew on 25/8/09
      if (imode .eq. 3) then
       write(6,99) 'i1,i2,i3,i4 Hqagg',i1,i2,i3,i4,
     .  1d0-Hqagg_1/Hqagg_2,1d0-Hqagg_ab_1/Hqagg_ab_2,
     .  1d0-Hqagg_ba_1/Hqagg_ba_2,1d0-Hqagg_sym_1/Hqagg_sym_2
      endif
      
      if (imode .eq. 2) then
        Hqagg=Hqagg_2
        Hqagg_ab=Hqagg_ab_2
        Hqagg_ba=Hqagg_ba_2
        Hqagg_sym=Hqagg_sym_2
      else
        Hqagg=Hqagg_1
        Hqagg_ab=Hqagg_ab_1
        Hqagg_ba=Hqagg_ba_1
        Hqagg_sym=Hqagg_sym_1
      endif
      
      return
      
   99 format(a18,4i3,4e12.4) 
      
      end
      

      subroutine HqarbLO(i1,i2,i3,i4,Hqarb)
c--- This routine is simply a wrapper to the four
c--- gluon Born routines. By changing "imode" below, one can
c--- switch between a squared calc. of Kauffman, Desai and Risal
c--- and one corresponding directly to the virtual amps (from H4pCode) 
      implicit none
      integer i1,i2,i3,i4,imode
      double precision Hqarb,Hqarb_1,Hqarb_2
      
      imode=2
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square using results from DFM
c--- imode=3   ! compare the two calculations 
      
      if ((imode .eq. 1) .or. (imode .eq. 3)) then
        call Ampsq_AQaq_nonid(i1,i2,i3,i4,Hqarb_1)
      endif
      
      if ((imode .eq. 2) .or. (imode .eq. 3)) then
        call H4qn(i1,i2,i3,i4,Hqarb_2)
      endif

c--- Checked that H4qn == Ampsq_AQaq_nonid on 23/10/09
      if (imode .eq. 3) then
       write(6,99) 'i1,i2,i3,i4 Hqarb',i1,i2,i3,i4,1d0-Hqarb_1/Hqarb_2
      endif
      
      if (imode .eq. 2) then
        Hqarb=Hqarb_2
      else
        Hqarb=Hqarb_1
      endif
      
      return
      
   99 format(a18,4i3,e12.4) 
      
      end
      

      subroutine HqaqaLO(i1,i2,i3,i4,Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i)
c--- This routine is simply a wrapper to the four
c--- gluon Born routines. By changing "imode" below, one can
c--- switch between a squared calc. of Kauffman, Desai and Risal
c--- and one corresponding directly to the virtual amps (from H4pCode) 
      implicit none
      integer i1,i2,i3,i4,imode
      double precision Hqaqa,Hqaqa_a,Hqaqa_b,Hqaqa_i,
     .                 Hqaqa_1,Hqaqa_a_1,Hqaqa_b_1,Hqaqa_i_1,
     .                 Hqaqa_2,Hqaqa_a_2,Hqaqa_b_2,Hqaqa_i_2
      
      imode=2
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square using results from DFM
c--- imode=3   ! compare the two calculations 
      
      if ((imode .eq. 1) .or. (imode .eq. 3)) then
        call Ampsq_AQaq_ident(i1,i2,i3,i4,
     .                        Hqaqa_1,Hqaqa_a_1,Hqaqa_b_1,Hqaqa_i_1)
      endif
      
      if ((imode .eq. 2) .or. (imode .eq. 3)) then
        call H4qi(i1,i2,i3,i4,Hqaqa_2,Hqaqa_a_2,Hqaqa_b_2,Hqaqa_i_2)
      endif

c--- Checked that H4qi == Ampsq_AQaq_ident on 24/10/09
      if (imode .eq. 3) then
       write(6,99) 'i1,i2,i3,i4 Hqaqa',i1,i2,i3,i4,
     .  1d0-Hqaqa_1/Hqaqa_2,1d0-Hqaqa_a_1/Hqaqa_a_2,
     .  1d0-Hqaqa_b_1/Hqaqa_b_2,1d0-Hqaqa_i_1/Hqaqa_i_2
      endif
      
      if (imode .eq. 2) then
        Hqaqa=Hqaqa_2
        Hqaqa_a=Hqaqa_a_2
        Hqaqa_b=Hqaqa_b_2
        Hqaqa_i=Hqaqa_i_2
      else
        Hqaqa=Hqaqa_1
        Hqaqa_a=Hqaqa_a_1
        Hqaqa_b=Hqaqa_b_1
        Hqaqa_i=Hqaqa_i_1
      endif
      
      return
      
   99 format(a18,4i3,4e12.4) 
      
      end
      
