      subroutine fdist(ih_call,xx,xxmu,ffx)
      implicit none
      include 'types.f'
      include 'pdlabel.f'
      include 'vanillafiles.f'
      real(dp),intent(in)::xx,xxmu
      integer,intent(in)::ih_call
      real(dp),intent(out)::ffx(-5:5)
      double precision fx(-5:5),xZ,xA,eks98r,xmu_safe,x,xmu
      double precision u_val,d_val,u_sea,d_sea,s_sea,c_sea,b_sea,gluon
      double precision Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf,Ctq5L,CT10Pdf,
     & CT14Pdf
      double precision fxnnpdf(-6:7)
      integer mode,Iprtn,ih,iZ,iA,Irt
      logical first,nucleon
c--- extra variables for MSTW08 implementation
      character*72 prefix,checkpath
      double precision str,sbar,chm,cbar,bot,bbar,photon
      data first/.true./
      save first
c---  ih1=+1 proton
c---  ih1=-1 pbar

c---  Extended 5/24/05 to calculate nucleon parton distributions
c---  for a nucleus of atomic number Z and mass A via:
C---   K.J. Eskola, V.J. Kolhinen and C.A. Salgado,
C---   "The scale dependent nuclear effects in parton distributions for
C---   practical applications", Eur. Phys. J. C9 (1999) 61,
C---   JYFL-8/98, US-FT/14-98, hep-ph/9807297.

      x=dble(xx)
      xmu=dble(xxmu)
      nucleon=.false.
      if     (ih_call .gt. 1000d0) then
c--- nucleon distribution functions
        ih=1
        nucleon=.true.
        iA=mod(ih_call,1000)
        iZ=(ih_call-iZ)/1000
        xA=real(iA,dp)
        xZ=real(iZ,dp)
        if (first) then
        write(6,*)
        write(6,*)'******************* Nucleon beam *******************'
        write(6,*)'*                                                  *'
        write(6,76) iZ,iA
        write(6,*)'****************************************************'
        first=.false.
        endif
      elseif (abs(ih_call) .ne. 1) then
c--- error
        write(6,*) 'Input beam ih=',ih,' does not make sense!'
        stop
      else
c--- parton distribution functions
        ih=ih_call
      endif

C---set to zero if x out of range
      if ((x .le. 0d0) .or. (x .ge. 1d0)) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif

      if      ((pdlabel(1:5) .eq. 'mstw8') .or.
     &         (pdlabel(1:4) .eq. 'MMHT') ) then
            if         (pdlabel .eq. 'mstw8lo') then
            prefix = checkpath('Pdfdata/mstw2008lo') ! LO grid
            elseif     (pdlabel .eq. 'mstw8nl') then
            prefix = checkpath('Pdfdata/mstw2008nlo') ! NLO grid
            elseif     (pdlabel .eq. 'mstw8nn') then
            prefix = checkpath('Pdfdata/mstw2008nnlo') ! NNLO grid
            elseif     (pdlabel .eq. 'MMHT_lo') then
            prefix = checkpath('Pdfdata/mmht2014lo135') ! LO grid
            elseif     (pdlabel .eq. 'MMHT_nl') then
            prefix = checkpath('Pdfdata/mmht2014nlo120') ! NLO grid
            elseif     (pdlabel .eq. 'MMHT_nn') then
            prefix = checkpath('Pdfdata/mmht2014nnlo118') ! NNLO grid
            endif
            if (index(prefix,'mstw2008') > 0) then
              call GetAllPDFs(prefix,0,x,xmu,u_val,d_val,u_sea,d_sea,
     &                        str,sbar,chm,cbar,bot,bbar,gluon,photon)
            else
              call MMHTGetAllPDFs(prefix,0,x,xmu,u_val,d_val,u_sea,d_sea,
     &                        str,sbar,chm,cbar,bot,bbar,gluon,photon)
            endif
c-----assign MSTW to standard grid
            fx(0)=gluon/x
            if (ih.eq.1) then
               fx(1)=(d_val+d_sea)/x
               fx(2)=(u_val+u_sea)/x
               fx(-1)=d_sea/x
               fx(-2)=u_sea/x
               fx(+3)=str/x
               fx(+4)=chm/x
               fx(+5)=bot/x
               fx(-3)=sbar/x
               fx(-4)=cbar/x
               fx(-5)=bbar/x
            elseif(ih.eq.-1) then
               fx(-1)=(d_val+d_sea)/x
               fx(-2)=(u_val+u_sea)/x
               fx(+1)=d_sea/x
               fx(+2)=u_sea/x
               fx(-3)=str/x
               fx(-4)=chm/x
               fx(-5)=bot/x
               fx(+3)=sbar/x
               fx(+4)=cbar/x
               fx(+5)=bbar/x
            endif

      elseif ((pdlabel(1:3) .eq. 'mrs')
     .   .or. (pdlabel(2:4) .eq. 'mrs')) then

             if     (pdlabel .eq. 'mrs4nf3') then
             mode=1
             call mrst2004f3(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs4lf3') then
             mode=2
             call mrst2004f3(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs4nf4') then
             mode=1
             call mrst2004f4(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs4lf4') then
             mode=2
             call mrst2004f4(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs04nl') then
             mode=1
             call mrst2004(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs04nn') then
             mode=2
             call mrst2004(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs02nl') then
             mode=1
             call mrst2002(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs02nn') then
             mode=2
             call mrst2002(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs0119') then
             mode=1
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs0117') then
             mode=2
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs0121') then
             mode=3
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs01_j') then
             mode=4
             call mrst2001(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs01lo') then
c             write(6,*) 'here'
             mode=1
             call mrstlo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_1') then
             mode=1
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_2') then
             mode=2
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_3') then
             mode=3
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_4') then
             mode=4
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_5') then
             mode=5
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_6') then
             mode=6
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_7') then
             mode=7
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_8') then
             mode=8
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs99_9') then
             mode=9
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9910') then
             mode=10
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9911') then
             mode=11
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs9912') then
             mode=12
             call mrs99(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z1') then
             mode=1
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z2') then
             mode=2
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z3') then
             mode=3
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z4') then
             mode=4
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98z5') then
             mode=5
             call mrs98(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98l1') then
             mode=1
             call mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98l2') then
             mode=2
             call mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98l3') then
             mode=3
             call mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98l4') then
             mode=4
             call mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98l5') then
             mode=5
             call mrs98lo(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs98ht') then
             mode=1
             call mrs98ht(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r1') then
             mode=1
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r2') then
             mode=2
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r3') then
             mode=3
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs96r4') then
             mode=4
             call mrs96(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'hmrs90e') then
             mode=1
             call mrsebh(x,xmu,mode,u_val,d_val,u_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             d_sea=u_sea
             elseif (pdlabel .eq. 'hmrs90b') then
             mode=2
             call mrsebh(x,xmu,mode,u_val,d_val,u_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             d_sea=u_sea
             elseif (pdlabel .eq. 'mrs95ap') then
             mode=20
             call mrseb(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             elseif (pdlabel .eq. 'mrs95_g') then
             mode=21
             call mrseb(x,xmu,mode,u_val,d_val,u_sea,d_sea,
     &                          s_sea,c_sea,b_sea,gluon)
             endif
c-----assign mrs to standard grid
            fx(-5)=b_sea/x
            fx(-4)=c_sea/x
            fx(-3)=s_sea/x
            fx( 0)=gluon/x
            fx(+3)=fx(-3)
            fx(+4)=fx(-4)
            fx(+5)=fx(-5)
            if (ih.eq.1) then
               fx(1)=(d_val+d_sea)/x
               fx(2)=(u_val+u_sea)/x
               fx(-1)=d_sea/x
               fx(-2)=u_sea/x
            elseif(ih.eq.-1) then
               fx(-1)=(d_val+d_sea)/x
               fx(-2)=(u_val+u_sea)/x
               fx(+1)=d_sea/x
               fx(+2)=u_sea/x
            endif

      elseif (pdlabel(1:5) .eq. 'cteq3') then
C   1      CTEQ3M   Standard MSbar scheme   0.116
C   3      CTEQ3L   Leading Order           0.116
C   2      CTEQ3D   Standard DIS scheme     0.116
          if (pdlabel .eq. 'cteq3_m') then
             mode=1
          elseif (pdlabel .eq. 'cteq3_l') then
             mode=2
          elseif (pdlabel .eq. 'cteq3_d') then
             mode=3
          endif
             fx(-5)=Ctq3df(mode,-5,x,xmu,Irt)/x
             fx(-4)=Ctq3df(mode,-4,x,xmu,Irt)/x
             fx(-3)=Ctq3df(mode,-3,x,xmu,Irt)/x
             fx(0)=Ctq3df(mode,0,x,xmu,Irt)/x

             fx(+3)=fx(-3)
             fx(+4)=fx(-4)
             fx(+5)=fx(-5)
             if (ih.eq.1) then
               fx(-1)=Ctq3df(mode,-2,x,xmu,Irt)/x
               fx(-2)=Ctq3df(mode,-1,x,xmu,Irt)/x
               fx(1)=Ctq3df(mode,+2,x,xmu,Irt)/x+fx(-1)
               fx(2)=Ctq3df(mode,+1,x,xmu,Irt)/x+fx(-2)
             elseif(ih.eq.-1) then
               fx(1)=Ctq3df(mode,-2,x,xmu,Irt)/x
               fx(2)=Ctq3df(mode,-1,x,xmu,Irt)/x
               fx(-1)=Ctq3df(mode,+2,x,xmu,Irt)/x+fx(1)
               fx(-2)=Ctq3df(mode,+1,x,xmu,Irt)/x+fx(2)
             endif

      elseif (pdlabel(1:5) .eq. 'cteq4') then
C   1      CTEQ4M   Standard MSbar scheme   0.116        1.6      cteq4m.tbl
C   2      CTEQ4D   Standard DIS scheme     0.116        1.6      cteq4d.tbl
C   3      CTEQ4L   Leading Order           0.116        1.6      cteq4l.tbl
C   4      CTEQ4A1  Alpha_s series          0.110        1.6      cteq4a1.tbl
C   5      CTEQ4A2  Alpha_s series          0.113        1.6      cteq4a2.tbl
C   6      CTEQ4A3  same as CTEQ4M          0.116        1.6      cteq4m.tbl
C   7      CTEQ4A4  Alpha_s series          0.119        1.6      cteq4a4.tbl
C   8      CTEQ4A5  Alpha_s series          0.122        1.6      cteq4a5.tbl
C   9      CTEQ4HJ  High Jet                0.116        1.6      cteq4hj.tbl
C   10     CTEQ4LQ  Low Q0                  0.114        0.7      cteq4lq.tbl

          if (pdlabel .eq. 'cteq4_m') then
             mode=1
          elseif (pdlabel .eq. 'cteq4_d') then
             mode=2
          elseif (pdlabel .eq. 'cteq4_l') then
             mode=3
          elseif (pdlabel .eq. 'cteq4a1') then
             mode=4
          elseif (pdlabel .eq. 'cteq4a2') then
             mode=5
          elseif (pdlabel .eq. 'cteq4a3') then
             mode=6
          elseif (pdlabel .eq. 'cteq4a4') then
             mode=7
          elseif (pdlabel .eq. 'cteq4a5') then
             mode=8
          elseif (pdlabel .eq. 'cteq4hj') then
             mode=9
          elseif (pdlabel .eq. 'cteq4lq') then
             mode=10
          endif

             fx(-5)=Ctq4Fn(mode,-5,x,xmu)
             fx(-4)=Ctq4Fn(mode,-4,x,xmu)
             fx(-3)=Ctq4Fn(mode,-3,x,xmu)

             fx(0)=Ctq4Fn(mode,0,x,xmu)

             fx(+3)=Ctq4Fn(mode,+3,x,xmu)
             fx(+4)=Ctq4Fn(mode,+4,x,xmu)
             fx(+5)=Ctq4Fn(mode,+5,x,xmu)
             if (ih.eq.1) then
               fx(1)=Ctq4Fn(mode,+2,x,xmu)
               fx(2)=Ctq4Fn(mode,+1,x,xmu)
               fx(-1)=Ctq4Fn(mode,-2,x,xmu)
               fx(-2)=Ctq4Fn(mode,-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=Ctq4Fn(mode,-2,x,xmu)
               fx(2)=Ctq4Fn(mode,-1,x,xmu)
               fx(-1)=Ctq4Fn(mode,+2,x,xmu)
               fx(-2)=Ctq4Fn(mode,+1,x,xmu)
             endif

      elseif (pdlabel .eq. 'cteq5l1') then
             fx(-5)=Ctq5L(-5,x,xmu)
             fx(-4)=Ctq5L(-4,x,xmu)
             fx(-3)=Ctq5L(-3,x,xmu)

             fx(0)=Ctq5L(0,x,xmu)

             fx(+3)=Ctq5L(+3,x,xmu)
             fx(+4)=Ctq5L(+4,x,xmu)
             fx(+5)=Ctq5L(+5,x,xmu)

             if (ih.eq.1) then
               fx(1)=Ctq5L(+2,x,xmu)
               fx(2)=Ctq5L(+1,x,xmu)
               fx(-1)=Ctq5L(-2,x,xmu)
               fx(-2)=Ctq5L(-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=Ctq5L(-2,x,xmu)
               fx(2)=Ctq5L(-1,x,xmu)
               fx(-1)=Ctq5L(+2,x,xmu)
               fx(-2)=Ctq5L(+1,x,xmu)
             endif

      elseif ((pdlabel(1:5) .eq. 'cteq5') .or.
     .        (pdlabel(1:4) .eq. 'ctq5')) then

             fx(-5)=Ctq5Pdf(-5,x,xmu)
             fx(-4)=Ctq5Pdf(-4,x,xmu)
             fx(-3)=Ctq5Pdf(-3,x,xmu)

             fx(0)=Ctq5Pdf(0,x,xmu)

             fx(+3)=Ctq5Pdf(+3,x,xmu)
             fx(+4)=Ctq5Pdf(+4,x,xmu)
             fx(+5)=Ctq5Pdf(+5,x,xmu)

             if (ih.eq.1) then
               fx(1)=Ctq5Pdf(+2,x,xmu)
               fx(2)=Ctq5Pdf(+1,x,xmu)
               fx(-1)=Ctq5Pdf(-2,x,xmu)
               fx(-2)=Ctq5Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=Ctq5Pdf(-2,x,xmu)
               fx(2)=Ctq5Pdf(-1,x,xmu)
               fx(-1)=Ctq5Pdf(+2,x,xmu)
               fx(-2)=Ctq5Pdf(+1,x,xmu)
             endif

      elseif (pdlabel(1:5) .eq. 'cteq6') then

             fx(-5)=Ctq6Pdf(-5,x,xmu)
             fx(-4)=Ctq6Pdf(-4,x,xmu)
             fx(-3)=Ctq6Pdf(-3,x,xmu)

             fx(0)=Ctq6Pdf(0,x,xmu)

             fx(+3)=Ctq6Pdf(+3,x,xmu)
             fx(+4)=Ctq6Pdf(+4,x,xmu)
             fx(+5)=Ctq6Pdf(+5,x,xmu)

             if (ih.eq.1) then
               fx(1)=Ctq6Pdf(+2,x,xmu)
               fx(2)=Ctq6Pdf(+1,x,xmu)
               fx(-1)=Ctq6Pdf(-2,x,xmu)
               fx(-2)=Ctq6Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=Ctq6Pdf(-2,x,xmu)
               fx(2)=Ctq6Pdf(-1,x,xmu)
               fx(-1)=Ctq6Pdf(+2,x,xmu)
               fx(-2)=Ctq6Pdf(+1,x,xmu)
             endif

      elseif (pdlabel(1:4) .eq. 'CT10') then

             fx(-5)=CT10Pdf(-5,x,xmu)
             fx(-4)=CT10Pdf(-4,x,xmu)
             fx(-3)=CT10Pdf(-3,x,xmu)

             fx(0)=CT10Pdf(0,x,xmu)

             fx(+3)=CT10Pdf(+3,x,xmu)
             fx(+4)=CT10Pdf(+4,x,xmu)
             fx(+5)=CT10Pdf(+5,x,xmu)

             if (ih.eq.1) then
               fx(1)=CT10Pdf(+2,x,xmu)
               fx(2)=CT10Pdf(+1,x,xmu)
               fx(-1)=CT10Pdf(-2,x,xmu)
               fx(-2)=CT10Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=CT10Pdf(-2,x,xmu)
               fx(2)=CT10Pdf(-1,x,xmu)
               fx(-1)=CT10Pdf(+2,x,xmu)
               fx(-2)=CT10Pdf(+1,x,xmu)
             endif

      elseif (pdlabel(1:4) .eq. 'CT14') then

             fx(-5)=CT14Pdf(-5,x,xmu)
             fx(-4)=CT14Pdf(-4,x,xmu)
             fx(-3)=CT14Pdf(-3,x,xmu)

             fx(0)=CT14Pdf(0,x,xmu)

             fx(+3)=CT14Pdf(+3,x,xmu)
             fx(+4)=CT14Pdf(+4,x,xmu)
             fx(+5)=CT14Pdf(+5,x,xmu)

             if (ih.eq.1) then
               fx(1)=CT14Pdf(+2,x,xmu)
               fx(2)=CT14Pdf(+1,x,xmu)
               fx(-1)=CT14Pdf(-2,x,xmu)
               fx(-2)=CT14Pdf(-1,x,xmu)
             elseif(ih.eq.-1) then
               fx(1)=CT14Pdf(-2,x,xmu)
               fx(2)=CT14Pdf(-1,x,xmu)
               fx(-1)=CT14Pdf(+2,x,xmu)
               fx(-2)=CT14Pdf(+1,x,xmu)
             endif

      elseif (pdlabel(1:2) .eq. 'NN') then

             call NNevolvePDF(x,xmu,fxnnpdf)
             fx(-5:5)=fxnnpdf(-5:5)/x
             if (ih == -1) then
               fx(1)=fxnnpdf(-1)/x
               fx(2)=fxnnpdf(-2)/x
               fx(-1)=fxnnpdf(1)/x
               fx(-2)=fxnnpdf(2)/x
             endif

c--- NEW ATTEMPT
      elseif (pdlabel(1:5) .eq. 'mtung') then
            if     (pdlabel .eq. 'mtungs1') then
              mode=1
            elseif (pdlabel .eq. 'mtunge1') then
              mode=2
            elseif (pdlabel .eq. 'mtungb1') then
              mode=3
            elseif (pdlabel .eq. 'mtungb2') then
              mode=4
            elseif (pdlabel .eq. 'mtungn1') then
              mode=5
            endif
            call mt(x,xmu,mode,u_val,d_val,
     .               u_sea,s_sea,c_sea,b_sea,gluon)
            d_sea=u_sea
c-----assign to standard grid
            fx(-5)=b_sea/x
            fx(-4)=c_sea/x
            fx(-3)=s_sea/x
            fx( 0)=gluon/x
            fx(+3)=fx(-3)
            fx(+4)=fx(-4)
            fx(+5)=fx(-5)
            if (ih.eq.1) then
               fx(1)=(d_val+d_sea)/x
               fx(2)=(u_val+u_sea)/x
               fx(-1)=d_sea/x
               fx(-2)=u_sea/x
            elseif(ih.eq.-1) then
               fx(-1)=(d_val+d_sea)/x
               fx(-2)=(u_val+u_sea)/x
               fx(+1)=d_sea/x
               fx(+2)=u_sea/x
            endif

      else
          write(6,*) 'Unimplemented pdf distribution'
          write(6,*) 'pdlabel= ',pdlabel
          write(6,*) 'Implemented are: ',
     . 'mrs4nf3,mrs4lf3,mrs4nf4,mrs4lf4,',
     . 'mrs04nl,mrs04nn,mrs02nl,mrs02nn,',
     . 'mrs0119,mrs0177,mrs0121,mrs01_j,',
     . 'mrs99_1,mrs99_2,mrs99_3,mrs99_4,mrs99_5,mrs99_6,',
     . 'mrs99_7,mrs99_8,mrs99_9,mrs9910,mrs9911,mrs9912,',
     . 'mrs98z1,',
     . 'mrs98z2,',
     . 'mrs98z3,',
     . 'mrs98z4,',
     . 'mrs98z5,',
     . 'mrs96r1,',
     . 'mrs96r2,',
     . 'mrs96r3,',
     . 'mrs96r4,',
     . 'hmrs90e,',
     . 'hmrs90b,',
     . 'cteq4_m,',
     . 'cteq5_m,',
     . 'cteq5_d,',
     . 'cteq5_l,',
     . 'cteq5l1,',
     . 'cteq5hj,',
     . 'cteq5hq,',
     . 'cteq5f3,',
     . 'cteq5f4,',
     . 'cteq6l1,',
     . 'cteq6_l,',
     . 'cteq6_d,',
     . 'cteq6_m,',
     . 'cteq6l1,',
     . 'cteq61m,',
     . 'cteq66m,',
     . 'mtungs1,',
     . 'mtunge1,',
     . 'mtungb1,',
     . 'mtungb2,',
     . 'mtungn1,',
     . 'CT10.00,',
     . 'CT14.LL,',
     . 'CT14.NL,',
     . 'NN2.3NL',
     . 'NN2.3NN,'

         stop
      endif

c--- now perform the corrections for a nucleon beam, if necessary
      if (nucleon) then
c        write(6,*) 'x,Q',x,xmu
c        write(6,*) 'before'
c        write(6,*) 'fx(+5)=',fx(+5)
c        write(6,*) 'fx(+4)=',fx(+4)
c        write(6,*) 'fx(+3)=',fx(+3)
c        write(6,*) 'fx(+2)=',fx(+2)
c        write(6,*) 'fx(+1)=',fx(+1)
c        write(6,*) 'fx( 0)=',fx( 0)
c        write(6,*) 'fx(-1)=',fx(-1)
c        write(6,*) 'fx(-2)=',fx(-2)
c        write(6,*) 'fx(-3)=',fx(-3)
c        write(6,*) 'fx(-4)=',fx(-4)
c        write(6,*) 'fx(-5)=',fx(-5)
c--- extract valence and sea components from the parton distributions
        u_sea=fx(-2)
        d_sea=fx(-1)
        u_val=fx(2)-u_sea
        d_val=fx(1)-d_sea
c--- apply scale factors from eks98r, with a maximum scale of 100 GeV
        xmu_safe=min(xmu,100d0)
        u_val=u_val*eks98r(x,xmu_safe,xA,1)
        d_val=d_val*eks98r(x,xmu_safe,xA,2)
        u_sea=u_sea*eks98r(x,xmu_safe,xA,3)
        d_sea=d_sea*eks98r(x,xmu_safe,xA,4)
        s_sea=fx(3)*eks98r(x,xmu_safe,xA,5)
        c_sea=fx(4)*eks98r(x,xmu_safe,xA,6)
        b_sea=fx(5)*eks98r(x,xmu_safe,xA,7)
        gluon=fx(0)*eks98r(x,xmu_safe,xA,8)
c--- write new nucleon distributions
        fx(1)=(xZ*(d_val+d_sea)+(xA-xZ)*(u_val+u_sea))/xA
        fx(2)=(xZ*(u_val+u_sea)+(xA-xZ)*(d_val+d_sea))/xA
        fx(-1)=(xZ*d_sea+(xA-xZ)*u_sea)/xA
        fx(-2)=(xZ*u_sea+(xA-xZ)*d_sea)/xA
        fx(-3)=s_sea
        fx(-4)=c_sea
        fx(-5)=b_sea
        fx(+3)=s_sea
        fx(+4)=c_sea
        fx(+5)=b_sea
        fx(0)=gluon
      endif

        ffx(:)=real(fx(:))

!        ffx(0)=0d0      ! DEBUG: remove gluon pdfs
!        ffx(1:5)=0d0    ! DEBUG: remove quark pdfs
!        ffx(-5:-1)=0d0  ! DEBUG: remove antiquark pdfs

      return

   76 format(' *    (Atomic number, mass) = (Z,A) = (',
     . i4,',',i4,')   *')

      end


      subroutine InitPDF(dummy)
      integer dummy

c--- this is a dummy routine that exists in LHAPDF only

      return
      end

