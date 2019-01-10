      double precision function Hqaqavsqanal(i1,i2,i3,i4)
c--- This routine is simply a wrapper to the identical)four
c--- quark virtual routines. By changing "imode" below, one can
c--- switch between a squared ME calculation and one using amplitudes
      implicit none
      integer i1,i2,i3,i4,imode
      double precision Hqaqavsq,Ampvirtsq_AQaq_ident,res_sq,res_amp
      
      imode=2
c--- imode=1   ! compute square using amplitudes from H4pCode      
c--- imode=2   ! compute square directly using results from EGZ
c--- imode=3   ! compare the two calculations 
      
      if ((imode .eq. 1) .or. (imode .eq. 3)) then
        res_amp=Ampvirtsq_AQaq_ident(i1,i2,i3,i4)
      endif
      
      if ((imode .eq. 2) .or. (imode .eq. 3)) then
        res_sq=Hqaqavsq(i1,i2,i3,i4)
      endif

c--- Checked that Hqaqavsq == Ampvirtsq_AQaq_ident on 27/10/09
      if (imode .eq. 3) then
       write(6,99) 'i1,i2,i3,i4 Hqaqavsq',i1,i2,i3,i4,1d0-res_sq/res_amp
      endif
      
      if (imode .eq. 2) then
        Hqaqavsqanal=res_sq
      else
        Hqaqavsqanal=res_amp
      endif
      
      return
      
   99 format(a20,4i3,e21.12) 
      
      end
      
