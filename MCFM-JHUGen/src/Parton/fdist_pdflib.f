      subroutine fdist(ih,x,xmu,fx)
      implicit none
      double precision fx(-5:5),x,xmu
      double precision u_val,d_val,u_sea,d_sea,s_sea,c_sea,
     . b_sea,t_sea,gluon
      double precision Ctq3df,Ctq4Fn,Ctq5Pdf,Ctq6Pdf,Ctq5L
      integer mode,Iprtn,ih,Irt
c---  ih1=+1 proton 
c---  ih1=-1 pbar 

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-5,5
             fx(Iprtn)=0d0
          enddo
          return
      endif
 
      call structm(x,xmu,u_val,d_val,u_sea,d_sea,
     &                   s_sea,c_sea,b_sea,t_sea,gluon)
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
                     
      return

      end

  
      subroutine InitPDF(dummy)
      integer dummy

c--- this is a dummy routine that exists in LHAPDF only
      
      return
      end
      

