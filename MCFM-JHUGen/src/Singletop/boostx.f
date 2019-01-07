      subroutine boostx(p_in,pt,ptt,p_out)
      implicit none
c--- Boost input vector p_in to output vector p_out using the same
c--- transformation as required to boost massive vector pt to ptt
      double precision p_in(4),pt(4),ptt(4),p_out(4),
     . p_tmp(4),beta(3),mass,gam,bdotp
      integer j
    
      mass=pt(4)**2-pt(1)**2-pt(2)**2-pt(3)**2  
      if (mass .lt. 0d0) then
        write(6,*) 'mass**2 .lt. 0 in boostx.f, mass**2=',mass
        stop
      endif
      mass=dsqrt(mass)

c--- boost to the rest frame of pt
      gam=pt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=-pt(j)/pt(4)
        bdotp=bdotp+beta(j)*p_in(j)
      enddo
      p_tmp(4)=gam*(p_in(4)+bdotp)
      do j=1,3
        p_tmp(j)=p_in(j)+gam*beta(j)/(1d0+gam)*(p_in(4)+p_tmp(4))
      enddo     

c--- boost from rest frame of pt to frame in which pt is identical
c--- with ptt, thus completing the transformation          
      gam=ptt(4)/mass

      bdotp=0d0
      do j=1,3
        beta(j)=+ptt(j)/ptt(4)
        bdotp=bdotp+beta(j)*p_tmp(j)
      enddo
      p_out(4)=gam*(p_tmp(4)+bdotp)
      do j=1,3
        p_out(j)=p_tmp(j)+gam*beta(j)/(1d0+gam)*(p_out(4)+p_tmp(4))
      enddo

      return
      end
      

