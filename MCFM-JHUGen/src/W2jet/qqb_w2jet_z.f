      subroutine qqb_w2jet_z(p,z)
************************************************************************
*     Author: J.M. Campbell                                            *
*     September, 2001.                                                   *
*     Additions Aug. 2001, for 4Q piece                                *
************************************************************************
*                                                                      *
*     The value of COLOURCHOICE determines which colour structures     *
*     are included in the terms for the QQGG piece                     *
*                                                                      *
************************************************************************
      implicit none
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'agq.f'
      include 'PR_cs_new.f'
      include 'lc.f'
      include 'flags.f'
      double precision z,p(mxpart,4),dot,xninv
      double precision xl12,xl15,xl16,xl25,xl26,xl56
      double precision 
     .                 ii_qq,ii_qg,ii_gq,ii_gg,
     .                 if_qq,if_gg,
     .                 fi_qq,fi_gg,
     .                 ff_qq,ff_gg
      double precision tempgq,tempqg
      integer is
      
      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file input.DAT'
        write(6,*) 'Gflag=',Gflag
        write(6,*) 'Qflag=',Qflag
        write(6,*) 'Failed in qqb_w2jet_z.f'
        stop
      endif

      xl12=dlog(+two*dot(p,1,2)/musq)
      xl15=dlog(-two*dot(p,1,5)/musq)
      xl16=dlog(-two*dot(p,1,6)/musq)
      xl25=dlog(-two*dot(p,2,5)/musq)
      xl26=dlog(-two*dot(p,2,6)/musq)
      xl56=dlog(+two*dot(p,5,6)/musq)

************************************************************************
*     Contributions from QQGG matrix elements                          *
************************************************************************            
      if (Gflag) then      
c--- QUARK-ANTIQUARK contributions
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,a,1,is)=ason4pi*xn*(if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R1(q,q,a,2,is)=ason4pi*xn*(if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R2(a,a,q,1,is)=ason4pi*xn*(if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R2(a,a,q,2,is)=ason4pi*xn*(if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      enddo
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,a,0,is)=R1(q,q,a,0,is)+ason4pi*xn*
     .  (if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2d0
     .  +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2d0)
      R1(q,q,a,1,is)=R1(q,q,a,1,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R1(q,q,a,2,is)=R1(q,q,a,2,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R2(a,a,q,0,is)=R2(a,a,q,0,is)+ason4pi*xn*
     .  (if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2d0
     .  +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2d0)
      R2(a,a,q,1,is)=R2(a,a,q,1,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R2(a,a,q,2,is)=R2(a,a,q,2,is)-ason4pi/xn*ii_qq(z,xl12,is)
      enddo
      endif

      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,a,0,is)=R1(q,q,a,0,is)-ason4pi*(xn+1d0/xn)
     .  *ii_qq(z,xl12,is)
      R2(a,a,q,0,is)=R2(a,a,q,0,is)-ason4pi*(xn+1d0/xn)
     .  *ii_qq(z,xl12,is)
      enddo
      endif

c--- ANTIQUARK-QUARK contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*ff_gg(z,xl56,2)*( (1) + (2) )
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(a,a,q,1,is)=ason4pi*xn*(if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R1(a,a,q,2,is)=ason4pi*xn*(if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R2(q,q,a,1,is)=ason4pi*xn*(if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      R2(q,q,a,2,is)=ason4pi*xn*(if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2d0
     .                          +ff_gg(z,xl56,is)/2d0)
      enddo
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(a,a,q,0,is)=R1(a,a,q,0,is)+ason4pi*xn*
     .  (if_qq(z,xl15,is)+fi_gg(z,xl15,is)/2d0
     .  +if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2d0)
      R1(a,a,q,1,is)=R1(a,a,q,1,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R1(a,a,q,2,is)=R1(a,a,q,2,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R2(q,q,a,0,is)=R2(q,q,a,0,is)+ason4pi*xn*
     .  (if_qq(z,xl25,is)+fi_gg(z,xl25,is)/2d0
     .  +if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2d0)
      R2(q,q,a,1,is)=R2(q,q,a,1,is)-ason4pi/xn*ii_qq(z,xl12,is)
      R2(q,q,a,2,is)=R2(q,q,a,2,is)-ason4pi/xn*ii_qq(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(a,a,q,0,is)=R1(a,a,q,0,is)-ason4pi*(xn+1d0/xn)
     .  *ii_qq(z,xl12,is)
      R2(q,q,a,0,is)=R2(q,q,a,0,is)-ason4pi*(xn+1d0/xn)
     .  *ii_qq(z,xl12,is)
      enddo
      endif

c--- GLUON-GLUON contributions
c--- no additional final-final pieces
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(g,g,g,1,is)=ason4pi*xn*(if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     .                           +ii_gg(z,xl12,is))
      R1(g,g,g,2,is)=ason4pi*xn*(if_gg(z,xl16,is)+fi_qq(z,xl16,is)
     .                           +ii_gg(z,xl12,is))
      R2(g,g,g,1,is)=ason4pi*xn*(if_gg(z,xl26,is)+fi_qq(z,xl26,is)
     .                           +ii_gg(z,xl12,is))
      R2(g,g,g,2,is)=ason4pi*xn*(if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     .                           +ii_gg(z,xl12,is))
      R1(q,g,g,1,is)=ason4pi*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R1(q,g,g,2,is)=ason4pi*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,1,is)=ason4pi*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,2,is)=ason4pi*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(g,g,g,0,is)=R1(g,g,g,0,is)+ason4pi*xn*
     .  (if_gg(z,xl15,is)+fi_qq(z,xl15,is)
     .  +if_gg(z,xl16,is)+fi_qq(z,xl16,is))
      R2(g,g,g,0,is)=R2(g,g,g,0,is)+ason4pi*xn*
     .  (if_gg(z,xl25,is)+fi_qq(z,xl25,is)
     .  +if_gg(z,xl26,is)+fi_qq(z,xl26,is))
      R1(g,g,g,1,is)=R1(g,g,g,1,is)-ason4pi/xn*ff_qq(z,xl56,is)
      R2(g,g,g,1,is)=R2(g,g,g,1,is)-ason4pi/xn*ff_qq(z,xl56,is)
      R1(g,g,g,2,is)=R1(g,g,g,2,is)-ason4pi/xn*ff_qq(z,xl56,is)
      R2(g,g,g,2,is)=R2(g,g,g,2,is)-ason4pi/xn*ff_qq(z,xl56,is)
      R1(q,g,g,1,is)=R1(q,g,g,1,is)-
     . ason4pi/xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R1(q,g,g,2,is)=R1(q,g,g,2,is)-
     . ason4pi/xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R1(q,g,g,0,is)=R1(q,g,g,0,is)+
     . ason4pi*2d0*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,1,is)=R2(q,g,g,1,is)-
     . ason4pi/xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,2,is)=R2(q,g,g,2,is)-
     . ason4pi/xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,0,is)=R2(q,g,g,0,is)+
     . ason4pi*2d0*xn*(avegg/aveqg)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(g,g,g,0,is)=R1(g,g,g,0,is)-ason4pi*(xn+1d0/xn)*ff_qq(z,xl56,is)
      R2(g,g,g,0,is)=R2(g,g,g,0,is)-ason4pi*(xn+1d0/xn)*ff_qq(z,xl56,is)
      R1(q,g,g,0,is)=R1(q,g,g,0,is)-
     . ason4pi*(xn+1d0/xn)*(avegg/aveqg)*ii_qg(z,xl12,is)
      R2(q,g,g,0,is)=R2(q,g,g,0,is)-
     . ason4pi*(xn+1d0/xn)*(avegg/aveqg)*ii_qg(z,xl12,is)
      enddo
      endif

c--- QUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,g,1,is)=ason4pi*half*xn*(ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     .              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R1(q,q,g,2,is)=ason4pi*xn*(if_qq(z,xl16,is)+fi_gg(z,xl16,is)/2d0)
      R2(g,g,q,1,is)=ason4pi*half*xn
     . *(2d0*if_gg(z,xl26,is)+fi_gg(z,xl26,is)
     .      +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     .      +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R2(g,g,q,2,is)=ason4pi*xn*(if_gg(z,xl26,is)+fi_gg(z,xl26,is)/2d0
     .                          +if_gg(z,xl25,is)+fi_qq(z,xl25,is))
      R2(a,g,q,1,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,2,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,g,1,is)=R1(q,q,g,1,is)
     .           -ason4pi/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R1(q,q,g,2,is)=R1(q,q,g,2,is)
     .           -ason4pi/xn*(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R1(q,q,g,0,is)=R1(q,q,g,0,is)+ason4pi*half*xn*
     .            (2d0*if_qq(z,xl16,is)+fi_gg(z,xl16,is)
     .            +ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     .            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R2(g,g,q,0,is)=R2(g,g,q,0,is)+ason4pi*half*xn*
     .            (2d0*if_gg(z,xl25,is)+2d0*fi_qq(z,xl25,is)
     .            +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     .            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R2(a,g,q,1,is)=R2(a,g,q,1,is)-
     . ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,2,is)=R2(a,g,q,2,is)-
     . ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R2(a,g,q,0,is)=R2(a,g,q,0,is)+
     . ason4pi*2d0*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(q,q,g,0,is)=R1(q,q,g,0,is)-ason4pi*(xn+1d0/xn)
     .                       *(if_qq(z,xl15,is)+fi_qq(z,xl15,is))
      R2(a,g,q,0,is)=R2(a,g,q,0,is)-
     . ason4pi*(xn+1d0/xn)*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif

c--- GLUON-QUARK contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      if ((colourchoice .eq. 1) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R1(g,g,q,1,is)=ason4pi*half*xn*(ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     .                           +2d0*if_gg(z,xl16,is)+fi_gg(z,xl16,is)
     .              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R1(g,g,q,2,is)=ason4pi*xn*(if_gg(z,xl16,is)+fi_gg(z,xl16,is)/2d0
     .                           +if_gg(z,xl15,is)+fi_qq(z,xl15,is))
      R2(q,q,g,1,is)=ason4pi*half*xn*(ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     .              +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R2(q,q,g,2,is)=ason4pi*xn*(if_qq(z,xl26,is)+fi_gg(z,xl26,is)/2d0)
      R1(a,g,q,1,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,2,is)=ason4pi*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 2) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R2(q,q,g,1,is)=R2(q,q,g,1,is)
     .           -ason4pi/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R2(q,q,g,2,is)=R2(q,q,g,2,is)
     .           -ason4pi/xn*(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R1(g,g,q,0,is)=R1(g,g,q,0,is)+ason4pi*half*xn*
     .            (2d0*if_gg(z,xl15,is)+2d0*fi_qq(z,xl15,is)
     .            +ii_gg(z,xl12,is)+ii_gg(z,xl12,is)
     .            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R2(q,q,g,0,is)=R2(q,q,g,0,is)+ason4pi*half*xn*
     .            (2d0*if_qq(z,xl26,is)+fi_gg(z,xl26,is)
     .            +ii_qq(z,xl12,is)+ii_qq(z,xl12,is)
     .            +ff_qq(z,xl56,is)+ff_gg(z,xl56,is)/2d0)
      R1(a,g,q,1,is)=R1(a,g,q,1,is)-
     . ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,2,is)=R1(a,g,q,2,is)-
     . ason4pi/xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      R1(a,g,q,0,is)=R1(a,g,q,0,is)+
     . ason4pi*2d0*xn*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif
      if ((colourchoice .eq. 3) .or. (colourchoice .eq. 0)) then
      do is=1,3
      R2(q,q,g,0,is)=R2(q,q,g,0,is)-ason4pi*(xn+1d0/xn)
     .                       *(if_qq(z,xl25,is)+fi_qq(z,xl25,is))
      R1(a,g,q,0,is)=R1(a,g,q,0,is)-
     . ason4pi*(xn+1d0/xn)*(aveqg/aveqq)*ii_qg(z,xl12,is)
      enddo
      endif

c--- GLUON-ANTIQUARK contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      do is=1,3
      R1(g,g,a,1,is)= R1(g,g,q,2,is)
      R1(g,g,a,2,is)= R1(g,g,q,1,is)
      R1(g,g,a,0,is)= R1(g,g,q,0,is)
      R2(a,a,g,1,is)= R2(q,q,g,2,is)    
      R2(a,a,g,2,is)= R2(q,q,g,1,is)
      R2(a,a,g,0,is)= R2(q,q,g,0,is)
      R1(q,g,a,0,is)= R1(a,g,q,0,is)
      R1(q,g,a,1,is)= R1(a,g,q,1,is)
      R1(q,g,a,2,is)= R1(a,g,q,2,is)
      enddo

c--- ANTIQUARK-GLUON contributions
c--- additional final-final pieces that are 0:
c--- ason4pi*xn*(ff_qq(z,xl56,2) + ff_gg(z,xl56,2)/2d0) * (1)
      do is=1,3
      R1(a,a,g,1,is)= R1(q,q,g,2,is)
      R1(a,a,g,2,is)= R1(q,q,g,1,is)
      R1(a,a,g,0,is)= R1(q,q,g,0,is)
      R2(g,g,a,1,is)= R2(g,g,q,2,is)
      R2(g,g,a,2,is)= R2(g,g,q,1,is)
      R2(g,g,a,0,is)= R2(g,g,q,0,is)
      R2(q,g,a,0,is)= R2(a,g,q,0,is)
      R2(q,g,a,1,is)= R2(a,g,q,1,is)
      R2(q,g,a,2,is)= R2(a,g,q,2,is)
      enddo
       
            
C--- off-diagonal Qflag -> Gflag pieces

      do is=1,3
      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
      R1(g,q,a,0,is)=tempgq
      R1(g,q,a,1,is)=tempgq
      R1(g,q,a,2,is)=tempgq
      R2(g,q,a,0,is)=tempgq
      R2(g,q,a,1,is)=tempgq
      R2(g,q,a,2,is)=tempgq

      R1(g,q,g,0,is)=tempgq
      R1(g,q,g,1,is)=tempgq
      R1(g,q,g,2,is)=tempgq
      R2(g,q,g,0,is)=tempgq
      R2(g,q,g,1,is)=tempgq
      R2(g,q,g,2,is)=tempgq

      R1(g,q,q,0,is)=tempgq
      R1(g,q,q,1,is)=tempgq
      R1(g,q,q,2,is)=tempgq
      R2(g,q,q,0,is)=tempgq
      R2(g,q,q,1,is)=tempgq
      R2(g,q,q,2,is)=tempgq

      R1(g,a,a,0,is)=tempgq
      R1(g,a,a,1,is)=tempgq
      R1(g,a,a,2,is)=tempgq
      R2(g,a,a,0,is)=tempgq
      R2(g,a,a,1,is)=tempgq
      R2(g,a,a,2,is)=tempgq

      R1(g,a,g,0,is)=tempgq
      R1(g,a,g,1,is)=tempgq
      R1(g,a,g,2,is)=tempgq
      R2(g,a,g,0,is)=tempgq
      R2(g,a,g,1,is)=tempgq
      R2(g,a,g,2,is)=tempgq

      R1(g,a,q,0,is)=tempgq
      R1(g,a,q,1,is)=tempgq
      R1(g,a,q,2,is)=tempgq
      R2(g,a,q,0,is)=tempgq
      R2(g,a,q,1,is)=tempgq
      R2(g,a,q,2,is)=tempgq

      enddo

      endif 

************************************************************************
*     Contributions from QQQQ matrix elements                          *
************************************************************************            
      if (Qflag) then

c--- implement leading color by defining color factors in which 1/xn=0
      if (colourchoice .eq. 1) then
        xninv=0d0
      else
        xninv=1d0/xn   
      endif
      
c--- QUARK-QUARK contributions
      do is=1,3
      R1(q,q,q,0,is)=ason4pi*(
     . -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*xninv
     . -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*xninv
     . +ii_qq(z,xl12,is)*(xn+xninv)
     . +ff_qq(z,xl56,is)*(xn+xninv))
      R1(q,q,q,1,is)=ason4pi*(
     . -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*xninv
     . +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn-two*xninv)
     . +ii_qq(z,xl12,is)*two*xninv
     . +ff_qq(z,xl56,is)*two*xninv)
      R1(q,q,q,2,is)=ason4pi*(
     . +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn-two*xninv)
     . -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*xninv
     . +ii_qq(z,xl12,is)*two*xninv
     . +ff_qq(z,xl56,is)*two*xninv)
      R2(q,q,q,0,is)=ason4pi*(
     . -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*xninv
     . -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*xninv
     . +ii_qq(z,xl12,is)*(xn+xninv)
     . +ff_qq(z,xl56,is)*(xn+xninv))
      R2(q,q,q,1,is)=ason4pi*(
     . +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn-two*xninv)
     . -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*xninv
     . +ii_qq(z,xl12,is)*two*xninv
     . +ff_qq(z,xl56,is)*two*xninv)
      R2(q,q,q,2,is)=ason4pi*(
     . -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*xninv
     . +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn-two*xninv)
     . +ii_qq(z,xl12,is)*two*xninv
     . +ff_qq(z,xl56,is)*two*xninv)

      enddo


c--- ANTIQUARK-ANTIQUARK contributions
      do is=1,3
      R1(a,a,a,0,is)=R1(q,q,q,0,is)
      R1(a,a,a,1,is)=R1(q,q,q,1,is)
      R1(a,a,a,2,is)=R1(q,q,q,2,is)
      R2(a,a,a,0,is)=R2(q,q,q,0,is)
      R2(a,a,a,1,is)=R2(q,q,q,1,is)
      R2(a,a,a,2,is)=R2(q,q,q,2,is)
      
      enddo

c--- QUARK-ANTIQUARK contributions
      do is=1,3
      R1(q,q,a,0,is)=ason4pi*(
     . -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*xninv
     . +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn+xninv)
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R1(q,q,a,1,is)=ason4pi*(
     . +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn-two*xninv)
     . +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two*xninv
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R1(q,q,a,2,is)=ason4pi*(
     . -(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*xninv
     . +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*two*xninv
     . +ii_qq(z,xl12,is)*(xn-two*xninv)
     . +ff_qq(z,xl56,is)*(xn-two*xninv))
      R2(a,a,q,0,is)=ason4pi*(
     . +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn+xninv)
     . -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*xninv
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R2(a,a,q,1,is)=ason4pi*(
     . +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two*xninv
     . +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn-two*xninv)
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R2(a,a,q,2,is)=ason4pi*(
     . +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*two*xninv
     . -(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*xninv
     . +ii_qq(z,xl12,is)*(xn-two*xninv)
     . +ff_qq(z,xl56,is)*(xn-two*xninv))

      enddo
      
c--- ANTIQUARK-QUARK contributions
      do is=1,3
      R1(a,a,q,0,is)=ason4pi*(
     . +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*(xn+xninv)
     . -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*xninv
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R1(a,a,q,1,is)=ason4pi*(
     . +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two*xninv
     . +(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*(xn-two*xninv)
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R1(a,a,q,2,is)=ason4pi*(
     . +(if_qq(z,xl15,is)+fi_qq(z,xl15,is))*two*xninv
     . -(if_qq(z,xl16,is)+fi_qq(z,xl16,is))*xninv
     . +ii_qq(z,xl12,is)*(xn-two*xninv)
     . +ff_qq(z,xl56,is)*(xn-two*xninv))
      R2(q,q,a,0,is)=ason4pi*(
     . -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*xninv
     . +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*(xn+xninv)
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R2(q,q,a,1,is)=ason4pi*(
     . +(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*(xn-two*xninv)
     . +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two*xninv
     . -ii_qq(z,xl12,is)*xninv
     . -ff_qq(z,xl56,is)*xninv)
      R2(q,q,a,2,is)=ason4pi*(
     . -(if_qq(z,xl25,is)+fi_qq(z,xl25,is))*xninv
     . +(if_qq(z,xl26,is)+fi_qq(z,xl26,is))*two*xninv
     . +ii_qq(z,xl12,is)*(xn-two*xninv)
     . +ff_qq(z,xl56,is)*(xn-two*xninv))

      enddo

      do is=1,3

c--- GLUON-QUARK contributions
      tempqg=ason4pi*ii_qg(z,xl12,is)
      R1(q,g,q,0,is)=tempqg
      R1(q,g,q,1,is)=tempqg
      R1(q,g,q,2,is)=tempqg

      R1(a,g,q,0,is)=tempqg
      R1(a,g,q,1,is)=tempqg
      R1(a,g,q,2,is)=tempqg

c--- GLUON-ANTIQUARK contributions
      R1(q,g,a,0,is)=tempqg
      R1(q,g,a,1,is)=tempqg
      R1(q,g,a,2,is)=tempqg

      R1(a,g,a,0,is)=tempqg
      R1(a,g,a,1,is)=tempqg
      R1(a,g,a,2,is)=tempqg


c--- QUARK-GLUON contributions
      R2(q,g,q,0,is)=tempqg
      R2(q,g,q,1,is)=tempqg
      R2(q,g,q,2,is)=tempqg

      R2(a,g,q,0,is)=tempqg
      R2(a,g,q,1,is)=tempqg
      R2(a,g,q,2,is)=tempqg

c--- ANTIQUARK-GLUON contributions
      R2(q,g,a,0,is)=tempqg
      R2(q,g,a,1,is)=tempqg
      R2(q,g,a,2,is)=tempqg

      R2(a,g,a,0,is)=tempqg
      R2(a,g,a,1,is)=tempqg
      R2(a,g,a,2,is)=tempqg

      enddo

      endif

      return
      end

