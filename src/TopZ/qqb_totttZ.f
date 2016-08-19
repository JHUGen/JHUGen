      subroutine qqb_totttZ(p,msq)
      implicit none
************************************************************************
*     Author: R.K. Ellis                                               *
*     May, 2013.                                                       *
*     calculate the element squared                                    *
*     for the process                                                  *
c----My notation                                                       *
C     q(-p1) +qbar(-p2)=Z(e-(p3)+e+(p4))+t(p5)+t(p6)                   *
C                                                                      *
************************************************************************
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'msq_cs.f'
      
      integer j,jud,etatb,etat
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),
     & wtgg1,wtgg2,wtgg0,facqq,facgg,
     & wtuub,wtubu,wtddb,wtdbd

C-----vector for demassifying t~ and t
      etatb=1
      etat=1

C----down quarks
      jud=1
      call ttqqZampsq(p,jud,6,5,2,1,3,4,wtddb)
      call ttqqZampsq(p,jud,6,5,1,2,3,4,wtdbd)
C----up quarks
      jud=2
      call ttqqZampsq(p,jud,6,5,2,1,3,4,wtuub)
      call ttqqZampsq(p,jud,6,5,1,2,3,4,wtubu)

      call ttggZdriver(p,etatb,etat,wtgg1,wtgg2,wtgg0)


      facqq=V/4d0*gsq**2*esq**2
      facgg=facqq*xn

C----set all elements to zero
      msq(:,:)=0d0

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if ((j .eq. 1).or.(j .eq. 3).or.(j .eq. 5)) then
          msq(j,-j)=aveqq*facqq*wtddb
      elseif ((j .eq. -1).or.(j .eq. -3).or.(j .eq. -5)) then
          msq(j,-j)=aveqq*facqq*wtdbd
      elseif ((j .eq. 2).or.(j .eq. 4)) then
          msq(j,-j)=aveqq*facqq*wtuub
      elseif ((j .eq. -2).or.(j .eq. -4)) then
          msq(j,-j)=aveqq*facqq*wtubu
      elseif (j .eq. 0) then
          msq(j,j)=avegg*facgg*(wtgg1+wtgg2+wtgg0)
          msq_cs(1,j,j)=avegg*facgg*wtgg1
          msq_cs(2,j,j)=avegg*facgg*wtgg2
          msq_cs(0,j,j)=avegg*facgg*wtgg0
      endif
      enddo
      return
      end
 
