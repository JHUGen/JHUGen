      subroutine qqbw2j_loop(i1,i2,i3,i4,i5,i6,
     .                       qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     .                       qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji)
      implicit none
      include 'constants.f'
      include 'zprods_com.f'
      integer i1,i2,i3,i4,i5,i6
      double complex a61_1234(2),a61_3412(2),a61_1432(2),a61_3214(2),
     .               a62_1234(2),a62_3412(2),a62_1432(2),a62_3214(2),
     .               atr_1234(2),atr_3412(2),atr_1432(2),atr_3214(2)
      double complex atrLLL,atrLRL,a61LLL,a61LRL,a62LLL,a62LRL
      double precision qqb_ijkk,qqb_iikl,qqb_ijkj,qqb_ijik,
     .                 qqb_ijii,qqb_ijjj,qqb_iiij,qqb_iiji
c--- this routine returns various pieces of the tree-loop interference
c--- for the process 0 -> q(1) qb(4) Q(3) Qb(2) + W (-> l(5) + lbar(6))
c--- the returned pieces are labelled according to the qqb -> qqb
c--- crossing and are calculated in the FORM program sq4q, with
c--- reference to the lowest order qqb_w2jetx.f

c--- there are 2 helicity amplitudes, we label LLL as 1, LRL as 2
      a61_1234(1)=a61LLL(i1,i2,i3,i4,i5,i6,za,zb)
      a61_1234(2)=a61LRL(i1,i2,i3,i4,i5,i6,za,zb)
 
      a61_3412(1)=a61LLL(i3,i4,i1,i2,i5,i6,za,zb)
      a61_3412(2)=a61LRL(i3,i4,i1,i2,i5,i6,za,zb)
 
      a61_1432(1)=a61LLL(i1,i4,i3,i2,i5,i6,za,zb)
      a61_1432(2)=a61LRL(i1,i4,i3,i2,i5,i6,za,zb)
 
      a61_3214(1)=a61LLL(i3,i2,i1,i4,i5,i6,za,zb)
      a61_3214(2)=a61LRL(i3,i2,i1,i4,i5,i6,za,zb)
 
      a62_1234(1)=a62LLL(i1,i2,i3,i4,i5,i6,za,zb)
      a62_1234(2)=a62LRL(i1,i2,i3,i4,i5,i6,za,zb)
 
      a62_3412(1)=a62LLL(i3,i4,i1,i2,i5,i6,za,zb)
      a62_3412(2)=a62LRL(i3,i4,i1,i2,i5,i6,za,zb)
 
      a62_1432(1)=a62LLL(i1,i4,i3,i2,i5,i6,za,zb)
      a62_1432(2)=a62LRL(i1,i4,i3,i2,i5,i6,za,zb)
 
      a62_3214(1)=a62LLL(i3,i2,i1,i4,i5,i6,za,zb)
      a62_3214(2)=a62LRL(i3,i2,i1,i4,i5,i6,za,zb)

      atr_1234(1)=Dconjg(atrLLL(i1,i2,i3,i4,i5,i6,za,zb))
      atr_1234(2)=Dconjg(atrLRL(i1,i2,i3,i4,i5,i6,za,zb))
 
      atr_3412(1)=Dconjg(atrLLL(i3,i4,i1,i2,i5,i6,za,zb))
      atr_3412(2)=Dconjg(atrLRL(i3,i4,i1,i2,i5,i6,za,zb))
 
      atr_1432(1)=Dconjg(atrLLL(i1,i4,i3,i2,i5,i6,za,zb))
      atr_1432(2)=Dconjg(atrLRL(i1,i4,i3,i2,i5,i6,za,zb))
 
      atr_3214(1)=Dconjg(atrLLL(i3,i2,i1,i4,i5,i6,za,zb))
      atr_3214(2)=Dconjg(atrLRL(i3,i2,i1,i4,i5,i6,za,zb))

* Note: in the expressions below, there are no LRL (2) interference
*       terms, since these would correspond to polarizations
*       that don't match, ie. (LL -> RR) x (LR -> LR)
      qqb_ijkk=dble(a61_1234(1)*atr_1234(1)+a61_1234(2)*atr_1234(2))
      qqb_iikl=dble(a61_3412(1)*atr_3412(1)+a61_3412(2)*atr_3412(2))
      qqb_ijkj=dble(a61_3214(1)*atr_3214(1)+a61_3214(2)*atr_3214(2))
      qqb_ijik=dble(a61_1432(1)*atr_1432(1)+a61_1432(2)*atr_1432(2))
      qqb_ijii=dble(a61_1234(1)*atr_1234(1)+a61_1234(2)*atr_1234(2)
     .             +a61_1432(1)*atr_1432(1)+a61_1432(2)*atr_1432(2)
     .            -(a62_1234(1)*atr_1432(1)+a62_1432(1)*atr_1234(1))/xn)
      qqb_ijjj=dble(a61_1234(1)*atr_1234(1)+a61_1234(2)*atr_1234(2)
     .             +a61_3214(1)*atr_3214(1)+a61_3214(2)*atr_3214(2)
     .            -(a62_1234(1)*atr_3214(1)+a62_3214(1)*atr_1234(1))/xn)
      qqb_iiij=dble(a61_1432(1)*atr_1432(1)+a61_1432(2)*atr_1432(2)
     .             +a61_3412(1)*atr_3412(1)+a61_3412(2)*atr_3412(2)
     .            -(a62_1432(1)*atr_3412(1)+a62_3412(1)*atr_1432(1))/xn)
      qqb_iiji=dble(a61_3214(1)*atr_3214(1)+a61_3214(2)*atr_3214(2)
     .             +a61_3412(1)*atr_3412(1)+a61_3412(2)*atr_3412(2)
     .            -(a62_3214(1)*atr_3412(1)+a62_3412(1)*atr_3214(1))/xn)
      
      return
      end
      
      
      
      
      
