      subroutine ALC_pp(k1,k2,k3,k4,k5,k6,coeff,ampALC)
c--- This is the leading colour amplitude in the notation
c---  (q1, Qb2, Q3, qb4)

c--- This routine just applies the appropriate interchange to get from
c---  the mm amplitude to the pp one, namely 1<->4, 2<->3, 5<->6

      implicit none
      include 'Wbbmlabels.f'
      integer k1,k2,k3,k4,k5,k6
      double complex ampALC(2,2)
      double complex coefftmp(0:4,20)
      
c--- zero out all integral coefficients      
      call clearcoeffs(coeff)
      
      call ALC_mm(k4,k3,k2,k1,k6,k5,coefftmp,ampALC)

      coeff(4,d2x3x4)=dconjg(coefftmp(4,d1x2x3))
      coeff(4,d1x2x3)=dconjg(coefftmp(4,d2x3x4))
      coeff(4,d1x23x4)=dconjg(coefftmp(4,d1x23x4))
      coeff(4,d1x2x34)=dconjg(coefftmp(4,d12x3x4))
      coeff(4,d12x3x4)=dconjg(coefftmp(4,d1x2x34))
      
      coeff(3,c23x4)=dconjg(coefftmp(3,c1x23))
      coeff(3,c1x23)=dconjg(coefftmp(3,c23x4))
      coeff(3,c2x3)=dconjg(coefftmp(3,c2x3))
      coeff(3,c2x34)=dconjg(coefftmp(3,c12x3))
      coeff(3,c3x4)=dconjg(coefftmp(3,c1x2))
      coeff(3,c12x3)=dconjg(coefftmp(3,c2x34))
      coeff(3,c1x2)=dconjg(coefftmp(3,c3x4))
      coeff(3,c1x234)=dconjg(coefftmp(3,c123x4))
      coeff(3,c123x4)=dconjg(coefftmp(3,c1x234))
      coeff(3,c12x34)=dconjg(coefftmp(3,c12x34))

      coeff(2,b123)=dconjg(coefftmp(2,b234))
      coeff(2,b234)=dconjg(coefftmp(2,b123))
      coeff(2,b23)=dconjg(coefftmp(2,b23))
      coeff(2,b1234)=dconjg(coefftmp(2,b1234))
      coeff(2,b12)=dconjg(coefftmp(2,b34))
      coeff(2,b34)=dconjg(coefftmp(2,b12))
      coeff(2,b2x1m)=dconjg(coefftmp(2,b2x1m))
      
      coeff(1,a0m)=dconjg(coefftmp(1,a0m))
      
      coeff(0,irat)=dconjg(coefftmp(0,irat))
      
      return
      end
      
