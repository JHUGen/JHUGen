      subroutine BBamps(p,ig,ia,ie,in,it,ib,BB)
      implicit none
      include 'constants.f'
      include 'masses.f'
      include 'zprods_decl.f'
      integer ig,ia,it,ib,ie,in,nu,j,jpart,i1,i2,i3,i4
      double precision p(mxpart,4),q(mxpart,4),tvec(4),bvec(4)
      double precision dot,mtsq,
     . aDg,aDb,aDt,bDg,nDt,bDe,bga2
      double complex BB(2,2,2,2),tga2,gDt
      double complex tzab(mxpart,mxpart),tzba(mxpart,mxpart)
      double complex bzab(mxpart,mxpart),bzba(mxpart,mxpart)
      mtsq=mt**2
      mbsq=mb**2
      aDb=dot(p,ia,ib)
      aDg=dot(p,ia,ig)
      aDt=dot(p,ia,it)
      bDg=dot(p,ib,ig)
      gDt=dcmplx(dot(p,ig,it),0d0)
      nDt=dot(p,in,it)
      bDe=dot(p,ib,ie)
      bga2=2d0*(bDg+aDb+aDg)
      tga2=2d0*gDt+dcmplx(2d0*(aDt+aDg),0d0)
      do j=1,4
      tvec(j)=p(it,j)
      bvec(j)=p(ib,j)
      enddo
      do nu=1,4
      do jpart=1,6
      q(jpart,nu)=p(jpart,nu)
      if (jpart.eq.it) q(it,nu)=p(it,nu)-0.5D0*mtsq/nDt*p(in,nu)
      if (jpart.eq.ib) q(ib,nu)=p(ib,nu)-0.5D0*mbsq/bDe*p(ie,nu)
      enddo
      enddo
      call spinoru(6,q,za,zb)
      call spinork(6,q,tzab,tzba,tvec)
      call spinork(6,q,bzab,bzba,bvec)

c--- Use the "non-overall scheme" to include the top width when necessary
c      if (dble(tga2+mt**2) .gt. 0d0) then
c        tga2=tga2+dcmplx(0d0,mt*twidth)
c        tga2=dcmplx(dsqrt(dble(tga2)**2+(mt*twidth)**2),0d0)
c      endif
c      if (dble(2d0*gDt+mt**2) .gt. 0d0) then
c        gDt=gDt+dcmplx(0d0,0.5d0*mt*twidth)
c      endif

      BB(1,1,1,1)= +one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & *mbsq*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie)*za(ie,
     &    ia)*zb(ig,in)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzab(ig,ig)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*tzab(ig,ig)*bzba(in,ie)*mtsq*gDt**(-1) * (  - za(ie
     &    ,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ig)*bzba(ia,ia)*tzba(in,ie)*
     & tga2*bga2**(-1)*aDb**(-1) * ( 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ig)*tzba(in,ie)*tga2*bga2**(-1)
     &  * (  - 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ia)*half*mtsq*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ia)*tzba(in,ig)*tzba(ia,ie)*
     & half*tga2*gDt**(-1)*aDb**(-1) * (  - 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(ia,ia)*tzba(in,ie)*tga2*bga2**(-1)
     & *aDb**(-1) * (  - za(ig,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ig,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*tzba(in,ie)*tga2*bga2**(-1) * ( za(ig,
     &    ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*mtsq * (  - za(ig,ie)*zb(ig,in) + za(ie
     &    ,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ie)*mtsq*gDt**(-1) * ( za(ig,ia
     &    )*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ie)*tzba(in,ig)*tzba(ia,ia)*
     & gDt**(-1) * ( 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ie)*tzba(in,ig) * ( 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ia)*tzba(ia,ie) * (  - 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ia)*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1)
     &  * ( 1.D0 )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1) * (  - za(
     &    ig,ia)*zb(ig,in) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*tzba(in,ig)*tzba(ia,ie)*mbsq*half*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(ig,in) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(in,it))*tzba(in,ie)*mbsq*tga2*bga2**(-1) * ( za(ig,ie)*zb(
     &    ig,in) - za(ie,ia)*zb(in,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(ig,ie)*bzba(in,ia)*tzba(in,ig)*half*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,ie) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ig)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*aDb**(-1) * ( za(
     &    ie,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*mtsq*gDt**(-1) * (  - za(ig,ia)*za(ie,ia)*zb(in,ia)
     &     )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*tzba(ig,ie)*tzba(in,ig)*gDt**(-1) * (  - za(ig,ia)
     &     )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*tzba(in,ig)*mtsq*gDt**(-1) * ( za(ie,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*tzba(in,ig)*tzba(ia,ia)*gDt**(-1) * (  - za(ie,ia)
     &     )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*tzba(in,ig) * (  - za(ie,ia) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(in,ie)*tzba(in,ia) * ( za(ig,ie) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *tzba(in,ig)*mbsq*half*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,ie)*
     &    za(ie,ia)*zb(ig,in) )
      BB(1,1,1,1) = BB(1,1,1,1)+one/(za(ib,ie))/(zb(in,it))*bzba(in,ie)
     & *tzba(in,ig)*gDt**(-1) * ( za(ig,ia)*za(ie,ia) )

      BB(1,1,1,2)= +one/(za(in,it))/(za(ib,ie))*bzba(in,ie)*mt*
     & gDt**(-1) * (  - za(ig,in)*za(ig,ia)*za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,in)*za(ig,ie
     &    )*za(ie,ia)*zb(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *tzab(in,ia)*bzba(in,ie)*mt*gDt**(-1) * (  - za(ig,ia)*za(ie,ia)
     &     )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(ig,ie)*bzba(in,ia)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - za(ig,in)*za(ig,ie) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(in,ig)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*
     &    za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(in,ie)*mt*mtsq*gDt**(-1) * (  - za(ig,in)*za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(in,ie)*mt * ( za(ig,in)*za(ie,ia) + za(ig,ie)*za(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(in,ie)*tzba(ig,ie)*mt*gDt**(-1) * ( za(ig,in)*za(ig,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(in,ie)*tzba(ia,ia)*mt*gDt**(-1) * ( za(ig,in)*za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*mbsq*mt*tga2*bga2**(-1) * (  - za(ig,ie)*za(ie,in)*
     &    zb(ig,in) + za(ie,in)*za(ie,ia)*zb(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzab(ig,ig)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &    za(ie,in)*za(ie,ia)*zb(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(ig,ig)*tzab(in,ia)*bzba(in,ie)*mt*gDt**(-1)
     &  * (  - za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(in,ig)*bzba(in,ie)*mt * ( za(ig,ie) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(in,ia)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ig,ie)*za(ie,ia)*zb(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(in,ia)*bzba(ig,ie)*bzba(in,ia)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * (  - za(ig,ie) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(in,ia)*bzba(in,ie)*mt * ( za(ie,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzab(in,ia)*bzba(in,ie)*tzba(ig,ie)*mt*gDt**(-1)
     &  * ( za(ig,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &     - za(ig,ia)*za(ie,in)*zb(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*mt*tga2*bga2**(-1) * (  - za(ig,ia)*za(
     &    ie,in)*zb(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(in,ig)*mt*tga2*bga2**(-1) * ( za(
     &    ie,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(in,ig)*bzba(ia,ia)*mt*tga2*
     & bga2**(-1)*aDb**(-1) * (  - za(ie,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(in,ia)*tzba(ia,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * ( za(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(ia,ia)*mt*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ig,ia)*za(ie,in)*zb(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(in,ie)*tzba(ig,ie)*mt * (  - za(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(in,ie)*tzba(ig,ie)*tzba(ia,ia)*mt*gDt**(-1)
     &  * (  - za(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(in,ie)*tzba(ia,ie)*mt * (  - za(in,ia) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(in,ia)*bzba(ia,ie)*mt*tga2*bga2**(-1) * (  - 
     &    za(ie,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*bzba(ia,ie)*mt*tga2*bga2**(-1) * ( za(ig,ia)*za(ie,
     &    in)*zb(ig,in) )
      BB(1,1,1,2) = BB(1,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(ig,ia))*tzba(ia,ie)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ig,in)*za(ie,ia)*zb(ig,in) )

      BB(1,1,2,1)= +one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *mb*mtsq * (  - za(ig,ie)*zb(ig,in)*zb(ie,in) + za(ie,ia)*zb(in,
     &    ia)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ig)*bzab(ia,ie)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * (  - zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(in,ie)*mb*tga2*bga2**(-1) * (  - 
     &    zb(ig,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ig,ie)*zb(ig,in)*zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(in,ig)*tzba(ia,ie)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * (  - zb(ig,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(in,ie)*mb*tga2*bga2**(-1) * (  - 
     &    zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzab(ig,ig)*mb*mtsq*gDt**(-1) * (  - za(ie,ia)*zb(
     &    in,ia)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*bzba(ia,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * (  - zb(ig,ie) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(in,ie)*mb*tga2*bga2**(-1) * ( zb(
     &    ig,ie) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ig,ie)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(in,ig)*tzba(ia,ie)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * ( zb(ig,ie) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(in,ie)*mb*tga2*bga2**(-1) * ( zb(
     &    ie,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(ia,ia)*tzba(in,ie)*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ig,ia)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*mb*mtsq*gDt**(-1) * ( za(ig,ia)*zb(in,
     &    ia)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(in,ig)*mb * ( zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(in,ig)*tzba(ia,ia)*mb*gDt**(-1)
     &  * ( zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ie)*mbsq*mb*tga2*bga2**(-1)*aDb**(-1) * ( 
     &     - za(ig,ia)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ie)*mb*tga2*bga2**(-1) * (  - za(ig,ia)*zb(
     &    ig,in)*zb(ie,ia) - za(ig,ia)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ia)*tzba(ia,ie)*mb * (  - zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *mb*mtsq*gDt**(-1) * (  - za(ig,ia)*za(ie,ia)*zb(in,ia)*zb(ie,in
     &    ) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ia,ie)*bzba(in,ig)*tzba(in,ie)*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( 1.D0 )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ia,ie)*tzba(in,ig)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    za(ig,ie)*zb(ig,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(in,ia)*tzba(in,ig)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - za(ig,ie)*zb(ig,ie) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(ig,ie)*tzba(in,ig)*mb*gDt**(-1) * (  - za(ig,ia)*zb(ie,in)
     &     )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ig)*mb*mtsq*gDt**(-1) * ( za(ie,ia)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ig)*mb * (  - za(ie,ia)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ig)*tzba(ia,ia)*mb*gDt**(-1) * (  - za(ie,ia)*zb(ie,in)
     &     )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ia)*mb * ( za(ig,ie)*zb(ie,in) )
      BB(1,1,2,1) = BB(1,1,2,1)+one/(zb(in,it))/(zb(ib,ie))*tzba(in,ig)
     & *mb*gDt**(-1) * ( za(ig,ia)*za(ie,ia)*zb(ie,in) )

      BB(1,1,2,2)= +one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))/(zb(ib,ie))
     & *mbsq*mt*mb*tga2*bga2**(-1)*aDb**(-1) * ( za(ig,ia)*za(ie,in)*
     &    zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*mt*mb*tga2*bga2**(-1) * ( za(ig,ia)*za(ie,in)*zb(ig
     &    ,in)*zb(ie,ia) + za(ig,ia)*za(ie,in)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ig,ig)*bzab(ia,ie)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ie,in)*zb(in,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ig,ie)*mt*mb*tga2*bga2**(-1) * ( za(ie,in)*zb(
     &    ig,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ia,ie)*mt*mb*tga2*bga2**(-1) * ( za(ie,in)*zb(
     &    in,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ia,ie)*tzab(in,ia)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(ig,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(ia,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,in)*zb(ig,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(ig,ig)*tzab(in,ia)*mt*mb*gDt**(-1) * (  - za(
     &    ie,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ig)*mt*mb * ( za(ig,ie)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ia)*mt*mb * ( za(ie,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ia)*bzba(in,ia)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,ie)*zb(ig,ie) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ia)*tzba(ig,ie)*mt*mb*gDt**(-1) * ( za(ig,
     &    ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ig)*mt*mb*tga2*bga2**(-1) * (  - za(ie,in)*
     &    zb(ig,ie) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ig)*bzba(ia,ia)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ie,in)*zb(ig,ie) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ia)*mt*mb*tga2*bga2**(-1) * (  - za(ie,in)*
     &    zb(ie,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(ia,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,in)*zb(ig,ie) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(ia,ia)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - 
     &    za(ig,ia)*za(ie,in)*zb(ig,ie)*zb(in,ia) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzba(ig,ie)*mt*mb * (  - za(ig,in)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(ia,ia)*mt*mb*gDt**(-1) * (  - za(
     &    ig,in)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ig,ia))
     & /(zb(ib,ie))*tzba(ia,ie)*mt*mb * (  - za(in,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *mt*mb*mtsq*gDt**(-1) * (  - za(ig,in)*za(ie,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *mt*mb * ( za(ig,in)*za(ie,ia)*zb(ie,in) + za(ig,ie)*za(in,ia)*
     &    zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzab(ia,ie)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,in
     &    )*za(ig,ie)*zb(ig,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzab(ia,ie)*bzba(in,ig)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - 
     &    za(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *tzab(in,ia)*mt*mb*gDt**(-1) * (  - za(ig,ia)*za(ie,ia)*zb(ie,in
     &    ) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzba(in,ia)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,in)*
     &    za(ig,ie)*zb(ig,ie) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *tzba(ig,ie)*mt*mb*gDt**(-1) * ( za(ig,in)*za(ig,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *tzba(ia,ia)*mt*mb*gDt**(-1) * ( za(ig,in)*za(ie,ia)*zb(ie,in) )
      BB(1,1,2,2) = BB(1,1,2,2)+one/(za(in,it))/(zb(ib,ie))*mt*mb*
     & gDt**(-1) * (  - za(ig,in)*za(ig,ia)*za(ie,ia)*zb(ie,in) )

      BB(1,2,1,1)= +one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *mbsq*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * ( za(ie,ia)**2*zb(ig,
     &    in)**2 )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzab(ia,ig)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(ig,in) )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*tzab(ia,ig)*bzba(in,ie)*mtsq*gDt**(-1) * (  - za(ie
     &    ,ia)*zb(ig,in) )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(ig,ia)*bzba(in,ia)*tzba(in,ie)*
     & tga2*bga2**(-1)*aDb**(-1) * (  - 1.D0 )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ia)*half*mtsq*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ie,ia)*zb(ig,in) )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ia)*tzba(ig,ie)*tzba(in,ia)*
     & half*tga2*gDt**(-1)*aDb**(-1) * ( 1.D0 )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ie)*tzba(ig,ia)*tzba(in,ia)*
     & gDt**(-1) * (  - 1.D0 )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*tzba(ig,ie)*tzba(in,ia)*mbsq*half*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ie,ia)*zb(ig,in) )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(in,ie)*tzba(ig,ia)*tzba(in,ia)*gDt**(-1) * ( za(ie,ia) )
      BB(1,2,1,1) = BB(1,2,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(ig,ie)*bzba(ig,ia)*tzba(in,ie)*tga2*bga2**(-1)*aDb**(-1)
     &  * ( zb(ig,in) )

      BB(1,2,1,2)= +one/(za(ig,ia))/(za(in,it))/(za(ib,ie))*bzba(in,ie)
     & *tzba(ig,ia)*mt*gDt**(-1) * ( za(in,ia)*za(ie,ia) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzab(ia,ig)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &    za(ie,in)*za(ie,ia)*zb(ig,in) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(in,ig)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ie,ia)**2*zb(ig,in) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(in,ig)*tzab(ia,ig)*bzba(in,ie)*mt*gDt**(-1)
     &  * ( za(ie,ia) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(in,ig)*bzba(ig,ie)*bzba(in,ia)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * (  - za(ie,ia) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(ig,ia)*bzba(in,ia)*mt*tga2*
     & bga2**(-1)*aDb**(-1) * ( za(ie,in) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(ig,ie)*bzba(in,ia)*tzba(ig,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * ( za(in,ia) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(in,ie)*tzba(ig,ie)*tzba(ig,ia)*mt*gDt**(-1)
     &  * (  - za(in,ia) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzba(ig,ie)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(in,ia)*za(ie,ia)*zb(ig,in) )
      BB(1,2,1,2) = BB(1,2,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(ig,ie)*bzba(ig,ia)*mt*tga2*bga2**(-1)*aDb**(-1) * (  - za(
     &    ie,in)*zb(ig,in) )

      BB(1,2,2,1)= +one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ia,ig)*bzab(ia,ie)*tzba(in,ie)*mb*tga2*bga2**(-1)*
     & aDb**(-1) * (  - zb(ig,in) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ie,ia)*zb(ig,in)**2 )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(ig,ie)*tzba(in,ia)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * ( zb(ig,in) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzab(ia,ig)*mb*mtsq*gDt**(-1) * (  - za(ie,ia)*zb(
     &    ig,in)*zb(ie,in) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(ig,ia)*bzba(in,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * ( zb(ig,ie) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ie,ia)*zb(ig,in)*zb(ig,ie) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(ig,ie)*tzba(in,ia)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * (  - zb(ig,ie) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(ig,ia)*tzba(in,ia)*mb*gDt**(-1)
     &  * (  - zb(ie,in) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(ig,ia)*tzba(in,ia)*mb*gDt**(-1) * ( za(ie,ia)*zb(ie,in) )
      BB(1,2,2,1) = BB(1,2,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(ig,ia)*tzba(in,ie)*mb*tga2*bga2**(-1)*aDb**(-1) * (  - zb(
     &    ig,in)*zb(ig,ie) )

      BB(1,2,2,2)= +one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzab(ia,ig)*bzab(ia,ie)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * ( za(
     &    ie,in)*zb(ig,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ia,ie)*tzab(in,ig)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(ig,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(ig,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(in,ia)*zb(ig,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ig)*tzab(ia,ig)*mt*mb*gDt**(-1) * ( za(ie,
     &    ia)*zb(ie,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ig)*bzba(in,ia)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ie,ia)*zb(ig,ie) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(ig,ia)*bzba(in,ia)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ie,in)*zb(ig,ie) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(ig,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(in,ia)*zb(ig,ie) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(ig,ia)*mt*mb*gDt**(-1) * (  - za(
     &    in,ia)*zb(ie,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *tzba(ig,ia)*mt*mb*gDt**(-1) * ( za(in,ia)*za(ie,ia)*zb(ie,in) )
      BB(1,2,2,2) = BB(1,2,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzba(ig,ia)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * ( za(ie,in)*zb(ig
     &    ,in)*zb(ig,ie) )

      BB(2,1,1,1)= +one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *mbsq*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,ie)**2*zb(in,
     &    ia)**2 )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzab(ig,ia)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*tzab(ig,ia)*bzba(in,ie)*mtsq*gDt**(-1) * (  - za(ig
     &    ,ie)*zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ig)*bzba(ia,ig)*bzba(ia,ie)*tzba(in,ie)*
     & tga2*bga2**(-1)*aDb**(-1) * (  - 1.D0 )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ig)*bzba(ia,ie)*half*mtsq*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,ie)*zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ig)*bzba(ia,ie)*tzba(in,ig)*tzba(ia,ie)*
     & half*tga2*gDt**(-1)*aDb**(-1) * ( 1.D0 )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ig)*tzba(ia,ig)*tzba(ia,ie)*
     & gDt**(-1) * (  - 1.D0 )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(ig,ia))
     & /(zb(in,it))*tzba(in,ig)*tzba(ia,ie)*mbsq*half*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,ie)*zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(in,ig)*bzba(ia,ie)*tzba(in,ig)*half*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,ie) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(in,ie)*tzba(in,ig)*tzba(ia,ig)*gDt**(-1) * ( za(ig,ie) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *tzba(in,ig)*mbsq*half*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie)
     &    **2*zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(ia,ig)*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1)*aDb**(-1)
     &  * ( zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(ia,ie)*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie)
     &    *zb(in,ia)**2 )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ib,ie))/(zb(ig,ia))/(zb(in,it))
     & *bzba(ia,ie)*tzba(in,ig)*tzba(ia,ie)*half*tga2*gDt**(-1)*
     & aDb**(-1) * (  - zb(in,ia) )
      BB(2,1,1,1) = BB(2,1,1,1)+one/(za(ib,ie))/(zb(in,it))*bzba(ia,ie)
     & *tzba(in,ig)*half*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,ie)*zb(in,
     &    ia) )

      BB(2,1,1,2)= +one/(za(ig,ia))/(za(in,it))/(za(ib,ie))*mbsq*half*
     & mt*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,in)*za(ig,ie)**2*zb(in,ia)
     &     )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(in,ig)*bzba(ia,ie)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    za(ig,in)*za(ig,ie) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(in,ie)*tzba(ia,ig)*mt*gDt**(-1) * (  - za(ig,in)*za(ig,ie)
     &     )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzab(ig,ia)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &    za(ig,ie)*za(ie,in)*zb(in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(ig,ia)*tzab(in,ia)*bzba(in,ie)*mt*gDt**(-1)
     &  * (  - za(ig,ie) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(in,ia)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ig,ie)**2*zb(in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzab(in,ia)*bzba(in,ig)*bzba(ia,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * ( za(ig,ie) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(in,ig)*bzba(ia,ig)*bzba(ia,ie)*mt*tga2*
     & bga2**(-1)*aDb**(-1) * ( za(ie,in) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(in,ig)*bzba(ia,ie)*tzba(ia,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * (  - za(ig,in) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*bzba(in,ie)*tzba(ia,ig)*tzba(ia,ie)*mt*gDt**(-1)
     &  * ( za(ig,in) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & /(zb(ig,ia))*tzba(ia,ie)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ig,in)*za(ig,ie)*zb(in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(in,it))/(za(ib,ie))*bzba(ia,ie)
     & *half*mt*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,in)*za(ig,ie)*zb(
     &    in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *tzab(in,ia)*bzba(ia,ie)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - za(ig,ie)*zb(in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(ia,ig)*bzba(ia,ie)*mt*tga2*bga2**(-1)*aDb**(-1) * (  - za(
     &    ie,in)*zb(in,ia) )
      BB(2,1,1,2) = BB(2,1,1,2)+one/(za(in,it))/(za(ib,ie))/(zb(ig,ia))
     & *bzba(ia,ie)*tzba(ia,ie)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    za(ig,in)*zb(in,ia) )

      BB(2,1,2,1)= +one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ig,ie)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,
     &    ie)*zb(in,ia)**2 )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*bzab(ig,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * ( zb(in,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(in,ig)*tzba(ia,ie)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * (  - zb(in,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzab(ig,ia)*mb*mtsq*gDt**(-1) * (  - za(ig,ie)*zb(
     &    in,ia)*zb(ie,in) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ig,ie)*zb(in,ia)*zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*bzba(ia,ig)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * (  - zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(in,ig)*tzba(ia,ie)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * ( zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ig)*tzba(ia,ig)*tzba(ia,ie)*mb*gDt**(-1)
     &  * (  - zb(ie,in) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ig,ie)*tzba(in,ig)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    za(ig,ie)*zb(in,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(in,ig)*tzba(in,ig)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - za(ig,ie)*zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ig)*tzba(ia,ig)*mb*gDt**(-1) * ( za(ig,ie)*zb(ie,in) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie)*zb(in,ia
     &    )**2*zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(ia,ig)*tzba(in,ie)*mb*tga2*bga2**(-1)*aDb**(-1) * ( zb(in,
     &    ia)*zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(zb(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ig)*tzba(ia,ie)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - zb(in,ia)*zb(ie,ia) )
      BB(2,1,2,1) = BB(2,1,2,1)+one/(zb(in,it))/(zb(ib,ie))*tzba(in,ig)
     & *half*mb*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,ie)*zb(in,ia)*zb(ie,
     &    ia) )

      BB(2,1,2,2)= +one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzab(ig,ie)*bzab(ig,ia)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - 
     &    za(ie,in)*zb(in,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ig,ie)*tzab(in,ia)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(in,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(ia,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,in)*zb(in,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(ig,ia)*tzab(in,ia)*mt*mb*gDt**(-1) * (  - za(
     &    ig,ie)*zb(ie,in) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzab(in,ia)*bzba(in,ig)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ig,ie)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ig)*bzba(ia,ig)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ie,in)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(ia,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,in)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ig,ia))
     & /(zb(ib,ie))*tzba(ia,ig)*tzba(ia,ie)*mt*mb*gDt**(-1) * ( za(ig,
     &    in)*zb(ie,in) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *bzab(ig,ie)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,in
     &    )*za(ig,ie)*zb(in,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *bzba(in,ig)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,in)*
     &    za(ig,ie)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *tzba(ia,ig)*mt*mb*gDt**(-1) * (  - za(ig,in)*za(ig,ie)*zb(ie,in
     &    ) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *tzab(in,ia)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie
     &    )*zb(in,ia)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *bzba(ia,ig)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*
     &    zb(in,ia)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(in,it))/(zb(ig,ia))/(zb(ib,ie))
     & *tzba(ia,ie)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * ( za(ig,in)*
     &    zb(in,ia)*zb(ie,ia) )
      BB(2,1,2,2) = BB(2,1,2,2)+one/(za(in,it))/(zb(ib,ie))*half*mt*mb*
     & tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,in)*za(ig,ie)*zb(in,ia)*
     &    zb(ie,ia) )

      BB(2,2,1,1)= +one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *mbsq*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * (  - za(ig,ie)*za(ie,
     &    ia)*zb(ig,in)*zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzab(ia,ia)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*tzab(ia,ia)*bzba(in,ie)*mtsq*gDt**(-1) * (  - za(ig
     &    ,ie)*zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(ig,ig)*bzba(in,ia)*bzba(ia,ie)*tzba(in,ie)*
     & tga2*bga2**(-1)*aDb**(-1) * ( 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(ig,ie)*bzba(in,ig)*tzba(in,ie)*tga2*bga2**(-1)
     &  * (  - 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ig)*bzba(ia,ie)*half*mtsq*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ig)*bzba(ia,ie)*tzba(ig,ie)*tzba(in,ia)*
     & half*tga2*gDt**(-1)*aDb**(-1) * (  - 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*mtsq * (  - za(ig,ie)*zb(ig,in) + za(ie
     &    ,ia)*zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ig)*tzba(in,ia)*gDt**(-1) * ( 
     &     - za(ig,ie)*zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ig)*tzba(in,ia)*tzba(ia,ie)*
     & gDt**(-1) * ( 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(ig,ie)*tzba(in,ig) * ( 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ig) * (  - za(ie,ia)*zb(ig,ia)
     &     )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ia)*mtsq*gDt**(-1) * ( za(ig,ie
     &    )*zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ia)*tzba(ia,ie) * (  - 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ie)*tzba(in,ia) * ( za(ig,ie)*zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ia)*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1)
     &  * ( 1.D0 )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*bzba(in,ia)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ig,ie)*zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*tzba(ig,ie)*tzba(in,ia)*mbsq*half*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ig,ie)*zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ig,ia))/(za(ib,ie))
     & /(zb(in,it))*tzba(in,ie)*mbsq*tga2*bga2**(-1) * ( za(ig,ie)*zb(
     &    ig,in) - za(ie,ia)*zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ig,ig)*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1)*aDb**(-1)
     &  * (  - zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ig,ie)*tzba(in,ie)*tga2*bga2**(-1) * ( zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(in,ie)*tzba(ia,ie)*mtsq*gDt**(-1) * ( zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(in,ia)*bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1)*aDb**(-1)
     &  * (  - zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ia,ie)*half*mtsq*tga2*gDt**(-1)*aDb**(-1) * ( za(ie,ia)*
     &    zb(ig,in)*zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ia,ie)*tzba(ig,ie)*tzba(in,ia)*half*tga2*gDt**(-1)*
     & aDb**(-1) * ( zb(in,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ia,ie)*tzba(in,ie)*mbsq*tga2*bga2**(-1)*aDb**(-1) * ( zb(
     &    ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *bzba(ia,ie)*tzba(in,ie)*tga2*bga2**(-1) * (  - zb(ig,in) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ig,ia))/(za(ib,ie))/(zb(in,it))
     & *tzba(in,ie)*mbsq*tga2*bga2**(-1)*aDb**(-1) * (  - za(ig,ie)*zb(
     &    ig,in)*zb(ig,ia) )
      BB(2,2,1,1) = BB(2,2,1,1)+one/(za(ib,ie))/(zb(in,it))*bzba(ia,ie)
     & *tzba(in,ie)*tga2*bga2**(-1)*aDb**(-1) * ( zb(ig,in)*zb(ig,ia) )

      BB(2,2,1,2)= +one/(za(ig,ia))/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *mbsq*mt*tga2*bga2**(-1) * (  - za(ig,ie)*za(ie,in)*zb(ig,in) + 
     &    za(ie,in)*za(ie,ia)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzab(ia,ia)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &    za(ig,ie)*za(ie,in)*zb(ig,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzab(in,ig)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ig,ie)*za(ie,ia)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzab(in,ig)*tzab(ia,ia)*bzba(in,ie)*mt*gDt**(-1)
     &  * ( za(ig,ie) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzab(in,ig)*bzba(in,ig)*bzba(ia,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * ( za(ie,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzab(in,ig)*bzba(in,ie)*mt * ( za(ig,ie) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzab(in,ia)*bzba(in,ie)*mt * ( za(ie,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(ig,ig)*bzba(in,ia)*bzba(ia,ie)*mt*tga2*
     & bga2**(-1)*aDb**(-1) * (  - za(ie,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(ig,ie)*bzba(in,ig)*mt*tga2*bga2**(-1) * ( za(
     &    ie,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ig)*bzba(ia,ie)*tzba(ig,ie)*half*mt*tga2*
     & gDt**(-1)*aDb**(-1) * (  - za(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*mt*mtsq*gDt**(-1) * ( za(ig,ie)*za(in,
     &    ia)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*mt * ( za(ig,in)*za(ie,ia)*zb(ig,ia) + 
     &    za(ig,ie)*za(in,ia)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*tzba(ig,ig)*mt*gDt**(-1) * (  - za(ig,
     &    ie)*za(in,ia)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*tzba(ig,ig)*tzba(ia,ie)*mt*gDt**(-1)
     &  * ( za(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*tzba(ig,ie)*mt * (  - za(ig,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ie)*tzba(ia,ie)*mt * (  - za(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ia)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( 
     &     - za(ig,ie)*za(ie,in)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*bzba(in,ia)*bzba(ia,ie)*mt*tga2*bga2**(-1) * (  - 
     &    za(ie,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(za(ib,ie))*tzba(ig,ie)*mbsq*half*mt*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ig,ie)*za(in,ia)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * ( za(ig,ie)*za(ie,in)*zb(ig
     &    ,in)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *tzab(in,ig)*bzba(in,ie)*tzba(ia,ie)*mt*gDt**(-1) * (  - 1.D0 )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *tzab(in,ig)*bzba(ia,ie)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &     - za(ie,ia)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(ig,ig)*bzba(ia,ie)*mt*tga2*bga2**(-1)*aDb**(-1) * ( za(ie,
     &    in)*zb(ig,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(ig,ie)*mt*tga2*bga2**(-1) * (  - za(ie,in)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(in,ia)*bzba(ia,ie)*mt*tga2*bga2**(-1)*aDb**(-1) * ( za(ie,
     &    in)*zb(ig,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(ia,ie)*mbsq*mt*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*
     &    zb(ig,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(ia,ie)*mt*tga2*bga2**(-1) * ( za(ie,in)*zb(ig,in) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(ig,ia))/(za(in,it))/(za(ib,ie))
     & *bzba(ia,ie)*tzba(ig,ie)*half*mt*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    za(in,ia)*zb(in,ia) )
      BB(2,2,1,2) = BB(2,2,1,2)+one/(za(in,it))/(za(ib,ie))*bzba(ia,ie)
     & *mt*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*zb(ig,in)*zb(ig,
     &    ia) )

      BB(2,2,2,1)= +one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *mb*mtsq * (  - za(ig,ie)*zb(ig,in)*zb(ie,in) + za(ie,ia)*zb(in,
     &    ia)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * ( za(ie,ia)*zb(ig,in)*zb(in,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*bzab(ia,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * ( zb(ig,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*bzba(in,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * (  - zb(ig,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(ig,ie)*tzba(in,ia)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * ( zb(in,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(in,ie)*mb*tga2*bga2**(-1) * (  - 
     &    zb(ig,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*tzba(in,ie)*mb*tga2*bga2**(-1) * (  - 
     &    zb(in,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzab(ia,ia)*mb*mtsq*gDt**(-1) * (  - za(ig,ie)*zb(
     &    ig,in)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(ig,ig)*bzba(in,ia)*tzba(in,ie)*mb*tga2*
     & bga2**(-1)*aDb**(-1) * ( zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1)
     &  * (  - za(ie,ia)*zb(ig,in)*zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(ig,ie)*tzba(in,ia)*half*mb*tga2*
     & gDt**(-1)*aDb**(-1) * (  - zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(in,ie)*mb*tga2*bga2**(-1) * ( zb(
     &    ig,ie) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*tzba(in,ie)*mb*tga2*bga2**(-1) * ( zb(
     &    ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ig)*tzba(in,ia)*mb*gDt**(-1) * (  - za(ig,
     &    ie)*zb(ig,ia)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ig)*tzba(in,ia)*tzba(ia,ie)*mb*gDt**(-1)
     &  * ( zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*tzba(in,ig)*mb * ( zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ig)*mb * (  - za(ie,ia)*zb(ig,ia)*zb(ie,in)
     &     )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ia)*mb*mtsq*gDt**(-1) * ( za(ig,ie)*zb(ig,
     &    ia)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ia)*mb * ( za(ig,ie)*zb(ig,ia)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(za(ig,ia))/(zb(in,it))
     & /(zb(ib,ie))*tzba(in,ia)*tzba(ia,ie)*mb * (  - zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *half*mb*mtsq*tga2*gDt**(-1)*aDb**(-1) * ( za(ie,ia)*zb(ig,in)*
     &    zb(in,ia)*zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzab(ig,ie)*tzba(in,ie)*mb*tga2*bga2**(-1)*aDb**(-1) * ( zb(ig,
     &    in)*zb(ig,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(ig,ig)*tzba(in,ie)*mb*tga2*bga2**(-1)*aDb**(-1) * (  - zb(
     &    ig,in)*zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *bzba(in,ia)*tzba(in,ie)*mb*tga2*bga2**(-1)*aDb**(-1) * (  - zb(
     &    ig,ia)*zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(ig,ie)*tzba(in,ia)*half*mb*tga2*gDt**(-1)*aDb**(-1) * ( 
     &    zb(in,ia)*zb(ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ie)*mbsq*mb*tga2*bga2**(-1)*aDb**(-1) * ( zb(ig,in)*zb(
     &    ie,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(in,ie)*mb*tga2*bga2**(-1) * (  - zb(ig,in)*zb(ie,ia) - zb(
     &    ig,ie)*zb(in,ia) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(za(ig,ia))/(zb(in,it))/(zb(ib,ie))
     & *tzba(ia,ie)*mb*mtsq*gDt**(-1) * ( zb(ig,in)*zb(ie,in) )
      BB(2,2,2,1) = BB(2,2,2,1)+one/(zb(in,it))/(zb(ib,ie))*tzba(in,ie)
     & *mb*tga2*bga2**(-1)*aDb**(-1) * ( zb(ig,in)*zb(ig,ia)*zb(ie,ia)
     &     )

      BB(2,2,2,2)= +one/(za(ig,ia))/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *mt*mb*mtsq*gDt**(-1) * ( za(ig,ie)*za(in,ia)*zb(ig,ia)*zb(ie,in
     &    ) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*mt*mb * ( za(ig,in)*za(ie,ia)*zb(ig,ia)*zb(ie,in)
     &     + za(ig,ie)*za(in,ia)*zb(ig,ia)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*mt*mb*tga2*bga2**(-1) * ( za(ie,in)*zb(
     &    ig,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*bzab(ia,ia)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ie,in)*zb(ig,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzab(in,ig)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(ie,ia)*zb(in,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*bzba(in,ia)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * ( za(ie,in)*zb(ig,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ig,ie)*tzba(ig,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(in,ia)*zb(in,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzab(ia,ie)*mt*mb*tga2*bga2**(-1) * ( za(ie,in)*zb(
     &    in,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzab(in,ig)*mt*mb * ( za(ig,ie)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzab(in,ig)*tzab(ia,ia)*mt*mb*gDt**(-1) * ( za(ig,
     &    ie)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzab(in,ig)*bzba(in,ig)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * ( za(ie,ia)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzab(in,ia)*mt*mb * ( za(ie,ia)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzba(ig,ig)*bzba(in,ia)*mt*mb*tga2*bga2**(-1)*
     & aDb**(-1) * (  - za(ie,in)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*mt*mb*tga2*bga2**(-1) * (  - za(ie,in)*
     &    zb(ig,ie) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzba(in,ig)*tzba(ig,ie)*half*mt*mb*tga2*gDt**(-1)*
     & aDb**(-1) * (  - za(in,ia)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*bzba(in,ia)*mt*mb*tga2*bga2**(-1) * (  - za(ie,in)*
     &    zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzba(ig,ig)*mt*mb*gDt**(-1) * (  - za(ig,ie)*za(in,
     &    ia)*zb(ig,ia)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzba(ig,ig)*tzba(ia,ie)*mt*mb*gDt**(-1) * ( za(in,
     &    ia)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzba(ig,ie)*mt*mb * (  - za(ig,in)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(ig,ia))/(za(in,it))
     & /(zb(ib,ie))*tzba(ia,ie)*mt*mb * (  - za(in,ia)*zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *mbsq*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*zb(ig,in)
     &    *zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *mt*mb*tga2*bga2**(-1) * ( za(ie,in)*zb(ig,in)*zb(ie,ia) + za(ie
     &    ,in)*zb(ig,ie)*zb(in,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *bzab(ig,ie)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * (  - za(ie,in)*
     &    zb(ig,in)*zb(ig,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *tzab(in,ig)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * (  - za(ie,ia
     &    )*zb(in,ia)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *tzab(in,ig)*tzba(ia,ie)*mt*mb*gDt**(-1) * (  - zb(ie,in) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *bzba(ig,ig)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * ( za(ie,in)*zb(ig
     &    ,in)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *bzba(in,ia)*mt*mb*tga2*bga2**(-1)*aDb**(-1) * ( za(ie,in)*zb(ig
     &    ,ia)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(ig,ia))/(za(in,it))/(zb(ib,ie))
     & *tzba(ig,ie)*half*mt*mb*tga2*gDt**(-1)*aDb**(-1) * ( za(in,ia)*
     &    zb(in,ia)*zb(ie,ia) )
      BB(2,2,2,2) = BB(2,2,2,2)+one/(za(in,it))/(zb(ib,ie))*mt*mb*tga2*
     & bga2**(-1)*aDb**(-1) * (  - za(ie,in)*zb(ig,in)*zb(ig,ia)*zb(ie,
     &    ia) )

c--- Use the "overall scheme" to include the top width when necessary
      if (dble(tga2)+mt**2 .gt. 0d0) then
        tga2=tga2+dcmplx(0d0,mt*twidth)
      endif

      do i1=1,2
      do i2=1,2
      do i3=1,2
      do i4=1,2
      BB(i1,i2,i3,i4)=BB(i1,i2,i3,i4)/tga2
      enddo
      enddo
      enddo
      enddo
     
      return
      end
