      subroutine qq4lggampf(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,msq)
      implicit none
      include 'constants.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zprods_decl.f'
      include 'interference.f'
      include 'pid_pdg.f'
      integer i1,i2,i3,i4,i5,i6,i7,i8,
     & jdu,h1,h3,h5,h7,h8
      double precision msq(2),colfac34_56
      double complex
     & a78(2,2,2,2,2,2),a87(2,2,2,2,2,2),aq(2,2,2,2,2,2),
     & amp78xy(2,2,2,2,2,2),amp87xy(2,2,2,2,2,2),
     & amp78yx(2,2,2,2,2,2),amp87yx(2,2,2,2,2,2),
     & tmp78(2,2,2,2,2,2),tmp87(2,2,2,2,2,2),
     & a78_swap(2,2,2,2,2,2),a87_swap(2,2,2,2,2,2),aq_swap(2,2,2,2,2,2),
     & amp78xy_swap(2,2,2,2,2,2),amp87xy_swap(2,2,2,2,2,2),
     & amp78yx_swap(2,2,2,2,2,2),amp87yx_swap(2,2,2,2,2,2),
     & tmp78_swap(2,2,2,2,2,2),tmp87_swap(2,2,2,2,2,2)

c---color factors for Z decays
      colfac34_56=1d0
      if (abs(pid_pdg(i3)).ge.0 .and. abs(pid_pdg(i3)).le.5) then
        colfac34_56=colfac34_56*xn
      endif
      if (abs(pid_pdg(i5)).ge.0 .and. abs(pid_pdg(i5)).le.5) then
        colfac34_56=colfac34_56*xn
      endif

      call qq4lggamp(i1,i2,i3,i4,i5,i6,i7,i8,za,zb,amp78xy,amp87xy)
      call qq4lggamp(i1,i2,i5,i6,i3,i4,i7,i8,za,zb,amp78yx,amp87yx)
      tmp78=amp78yx
      tmp87=amp87yx
      do h3=1,2
      do h5=1,2
      amp78yx(:,:,h3,h5,:,:)=tmp78(:,:,h5,h3,:,:)
      amp87yx(:,:,h3,h5,:,:)=tmp87(:,:,h5,h3,:,:)
      enddo
      enddo
      a78(:,:,:,:,:,:)=amp78xy(:,:,:,:,:,:)+amp78yx(:,:,:,:,:,:)
      a87(:,:,:,:,:,:)=amp87xy(:,:,:,:,:,:)+amp87yx(:,:,:,:,:,:)
      aq(:,:,:,:,:,:)=a78(:,:,:,:,:,:)+a87(:,:,:,:,:,:)

      if (interference) then
         call qq4lggamp(i1,i2,i3,i6,i5,i4,i7,i8,za,zb,
     &    amp78xy_swap,amp87xy_swap)
         call qq4lggamp(i1,i2,i5,i4,i3,i6,i7,i8,za,zb,
     &    amp78yx_swap,amp87yx_swap)
         tmp78_swap=amp78yx_swap
         tmp87_swap=amp87yx_swap
         do h3=1,2
         do h5=1,2
         amp78yx_swap(:,:,h3,h5,:,:)=tmp78_swap(:,:,h5,h3,:,:)
         amp87yx_swap(:,:,h3,h5,:,:)=tmp87_swap(:,:,h5,h3,:,:)
         enddo
         enddo
         a78_swap(:,:,:,:,:,:)=
     &    amp78xy_swap(:,:,:,:,:,:)+amp78yx_swap(:,:,:,:,:,:)
         a87_swap(:,:,:,:,:,:)=
     &    amp87xy_swap(:,:,:,:,:,:)+amp87yx_swap(:,:,:,:,:,:)
         aq_swap(:,:,:,:,:,:)=
     &    a78_swap(:,:,:,:,:,:)+a87_swap(:,:,:,:,:,:)
      endif

      msq(:)=zip
      do jdu=1,2
      do h1=1,2
      do h3=1,2
      do h5=1,2
      do h7=1,2
      do h8=1,2
C23456789012345678901234567890123456789012345678901234567890123456789012
      msq(jdu)=msq(jdu)+gsq**2*esq**4*V*xn*(
     & +dble(a78(jdu,h1,h3,h5,h7,h8)*Dconjg(a78(jdu,h1,h3,h5,h7,h8)))
     & +dble(a87(jdu,h1,h3,h5,h7,h8)*Dconjg(a87(jdu,h1,h3,h5,h7,h8)))
     & -dble(aq(jdu,h1,h3,h5,h7,h8)*Dconjg(aq(jdu,h1,h3,h5,h7,h8)))
     &  /xn**2)

         if(interference) then
      msq(jdu)=msq(jdu)+gsq**2*esq**4*V*xn*(
     & +dble(a78_swap(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(a78_swap(jdu,h1,h3,h5,h7,h8)))
     & +dble(a87_swap(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(a87_swap(jdu,h1,h3,h5,h7,h8)))
     & -dble(aq_swap(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(aq_swap(jdu,h1,h3,h5,h7,h8)))
     &  /xn**2)

      if (h3.eq.h5) then
      msq(jdu)=msq(jdu)-2d0*gsq**2*esq**4*V*xn/sqrt(colfac34_56)*(
     & +dble(a78(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(a78_swap(jdu,h1,h3,h5,h7,h8)))
     & +dble(a87(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(a87_swap(jdu,h1,h3,h5,h7,h8)))
     & -dble(aq(jdu,h1,h3,h5,h7,h8)
     &    *Dconjg(aq_swap(jdu,h1,h3,h5,h7,h8)))
     &  /xn**2)
      endif
         endif

      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

      msq(:) = msq(:)*colfac34_56*vsymfact

      return
      end
