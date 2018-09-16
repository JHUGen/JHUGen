      subroutine msq_qq(ja,ka,amp1_a,amp1_b,amp2_a,amp2_b,prop,
     .                  msq0,msq1,msq2)
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer ja,ka
      double complex a111,a112,a121,a211,a122,a212,a221,a222
      double complex b111,b112,b121,b211,b122,b212,b221,b222
      double complex prop
      double complex amp1_a(2,2,2),amp1_b(2,2,2)
      double complex amp2_a(2,2,2),amp2_b(2,2,2)
      double precision msq0,msq1,msq2

      if (ja .ne. ka) then
         a111=(Q(ja)*q1+L(ja)*l1*prop)*amp1_a(1,1,1)
     .       +(Q(ka)*q1+L(ka)*l1*prop)*amp1_b(1,1,1)
         a121=(Q(ja)*q1+L(ja)*l1*prop)*amp1_a(1,2,1)
     .       +(Q(ka)*q1+R(ka)*l1*prop)*amp1_b(1,2,1)
         a112=(Q(ja)*q1+L(ja)*r1*prop)*amp1_a(1,1,2)
     .       +(Q(ka)*q1+L(ka)*r1*prop)*amp1_b(1,1,2)
         a122=(Q(ja)*q1+L(ja)*r1*prop)*amp1_a(1,2,2)
     .       +(Q(ka)*q1+R(ka)*r1*prop)*amp1_b(1,2,2)
         a211=(Q(ja)*q1+R(ja)*l1*prop)*amp1_a(2,1,1)
     .       +(Q(ka)*q1+L(ka)*l1*prop)*amp1_b(2,1,1)
         a221=(Q(ja)*q1+R(ja)*l1*prop)*amp1_a(2,2,1)
     .       +(Q(ka)*q1+R(ka)*l1*prop)*amp1_b(2,2,1)
         a212=(Q(ja)*q1+R(ja)*r1*prop)*amp1_a(2,1,2)
     .       +(Q(ka)*q1+L(ka)*r1*prop)*amp1_b(2,1,2)
         a222=(Q(ja)*q1+R(ja)*r1*prop)*amp1_a(2,2,2)
     .       +(Q(ka)*q1+R(ka)*r1*prop)*amp1_b(2,2,2)

         msq0=zip
         msq1=(abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .        +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
         msq2=zip
      
      elseif (ja .eq. ka) then
         a111=(Q(ja)*q1+L(ja)*l1*prop)*(amp1_a(1,1,1)+amp1_b(1,1,1))
         b111=(Q(ja)*q1+L(ja)*l1*prop)*(amp2_a(1,1,1)+amp2_b(1,1,1))
         a112=(Q(ja)*q1+L(ja)*r1*prop)*(amp1_a(1,1,2)+amp1_b(1,1,2))
         b112=(Q(ja)*q1+L(ja)*r1*prop)*(amp2_a(1,1,2)+amp2_b(1,1,2))
         a221=(Q(ja)*q1+R(ja)*l1*prop)*(amp1_a(2,2,1)+amp1_b(2,2,1))
         b221=(Q(ja)*q1+R(ja)*l1*prop)*(amp2_a(2,2,1)+amp2_b(2,2,1))
         a222=(Q(ja)*q1+R(ja)*r1*prop)*(amp1_a(2,2,2)+amp1_b(2,2,2))
         b222=(Q(ja)*q1+R(ja)*r1*prop)*(amp2_a(2,2,2)+amp2_b(2,2,2))

         a121=(Q(ja)*q1+L(ja)*l1*prop)*amp1_a(1,2,1)
     .       +(Q(ka)*q1+R(ka)*l1*prop)*amp1_b(1,2,1)
         b121=(Q(ja)*q1+L(ja)*l1*prop)*amp2_a(1,2,1)
     .       +(Q(ka)*q1+R(ka)*l1*prop)*amp2_b(1,2,1)
         a122=(Q(ja)*q1+L(ja)*r1*prop)*amp1_a(1,2,2)
     .       +(Q(ka)*q1+R(ka)*r1*prop)*amp1_b(1,2,2)
         b122=(Q(ja)*q1+L(ja)*r1*prop)*amp2_a(1,2,2)
     .       +(Q(ka)*q1+R(ka)*r1*prop)*amp2_b(1,2,2)
         a211=(Q(ja)*q1+R(ja)*l1*prop)*amp1_a(2,1,1)
     .       +(Q(ka)*q1+L(ka)*l1*prop)*amp1_b(2,1,1)
         b211=(Q(ja)*q1+R(ja)*l1*prop)*amp2_a(2,1,1)
     .       +(Q(ka)*q1+L(ka)*l1*prop)*amp2_b(2,1,1)
         a212=(Q(ja)*q1+R(ja)*r1*prop)*amp1_a(2,1,2)
     .       +(Q(ka)*q1+L(ka)*r1*prop)*amp1_b(2,1,2)
         b212=(Q(ja)*q1+R(ja)*r1*prop)*amp2_a(2,1,2)
     .       +(Q(ka)*q1+L(ka)*r1*prop)*amp2_b(2,1,2)

         msq0=half*(
     .    +Dble(a111*dconjg(b111))+Dble(a112*dconjg(b112))
     .    +Dble(a221*dconjg(b221))+Dble(a222*dconjg(b222)))*two/xn
         msq1=half*
     .    (abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .    +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)
         msq2=half*(
     .    +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .    +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
      endif

      return
      end
      
      subroutine msq_qqb(ja,ka,amp1_a,amp1_b,amp2_a,amp2_b,prop,
     .                  msq0,msq1,msq2,msq_up,msq_down)
      implicit none
      include 'constants.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      integer ja,ka
      double complex a111,a112,a121,a211,a122,a212,a221,a222
      double complex b111,b112,b121,b211,b122,b212,b221,b222
      double complex prop
      double complex amp1_a(2,2,2),amp1_b(2,2,2)
      double complex amp2_a(2,2,2),amp2_b(2,2,2)
      double precision msq0,msq1,msq2,msq_up,msq_down

      if (ja .ne. -ka) then 
         a111=(Q(+ja)*q1+L(+ja)*l1*prop)*amp1_a(1,1,1)
     .       +(Q(-ka)*q1+L(-ka)*l1*prop)*amp1_b(1,1,1)
         a112=(Q(+ja)*q1+L(+ja)*r1*prop)*amp1_a(1,1,2)
     .       +(Q(-ka)*q1+L(-ka)*r1*prop)*amp1_b(1,1,2)
         a221=(Q(+ja)*q1+R(+ja)*l1*prop)*amp1_a(2,2,1)
     .       +(Q(-ka)*q1+R(-ka)*l1*prop)*amp1_b(2,2,1)
         a222=(Q(+ja)*q1+R(+ja)*r1*prop)*amp1_a(2,2,2)
     .       +(Q(-ka)*q1+R(-ka)*r1*prop)*amp1_b(2,2,2)

         a121=(Q(+ja)*q1+L(+ja)*l1*prop)*amp1_a(1,2,1)
     .       +(Q(-ka)*q1+R(-ka)*l1*prop)*amp1_b(1,2,1)
         a122=(Q(+ja)*q1+L(+ja)*r1*prop)*amp1_a(1,2,2)
     .       +(Q(-ka)*q1+R(-ka)*r1*prop)*amp1_b(1,2,2)
         a211=(Q(+ja)*q1+R(+ja)*l1*prop)*amp1_a(2,1,1)
     .       +(Q(-ka)*q1+L(-ka)*l1*prop)*amp1_b(2,1,1)
         a212=(Q(+ja)*q1+R(+ja)*r1*prop)*amp1_a(2,1,2)
     .       +(Q(-ka)*q1+L(-ka)*r1*prop)*amp1_b(2,1,2)

         msq0=zip
         msq1=zip
         msq2=(
     .   +abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .   +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

         msq_up=0d0
         msq_down=0d0

      elseif (ja .eq. -ka) then
c--case where final state from annihilation diagrams is the same quark
         a111=(Q(ja)*q1+L(ja)*l1*prop)*(amp1_a(1,1,1)+amp1_b(1,1,1))
         b111=(Q(ja)*q1+L(ja)*l1*prop)*(amp2_a(1,1,1)+amp2_b(1,1,1))
 
         a112=(Q(ja)*q1+L(ja)*r1*prop)*(amp1_a(1,1,2)+amp1_b(1,1,2))
         b112=(Q(ja)*q1+L(ja)*r1*prop)*(amp2_a(1,1,2)+amp2_b(1,1,2))
 
         a221=(Q(ja)*q1+R(ja)*l1*prop)*(amp1_a(2,2,1)+amp1_b(2,2,1))
         b221=(Q(ja)*q1+R(ja)*l1*prop)*(amp2_a(2,2,1)+amp2_b(2,2,1))

         a222=(Q(ja)*q1+R(ja)*r1*prop)*(amp1_a(2,2,2)+amp1_b(2,2,2))
         b222=(Q(ja)*q1+R(ja)*r1*prop)*(amp2_a(2,2,2)+amp2_b(2,2,2))
 
         a121=(Q(+ja)*q1+L(+ja)*l1*prop)*amp1_a(1,2,1)
     .       +(Q(-ka)*q1+R(-ka)*l1*prop)*amp1_b(1,2,1)
         a122=(Q(+ja)*q1+L(+ja)*r1*prop)*amp1_a(1,2,2)
     .       +(Q(-ka)*q1+R(-ka)*r1*prop)*amp1_b(1,2,2)
         a211=(Q(+ja)*q1+R(+ja)*l1*prop)*amp1_a(2,1,1)
     .       +(Q(-ka)*q1+L(-ka)*l1*prop)*amp1_b(2,1,1)
         a212=(Q(+ja)*q1+R(+ja)*r1*prop)*amp1_a(2,1,2)
     .       +(Q(-ka)*q1+L(-ka)*r1*prop)*amp1_b(2,1,2)
 
         b121=(Q(+ja)*q1+L(+ja)*l1*prop)*amp2_a(1,2,1)
     .       +(Q(-ka)*q1+R(-ka)*l1*prop)*amp2_b(1,2,1)
         b122=(Q(+ja)*q1+L(+ja)*r1*prop)*amp2_a(1,2,2)
     .       +(Q(-ka)*q1+R(-ka)*r1*prop)*amp2_b(1,2,2)
         b211=(Q(+ja)*q1+R(+ja)*l1*prop)*amp2_a(2,1,1)
     .       +(Q(-ka)*q1+L(-ka)*l1*prop)*amp2_b(2,1,1)
         b212=(Q(+ja)*q1+R(+ja)*r1*prop)*amp2_a(2,1,2)
     .       +(Q(-ka)*q1+L(-ka)*r1*prop)*amp2_b(2,1,2)

         msq0=(
     .   +Dble(a111*Dconjg(b111))+Dble(a112*Dconjg(b112))
     .   +Dble(a221*Dconjg(b221))+Dble(a222*Dconjg(b222)))*two/xn
         msq1=(
     .   +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .   +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
         msq2=(
     .   +abs(a111)**2+abs(a112)**2+abs(a221)**2+abs(a222)**2
     .   +abs(a122)**2+abs(a212)**2+abs(a121)**2+abs(a211)**2)

         b111=(Q(+ja)*q1+L(+ja)*l1*prop)*amp2_a(1,1,1)
     .       +(Q(+1)*q1+L(+1)*l1*prop)*amp2_b(1,1,1)
         b112=(Q(+ja)*q1+L(+ja)*r1*prop)*amp2_a(1,1,2)
     .       +(Q(+1)*q1+L(+1)*r1*prop)*amp2_b(1,1,2)
         b221=(Q(+ja)*q1+R(+ja)*l1*prop)*amp2_a(2,2,1)
     .       +(Q(+1)*q1+R(+1)*l1*prop)*amp2_b(2,2,1)
         b222=(Q(+ja)*q1+R(+ja)*r1*prop)*amp2_a(2,2,2)
     .       +(Q(+1)*q1+R(+1)*r1*prop)*amp2_b(2,2,2)
         b121=(Q(+ja)*q1+L(+ja)*l1*prop)*amp2_a(1,2,1)
     .       +(Q(+1)*q1+R(+1)*l1*prop)*amp2_b(1,2,1)
         b122=(Q(+ja)*q1+L(+ja)*r1*prop)*amp2_a(1,2,2)
     .       +(Q(+1)*q1+R(+1)*r1*prop)*amp2_b(1,2,2)
         b211=(Q(+ja)*q1+R(+ja)*l1*prop)*amp2_a(2,1,1)
     .       +(Q(+1)*q1+L(+1)*l1*prop)*amp2_b(2,1,1)
         b212=(Q(+ja)*q1+R(+ja)*r1*prop)*amp2_a(2,1,2)
     .       +(Q(+1)*q1+L(+1)*r1*prop)*amp2_b(2,1,2)
       
         msq_down=(
     .   +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .   +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)
 
         b111=(Q(+ja)*q1+L(+ja)*l1*prop)*amp2_a(1,1,1)
     .       +(Q(+2)*q1+L(+2)*l1*prop)*amp2_b(1,1,1)
         b112=(Q(+ja)*q1+L(+ja)*r1*prop)*amp2_a(1,1,2)
     .       +(Q(+2)*q1+L(+2)*r1*prop)*amp2_b(1,1,2)
         b221=(Q(+ja)*q1+R(+ja)*l1*prop)*amp2_a(2,2,1)
     .       +(Q(+2)*q1+R(+2)*l1*prop)*amp2_b(2,2,1)
         b222=(Q(+ja)*q1+R(+ja)*r1*prop)*amp2_a(2,2,2)
     .       +(Q(+2)*q1+R(+2)*r1*prop)*amp2_b(2,2,2)
         b121=(Q(+ja)*q1+L(+ja)*l1*prop)*amp2_a(1,2,1)
     .       +(Q(+2)*q1+R(+2)*l1*prop)*amp2_b(1,2,1)
         b122=(Q(+ja)*q1+L(+ja)*r1*prop)*amp2_a(1,2,2)
     .       +(Q(+2)*q1+R(+2)*r1*prop)*amp2_b(1,2,2)
         b211=(Q(+ja)*q1+R(+ja)*l1*prop)*amp2_a(2,1,1)
     .    +(Q(+2)*q1+L(+2)*l1*prop)*amp2_b(2,1,1)
         b212=(Q(+ja)*q1+R(+ja)*r1*prop)*amp2_a(2,1,2)
     .       +(Q(+2)*q1+L(+2)*r1*prop)*amp2_b(2,1,2)
 
         msq_up=(
     .   +abs(b111)**2+abs(b112)**2+abs(b221)**2+abs(b222)**2
     .   +abs(b122)**2+abs(b212)**2+abs(b121)**2+abs(b211)**2)

      endif
      
      return
      end
