!      - comments start with "!" rather than "C"
!      - continuation lines are avoided
   
      real(dp),parameter:: 
     & zip=0._dp,
     & zero=0._dp,
     & one=1._dp,
     & two=2._dp,
     & three=3._dp,
     & four=4._dp,
     & five=5._dp,
     & six=6._dp,
     & seven=7._dp,
     & eight=8._dp,
     & nine=9._dp,
     & ten=10._dp,
     & eleven=11._dp,
     & twelve=12._dp,
     & sixteen=16._dp,
     & half=0.5_dp,
     & quarter=0.25_dp


      real(dp),parameter::
     & pi=3.14159265358979311599796346854418516_dp,
     & pisq=pi*pi,
     & pisqo6=pisq/six,
     & twopi=two*pi,
     & fourpi=four*pi,
     & pion4=pi/four,  
     & pion10=pi/ten,
     & pisqm8=pisq-eight

      real(dp),parameter:: 
     & rt2=1.41421356237309504880168872420969798_dp,
     & twort2=two*rt2,
     & fourrt2=four*rt2,
     & rt2onpi=0.797884560802865371431347722487877591_dp 
! sqrt(two/pi)
!-----------------------------------------------------
!-----------------------------------------------------
      complex(dp),parameter:: 
     & im=(zip,one),
     & impi=(zip,pi),
     & czip=(zip,zip),
     & cone=(one,zip),
     & ctwo=(two,zip),
     & chalf=(half,zip)
!-----------------------------------------------------


      real(dp),parameter::
     & cf=four/three,
     & ca=three,
     & xn=three,
     & Nc=three,
     & Ncinv=one/three, 
     & xnsq=nine,
     & v=eight,
     & tr=half,
     & Von4=two,
     & ninth=one/nine,
     & xn4=xnsq-four,
     & qu=two/three,
     & qd=-one/three,
     & qe=-one,
     & spinave=one/four,
     & aveqq=spinave/xnsq,
     & aveqg=spinave/xn/v,
     & avegg=spinave/v**2,
     & aem=one/137.035989_dp

!-----------------------------------------------------

      real(dp),parameter:: 
     & nbGeV2=0.389379e6_dp,
     & pbGeV2=0.389379e9_dp,
     & fbGeV2=0.389379e12_dp,
!----decifemtobarns
     & dfbGeV2=0.389379e13_dp,
     & overa=pbGeV2/xn/256._dp/pi
      integer,parameter:: nloop=2,fn=-5

