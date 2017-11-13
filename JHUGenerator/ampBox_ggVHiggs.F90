function calc_MassiveTensorBox(MV2,MT2,MH2,shat,that,uhat,smT,y,LC,SP,DD0,DD1,DD2,DD3,kap,kapPr)
#if useCollier==1      
   use COLLIER
#endif
implicit none
real(8),parameter :: E=dexp(1d0)
real(8) :: shat,that,uhat,smT,y,MV2,MT2,MH2
complex(8) :: calc_MassiveTensorBox, LC(1:15),DD0(2:4),SP(1:7),kap,kapPr
integer,parameter :: rank=3
double complex :: DD1(0:rank/2,0:rank,0:rank,0:rank),DD1uv(0:rank/2,0:rank,0:rank,0:rank)
double complex :: DD2(0:rank/2,0:rank,0:rank,0:rank),DD2uv(0:rank/2,0:rank,0:rank,0:rank)
double complex :: DD3(0:rank/2,0:rank,0:rank,0:rank),DD3uv(0:rank/2,0:rank,0:rank,0:rank)



  calc_MassiveTensorBox = kap*(-48*DD1(1,0,0,1)*LC(1) - 48*DD2(1,0,0,1)*LC(1) + &
           48*DD1(1,1,0,0)*LC(11) + 48*DD2(1,1,0,0)*LC(11) + &
           16*DD2(1,0,0,0)*(LC(11) - LC(15)) + &
           48*DD2(1,0,1,0)*(LC(11) - LC(15)) + &
           16*DD1(1,0,0,0)*(-LC(1) + LC(15)) + &
           48*DD1(1,0,1,0)*(-LC(1) + LC(15)) + &
           16*DD1(0,1,0,2)*LC(3)*SP(3) + 16*DD2(0,1,0,2)*LC(3)*SP(3) + &
           16*DD2(0,0,1,2)*(LC(3) - LC(9))*SP(3) + &
           16*DD1(0,0,1,2)*LC(9)*SP(3) + &
           DD2(0,0,0,2)*(-4*shat*LC(1) + (4*smT*LC(1))/E**y + &
              8*(LC(3) - LC(9))*SP(3)) + &
           DD3(0,1,0,2)*((-8*smT*LC(1))/E**y + 16*LC(9)*SP(3)) + &
           8*DD3(0,0,1,2)*((shat - smT/E**y)*LC(1) + &
              2*(-LC(3) + LC(9))*SP(3)) + &
           DD1(0,0,0,2)*(4*(shat - smT/E**y)*LC(1) + &
              8*(-LC(3) + LC(9))*SP(3)) + &
           DD3(0,0,0,2)*(4*(shat - smT/E**y)*LC(1) + &
              8*(-LC(3) + LC(9))*SP(3)) + 16*DD1(0,2,0,1)*LC(3)*SP(4) + &
           16*DD2(0,2,0,1)*LC(3)*SP(4) - 16*DD2(0,2,1,0)*LC(13)*SP(4) + &
           16*DD1(0,2,1,0)*(LC(3) + LC(13))*SP(4) + &
           DD2(0,2,0,0)*(4*E**y*smT*LC(11) - 8*LC(13)*SP(4)) + &
           DD1(0,2,0,0)*(-4*E**y*smT*LC(11) + 8*LC(13)*SP(4)) + &
           DD2(0,1,0,0)*(4*(-2*MT2 + E**y*smT)*LC(11) - &
              8*(LC(3) + LC(13))*SP(4)) + &
           DD0(3)*((-2*smT*LC(11))/E**y + &
              4*MT2*(LC(1) - 2*LC(11) + LC(15)) + 4*LC(13)*SP(3) - &
              8*LC(3)*SP(4) - 4*LC(5)*SP(5)) + &
           DD0(2)*(2*E**y*smT*LC(1) + &
              MT2*(8*LC(1) - 4*(LC(11) + LC(15))) - 8*LC(3)*SP(3) - &
              4*LC(9)*SP(4) + 4*LC(5)*SP(5)) + &
           4*DD2(0,0,0,1)*((2*MT2 - shat + smT/E**y)*LC(1) + &
              LC(6)*SP(2) + (2*LC(3) - 2*LC(9) + LC(13))*SP(3) - &
              LC(9)*SP(4) - LC(5)*SP(5) + LC(4)*SP(6)) + &
           4*DD1(0,0,0,1)*((2*MT2 + shat - smT/E**y)*LC(1) - &
              LC(7)*SP(1) - 4*LC(3)*SP(3) + (2*LC(9) + LC(13))*SP(3) - &
              LC(9)*SP(4) + LC(5)*SP(5) - LC(2)*SP(7)) + &
           4*DD3(0,0,0,1)*((shat - smT/E**y)*LC(1) + LC(7)*SP(1) - &
              LC(6)*SP(2) + LC(4)*SP(6) - LC(2)*SP(7)) - &
           4*DD1(0,1,0,1)*(shat*LC(15) + LC(7)*SP(1) - &
              (4*LC(3) + LC(13))*SP(3) + (-2*LC(3) + LC(9))*SP(4) - &
              LC(5)*SP(5) + LC(2)*SP(7)) + &
           4*DD1(0,1,1,1)*((-3*smT*LC(11))/E**y - 3*shat*LC(15) + &
              8*LC(3)*SP(3) - 2*LC(9)*SP(4) - 2*LC(4)*SP(6) + &
              2*LC(2)*SP(7)) + &
           4*DD3(0,0,1,0)*(-2*MV2*LC(1) + (smT*LC(11))/E**y + &
              shat*LC(15) + LC(7)*SP(1) + LC(6)*SP(2) + 2*LC(3)*SP(4) + &
              3*LC(4)*SP(6) - 2*LC(14)*SP(6) + 3*(LC(2) - 2*LC(8))*SP(7) &
              ) + 4*DD2(0,0,1,0)* &
            (E**y*smT*LC(11) + 2*MT2*(-LC(11) + LC(15)) - LC(6)*SP(2) + &
              LC(13)*(SP(3) - 2*SP(4)) + (-2*LC(3) + LC(9))*SP(4) - &
              LC(5)*SP(5) + (LC(4) - 2*LC(14))*SP(6) + &
              2*(LC(2) - LC(8))*SP(7)) + &
           DD3(0,1,0,0)*(-8*MV2*LC(1) + 4*shat*LC(15) + 8*LC(7)*SP(1) + &
              8*LC(4)*SP(6) - 8*LC(14)*SP(6) - 24*LC(8)*SP(7)) + &
           DD0(4)*(4*(MT2 - MV2)*LC(1) - 4*MT2*LC(11) + 2*shat*LC(15) + &
              4*LC(7)*SP(1) + 4*LC(4)*SP(6) - 8*LC(8)*SP(7)) + &
           8*DD1(0,0,3,0)*(-3*MV2*LC(1) + MV2*LC(15) + 2*LC(9)*SP(3) + &
              2*LC(10)*SP(6) - 2*LC(8)*SP(7)) + &
           8*DD3(0,2,0,1)*(2*MV2*LC(1) - (smT*LC(15))/E**y - &
              2*LC(10)*SP(6) + 2*LC(8)*SP(7)) - &
           8*DD1(0,0,2,1)*(3*MV2*LC(1) - 4*LC(9)*SP(3) - &
              2*LC(10)*SP(6) + 2*LC(8)*SP(7)) + &
           4*DD3(0,0,1,1)*(2*MV2*LC(1) + shat*LC(1) - &
              (smT*(LC(1) - LC(11) + LC(15)))/E**y + 2*LC(3)*SP(4) + &
              4*LC(4)*SP(6) - 2*LC(10)*SP(6) + 2*LC(8)*SP(7)) - &
           4*DD3(0,2,0,0)*(MV2*LC(1) - MV2*LC(11) + 4*LC(14)*SP(6) + &
              4*LC(8)*SP(7)) + &
           2*DD3(0,1,0,1)*(4*MV2*LC(1) + &
              (smT*(-2*LC(1) + LC(11) - 2*LC(15)))/E**y + shat*LC(15) + &
              4*LC(7)*SP(1) - 4*LC(6)*SP(2) + 2*LC(9)*SP(4) + &
              10*LC(4)*SP(6) - 4*LC(10)*SP(6) - 2*LC(2)*SP(7) + &
              4*LC(8)*SP(7)) + &
           4*DD3(0,1,1,1)*(8*MV2*LC(1) + 2*E**y*smT*LC(1) - &
              shat*LC(15) - (smT*(LC(11) + 4*LC(15)))/E**y - &
              2*(LC(9)*SP(4) + (LC(4) + 4*LC(10))*SP(6)) + &
              2*(LC(2) + 4*LC(8))*SP(7)) + &
           8*DD3(0,2,1,0)*(-2*MV2*LC(11) + E**y*smT*LC(15) + &
              2*LC(14)*SP(6) - 2*LC(12)*SP(7)) + &
           8*DD2(0,1,2,0)*(3*MV2*LC(11) - 4*LC(13)*SP(4) - &
              2*LC(14)*SP(6) + 2*LC(12)*SP(7)) + &
           8*DD2(0,0,3,0)*(3*MV2*LC(11) - MV2*LC(15) - 2*LC(13)*SP(4) - &
              2*LC(14)*SP(6) + 2*LC(12)*SP(7)) - &
           8*DD3(0,1,2,0)*(4*MV2*LC(11) + &
              E**y*smT*(LC(11) - 2*LC(15)) - 2*LC(13)*SP(4) - &
              4*LC(14)*SP(6) + 4*LC(12)*SP(7)) + &
           4*DD1(0,1,2,0)*(6*MV2*LC(11) - (3*smT*LC(11))/E**y - &
              3*shat*LC(15) + 4*LC(3)*SP(3) - &
              2*(LC(9)*SP(4) + (LC(4) + 2*LC(14))*SP(6)) + &
              2*(LC(2) + 2*LC(12))*SP(7)) + &
           DD2(0,1,0,1)*(-2*E**y*smT*LC(1) + (2*smT*LC(11))/E**y + &
              4*shat*(-LC(1) + LC(11) + LC(15)) + &
              4*(-(LC(7)*SP(1)) + 4*LC(3)*SP(3) + 6*LC(3)*SP(4) + &
                 LC(5)*SP(5) - LC(2)*SP(7))) + &
           4*DD2(0,1,1,1)*((3*smT*LC(11))/E**y + 3*shat*LC(15) + &
              2*((4*LC(3) + LC(9))*SP(4) + LC(4)*SP(6) - LC(2)*SP(7))) &
            + 4*DD2(0,0,2,1)* &
            (-6*MV2*LC(1) + (3*smT*LC(11))/E**y + 3*shat*LC(15) + &
              2*((2*LC(3) + LC(9))*SP(4) + (LC(4) + 2*LC(10))*SP(6) - &
                 (LC(2) + 2*LC(8))*SP(7))) + &
           2*DD2(0,0,1,1)*(-4*MV2*LC(1) - 2*shat*LC(1) + &
              2*shat*LC(11) + (5*smT*LC(11))/E**y + 3*shat*LC(15) - &
              4*smT*LC(1)*Sinh(y) + &
              2*(4*(LC(3) - LC(9))*SP(3) + (6*LC(3) + LC(9))*SP(4) + &
                 (5*LC(4) + 2*LC(10))*SP(6) - (LC(2) + 2*LC(8))*SP(7))) &
            + (DD1(0,0,1,1)*(-4*smT*(2*LC(1) + LC(11)) - &
                4*E**y*(2*(MV2 - shat)*LC(1) + LC(7)*SP(1) + &
                   (4*LC(3) - 8*LC(9) - LC(13))*SP(3) + LC(9)*SP(4) - &
                   LC(5)*SP(5) + 4*LC(4)*SP(6) - 2*LC(10)*SP(6) + &
                   (LC(2) + 2*LC(8))*SP(7))))/E**y + &
           DD1(0,0,1,0)*(4*(2*MT2 + shat)*LC(1) + 2*E**y*smT*LC(1) - &
              (2*smT*(2*LC(1) + LC(11)))/E**y - 8*MT2*LC(15) - &
              4*(LC(7)*SP(1) - LC(6)*SP(2) + &
                 (4*LC(3) - 2*LC(9) - LC(13))*SP(3) + 3*LC(9)*SP(4) - &
                 2*LC(5)*SP(5) + 3*LC(4)*SP(6) + 2*LC(14)*SP(6) + &
                 (LC(2) + 2*LC(8))*SP(7))) + &
           (DD3(0,0,2,1)*(8*E**(2*y)*smT*LC(1) - &
                4*smT*(LC(11) + 2*LC(15)) + &
                4*E**y*(4*MV2*LC(1) - shat*(2*LC(11) + LC(15)) - &
                   2*(2*LC(3) + LC(9))*SP(4) - &
                   2*(LC(4) + 2*LC(10))*SP(6) + &
                   2*(LC(2) + 2*LC(8))*SP(7))))/E**y + &
           (4*DD1(0,0,2,0)*(-(smT*(LC(1) + LC(11))) + &
                E**y*(shat*LC(1) + MV2*(-4*LC(1) + LC(11) + LC(15)) - &
                   LC(7)*SP(1) - 2*LC(3)*SP(3) + &
                   (6*LC(9) + LC(13))*SP(3) - LC(9)*SP(4) + &
                   LC(5)*SP(5) + &
                   2*(-2*LC(4) + LC(10) - 2*LC(14))*SP(6) - &
                   (LC(2) + 6*LC(8))*SP(7))))/E**y + &
           2*DD3(0,1,1,0)*(-4*MV2*LC(1) + (smT*LC(11))/E**y - &
              2*E**y*smT*(LC(11) - LC(15)) + shat*LC(15) + &
              2*(LC(9)*SP(4) + (LC(4) - 6*LC(14))*SP(6) + &
                 (3*LC(2) - 8*LC(8) - 2*LC(12))*SP(7))) - &
           8*DD3(0,0,3,0)*(2*MV2*LC(11) + E**y*smT*(LC(11) - LC(15)) - &
              2*(LC(13)*SP(4) + LC(14)*SP(6) - LC(12)*SP(7))) + &
           2*DD1(0,1,1,0)*(2*E**y*smT*(LC(1) - LC(11)) + 4*MV2*LC(11) - &
              (3*smT*LC(11))/E**y - 5*shat*LC(15) + &
              2*(-(LC(7)*SP(1)) + (4*LC(3) + LC(13))*SP(3) + &
                 2*(LC(3) - 2*LC(9) + 2*LC(13))*SP(4) + LC(5)*SP(5) - &
                 (3*LC(4) + 2*LC(14))*SP(6) + 2*(-LC(2) + LC(12))*SP(7)) &
              ) + 2*DD2(0,1,1,0)* &
            (4*MV2*LC(11) + shat*LC(15) + &
              smT*LC(11)*(5*Cosh(y) + 3*Sinh(y)) + &
              2*((LC(9) - 8*LC(13))*SP(4) + (LC(4) - 2*LC(14))*SP(6) + &
                 (3*LC(2) + 2*LC(12))*SP(7))) + &
           2*DD2(0,0,2,0)*(shat*LC(15) - &
              2*MV2*(LC(1) - 4*LC(11) + LC(15)) + &
              smT*LC(11)*(3*Cosh(y) + Sinh(y)) + &
              2*((LC(9) - 6*LC(13))*SP(4) + (LC(4) - 6*LC(14))*SP(6) + &
                 (3*LC(2) - 4*LC(8) + 2*LC(12))*SP(7))) + &
           DD1(0,1,0,0)*(2*E**y*smT*(LC(1) - 2*LC(11)) - &
              2*(4*MT2*LC(11) + shat*LC(15) + &
                 2*(LC(6)*SP(2) + (LC(9) - 2*LC(13))*SP(4) - &
                    LC(5)*SP(5) + LC(2)*SP(7)))) + &
           (2*DD3(0,0,2,0)*(smT*LC(11) + &
                2*E**(2*y)*smT*(-2*LC(11) + LC(15)) + &
                E**y*(-2*MV2*(LC(1) + LC(11)) + shat*LC(15) + &
                   2*((LC(9) + 2*LC(13))*SP(4) + &
                      (LC(4) - 2*LC(14))*SP(6) + &
                      (3*LC(2) - 4*LC(8) - 2*LC(12))*SP(7)))))/E**y) + &
        kapPr*(-4*shat*DD1(0,0,0,1)*SP(3)*SP(5) + &
           4*(shat - smT/E**y)*DD3(0,0,0,1)*SP(3)*SP(5) - &
           4*shat*DD2(0,1,0,0)*SP(4)*SP(5) + &
           (8*smT*DD1(0,0,2,1)*SP(2)*SP(6))/E**y - &
           8*(MV2 - smT/E**y)*DD1(0,0,3,0)*SP(2)*SP(6) + &
           8*shat*DD1(0,1,1,1)*SP(2)*SP(6) + &
           8*(shat - E**y*smT)*DD1(0,1,2,0)*SP(2)*SP(6) - &
           48*DD1(1,0,1,0)*SP(2)*SP(6) + &
           8*(shat - smT/E**y)*DD2(0,0,2,1)*SP(1)*SP(7) - &
           8*(MV2 - E**y*smT)*DD2(0,0,3,0)*SP(1)*SP(7) + &
           8*shat*DD2(0,1,1,1)*SP(1)*SP(7) + &
           8*E**y*smT*DD2(0,1,2,0)*SP(1)*SP(7) - &
           48*DD2(1,0,1,0)*SP(1)*SP(7) + &
           16*DD3(1,0,0,0)*(SP(2)*SP(6) - 2*SP(1)*SP(7)) - &
           (8*(E**y*shat - smT)*DD3(0,0,2,1)* &
              (SP(2)*SP(6) - SP(1)*SP(7)))/E**y + &
           8*(MV2 - E**y*smT)*DD3(0,0,3,0)* &
            (SP(2)*SP(6) - SP(1)*SP(7)) - &
           (8*(E**y*shat - 2*smT)*DD3(0,1,1,1)* &
              (SP(2)*SP(6) - SP(1)*SP(7)))/E**y + &
           8*(3*MV2 - 2*E**y*smT)*DD3(0,1,2,0)* &
            (SP(2)*SP(6) - SP(1)*SP(7)) + &
           (8*smT*DD3(0,2,0,1)*(SP(2)*SP(6) - SP(1)*SP(7)))/E**y + &
           8*(3*MV2 - E**y*smT)*DD3(0,2,1,0)* &
            (SP(2)*SP(6) - SP(1)*SP(7)) + &
           8*MV2*DD3(0,3,0,0)*(SP(2)*SP(6) - SP(1)*SP(7)) + &
           48*DD3(1,0,1,0)*(SP(2)*SP(6) - SP(1)*SP(7)) + &
           48*DD3(1,1,0,0)*(SP(2)*SP(6) - SP(1)*SP(7)) + &
           4*shat*DD1(0,1,0,1)* &
            (SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7)) + &
           4*(shat - E**y*smT)*DD1(0,1,1,0)* &
            (SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7)) + &
           16*DD2(1,0,0,0)*(SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7)) + &
           16*DD1(1,0,0,0)*(SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7)) + &
           (4*(E**y*shat - smT)*DD2(0,0,1,1)* &
              (SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7)))/E**y + &
           4*shat*DD2(0,1,0,1)* &
            (SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7)) + &
           (4*(E**y*shat - smT)*DD3(0,0,1,1)* &
              ((SP(3) + SP(4))*SP(5) - 2*SP(2)*SP(6) + 2*SP(1)*SP(7)))/ &
            E**y + DD0(2)*(2*E**y*smT*SP(3)*SP(5) - &
              (2*smT*SP(4)*SP(5))/E**y - &
              2*shat*(2*SP(3)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7)) + &
              4*MT2*(SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7))) + &
           DD1(0,0,1,1)*(8*shat*SP(2)*SP(6) + &
              (4*smT*(SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7)))/E**y) + &
           DD2(0,1,1,0)*(8*shat*SP(1)*SP(7) + &
              4*E**y*smT*(SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7))) + &
           DD0(3)*(-2*E**y*smT*SP(3)*SP(5) + (2*smT*SP(4)*SP(5))/E**y + &
              4*MT2*(SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7)) - &
              2*shat*(2*SP(4)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7))) + &
           (DD3(0,1,0,1)*(-8*E**y*shat*SP(2)*SP(6) - &
                4*smT*((SP(3) + SP(4))*SP(5) - 2*SP(2)*SP(6) + &
                   2*SP(1)*SP(7))))/E**y + &
           4*DD3(0,0,1,0)*(shat*SP(4)*SP(5) - 2*MT2*SP(2)*SP(6) - &
              shat*SP(2)*SP(6) + &
              E**y*smT*(2*SP(3)*SP(5) - SP(2)*SP(6)) + &
              MV2*(-2*SP(3)*SP(5) + 2*SP(2)*SP(6)) + &
              ((2*MT2 + 3*shat)*SP(1) + 2*(3*SP(3) + SP(4))*SP(6))* &
               SP(7) - (smT*(SP(4)*SP(5) + 3*SP(1)*SP(7)))/E**y) + &
           (DD0(4)*(2*E**(2*y)*smT*SP(3)*SP(5) - &
                2*smT*(SP(4)*SP(5) + 2*SP(1)*SP(7)) + &
                2*E**y*(2*(-(MV2*SP(3)) + MT2*(SP(3) + SP(4)))*SP(5) - &
                   shat*SP(2)*SP(6) + &
                   ((4*MT2 + shat)*SP(1) + 4*SP(3)*SP(6))*SP(7))))/E**y &
            + (4*DD2(0,0,1,0)* &
              (-(E**(2*y)*smT*SP(2)*SP(6)) + &
                smT*(SP(4)*SP(5) - SP(1)*SP(7)) + &
                E**y*(shat*(-(SP(4)*SP(5)) + SP(2)*SP(6)) + &
                   ((2*MT2 + shat)*SP(1) + 2*(SP(3) + SP(4))*SP(6))* &
                    SP(7))))/E**y + &
           (DD3(0,0,2,0)*(-8*smT*SP(1)*SP(7) + &
                4*E**(2*y)*smT* &
                 ((SP(3) + SP(4))*SP(5) - 4*SP(2)*SP(6) + 2*SP(1)*SP(7)) &
                  - 4*E**y*(MV2*(SP(3) + SP(4))*SP(5) - &
                   4*MV2*SP(2)*SP(6) - &
                   2*((-MV2 + shat)*SP(1) + 2*(SP(3) + SP(4))*SP(6))* &
                    SP(7))))/E**y + &
           (4*DD2(0,0,2,0)*(-2*smT*SP(1)*SP(7) + &
                E**(2*y)*smT* &
                 (SP(3)*SP(5) - SP(2)*SP(6) + SP(1)*SP(7)) + &
                E**y*(MV2*(-(SP(3)*SP(5)) + SP(2)*SP(6)) + &
                   (-(MV2*SP(1)) + 2*shat*SP(1) + &
                      4*(SP(3) + SP(4))*SP(6))*SP(7))))/E**y + &
           (DD3(0,1,1,0)*(-16*smT*SP(1)*SP(7) + &
                4*E**(2*y)*smT* &
                 ((SP(3) + SP(4))*SP(5) - 6*SP(2)*SP(6) + 2*SP(1)*SP(7)) &
                  + 8*E**y*(-(MV2*(SP(3) + SP(4))*SP(5)) + &
                   4*MV2*SP(2)*SP(6) + &
                   ((-2*MV2 + shat)*SP(1) + 4*(SP(3) + SP(4))*SP(6))* &
                    SP(7))))/E**y + &
           (DD3(0,1,0,0)*(4*E**(2*y)*smT*(SP(3)*SP(5) - SP(2)*SP(6)) - &
                4*smT*(SP(4)*SP(5) + 3*SP(1)*SP(7)) + &
                E**y*(-8*MV2*SP(3)*SP(5) - &
                   4*(2*MT2 - 2*MV2 + shat)*SP(2)*SP(6) + &
                   4*((2*MT2 + shat)*SP(1) + 2*(3*SP(3) + SP(4))*SP(6))* &
                    SP(7))))/E**y - &
           (4*DD3(0,2,0,0)*(2*E**(2*y)*smT*SP(2)*SP(6) + &
                2*smT*SP(1)*SP(7) + &
                E**y*(-4*(SP(3) + SP(4))*SP(6)*SP(7) + &
                   MV2*((SP(3) + SP(4))*SP(5) - 4*SP(2)*SP(6) + &
                      2*SP(1)*SP(7)))))/E**y + &
           (4*DD1(0,0,1,0)*(E**(2*y)*smT*(SP(3)*SP(5) - SP(2)*SP(6)) - &
                smT*SP(1)*SP(7) + &
                E**y*(shat*(-(SP(3)*SP(5)) + SP(2)*SP(6) + &
                      SP(1)*SP(7)) + &
                   2*SP(6)*(MT2*SP(2) + (SP(3) + SP(4))*SP(7)))))/E**y &
            + 4*DD1(0,0,2,0)* &
            (-(MV2*(SP(4)*SP(5) + SP(2)*SP(6) - SP(1)*SP(7))) + &
              (smT*(SP(4)*SP(5) + (1 - 2*E**(2*y))*SP(2)*SP(6) - &
                    SP(1)*SP(7)) + &
                 2*E**y*SP(6)*(shat*SP(2) + 2*(SP(3) + SP(4))*SP(7)))/ &
               E**y))

  
  
  
  
return
end function




#if 0==1 



function calc_MassiveBox(MV2,MH2,MT2,shat,smT,y,LC,SP,SI,kap,kapPr)
implicit none
real(8),parameter :: E=dexp(1d0)
real(8) :: shat,smT,y,MV2,MH2,MT2,smT2
complex(8) :: calc_MassiveBox, LC(1:15),SI(0:16),SP(1:9),kap,kapPr

     smT2 = smT**2

! print *,"DB",shat,smT,y,MV2,MH2,MT2,smT2,SI(0:16),LC(1:15)

     calc_MassiveBox =  kap*SI(2)*((-8*LC(1))/shat - &
           (8*E**y*smT*LC(11))/(shat**2 - E**y*shat*smT) + &
           (16*E**y*LC(3)*SP(3))/(E**y*shat**2 - shat*smT) - &
           (16*E**y*LC(9)*SP(3))/(E**y*shat**2 - shat*smT) + &
           (16*LC(13)*SP(4))/(shat**2 - E**y*shat*smT)) + &
        kap*SI(7)*((2*E**y*(-MV2 + E**y*smT)* &
              (16*MV2**2*shat**2 + &
                MV2*shat*(3*E**y*shat - 37*smT)*smT + &
                9*smT**3*(E**y*shat + smT))*LC(1))/ &
            (shat*(E**y*shat - smT)*(-(MV2*shat) + smT2)**2) + &
           (2*(-(smT2/(E**y*shat - smT)) + &
                (2*smT**3*(MV2*shat + 2*smT2))/ &
                 (-(MV2*shat) + smT2)**2 + &
                (MV2*(-8*MV2**2*shat**2 + 7*MV2*shat*smT2 - 5*smT2**2))/ &
                 (E**y*(-(MV2*shat) + smT2)**2))*LC(11))/(shat*smT) - &
           (2*(2*MV2*smT*(3*MV2*shat + 2*smT2) + &
                E**(2*y)*shat*smT*(MV2*shat + 9*smT2) - &
                E**y*(5*MV2**2*shat**2 + 7*MV2*shat*smT2 + 8*smT2**2))* &
              LC(15))/((E**y*shat - smT)*(-(MV2*shat) + smT2)**2) - &
           (4*(MV2 - E**y*smT)*LC(7)*SP(1))/(shat*(MV2*shat - smT2)) + &
           (8*(MV2 - E**y*smT)*LC(6)*SP(2))/(shat*(MV2*shat - smT2)) + &
           (16*MV2*(MV2 - E**y*smT)*LC(13)*SP(4))/ &
            (E**(2*y)*(-(MV2*shat) + smT2)**2) - &
           (4*(-MV2 + E**y*smT)*smT2*LC(9)* &
              (4*E**(3*y)*shat*smT*(MV2*shat - 3*smT2)*SP(3) + &
                4*E**(4*y)*shat**2*smT2*SP(3) + 2*MV2*shat*smT2*SP(4) + &
                E**(2*y)*(-4* &
                    (MV2*shat**2*(3*MV2 + shat) - &
                      shat*(5*MV2 + shat)*smT2 + smT2**2)*SP(3) + &
                   shat**2*(MV2*shat + smT2)*SP(4)) + &
                E**y*shat*smT* &
                 (4*(MV2*shat - smT2)*SP(3) - (3*MV2*shat + smT2)*SP(4)) &
                ))/ &
            (shat*smT**2*(-(E**y*shat) + smT)**2* &
              (-(MV2*shat) + smT2)**2) + &
           (16*(-MV2 + E**y*smT)*LC(3)* &
              (E**(4*y)*shat**2*smT*(2*MV2*shat - smT2)*smT2*SP(3) - &
                E**(3*y)*shat*smT2* &
                 (4*MV2**2*shat**2 - 5*MV2*shat*smT2 + 3*smT2**2)*SP(3) &
                 + smT**5*(-2*MV2*shat + smT2)*SP(4) + &
                E**(2*y)*smT*smT2* &
                 ((MV2*(MV2 - shat)*shat**2 + shat*(-MV2 + shat)*smT2 + &
                      smT2**2)*SP(3) + &
                   shat**2*(-2*MV2*shat + smT2)*SP(4)) + &
                E**y*shat*smT2**2* &
                 (-(smT2*(SP(3) + 2*SP(4))) + &
                   MV2*shat*(SP(3) + 4*SP(4)))))/ &
            (E**y*shat**2*smT**2*(-(E**y*shat) + smT)**2* &
              (-(MV2*shat) + smT2)**2) - &
           (4*(MV2 - E**y*smT)*LC(5)*SP(5))/(shat*(MV2*shat - smT2)) + &
           (4*(MV2 - E**y*smT)*(-2*MV2*smT + E**y*(MV2*shat + smT2))* &
              LC(4)*SP(6))/((E**y*shat - smT)*(-(MV2*shat) + smT2)**2) &
            + (16*E**y*smT*(MV2 - E**y*smT)*LC(10)*SP(6))/ &
            (-(MV2*shat) + smT2)**2 + &
           (16*(MV2**2*shat - E**y*smT**3)*LC(14)*SP(6))/ &
            (E**y*smT*(-(MV2*shat) + smT2)**2) + &
           (4*smT*(-MV2 + E**y*smT)* &
              (MV2*shat + smT*(-2*E**y*shat + smT))*LC(2)*SP(7))/ &
            (shat*(-(E**y*shat) + smT)*(-(MV2*shat) + smT2)**2) + &
           (16*E**y*smT*(-MV2 + E**y*smT)*LC(8)*SP(7))/ &
            (-(MV2*shat) + smT2)**2 + &
           (16*(-(MV2**2*shat) + E**y*smT**3)*LC(12)*SP(7))/ &
            (E**y*smT*(-(MV2*shat) + smT2)**2)) + &
        kap*SI(6)*((4*smT*(-(E**y*MV2) + smT)* &
              (3*E**y*smT*(MV2*shat + smT2) - shat*(5*MV2*shat + smT2))* &
              LC(1))/(shat*(shat - E**y*smT)*(-(MV2*shat) + smT2)**2) + &
           ((6*MV2)/(MV2*shat - smT2) + &
              (-30*MV2**2*shat**2*smT + 56*MV2*shat*smT**3 - 14*smT**5)/ &
               (E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
              (16*(-(MV2*shat) + smT2))/(shat*(shat - E**y*smT)**2) - &
              (2*smT2*(MV2*shat + 5*smT2))/ &
               (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
              (2*(4*shat*(-6*MV2 + shat) + 7*smT2))/ &
               (shat**2*(shat - E**y*smT)))*LC(11) + &
           2*(-5/(shat - E**y*smT) + &
              (2*E**y*MV2*smT*(MV2*shat - 6*smT2) + &
                 E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
                 smT2*(3*MV2*shat + 7*smT2))/ &
               (E**y*smT*(-(MV2*shat) + smT2)**2))*LC(15) - &
           (4*(E**(2*y)*MV2*shat*smT + smT**3 - &
                E**y*MV2*(shat**2 + smT2))*LC(7)*SP(1))/ &
            (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) + &
           (8*(E**(2*y)*MV2*shat*smT + smT**3 - &
                E**y*MV2*(shat**2 + smT2))*LC(6)*SP(2))/ &
            (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) - &
           (4*(E**y*MV2 - smT)*LC(9)* &
              (-4*E**(2*y)*MV2*shat**2*SP(3) + &
                4*E**(3*y)*MV2*shat*smT*SP(3) + &
                shat*(-3*MV2*shat + smT2)*SP(4) + &
                E**y*smT*(MV2*shat + smT2)*SP(4)))/ &
            (E**y*shat*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2) + &
           (16*(E**y*MV2 - smT)*LC(3)* &
              (E**(4*y)*smT**3*(2*MV2*shat - smT2)*SP(3) - &
                2*E**(3*y)*shat*(2*MV2*shat - smT2)*smT2*SP(3) - &
                MV2*shat**3*smT*SP(4) + &
                2*E**y*shat*(MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)* &
                 SP(4) - E**(2*y)*smT* &
                 (shat**2*(-2*MV2*shat + smT2)*SP(3) + &
                   (MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4))))/ &
            (E**(2*y)*shat**2*(shat - E**y*smT)**2* &
              (-(MV2*shat) + smT2)**2) + &
           4*LC(13)*((E**y*MV2*SP(3))/(MV2*shat*smT - smT**3) - &
              (MV2*SP(3))/(MV2*shat**2 - shat*smT2) - &
              (4*smT*(4*MV2*shat + shat**2 - 4*smT2)*SP(4))/ &
               (E**y*shat**3*(MV2*shat - smT2)) + &
              (12*(MV2*shat - smT2)*SP(4))/ &
               (shat**2*(shat - E**y*smT)**2) + &
              (4*smT**3*SP(4))/ &
               (E**(3*y)*shat*(-(MV2*shat) + smT2)**2) - &
              (4*smT2**2*SP(4))/ &
               (E**(2*y)*shat**2*(-(MV2*shat) + smT2)**2) + &
              (-16*smT2*SP(4) + shat**2*(SP(3) + 4*SP(4)))/ &
               (shat**3*(shat - E**y*smT))) - &
           (4*(E**(2*y)*MV2*shat*smT + smT**3 - &
                E**y*MV2*(shat**2 + smT2))*LC(5)*SP(5))/ &
            (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) + &
           4*(1/(-shat**2 + E**y*shat*smT) + &
              (-2*E**y*MV2**2*shat*smT + &
                 E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
                 (3*MV2*shat - smT2)*smT2)/ &
               (E**y*shat*smT*(-(MV2*shat) + smT2)**2))*LC(4)*SP(6) + &
           (16*smT*(-(E**y*MV2) + smT)*LC(10)*SP(6))/ &
            (-(MV2*shat) + smT2)**2 + &
           (16*(E**y*MV2 - smT)*(shat*smT + E**y*(MV2*shat - 2*smT2))* &
              LC(14)*SP(6))/ &
            (E**(2*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2) - &
           (4*(2*E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
                3*MV2*shat*smT2 - smT2**2 + &
                E**y*MV2*smT*(-3*MV2*shat + smT2))*LC(2)*SP(7))/ &
            (E**y*shat*smT*(-(MV2*shat) + smT2)**2) - &
           (16*smT*(-(E**y*MV2) + smT)*LC(8)*SP(7))/ &
            (-(MV2*shat) + smT2)**2 + &
           (16*(-(E**y*MV2) + smT)* &
              (shat*smT + E**y*(MV2*shat - 2*smT2))*LC(12)*SP(7))/ &
            (E**(2*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2)) + &
        kap*SI(4)*((MV2*(-16*MV2**2*shat**2*smT**3 + &
                37*MV2*shat*smT**5 - 9*smT**7 - &
                6*E**(8*y)*smT**5*(MV2*shat + smT2) + &
                E**(7*y)*shat*smT2**2*(MV2*shat + 23*smT2) - &
                E**(5*y)*shat*smT2* &
                 (20*MV2**2*shat**2 + 131*MV2*shat*smT2 - 79*smT2**2) - &
                E**(4*y)*smT*(MV2*shat - smT2)* &
                 (16*MV2**2*shat**2 + 32*MV2*shat*smT2 - 45*smT2**2) + &
                E**(6*y)*smT**3* &
                 (46*MV2**2*shat**2 + 5*MV2*shat*smT2 - 27*smT2**2) + &
                E**y*shat*smT2* &
                 (56*MV2**2*shat**2 - 113*MV2*shat*smT2 + 33*smT2**2) - &
                E**(3*y)*shat* &
                 (80*MV2**3*shat**3 - 348*MV2**2*shat**2*smT2 + &
                   285*MV2*shat*smT2**2 - 89*smT2**3) - &
                E**(2*y)*smT* &
                 (8*MV2**3*shat**3 + 86*MV2**2*shat**2*smT2 - &
                   103*MV2*shat*smT2**2 + 33*smT2**3))*LC(1))/ &
            (2d0*E**(3*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (MV2*(3*E**y*shat*smT**5*(MV2*shat - 5*smT2) - &
                3*E**(9*y)*shat*smT**5*(MV2*shat - smT2) + &
                8*MV2**2*shat**2*smT2**2 - 7*MV2*shat*smT2**3 + &
                5*smT2**4 + 2*E**(8*y)*smT2**2* &
                 (12*MV2**2*shat**2 - 20*MV2*shat*smT2 + 5*smT2**2) - &
                2*E**(7*y)*shat*smT**3* &
                 (16*MV2**2*shat**2 - 35*MV2*shat*smT2 + 13*smT2**2) - &
                2*E**(3*y)*shat*smT**3* &
                 (10*MV2**2*shat**2 - 59*MV2*shat*smT2 + 31*smT2**2) + &
                4*E**(5*y)*shat*smT* &
                 (8*MV2**3*shat**3 - 47*MV2**2*shat**2*smT2 + &
                   49*MV2*shat*smT2**2 - 19*smT2**3) + &
                E**(2*y)*smT2* &
                 (-64*MV2**3*shat**3 + 88*MV2**2*shat**2*smT2 - &
                   61*MV2*shat*smT2**2 + 25*smT2**3) + &
                E**(6*y)*smT2* &
                 (-56*MV2**3*shat**3 + 144*MV2**2*shat**2*smT2 - &
                   111*MV2*shat*smT2**2 + 35*smT2**3) + &
                E**(4*y)*(-(MV2*shat) + smT2)* &
                 (-128*MV2**3*shat**3 + 112*MV2**2*shat**2*smT2 - &
                   80*MV2*shat*smT2**2 + 45*smT2**3))*LC(11))/ &
            (2d0*E**(5*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (MV2*(-6*MV2*shat*smT**5 - 4*smT**7 + &
                E**(6*y)*(-6*MV2**2*shat**2*smT**3 - &
                   46*MV2*shat*smT**5 + 32*smT**7) - &
                2*E**(8*y)*smT**5*(MV2*shat - 6*smT2) + &
                10*E**(2*y)*MV2*shat*smT**3*(7*MV2*shat - 5*smT2) - &
                2*E**y*shat*(3*MV2*shat - 13*smT2)*smT2**2 - &
                E**(9*y)*shat*(MV2*shat - smT2)*smT2**2 - &
                E**(3*y)*shat*smT2* &
                 (32*MV2**2*shat**2 + 57*MV2*shat*smT2 - 29*smT2**2) + &
                E**(7*y)*shat*smT2* &
                 (8*MV2**2*shat**2 - 7*MV2*shat*smT2 - 21*smT2**2) - &
                E**(5*y)*shat*(4*MV2*shat - 19*smT2)* &
                 (4*MV2**2*shat**2 + MV2*shat*smT2 - smT2**2) + &
                8*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
                 (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2))* &
              LC(15))/ &
            (2d0*E**(4*y)*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
           (4*MV2*(-(E**y*shat) + smT)*LC(7)*SP(1))/ &
            (shat*smT*(-(MV2*shat) + smT2)) + &
           (LC(6)*(8*E**y*MV2*shat*SP(2) - 8*MV2*smT*SP(2)))/ &
            (MV2*shat**2*smT - shat*smT**3) + &
           (MV2*LC(9)*(-8*E**(9*y)*shat*smT**7*SP(3) + &
                4*E**(10*y)*MV2*shat*smT2**3*SP(3) - &
                2*MV2*shat*smT2**3*SP(4) + &
                2*E**(7*y)*shat*smT**5* &
                 (6*(5*MV2*shat - 3*smT2)*SP(3) - &
                   (MV2*shat + smT2)*SP(4)) - &
                2*E**y*shat*smT**5* &
                 (2*(MV2*shat - smT2)*SP(3) - (MV2*shat + smT2)*SP(4)) &
                 + E**(8*y)*smT2**2* &
                 (4*(-8*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2)* &
                    SP(3) + smT2*(MV2*shat + smT2)*SP(4)) - &
                2*E**(3*y)*shat*smT**3*(2*MV2*shat - smT2)* &
                 (2*(MV2*shat - 3*smT2)*SP(3) + &
                   (5*MV2*shat + smT2)*SP(4)) + &
                E**(4*y)*smT2* &
                 (4*smT2*(4*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
                      3*smT2**2)*SP(3) + &
                   (-(MV2*shat) + smT2)* &
                    (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(4)) - &
                2*E**(5*y)*shat*smT**3* &
                 (2*MV2**2*shat**2*(22*SP(3) - 5*SP(4)) + &
                   3*MV2*shat*smT2*(-18*SP(3) + SP(4)) + &
                   smT2**2*(22*SP(3) + SP(4))) + &
                E**(6*y)*smT2* &
                 (64*MV2**3*shat**3*SP(3) + &
                   3*smT2**3*(4*SP(3) + SP(4)) + &
                   MV2*shat*smT2**2*(16*SP(3) + SP(4)) - &
                   4*MV2**2*shat**2*smT2*(23*SP(3) + 2*SP(4))) + &
                E**(2*y)*smT2**2* &
                 (smT2**2*(4*SP(3) + SP(4)) + &
                   4*MV2**2*shat**2*(3*SP(3) + 4*SP(4)) - &
                   MV2*shat*smT2*(20*SP(3) + 13*SP(4)))))/ &
            (E**(4*y)*shat*smT**2*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*MV2*smT2*LC(3)* &
              (E**(9*y)*shat*smT**3*(3*MV2*shat - smT2)*smT2*SP(3) + &
                E**(10*y)*smT2**3*(-2*MV2*shat + smT2)*SP(3) + &
                smT2**3*(-2*MV2*shat + smT2)*SP(4) - &
                2*E**(3*y)*shat*smT**3* &
                 ((4*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
                    SP(3) + MV2*shat*(7*MV2*shat - 4*smT2)*SP(4)) + &
                E**(4*y)*smT2* &
                 (2*(-5*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 - &
                      7*MV2*shat*smT2**2 + 2*smT2**3)*SP(3) - &
                   3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
                    (2*MV2*shat - smT2)*SP(4)) + &
                2*E**(5*y)*shat*smT* &
                 ((2*MV2*shat - smT2)* &
                    (7*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(3) + MV2*shat* &
                    (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
                    SP(4)) + &
                E**(6*y)*smT2* &
                 (-3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
                    (2*MV2*shat - smT2)*SP(3) + &
                   2*(-5*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 - &
                      7*MV2*shat*smT2**2 + 2*smT2**3)*SP(4)) - &
                2*E**(7*y)*shat*smT*smT2* &
                 (11*MV2**2*shat**2*SP(3) + 2*smT2**2*SP(3) + &
                   MV2*shat*smT2*(-10*SP(3) + SP(4))) + &
                E**(8*y)*smT2**2* &
                 (smT2**2*(4*SP(3) + SP(4)) + &
                   MV2**2*shat**2*(14*SP(3) + SP(4)) - &
                   MV2*shat*smT2*(16*SP(3) + SP(4))) + &
                E**y*shat*smT**5* &
                 (-(smT2*SP(3)) + MV2*shat*(SP(3) + 2*SP(4))) + &
                E**(2*y)*smT2**2* &
                 (smT2**2*(SP(3) + 4*SP(4)) + &
                   MV2**2*shat**2*(SP(3) + 14*SP(4)) - &
                   MV2*shat*smT2*(SP(3) + 16*SP(4)))))/ &
            (E**(5*y)*shat**2*smT**3*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (MV2*LC(13)*(E**(10*y)*smT**5*(MV2*shat - smT2)*SP(3) - &
                E**(11*y)*shat*(MV2*shat - smT2)*smT2**2*SP(3) - &
                4*MV2*shat*smT**5*SP(4) - &
                12*E**(3*y)*shat*(5*MV2*shat - 3*smT2)*smT2**2*SP(4) + &
                8*E**y*shat*smT2**3*SP(4) + &
                4*E**(2*y)*smT**3* &
                 (8*MV2**2*shat**2 - 5*MV2*shat*smT2 - smT2**2)*SP(4) + &
                E**(9*y)*shat*(MV2*shat - smT2)*smT2* &
                 (8*MV2*shat*SP(3) - 5*smT2*SP(3) + 4*smT2*SP(4)) - &
                E**(7*y)*shat* &
                 ((MV2*shat - smT2)* &
                    (16*MV2**2*shat**2 - 20*MV2*shat*smT2 + 7*smT2**2)* &
                    SP(3) - 4*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
                    smT2*SP(4)) - &
                E**(8*y)*smT**3* &
                 (3*(MV2*shat - smT2)*(2*MV2*shat - smT2)*SP(3) + &
                   4*(3*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
                    SP(4)) - &
                E**(4*y)*smT*(MV2*shat - smT2)* &
                 (2*MV2*shat*smT2*(SP(3) - 14*SP(4)) + &
                   64*MV2**2*shat**2*SP(4) - smT2**2*(SP(3) + 12*SP(4))) &
                  + E**(6*y)*smT* &
                 (8*MV2**3*shat**3*SP(3) - &
                   16*MV2**2*shat**2*smT2*(SP(3) + SP(4)) - &
                   3*smT2**3*(SP(3) + 4*SP(4)) + &
                   MV2*shat*smT2**2*(11*SP(3) + 20*SP(4))) + &
                E**(5*y)*shat*smT2* &
                 (4*MV2**2*shat**2*(SP(3) + 22*SP(4)) + &
                   smT2**2*(3*SP(3) + 44*SP(4)) - &
                   MV2*shat*smT2*(7*SP(3) + 108*SP(4)))))/ &
            (E**(6*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
           (4*MV2*(-(E**y*shat) + smT)*LC(5)*SP(5))/ &
            (shat*smT*(-(MV2*shat) + smT2)) + &
           (MV2*(-2*MV2*smT**5 + 2*E**(8*y)*MV2*smT**5 + &
                2*E**(2*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) - &
                2*E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) + &
                E**(9*y)*smT2**2*(-(MV2*shat) + smT2) + &
                2*E**y*smT2**2*(MV2*shat + smT2) + &
                E**(7*y)*smT2* &
                 (8*MV2**2*shat**2 - 15*MV2*shat*smT2 + 3*smT2**2) + &
                E**(5*y)*(-4*MV2*shat + smT2)* &
                 (4*MV2**2*shat**2 - 13*MV2*shat*smT2 + 5*smT2**2) + &
                E**(3*y)*smT2* &
                 (-16*MV2**2*shat**2 - MV2*shat*smT2 + 5*smT2**2))* &
              LC(4)*SP(6))/ &
            (E**(4*y)*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-(MV2*smT**5) + E**(8*y)*MV2*smT**5 - &
                E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) + &
                6*E**(4*y)*MV2*smT*(-(MV2*shat) + smT2)**2 - &
                E**(7*y)*smT2**2*(MV2*shat + smT2) + &
                E**(5*y)*smT2* &
                 (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
                E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
                   smT2**2) + &
                E**(2*y)*MV2*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
                E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(10)*SP(6))/ &
            (E**(3*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-2*E**(4*y)*MV2**2*shat*(8*MV2*shat - 5*smT2)* &
                 (MV2*shat - smT2) + &
                E**(6*y)*MV2*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
                 smT2 - MV2**2*shat*smT2**2 - &
                E**(8*y)*MV2*(MV2*shat - 2*smT2)*smT2**2 + &
                E**y*smT**5*(MV2*shat + smT2) - &
                E**(3*y)*smT**3* &
                 (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) + &
                E**(2*y)*MV2*smT2* &
                 (8*MV2**2*shat**2 - 5*MV2*shat*smT2 - smT2**2) + &
                E**(7*y)*smT**3* &
                 (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(5*y)*smT* &
                 (4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(14)*SP(6))/ &
            (E**(5*y)*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (MV2*(MV2*shat*smT**5 + smT**7 + &
                E**(8*y)*(-3*MV2*shat*smT**5 + smT**7) + &
                E**(6*y)*(22*MV2**2*shat**2*smT**3 - &
                   22*MV2*shat*smT**5 + 4*smT**7) + &
                2*E**(9*y)*shat*(MV2*shat - smT2)*smT2**2 + &
                2*E**(2*y)*smT**3*(MV2*shat + smT2)* &
                 (-3*MV2*shat + 2*smT2) - &
                E**y*shat*smT2**2*(MV2*shat + 3*smT2) + &
                E**(3*y)*shat*smT2* &
                 (8*MV2**2*shat**2 + 13*MV2*shat*smT2 - 9*smT2**2) + &
                2*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
                 (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2) - &
                E**(7*y)*shat*smT2* &
                 (16*MV2**2*shat**2 - 27*MV2*shat*smT2 + 7*smT2**2) + &
                E**(5*y)*shat* &
                 (32*MV2**3*shat**3 - 88*MV2**2*shat**2*smT2 + &
                   55*MV2*shat*smT2**2 - 11*smT2**3))*LC(2)*SP(7))/ &
            (E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(MV2*smT**5 - E**(8*y)*MV2*smT**5 + &
                E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) - &
                6*E**(4*y)*MV2*smT*(-(MV2*shat) + smT2)**2 + &
                E**(7*y)*smT2**2*(MV2*shat + smT2) + &
                E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
                   smT2**2) - &
                E**(2*y)*MV2*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(5*y)*smT2* &
                 (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
                E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(8)*SP(7))/ &
            (E**(3*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(2*E**(4*y)*MV2**2*shat*(8*MV2*shat - 5*smT2)* &
                 (MV2*shat - smT2) - &
                E**(6*y)*MV2*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
                 smT2 + MV2**2*shat*smT2**2 + &
                E**(8*y)*MV2*(MV2*shat - 2*smT2)*smT2**2 - &
                E**y*smT**5*(MV2*shat + smT2) + &
                E**(3*y)*smT**3* &
                 (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
                E**(7*y)*smT**3* &
                 (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(2*y)*MV2*smT2* &
                 (-8*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2) - &
                E**(5*y)*smT* &
                 (4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(12)*SP(7))/ &
            (E**(5*y)*smT*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
        kap*SI(3)*(((-2*shat*smT**7*(5*MV2*shat + smT2) + &
                3*E**(12*y)*shat*smT**7*(MV2*shat + 3*smT2) + &
                E**(11*y)*smT2**3* &
                 (MV2*(13*MV2 - 6*shat)*shat**2 - &
                   2*shat*(23*MV2 + 9*shat)*smT2 + 9*smT2**2) + &
                E**(10*y)*smT**5* &
                 (MV2*shat**2*(-16*MV2**2 - 41*MV2*shat + 3*shat**2) + &
                   shat*(37*MV2**2 + 35*MV2*shat + 9*shat**2)*smT2 + &
                   3*(-3*MV2 + 10*shat)*smT2**2) + &
                2*E**y*smT2**2* &
                 (-8*MV2**2*shat**4 + &
                   MV2*shat**2*(5*MV2 + 26*shat)*smT2 + &
                   2*(2*MV2 - 3*shat)*shat*smT2**2 + 3*smT2**3) + &
                2*E**(9*y)*smT2**2* &
                 (3*MV2**2*shat**3*(-15*MV2 + 7*shat) + &
                   MV2*shat**2*(184*MV2 + 33*shat)*smT2 - &
                   2*shat*(74*MV2 + 21*shat)*smT2**2 + 21*smT2**3) + &
                E**(8*y)*smT**3* &
                 (2*MV2**2*shat**3* &
                    (64*MV2**2 + 77*MV2*shat - 9*shat**2) - &
                   MV2*shat**2*(358*MV2**2 + 595*MV2*shat + 47*shat**2)* &
                    smT2 + shat* &
                    (239*MV2**2 + 324*MV2*shat + 41*shat**2)*smT2**2 + &
                   33*(-MV2 + shat)*smT2**3) + &
                E**(2*y)*smT**3* &
                 (16*MV2**2*shat**4*(8*MV2 + shat) - &
                   2*MV2*shat**3*(94*MV2 + 21*shat)*smT2 + &
                   shat*(-6*MV2**2 + 45*MV2*shat + 14*shat**2)* &
                    smT2**2 - 3*(2*MV2 + 3*shat)*smT2**3) + &
                E**(6*y)*smT*(-(MV2*shat) + smT2)* &
                 (16*MV2**2*shat**3* &
                    (16*MV2**2 + 32*MV2*shat + 3*shat**2) - &
                   32*MV2*shat**2*(16*MV2**2 + 34*MV2*shat + 3*shat**2)* &
                    smT2 + shat* &
                    (304*MV2**2 + 607*MV2*shat + 69*shat**2)*smT2**2 + &
                   (-45*MV2 + 7*shat)*smT2**3) + &
                E**(3*y)*smT2* &
                 (-8*MV2**3*shat**4*(17*MV2 + 15*shat) + &
                   6*MV2**2*shat**3*(11*MV2 + shat)*smT2 + &
                   MV2*shat**2*(139*MV2 + 156*shat)*smT2**2 - &
                   6*shat*(21*MV2 + 11*shat)*smT2**3 + 33*smT2**4) + &
                2*E**(7*y)*smT2* &
                 (16*MV2**3*shat**4*(5*MV2 + shat) - &
                   MV2**2*shat**3*(301*MV2 + 15*shat)*smT2 + &
                   MV2*shat**2*(504*MV2 + 101*shat)*smT2**2 - &
                   2*shat*(149*MV2 + 39*shat)*smT2**3 + 39*smT2**4) - &
                E**(4*y)*smT* &
                 (8*MV2**3*(19*MV2 - shat)*shat**5 - &
                   2*MV2**2*shat**3* &
                    (68*MV2**2 + 459*MV2*shat + 61*shat**2)*smT2 + &
                   MV2*shat**2* &
                    (226*MV2**2 + 1033*MV2*shat + 157*shat**2)*smT2**2 &
                    - 3*shat*(47*MV2**2 + 121*MV2*shat + 17*shat**2)* &
                    smT2**3 + 3*(9*MV2 + 4*shat)*smT2**4) + &
                2*E**(5*y)*(16*MV2**4*shat**5*(8*MV2 + 5*shat) - &
                   4*MV2**3*shat**4*(69*MV2 + 61*shat)*smT2 + &
                   MV2**2*shat**3*(-71*MV2 + 95*shat)*smT2**2 + &
                   3*MV2*shat**2*(133*MV2 + 39*shat)*smT2**3 - &
                   24*shat*(10*MV2 + 3*shat)*smT2**4 + 36*smT2**5))* &
              LC(1))/ &
            (2d0*E**(4*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)* &
              (-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           ((-2*E**(13*y)*shat*smT2**4*(MV2*shat + 2*smT2) - &
                shat**2*smT**7*(MV2*shat + 5*smT2) + &
                E**y*shat*smT2**3* &
                 (MV2*shat**2*(-15*MV2 + 2*shat) + &
                   10*shat*(3*MV2 + shat)*smT2 + 3*smT2**2) + &
                E**(12*y)*smT**5* &
                 (8*MV2**2*(MV2 - 3*shat)*shat**3 + &
                   3*MV2*shat**2*(-2*MV2 + 17*shat)*smT2 + &
                   (5*MV2 - 9*shat)*shat*smT2**2 + 5*smT2**3) + &
                E**(2*y)*smT**5* &
                 (MV2*shat**3*(16*MV2**2 + 27*MV2*shat - shat**2) - &
                   shat**2*(MV2**2 + 15*MV2*shat + 5*shat**2)*smT2 - &
                   shat*(43*MV2 + 36*shat)*smT2**2 + 10*smT2**3) + &
                E**(3*y)*smT2**2* &
                 (MV2**2*(104*MV2 - 21*shat)*shat**4 - &
                   MV2*shat**2*(24*MV2**2 + 308*MV2*shat + 43*shat**2)* &
                    smT2 + shat*(40*MV2**2 + 199*MV2*shat + 58*shat**2)* &
                    smT2**2 + 5*(-2*MV2 + shat)*smT2**3) + &
                E**(11*y)*smT2**2* &
                 (8*MV2**2*shat**4*(16*MV2 + 5*shat) - &
                   8*MV2*shat**2*(MV2**2 + 30*MV2*shat + 10*shat**2)* &
                    smT2 + shat*(7*MV2**2 + 107*MV2*shat + 22*shat**2)* &
                    smT2**2 - (5*MV2 + 19*shat)*smT2**3) + &
                E**(10*y)*smT**3* &
                 (-8*MV2**2*shat**4* &
                    (25*MV2**2 + 23*MV2*shat + 2*shat**2) + &
                   MV2*shat**3*(224*MV2**2 + 201*MV2*shat + 31*shat**2)* &
                    smT2 + shat**2* &
                    (85*MV2**2 + 39*MV2*shat - 9*shat**2)*smT2**2 - &
                   shat*(145*MV2 + 56*shat)*smT2**3 + 30*smT2**4) + &
                E**(4*y)*smT**3* &
                 (-2*MV2**2*shat**4* &
                    (64*MV2**2 + 82*MV2*shat - 5*shat**2) + &
                   2*MV2*shat**3* &
                    (48*MV2**2 + 205*MV2*shat + 13*shat**2)*smT2 + &
                   shat**2*(298*MV2**2 - 85*MV2*shat - 24*shat**2)* &
                    smT2**2 - shat*(287*MV2 + 107*shat)*smT2**3 + &
                   45*smT2**4) - &
                E**(5*y)*smT2* &
                 (4*MV2**3*(28*MV2 - 15*shat)*shat**5 + &
                   MV2**2*shat**3* &
                    (-192*MV2**2 - 720*MV2*shat + 83*shat**2)*smT2 + &
                   MV2*shat**2* &
                    (416*MV2**2 + 1326*MV2*shat + 149*shat**2)*smT2**2 &
                    - shat*(247*MV2**2 + 662*MV2*shat + 136*shat**2)* &
                    smT2**3 + 5*(7*MV2 + 2*shat)*smT2**4) - &
                E**(9*y)*smT2* &
                 (16*MV2**3*(MV2 - 4*shat)*shat**5 - &
                   MV2**2*shat**3* &
                    (200*MV2**2 + 780*MV2*shat + 63*shat**2)*smT2 + &
                   MV2*shat**2* &
                    (360*MV2**2 + 1267*MV2*shat + 197*shat**2)*smT2**2 &
                    - shat*(197*MV2**2 + 594*MV2*shat + 94*shat**2)* &
                    smT2**3 + (25*MV2 + 37*shat)*smT2**4) + &
                2*E**(6*y)*smT* &
                 (4*MV2**3*shat**5* &
                    (32*MV2**2 + 36*MV2*shat - 3*shat**2) - &
                   2*MV2**3*shat**4*(116*MV2 + 257*shat)*smT2 + &
                   3*MV2*shat**3* &
                    (-26*MV2**2 + 190*MV2*shat + 11*shat**2)*smT2**2 + &
                   shat**2*(454*MV2**2 - 106*MV2*shat - 21*shat**2)* &
                    smT2**3 - 2*shat*(150*MV2 + 41*shat)*smT2**4 + &
                   40*smT2**5) + &
                E**(8*y)*smT* &
                 (8*MV2**3*shat**5* &
                    (48*MV2**2 + 23*MV2*shat - shat**2) - &
                   2*MV2**2*shat**4* &
                    (496*MV2**2 + 484*MV2*shat + 21*shat**2)*smT2 + &
                   2*MV2*shat**3*(14*MV2 + shat)*(13*MV2 + 35*shat)* &
                    smT2**2 + &
                   shat**2*(644*MV2**2 - 153*MV2*shat - 32*shat**2)* &
                    smT2**3 - shat*(506*MV2 + 135*shat)*smT2**4 + &
                   70*smT2**5) - &
                E**(7*y)*(32*MV2**4*shat**6*(8*MV2 + shat) + &
                   4*MV2**3*shat**4* &
                    (96*MV2**2 - 14*MV2*shat - 53*shat**2)*smT2 + &
                   MV2**2*shat**3* &
                    (-1024*MV2**2 - 1556*MV2*shat + 95*shat**2)*smT2**2 &
                    + MV2*shat**2* &
                    (992*MV2**2 + 2276*MV2*shat + 221*shat**2)*smT2**3 &
                    - shat*(397*MV2**2 + 978*MV2*shat + 160*shat**2)* &
                    smT2**4 + (45*MV2 + 34*shat)*smT2**5))*LC(11))/ &
            (2d0*E**(6*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)**2* &
              (-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           ((-(shat*smT2**3*(3*MV2*shat + 7*smT2)) + &
                E**(12*y)*shat*smT2**3*(MV2*shat + 9*smT2) - &
                E**(11*y)*smT**5* &
                 (MV2*shat**2*(5*MV2 + 3*shat) + &
                   shat*(7*MV2 + 17*shat)*smT2 + 8*smT2**2) + &
                E**y*smT**5*(3*MV2*shat**2*(MV2 + 2*shat) + &
                   shat*(5*MV2 + 14*shat)*smT2 + 12*smT2**2) + &
                E**(10*y)*smT2**2* &
                 (MV2*shat**3*(-19*MV2 + 2*shat) + &
                   shat*(6*MV2**2 - MV2*shat + 8*shat**2)*smT2 + &
                   4*(MV2 + 10*shat)*smT2**2) - &
                E**(8*y)*shat*smT2* &
                 (2*MV2**2*shat**3*(41*MV2 + 18*shat) + &
                   MV2*shat*(70*MV2**2 + 68*MV2*shat + shat**2)*smT2 - &
                   (50*MV2**2 + 32*MV2*shat + 17*shat**2)*smT2**2 - &
                   48*smT2**3) + &
                E**(9*y)*smT**3* &
                 (MV2**2*shat**3*(62*MV2 + 55*shat) + &
                   MV2*shat**2*(-12*MV2 + 19*shat)*smT2 - &
                   2*shat*(5*MV2 + 27*shat)*smT2**2 - 20*smT2**3) + &
                E**(3*y)*smT**3* &
                 (-(MV2**2*shat**3*(46*MV2 + 79*shat)) + &
                   2*MV2*shat**2*(9*MV2 + 10*shat)*smT2 + &
                   shat*(-32*MV2 + 39*shat)*smT2**2 + 40*smT2**3) - &
                E**(6*y)*(-(MV2*shat) + smT2)* &
                 (-8*MV2**2*shat**4*(4*MV2 + shat) + &
                   8*MV2*shat**2*(8*MV2**2 + 11*MV2*shat + shat**2)* &
                    smT2 - shat*(16*MV2 + shat)*(4*MV2 + 3*shat)* &
                    smT2**2 + (24*MV2 + 13*shat)*smT2**3) + &
                E**(4*y)*smT2* &
                 (2*MV2**2*shat**4*(41*MV2 + 22*shat) + &
                   MV2*shat**2*(6*MV2**2 + 28*MV2*shat - 11*shat**2)* &
                    smT2 + shat*(46*MV2**2 + 22*MV2*shat - 13*shat**2)* &
                    smT2**2 - 2*(16*MV2 + 31*shat)*smT2**3) + &
                E**(7*y)*shat*smT* &
                 (4*MV2**3*shat**3*(-12*MV2 + 17*shat) + &
                   MV2**2*(170*MV2 - 37*shat)*shat**2*smT2 + &
                   MV2*shat*(-41*MV2 + 55*shat)*smT2**2 - &
                   (41*MV2 + 46*shat)*smT2**3) + &
                E**(5*y)*smT* &
                 (4*MV2**3*(4*MV2 - 5*shat)*shat**4 - &
                   MV2**2*shat**3*(154*MV2 + 83*shat)*smT2 + &
                   MV2*shat**2*(133*MV2 + 47*shat)*smT2**2 + &
                   shat*(-75*MV2 + 16*shat)*smT2**3 + 40*smT2**4) - &
                E**(2*y)*smT2**2* &
                 (3*MV2*(shat**2 + 2*smT2)**2 + &
                   shat*smT2*(7*shat**2 + 39*smT2) - &
                   MV2**2*(31*shat**3 + 2*shat*smT2)))*LC(15))/ &
            (2d0*E**(5*y)*(E**y*shat - smT)*(-shat + E**y*smT)* &
              (-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(7)*SP(1))/ &
            (shat*(-shat + E**y*smT)*(MV2*shat - smT2)) + &
           (8*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(6)*SP(2))/ &
            (shat*(shat - E**y*smT)*(MV2*shat - smT2)) + &
           (LC(9)*(4*E**(15*y)*shat**2*smT2**4*SP(3) - &
                4*E**(14*y)*shat*smT**7*(2*shat**2 + 3*smT2)*SP(3) + &
                2*shat*smT**7*(2*MV2*shat - smT2)*SP(4) - &
                E**y*smT2**3* &
                 (MV2*shat**2*(-MH2 + 4*MV2 + 11*shat) + &
                   (MH2 - MV2 - 5*shat)*shat*smT2 + smT2**2)*SP(4) - &
                2*E**(13*y)*smT2**3* &
                 (2*(shat**2*(4*MV2**2 + 7*MV2*shat - shat**2) - &
                      shat*(8*MV2 + 11*shat)*smT2 + smT2**2)*SP(3) - &
                   MV2*shat**3*SP(4)) + &
                E**(12*y)*smT**5* &
                 (4*(MV2*shat**2*(3*MV2**2 + 5*MV2*shat + 14*shat**2) + &
                      shat*(-5*MV2**2 + 9*MV2*shat - 13*shat**2)*smT2 + &
                      (MV2 - 14*shat)*smT2**2)*SP(3) - &
                   shat**2*(MV2*shat*(-MH2 + 2*MV2 + 3*shat) + &
                      (MH2 + 4*MV2 + shat)*smT2)*SP(4)) + &
                E**(2*y)*smT**5* &
                 (4*MV2*shat**2*smT2*SP(3) + &
                   (-(MV2*(MH2 + 21*MV2 - 10*shat)*shat**3) + &
                      shat*(MV2**2 + 36*MV2*shat + (MH2 - 4*shat)*shat)* &
                       smT2 + (MV2 - 7*shat)*smT2**2)*SP(4)) + &
                E**(11*y)*smT2**2* &
                 (4*(MV2*shat**3*(29*MV2**2 + 9*MV2*shat - 7*shat**2) + &
                      shat**2*(-75*MV2**2 - 57*MV2*shat + 5*shat**2)* &
                       smT2 + shat*(46*MV2 + 39*shat)*smT2**2 - &
                      4*smT2**3)*SP(3) + &
                   shat*(MV2*shat**3*(-14*MV2 + shat) + &
                      shat*(-(MH2*MV2) + 4*MV2**2 + 21*MV2*shat + &
                         shat**2)*smT2 + (MH2 + 2*MV2 + shat)*smT2**2)* &
                    SP(4)) + &
                E**(3*y)*smT2**2* &
                 (4*shat*(5*MV2**2*shat**3 - &
                      MV2*shat*(MV2 + 12*shat)*smT2 - &
                      (MV2 - 4*shat)*smT2**2)*SP(3) + &
                   (MV2*shat**3* &
                       (-8*(MH2 - 4*MV2)*MV2 + 76*MV2*shat - 3*shat**2) &
                       + shat**2* &
                       (11*MH2*MV2 - 34*MV2**2 - 97*MV2*shat + shat**2)* &
                       smT2 + shat*(-3*MH2 + 8*MV2 + 21*shat)*smT2**2 - &
                      4*smT2**3)*SP(4)) + &
                E**(10*y)*smT**3* &
                 (-4*(MV2**2*shat**3* &
                       (24*MV2**2 + 41*MV2*shat + 19*shat**2) - &
                      2*MV2*shat**2* &
                       (26*MV2**2 + 43*MV2*shat + 29*shat**2)*smT2 + &
                      shat*(29*MV2**2 + 7*MV2*shat + 33*shat**2)* &
                       smT2**2 - 3*(MV2 - 9*shat)*smT2**3)*SP(3) + &
                   shat*(MV2**2*shat**3*(-8*MH2 + 16*MV2 + 19*shat) + &
                      MV2*(11*MH2 + 15*MV2 - 9*shat)*shat**2*smT2 - &
                      (2*MV2**2 + 32*MV2*shat + shat*(3*MH2 + 8*shat))* &
                       smT2**2 + smT2**3)*SP(4)) - &
                E**(4*y)*smT**3* &
                 (4*shat*(2*MV2**2*shat**3*(13*MV2 + 4*shat) - &
                      MV2*shat**2*(41*MV2 + 17*shat)*smT2 + &
                      (-MV2**2 + 8*MV2*shat + 6*shat**2)*smT2**2 + &
                      3*smT2**3)*SP(3) + &
                   (MV2**2*shat**4*(-8*MH2 + 36*MV2 + 73*shat) + &
                      MV2*shat**2* &
                       (8*MV2**2 + 30*MV2*shat + &
                         (11*MH2 - 82*shat)*shat)*smT2 + &
                      shat*(-MV2**2 - 56*MV2*shat + &
                         shat*(-3*MH2 + 17*shat))*smT2**2 + &
                      (-3*MV2 + 8*shat)*smT2**3)*SP(4)) + &
                E**(9*y)*smT2* &
                 (-4*(MV2**2*shat**4* &
                       (40*MV2**2 - 5*MV2*shat - 9*shat**2) - &
                      2*MV2*shat**3* &
                       (71*MV2**2 + 19*MV2*shat - 9*shat**2)*smT2 + &
                      shat**2*(189*MV2**2 + 112*MV2*shat - 9*shat**2)* &
                       smT2**2 - shat*(89*MV2 + 65*shat)*smT2**3 + &
                      6*smT2**4)*SP(3) + &
                   (2*MV2**2*(2*MV2 - 3*shat)*shat**5 - &
                      2*MV2*shat**3* &
                       (-4*(MH2 - 4*MV2)*MV2 + 39*MV2*shat + shat**2)* &
                       smT2 + &
                      shat**2* &
                       (-11*MH2*MV2 + 12*MV2**2 + 47*MV2*shat + &
                         4*shat**2)*smT2**2 + &
                      shat*(3*MH2 + 13*MV2 + 9*shat)*smT2**3 - smT2**4)* &
                    SP(4)) - &
                2*E**(5*y)*smT2* &
                 (2*(-3*MV2**2*shat**4*(MV2 + shat)*(8*MV2 + shat) + &
                      2*MV2*shat**3*(7*MV2**2 + 3*MV2*shat + 3*shat**2)* &
                       smT2 + &
                      shat**2*(27*MV2**2 + 45*MV2*shat - 2*shat**2)* &
                       smT2**2 - shat*(19*MV2 + 24*shat)*smT2**3 + &
                      smT2**4)*SP(3) + &
                   (MV2**2*shat**4* &
                       (-8*MH2*MV2 + 20*MV2**2 + 20*MV2*shat - &
                         11*shat**2) + &
                      MV2*shat**3* &
                       (12*MH2*MV2 - 58*MV2**2 - 95*MV2*shat + &
                         11*shat**2)*smT2 + &
                      shat**2* &
                       (-5*(MH2 - 8*MV2)*MV2 + 81*MV2*shat - 2*shat**2)* &
                       smT2**2 + &
                      (MH2 - 11*MV2 - 17*shat)*shat*smT2**3 + 3*smT2**4) &
                     *SP(4)) + &
                E**(6*y)*smT* &
                 (4*(4*MV2**3*(8*MV2 - shat)*shat**5 - &
                      MV2**2*shat**3* &
                       (24*MV2**2 + 175*MV2*shat + 37*shat**2)*smT2 + &
                      2*MV2*shat**2* &
                       (20*MV2**2 + 97*MV2*shat + 31*shat**2)*smT2**2 - &
                      shat*(19*MV2**2 + 46*MV2*shat + 25*shat**2)* &
                       smT2**3 + (MV2 - 14*shat)*smT2**4)*SP(3) + &
                   (4*MV2**3*shat**5*(-4*MH2 + 20*MV2 + 15*shat) - &
                      MV2**2*shat**3* &
                       (8*MV2**2 + 188*MV2*shat + &
                         3*shat*(-8*MH2 + 53*shat))*smT2 + &
                      MV2*shat**2* &
                       (16*MV2**2 + 118*MV2*shat + &
                         shat*(-10*MH2 + 119*shat))*smT2**2 - &
                      shat*(11*MV2**2 + 16*MV2*shat - &
                         2*(MH2 - 14*shat)*shat)*smT2**3 + &
                      (3*MV2 - 2*shat)*smT2**4)*SP(4)) + &
                E**(8*y)*smT* &
                 (4*(4*MV2**3*shat**4* &
                       (12*MV2**2 + 24*MV2*shat + shat**2) - &
                      2*MV2**2*shat**3* &
                       (64*MV2**2 + 147*MV2*shat + 28*shat**2)*smT2 + &
                      MV2*shat**2* &
                       (121*MV2**2 + 274*MV2*shat + 89*shat**2)*smT2**2 &
                       - shat*(44*MV2**2 + 53*MV2*shat + 41*shat**2)* &
                       smT2**3 + 3*(MV2 - 9*shat)*smT2**4)*SP(3) + &
                   (4*MV2**3*shat**5*(4*MH2 - 2*MV2 + 5*shat) + &
                      MV2**2*(-24*MH2 + 8*MV2 - 27*shat)*shat**4*smT2 + &
                      MV2*shat**2* &
                       (16*MV2**2 + 88*MV2*shat + &
                         shat*(10*MH2 + 41*shat))*smT2**2 - &
                      shat*(13*MV2**2 + 68*MV2*shat + &
                         2*shat*(MH2 + 11*shat))*smT2**3 + &
                      (MV2 + 2*shat)*smT2**4)*SP(4)) - &
                2*E**(7*y)*(2* &
                    (16*MV2**4*shat**5*(3*MV2 + 2*shat) - &
                      4*MV2**2*shat**4* &
                       (24*MV2**2 + 26*MV2*shat + 3*shat**2)*smT2 + &
                      MV2*shat**3* &
                       (-27*MV2**2 + 22*MV2*shat + 17*shat**2)*smT2**2 &
                       + shat**2*(136*MV2**2 + 95*MV2*shat - 7*shat**2)* &
                       smT2**3 - shat*(71*MV2 + 56*shat)*smT2**4 + &
                      4*smT2**5)*SP(3) + &
                   (4*MV2**3*shat**6*(5*MV2 + 2*shat) - &
                      2*MV2**2*shat**4* &
                       (-4*MH2*MV2 + 4*MV2**2 + 15*MV2*shat + &
                         8*shat**2)*smT2 + &
                      MV2*shat**3* &
                       (-12*MH2*MV2 + 14*MV2**2 + 3*MV2*shat + &
                         11*shat**2)*smT2**2 + &
                      shat**2* &
                       (5*MH2*MV2 + 9*MV2**2 + 24*MV2*shat - 3*shat**2)* &
                       smT2**3 - shat*(MH2 + 13*(MV2 + shat))*smT2**4 + &
                      2*smT2**5)*SP(4))))/ &
            (E**(5*y)*shat*(-(E**y*shat) + smT)**2*(-shat + E**y*smT)* &
              (-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (LC(13)*(-(E**(13*y)*shat*(MV2*shat - smT2)*smT2**3*SP(3)) - &
                4*shat**2*smT**7*SP(4) + &
                4*E**y*shat*smT2**3*(shat**2 + 3*smT2)*SP(4) + &
                4*E**(2*y)*smT**5* &
                 (MV2*shat**2*(4*MV2 + 7*shat) - &
                   8*shat*(MV2 + shat)*smT2 + smT2**2)*SP(4) + &
                E**(12*y)*smT**5* &
                 ((MV2*shat - smT2)*(3*MV2*shat + 2*shat**2 - smT2)* &
                    SP(3) + 4*MV2*shat*smT2*SP(4)) - &
                E**(4*y)*smT**3* &
                 ((MV2*shat - smT2)* &
                    (smT2**2 + MV2*shat*(2*shat**2 + smT2))*SP(3) + &
                   4*(2*MV2**2*shat**3*(16*MV2 + 5*shat) - &
                      40*MV2*shat**2*(2*MV2 + shat)*smT2 + &
                      shat*(47*MV2 + 24*shat)*smT2**2 - 4*smT2**3)*SP(4) &
                   ) - 2*E**(10*y)*smT**3* &
                 ((MV2*shat - smT2)* &
                    (MV2*shat**2*(6*MV2 + 5*shat) - &
                      3*shat*(2*MV2 + shat)*smT2 + 2*smT2**2)*SP(3) + &
                   2*(2*MV2**2*shat**3*(13*MV2 + 4*shat) - &
                      MV2*shat**2*(46*MV2 + 17*shat)*smT2 + &
                      2*shat*(10*MV2 + 3*shat)*smT2**2 - smT2**3)*SP(4)) &
                  + E**(9*y)*smT2* &
                 ((MV2*shat - smT2)* &
                    (2*MV2*shat**2*(MV2 + shat)*(4*MV2 + 3*shat) - &
                      shat*(8*MV2**2 + 7*MV2*shat + 3*shat**2)*smT2 + &
                      (3*MV2 + shat)*smT2**2)*SP(3) + &
                   4*(3*MV2**2*shat**3*(MV2 + shat)*(8*MV2 + shat) - &
                      2*MV2*shat**2* &
                       (20*MV2**2 + 7*MV2*shat + 3*shat**2)*smT2 + &
                      shat*(19*MV2**2 - 28*MV2*shat + 2*shat**2)* &
                       smT2**2 - (MV2 - 18*shat)*smT2**3)*SP(4)) + &
                E**(5*y)*smT2* &
                 ((MV2*shat - smT2)* &
                    (2*MV2*shat**3*(MV2 + shat) - &
                      shat**2*(5*MV2 + shat)*smT2 + &
                      (MV2 + 4*shat)*smT2**2)*SP(3) + &
                   4*(3*MV2**2*shat**3* &
                       (8*MV2**2 + 3*MV2*shat + 3*shat**2) - &
                      2*MV2*shat**2* &
                       (26*MV2**2 + 3*MV2*shat + 9*shat**2)*smT2 + &
                      shat*(29*MV2**2 - 40*MV2*shat + 9*shat**2)* &
                       smT2**2 + (-3*MV2 + 31*shat)*smT2**3)*SP(4)) - &
                2*E**(8*y)*smT* &
                 (-((shat - smT)*(shat + smT)*(4*MV2*shat - 3*smT2)* &
                      (-(MV2*shat) + smT2)**2*SP(3)) + &
                   2*(4*MV2**3*shat**4*(-8*MV2 + shat) + &
                      2*MV2**2*shat**3*(74*MV2 + 17*shat)*smT2 - &
                      4*MV2*shat**2*(45*MV2 + 14*shat)*smT2**2 + &
                      shat*(74*MV2 + 23*shat)*smT2**3 - 4*smT2**4)*SP(4) &
                   ) + 2*E**(6*y)*smT* &
                 ((MV2*shat - smT2)* &
                    (4*MV2**2*shat**4 + &
                      MV2*(2*MV2 - 3*shat)*shat**2*smT2 + &
                      shat*(2*MV2 + shat)*smT2**2 - 2*smT2**3)*SP(3) + &
                   4*(2*MV2**3*shat**4*(16*MV2 + shat) - &
                      MV2**2*shat**3*(97*MV2 + 22*shat)*smT2 + &
                      MV2*shat**2*(109*MV2 + 36*shat)*smT2**2 - &
                      shat*(46*MV2 + 17*shat)*smT2**3 + 3*smT2**4)*SP(4) &
                   ) - E**(7*y)* &
                 ((MV2*shat - smT2)* &
                    (8*MV2**2*shat**4*(MV2 + shat) - &
                      8*MV2*shat**3*(MV2 + shat)*smT2 + &
                      shat*(6*MV2**2 + 13*MV2*shat + 3*shat**2)* &
                       smT2**2 - (3*MV2 + 5*shat)*smT2**3)*SP(3) + &
                   4*(16*MV2**4*shat**4*(3*MV2 + 2*shat) - &
                      4*MV2**2*shat**3* &
                       (32*MV2**2 + 25*MV2*shat + 3*shat**2)*smT2 + &
                      MV2*shat**2* &
                       (121*MV2**2 + 56*MV2*shat + 17*shat**2)*smT2**2 &
                       - (11*MV2 - 7*shat)*(4*MV2 - shat)*shat* &
                       smT2**3 + 3*(MV2 - 11*shat)*smT2**4)*SP(4)) - &
                E**(11*y)*smT2**2* &
                 ((MV2*shat - smT2)* &
                    (2*MV2**2*shat + shat**3 + 2*shat*smT2 - &
                      MV2*(shat**2 + smT2))*SP(3) + &
                   4*shat*(12*MV2*shat*smT2 - 4*smT2**2 + &
                      MV2**2*(-5*shat**2 + smT2))*SP(4)) - &
                E**(3*y)*smT2**2* &
                 (shat*smT2*(-(MV2*shat) + smT2)*SP(3) + &
                   4*(3*MV2**3*shat**2 - &
                      5*shat*smT2*(shat**2 + 3*smT2) + &
                      MV2**2*(shat**3 - 5*shat*smT2) + &
                      MV2*(7*shat**4 + 17*shat**2*smT2 + smT2**2))*SP(4) &
                   )))/ &
            (E**(7*y)*shat*(shat - E**y*smT)**2*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*LC(3)*(-(E**(16*y)*shat**2*(2*MV2*shat - smT2)*smT2**4* &
                   SP(3)) + 3*E**(15*y)*shat*smT**7* &
                 (2*MV2*shat**2*(MV2 + shat) - &
                   shat*(2*MV2 + shat)*smT2 + smT2**2)*SP(3) - &
                MV2*shat**3*smT2**4*SP(4) + &
                E**y*shat*smT**7* &
                 (3*MV2*shat**2*(MV2 + shat) - 2*MV2*shat*smT2 + &
                   2*smT2**2)*SP(4) - &
                E**(14*y)*smT2**3* &
                 ((2*MV2*shat**3*(2*MV2**2 + MV2*shat + 3*shat**2) + &
                      shat**2*(-4*MV2**2 + 2*MV2*shat - 3*shat**2)* &
                       smT2 + shat*(2*MV2 + 3*shat)*smT2**2 + smT2**3)* &
                    SP(3) + shat**2*smT2*(-2*MV2*shat + smT2)*SP(4)) - &
                E**(2*y)*smT2**3* &
                 (shat**2*smT2*(-2*MV2*shat + smT2)*SP(3) + &
                   (MV2*shat**3*(2*MV2**2 + 2*MV2*shat + 3*shat**2) - &
                      MV2**2*shat**2*smT2 + &
                      shat*(MV2 + 5*shat)*smT2**2 + smT2**3)*SP(4)) + &
                E**(3*y)*smT**5* &
                 (shat*(MV2**2*shat**4 - &
                      MV2*shat**2*(2*MV2 + 7*shat)*smT2 + &
                      3*shat*(-MV2 + shat)*smT2**2 + 2*smT2**3)*SP(3) + &
                   (MV2*shat**4*(-18*MV2**2 - 11*MV2*shat + shat**2) + &
                      MV2*shat**2*(MV2**2 + 34*MV2*shat + 10*shat**2)* &
                       smT2 - &
                      shat*(MV2**2 + 25*MV2*shat - 4*shat**2)*smT2**2 + &
                      (MV2 + 12*shat)*smT2**3)*SP(4)) - &
                E**(4*y)*smT2**2* &
                 ((2*MV2**3*shat**5 + &
                      4*MV2*(2*MV2 - shat)*shat**4*smT2 + &
                      shat**2*(-4*MV2**2 - 28*MV2*shat + shat**2)* &
                       smT2**2 + 11*shat**2*smT2**3 + smT2**4)*SP(3) + &
                   (-2*MV2**2*shat**4* &
                       (8*MV2**2 + 27*MV2*shat + 8*shat**2) + &
                      MV2*shat**3* &
                       (22*MV2**2 + 75*MV2*shat + 12*shat**2)*smT2 + &
                      shat**2*(-10*MV2**2 - 55*MV2*shat + shat**2)* &
                       smT2**2 + shat*(-4*MV2 + 25*shat)*smT2**3 + &
                      5*smT2**4)*SP(4)) + &
                E**(11*y)*smT**3* &
                 ((2*MV2**2*shat**5* &
                       (9*MV2**2 - 26*MV2*shat - 7*shat**2) + &
                      MV2*shat**3* &
                       (-8*MV2**3 - 150*MV2**2*shat + 35*MV2*shat**2 + &
                         16*shat**3)*smT2 + &
                      shat**2* &
                       (14*MV2**3 + 212*MV2**2*shat - 3*MV2*shat**2 - &
                         4*shat**3)*smT2**2 + &
                      shat*(-12*MV2**2 - 135*MV2*shat + 5*shat**2)* &
                       smT2**3 + 4*(MV2 + 10*shat)*smT2**4)*SP(3) + &
                   (MV2**2*shat**5*(2*MV2**2 + 9*MV2*shat - shat**2) + &
                      2*MV2**2*shat**4*(5*MV2 + 6*shat)*smT2 + &
                      3*MV2*(MV2 - 12*shat)*shat**3*smT2**2 - &
                      2*shat*(MV2**2 + 15*MV2*shat - 6*shat**2)* &
                       smT2**3 + (MV2 + 12*shat)*smT2**4)*SP(4)) + &
                E**(5*y)*smT**3* &
                 ((MV2**2*shat**5*(2*MV2**2 - 7*MV2*shat - shat**2) + &
                      MV2*shat**4*(14*MV2**2 + 51*MV2*shat + shat**2)* &
                       smT2 - &
                      shat**3*(MV2**2 + 64*MV2*shat + shat**2)* &
                       smT2**2 + &
                      shat*(-2*MV2**2 - 31*MV2*shat + 17*shat**2)* &
                       smT2**3 + (MV2 + 13*shat)*smT2**4)*SP(3) + &
                   (2*MV2**2*shat**5* &
                       (MV2**2 - 22*MV2*shat - 3*shat**2) - &
                      MV2*(2*MV2 - shat)*shat**3* &
                       (4*MV2**2 + 45*MV2*shat + 4*shat**2)*smT2 + &
                      2*MV2*shat**2* &
                       (7*MV2**2 + 61*MV2*shat - 13*shat**2)*smT2**2 + &
                      shat*(-12*MV2**2 - 83*MV2*shat + 18*shat**2)* &
                       smT2**3 + 2*(2*MV2 + 15*shat)*smT2**4)*SP(4)) - &
                E**(10*y)*smT2* &
                 ((4*MV2**3*shat**5* &
                       (16*MV2**2 + 47*MV2*shat + 6*shat**2) - &
                      2*MV2**2*shat**4* &
                       (86*MV2**2 + 285*MV2*shat + 41*shat**2)*smT2 + &
                      2*MV2*shat**3* &
                       (69*MV2**2 + 290*MV2*shat + 33*shat**2)*smT2**2 &
                       - 2*shat**2* &
                       (17*MV2**2 + 140*MV2*shat + 7*shat**2)*smT2**3 + &
                      shat*(-12*MV2 + 67*shat)*smT2**4 + 10*smT2**5)* &
                    SP(3) + (4*MV2**3*shat**6*(7*MV2 + 5*shat) + &
                      2*MV2**2*shat**4* &
                       (2*MV2**2 - 11*MV2*shat - 5*shat**2)*smT2 + &
                      MV2*shat**3* &
                       (20*MV2**2 + 70*MV2*shat - 11*shat**2)*smT2**2 + &
                      shat**2*(-19*MV2**2 - 98*MV2*shat + 4*shat**2)* &
                       smT2**3 + shat*(-7*MV2 + 31*shat)*smT2**4 + &
                      5*smT2**5)*SP(4)) - &
                E**(6*y)*smT2* &
                 ((4*MV2**3*shat**6*(3*MV2 + shat) + &
                      2*MV2**2*shat**4* &
                       (2*MV2**2 - 11*MV2*shat + 7*shat**2)*smT2 + &
                      2*MV2*shat**3* &
                       (11*MV2**2 + 61*MV2*shat - 8*shat**2)*smT2**2 + &
                      shat**2*(-22*MV2**2 - 144*MV2*shat + shat**2)* &
                       smT2**3 + shat*(-6*MV2 + 41*shat)*smT2**4 + &
                      5*smT2**5)*SP(3) + &
                   (4*MV2**3*shat**5* &
                       (8*MV2**2 + 19*MV2*shat - 2*shat**2) - &
                      2*MV2**2*shat**4* &
                       (46*MV2**2 + 141*MV2*shat + 5*shat**2)*smT2 + &
                      2*MV2*shat**3* &
                       (31*MV2**2 + 152*MV2*shat + 4*shat**2)*smT2**2 + &
                      4*shat**3*(-41*MV2 + shat)*smT2**3 + &
                      3*shat*(-6*MV2 + 17*shat)*smT2**4 + 10*smT2**5)* &
                    SP(4)) + &
                E**(13*y)*smT**5* &
                 ((2*MV2*shat**4*(-19*MV2**2 - 14*MV2*shat + shat**2) + &
                      shat**2* &
                       (MV2**3 + 74*MV2**2*shat + 38*MV2*shat**2 - &
                         shat**3)*smT2 - &
                      shat*(MV2**2 + 50*MV2*shat + 7*shat**2)*smT2**2 + &
                      (MV2 + 17*shat)*smT2**3)*SP(3) - &
                   shat*(-2*smT2**2*(shat**2 + smT2) + &
                      MV2*shat*smT2*(4*shat**2 + 3*smT2) + &
                      MV2**2*(shat**4 + 2*shat**2*smT2))*SP(4)) + &
                E**(12*y)*smT2**2* &
                 ((2*MV2**2*shat**4* &
                       (16*MV2**2 + 50*MV2*shat + 19*shat**2) - &
                      2*MV2*shat**3*(3*MV2 + shat)*(9*MV2 + 23*shat)* &
                       smT2 + &
                      shat**2*(30*MV2**2 + 94*MV2*shat + 11*shat**2)* &
                       smT2**2 - 29*shat**2*smT2**3 - 5*smT2**4)*SP(3) &
                    + (MV2*shat**3*smT2*(2*shat**2 + 23*smT2) - &
                      smT2**2*(shat**4 + 9*shat**2*smT2 + smT2**2) + &
                      MV2**2* &
                       (2*shat**6 - 7*shat**4*smT2 + 4*shat**2*smT2**2)) &
                     *SP(4)) + &
                E**(9*y)*smT* &
                 (MV2*smT2**2* &
                    ((27*shat**6 - 146*shat**4*smT2 - &
                         169*shat**2*smT2**2 + 6*smT2**3)*SP(3) + &
                      (2*shat**6 - 97*shat**4*smT2 - &
                         90*shat**2*smT2**2 + 4*smT2**3)*SP(4)) + &
                   16*MV2**5* &
                    (shat**4*smT2*SP(3) + shat**6*(8*SP(3) + SP(4))) + &
                   2*shat*smT2**3* &
                    (-3*shat**4*SP(3) + 5*smT2**2*(5*SP(3) + 3*SP(4)) + &
                      shat**2*smT2*(15*SP(3) + 14*SP(4))) + &
                   2*MV2**4*shat**3* &
                    (5*shat**2*smT2*(-23*SP(3) + SP(4)) + &
                      smT2**2*(-23*SP(3) + SP(4)) + &
                      shat**4*(97*SP(3) + 29*SP(4))) + &
                   MV2**3*shat**2* &
                    (2*shat**2*smT2**2*(SP(3) - 12*SP(4)) + &
                      2*shat**6*(9*SP(3) + 5*SP(4)) + &
                      smT2**3*(47*SP(3) + 10*SP(4)) - &
                      shat**4*smT2*(469*SP(3) + 88*SP(4))) + &
                   MV2**2*shat*smT2* &
                    (73*shat**2*smT2**2*(3*SP(3) + SP(4)) - &
                      shat**6*(39*SP(3) + 10*SP(4)) - &
                      smT2**3*(23*SP(3) + 14*SP(4)) + &
                      shat**4*smT2*(391*SP(3) + 114*SP(4)))) - &
                E**(8*y)*(32*MV2**5*shat**5* &
                    (smT2*(SP(3) + SP(4)) + shat**2*(2*SP(3) + SP(4))) &
                    + 8*MV2**4*shat**4*(shat - smT)*(shat + smT)* &
                    (shat**2*(7*SP(3) + 3*SP(4)) + &
                      2*smT2*(5*SP(3) + 4*SP(4))) + &
                   smT2**4*(6*shat**4*(-SP(3) + SP(4)) + &
                      10*smT2**2*(SP(3) + SP(4)) + &
                      shat**2*smT2*(73*SP(3) + 54*SP(4))) - &
                   2*MV2**3*shat**3*smT2* &
                    (2*shat**4*(23*SP(3) + 7*SP(4)) - &
                      smT2**2*(43*SP(3) + 27*SP(4)) + &
                      shat**2*smT2*(175*SP(3) + 101*SP(4))) + &
                   2*MV2**2*shat**2*smT2**2* &
                    (-(smT2**2*(13*SP(3) + 3*SP(4))) + &
                      shat**4*(17*SP(3) + 7*SP(4)) + &
                      shat**2*smT2*(254*SP(3) + 139*SP(4))) + &
                   MV2*shat*smT2**3* &
                    (2*shat**4*(7*SP(3) - 5*SP(4)) - &
                      4*smT2**2*(4*SP(3) + 5*SP(4)) - &
                      shat**2*smT2*(304*SP(3) + 187*SP(4)))) + &
                E**(7*y)*smT* &
                 (shat*smT2**3* &
                    ((-4*shat**4 + 35*shat**2*smT2 + 35*smT2**2)* &
                       SP(3) + 8*smT2*(4*shat**2 + 5*smT2)*SP(4)) + &
                   MV2**2*shat*smT2* &
                    (2*(-9*shat**6 + 141*shat**4*smT2 + &
                         34*shat**2*smT2**2 - 7*smT2**3)*SP(3) + &
                      (-7*shat**6 + 197*shat**4*smT2 + &
                         167*shat**2*smT2**2 - 23*smT2**3)*SP(4)) + &
                   MV2*smT2**2* &
                    (2*(7*shat**6 - 84*shat**4*smT2 - &
                         53*shat**2*smT2**2 + 2*smT2**3)*SP(3) + &
                      (5*shat**6 - 98*shat**4*smT2 - &
                         123*shat**2*smT2**2 + 6*smT2**3)*SP(4)) + &
                   16*MV2**5* &
                    (shat**4*smT2*SP(4) + shat**6*(SP(3) + 4*SP(4))) + &
                   2*MV2**4*shat**3* &
                    (smT2**2*(SP(3) - 23*SP(4)) + &
                      shat**4*(45*SP(3) + 41*SP(4)) - &
                      shat**2*smT2*(11*SP(3) + 51*SP(4))) + &
                   MV2**3*shat**2* &
                    (2*shat**2*smT2**2*(12*SP(3) - 23*SP(4)) + &
                      2*shat**6*(5*SP(3) + SP(4)) + &
                      smT2**3*(10*SP(3) + 47*SP(4)) - &
                      shat**4*smT2*(224*SP(3) + 213*SP(4))))))/ &
            (E**(6*y)*shat**2*(-(E**y*shat) + smT)**2* &
              (shat - E**y*smT)**2*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(5)*SP(5))/ &
            (shat*(-shat + E**y*smT)*(MV2*shat - smT2)) + &
           ((-3*MV2*shat*smT2**3 + smT2**4 + &
                E**(12*y)*smT2**3*(MV2*shat + smT2) + &
                E**y*smT**5*(3*MV2*shat*(MV2 + 2*shat) + &
                   (MV2 - 2*shat)*smT2) - &
                E**(11*y)*smT**5* &
                 (MV2*shat*(MV2 + 3*shat) + (3*MV2 + shat)*smT2) + &
                E**(10*y)*smT2**2* &
                 (MV2*shat**2*(-3*MV2 + 2*shat) + &
                   MV2*(2*MV2 + 3*shat)*smT2 + 4*smT2**2) + &
                E**(3*y)*smT**3* &
                 (-(MV2**2*shat**2*(22*MV2 + 39*shat)) + &
                   2*MV2*shat*(7*MV2 + 22*shat)*smT2 + &
                   (4*MV2 - 9*shat)*smT2**2) + &
                E**(9*y)*smT**3* &
                 (3*MV2**2*shat**2*(2*MV2 + 5*shat) + &
                   MV2*(12*MV2 - 5*shat)*shat*smT2 - &
                   2*(7*MV2 + 3*shat)*smT2**2) + &
                E**(8*y)*smT2* &
                 (-2*MV2**2*shat**3*(13*MV2 + 6*shat) + &
                   MV2*shat*(-14*MV2**2 - 4*MV2*shat + 7*shat**2)* &
                    smT2 + (10*MV2**2 + 8*MV2*shat + shat**2)*smT2**2 + &
                   8*smT2**3) + &
                E**(4*y)*smT2* &
                 (2*MV2**2*shat**3*(13*MV2 + 10*shat) + &
                   MV2*shat*(14*MV2**2 + 28*MV2*shat - 19*shat**2)* &
                    smT2 + (-10*MV2**2 - 50*MV2*shat + 3*shat**2)* &
                    smT2**2 + 10*smT2**3) + &
                E**(6*y)*(-(MV2*shat) + smT2)* &
                 (8*MV2**2*shat**3*(4*MV2 + shat) - &
                   8*MV2*shat**2*(3*MV2 + shat)*smT2 + &
                   3*shat*(-4*MV2 + shat)*smT2**2 + 11*smT2**3) + &
                E**(7*y)*smT* &
                 (4*MV2**3*shat**3*(4*MV2 + 9*shat) - &
                   MV2**2*shat**2*(14*MV2 + 45*shat)*smT2 + &
                   MV2*shat*(23*MV2 + 31*shat)*smT2**2 - &
                   (17*MV2 + 14*shat)*smT2**3) + &
                E**(5*y)*smT* &
                 (4*MV2**3*shat**3*(4*MV2 + 3*shat) - &
                   3*MV2**2*shat**2*(22*MV2 + 25*shat)*smT2 + &
                   MV2*shat*(45*MV2 + 71*shat)*smT2**2 - &
                   (3*MV2 + 16*shat)*smT2**3) + &
                E**(2*y)*smT2**2* &
                 (MV2**2*(15*shat**2 - 2*smT2) + &
                   smT2*(shat**2 + 5*smT2) - &
                   3*MV2*(shat**3 + 8*shat*smT2)))*LC(4)*SP(6))/ &
            (E**(5*y)*(E**y*shat - smT)*(-shat + E**y*smT)* &
              (-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(E**y*(MV2 + shat)*smT**5 - &
                E**(9*y)*(MV2 + shat)*smT**5 + &
                E**(7*y)*(MV2 + shat)*smT**3*(7*MV2*shat - 5*smT2) - &
                2*E**(6*y)*(5*MV2*shat - 4*smT2)*smT2**2 - smT2**3 + &
                E**(10*y)*smT2**3 - &
                E**(3*y)*(MV2 + shat)*smT*(-3*MV2*shat + smT2)* &
                 (-2*MV2*shat + smT2) - &
                6*E**(5*y)*(MV2 + shat)*smT*(-(MV2*shat) + smT2)**2 + &
                E**(8*y)*smT2**2*(-6*MV2*shat + 5*smT2) - &
                E**(2*y)*smT2* &
                 (-4*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) + &
                E**(4*y)*(8*MV2**3*shat**3 - 4*MV2**2*shat**2*smT2 - &
                   6*MV2*shat*smT2**2 + 4*smT2**3))*LC(10)*SP(6))/ &
            (E**(4*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(E**(11*y)*smT**7 + shat*smT2**3 - &
                E**(10*y)*shat*smT2**2*(MV2**2 + 2*smT2) - &
                E**y*smT**5*(shat**2 + 2*smT2) + &
                E**(9*y)*shat*smT**3* &
                 (-4*MV2**2*shat + 4*MV2*smT2 + shat*smT2) + &
                E**(2*y)*smT2**2* &
                 (-(MV2*shat*(MV2 + 6*shat)) + 2*(MV2 + 3*shat)*smT2) + &
                E**(8*y)*MV2*smT2* &
                 (2*MV2*shat**2*(7*MV2 + 5*shat) - &
                   shat*(17*MV2 + 7*shat)*smT2 + 5*smT2**2) - &
                E**(3*y)*smT**3* &
                 (MV2*(MV2 - 7*shat)*shat**2 + &
                   shat*(-13*MV2 + 5*shat)*smT2 + 9*smT2**2) + &
                E**(4*y)*smT2* &
                 (8*MV2**3*shat**2 - MV2*shat*(19*MV2 + 17*shat)*smT2 + &
                   (9*MV2 + 13*shat)*smT2**2) - &
                2*E**(5*y)*smT* &
                 (MV2**2*shat**3*(-5*MV2 + 3*shat) + &
                   2*MV2*(5*MV2 - 3*shat)*shat**2*smT2 + &
                   shat*(-13*MV2 + 3*shat)*smT2**2 + 7*smT2**3) - &
                E**(7*y)*smT* &
                 (2*MV2**2*shat**3*(7*MV2 + 3*shat) - &
                   MV2*shat**2*(MV2 + 5*shat)*smT2 + &
                   shat*(-17*MV2 + shat)*smT2**2 + 8*smT2**3) + &
                2*E**(6*y)*(4*MV2**3*shat**3*(-2*MV2 + shat) + &
                   MV2**2*shat**2*(19*MV2 + shat)*smT2 - &
                   MV2*shat*(17*MV2 + 9*shat)*smT2**2 + &
                   (6*MV2 + 5*shat)*smT2**3))*LC(14)*SP(6))/ &
            (E**(6*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           ((-2*E**(11*y)*shat*smT**7 + &
                E**y*(MV2 + 2*shat)*smT**5*(3*MV2*shat - smT2) + &
                smT2**3*(-3*MV2*shat + smT2) + &
                E**(2*y)*smT2**2*(-3*MV2*shat + smT2)* &
                 (-5*MV2*shat + shat**2 + 5*smT2) + &
                E**(10*y)*smT2**2* &
                 (2*MV2*shat**3 + 3*MV2*shat*smT2 + smT2**2) - &
                E**(9*y)*smT**3* &
                 (5*MV2**2*shat**3 + MV2*(MV2 - 14*shat)*shat*smT2 + &
                   (MV2 + 9*shat)*smT2**2) - &
                E**(3*y)*smT**3* &
                 (MV2**2*shat**2*(22*MV2 + 39*shat) - &
                   2*MV2*shat*(11*MV2 + 21*shat)*smT2 + &
                   (4*MV2 + 9*shat)*smT2**2) + &
                E**(8*y)*smT2* &
                 (2*MV2**2*(MV2 - 6*shat)*shat**3 + &
                   MV2*shat**2*(-17*MV2 + 7*shat)*smT2 + &
                   shat*(4*MV2 + shat)*smT2**2 + 5*smT2**3) + &
                E**(4*y)*smT2* &
                 (2*MV2**2*shat**3*(13*MV2 + 10*shat) + &
                   MV2*(5*MV2 - 19*shat)*shat**2*smT2 + &
                   3*shat*(-11*MV2 + shat)*smT2**2 + 10*smT2**3) + &
                E**(7*y)*smT* &
                 (44*MV2**3*shat**4 + &
                   MV2**2*(6*MV2 - 61*shat)*shat**2*smT2 + &
                   2*MV2*shat*(MV2 + 21*shat)*smT2**2 - &
                   (4*MV2 + 17*shat)*smT2**3) + &
                E**(5*y)*smT* &
                 (4*MV2**3*shat**3*(4*MV2 + 3*shat) - &
                   MV2**2*shat**2*(32*MV2 + 63*shat)*smT2 + &
                   2*MV2*shat*(11*MV2 + 32*shat)*smT2**2 - &
                   (6*MV2 + 17*shat)*smT2**3) + &
                E**(6*y)*(-8*MV2**3*shat**4*(4*MV2 + shat) + &
                   4*MV2**2*shat**3*(9*MV2 + 4*shat)*smT2 - &
                   MV2*shat**2*(3*MV2 + 11*shat)*smT2**2 + &
                   3*shat*(-5*MV2 + shat)*smT2**3 + 10*smT2**4))*LC(2)* &
              SP(7))/ &
            (E**(5*y)*shat*(E**y*shat - smT)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-(E**y*(MV2 + shat)*smT**5) + &
                E**(9*y)*(MV2 + shat)*smT**5 + &
                E**(8*y)*(6*MV2*shat - 5*smT2)*smT2**2 + &
                2*E**(6*y)*(5*MV2*shat - 4*smT2)*smT2**2 + smT2**3 - &
                E**(10*y)*smT2**3 + &
                E**(3*y)*(MV2 + shat)*smT*(-3*MV2*shat + smT2)* &
                 (-2*MV2*shat + smT2) + &
                6*E**(5*y)*(MV2 + shat)*smT*(-(MV2*shat) + smT2)**2 + &
                E**(7*y)*(MV2 + shat)*smT**3*(-7*MV2*shat + 5*smT2) + &
                E**(2*y)*smT2* &
                 (-4*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) + &
                E**(4*y)*(-8*MV2**3*shat**3 + 4*MV2**2*shat**2*smT2 + &
                   6*MV2*shat*smT2**2 - 4*smT2**3))*LC(8)*SP(7))/ &
            (E**(4*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-(E**(11*y)*smT**7) - shat*smT2**3 + &
                E**(10*y)*shat*smT2**2*(MV2**2 + 2*smT2) + &
                E**y*smT**5*(shat**2 + 2*smT2) + &
                E**(9*y)*shat*smT**3* &
                 (4*MV2**2*shat - (4*MV2 + shat)*smT2) + &
                E**(2*y)*smT2**2* &
                 (MV2*shat*(MV2 + 6*shat) - 2*(MV2 + 3*shat)*smT2) - &
                E**(8*y)*MV2*smT2* &
                 (2*MV2*shat**2*(7*MV2 + 5*shat) - &
                   shat*(17*MV2 + 7*shat)*smT2 + 5*smT2**2) + &
                E**(3*y)*smT**3* &
                 (MV2*(MV2 - 7*shat)*shat**2 + &
                   shat*(-13*MV2 + 5*shat)*smT2 + 9*smT2**2) - &
                E**(4*y)*smT2* &
                 (8*MV2**3*shat**2 - MV2*shat*(19*MV2 + 17*shat)*smT2 + &
                   (9*MV2 + 13*shat)*smT2**2) + &
                2*E**(5*y)*smT* &
                 (MV2**2*shat**3*(-5*MV2 + 3*shat) + &
                   2*MV2*(5*MV2 - 3*shat)*shat**2*smT2 + &
                   shat*(-13*MV2 + 3*shat)*smT2**2 + 7*smT2**3) + &
                E**(7*y)*smT* &
                 (2*MV2**2*shat**3*(7*MV2 + 3*shat) - &
                   MV2*shat**2*(MV2 + 5*shat)*smT2 + &
                   shat*(-17*MV2 + shat)*smT2**2 + 8*smT2**3) - &
                2*E**(6*y)*(4*MV2**3*shat**3*(-2*MV2 + shat) + &
                   MV2**2*shat**2*(19*MV2 + shat)*smT2 - &
                   MV2*shat*(17*MV2 + 9*shat)*smT2**2 + &
                   (6*MV2 + 5*shat)*smT2**3))*LC(12)*SP(7))/ &
            (E**(6*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
        kap*SI(0)*((12*MV2*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(1))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (6*MV2*(smT2 + E**(2*y)*(-4*MV2*shat + 3*smT2))*LC(11))/ &
            (E**(2*y)*(MV2*shat - smT2)* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (20*MV2*shat*smT*LC(15)*Sinh(y))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) + &
           (8*(-smT**3 - E**(3*y)*MV2*(4*MV2*shat - 3*smT2) + &
                E**(2*y)*smT*(3*MV2*shat - 2*smT2) + &
                E**(4*y)*smT*(MV2*shat - smT2) + E**y*MV2*smT2)*LC(13)* &
              SP(4))/ &
            (E**(3*y)*(-shat + E**y*smT)*(MV2*shat - smT2)* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (8*LC(3)*(E**(6*y)*smT**5*SP(3) + smT**5*SP(4) - &
                E**(3*y)*(-2*MV2*shat + smT2)* &
                 (2*shat*smT2 + MV2*(-shat**2 + smT2))*(SP(3) + SP(4)) &
                 - E**y*smT2* &
                 (shat*(-(MV2*shat) + smT2)*SP(3) + &
                   (MV2 + shat)*smT2*SP(4)) - &
                E**(5*y)*smT2* &
                 ((MV2 + shat)*smT2*SP(3) + &
                   shat*(-(MV2*shat) + smT2)*SP(4)) + &
                E**(2*y)*smT* &
                 ((-2*MV2**2*shat**2 + smT2**2)*SP(3) + &
                   2*(-(MV2*shat) + smT2)**2*SP(4)) + &
                E**(4*y)*smT* &
                 (2*(-(MV2*shat) + smT2)**2*SP(3) + &
                   (-2*MV2**2*shat**2 + smT2**2)*SP(4))))/ &
            (E**(2*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)* &
              (MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
            + (4*LC(9)*(-2*E**(5*y)*smT**3*SP(3) + &
                2*E**(4*y)*MV2*smT2*SP(3) - MV2*smT2*SP(4) - &
                E**(3*y)*smT* &
                 (4*smT2*SP(3) + MV2*shat*(-6*SP(3) + SP(4))) + &
                E**y*smT*(-2*smT2*SP(3) + MV2*shat*(2*SP(3) + SP(4))) + &
                E**(2*y)*MV2* &
                 (-8*MV2*shat*SP(3) + smT2*(6*SP(3) + SP(4)))))/ &
            (E**y*(E**y*shat - smT)*(MV2*shat - smT2)* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (8*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(10)*SP(6))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (8*(2*MV2*shat - (2*smT2*Cosh(y))/E**y)*LC(14)*SP(6))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) + &
           (8*MV2*smT*LC(4)*Sinh(y)*SP(6))/ &
            ((-(MV2*shat) + smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
             + (8*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(8)*SP(7))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (8*(-2*MV2*shat + smT2 + smT2/E**(2*y))*LC(12)*SP(7))/ &
            ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
           (8*MV2*smT*LC(2)*Sinh(y)*SP(7))/ &
            ((-(MV2*shat) + smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
           ) + kap*SI(5)*(((4*smT2**2*(-(MV2*shat) + smT2)**2 - &
                3*E**(9*y)*shat*smT**5*(MV2*shat + 3*smT2) + &
                E**(7*y)*shat*smT**3* &
                 (18*MV2**2*shat**2 + 47*MV2*shat*smT2 - 41*smT2**2) + &
                E**(8*y)*smT2**2* &
                 (3*MV2**2*shat**2 + 17*MV2*shat*smT2 + 4*smT2**2) - &
                2*E**y*shat*smT**3* &
                 (8*MV2**2*shat**2 - 21*MV2*shat*smT2 + 7*smT2**2) + &
                3*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
                 (16*MV2**2*shat**2 - 32*MV2*shat*smT2 + 23*smT2**2) + &
                E**(2*y)*smT2* &
                 (24*MV2**3*shat**3 - 31*MV2**2*shat**2*smT2 - &
                   33*MV2*shat*smT2**2 + 16*smT2**3) + &
                E**(6*y)*smT2* &
                 (-44*MV2**3*shat**3 - 61*MV2**2*shat**2*smT2 + &
                   17*MV2*shat*smT2**2 + 16*smT2**3) - &
                E**(3*y)*shat*smT* &
                 (8*MV2**3*shat**3 + 122*MV2**2*shat**2*smT2 - &
                   157*MV2*shat*smT2**2 + 51*smT2**3) + &
                E**(4*y)*(-16*MV2**4*shat**4 + &
                   148*MV2**3*shat**3*smT2 - &
                   59*MV2**2*shat**2*smT2**2 - 25*MV2*shat*smT2**3 + &
                   24*smT2**4))*LC(1))/ &
            (2d0*E**(4*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           ((-4*E**(9*y)*smT2**2*(-(MV2*shat) + smT2)**2 + &
                shat*smT**5*(MV2*shat + 5*smT2) - &
                2*E**(2*y)*shat*smT**3* &
                 (5*MV2**2*shat**2 + 13*MV2*shat*smT2 - 12*smT2**2) + &
                E**y*smT2**2* &
                 (MV2**2*shat**2 - 9*MV2*shat*smT2 - 4*smT2**2) + &
                E**(8*y)*shat*smT**3* &
                 (16*MV2**2*shat**2 - 31*MV2*shat*smT2 + 9*smT2**2) + &
                6*E**(4*y)*shat*smT* &
                 (4*MV2**3*shat**3 - 11*MV2*shat*smT2**2 + 7*smT2**3) + &
                2*E**(6*y)*shat*smT* &
                 (4*MV2**3*shat**3 + 21*MV2**2*shat**2*smT2 - &
                   35*MV2*shat*smT2**2 + 16*smT2**3) - &
                E**(7*y)*smT2* &
                 (24*MV2**3*shat**3 - 27*MV2**2*shat**2*smT2 - &
                   25*MV2*shat*smT2**2 + 16*smT2**3) - &
                E**(3*y)*smT2* &
                 (8*MV2**3*shat**3 - 69*MV2**2*shat**2*smT2 + &
                   9*MV2*shat*smT2**2 + 16*smT2**3) + &
                E**(5*y)*(16*MV2**4*shat**4 - 96*MV2**3*shat**3*smT2 + &
                   51*MV2**2*shat**2*smT2**2 + 17*MV2*shat*smT2**3 - &
                   24*smT2**4))*LC(11))/ &
            (2d0*E**(5*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (shat*LC(15)*(-((MV2*shat - smT2)* &
                   (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)) + &
                smT*(4*smT*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
                    Cosh(2*y) + &
                   (-(MV2*shat*smT**3) + smT**5)*Cosh(4*y) + &
                   4*MV2*(18*MV2**2*shat**2 + 25*MV2*shat*smT2 - &
                      13*smT2**2)*Sinh(y) + &
                   10*smT*(-8*MV2**2*shat**2 + MV2*shat*smT2 + &
                      3*smT2**2)*Sinh(2*y) + &
                   4*MV2*(3*MV2*shat - 13*smT2)*smT2*Sinh(3*y) + &
                   5*smT**3*(MV2*shat + 3*smT2)*Sinh(4*y))))/ &
            (2d0*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*LC(3)*(-(E**(10*y)*smT**5*(2*MV2*shat - smT2)*SP(3)) + &
                2*E**(9*y)*MV2**2*shat*smT2**2*SP(3) - &
                MV2*shat*smT**5*SP(4) + &
                2*E**y*MV2**2*shat*smT2**2*SP(4) + &
                2*E**(5*y)*MV2**2*shat* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
                 (SP(3) + SP(4)) - &
                2*E**(3*y)*MV2**2*shat*smT2* &
                 (smT2*(SP(3) - 4*SP(4)) + 7*MV2*shat*SP(4)) - &
                E**(4*y)*smT* &
                 (2*(5*MV2**3*shat**3 - 9*MV2**2*shat**2*smT2 + &
                      7*MV2*shat*smT2**2 - 2*smT2**3)*SP(3) + &
                   MV2*shat*(2*MV2*shat - 5*smT2)*(MV2*shat - smT2)* &
                    SP(4)) - &
                E**(6*y)*smT* &
                 (3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
                    (2*MV2*shat - smT2)*SP(3) + &
                   2*MV2*shat* &
                    (5*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)*SP(4) &
                   ) - 2*E**(7*y)*MV2**2*shat*smT2* &
                 (7*MV2*shat*SP(3) + smT2*(-4*SP(3) + SP(4))) + &
                E**(8*y)*smT**3* &
                 (-16*MV2*shat*smT2*SP(3) + 4*smT2**2*SP(3) + &
                   MV2**2*shat**2*(14*SP(3) + SP(4))) + &
                E**(2*y)*smT**3* &
                 (smT2**2*SP(3) - MV2*shat*smT2*(SP(3) + 4*SP(4)) + &
                   MV2**2*shat**2*(SP(3) + 6*SP(4)))))/ &
            (E**(5*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (LC(9)*(4*E**(10*y)*shat*smT2**3*SP(3) - &
                4*E**(9*y)*smT**5*(MV2*shat + smT2)*SP(3) + &
                shat*smT2**2*(-3*MV2*shat + smT2)*SP(4) - &
                E**(8*y)*shat*smT2**2* &
                 (4*(7*MV2*shat - 5*smT2)*SP(3) - &
                   (MV2*shat + smT2)*SP(4)) + &
                2*E**(6*y)*shat*smT2* &
                 (18*(-(MV2*shat) + smT2)**2*SP(3) - &
                   (3*MV2*shat - 2*smT2)*(MV2*shat + smT2)*SP(4)) - &
                E**y*smT**3*(4*(-(MV2*shat) + smT2)**2*SP(3) - &
                   MV2*shat*(MV2*shat + 3*smT2)*SP(4)) + &
                E**(7*y)*smT**3* &
                 (4*(7*MV2*shat - 4*smT2)*(MV2*shat + smT2)*SP(3) - &
                   MV2*shat*(MV2*shat + 3*smT2)*SP(4)) - &
                2*E**(4*y)*shat* &
                 (-2*smT2*(12*MV2**2*shat**2 - 17*MV2*shat*smT2 + &
                      7*smT2**2)*SP(3) + &
                   (MV2*shat - smT2)* &
                    (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(4)) - &
                E**(5*y)*smT* &
                 (4*(MV2*shat + 2*smT2)* &
                    (6*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(3) + MV2*shat* &
                    (-16*MV2**2*shat**2 + MV2*shat*smT2 + 3*smT2**2)* &
                    SP(4)) + &
                E**(3*y)*(-16*smT**7*SP(3) + &
                   MV2**2*shat**2*smT**3*(-4*SP(3) + SP(4)) + &
                   3*MV2*shat*smT**5*(12*SP(3) + SP(4)) - &
                   8*MV2**3*shat**3*smT*(SP(3) + 2*SP(4))) + &
                2*E**(2*y)*shat*smT2* &
                 (2*smT2**2*(2*SP(3) + SP(4)) + &
                   MV2**2*shat**2*(6*SP(3) + 11*SP(4)) - &
                   MV2*shat*smT2*(12*SP(3) + 11*SP(4)))))/ &
            (E**(4*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (LC(13)*(E**(10*y)*shat*(MV2*shat - smT2)*smT2**2*SP(3) + &
                4*E**(2*y)*shat*(7*MV2*shat - 5*smT2)*smT2**2*SP(4) - &
                4*shat*smT2**3*SP(4) + &
                4*E**y*smT**5*(MV2*shat + smT2)*SP(4) - &
                E**(9*y)*smT**3*(MV2*shat - smT2)* &
                 (MV2*shat*(SP(3) - 4*SP(4)) + 4*smT2*SP(4)) + &
                E**(3*y)*smT**3* &
                 (MV2*shat*(MV2*shat - smT2)*SP(3) - &
                   4*(7*MV2*shat - 4*smT2)*(MV2*shat + smT2)*SP(4)) + &
                E**(7*y)*smT* &
                 (MV2*shat*(MV2*shat - smT2)*(4*MV2*shat - smT2)* &
                    SP(3) + 4*(2*MV2*shat - smT2)* &
                    (MV2**2*shat**2 + MV2*shat*smT2 - 4*smT2**2)*SP(4)) &
                 - E**(5*y)*smT* &
                 (MV2*shat*(MV2*shat - smT2)*(4*MV2*shat - smT2)* &
                    SP(3) - 4*(MV2*shat + 2*smT2)* &
                    (6*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(4)) + &
                E**(6*y)*shat* &
                 ((MV2*shat - smT2)* &
                    (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
                    SP(3) - 4*smT2* &
                    (12*MV2**2*shat**2 - 17*MV2*shat*smT2 + 7*smT2**2)* &
                    SP(4)) - &
                E**(8*y)*shat*smT2* &
                 (6*MV2**2*shat**2*(SP(3) + 2*SP(4)) - &
                   3*MV2*shat*smT2*(3*SP(3) + 8*SP(4)) + &
                   smT2**2*(3*SP(3) + 8*SP(4))) - &
                E**(4*y)*shat*(MV2*shat - smT2)*smT2* &
                 (2*MV2*shat*(SP(3) + 18*SP(4)) - &
                   smT2*(SP(3) + 36*SP(4)))))/ &
            (E**(6*y)*shat*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 - &
                E**(6*y)*shat*smT**3*(7*MV2*shat - 5*smT2) + &
                6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 - &
                E**(7*y)*smT2**2*(MV2*shat + smT2) + &
                E**(5*y)*smT2* &
                 (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
                E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
                   smT2**2) + &
                E**(2*y)*shat*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
                E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(10)*SP(6))/ &
            (E**(3*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 + &
                E**(2*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
                6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
                E**y*smT2**2*(MV2*shat + smT2) + &
                E**(7*y)*smT2* &
                 (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
                E**(6*y)*shat*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(3*y)*smT2* &
                 (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
                E**(5*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(14)*SP(6))/ &
            (E**(5*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (2*LC(4)*(-(MV2*shat) + smT2 + &
                (12*MV2*shat*smT*(MV2*shat - smT2)*(MV2 - smT*Cosh(y))* &
                   Sinh(y))/(-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2 - &
                (2*smT*(2*MV2*(MV2*shat + smT2) + &
                     (-5*MV2*shat*smT + smT**3)*Cosh(y))*Sinh(y))/ &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y)))*SP(6))/ &
            (-(MV2*shat) + smT2)**2 + &
           (4*(shat*smT**5 - E**(8*y)*shat*smT**5 + &
                E**(6*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
                6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
                E**(7*y)*smT2**2*(MV2*shat + smT2) + &
                E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
                   smT2**2) - &
                E**(2*y)*shat*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(5*y)*smT2* &
                 (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
                E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(8)*SP(7))/ &
            (E**(3*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
           (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 + &
                E**(2*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
                6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
                E**y*smT2**2*(MV2*shat + smT2) + &
                E**(7*y)*smT2* &
                 (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
                E**(6*y)*shat*smT* &
                 (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
                E**(3*y)*smT2* &
                 (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
                E**(5*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                   7*MV2*shat*smT2**2 + 3*smT2**3))*LC(12)*SP(7))/ &
            (E**(5*y)*(-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
           (LC(2)*((MV2*shat - smT2)* &
                 (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2) + &
                smT*(-4*smT*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
                    Cosh(2*y) + smT**3*(MV2*shat - smT2)*Cosh(4*y) + &
                   4*MV2*(-2*MV2*shat + smT2)*(5*MV2*shat + smT2)* &
                    Sinh(y) + &
                   2*smT*(16*MV2**2*shat**2 - 13*MV2*shat*smT2 + &
                      smT2**2)*Sinh(2*y) + &
                   4*MV2*smT2*(MV2*shat + smT2)*Sinh(3*y) + &
                   smT**3*(-5*MV2*shat + smT2)*Sinh(4*y)))*SP(7))/ &
            ((-(MV2*shat) + smT2)**2* &
              (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
        SI(13)*(kap*(((-12*MV2*shat**2*smT2**4 + &
                   E**(14*y)*shat*smT2**4*(5*MV2*shat + 7*smT2) - &
                   2*E**(13*y)*smT**7* &
                    (2*MV2**2*shat**2 + 9*MV2*shat*smT2 + smT2**2) + &
                   E**(12*y)*smT2**3* &
                    (MV2**2*(-8*MT2 + 3*MV2 - 51*shat)*shat**2 + &
                      MV2*(-8*MT2 + 17*MV2 - 47*shat)*shat*smT2 + &
                      2*(8*MT2 + 2*MV2 + 19*shat)*smT2**2) + &
                   E**(11*y)*smT**5* &
                    (42*MV2**3*shat**3 + 143*MV2**2*shat**2*smT2 - &
                      76*MV2*shat*smT2**2 - 13*smT2**3) + &
                   E**y*smT**5* &
                    (4*MV2**3*shat**3 - MV2**2*shat**2*smT2 + &
                      24*MV2*shat*smT2**2 - 3*smT2**3) + &
                   E**(2*y)*smT2**2* &
                    (-8*MV2**3*shat**3*(4*MT2 + MV2 + shat) + &
                      MV2**2*shat**2*(104*MT2 + 21*MV2 + 145*shat)* &
                       smT2 - &
                      MV2*shat*(88*MT2 + 41*(MV2 + 2*shat))*smT2**2 + &
                      (16*MT2 + 4*MV2 + 5*shat)*smT2**3) + &
                   E**(10*y)*smT2**2* &
                    (2*MV2**3*shat**3*(40*MT2 - 15*MV2 + 79*shat) + &
                      MV2**2*shat**2*(40*MT2 - 155*MV2 + 81*shat)* &
                       smT2 + &
                      MV2*(-200*MT2 + 45*MV2 - 194*shat)*shat*smT2**2 + &
                      (80*MT2 + 20*MV2 + 87*shat)*smT2**3) + &
                   2*E**(5*y)*smT**3* &
                    (292*MV2**4*shat**4 - 356*MV2**3*shat**3*smT2 + &
                      101*MV2**2*shat**2*smT2**2 + &
                      43*MV2*shat*smT2**3 - 20*smT2**4) + &
                   E**(3*y)*smT**3* &
                    (56*MV2**4*shat**4 - 254*MV2**3*shat**3*smT2 + &
                      35*MV2**2*shat**2*smT2**2 + 84*MV2*shat*smT2**3 - &
                      17*smT2**4) - &
                   2*E**(7*y)*smT*(-(MV2*shat) + smT2)* &
                    (-48*MV2**4*shat**4 + 108*MV2**3*shat**3*smT2 - &
                      162*MV2**2*shat**2*smT2**2 + &
                      37*MV2*shat*smT2**3 + 25*smT2**4) - &
                   E**(9*y)*smT**3* &
                    (152*MV2**4*shat**4 + 204*MV2**3*shat**3*smT2 - &
                      379*MV2**2*shat**2*smT2**2 + &
                      108*MV2*shat*smT2**3 + 35*smT2**4) + &
                   E**(4*y)*smT2* &
                    (-16*MV2**4*shat**4*(-10*MT2 + 3*(MV2 + shat)) + &
                      2*MV2**3*(-360*MT2 + 67*MV2 - 139*shat)*shat**3* &
                       smT2 + &
                      MV2**2*shat**2*(920*MT2 + 131*MV2 + 357*shat)* &
                       smT2**2 - &
                      MV2*shat*(440*MT2 + 117*MV2 + 193*shat)*smT2**3 + &
                      10*(8*MT2 + 2*MV2 + 3*shat)*smT2**4) + &
                   2*E**(8*y)*smT2* &
                    (4*MV2**4*(-20*MT2 + 17*MV2 - 9*shat)*shat**4 + &
                      MV2**3*(-120*MT2 + 141*MV2 - 101*shat)*shat**3* &
                       smT2 + &
                      MV2**2*shat**2*(440*MT2 - 111*MV2 + 139*shat)* &
                       smT2**2 + &
                      2*MV2*(-160*MT2 + MV2 - 73*shat)*shat*smT2**3 + &
                      2*(40*MT2 + 10*MV2 + 27*shat)*smT2**4) + &
                   E**(6*y)*(-32*MV2**11 + &
                      32*MV2**5* &
                       (MV2**6 + (-4*MT2 + MV2)*shat**5 + shat**6) + &
                      8*MV2**4*shat**4*(128*MT2 - 47*MV2 + 7*shat)* &
                       smT2 + &
                      2*MV2**3*(-952*MT2 + 83*MV2 - 11*shat)*shat**3* &
                       smT2**2 - &
                      2*MV2**2*(-824*MT2 + MV2 - 155*shat)*shat**2* &
                       smT2**3 - &
                      MV2*shat*(100*(8*MT2 + MV2) + 273*shat)*smT2**4 + &
                      (40*(4*MT2 + MV2) + 77*shat)*smT2**5))*LC(1))/ &
               (4d0*E**(6*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((-6*shat*smT2**5 - &
                   E**(14*y)*shat*(MV2*shat - 3*smT2)*smT2**3* &
                    (2*MV2*shat + smT2) + &
                   2*E**y*smT**9*(5*MV2*shat + smT2) + &
                   E**(13*y)*smT**7* &
                    (5*MV2**2*shat**2 - 19*MV2*shat*smT2 + 2*smT2**2) - &
                   E**(2*y)*smT2**3* &
                    (-(MV2**2*shat**2*(-4*MT2 + MV2 + shat)) + &
                      MV2*(-20*MT2 + 9*MV2 - 59*shat)*shat*smT2 + &
                      2*(8*MT2 + 2*MV2 + 15*shat)*smT2**2) + &
                   E**(3*y)*smT**5* &
                    (-4*MV2**3*shat**3 - 93*MV2**2*shat**2*smT2 + &
                      37*MV2*shat*smT2**2 + 12*smT2**3) + &
                   E**(12*y)*smT2**2* &
                    (8*MV2**3*shat**3*(4*MT2 + MV2 + 3*shat) - &
                      MV2**2*shat**2*(92*MT2 + 25*MV2 + 73*shat)*smT2 + &
                      MV2*shat*(76*MT2 + 33*MV2 + 8*shat)*smT2**2 + &
                      (-4*(4*MT2 + MV2) + 11*shat)*smT2**3) - &
                   E**(4*y)*smT2**2* &
                    (2*MV2**3*shat**3*(-20*MT2 + 5*MV2 + 6*shat) + &
                      MV2**2*shat**2*(220*MT2 - 95*MV2 + 157*shat)* &
                       smT2 + &
                      MV2*(-260*MT2 + 5*MV2 - 164*shat)*shat*smT2**2 + &
                      (80*MT2 + 20*MV2 + 61*shat)*smT2**3) + &
                   E**(11*y)*smT**3* &
                    (-88*MV2**4*shat**4 + 210*MV2**3*shat**3*smT2 - &
                      37*MV2**2*shat**2*smT2**2 - 49*MV2*shat*smT2**3 + &
                      12*smT2**4) + &
                   2*E**(7*y)*smT*(-(MV2*shat) + smT2)* &
                    (40*MV2**4*shat**4 + 28*MV2**3*shat**3*smT2 - &
                      121*MV2**2*shat**2*smT2**2 + &
                      38*MV2*shat*smT2**3 + 20*smT2**4) + &
                   E**(5*y)*smT**3* &
                    (36*MV2**4*shat**4 + 222*MV2**3*shat**3*smT2 - &
                      287*MV2**2*shat**2*smT2**2 + &
                      59*MV2*shat*smT2**3 + 30*smT2**4) + &
                   E**(10*y)*smT2* &
                    (16*MV2**4*shat**4*(-10*MT2 + 3*MV2 + shat) + &
                      2*MV2**3*shat**3*(300*MT2 - 47*MV2 + 10*shat)* &
                       smT2 + &
                      MV2**2*shat**2*(-740*MT2 - 71*MV2 + 49*shat)* &
                       smT2**2 + &
                      MV2*(380*MT2 + 77*MV2 - 27*shat)*shat*smT2**3 - &
                      4*(20*MT2 + 5*MV2 - 2*shat)*smT2**4) - &
                   E**(6*y)*smT2* &
                    (-8*MV2**4*shat**4*(-16*MT2 + 4*MV2 + 5*shat) + &
                      2*MV2**3*shat**3*(-372*MT2 + 153*MV2 + 2*shat)* &
                       smT2 + &
                      2*MV2**2*shat**2*(572*MT2 - 83*MV2 + 4*shat)* &
                       smT2**2 - &
                      MV2*shat*(688*MT2 + 28*MV2 + 123*shat)*smT2**3 + &
                      (40*(4*MT2 + MV2) + 61*shat)*smT2**4) + &
                   2*E**(9*y)*smT* &
                    (32*MV2**5*shat**5 - 178*MV2**4*shat**4*smT2 + &
                      197*MV2**3*shat**3*smT2**2 - &
                      83*MV2**2*shat**2*smT2**3 - 13*MV2*shat*smT2**4 + &
                      15*smT2**5) - &
                   2*E**(8*y)* &
                    (16*MV2**5*shat**5*(-4*MT2 + MV2 + shat) - &
                      4*MV2**4*shat**4*(-92*MT2 + 26*MV2 + 25*shat)* &
                       smT2 + &
                      MV2**3*shat**3*(-700*MT2 + 71*MV2 + 222*shat)* &
                       smT2**2 + &
                      MV2**2*(692*MT2 - 29*MV2 - 111*shat)*shat**2* &
                       smT2**3 - &
                      2*MV2*(188*MT2 + 17*MV2 - 3*shat)*shat*smT2**4 + &
                      4*(20*MT2 + 5*MV2 + 3*shat)*smT2**5))*LC(11))/ &
               (4d0*E**(8*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
              (LC(15)*(2*(MV2*shat - smT2)* &
                    (8*MV2**4*shat**4 - 32*MV2**3*shat**3*smT2 + &
                      35*MV2**2*shat**2*smT2**2 - 21*MV2*shat*smT2**3 + &
                      5*smT2**4) + &
                   smT*(2*shat*(-(MV2*shat) + smT2)* &
                       (-32*MV2**3*shat**3 + 64*MV2**2*shat**2*smT2 - &
                         52*MV2*shat*smT2**2 + 15*smT2**3)*Cosh(y) - &
                      smT*(-(MV2*shat) + smT2)* &
                       (-16*MV2**3*shat**3 + 72*MV2**2*shat**2*smT2 - &
                         56*MV2*shat*smT2**2 + 15*smT2**3)*Cosh(2*y) + &
                      2*shat*smT2*(-(MV2*shat) + smT2)* &
                       (32*MV2**2*shat**2 - 34*MV2*shat*smT2 + &
                         11*smT2**2)*Cosh(3*y) + &
                      2*smT**3*(MV2*shat - smT2)* &
                       (MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
                       Cosh(4*y) + &
                      10*shat*smT2**2*(-2*MV2*shat + smT2)* &
                       (-(MV2*shat) + smT2)*Cosh(5*y) + &
                      smT**7*(MV2*shat - smT2)*Cosh(6*y) + &
                      2*shat*smT2**3*(-(MV2*shat) + smT2)*Cosh(7*y) + &
                      2*(8*MV2**4*(-28*MT2 + 13*MV2 - 2*shat)*shat**4 + &
                         2*MV2**3*(100*MT2 + 86*MV2 - 67*shat)*shat**3* &
                          smT2 + &
                         4*MV2**2*shat**2*(26*MT2 - 25*MV2 + 12*shat)* &
                          smT2**2 + &
                         MV2*(-112*MT2 + 8*MV2 - 85*shat)*shat* &
                          smT2**3 + (32*MT2 + 16*MV2 + 37*shat)*smT2**4) &
                        *Sinh(y) - &
                      smT*(328*MV2**4*shat**4 + &
                         196*MV2**3*shat**3*smT2 - &
                         546*MV2**2*shat**2*smT2**2 + &
                         187*MV2*shat*smT2**3 + 35*smT2**4)*Sinh(2*y) + &
                      4*smT2* &
                       (MV2**3*shat**3*(36*MT2 - 2*MV2 + 53*shat) + &
                         MV2**2*shat**2*(2*MT2 - 69*MV2 + 41*shat)* &
                          smT2 + &
                         MV2*(-62*MT2 + 9*(MV2 - 8*shat))*shat* &
                          smT2**2 + 3*(8*MT2 + 4*MV2 + 11*shat)*smT2**3) &
                        *Sinh(3*y) + &
                      2*smT**3* &
                       (23*MV2**3*shat**3 + 138*MV2**2*shat**2*smT2 - &
                         67*MV2*shat*smT2**2 - 14*smT2**3)*Sinh(4*y) + &
                      4*smT2**2* &
                       (MV2**2*(-2*MT2 + MV2 - 15*shat)*shat**2 + &
                         MV2*(-6*MT2 + 5*MV2 - 28*shat)*shat*smT2 + &
                         2*(4*MT2 + 2*MV2 + 9*shat)*smT2**2)*Sinh(5*y) &
                       - smT**5* &
                       (6*MV2**2*shat**2 + 27*MV2*shat*smT2 + &
                         7*smT2**2)*Sinh(6*y) + &
                      2*shat*smT2**3*(3*MV2*shat + 7*smT2)*Sinh(7*y))))/ &
               (4d0*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
              (2*(-2*E**(3*y)*MV2*shat**2*smT + E**y*shat*smT**3 + &
                   E**(5*y)*shat*smT**3 - MV2*shat*smT2 + &
                   E**(4*y)*smT2*(-2*MV2*shat + smT2) + &
                   E**(2*y)*(4*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
                      smT2**2))*LC(7)*SP(1))/ &
               (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
              (8*smT*(shat + smT*Cosh(y))*LC(6)*Sinh(y)*SP(2))/ &
               (shat*(MV2*shat - smT2)) + &
              (LC(13)*(-(E**(15*y)*smT2**3*(-(MV2*shat) + smT2)* &
                      (-(MV2*shat) + 2*smT2)*SP(3)) + &
                   4*shat*smT**9*SP(4) - &
                   2*E**y*smT2**4*(3*MV2*shat + smT2)*SP(4) - &
                   4*E**(2*y)*smT**7* &
                    (MV2*shat*(2*MT2 - MV2 + 11*shat) - &
                      (2*MT2 + MV2 + 6*shat)*smT2)*SP(4) + &
                   E**(14*y)*shat*smT**5* &
                    (-((MV2*shat - smT2)*(MV2**2 - MV2*shat + 2*smT2)* &
                         SP(3)) - 4*MV2*shat*smT2*SP(4)) - &
                   E**(3*y)*smT2**3* &
                    (MV2*shat*(MV2*shat - smT2)*SP(3) - &
                      4*(5*MV2*shat - 3*smT2)*(3*MV2*shat + smT2)*SP(4)) &
                     - E**(7*y)*smT2*(-(MV2*shat) + smT2)* &
                    ((-56*MV2**3*shat**3 + 92*MV2**2*shat**2*smT2 - &
                         51*MV2*shat*smT2**2 + 10*smT2**3)*SP(3) + &
                      8*(-(MV2*shat) + smT2)**2*(12*MV2*shat + 5*smT2)* &
                       SP(4)) - &
                   E**(11*y)*smT2* &
                    (5*(MV2*shat - smT2)* &
                       (8*MV2**3*shat**3 - 20*MV2**2*shat**2*smT2 + &
                         15*MV2*shat*smT2**2 - 4*smT2**3)*SP(3) - &
                      4*(2*MV2*shat - smT2)*(MV2*shat + smT2)* &
                       (4*MV2**2*shat**2 - 11*MV2*shat*smT2 + &
                         3*smT2**2)*SP(4)) - &
                   2*E**(5*y)*smT2**2* &
                    (-((MV2*shat - smT2)* &
                         (7*MV2**2*shat**2 - 6*MV2*shat*smT2 + smT2**2)* &
                         SP(3)) + &
                      5*(3*MV2*shat + smT2)* &
                       (6*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
                       SP(4)) + &
                   E**(4*y)*smT**5* &
                    (-(shat*(MV2*(-MV2 + shat) - 2*smT2)* &
                         (MV2*shat - smT2)*SP(3)) + &
                      4*(-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
                         MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
                         5*(2*MT2 + MV2 + 3*shat)*smT2**2)*SP(4)) + &
                   E**(12*y)*smT**3* &
                    (-(shat*(MV2*shat - smT2)* &
                         (2*MV2**2*shat*(-5*MV2 + 3*shat) + &
                           MV2*(5*MV2 - 17*shat)*smT2 + 6*smT2**2)*SP(3) &
                         ) + &
                      4*(-2*MV2**3*shat**3*(2*MT2 + MV2 + shat) + &
                         2*MV2**2*shat**2*(6*MT2 + 3*MV2 + 8*shat)* &
                          smT2 - &
                         MV2*shat*(10*MT2 + 7*MV2 + 10*shat)*smT2**2 + &
                         (2*MT2 + MV2 + shat)*smT2**3)*SP(4)) + &
                   E**(6*y)*smT**3* &
                    (shat*(MV2*shat - smT2)* &
                       (2*MV2**2*shat*(-5*MV2 + 3*shat) + &
                         MV2*(5*MV2 - 17*shat)*smT2 + 6*smT2**2)*SP(3) &
                       - 20*(2*MV2**3*shat**3* &
                          (6*MT2 - 3*MV2 + 5*shat) - &
                         8*MV2**2*shat**2*(3*MT2 + 2*shat)*smT2 + &
                         MV2*shat*(16*MT2 + 4*MV2 + 13*shat)*smT2**2 - &
                         2*(2*MT2 + MV2 + 2*shat)*smT2**3)*SP(4)) + &
                   2*E**(9*y)* &
                    (2*(2*MV2*shat - smT2)*(-(MV2*shat) + smT2)**2* &
                       (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 5*smT2**2)* &
                       SP(3) + &
                      smT2*(64*MV2**4*shat**4 - &
                         72*MV2**3*shat**3*smT2 + &
                         12*MV2**2*shat**2*smT2**2 + &
                         31*MV2*shat*smT2**3 - 15*smT2**4)*SP(4)) + &
                   2*E**(8*y)*smT* &
                    (-(shat*(MV2*shat - smT2)* &
                         (4*MV2**3*shat**2*(-3*MV2 + shat) + &
                           3*MV2**2*(3*MV2 - 5*shat)*shat*smT2 - &
                           2*MV2*(MV2 - 4*shat)*smT2**2 - 2*smT2**3)* &
                         SP(3)) + &
                      2*(4*MV2**4*shat**4*(12*MT2 - 3*MV2 + 5*shat) - &
                         2*MV2**3*shat**3*(62*MT2 + 7*MV2 + 15*shat)* &
                          smT2 + &
                         12*MV2**2*shat**2*(12*MT2 + 2*MV2 + 5*shat)* &
                          smT2**2 - &
                         2*MV2*shat*(44*MT2 + 14*MV2 + 25*shat)* &
                          smT2**3 + 5*(4*MT2 + 2*MV2 + 3*shat)*smT2**4)* &
                       SP(4)) - &
                   2*E**(10*y)*smT* &
                    (-(shat*(MV2*shat - smT2)* &
                         (4*MV2**3*shat**2*(-3*MV2 + shat) + &
                           3*MV2**2*(3*MV2 - 5*shat)*shat*smT2 - &
                           2*MV2*(MV2 - 4*shat)*smT2**2 - 2*smT2**3)* &
                         SP(3)) - &
                      2*(-4*MV2**4*shat**4*(-4*MT2 + MV2 + shat) + &
                         2*MV2**3*(-34*MT2 + MV2 - 15*shat)*shat**3* &
                          smT2 + &
                         4*MV2**2*shat**2*(22*MT2 + 7*MV2 + 11*shat)* &
                          smT2**2 - &
                         MV2*shat*(46*MT2 + 21*MV2 + 27*shat)*smT2**3 + &
                         (10*MT2 + 5*MV2 + 6*shat)*smT2**4)*SP(4)) - &
                   2*E**(13*y)*smT2**2* &
                    (smT2**3*(5*SP(3) + SP(4)) - &
                      MV2**3*shat**3*(5*SP(3) + 2*SP(4)) + &
                      MV2**2*shat**2*smT2*(17*SP(3) + 6*SP(4)) - &
                      MV2*shat*smT2**2*(17*SP(3) + 9*SP(4)))))/ &
               (2d0*E**(9*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              (LC(9)*(-4*E**(16*y)*shat*smT**9*SP(3) + &
                   2*E**(15*y)*smT2**4*(3*MV2*shat + smT2)*SP(3) + &
                   2*shat*smT**7*(2*MV2*shat - smT2)*SP(4) + &
                   E**y*smT2**3* &
                    (-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2)* &
                    SP(4) + 2*E**(14*y)*smT**7* &
                    (2*(MV2*shat*(2*MT2 - MV2 + 11*shat) - &
                         (2*MT2 + MV2 + 6*shat)*smT2)*SP(3) - &
                      MV2*shat**2*SP(4)) + &
                   E**(2*y)*shat*smT**5* &
                    (4*MV2*shat*smT2*SP(3) + &
                      (MV2**2*(-4*MT2 + MV2 - 43*shat)*shat + &
                         MV2*(4*MT2 + 3*MV2 + 44*shat)*smT2 - 11*smT2**2 &
                         )*SP(4)) - &
                   E**(12*y)*smT**5* &
                    (4*(-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
                         MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
                         5*(2*MT2 + MV2 + 3*shat)*smT2**2)*SP(3) + &
                      shat*(MV2**2*(-4*MT2 + MV2 - 23*shat)*shat + &
                         MV2*(4*MT2 + 3*MV2 + 14*shat)*smT2 - smT2**2)* &
                       SP(4)) + &
                   2*E**(13*y)*smT2**3* &
                    (-2*(5*MV2*shat - 3*smT2)*(3*MV2*shat + smT2)* &
                       SP(3) + &
                      (2*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) &
                      ) + E**(4*y)*smT**3* &
                    (4*(2*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
                         2*MV2**2*shat**2*(6*MT2 + 3*MV2 + 8*shat)* &
                          smT2 + &
                         MV2*shat*(10*MT2 + 7*MV2 + 10*shat)*smT2**2 - &
                         (2*MT2 + MV2 + shat)*smT2**3)*SP(3) + &
                      shat*(2*MV2**3*shat**2* &
                          (20*MT2 - 5*MV2 + 74*shat) - &
                         MV2**2*shat*(60*MT2 + 25*MV2 + 229*shat)* &
                          smT2 + &
                         MV2*(20*MT2 + 15*MV2 + 126*shat)*smT2**2 - &
                         23*smT2**3)*SP(4)) + &
                   E**(10*y)*smT**3* &
                    (20*(2*MV2**3*shat**3*(6*MT2 - 3*MV2 + 5*shat) - &
                         8*MV2**2*shat**2*(3*MT2 + 2*shat)*smT2 + &
                         MV2*shat*(16*MT2 + 4*MV2 + 13*shat)*smT2**2 - &
                         2*(2*MT2 + MV2 + 2*shat)*smT2**3)*SP(3) + &
                      shat*(2*MV2**3*(-20*MT2 + 5*MV2 - 42*shat)* &
                          shat**2 + &
                         MV2**2*shat*(60*MT2 + 25*MV2 + 97*shat)*smT2 - &
                         MV2*(20*MT2 + 15*MV2 + 36*shat)*smT2**2 + &
                         smT2**3)*SP(4)) + &
                   E**(11*y)*smT2**2* &
                    (10*(3*MV2*shat + smT2)* &
                       (6*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
                       SP(3) + &
                      (-42*MV2**3*shat**3 + 48*MV2**2*shat**2*smT2 - &
                         33*MV2*shat*smT2**2 + 11*smT2**3)*SP(4)) - &
                   2*E**(8*y)*smT* &
                    (2*(4*MV2**4*shat**4*(12*MT2 - 3*MV2 + 5*shat) - &
                         2*MV2**3*shat**3*(62*MT2 + 7*MV2 + 15*shat)* &
                          smT2 + &
                         12*MV2**2*shat**2*(12*MT2 + 2*MV2 + 5*shat)* &
                          smT2**2 - &
                         2*MV2*shat*(44*MT2 + 14*MV2 + 25*shat)* &
                          smT2**3 + 5*(4*MT2 + 2*MV2 + 3*shat)*smT2**4)* &
                       SP(3) + &
                      shat*(12*MV2**4*(-4*MT2 + 2*MV2 - 3*shat)* &
                          shat**3 + &
                         MV2**3*shat**2*(84*MT2 + 3*MV2 + 34*shat)* &
                          smT2 - &
                         MV2**2*shat*(44*MT2 + 13*MV2 + 11*shat)* &
                          smT2**2 + &
                         2*MV2*(4*MT2 + 3*MV2 - 3*shat)*smT2**3 + &
                         4*smT2**4)*SP(4)) - &
                   2*E**(6*y)*smT* &
                    (2*(-4*MV2**4*shat**4*(-4*MT2 + MV2 + shat) + &
                         2*MV2**3*(-34*MT2 + MV2 - 15*shat)*shat**3* &
                          smT2 + &
                         4*MV2**2*shat**2*(22*MT2 + 7*MV2 + 11*shat)* &
                          smT2**2 - &
                         MV2*shat*(46*MT2 + 21*MV2 + 27*shat)*smT2**3 + &
                         (10*MT2 + 5*MV2 + 6*shat)*smT2**4)*SP(3) + &
                      shat*(4*MV2**4*shat**3* &
                          (12*MT2 - 6*MV2 + 17*shat) - &
                         MV2**3*shat**2*(84*MT2 + 3*MV2 + 130*shat)* &
                          smT2 + &
                         MV2**2*shat*(44*MT2 + 13*MV2 + 127*shat)* &
                          smT2**2 - &
                         MV2*(8*MT2 + 6*MV2 + 61*shat)*smT2**3 + &
                         11*smT2**4)*SP(4)) + &
                   E**(9*y)*smT2* &
                    (-8*(MV2*shat - smT2)**3*(12*MV2*shat + 5*smT2)* &
                       SP(3) + &
                      (148*MV2**4*shat**4 - 294*MV2**3*shat**3*smT2 + &
                         254*MV2**2*shat**2*smT2**2 - &
                         113*MV2*shat*smT2**3 + 25*smT2**4)*SP(4)) + &
                   E**(3*y)*smT2**2* &
                    (-4*MV2**3*shat**3*(SP(3) - 4*SP(4)) + &
                      smT2**3*(2*SP(3) + 7*SP(4)) + &
                      4*MV2**2*shat**2*smT2*(3*SP(3) + 7*SP(4)) - &
                      MV2*shat*smT2**2*(18*SP(3) + 35*SP(4))) + &
                   2*E**(5*y)*smT2* &
                    (MV2**3*shat**3*smT2*(36*SP(3) - 41*SP(4)) + &
                      2*smT2**4*(3*SP(3) + 5*SP(4)) - &
                      2*MV2**4*shat**4*(8*SP(3) + 9*SP(4)) + &
                      6*MV2**2*shat**2*smT2**2*(3*SP(3) + 16*SP(4)) - &
                      MV2*shat*smT2**3*(28*SP(3) + 57*SP(4))) - &
                   2*E**(7*y)* &
                    (6*MV2**2*shat**2*smT2**3*(2*SP(3) - 31*SP(4)) + &
                      8*MV2**4*shat**4*smT2*(8*SP(3) - 19*SP(4)) + &
                      40*MV2**5*shat**5*SP(4) - &
                      15*smT2**5*(SP(3) + SP(4)) + &
                      3*MV2**3*shat**3*smT2**2*(-24*SP(3) + 77*SP(4)) + &
                      MV2*shat*smT2**4*(31*SP(3) + 82*SP(4)))))/ &
               (2d0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              (LC(3)*(-2*E**(15*y)*MV2*shat*smT**7*(3*MV2*shat - smT2)* &
                    SP(3) + E**(16*y)*shat*(3*MV2*shat - smT2)*smT2**4* &
                    SP(3) - 4*E**y*MV2**2*shat**2*smT**7*SP(4) + &
                   2*MV2*shat**2*smT2**4*SP(4) - &
                   E**(14*y)*smT2**3* &
                    ((MV2**2*shat**2*(8*MT2 - 4*MV2 + 31*shat) - &
                         2*MV2*shat*(6*MT2 + 13*shat)*smT2 + &
                         (4*MT2 + 5*shat)*smT2**2)*SP(3) + &
                      MV2*shat**2*(3*MV2*shat - smT2)*SP(4)) + &
                   2*E**(3*y)*MV2**2*shat**2*smT**5* &
                    (2*smT2*SP(3) + 19*MV2*shat*SP(4) - 11*smT2*SP(4)) &
                    + 2*E**(7*y)*MV2*shat*smT* &
                    ((32*MV2**4*shat**4 - 46*MV2**3*shat**3*smT2 + &
                         36*MV2**2*shat**2*smT2**2 - &
                         17*MV2*shat*smT2**3 + 5*smT2**4)*SP(3) - &
                      10*MV2*shat*(3*MV2*shat - 2*smT2)* &
                       (MV2*shat - smT2)*smT2*SP(4)) + &
                   E**(2*y)*smT2**3* &
                    (-2*MV2**2*shat**3*SP(3) + &
                      (4*MV2**2*(-2*MT2 + MV2 - 5*shat)*shat**2 + &
                         3*MV2*shat*(4*MT2 + 3*shat)*smT2 + &
                         (-4*MT2 + shat)*smT2**2)*SP(4)) - &
                   E**(10*y)*smT2* &
                    (2*(4*MV2**4*shat**4*(30*MT2 - 15*MV2 + 11*shat) + &
                         12*MV2**3*(-25*MT2 + 5*MV2 - 6*shat)*shat**3* &
                          smT2 + &
                         4*MV2**2*shat**2*(70*MT2 - 5*MV2 + 18*shat)* &
                          smT2**2 - &
                         2*MV2*shat*(60*MT2 + 17*shat)*smT2**3 + &
                         5*(4*MT2 + shat)*smT2**4)*SP(3) + &
                      (80*MV2**4*shat**4*(MT2 + shat) - &
                         40*MV2**3*shat**3*(5*MT2 + MV2 + 2*shat)* &
                          smT2 + &
                         10*MV2**2*shat**2*(20*MT2 + 2*MV2 + shat)* &
                          smT2**2 + &
                         MV2*shat*(-100*MT2 + 17*shat)*smT2**3 - &
                         5*(-4*MT2 + shat)*smT2**4)*SP(4)) - &
                   E**(6*y)*smT2* &
                    ((16*MV2**4*shat**4*(5*MT2 + 4*shat) - &
                         8*MV2**3*shat**3*(25*MT2 + 5*MV2 + 8*shat)* &
                          smT2 + &
                         MV2**2*shat**2*(200*MT2 + 20*MV2 + 23*shat)* &
                          smT2**2 - &
                         2*MV2*shat*(50*MT2 + shat)*smT2**3 + &
                         (20*MT2 + shat)*smT2**4)*SP(3) + &
                      5*(8*MV2**4*shat**4*(6*MT2 - 3*MV2 + shat) + &
                         24*MV2**3*(-5*MT2 + MV2)*shat**3*smT2 + &
                         MV2**2*(112*MT2 - 8*MV2 - 5*shat)*shat**2* &
                          smT2**2 + &
                         MV2*shat*(-48*MT2 + 5*shat)*smT2**3 - &
                         2*(-4*MT2 + shat)*smT2**4)*SP(4)) - &
                   2*E**(11*y)*MV2*shat*smT**3* &
                    (-10*smT2**3*SP(3) + &
                      7*MV2*shat*smT2**2*(10*SP(3) - SP(4)) + &
                      5*MV2**2*shat**2*smT2*(-28*SP(3) + SP(4)) + &
                      10*MV2**3*shat**3*(9*SP(3) + SP(4))) + &
                   2*E**(13*y)*MV2*shat*smT**5* &
                    (5*smT2**2*SP(3) + &
                      MV2*shat*smT2*(-27*SP(3) + SP(4)) + &
                      MV2**2*shat**2*(30*SP(3) + SP(4))) - &
                   2*E**(5*y)*MV2*shat*smT**3* &
                    (-(smT2**3*SP(3)) + &
                      MV2**2*shat**2*smT2*(14*SP(3) - 65*SP(4)) + &
                      MV2*shat*smT2**2*(-7*SP(3) + 25*SP(4)) + &
                      2*MV2**3*shat**3*(SP(3) + 25*SP(4))) + &
                   4*E**(9*y)*MV2*shat*smT* &
                    (5*smT2**4*SP(3) + &
                      MV2*shat*smT2**3*(-34*SP(3) + SP(4)) + &
                      24*MV2**4*shat**4*(SP(3) + SP(4)) + &
                      MV2**2*shat**2*smT2**2*(72*SP(3) + 7*SP(4)) - &
                      MV2**3*shat**3*smT2*(67*SP(3) + 27*SP(4))) + &
                   E**(12*y)*smT2**2* &
                    (shat*(20*MV2**3*shat**2*(-2*MV2 + 5*shat) + &
                         MV2**2*(20*MV2 - 133*shat)*shat*smT2 + &
                         65*MV2*shat*smT2**2 - 10*smT2**3)*SP(3) + &
                      shat*(28*MV2**3*shat**3 - &
                         MV2**2*shat*(4*MV2 + 19*shat)*smT2 + smT2**3)* &
                       SP(4) + &
                      4*MT2*(MV2*shat - smT2)* &
                       (5*(-2*MV2*shat + smT2)**2*SP(3) + &
                         (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
                          SP(4))) + &
                   E**(4*y)*smT2**2* &
                    (4*MT2*(MV2*shat - smT2)* &
                       ((2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
                          SP(3) + 5*(-2*MV2*shat + smT2)**2*SP(4)) + &
                      shat*(MV2*shat* &
                          (4*MV2**2*(5*shat**2 - smT2) - &
                            9*MV2*shat*smT2 - smT2**2)*SP(3) + &
                         (20*MV2**3*shat**2*(-2*MV2 + 3*shat) + &
                            MV2**2*(20*MV2 - 47*shat)*shat*smT2 + &
                            4*MV2*shat*smT2**2 + 5*smT2**3)*SP(4))) + &
                   E**(8*y)*(8*MT2*(MV2*shat - smT2)* &
                       (24*MV2**4*shat**4 - 48*MV2**3*shat**3*smT2 + &
                         49*MV2**2*shat**2*smT2**2 - &
                         25*MV2*shat*smT2**3 + 5*smT2**4)* &
                       (SP(3) + SP(4)) + &
                      shat*(MV2*shat*smT2**4*(29*SP(3) - 38*SP(4)) - &
                         5*smT2**5*(SP(3) - 2*SP(4)) + &
                         58*MV2**2*shat**2*smT2**3*(-SP(3) + SP(4)) - &
                         48*MV2**6*shat**4*(SP(3) + SP(4)) + &
                         8*MV2**5*shat**3*(2*shat**2 + 3*smT2)* &
                          (SP(3) + SP(4)) + &
                         8*MV2**3*shat*smT2**2* &
                          (shat**2*(5*SP(3) - 7*SP(4)) + &
                            smT2*(SP(3) + SP(4))) + &
                         8*MV2**4*shat**2*smT2* &
                          (-3*smT2*(SP(3) + SP(4)) + &
                            shat**2*(SP(3) + 5*SP(4)))))))/ &
               (E**(8*y)*shat**2*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              (2*(E**(5*y)*shat*smT**3 + &
                   E**y*shat*smT*(2*MV2*shat - smT2) - &
                   2*E**(3*y)*shat*smT*(2*MV2*shat - smT2) - &
                   E**(4*y)*MV2*shat*smT2 + smT2*(-2*MV2*shat + smT2) + &
                   E**(2*y)*(4*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
                      smT2**2))*LC(5)*SP(5))/ &
               (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
              ((-2*E**(14*y)*shat*smT**7*(3*MV2*shat - 2*smT2) + &
                   2*shat*smT**7*(2*MV2*shat - smT2) - &
                   4*E**y*MV2**2*shat**2*smT2**3 + &
                   E**(13*y)*smT2**3* &
                    (2*MV2**2*shat**2 + 3*MV2*shat*smT2 - smT2**2) - &
                   E**(2*y)*shat*smT**5* &
                    (-2*MV2**2*(-2*MT2 + MV2 - 22*shat)*shat - &
                      MV2*(4*MT2 + 2*MV2 + 43*shat)*smT2 + 9*smT2**2) + &
                   E**(12*y)*shat*smT**5* &
                    (-2*MV2**2*(-2*MT2 + MV2 - 32*shat)*shat - &
                      MV2*(4*MT2 + 2*MV2 + 73*shat)*smT2 + 19*smT2**2) &
                    - E**(10*y)*shat*smT**3* &
                    (2*MV2**3*shat**2*(20*MT2 - 10*MV2 + 109*shat) - &
                      2*MV2**2*shat*(30*MT2 + 5*MV2 + 172*shat)*smT2 + &
                      MV2*(20*MT2 + 10*MV2 + 183*shat)*smT2**2 - &
                      35*smT2**3) + &
                   E**(4*y)*shat*smT**3* &
                    (2*MV2**3*shat**2*(20*MT2 - 10*MV2 + 77*shat) - &
                      2*MV2**2*shat*(30*MT2 + 5*MV2 + 106*shat)*smT2 + &
                      MV2*(20*MT2 + 10*MV2 + 93*shat)*smT2**2 - &
                      13*smT2**3) - &
                   E**(3*y)*smT2**2* &
                    (-42*MV2**3*shat**3 + 28*MV2**2*shat**2*smT2 - &
                      3*MV2*shat*smT2**2 + smT2**3) - &
                   E**(11*y)*smT2**2* &
                    (20*MV2**3*shat**3 + 16*MV2**2*shat**2*smT2 - &
                      25*MV2*shat*smT2**2 + 5*smT2**3) - &
                   2*E**(6*y)*shat*smT* &
                    (12*MV2**4*shat**3*(4*MT2 - 3*MV2 + 6*shat) + &
                      MV2**3*(-84*MT2 + 18*MV2 - 85*shat)*shat**2* &
                       smT2 + &
                      2*MV2**2*shat*(22*MT2 + MV2 + 19*shat)*smT2**2 - &
                      MV2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3 - smT2**4) - &
                   10*E**(9*y)*smT2* &
                    (-6*MV2**4*shat**4 + MV2**3*shat**3*smT2 + &
                      8*MV2**2*shat**2*smT2**2 - 6*MV2*shat*smT2**3 + &
                      smT2**4) - &
                   5*E**(5*y)*smT2* &
                    (28*MV2**4*shat**4 - 42*MV2**3*shat**3*smT2 + &
                      22*MV2**2*shat**2*smT2**2 - 5*MV2*shat*smT2**3 + &
                      smT2**4) + &
                   2*E**(7*y)*(MV2*shat - smT2)* &
                    (24*MV2**4*shat**4 - 48*MV2**3*shat**3*smT2 + &
                      49*MV2**2*shat**2*smT2**2 - 25*MV2*shat*smT2**3 + &
                      5*smT2**4) + &
                   2*E**(8*y)*shat*smT* &
                    (4*MV2**4*shat**3*(12*MT2 - 9*MV2 + 26*shat) + &
                      MV2**3*(-84*MT2 + 18*MV2 - 181*shat)*shat**2* &
                       smT2 + &
                      2*MV2**2*shat*(22*MT2 + MV2 + 77*shat)*smT2**2 - &
                      4*MV2*(2*MT2 + MV2 + 19*shat)*smT2**3 + 14*smT2**4 &
                      ))*LC(4)*SP(6))/ &
               (2d0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((E**(5*y)*(-32*MV2**4*shat**4*smT - &
                      12*MV2**3*shat**3*smT**3 + &
                      48*MV2**2*shat**2*smT**5 - 27*MV2*shat*smT**7 + &
                      3*smT**9) + 2*shat*smT2**4 - &
                   2*E**(14*y)*shat*smT2**4 - &
                   E**y*smT**7*(3*MV2*shat + smT2) + &
                   E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
                   2*E**(11*y)*smT**5*(3*MV2*shat + smT2)* &
                    (-5*MV2*shat + 3*smT2) + &
                   4*E**(7*y)*smT*(-(MV2*shat) + smT2)**3* &
                    (8*MV2*shat + 3*smT2) + &
                   2*E**(2*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
                      (2*MT2 + MV2 + 5*shat)*smT2) - &
                   2*E**(12*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
                      (2*MT2 + MV2 + 5*shat)*smT2) - &
                   2*E**(10*y)*smT2**2* &
                    (-10*MV2**2*(-2*MT2 + MV2 - 3*shat)*shat**2 - &
                      5*MV2*shat*(6*MT2 + MV2 + 6*shat)*smT2 + &
                      (10*MT2 + 5*MV2 + 11*shat)*smT2**2) - &
                   2*E**(3*y)*smT**3* &
                    (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
                      10*MV2*shat*smT2**2 + smT2**3) + &
                   E**(9*y)*smT**3* &
                    (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
                      11*MV2*shat*smT2**2 + 13*smT2**3) + &
                   2*E**(4*y)*smT2* &
                    (4*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
                      2*MV2**2*(2*MT2 + 11*MV2 - 9*shat)*shat**2*smT2 + &
                      MV2*(-6*MT2 + 7*MV2 - 18*shat)*shat*smT2**2 + &
                      (2*MT2 + MV2 + 7*shat)*smT2**3) + &
                   2*E**(8*y)*smT2* &
                    (8*MV2**3*shat**3*(7*MT2 - 4*MV2 + 2*shat) + &
                      6*MV2**2*(-18*MT2 + MV2 - 3*shat)*shat**2*smT2 + &
                      2*MV2*shat*(34*MT2 + 7*(MV2 + 2*shat))*smT2**2 - &
                      (16*MT2 + 8*MV2 + 11*shat)*smT2**3) - &
                   2*E**(6*y)* &
                    (-8*MV2**4*shat**4*(-4*MT2 + MV2 + shat) - &
                      4*MV2**3*(16*MT2 + 5*MV2 - 7*shat)*shat**3*smT2 + &
                      6*MV2**2*(10*MT2 + 3*MV2 - shat)*shat**2* &
                       smT2**2 - 2*MV2*(18*MT2 + 7*MV2)*shat*smT2**3 + &
                      (8*MT2 + 4*MV2 + shat)*smT2**4))*LC(10)*SP(6))/ &
               (E**(6*y)*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((-2*E**(14*y)*MV2*shat**2*smT2**3 + 2*shat*smT2**4 - &
                   E**y*smT**7*(3*MV2*shat + smT2) + &
                   E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
                   2*E**(3*y)*smT**5*(5*MV2*shat - 3*smT2)* &
                    (3*MV2*shat + smT2) + &
                   4*E**(7*y)*smT*(MV2*shat - smT2)**3* &
                    (8*MV2*shat + 3*smT2) + &
                   2*E**(2*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 11*shat)*shat + &
                      (2*MT2 + MV2 + 6*shat)*smT2) - &
                   2*E**(12*y)*smT2**2* &
                    (-10*MV2**2*shat**3 + &
                      MV2*shat*(-2*MT2 + MV2 + 4*shat)*smT2 + &
                      (2*MT2 + MV2 + shat)*smT2**2) + &
                   2*E**(4*y)*smT2**2* &
                    (-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
                      MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
                      5*(2*MT2 + MV2 + 3*shat)*smT2**2) + &
                   2*E**(11*y)*smT**3* &
                    (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
                      10*MV2*shat*smT2**2 + smT2**3) - &
                   E**(5*y)*smT**3* &
                    (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
                      11*MV2*shat*smT2**2 + 13*smT2**3) - &
                   2*E**(10*y)*smT2* &
                    (4*MV2**3*shat**3*(2*MT2 + MV2 + 9*shat) - &
                      2*MV2**2*shat**2*(2*MT2 + 11*MV2 + 19*shat)* &
                       smT2 + &
                      MV2*shat*(-6*MT2 + 7*MV2 + 13*shat)*smT2**2 + &
                      (2*MT2 + MV2)*smT2**3) + &
                   2*E**(6*y)*smT2* &
                    (8*MV2**3*(-7*MT2 + 4*MV2 - 6*shat)*shat**3 + &
                      2*MV2**2*shat**2*(54*MT2 - 3*MV2 + 37*shat)* &
                       smT2 - &
                      MV2*shat*(68*MT2 + 14*MV2 + 59*shat)*smT2**2 + &
                      2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3) + &
                   E**(9*y)*smT* &
                    (32*MV2**4*shat**4 + 12*MV2**3*shat**3*smT2 - &
                      48*MV2**2*shat**2*smT2**2 + 27*MV2*shat*smT2**3 - &
                      3*smT2**4) + &
                   2*E**(8*y)* &
                    (-8*MV2**4*(-4*MT2 + MV2 - 3*shat)*shat**4 - &
                      4*MV2**3*shat**3*(16*MT2 + 5*MV2 + 9*shat)*smT2 + &
                      6*MV2**2*shat**2*(10*MT2 + 3*MV2 + 9*shat)* &
                       smT2**2 - &
                      2*MV2*shat*(18*MT2 + 7*MV2 + 18*shat)*smT2**3 + &
                      (8*MT2 + 4*MV2 + 9*shat)*smT2**4))*LC(14)*SP(6))/ &
               (E**(8*y)*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((-2*shat*smT**9 + &
                   E**(14*y)*shat*smT**7*(MV2*shat + smT2) + &
                   E**(13*y)*smT2**3* &
                    (-(MV2**2*shat**2) - 4*MV2*shat*smT2 + smT2**2) + &
                   E**y*smT2**3* &
                    (6*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2) - &
                   2*E**(2*y)*shat*smT**5* &
                    (MV2**2*shat*(-2*MT2 + MV2 + shat) + &
                      MV2*(2*MT2 + MV2 - 10*shat)*smT2 + 4*smT2**2) + &
                   E**(4*y)*shat*smT**3* &
                    (2*MV2**3*shat**2*(-20*MT2 + 10*MV2 + 11*shat) + &
                      2*MV2**2*(30*MT2 + 5*MV2 - 38*shat)*shat*smT2 - &
                      5*MV2*(4*MT2 + 2*MV2 - 9*shat)*smT2**2 - &
                      13*smT2**3) + &
                   2*E**(10*y)*shat*smT**3* &
                    (5*MV2**3*shat**2*(4*MT2 - 2*MV2 + shat) - &
                      5*MV2**2*(6*MT2 + MV2 - shat)*shat*smT2 + &
                      5*MV2*(2*MT2 + MV2)*smT2**2 + smT2**3) + &
                   E**(11*y)*smT2**2* &
                    (8*MV2**3*shat**3 + 35*MV2**2*shat**2*smT2 - &
                      33*MV2*shat*smT2**2 + 6*smT2**3) + &
                   E**(3*y)*smT2**2* &
                    (-62*MV2**3*shat**3 + 67*MV2**2*shat**2*smT2 - &
                      27*MV2*shat*smT2**2 + 6*smT2**3) + &
                   E**(6*y)*shat*smT* &
                    (-8*MV2**4*shat**3*(-12*MT2 + 9*MV2 + 10*shat) + &
                      2*MV2**3*shat**2*(-84*MT2 + 18*MV2 + 83*shat)* &
                       smT2 + &
                      4*MV2**2*(22*MT2 + MV2 - 22*shat)*shat*smT2**2 + &
                      MV2*(-8*(2*MT2 + MV2) + 43*shat)*smT2**3 - &
                      11*smT2**4) + &
                   2*E**(8*y)*shat*smT* &
                    (12*MV2**4*shat**3*(-4*MT2 + 3*MV2 + 2*shat) + &
                      MV2**3*(84*MT2 - 18*MV2 - 35*shat)*shat**2*smT2 - &
                      2*MV2**2*shat*(22*MT2 + MV2 + 7*shat)*smT2**2 + &
                      4*MV2*(2*MT2 + MV2 + 3*shat)*smT2**3 - 2*smT2**4) &
                    + 2*E**(7*y)*(-(MV2*shat) + smT2)* &
                    (56*MV2**4*shat**4 - 104*MV2**3*shat**3*smT2 + &
                      105*MV2**2*shat**2*smT2**2 - &
                      52*MV2*shat*smT2**3 + 10*smT2**4) + &
                   E**(9*y)*smT2* &
                    (-12*MV2**4*shat**4 - 94*MV2**3*shat**3*smT2 + &
                      164*MV2**2*shat**2*smT2**2 - &
                      93*MV2*shat*smT2**3 + 15*smT2**4) + &
                   E**(5*y)*smT2* &
                    (204*MV2**4*shat**4 - 362*MV2**3*shat**3*smT2 + &
                      247*MV2**2*shat**2*smT2**2 - &
                      84*MV2*shat*smT2**3 + 15*smT2**4) - &
                   E**(12*y)*shat*smT**5* &
                    (4*MT2*MV2*(MV2*shat - smT2) - &
                      (MV2*shat + smT2)*(2*MV2*(MV2 - 4*shat) + 3*smT2)) &
                   )*LC(2)*SP(7))/ &
               (2d0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((-2*MV2*shat**2*smT2**3 + 2*E**(14*y)*shat*smT2**4 + &
                   E**y*smT**7*(3*MV2*shat + smT2) - &
                   E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
                   2*E**(11*y)*smT**5*(5*MV2*shat - 3*smT2)* &
                    (3*MV2*shat + smT2) + &
                   4*E**(7*y)*smT*(MV2*shat - smT2)**3* &
                    (8*MV2*shat + 3*smT2) + &
                   2*E**(12*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 11*shat)*shat + &
                      (2*MT2 + MV2 + 6*shat)*smT2) - &
                   2*E**(2*y)*smT2**2* &
                    (-10*MV2**2*shat**3 + &
                      MV2*shat*(-2*MT2 + MV2 + 4*shat)*smT2 + &
                      (2*MT2 + MV2 + shat)*smT2**2) + &
                   2*E**(10*y)*smT2**2* &
                    (-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
                      MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
                      5*(2*MT2 + MV2 + 3*shat)*smT2**2) + &
                   2*E**(3*y)*smT**3* &
                    (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
                      10*MV2*shat*smT2**2 + smT2**3) - &
                   E**(9*y)*smT**3* &
                    (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
                      11*MV2*shat*smT2**2 + 13*smT2**3) - &
                   2*E**(4*y)*smT2* &
                    (4*MV2**3*shat**3*(2*MT2 + MV2 + 9*shat) - &
                      2*MV2**2*shat**2*(2*MT2 + 11*MV2 + 19*shat)* &
                       smT2 + &
                      MV2*shat*(-6*MT2 + 7*MV2 + 13*shat)*smT2**2 + &
                      (2*MT2 + MV2)*smT2**3) + &
                   2*E**(8*y)*smT2* &
                    (8*MV2**3*(-7*MT2 + 4*MV2 - 6*shat)*shat**3 + &
                      2*MV2**2*shat**2*(54*MT2 - 3*MV2 + 37*shat)* &
                       smT2 - &
                      MV2*shat*(68*MT2 + 14*MV2 + 59*shat)*smT2**2 + &
                      2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3) + &
                   E**(5*y)*smT* &
                    (32*MV2**4*shat**4 + 12*MV2**3*shat**3*smT2 - &
                      48*MV2**2*shat**2*smT2**2 + 27*MV2*shat*smT2**3 - &
                      3*smT2**4) + &
                   2*E**(6*y)* &
                    (-8*MV2**4*(-4*MT2 + MV2 - 3*shat)*shat**4 - &
                      4*MV2**3*shat**3*(16*MT2 + 5*MV2 + 9*shat)*smT2 + &
                      6*MV2**2*shat**2*(10*MT2 + 3*MV2 + 9*shat)* &
                       smT2**2 - &
                      2*MV2*shat*(18*MT2 + 7*MV2 + 18*shat)*smT2**3 + &
                      (8*MT2 + 4*MV2 + 9*shat)*smT2**4))*LC(8)*SP(7))/ &
               (E**(6*y)*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
              ((E**(9*y)*(-32*MV2**4*shat**4*smT - &
                      12*MV2**3*shat**3*smT**3 + &
                      48*MV2**2*shat**2*smT**5 - 27*MV2*shat*smT**7 + &
                      3*smT**9) - 2*shat*smT2**4 + &
                   2*E**(14*y)*shat*smT2**4 + &
                   E**y*smT**7*(3*MV2*shat + smT2) - &
                   E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
                   2*E**(3*y)*smT**5*(3*MV2*shat + smT2)* &
                    (-5*MV2*shat + 3*smT2) + &
                   4*E**(7*y)*smT*(-(MV2*shat) + smT2)**3* &
                    (8*MV2*shat + 3*smT2) - &
                   2*E**(2*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
                      (2*MT2 + MV2 + 5*shat)*smT2) + &
                   2*E**(12*y)*smT2**3* &
                    (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
                      (2*MT2 + MV2 + 5*shat)*smT2) - &
                   2*E**(4*y)*smT2**2* &
                    (-10*MV2**2*(-2*MT2 + MV2 - 3*shat)*shat**2 - &
                      5*MV2*shat*(6*MT2 + MV2 + 6*shat)*smT2 + &
                      (10*MT2 + 5*MV2 + 11*shat)*smT2**2) - &
                   2*E**(11*y)*smT**3* &
                    (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
                      10*MV2*shat*smT2**2 + smT2**3) + &
                   E**(5*y)*smT**3* &
                    (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
                      11*MV2*shat*smT2**2 + 13*smT2**3) + &
                   2*E**(10*y)*smT2* &
                    (4*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
                      2*MV2**2*(2*MT2 + 11*MV2 - 9*shat)*shat**2*smT2 + &
                      MV2*(-6*MT2 + 7*MV2 - 18*shat)*shat*smT2**2 + &
                      (2*MT2 + MV2 + 7*shat)*smT2**3) + &
                   2*E**(6*y)*smT2* &
                    (8*MV2**3*shat**3*(7*MT2 - 4*MV2 + 2*shat) + &
                      6*MV2**2*(-18*MT2 + MV2 - 3*shat)*shat**2*smT2 + &
                      2*MV2*shat*(34*MT2 + 7*(MV2 + 2*shat))*smT2**2 - &
                      (16*MT2 + 8*MV2 + 11*shat)*smT2**3) - &
                   2*E**(8*y)* &
                    (-8*MV2**4*shat**4*(-4*MT2 + MV2 + shat) - &
                      4*MV2**3*(16*MT2 + 5*MV2 - 7*shat)*shat**3*smT2 + &
                      6*MV2**2*(10*MT2 + 3*MV2 - shat)*shat**2* &
                       smT2**2 - 2*MV2*(18*MT2 + 7*MV2)*shat*smT2**3 + &
                      (8*MT2 + 4*MV2 + shat)*smT2**4))*LC(12)*SP(7))/ &
               (E**(8*y)*(MV2*shat - smT2)**3* &
                 (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
           (8*kapPr*(MV2*shat - smT2*Cosh(2*y))*(SP(3) + SP(4))* &
              (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)))/ &
            (-(MV2*shat) + smT2)**2) + &
        SI(16)*(kap*(((-2*E**(2*y)*shat*smT2* &
                    (2*MT2*shat - MV2*shat + smT2) + &
                   2*E**y*smT*(-(MV2*shat) + smT2)* &
                    (shat*(12*MT2 - 3*MV2 + shat) + 4*smT2) + &
                   4*shat*(MT2*shat*(4*MV2*shat - 3*smT2) - &
                      (-(MV2*shat) + smT2)**2))*LC(1))/ &
               (shat**2*(MV2*shat - smT2)) + &
              ((shat*smT2*(-((2*MT2 + MV2)*shat) + smT2) - &
                   2*E**y*smT* &
                    (MV2*(-4*MT2 + MV2)*shat**2 + &
                      (12*MT2 - 5*MV2)*shat*smT2 + 4*smT2**2) - &
                   2*E**(2*y)*shat* &
                    (MT2*shat*(4*MV2*shat - 5*smT2) - &
                      (-(MV2*shat) + smT2)**2))*LC(11))/ &
               (E**(2*y)*shat**2*(MV2*shat - smT2)) + &
              smT*LC(15)*(E**(-y) + &
                 (8*MT2*shat*smT - 2*MV2*shat*smT + 2*smT**3 + &
                    4*MT2*shat**2*Sinh(y))/(MV2*shat**2 - shat*smT2)) + &
              (4*LC(7)*(2*MT2*shat*(-2*MV2*shat + smT2) + &
                   (-(MV2*shat) + smT2)**2 + &
                   shat*smT*(2*MT2*shat*Cosh(y) + &
                      (2*MT2*shat - MV2*shat + smT2)*Sinh(y)))*SP(1))/ &
               (shat**2*(MV2*shat - smT2)) - &
              (4*LC(6)*(2*MT2*shat*(-2*MV2*shat + smT2) + &
                   (-(MV2*shat) + smT2)**2 + &
                   shat*smT*(2*MT2*shat*Cosh(y) + &
                      (2*MT2*shat - MV2*shat + smT2)*Sinh(y)))*SP(2))/ &
               (shat**2*(MV2*shat - smT2)) + &
              (4*smT2*(4*MT2*shat - MV2*shat + smT2)*LC(13)*SP(4))/ &
               (E**(2*y)*shat**2*(MV2*shat - smT2)) - &
              (2*smT*LC(9)*(2*E**(3*y)*smT* &
                    (4*MT2*shat - MV2*shat + smT2)*SP(3) - &
                   2*E**(2*y)*MT2*shat**2*SP(4) + &
                   shat*(2*MT2*shat - MV2*shat + smT2)*SP(4)))/ &
               (E**y*shat**2*(MV2*shat - smT2)) + &
              (4*LC(3)*(E**(3*y)*shat*smT2* &
                    (2*MT2*shat - MV2*shat + smT2)*SP(3) - &
                   2*E**(2*y)*smT* &
                    ((-(MV2*shat) + smT2)**2 + &
                      MT2*shat*(-4*MV2*shat + 3*smT2))*SP(3) + &
                   smT*(-(MV2*shat**2*(-8*MT2 + 2*MV2 + shat)) + &
                      shat*(-6*MT2 + 4*MV2 + shat)*smT2 - 2*smT2**2)* &
                    SP(4) + E**y*shat* &
                    (2*MT2*shat*(2*MV2*shat - smT2) - &
                      (-(MV2*shat) + smT2)**2)*SP(4)))/ &
               (E**y*shat**3*(MV2*shat - smT2)) + &
              (2*smT*(-2*MT2*shat + &
                   E**(2*y)*(10*MT2*shat - MV2*shat + smT2))*LC(4)*SP(6) &
                 )/(E**y*shat*(MV2*shat - smT2)) + &
              (4*(E**y*shat + smT)*(4*MT2*shat - MV2*shat + smT2)* &
                 LC(14)*SP(6))/(E**y*shat*(MV2*shat - smT2)) + &
              (2*smT*(-6*MT2*shat - &
                   E**(2*y)*(2*MT2*shat - MV2*shat + smT2))*LC(2)*SP(7)) &
                /(E**y*shat*(MV2*shat - smT2)) + &
              (4*(4*MT2*shat - MV2*shat + smT2)*LC(8)*SP(7))/ &
               (MV2*shat - smT2) - &
              (4*smT*(4*MT2*shat - MV2*shat + smT2)*LC(12)*SP(7))/ &
               (E**y*shat*(MV2*shat - smT2))) + &
           kapPr*(-2*E**y*smT*(SP(3)*SP(5) + SP(2)*SP(6)) - &
              (2*smT*(SP(4)*SP(5) + SP(1)*SP(7)))/E**y + &
              (4*(4*MT2*shat - MV2*shat + smT2)*(SP(3) + SP(4))* &
                 (MV2*shat*SP(5) - smT2*SP(5) - shat*SP(6)*SP(7)))/ &
               (shat*(MV2*shat - smT2)))) + &
        SI(12)*(kap*(((-2*shat*smT*(2*MV2*shat - 3*smT2)* &
                    (MV2*shat - smT2)*(2*MV2*shat - smT2) - &
                   E**(4*y)*shat**3*smT**3*(5*MV2*shat + 7*smT2) + &
                   E**(3*y)*shat**2*smT2* &
                    (8*MV2**2*shat**2 + 17*MV2*shat*smT2 + 11*smT2**2) &
                    + E**(2*y)*shat*smT* &
                    (MV2**2*(8*MT2 - 19*MV2 - 4*shat)*shat**3 + &
                      MV2*shat**2*(8*MT2 + 35*MV2 + 10*shat)*smT2 - &
                      2*shat*(8*MT2 + 32*MV2 + 3*shat)*smT2**2 + &
                      12*smT2**3) + &
                   E**y*(8*MV2**3*shat**4*(4*MT2 + MV2 + shat) - &
                      MV2**2*shat**3*(13*(8*MT2 + MV2) + 20*shat)* &
                       smT2 + &
                      MV2*shat**2*(88*MT2 + 5*MV2 + 12*shat)*smT2**2 + &
                      4*(-4*MT2 + 5*MV2)*shat*smT2**3 - 8*smT2**4))* &
                 LC(1))/(E**y*shat**2*(MV2*shat - smT2)**3) + &
              ((E**y*shat - smT)* &
                 (shat*(MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2 + &
                   E**(3*y)*shat**2*smT*(MV2*shat - 3*smT2)* &
                    (2*MV2*shat + smT2) + &
                   E**(2*y)*shat* &
                    (-4*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 + &
                      5*MV2*shat*smT2**2 + 2*smT2**3) - &
                   E**y*smT*(2*(2*MT2 - MV2)*MV2**2*shat**3 + &
                      MV2*shat**2*(-20*MT2 + 16*MV2 + shat)*smT2 + &
                      (16*MT2 - 16*MV2 - shat)*shat*smT2**2 + 8*smT2**3) &
                   )*LC(11))/(E**(3*y)*shat**2*(MV2*shat - smT2)**3) + &
              ((E**y*shat - smT)* &
                 (shat*smT**3*(MV2*shat - smT2) + &
                   2*E**(3*y)*shat**2*smT2*(MV2*shat + 4*smT2) - &
                   E**(2*y)*shat*smT*(3*MV2*shat + smT2)* &
                    (MV2*shat + 4*smT2) + &
                   E**y*(MV2**2*shat**3*(-4*MT2 + MV2 + 2*shat) + &
                      MV2*(-12*MT2 + 7*MV2 - 5*shat)*shat**2*smT2 + &
                      shat*(16*MT2 + 3*shat)*smT2**2 + 2*smT2**3))* &
                 LC(15))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(-3*E**(3*y)*MV2*shat**3*smT + E**(4*y)*shat**3*smT2 + &
                   shat*(MV2*shat - smT2)*smT2 + &
                   E**(2*y)*shat* &
                    (MV2*shat**2*(MV2 + 2*shat) + &
                      (3*MV2 - 2*shat)*shat*smT2 - smT2**2) + &
                   E**y*smT*(MV2*(MV2 - 3*shat)*shat**2 + &
                      shat*(-4*MV2 + 3*shat)*smT2 + 2*smT2**2))*LC(7)* &
                 SP(1))/(E**(2*y)*shat**2*(-(MV2*shat) + smT2)**2) + &
              (2*(E**(2*y)*shat**2*(-2*MV2 + shat) + &
                   2*E**(3*y)*shat**2*smT - shat*smT2 + &
                   2*E**y*smT*(-(MV2*shat) + smT2))*LC(6)*SP(2))/ &
               (E**(2*y)*shat**2*(MV2*shat - smT2)) + &
              (LC(13)*(2*E**(4*y)*shat**2*smT*(-(MV2*shat) + smT2)**2* &
                    SP(3) - 4*smT*smT2* &
                    (2*MV2*(-2*MT2 + MV2)*shat**2 + &
                      (4*MT2 - MV2)*shat*smT2 + smT2**2)*SP(4) + &
                   4*E**y*shat*smT2* &
                    (MV2*(-4*MT2 + 3*MV2)*shat**2 + &
                      2*(2*MT2 + MV2)*shat*smT2 + smT2**2)*SP(4) - &
                   2*E**(3*y)*shat* &
                    ((2*MV2*shat + shat**2 - smT2)* &
                       (-(MV2*shat) + smT2)**2*SP(3) - &
                      4*MV2*shat**3*smT2*SP(4)) + &
                   2*E**(2*y)*shat**2*smT* &
                    ((-(MV2*shat) + smT2)**2*SP(3) - &
                      2*MV2*shat*(MV2*shat + 5*smT2)*SP(4))))/ &
               (E**(3*y)*shat**2*(MV2*shat - smT2)**3) + &
              (LC(9)*(8*E**(7*y)*shat**4*smT2**2*SP(3) - &
                   4*E**(6*y)*shat**3*smT**3*(3*MV2*shat + 5*smT2)* &
                    SP(3) + 2*shat*smT**5*(MV2*shat - smT2)*SP(4) - &
                   2*E**y*shat**2*smT2* &
                    (MV2**2*shat*(4*MT2 - MV2 + shat) + &
                      MV2*(-4*MT2 - MV2 + shat)*smT2 - 2*smT2**2)*SP(4) &
                    - 4*E**(5*y)*shat**2*smT2* &
                    ((2*MV2*shat**2*(2*MT2 - MV2 + shat) - &
                         shat*(4*MT2 + 5*MV2 + 2*shat)*smT2 - 5*smT2**2) &
                        *SP(3) - MV2*shat**3*SP(4)) - &
                   2*E**(4*y)*shat*smT* &
                    (-2*(MV2**2*shat**3*(-6*MV2 + shat) + &
                         2*MV2*shat**2*(4*MT2 + 7*MV2 + 2*shat)*smT2 - &
                         shat*(8*MT2 + 19*MV2 + 5*shat)*smT2**2 + &
                         3*smT2**3)*SP(3) + &
                      shat**2* &
                       (4*MV2**2*shat**2 + 3*MV2*shat*smT2 + smT2**2)* &
                       SP(4)) + &
                   2*E**(3*y)* &
                    (2*(2*MV2**3*shat**4*(4*MT2 + 2*MV2 + shat) - &
                         2*MV2**2*shat**3*(12*MT2 + 5*MV2 + 4*shat)* &
                          smT2 + &
                         4*MV2*shat**2*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
                         shat*(-4*MT2 + MV2 + 2*shat)*smT2**3 - smT2**4) &
                        *SP(3) + &
                      MV2*shat**3* &
                       (-(MV2*shat**2*(4*MT2 - 3*MV2 + shat)) + &
                         shat*(4*MT2 + 3*MV2 + shat)*smT2 + 6*smT2**2)* &
                       SP(4)) - &
                   2*E**(2*y)*shat*smT* &
                    (2*(MV2*shat - smT2)* &
                       (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
                       SP(3) + &
                      (-2*MV2**2*shat**3*(4*MT2 - 2*MV2 + shat) + &
                         MV2*shat**2*(8*MT2 + shat)*smT2 + &
                         shat*(5*MV2 + shat)*smT2**2 - smT2**3)*SP(4)))) &
                /(E**(2*y)*shat**2*(E**y*shat - smT)* &
                 (MV2*shat - smT2)**3) + &
              (LC(3)*(-4*E**(6*y)*shat**4*smT**3*(3*MV2*shat - smT2)* &
                    SP(3) + 4*E**(5*y)*shat**3*smT2* &
                    (4*MV2**2*shat**2 + 7*MV2*shat*smT2 - 3*smT2**2)* &
                    SP(3) - 8*smT**3* &
                    (-(MV2**2*shat**3*(-4*MT2 + 2*MV2 + shat)) + &
                      MV2*shat**2*(-6*MT2 + 3*MV2 + 2*shat)*smT2 - &
                      shat*(-2*MT2 + 3*MV2 + shat)*smT2**2 + smT2**3)* &
                    SP(4) + 4*E**(4*y)*shat**2*smT* &
                    ((MV2**2*shat**3*(8*MT2 + 6*MV2 + 3*shat) - &
                         4*MV2*shat**2*(3*MT2 + 8*MV2 + shat)*smT2 + &
                         shat*(4*MT2 + 19*MV2 + shat)*smT2**2 - &
                         5*smT2**3)*SP(3) + &
                      MV2*shat**3*(3*MV2*shat - smT2)*SP(4)) - &
                   4*E**(2*y)*smT* &
                    ((-2*MV2**3*shat**4*(4*MT2 + 3*shat) + &
                         MV2**2*shat**3*(16*MT2 + 7*shat)*smT2 - &
                         2*MV2*shat**2*(6*MT2 - 3*MV2 + shat)*smT2**2 + &
                         shat*(4*MT2 - 6*MV2 + shat)*smT2**3 + 2*smT2**4 &
                         )*SP(3) + &
                      shat**3* &
                       (-2*MV2**2*shat**2*(-4*MT2 + 4*MV2 + shat) + &
                         MV2*shat*(-12*MT2 - 5*MV2 + 4*shat)*smT2 + &
                         (4*MT2 + MV2 - 2*shat)*smT2**2)*SP(4)) - &
                   4*E**y*shat*smT2* &
                    (2*MV2**2*shat**2*(MV2*shat - smT2)*SP(3) + &
                      (2*MV2**2*shat**3*(-8*MT2 + 5*MV2 + 2*shat) - &
                         8*MV2*shat**2*(-3*MT2 + MV2 + shat)*smT2 + &
                         shat*(-8*MT2 + 9*MV2 + 4*shat)*smT2**2 - &
                         3*smT2**3)*SP(4)) - &
                   4*E**(3*y)*shat* &
                    ((4*MV2**3*shat**4*(4*MT2 + MV2 + shat) - &
                         2*MV2**2*shat**3*(16*MT2 + 5*MV2 + shat)* &
                          smT2 - &
                         2*MV2*shat**2*(-12*MT2 + 5*MV2 + shat)* &
                          smT2**2 + (-8*MT2 + 13*MV2)*shat*smT2**3 - &
                         5*smT2**4)*SP(3) + &
                      shat**2* &
                       (2*MV2**3*shat**3 + 10*MV2**2*shat**2*smT2 - &
                         5*MV2*shat*smT2**2 + smT2**3)*SP(4))))/ &
               (E**(2*y)*shat**3*(E**y*shat - smT)*(MV2*shat - smT2)**3) &
                + (2*(-(E**(3*y)*shat**2*smT2) + &
                   shat*smT*(MV2**2 - MV2*shat + smT2) + &
                   E**(2*y)*shat*smT*(MV2*shat + 2*smT2) + &
                   E**y*(MV2*shat**2*(MV2 + shat) - &
                      shat*(5*MV2 + shat)*smT2 + smT2**2))*LC(5)*SP(5))/ &
               (E**y*shat*(-(MV2*shat) + smT2)**2) + &
              (2*(E**y*shat - smT)* &
                 (-(E**(2*y)*MV2*smT*(7*MV2*shat - 3*smT2)) + &
                   MV2*smT*(MV2*shat - smT2) + &
                   2*E**(3*y)*(3*MV2*shat - 2*smT2)*smT2 + &
                   E**y*(MV2**2*shat*(-4*MT2 + MV2 + shat) + &
                      MV2*(4*MT2 + MV2 - 3*shat)*smT2 + 2*smT2**2))* &
                 LC(4)*SP(6))/(E**(2*y)*(MV2*shat - smT2)**3) - &
              (4*(E**y*shat - smT)*smT* &
                 ((MV2 - E**y*smT)* &
                    (MV2*shat + smT*(-2*E**y*shat + smT)) + &
                   4*MT2*(-(MV2*shat) + smT2))*LC(10)*SP(6))/ &
               (-(MV2*shat) + smT2)**3 + &
              (4*(E**y*shat - smT)*smT* &
                 (2*MV2**2*shat**2 + smT**3*(E**y*shat + smT) + &
                   MV2*shat*(E**y*shat*(2*E**y*shat - 5*smT) - smT2) + &
                   4*MT2*shat*(-(MV2*shat) + smT2))*LC(14)*SP(6))/ &
               (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(3*MV2*shat*smT2*(-(MV2*shat) + smT2) - &
                   E**(4*y)*shat**2*smT2*(MV2*shat + smT2) + &
                   E**y*shat*smT* &
                    (MV2**2*shat*(-4*MT2 + 5*shat) + &
                      2*MV2*(2*MT2 + MV2 - 2*shat)*smT2 - smT2**2) + &
                   2*E**(3*y)*shat*smT* &
                    (MV2**2*shat**2 + MV2*shat*smT2 + smT2**2) + &
                   E**(2*y)*(-2*MV2**2*shat**3*(-2*MT2 + MV2 + shat) + &
                      MV2*shat**2*(-4*MT2 + 2*MV2 + shat)*smT2 + &
                      shat*(-7*MV2 + shat)*smT2**2 + smT2**3))*LC(2)* &
                 SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) - &
              (4*(E**y*shat - smT)*smT* &
                 (2*MV2**2*shat**2 - MV2*shat*smT*(3*E**y*shat + smT) + &
                   4*MT2*shat*(-(MV2*shat) + smT2) + &
                   smT2*(E**y*shat*(2*E**y*shat - smT) + smT2))*LC(12)* &
                 SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) - &
              (4*(E**y*shat - smT)*LC(8)* &
                 (2*MV2**2*shat**2 - &
                   E**(2*y)*smT**3*(smT - 4*shat*Cosh(y)) + &
                   E**y*smT*(MV2*(-4*MT2 + MV2 - 2*shat)*shat + &
                      (4*MT2 + MV2)*smT2 - &
                      MV2*shat*smT*(5*Cosh(y) + Sinh(y))))*SP(7))/ &
               (E**y*(MV2*shat - smT2)**3)) + &
           (4*kapPr*(E**y*shat - smT)* &
              (-((MV2*shat - smT2)* &
                   (smT2*(SP(3) + SP(4))*SP(5) + shat**2*SP(2)*SP(6))) &
                 + 2*shat*(MV2*shat*SP(3) + smT2*SP(4))*SP(6)*SP(7) + &
                E**y*shat*smT* &
                 ((MV2*shat - smT2)* &
                    ((SP(3) + SP(4))*SP(5) + SP(2)*SP(6)) - &
                   2*shat*(SP(3) + SP(4))*SP(6)*SP(7))))/ &
            (E**y*shat*(-(MV2*shat) + smT2)**2)) + &
        SI(9)*(kap*(((12*MV2*shat**3*smT2**2 + &
                   2*E**(3*y)*shat*smT*(MV2*shat - smT2)* &
                    (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2) + &
                   E**y*shat*smT* &
                    (4*MV2**3*shat**3 - 23*MV2**2*shat**2*smT2 - &
                      2*MV2*shat*smT2**2 - 3*smT2**3) + &
                   E**(2*y)*(-4*MV2**3*shat**5 + &
                      MV2**2*shat**3*(-8*MV2 + 7*shat)*smT2 + &
                      4*MV2*(9*MV2 - shat)*shat**2*smT2**2 + &
                      shat*(-24*MV2 + shat)*smT2**3 + 8*smT2**4))*LC(1)) &
                /(E**(2*y)*shat**2*(MV2*shat - smT2)**3) + &
              (smT2*(2*E**(3*y)*MV2*shat**2*smT*(MV2*shat - smT2) + &
                   6*shat**2*smT2**2 - &
                   E**y*shat*smT* &
                    (-2*MV2**2*shat**2 + 13*MV2*shat*smT2 + smT2**2) + &
                   E**(2*y)*(MV2**2*(MV2 - shat)*shat**3 - &
                      7*MV2**2*shat**2*smT2 + &
                      shat*(20*MV2 + shat)*smT2**2 - 8*smT2**3))*LC(11)) &
                /(E**(4*y)*shat**2*(MV2*shat - smT2)**3) + &
              (smT*(-2*shat**2*smT2*(2*MV2*shat + 3*smT2) + &
                   E**y*shat*smT* &
                    (3*MV2**2*shat**2 + 15*MV2*shat*smT2 + 2*smT2**2) - &
                   E**(2*y)*(MV2**2*(MV2 - 2*shat)*shat**3 + &
                      MV2*shat**2*(3*MV2 + 4*shat)*smT2 + &
                      2*(4*MV2 - shat)*shat*smT2**2 - 2*smT2**3))*LC(15) &
                 )/(E**(3*y)*shat*(MV2*shat - smT2)**3) - &
              (2*smT*LC(7)*(MV2*(MV2 - 2*shat)*shat**2 + &
                   shat*(-4*MV2 + shat)*smT2 + 2*smT2**2 + &
                   shat*smT*(3*MV2*shat - smT2)*Cosh(y) + &
                   shat*smT*(-(MV2*shat) + smT2)*Sinh(y))*SP(1))/ &
               (E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
              (2*(E**y*shat*smT*(2*MV2*shat - smT2) - &
                   smT2*(-2*MV2*shat + shat**2 + 2*smT2))*LC(6)*SP(2))/ &
               (E**y*shat**2*(MV2*shat*smT - smT**3)) + &
              (LC(13)*(-2*E**(4*y)*shat**2*(MV2*shat - smT2)*smT2* &
                    (MV2**2 + smT2)*SP(3) - &
                   2*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
                    (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2)*SP(3) &
                     - 8*shat**2*smT2**3*SP(4) + &
                   4*E**y*shat*smT**3*smT2*(3*MV2*shat + smT2)*SP(4) + &
                   4*E**(2*y)*smT2**2* &
                    (2*MV2*shat**3 - shat*(3*MV2 + 2*shat)*smT2 + &
                      smT2**2)*SP(4) + &
                   2*E**(3*y)*shat*smT*(MV2*shat - smT2)*smT2* &
                    (MV2*shat*SP(3) - 2*(MV2*shat + smT2)*SP(4))))/ &
               (E**(5*y)*shat**2*smT*(MV2*shat - smT2)**3) + &
              (LC(9)*(-4*E**(4*y)*smT2* &
                    (-2*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                      3*MV2*shat*smT2**2 + smT2**3)*SP(3) + &
                   4*shat**2*smT2**2*(-2*MV2*shat + smT2)*SP(4) + &
                   4*E**y*shat*smT*smT2* &
                    (MV2**2*shat**2 + 2*MV2*shat*smT2 - smT2**2)*SP(4) &
                    + 2*E**(3*y)*shat*smT* &
                    (2*MV2*shat*smT2*(MV2*shat + 3*smT2)*SP(3) + &
                      (MV2*shat - smT2)* &
                       (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2)* &
                       SP(4)) + &
                   2*E**(2*y)*shat**2*smT2* &
                    (-4*MV2*shat*smT2*SP(3) + &
                      (2*MV2**2*(shat**2 - smT2) - 3*MV2*shat*smT2 + &
                         smT2**2)*SP(4))))/ &
               (E**(3*y)*shat**2*smT*(MV2*shat - smT2)**3) + &
              (4*LC(3)*(-(E**(5*y)*shat*(MV2*shat - smT2)*smT2* &
                      (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
                      SP(3)) + &
                   E**(4*y)*smT* &
                    (2*MV2**3*shat**5 + &
                      4*MV2**2*(MV2 - shat)*shat**3*smT2 + &
                      3*MV2*shat**2*(-2*MV2 + shat)*smT2**2 + &
                      (6*MV2 - shat)*shat*smT2**3 - 2*smT2**4)*SP(3) - &
                   2*MV2*shat**3*smT**5*SP(4) + &
                   4*E**y*MV2**2*shat**3*smT2**2*SP(4) - &
                   E**(3*y)*shat*smT2* &
                    (2*MV2**2*shat**2*(MV2*shat + smT2)*SP(3) + &
                      (MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2*SP(4)) &
                    + E**(2*y)*smT*smT2* &
                    (2*MV2**2*shat**4*SP(3) + &
                      (2*MV2**2*shat**4 - &
                         3*MV2*shat**2*(2*MV2 + shat)*smT2 + &
                         shat*(6*MV2 + shat)*smT2**2 - 2*smT2**3)*SP(4)) &
                   ))/(E**(4*y)*shat**3*smT*(MV2*shat - smT2)**3) - &
              (2*(E**y*MV2*shat*(MV2 + shat)*smT + &
                   smT2*(-2*MV2*shat + smT2) + &
                   E**(2*y)*(2*MV2**2*shat**2 - 4*MV2*shat*smT2 + &
                      smT2**2))*LC(5)*SP(5))/ &
               (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
              (2*(E**y*MV2*shat*(5*MV2*shat - smT2)*smT2 + &
                   2*shat*smT**3*(-2*MV2*shat + smT2) + &
                   E**(2*y)*shat*smT* &
                    (-(MV2**2*(MV2 - 5*shat)*shat) - &
                      MV2*(MV2 + 6*shat)*smT2 + smT2**2) + &
                   E**(3*y)*(-2*MV2**3*shat**3 + MV2**2*shat**2*smT2 + &
                      smT2**3))*LC(4)*SP(6))/ &
               (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
              (4*(E**y*MV2 - smT)*smT2* &
                 (-2*shat*smT + E**y*(MV2*shat + smT2))*LC(10)*SP(6))/ &
               (E**(2*y)*(-(MV2*shat) + smT2)**3) + &
              (4*(-shat + E**y*smT)*smT2*LC(14)* &
                 (-3*MV2*shat*smT + smT**3 - &
                   2*shat*((MV2*shat - 2*smT2)*Cosh(y) + &
                      MV2*shat*Sinh(y)))*SP(6))/ &
               (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(2*shat*smT**5 + &
                   E**(2*y)*shat*smT* &
                    (2*MV2**3*shat + MV2*shat*smT2 - smT2**2) - &
                   E**y*smT2*(3*MV2**2*shat**2 + smT2**2) + &
                   E**(3*y)*(MV2*shat - smT2)* &
                    (4*MV2**2*shat**2 - 7*MV2*shat*smT2 + 2*smT2**2))* &
                 LC(2)*SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3) - &
              (4*(E**y*MV2 - smT)*smT* &
                 (-2*MV2*shat**2 + E**y*smT*(MV2*shat + smT2))*LC(8)* &
                 SP(7))/(E**(2*y)*(-(MV2*shat) + smT2)**3) + &
              (4*smT**3*(shat - E**y*smT)* &
                 (2*shat*smT + E**y*(-3*MV2*shat + smT2))*LC(12)*SP(7))/ &
               (E**(4*y)*shat*(MV2*shat - smT2)**3)) + &
           (4*kapPr*(E**(2*y)*shat*smT**3*(MV2*shat - smT2)*SP(2)* &
                 SP(6) + shat*smT*smT2*(SP(3) + SP(4))* &
                 (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)) + &
                E**y*smT2*(-((MV2*shat - smT2)* &
                      (smT2*(SP(3) + SP(4))*SP(5) + shat**2*SP(2)*SP(6)) &
                      ) + 2*shat*(MV2*shat*SP(3) + smT2*SP(4))*SP(6)* &
                    SP(7))))/(E**(2*y)*shat*smT*(-(MV2*shat) + smT2)**2) &
           ) + SI(8)*(kap*((smT* &
                 (12*MV2*shat**2*smT2 - &
                   E**(2*y)*shat*(-(MV2*shat) + smT2)* &
                    (MV2**2 - 5*MV2*shat + 5*smT2) - &
                   E**(4*y)*shat*smT2*(5*MV2*shat + 7*smT2) + &
                   2*E**(3*y)*smT* &
                    (4*MV2**2*shat**2 + 7*MV2*shat*smT2 + smT2**2) + &
                   E**y*smT*(-11*MV2**2*shat**2 - 16*MV2*shat*smT2 + &
                      3*smT2**2))*LC(1))/(E**y*(-(MV2*shat) + smT2)**3) &
               + (smT*(6*shat*smT2**2 + &
                   E**(4*y)*shat*(MV2*shat - 3*smT2)* &
                    (2*MV2*shat + smT2) + &
                   E**(2*y)*shat*(-(MV2*shat) + smT2)* &
                    (MV2*(MV2 + shat) + 2*smT2) - &
                   2*E**y*smT* &
                    (-(MV2**2*shat**2) + 6*MV2*shat*smT2 + smT2**2) - &
                   E**(3*y)*smT* &
                    (MV2**2*shat**2 - 15*MV2*shat*smT2 + 2*smT2**2))* &
                 LC(11))/(E**(3*y)*(-(MV2*shat) + smT2)**3) + &
              (shat*smT*(E**(2*y)*shat*smT*(MV2*shat - smT2) + &
                   2*shat*smT*(2*MV2*shat + 3*smT2) - &
                   2*E**(4*y)*shat*smT*(MV2*shat + 4*smT2) + &
                   E**(3*y)*(3*MV2*shat + smT2)*(MV2*shat + 4*smT2) - &
                   E**y*(5*MV2**2*shat**2 + 12*MV2*shat*smT2 + &
                      3*smT2**2))*LC(15))/ &
               (E**(2*y)*(MV2*shat - smT2)**3) + &
              (2*(MV2*shat*smT - E**(3*y)*shat*smT2 - &
                   E**y*shat*(2*MV2**2 + smT2) + &
                   E**(2*y)*smT*(2*MV2*shat + smT2))*LC(7)*SP(1))/ &
               (E**y*(-(MV2*shat) + smT2)**2) + &
              (4*smT*LC(6)*Sinh(y)*SP(2))/(-(MV2*shat) + smT2) + &
              (LC(9)*(-8*E**(6*y)*shat*smT2**2*SP(3) + &
                   4*E**(5*y)*smT**3*(3*MV2*shat + smT2)*SP(3) + &
                   4*E**(4*y)*shat*smT2* &
                    (-2*smT2*SP(3) + MV2*shat*(2*SP(3) - SP(4))) + &
                   4*shat*(2*MV2*shat - smT2)*smT2*SP(4) + &
                   2*E**y*smT* &
                    (-4*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) &
                    - 2*E**(2*y)*shat* &
                    (-4*MV2*shat*smT2*SP(3) + &
                      (MV2*shat - smT2)*(MV2*(MV2 + shat) - smT2)*SP(4)) &
                     + 4*E**(3*y)*smT* &
                    (-3*MV2*shat*smT2*SP(3) + smT2**2*SP(3) + &
                      2*MV2**2*shat**2*(-SP(3) + SP(4)))))/ &
               (E**(2*y)*(MV2*shat - smT2)**3) + &
              (4*LC(3)*(E**(6*y)*smT**3*(3*MV2*shat - smT2)*SP(3) - &
                   4*E**(5*y)*MV2**2*shat*smT2*SP(3) + &
                   2*MV2*shat*smT**3*SP(4) - &
                   4*E**y*MV2**2*shat*smT2*SP(4) + &
                   4*E**(3*y)*MV2**3*shat**2*(SP(3) + SP(4)) - &
                   E**(2*y)*smT* &
                    (-3*MV2*shat*smT2*SP(4) + smT2**2*SP(4) + &
                      2*MV2**2*shat**2*(SP(3) + SP(4))) + &
                   E**(4*y)*(-2*smT**5*SP(3) + &
                      MV2*shat*smT**3*(7*SP(3) + SP(4)) - &
                      MV2**2*shat**2*smT*(5*SP(3) + 3*SP(4)))))/ &
               (E**(3*y)*(MV2*shat - smT2)**3) + &
              (LC(13)*(-2*E**(5*y)*MV2*shat*smT*(MV2*shat - smT2)* &
                    SP(3) - 8*E**(2*y)*shat*(MV2*shat - smT2)*smT2* &
                    SP(4) + 8*shat*smT2**2*SP(4) - &
                   4*E**y*smT**3*(3*MV2*shat + smT2)*SP(4) + &
                   2*E**(4*y)*MV2*shat* &
                    ((MV2 + shat)*(MV2*shat - smT2)*SP(3) - &
                      4*shat*smT2*SP(4)) - &
                   2*E**(3*y)*smT* &
                    (MV2**2*shat**2*(SP(3) - 4*SP(4)) + &
                      2*smT2**2*SP(4) - MV2*shat*smT2*(SP(3) + 6*SP(4))) &
                   ))/(E**(4*y)*(MV2*shat - smT2)**3) + &
              (2*(-2*MV2*shat*smT + smT**3 + E**(3*y)*shat*smT2 + &
                   E**y*shat*(2*MV2**2 + smT2) - &
                   E**(2*y)*smT*(MV2*shat + 2*smT2))*LC(5)*SP(5))/ &
               (E**y*(-(MV2*shat) + smT2)**2) - &
              (2*(E**(3*y)*(-10*MV2**2*shat**2*smT + &
                      5*MV2*shat*smT**3 + smT**5) + &
                   2*E**y*MV2*shat*smT*(3*MV2*shat - smT2) + &
                   2*E**(4*y)*shat*(3*MV2*shat - 2*smT2)*smT2 + &
                   2*shat*smT2*(-2*MV2*shat + smT2) + &
                   E**(2*y)*shat*(MV2*shat - smT2)*(2*MV2*shat + smT2))* &
                 LC(4)*SP(6))/(E**(2*y)*(MV2*shat - smT2)**3) - &
              (4*shat*(2*E**(4*y)*MV2*shat**2*smT - 2*shat*smT**3 + &
                   2*E**(2*y)*shat*smT*(MV2*shat - smT2) + &
                   E**y*smT2*(3*MV2*shat + smT2) + &
                   E**(3*y)*(-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
                      smT2**2))*LC(14)*SP(6))/ &
               (E**(3*y)*(MV2*shat - smT2)**3) + &
              (8*E**y*shat*smT2* &
                 (3*MV2*shat + smT2 - 4*shat*smT*Cosh(y))*LC(10)* &
                 Sinh(y)*SP(6))/(MV2*shat - smT2)**3 + &
              (2*(2*E**(2*y)*MV2*shat*(-MV2 + shat)*(MV2*shat - smT2) - &
                   2*shat*smT2**2 + &
                   E**(4*y)*shat*smT2*(MV2*shat + smT2) - &
                   E**(3*y)*smT*(MV2*shat + smT2)**2 + &
                   E**y*smT**3*(3*MV2*shat + smT2))*LC(2)*SP(7))/ &
               (E**(2*y)*(MV2*shat - smT2)**3) - &
              (4*shat*(2*MV2*shat**2*smT - 2*E**(4*y)*shat*smT**3 + &
                   2*E**(2*y)*shat*smT*(MV2*shat - smT2) + &
                   E**(3*y)*smT2*(3*MV2*shat + smT2) + &
                   E**y*(-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2)) &
                  *LC(8)*SP(7))/(E**y*(MV2*shat - smT2)**3) - &
              (8*shat*smT2*(3*MV2*shat + smT2 - 4*shat*smT*Cosh(y))* &
                 LC(12)*Sinh(y)*SP(7))/(E**y*(MV2*shat - smT2)**3)) + &
           (4*kapPr*shat*(smT*(-(MV2*shat) + smT2)*(SP(3) + SP(4))* &
                 SP(5) + smT* &
                 ((-(MV2*shat) + smT2)*SP(1) + &
                   2*shat*(SP(3) + SP(4))*SP(6))*SP(7) + &
                E**(2*y)*smT* &
                 (-((MV2*shat - smT2)* &
                      ((SP(3) + SP(4))*SP(5) + SP(2)*SP(6))) + &
                   2*shat*(SP(3) + SP(4))*SP(6)*SP(7)) + &
                E**y*((MV2*shat - smT2)* &
                    (2*MV2*(SP(3) + SP(4))*SP(5) + shat*SP(2)*SP(6)) + &
                   (shat*(MV2*shat - smT2)*SP(1) - &
                      2*(MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7))) &
              )/(E**y*(-(MV2*shat) + smT2)**2)) + &
        SI(15)*(kap*((smT*(-(E**(3*y)*shat*smT2*(11*MV2*shat + smT2)* &
                      (MV2*shat + 2*smT2)) + &
                   E**(4*y)*shat**2*smT**3*(5*MV2*shat + 7*smT2) + &
                   E**(2*y)*shat*smT* &
                    (MV2**2*shat**2*(-14*MT2 + 7*MV2 + 2*shat) + &
                      MV2*(-8*MT2 + 25*MV2 - 6*shat)*shat*smT2 + &
                      2*(11*MT2 + 2*(MV2 + shat))*smT2**2) - &
                   E**y*(MV2**3*shat**3*(-4*MT2 + MV2 + 2*shat) + &
                      MV2**2*(-26*MT2 + 9*MV2 - 6*shat)*shat**2*smT2 + &
                      2*MV2*shat*(14*MT2 + MV2 + shat)*smT2**2 + &
                      2*(MT2 + shat)*smT2**3) + &
                   2*shat*smT*(-(MV2*shat) + smT2)* &
                    (MV2*smT2 + 2*MT2*(-(MV2*shat) + smT2)))*LC(1))/ &
               (shat*(MV2*shat - smT2)**3) + &
              ((-(E**(4*y)*shat**2*(MV2*shat - 3*smT2)*smT2* &
                      (2*MV2*shat + smT2)) + &
                   shat*(MV2*shat - smT2)*smT2* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
                   2*E**(3*y)*MV2*shat**2*smT* &
                    (2*MV2**2*shat**2 - 4*MV2*shat*smT2 - 7*smT2**2) - &
                   2*E**y*smT* &
                    (MV2**3*shat**4 + &
                      2*MV2**2*(-3*MT2 + MV2 - shat)*shat**2*smT2 + &
                      MV2*shat*(3*MT2 + MV2 + 2*shat)*smT2**2 + &
                      (3*MT2 - shat)*smT2**3) + &
                   E**(2*y)*shat* &
                    (2*(4*MT2 - MV2)*MV2**3*shat**3 + &
                      3*MV2**2*(-8*MT2 + MV2)*shat**2*smT2 + &
                      MV2*shat*(6*MT2 + 15*MV2 + shat)*smT2**2 + &
                      (10*MT2 + 2*MV2 - shat)*smT2**3))*LC(11))/ &
               (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              ((-2*E**(4*y)*shat**2*smT**3*(MV2*shat + 4*smT2) + &
                   E**(3*y)*shat*smT2*(5*MV2*shat + smT2)* &
                    (MV2*shat + 4*smT2) + &
                   shat*smT*(MV2*shat - smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
                   E**(2*y)*shat*smT* &
                    (2*(3*MT2 - 2*MV2)*MV2**2*shat**2 + &
                      MV2*shat*(18*MT2 - 18*MV2 + shat)*smT2 - &
                      (24*MT2 + 8*MV2 + shat)*smT2**2) + &
                   E**y*(MV2**3*(-4*MT2 + MV2)*shat**3 + &
                      MV2**2*(-14*MT2 + 5*MV2 - 2*shat)*shat**2*smT2 + &
                      2*MV2*shat*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
                      8*MT2*smT2**3))*LC(15))/ &
               (E**y*(MV2*shat - smT2)**3) + &
              (2*(shat**2*(MV2 - E**y*smT)**3 + &
                   2*MT2*(-2*MV2*shat + smT*(E**y*shat + smT))* &
                    (MV2*shat - smT2))*LC(7)*SP(1))/ &
               (shat*(-(MV2*shat) + smT2)**2) + &
              (2*(MV2 - E**y*smT)*(shat + E**y*smT)*LC(6)*SP(2))/ &
               (MV2*shat - smT2) + &
              (LC(13)*(-2*E**(4*y)*shat*(-(MV2*shat*smT) + smT**3)**2* &
                    SP(3) - 4*E**y*shat*smT* &
                    (MV2**2*(-4*MT2 + MV2)*shat**2 + &
                      2*MV2*(MT2 + 2*MV2)*shat*smT2 + &
                      (2*MT2 + MV2)*smT2**2)*SP(4) + &
                   4*smT2*(MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(4) &
                     + 2*E**(3*y)*shat*smT* &
                    ((MV2 + shat)*(-(MV2*shat) + smT2)**2*SP(3) - &
                      4*MV2*shat**2*smT2*SP(4)) + &
                   2*E**(2*y)*shat* &
                    ((MV2*shat - 2*smT2)*(-(MV2*shat) + smT2)**2* &
                       SP(3) + 6*MV2*shat*smT2*(MV2*shat + smT2)*SP(4))) &
                 )/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (2*LC(9)*(-4*E**(6*y)*shat**2*smT**5*SP(3) + &
                   2*E**(5*y)*shat*smT2**2*(5*MV2*shat + smT2)*SP(3) + &
                   shat*smT*(MV2*shat - smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2)*SP(4) + &
                   E**y*MV2*shat* &
                    (shat*(MV2**2 - MV2*shat + smT2)* &
                       (MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(4) &
                     + 2*E**(4*y)*shat*smT**3* &
                    (2*(MV2*shat*(3*MT2 - 2*MV2 + shat) - &
                         (3*MT2 + MV2 + shat)*smT2)*SP(3) - &
                      MV2*shat**2*SP(4)) + &
                   E**(3*y)*smT2* &
                    (2*(MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
                         MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
                         (2*MT2 + shat)*smT2**2)*SP(3) + &
                      shat*(6*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)* &
                       SP(4)) - &
                   E**(2*y)*shat*smT* &
                    (-2*(MV2*shat - smT2)* &
                       (4*MT2*(-(MV2*shat) + smT2) + &
                         MV2*(MV2*shat + smT2))*SP(3) + &
                      MV2*(-(MV2*shat**2*(6*MT2 - 5*MV2 + shat)) + &
                         shat*(6*MT2 + shat)*smT2 + smT2**2)*SP(4))))/ &
               (E**y*shat*(MV2*shat - smT2)**3) + &
              (4*LC(3)*(-2*E**(4*y)*MV2*shat*smT**3*(4*MV2*shat - smT2)* &
                    SP(3) + E**(5*y)*shat*(3*MV2*shat - smT2)*smT2**2* &
                    SP(3) + MV2*smT* &
                    (MV2**2*shat**2*(-8*MT2 + 2*MV2 + shat) - &
                      2*MV2*shat*(-5*MT2 + shat)*smT2 + &
                      (-2*MT2 + shat)*smT2**2)*SP(4) + &
                   E**(2*y)*MV2*smT* &
                    ((MV2**2*shat**2*(8*MT2 - 2*MV2 + 5*shat) - &
                         2*MV2*shat*(5*MT2 + 3*shat)*smT2 + &
                         (2*MT2 + shat)*smT2**2)*SP(3) + &
                      6*MV2**2*shat**3*SP(4)) - &
                   E**(3*y)*smT2* &
                    ((MV2**2*shat**2*(10*MT2 - 7*MV2 + 3*shat) + &
                         MV2*(-14*MT2 + MV2 - 4*shat)*shat*smT2 + &
                         (4*MT2 + shat)*smT2**2)*SP(3) + &
                      MV2*shat**2*(3*MV2*shat - smT2)*SP(4)) + &
                   E**y*(-2*(MV2*shat - smT2)* &
                       (MV2**2*(-4*MT2 + MV2)*shat**2 + &
                         6*MT2*MV2*shat*smT2 - 2*MT2*smT2**2)*SP(3) + &
                      (MV2**3*(12*MT2 - 3*MV2 - 2*shat)*shat**3 + &
                         MV2**2*shat**2*(-22*MT2 - 3*MV2 + 5*shat)* &
                          smT2 + 2*MV2*(7*MT2 - 2*shat)*shat*smT2**2 + &
                         (-4*MT2 + shat)*smT2**3)*SP(4))))/ &
               (E**y*shat*(MV2*shat - smT2)**3) + &
              (2*LC(5)*(-2*MT2*(-2*MV2*shat + smT*(E**y*shat + smT))* &
                    (MV2*shat - smT2) + &
                   shat*(-(MV2**3*shat) - &
                      MV2**2*shat*(shat - 2*E**y*smT) + &
                      smT**3* &
                       (-2*smT + E**y*(shat + E**(2*y)*shat - E**y*smT)) &
                        + E**y*MV2*smT* &
                       (-shat**2 + smT2 + &
                         shat*smT*(Cosh(y) - 5*Sinh(y)))))*SP(5))/ &
               (shat*(-(MV2*shat) + smT2)**2) + &
              (2*(-2*E**(4*y)*shat*smT**3*(3*MV2*shat - 2*smT2) + &
                   smT*(MV2*shat - smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) - &
                   E**(3*y)*smT2* &
                    (-12*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2) + &
                   E**(2*y)*smT* &
                    (MV2**2*shat**2*(14*MT2 - 7*MV2 + shat) - &
                      MV2*shat*(22*MT2 + shat)*smT2 + &
                      (8*MT2 + MV2)*smT2**2) + &
                   E**y*MV2*(shat*(MV2**2 - MV2*shat + smT2)* &
                       (MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
                 LC(4)*SP(6))/(E**y*(MV2*shat - smT2)**3) + &
              (4*E**y*smT*(-(shat*(MV2 - E**y*smT)**2* &
                      (MV2*shat + smT*(-2*E**y*shat + smT))) + &
                   2*MT2*(2*MV2*shat + smT*(-3*E**y*shat + smT))* &
                    (MV2*shat - smT2))*LC(10)*SP(6))/ &
               (-(MV2*shat) + smT2)**3 - &
              (4*(E**y*shat - smT)* &
                 (MV2*shat*(MV2 - E**y*smT)* &
                    (MV2*shat + smT*(-2*E**y*shat + smT)) - &
                   2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*LC(14)* &
                 SP(6))/(E**y*(MV2*shat - smT2)**3) + &
              (2*(E**(4*y)*shat**2*smT**3*(MV2*shat + smT2) - &
                   2*E**(3*y)*MV2*shat**2*smT2*(MV2*shat + 2*smT2) + &
                   3*shat*smT*(MV2*shat - smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
                   E**(2*y)*shat*smT* &
                    (MV2**2*shat**2*(-4*MT2 + MV2 + 2*shat) + &
                      MV2*(2*MT2 + 5*MV2 - shat)*shat*smT2 + &
                      (2*MT2 - shat)*smT2**2) - &
                   2*E**y*(MV2**3*shat**4 + &
                      MV2**2*shat**2*(-4*MT2 + MV2 + shat)*smT2 + &
                      MV2*(5*MT2 - 2*shat)*shat*smT2**2 - MT2*smT2**3))* &
                 LC(2)*SP(7))/(E**y*shat*(MV2*shat - smT2)**3) + &
              (4*(shat*(-MV2 + E**y*smT)* &
                    (MV2*shat + smT*(-2*E**y*shat + smT))* &
                    (MV2*(shat + E**y*smT) - (1 + E**(2*y))*smT2) + &
                   2*MT2*(MV2*shat - smT2)* &
                    (2*MV2*shat*(shat + E**y*smT) + &
                      (-((2 + 3*E**(2*y))*shat) + E**y*smT)*smT2))* &
                 LC(8)*SP(7))/(MV2*shat - smT2)**3 + &
              (4*smT*(shat*(MV2 - E**y*smT)**2* &
                    (MV2*shat + smT*(-2*E**y*shat + smT)) - &
                   2*MT2*(2*MV2*shat + smT*(-3*E**y*shat + smT))* &
                    (MV2*shat - smT2))*LC(12)*SP(7))/ &
               (E**y*(-(MV2*shat) + smT2)**3)) - &
           (2*kapPr*(smT*(-(MV2*shat) + smT2)* &
                 ((MV2*shat - smT2)*SP(4)*SP(5) + MV2*shat*SP(1)*SP(7)) &
                 + E**(3*y)*shat*smT2* &
                 ((MV2*shat - smT2)* &
                    (2*(SP(3) + SP(4))*SP(5) + SP(2)*SP(6)) - &
                   4*shat*(SP(3) + SP(4))*SP(6)*SP(7)) + &
                E**(2*y)*smT* &
                 (-((MV2*shat - smT2)* &
                      (((3*MV2*shat + smT2)*SP(3) + 4*MV2*shat*SP(4))* &
                         SP(5) + shat*(MV2 + shat)*SP(2)*SP(6))) + &
                   shat*(shat*(-(MV2*shat) + smT2)*SP(1) + &
                      2*(3*MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7) &
                   ) - E**y*((MV2*shat - smT2)* &
                    (8*MT2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) - &
                      2*shat* &
                       (MV2**2*SP(3) + (MV2*(MV2 + shat) - smT2)*SP(4))* &
                       SP(5) - shat*smT2*SP(2)*SP(6)) + &
                   2*shat*(MV2*shat*(-(MV2*shat) + smT2)*SP(1) + &
                      (4*MT2*(-(MV2*shat) + smT2) + &
                         MV2*(MV2*shat + smT2))*(SP(3) + SP(4))*SP(6))* &
                    SP(7))))/(E**y*(-(MV2*shat) + smT2)**2)) + &
        SI(14)*(kap*((smT*(-12*MV2*shat**3*smT**3 + &
                   2*E**(4*y)*shat*smT*(MV2*shat - smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
                   E**y*shat*smT2* &
                    (23*MV2**2*shat**2 + 14*MV2*shat*smT2 - smT2**2) + &
                   E**(2*y)*shat*smT* &
                    (MV2**2*shat**2*(28*MT2 - 13*MV2 + 3*shat) - &
                      2*MV2*shat*(10*MT2 + 11*MV2 + 2*shat)*smT2 + &
                      (-8*MT2 - MV2 + shat)*smT2**2) + &
                   E**(3*y)*(-(MV2**3*shat**3*(8*MT2 - 2*MV2 + shat)) - &
                      4*MV2**2*shat**2*(4*MT2 - 2*MV2 + shat)*smT2 + &
                      MV2*shat*(20*MT2 + 2*MV2 + 7*shat)*smT2**2 + &
                      2*(2*MT2 - shat)*smT2**3))*LC(1))/ &
               (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              ((-6*shat**2*smT2**3 - &
                   E**y*shat*smT**3* &
                    (MV2**2*shat**2 - 17*MV2*shat*smT2 - 2*smT2**2) + &
                   2*E**(4*y)*shat*(MV2*shat - smT2)* &
                    (MV2**2*(-4*MT2 + MV2)*shat**2 + &
                      (6*MT2 - MV2)*MV2*shat*smT2 + &
                      (-2*MT2 + MV2)*smT2**2) + &
                   E**(2*y)*shat*smT2* &
                    (MV2**2*shat**2*(-2*MT2 + 2*MV2 + shat) + &
                      2*(11*MT2 - 8*MV2)*MV2*shat*smT2 - &
                      (20*MT2 + 4*MV2 + shat)*smT2**2) + &
                   E**(3*y)*smT* &
                    (-(MV2**3*shat**3*(-4*MT2 + MV2 + 3*shat)) + &
                      MV2**2*shat**2*(-22*MT2 + 5*MV2 + 4*shat)*smT2 + &
                      MV2*(14*MT2 + 2*MV2 - 3*shat)*shat*smT2**2 + &
                      2*(2*MT2 + shat)*smT2**3))*LC(11))/ &
               (E**(4*y)*shat*(MV2*shat - smT2)**3) + &
              ((2*shat**2*smT**3*(2*MV2*shat + 3*smT2) - &
                   E**y*shat*smT2* &
                    (8*MV2**2*shat**2 + 19*MV2*shat*smT2 + 3*smT2**2) + &
                   E**(2*y)*shat*smT* &
                    (5*MV2**2*(-2*MT2 + MV2)*shat**2 + &
                      2*MV2*(-5*MT2 + 9*MV2)*shat*smT2 + &
                      (20*MT2 + 7*MV2)*smT2**2) - &
                   E**(3*y)*(MV2*shat + 4*smT2)* &
                    (MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
                 LC(15))/(E**(3*y)*(MV2*shat - smT2)**3) + &
              (2*(MV2*shat**2*smT2 - &
                   E**y*shat**2*smT*(2*MV2**2 + smT2) + &
                   E**(3*y)*shat*smT* &
                    (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2) + &
                   E**(2*y)*(MV2**2*(-4*MT2 + MV2)*shat**2 + &
                      MV2*shat*(6*MT2 + shat)*smT2 + &
                      (-2*MT2 + shat)*smT2**2))*LC(7)*SP(1))/ &
               (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
              (2*(E**y*MV2 - smT)*(E**y*shat + smT)*LC(6)*SP(2))/ &
               (E**(2*y)*(MV2*shat - smT2)) + &
              (4*LC(3)*(E**(6*y)*smT2*(-(MV2*shat) + smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2)*SP(3) + &
                   2*E**(5*y)*MV2*smT* &
                    (MV2**2*shat**2*(-4*MT2 + MV2 + shat) + &
                      MV2*(5*MT2 - shat)*shat*smT2 - MT2*smT2**2)*SP(3) &
                    - 6*E**y*MV2**2*shat**2*smT**3*SP(4) + &
                   2*MV2*shat**2*smT2**2*SP(4) + &
                   E**(2*y)*smT2* &
                    (-2*MV2**2*shat**3*SP(3) + &
                      (-(MV2**2*shat**2*(8*MT2 - 6*MV2 + shat)) + &
                         MV2*shat*(10*MT2 + shat)*smT2 - 2*MT2*smT2**2)* &
                       SP(4)) + &
                   2*E**(3*y)*MV2*smT* &
                    (MV2*shat**2*(2*MV2*shat + smT2)*SP(3) + &
                      (MV2**2*shat**2*(4*MT2 - MV2 + shat) - &
                         MV2*shat*(5*MT2 + shat)*smT2 + MT2*smT2**2)* &
                       SP(4)) + &
                   E**(4*y)*((-2*MV2**3*shat**3*(-4*MT2 + MV2 + shat) - &
                         4*MV2**2*(3*MT2 + MV2 - shat)*shat**2*smT2 - &
                         3*MV2*shat*(-2*MT2 + shat)*smT2**2 + &
                         (-2*MT2 + shat)*smT2**3)*SP(3) + &
                      (MV2*shat - smT2)* &
                       ((4*MT2 - MV2)*MV2**2*shat**2 - &
                         6*MT2*MV2*shat*smT2 + 2*MT2*smT2**2)*SP(4))))/ &
               (E**(4*y)*shat*(MV2*shat - smT2)**3) + &
              (2*smT*LC(9)*(-2*E**(5*y)*smT* &
                    (MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(3) &
                     - 2*shat**2*smT2*(-2*MV2*shat + smT2)*SP(4) + &
                   E**y*shat*smT*(-7*MV2**2*shat**2 + smT2**2)*SP(4) + &
                   E**(4*y)*shat* &
                    (2*(MV2**2*(-4*MT2 + MV2)*shat**2 + &
                         2*MV2*(MT2 + 2*MV2)*shat*smT2 + &
                         (2*MT2 + MV2)*smT2**2)*SP(3) + &
                      (MV2*shat - smT2)* &
                       (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2)*SP(4)) + &
                   E**(3*y)*smT* &
                    (-6*MV2*shat**2*(MV2*shat + smT2)*SP(3) + &
                      (MV2**2*shat**2*(8*MT2 - 2*MV2 + 3*shat) - &
                         2*MV2*shat*(5*MT2 + 2*shat)*smT2 + &
                         (2*MT2 + shat)*smT2**2)*SP(4)) - &
                   E**(2*y)*shat* &
                    (-4*MV2*shat**2*smT2*SP(3) + &
                      (MV2**2*shat**2*(10*MT2 - 3*MV2 + 2*shat) - &
                         MV2*shat*(14*MT2 + 4*MV2 + 3*shat)*smT2 + &
                         (4*MT2 + MV2 + shat)*smT2**2)*SP(4))))/ &
               (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
              (LC(13)*(-2*E**(6*y)*shat*smT*(MV2*shat - smT2)* &
                    (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2)*SP(3) + &
                   2*E**(5*y)*(-(MV2*shat) + smT2)* &
                    (MV2**2*(-4*MT2 + MV2)*shat**2 + &
                      MV2*shat*(6*MT2 + shat)*smT2 + &
                      (-2*MT2 + shat)*smT2**2)*SP(3) + &
                   8*shat**2*smT**5*SP(4) - &
                   4*E**y*shat*smT2**2*(5*MV2*shat + smT2)*SP(4) - &
                   8*E**(2*y)*shat*smT**3* &
                    (MV2*shat*(3*MT2 - 2*MV2 + shat) - &
                      (3*MT2 + MV2 + shat)*smT2)*SP(4) - &
                   2*E**(3*y)*smT2* &
                    (MV2*shat**2*(MV2*shat - smT2)*SP(3) + &
                      2*(MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
                         MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
                         (2*MT2 + shat)*smT2**2)*SP(4)) + &
                   2*E**(4*y)*shat*smT*(MV2*shat - smT2)* &
                    (shat*(2*MV2**2 + smT2)*SP(3) - &
                      2*(4*MT2*(-(MV2*shat) + smT2) + &
                         MV2*(MV2*shat + smT2))*SP(4))))/ &
               (E**(5*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(E**y*MV2*shat*smT*(3*MV2*shat + shat**2 - smT2) + &
                   shat*smT2*(-2*MV2*shat + smT2) + &
                   E**(2*y)*(MV2**2*shat**2*(4*MT2 - MV2 + shat) - &
                      2*MV2*shat*(3*MT2 + 2*shat)*smT2 + &
                      (2*MT2 + shat)*smT2**2) + &
                   E**(3*y)*shat*smT* &
                    (MV2*smT2 + 2*MT2*(-(MV2*shat) + smT2)))*LC(5)*SP(5) &
                 )/(E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
              (2*(E**(2*y)*(MV2**2*(-10*MT2 + 6*MV2 - 3*shat)*shat**2* &
                       smT + 2*MV2*shat*(7*MT2 + shat)*smT**3 + &
                      (-4*MT2 + shat)*smT**5) + &
                   2*shat*smT**3*(2*MV2*shat - smT2) - &
                   3*E**y*MV2*shat*(3*MV2*shat - smT2)*smT2 + &
                   4*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
                    (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) - &
                   E**(3*y)*MV2* &
                    (MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
                      MV2*(2*MT2 + MV2 - 2*shat)*shat*smT2 + &
                      (2*MT2 + 5*shat)*smT2**2))*LC(4)*SP(6))/ &
               (E**(3*y)*(MV2*shat - smT2)**3) + &
              (4*smT*(-2*shat**2*smT**3 + &
                   E**y*shat*smT2*(5*MV2*shat + smT2) - &
                   2*E**(2*y)*shat*smT* &
                    (3*MT2*(-(MV2*shat) + smT2) + &
                      MV2*(2*MV2*shat + smT2)) + &
                   E**(3*y)*(MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
                 LC(10)*SP(6))/(E**(2*y)*(-(MV2*shat) + smT2)**3) + &
              (4*(2*shat**2*smT2**2 - &
                   E**y*shat*smT**3*(5*MV2*shat + smT2) - &
                   2*E**(2*y)*shat*smT2* &
                    (MV2*shat*(3*MT2 - 2*MV2 + shat) - &
                      (3*MT2 + MV2 + shat)*smT2) - &
                   E**(3*y)*smT* &
                    (MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
                      MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
                      (2*MT2 + shat)*smT2**2) - &
                   E**(4*y)*shat*(MV2*shat - smT2)* &
                    (4*MT2*(-(MV2*shat) + smT2) + MV2*(MV2*shat + smT2)) &
                   )*LC(14)*SP(6))/(E**(4*y)*(MV2*shat - smT2)**3) + &
              (2*(-2*shat**2*smT**5 + &
                   E**(4*y)*shat*smT*(MV2*shat - smT2)* &
                    (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2) + &
                   E**y*shat*smT2* &
                    (3*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) - &
                   E**(2*y)*shat*smT* &
                    (MV2**2*(-2*MT2 + 5*MV2)*shat**2 + &
                      MV2*shat*(-2*MT2 + shat)*smT2 + &
                      (4*MT2 + MV2 - shat)*smT2**2) + &
                   E**(3*y)*(2*MV2**3*(-4*MT2 + MV2)*shat**3 + &
                      MV2**2*shat**2*(12*MT2 + shat)*smT2 - &
                      6*MT2*MV2*shat*smT2**2 + (2*MT2 - shat)*smT2**3))* &
                 LC(2)*SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3) + &
              (4*(-shat + E**y*smT)* &
                 (2*MV2*shat**2*smT2 - &
                   E**y*MV2*shat*smT*(3*MV2*shat + smT2) + &
                   E**(2*y)*(MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
                 LC(8)*SP(7))/(E**(2*y)*(MV2*shat - smT2)**3) + &
              (4*smT*(-2*shat**2*smT**3 + &
                   E**y*shat*smT2*(5*MV2*shat + smT2) - &
                   2*E**(2*y)*shat*smT* &
                    (3*MT2*(-(MV2*shat) + smT2) + &
                      MV2*(2*MV2*shat + smT2)) + &
                   E**(3*y)*(MV2**2*shat*(MV2*shat + smT2) - &
                      2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
                 LC(12)*SP(7))/(E**(4*y)*(MV2*shat - smT2)**3)) + &
           (2*kapPr*(E**(3*y)*smT*(MV2*shat - smT2)* &
                 ((MV2*shat - smT2)*SP(3)*SP(5) + MV2*shat*SP(2)*SP(6)) &
                 + shat*smT2* &
                 (-2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) + &
                   ((-(MV2*shat) + smT2)*SP(1) + &
                      4*shat*(SP(3) + SP(4))*SP(6))*SP(7)) + &
                E**y*smT*((MV2*shat - smT2)* &
                    ((4*MV2*shat*SP(3) + (3*MV2*shat + smT2)*SP(4))* &
                       SP(5) + shat**2*SP(2)*SP(6)) + &
                   shat*((MV2 + shat)*(MV2*shat - smT2)*SP(1) - &
                      2*(3*MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7) &
                   ) + E**(2*y)* &
                 (2*(MV2*shat - smT2)* &
                    (4*MT2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) + &
                      shat*(-(((MV2*(MV2 + shat) - smT2)*SP(3) + &
                              MV2**2*SP(4))*SP(5)) - &
                         MV2*shat*SP(2)*SP(6))) + &
                   shat*(smT2*(-(MV2*shat) + smT2)*SP(1) + &
                      2*(4*MT2*(-(MV2*shat) + smT2) + &
                         MV2*(MV2*shat + smT2))*(SP(3) + SP(4))*SP(6))* &
                    SP(7))))/(E**(2*y)*(-(MV2*shat) + smT2)**2)) + &
        SI(11)*(kap*(((12*MV2*shat**4*smT**3 - &
                   2*E**(4*y)*shat*smT**5*(MV2*shat - smT2) - &
                   E**(2*y)*shat*smT* &
                    (MV2**2*shat**3*(8*MT2 + 10*MV2 + shat) + &
                      (8*MT2 - 57*MV2)*MV2*shat**2*smT2 - &
                      shat*(16*MT2 - 24*MV2 + shat)*smT2**2 - 13*smT2**3 &
                      ) + E**(3*y)*smT2* &
                    (MV2**2*shat**3*(8*MT2 + 4*MV2 + shat) + &
                      2*MV2*shat**2*(4*MT2 - 14*MV2 + shat)*smT2 + &
                      (-16*MT2 + 20*MV2 - 3*shat)*shat*smT2**2 - &
                      8*smT2**3) + &
                   E**y*shat**2* &
                    (4*MV2**3*shat**3 - 23*MV2**2*shat**2*smT2 - &
                      14*MV2*shat*smT2**2 - 3*smT2**3))*LC(1))/ &
               (E**y*shat**2*(MV2*shat - smT2)**3) + &
              ((-6*shat**4*smT**5 + &
                   2*E**(5*y)*shat*(MV2*shat - 2*smT2)* &
                    (MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2 + &
                   3*E**y*shat**3*smT2**2*(3*MV2*shat + 5*smT2) - &
                   E**(2*y)*shat**2*smT* &
                    (MV2**2*(4*MT2 + 11*MV2 - shat)*shat**3 - &
                      5*MV2*(4*MT2 + 5*MV2)*shat**2*smT2 + &
                      shat*(16*MT2 + 44*MV2 + shat)*smT2**2 + 6*smT2**3) &
                     + E**(3*y)*shat* &
                    (4*MV2**3*(8*MT2 + 2*MV2 - shat)*shat**4 + &
                      8*MV2**2*shat**3*(-11*MT2 + MV2 + shat)*smT2 - &
                      2*MV2*shat**2*(-28*MT2 + 28*MV2 + 5*shat)* &
                       smT2**2 + shat*(83*MV2 + 6*shat)*smT2**3 - &
                      19*smT2**4) - &
                   E**(4*y)*smT* &
                    (16*MV2**3*(2*MT2 + MV2)*shat**4 - &
                      MV2**2*shat**3*(92*MT2 + 5*(7*MV2 + shat))*smT2 + &
                      MV2*shat**2*(76*MT2 + 17*MV2 + 4*shat)*smT2**2 + &
                      shat*(-16*MT2 + 16*MV2 + shat)*smT2**3 - 8*smT2**4 &
                      ))*LC(11))/ &
               (E**(3*y)*shat**2*(-shat + E**y*smT)* &
                 (MV2*shat - smT2)**3) + &
              ((5*E**y*shat**2*smT*(MV2*shat + smT2)* &
                    (MV2*shat + 2*smT2) - &
                   2*shat**3*smT2*(2*MV2*shat + 3*smT2) + &
                   E**(3*y)*smT* &
                    (MV2**2*shat**3*(-4*MT2 + 5*MV2 + 2*shat) - &
                      MV2*shat**2*(12*MT2 + 5*MV2 + 4*shat)*smT2 + &
                      2*shat*(8*MT2 + 6*MV2 + shat)*smT2**2 - 2*smT2**3) &
                     - E**(2*y)*shat* &
                    (MV2**2*shat**3*(-4*MT2 + MV2 + 2*shat) - &
                      4*MV2*shat**2*(3*MT2 - 3*MV2 + shat)*smT2 + &
                      shat*(16*MT2 + 11*MV2 + 2*shat)*smT2**2 + &
                      6*smT2**3))*LC(15))/ &
               (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(-(shat**2*smT**3) - &
                   E**(3*y)*shat*(MV2*shat - 2*smT2)*smT2 + &
                   E**(2*y)*smT* &
                    (MV2*shat**2*(MV2 + 3*shat) - &
                      shat*(4*MV2 + 5*shat)*smT2 + 2*smT2**2) - &
                   E**y*shat* &
                    (MV2*shat**2*(3*MV2 + 2*shat) - &
                      shat*(8*MV2 + 3*shat)*smT2 + 3*smT2**2))*LC(7)* &
                 SP(1))/(E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
              (2*(-2*shat**2*smT + E**(3*y)*shat*smT2 + &
                   2*E**(2*y)*smT*(-(MV2*shat) + smT2) - &
                   E**y*shat*(-6*MV2*shat + shat**2 + 4*smT2))*LC(6)* &
                 SP(2))/(E**y*shat**2*(MV2*shat - smT2)) + &
              (2*(shat - E**y*smT)*LC(9)* &
                 (-2*E**(4*y)*smT2* &
                    (2*MV2*(-2*MT2 + MV2)*shat**2 + &
                      (4*MT2 - MV2)*shat*smT2 + smT2**2)*SP(3) + &
                   4*E**y*MV2**2*shat**3*smT*SP(4) + &
                   2*shat**2*smT2*(-2*MV2*shat + smT2)*SP(4) + &
                   E**(3*y)*shat*smT* &
                    (2*MV2*shat*(MV2*shat + 3*smT2)*SP(3) + &
                      smT2*(-(MV2*shat) + smT2)*SP(4)) + &
                   E**(2*y)*shat**2* &
                    (-4*MV2*shat*smT2*SP(3) + &
                      (2*MV2**2*shat*(2*MT2 + shat) - &
                         MV2*(4*MT2 + 2*MV2 + 3*shat)*smT2 + smT2**2)* &
                       SP(4))))/(E**(2*y)*shat**2*(MV2*shat - smT2)**3) &
               + (LC(13)*(-2*E**(7*y)*shat*smT**5*(MV2*shat - smT2)* &
                    SP(3) + 2*E**(6*y)*shat**2*(MV2*shat - smT2)*smT2* &
                    (MV2**2 + 3*smT2)*SP(3) + 8*shat**4*smT2**2*SP(4) - &
                   4*E**y*shat**3*smT**3*(3*MV2*shat + 5*smT2)*SP(4) - &
                   4*E**(2*y)*shat**2*smT2* &
                    (2*MV2*shat**2*(2*MT2 - MV2 + shat) - &
                      shat*(4*MT2 + 5*MV2 + 2*shat)*smT2 - 5*smT2**2)* &
                    SP(4) - 2*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
                    (shat*(2*MV2**2*shat + MV2*smT2 + 3*shat*smT2)* &
                       SP(3) + &
                      2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
                       SP(4)) - &
                   2*E**(3*y)*shat*smT* &
                    (MV2*shat**3*(MV2*shat - smT2)*SP(3) - &
                      2*(MV2**2*shat**3*(-6*MV2 + shat) + &
                         2*MV2*shat**2*(4*MT2 + 7*MV2 + 2*shat)*smT2 - &
                         shat*(8*MT2 + 19*MV2 + 5*shat)*smT2**2 + &
                         3*smT2**3)*SP(4)) + &
                   2*E**(4*y)* &
                    (shat**3*(MV2*shat - smT2)* &
                       (MV2**2*shat + 2*MV2*smT2 + shat*smT2)*SP(3) + &
                      2*(2*MV2**3*shat**4*(4*MT2 + 2*MV2 + shat) - &
                         2*MV2**2*shat**3*(12*MT2 + 5*MV2 + 4*shat)* &
                          smT2 + &
                         4*MV2*shat**2*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
                         shat*(-4*MT2 + MV2 + 2*shat)*smT2**3 - smT2**4) &
                        *SP(4))))/ &
               (E**(4*y)*shat**2*(-shat + E**y*smT)* &
                 (MV2*shat - smT2)**3) + &
              (4*LC(3)*(E**(7*y)*shat*(MV2*shat - smT2)*smT2**3*SP(3) + &
                   E**(6*y)*smT**3* &
                    (-2*MV2**2*shat**3*(-4*MT2 + 2*MV2 + shat) + &
                      MV2*shat**2*(-12*MT2 + 6*MV2 + shat)*smT2 + &
                      shat*(4*MT2 - 6*MV2 + shat)*smT2**2 + 2*smT2**3)* &
                    SP(3) + 2*MV2*shat**5*smT**3*SP(4) - &
                   4*E**y*MV2*shat**4*smT2*(MV2*shat + smT2)*SP(4) + &
                   E**(5*y)*shat*smT2* &
                    ((2*MV2**2*shat**3*(-8*MT2 + 5*MV2 + 2*shat) + &
                         MV2*shat**2*(24*MT2 - 5*(2*MV2 + shat))*smT2 + &
                         shat*(-8*MT2 + 12*MV2 + shat)*smT2**2 - &
                         4*smT2**3)*SP(3) + &
                      smT2*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
                       SP(4)) - &
                   E**(2*y)*shat**2*smT* &
                    (2*MV2**2*shat**4*SP(3) + &
                      (8*MT2*MV2**2*shat**3 + &
                         MV2*shat**2*(-12*MT2 - 14*MV2 + shat)*smT2 + &
                         (4*(MT2 + MV2) - shat)*shat*smT2**2 - 2*smT2**3 &
                         )*SP(4)) + &
                   E**(3*y)*shat* &
                    (2*MV2**2*shat**4*(MV2*shat + 3*smT2)*SP(3) + &
                      (2*MV2**3*(8*MT2 - shat)*shat**4 + &
                         8*MV2**2*shat**3*(-4*MT2 + shat)*smT2 + &
                         MV2*(24*MT2 - 16*MV2 - 7*shat)*shat**2* &
                          smT2**2 + &
                         shat*(-8*MT2 + 12*MV2 + shat)*smT2**3 - &
                         4*smT2**4)*SP(4)) + &
                   E**(4*y)*smT* &
                    (-(shat**2* &
                         (2*MV2**2*shat**3*(-4*MT2 + 4*MV2 + shat) - &
                           3*MV2*shat**2*(-4*MT2 + shat)*smT2 + &
                           shat*(-4*MT2 + 6*MV2 + shat)*smT2**2 - &
                           2*smT2**3)*SP(3)) + &
                      (2*MV2**3*shat**4*(-4*MT2 + shat) + &
                         2*MV2**2*(8*MT2 - 5*shat)*shat**3*smT2 + &
                         MV2*shat**2*(-12*MT2 + 6*MV2 + 11*shat)* &
                          smT2**2 + &
                         shat*(4*MT2 - 3*(2*MV2 + shat))*smT2**3 + &
                         2*smT2**4)*SP(4))))/ &
               (E**(3*y)*shat**3*(-shat + E**y*smT)* &
                 (MV2*shat - smT2)**3) + &
              (2*(shat*smT*(2*MV2*shat - smT2) - E**(3*y)*smT2**2 + &
                   E**(2*y)*shat*smT*(MV2*(MV2 + shat) + smT2) - &
                   E**y*(MV2*shat**2*(3*MV2 + shat) - 2*MV2*shat*smT2 + &
                      smT2**2))*LC(5)*SP(5))/ &
               (E**y*shat*(-(MV2*shat) + smT2)**2) + &
              (2*(shat - E**y*smT)* &
                 (E**y*MV2*shat*smT*(5*MV2*shat - smT2) + &
                   2*shat*smT2*(-2*MV2*shat + smT2) + &
                   E**(3*y)*smT*(-(MV2*shat) + smT2)* &
                    (3*MV2*shat + smT2) + &
                   E**(2*y)*shat* &
                    (MV2**2*shat*(4*MT2 - MV2 + shat) - &
                      MV2*(4*MT2 + MV2 - 2*shat)*smT2 - 3*smT2**2))* &
                 LC(4)*SP(6))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (4*smT*(shat - E**y*smT)* &
                 (2*shat*smT2 - E**y*smT*(3*MV2*shat + smT2) + &
                   E**(2*y)*(4*MT2*(-(MV2*shat) + smT2) + &
                      MV2*(MV2*shat + smT2)))*LC(10)*SP(6))/ &
               (E**y*(-(MV2*shat) + smT2)**3) + &
              (4*(-2*shat**3*smT**3 + &
                   2*E**(4*y)*MV2*shat**2*smT*(MV2*shat - smT2) + &
                   3*E**y*shat**2*smT2*(MV2*shat + smT2) - &
                   2*E**(2*y)*shat*smT* &
                    (MV2*(-2*MT2 + MV2 - shat)*shat**2 + &
                      shat*(2*MT2 + MV2 + shat)*smT2 + smT2**2) + &
                   E**(3*y)*(2*MV2**2*(MV2 - shat)*shat**3 - &
                      4*MV2*(MT2 + MV2)*shat**2*smT2 + &
                      shat*(4*MT2 + 5*MV2 + 2*shat)*smT2**2 - smT2**3))* &
                 LC(14)*SP(6))/(E**(3*y)*shat*(MV2*shat - smT2)**3) + &
              (2*(E**(4*y)*MV2*shat*(MV2*shat - smT2)*smT2 + &
                   2*shat**2*smT2**2 - &
                   3*E**y*shat*smT*(MV2**2*shat**2 + smT2**2) + &
                   E**(3*y)*shat*smT* &
                    (-(MV2**2*shat*(-4*MT2 + 2*MV2 + shat)) - &
                      4*MT2*MV2*smT2 + smT2**2) - &
                   E**(2*y)*(4*(MT2 - MV2)*MV2**2*shat**3 + &
                      MV2*(-4*MT2 + 3*MV2 - shat)*shat**2*smT2 + &
                      shat*(-6*MV2 + shat)*smT2**2 + smT2**3))*LC(2)* &
                 SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
              (4*smT*(-shat + E**y*smT)* &
                 (2*MV2*shat**2 + E**y*(-5*MV2*shat*smT + smT**3) + &
                   E**(2*y)*(4*MT2*(-(MV2*shat) + smT2) + &
                      MV2*(MV2*shat + smT2)))*LC(8)*SP(7))/ &
               (E**y*(-(MV2*shat) + smT2)**3) + &
              (4*(2*shat**3*smT**3 - &
                   3*E**y*shat**2*smT2*(MV2*shat + smT2) + &
                   2*E**(2*y)*shat*smT* &
                    (MV2*(-2*MT2 + MV2)*shat**2 + &
                      (2*MT2 + MV2)*shat*smT2 + smT2**2) + &
                   E**(3*y)*(-2*MV2**3*shat**3 + &
                      4*MV2*(MT2 + MV2)*shat**2*smT2 - &
                      (4*MT2 + 5*MV2)*shat*smT2**2 + smT2**3))*LC(12)* &
                 SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3)) - &
           (4*kapPr*(shat - E**y*smT)* &
              (shat*smT*(-((MV2*shat - smT2)*(SP(3) + SP(4))*SP(5)) + &
                   ((-(MV2*shat) + smT2)*SP(1) + &
                      2*shat*(SP(3) + SP(4))*SP(6))*SP(7)) + &
                E**y*((MV2*shat - smT2)*smT2*(SP(3) + SP(4))*SP(5) + &
                   shat*(shat*(MV2*shat - smT2)*SP(1) - &
                      2*(smT2*SP(3) + MV2*shat*SP(4))*SP(6))*SP(7))))/ &
            (E**y*shat*(-(MV2*shat) + smT2)**2)) + &
        SI(10)*(kap*((E**y*(4*shat*smT*(MV2*shat - smT2)**3 + &
                   2*shat*smT**3*smT2*(-(MV2*shat) + smT2) - &
                   E**(3*y)*shat**2*smT2**2*(5*MV2*shat + 7*smT2) - &
                   2*E**(2*y)*shat*smT* &
                    ((-(MV2*shat*smT) + smT**3)**2 - &
                      smT2*(3*MV2**2*shat**2 + 8*MV2*shat*smT2 + &
                         smT2**2)) + &
                   E**y*(2*(-(MV2*shat*smT) + smT**3)**2* &
                       (-3*MV2*shat + shat**2 + 4*smT2) - &
                      shat*smT2* &
                       (MV2**2*shat**2*(MV2 + 2*shat) + &
                         3*MV2*(3*MV2 - 2*shat)*shat*smT2 + &
                         2*(MV2 + 2*shat)*smT2**2)))*LC(1))/ &
               (shat**2*(MV2*shat - smT2)**3) + &
              ((E**(3*y)*shat**2*(MV2*shat - 3*smT2)*smT2* &
                    (2*MV2*shat + smT2) - &
                   shat*smT*(MV2*shat - smT2)* &
                    (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2) + &
                   E**(2*y)*shat*smT* &
                    (-4*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 + &
                      5*MV2*shat*smT2**2 + 2*smT2**3) + &
                   E**y*(8*MV2**4*shat**4 - 10*MV2**3*shat**3*smT2 - &
                      MV2*shat**2*(16*MV2 + shat)*smT2**2 + &
                      shat*(20*MV2 + shat)*smT2**3 - 8*smT2**4))*LC(11)) &
                /(E**y*shat**2*(MV2*shat - smT2)**3) + &
              ((2*E**(3*y)*shat**2*smT**3*(MV2*shat + 4*smT2) - &
                   E**(2*y)*shat*smT2*(3*MV2*shat + smT2)* &
                    (MV2*shat + 4*smT2) + &
                   shat*(MV2*shat - smT2)* &
                    (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2) + &
                   E**y*smT*(-(MV2**2*shat**3*(3*MV2 + 2*shat)) + &
                      3*MV2*shat**2*(5*MV2 + shat)*smT2 - &
                      shat*(4*MV2 + shat)*smT2**2 + 2*smT2**3))*LC(15))/ &
               (shat*(MV2*shat - smT2)**3) + &
              (2*(E**(3*y)*shat**2*smT**3 + &
                   shat*(MV2*shat - smT2)*(2*MV2*shat - smT2) - &
                   E**(2*y)*shat*smT2*(MV2*shat + smT2) - &
                   E**y*smT*(MV2*shat**2*(MV2 + 2*shat) - &
                      2*shat*(2*MV2 + shat)*smT2 + 2*smT2**2))*LC(7)* &
                 SP(1))/(shat**2*(-(MV2*shat) + smT2)**2) + &
              (2*(-2*MV2*shat**2 + &
                   E**y*smT*(2*MV2*shat + shat**2 - 2*smT2) + shat*smT2) &
                  *LC(6)*SP(2))/(shat**2*(MV2*shat - smT2)) + &
              (LC(13)*(2*E**(3*y)*shat*smT*smT2*(-(MV2*shat) + smT2)**2* &
                    SP(3) - 4*E**y*MV2*shat**2*smT**3* &
                    (MV2*shat + 3*smT2)*SP(4) + &
                   4*smT2*(-2*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
                      3*MV2*shat*smT2**2 + smT2**3)*SP(4) - &
                   2*E**(2*y)*shat**2*smT2* &
                    ((-(MV2*shat) + smT2)**2*SP(3) - &
                      4*MV2*shat*smT2*SP(4))))/ &
               (E**y*shat**2*smT*(MV2*shat - smT2)**3) - &
              (4*LC(3)*(E**(4*y)*shat**2*smT**3*(3*MV2*shat - smT2)* &
                    smT2*SP(3) - &
                   E**(3*y)*shat*smT2* &
                    ((-(MV2*shat*smT) + smT**3)**2 + &
                      MV2*shat*(5*MV2*shat - smT2)*smT2)*SP(3) + &
                   2*smT**3*(-2*MV2*shat + smT2)* &
                    (MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) + &
                   E**(2*y)*smT* &
                    (smT2*(-3*MV2**2*shat**4 + &
                         2*MV2*shat**2*(3*MV2 + 2*shat)*smT2 - &
                         shat*(6*MV2 + shat)*smT2**2 + 2*smT2**3)*SP(3) &
                       + MV2*shat**3*smT2*(-3*MV2*shat + smT2)*SP(4)) + &
                   E**y*shat*smT2* &
                    (-2*MV2**2*shat**2*smT2*(SP(3) - 2*SP(4)) - &
                      3*MV2*shat*smT2**2*SP(4) + smT2**3*SP(4) + &
                      2*MV2**3*shat**3*(SP(3) + SP(4)))))/ &
               (shat**3*smT*(MV2*shat - smT2)**3) + &
              (2*LC(9)*(4*E**(5*y)*shat**2*smT2**3*SP(3) - &
                   2*E**(4*y)*shat*smT**3*smT2*(3*MV2*shat + smT2)* &
                    SP(3) + E**y*MV2*shat**2*smT2* &
                    (MV2*(MV2 - shat)*shat + (MV2 + shat)*smT2)*SP(4) + &
                   shat*smT*(MV2*shat - smT2)* &
                    (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)*SP(4) &
                     - 2*E**(3*y)*smT2**2* &
                    ((2*MV2*shat**3 - shat*(3*MV2 + 2*shat)*smT2 + &
                         smT2**2)*SP(3) - MV2*shat**3*SP(4)) + &
                   E**(2*y)*shat*smT*smT2* &
                    (2*MV2**2*shat**2*(SP(3) - 2*SP(4)) + &
                      MV2*shat*smT2*SP(4) - smT2**2*(2*SP(3) + SP(4))))) &
                /(shat**2*smT*(MV2*shat - smT2)**3) + &
              (2*E**y*smT*(MV2*shat*(-MV2 + shat) + &
                   E**y*smT*(MV2*shat + smT2 - 2*shat*smT*Cosh(y)))* &
                 LC(5)*SP(5))/(shat*(-(MV2*shat) + smT2)**2) + &
              (2*smT*(MV2*shat*smT*(MV2*shat - smT2) + &
                   2*E**(3*y)*shat*(3*MV2*shat - 2*smT2)*smT2 - &
                   E**(2*y)*smT* &
                    (5*MV2**2*shat**2 + MV2*shat*smT2 - 2*smT2**2) + &
                   E**y*shat* &
                    (MV2**2*(MV2 - 3*shat)*shat + &
                      MV2*(MV2 + 5*shat)*smT2 - 2*smT2**2))*LC(4)*SP(6)) &
                /(shat*(MV2*shat - smT2)**3) + &
              (4*E**(2*y)*(-MV2 + E**y*smT)* &
                 (MV2*shat + smT*(-2*E**y*shat + smT))*smT2*LC(10)*SP(6) &
                 )/(-(MV2*shat) + smT2)**3 + &
              (4*(-2*MV2**3*shat**3 - &
                   2*MV2**2*shat**2*(E**y*shat - 3*smT)*smT + &
                   smT**5*(-(E**y*shat) + smT) + &
                   MV2*shat*(2*E**y*shat - 3*smT)*(E**y*shat + smT)*smT2 &
                   )*LC(14)*SP(6))/(shat*(MV2*shat - smT2)**3) - &
              (2*smT*(E**(2*y)*(-5*MV2*shat*smT**3 + smT**5) + &
                   3*MV2*shat*smT*(-(MV2*shat) + smT2) + &
                   E**(3*y)*shat*smT2*(MV2*shat + smT2) + &
                   E**y*shat* &
                    (-(MV2*shat*smT2) - smT2**2 + &
                      2*MV2**2*(shat**2 + smT2)))*LC(2)*SP(7))/ &
               (shat*(MV2*shat - smT2)**3) + &
              (4*E**y*smT2*(2*MV2*shat*smT - 2*smT**3 + &
                   2*E**(3*y)*shat*smT2 - &
                   E**(2*y)*smT*(3*MV2*shat + smT2) + &
                   E**y*(MV2*(MV2 - 2*shat)*shat + (MV2 + 2*shat)*smT2)) &
                  *LC(8)*SP(7))/(-(MV2*shat) + smT2)**3 + &
              (4*(2*MV2**3*shat**3 + &
                   3*MV2*shat*smT**3*(E**y*shat + smT) - &
                   6*MV2**2*shat**2*smT2 - &
                   smT2**2*(E**y*shat*(2*E**y*shat - smT) + smT2))* &
                 LC(12)*SP(7))/(shat*(MV2*shat - smT2)**3)) + &
           (4*kapPr*(shat*smT**3*(MV2*shat - smT2)*SP(1)*SP(7) + &
                E**(2*y)*shat*smT*smT2*(SP(3) + SP(4))* &
                 (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)) + &
                E**y*smT2*(smT2*(-(MV2*shat) + smT2)*(SP(3) + SP(4))* &
                    SP(5) + shat* &
                    (shat*(-(MV2*shat) + smT2)*SP(1) + &
                      2*(smT2*SP(3) + MV2*shat*SP(4))*SP(6))*SP(7))))/ &
            (shat*smT*(-(MV2*shat) + smT2)**2))

     
     
return              
end function               



#endif 



! function calc_MassiveBox_QP(MV2_in,MH2_in,MT2_in,shat_in,smT_in,y_in,LC_in,SP_in,SI_in,kap_in,kapPr_in)
! implicit none
! real(16),parameter :: E=exp(1q0)
! real(8) :: shat_in,smT_in,y_in,MV2_in,MH2_in,MT2_in
! complex(8) :: calc_MassiveBox_QP, LC_in(1:15),SI_in(0:16),SP_in(1:9),kap_in,kapPr_in
! real(16) :: shat,smT,y,MV2,MH2,MT2,smT2
! complex(16) :: calc_MassiveBox, LC(1:15),SI(0:16),SP(1:9),kap,kapPr
! 
! 
! shat=real(shat_in,kind=16)
! smT=real(smT_in,kind=16)
! y=real(y_in,kind=16)
! MV2=real(MV2_in,kind=16)
! MH2=real(MH2_in,kind=16)
! MT2=real(MT2_in,kind=16)
! LC(1:15)=cmplx(LC_in(1:15),kind=16)
! SI(0:16)=cmplx(SI_in(0:16),kind=16)
! SP(1:9)=cmplx(SP_in(1:9),kind=16)
! kap=cmplx(kap_in,kind=16)
! kapPr=cmplx(kapPr_in,kind=16)
! 
! 
!      smT2 = smT**2
! ! print *,"QP",shat,smT,y,MV2,MH2,MT2,smT2,SI(0:16),LC(1:15)
! 
!      
!      calc_MassiveBox =  kap*SI(2)*((-8*LC(1))/shat - &
!            (8*E**y*smT*LC(11))/(shat**2 - E**y*shat*smT) + &
!            (16*E**y*LC(3)*SP(3))/(E**y*shat**2 - shat*smT) - &
!            (16*E**y*LC(9)*SP(3))/(E**y*shat**2 - shat*smT) + &
!            (16*LC(13)*SP(4))/(shat**2 - E**y*shat*smT)) + &
!         kap*SI(7)*((2*E**y*(-MV2 + E**y*smT)* &
!               (16*MV2**2*shat**2 + &
!                 MV2*shat*(3*E**y*shat - 37*smT)*smT + &
!                 9*smT**3*(E**y*shat + smT))*LC(1))/ &
!             (shat*(E**y*shat - smT)*(-(MV2*shat) + smT2)**2) + &
!            (2*(-(smT2/(E**y*shat - smT)) + &
!                 (2*smT**3*(MV2*shat + 2*smT2))/ &
!                  (-(MV2*shat) + smT2)**2 + &
!                 (MV2*(-8*MV2**2*shat**2 + 7*MV2*shat*smT2 - 5*smT2**2))/ &
!                  (E**y*(-(MV2*shat) + smT2)**2))*LC(11))/(shat*smT) - &
!            (2*(2*MV2*smT*(3*MV2*shat + 2*smT2) + &
!                 E**(2*y)*shat*smT*(MV2*shat + 9*smT2) - &
!                 E**y*(5*MV2**2*shat**2 + 7*MV2*shat*smT2 + 8*smT2**2))* &
!               LC(15))/((E**y*shat - smT)*(-(MV2*shat) + smT2)**2) - &
!            (4*(MV2 - E**y*smT)*LC(7)*SP(1))/(shat*(MV2*shat - smT2)) + &
!            (8*(MV2 - E**y*smT)*LC(6)*SP(2))/(shat*(MV2*shat - smT2)) + &
!            (16*MV2*(MV2 - E**y*smT)*LC(13)*SP(4))/ &
!             (E**(2*y)*(-(MV2*shat) + smT2)**2) - &
!            (4*(-MV2 + E**y*smT)*smT2*LC(9)* &
!               (4*E**(3*y)*shat*smT*(MV2*shat - 3*smT2)*SP(3) + &
!                 4*E**(4*y)*shat**2*smT2*SP(3) + 2*MV2*shat*smT2*SP(4) + &
!                 E**(2*y)*(-4* &
!                     (MV2*shat**2*(3*MV2 + shat) - &
!                       shat*(5*MV2 + shat)*smT2 + smT2**2)*SP(3) + &
!                    shat**2*(MV2*shat + smT2)*SP(4)) + &
!                 E**y*shat*smT* &
!                  (4*(MV2*shat - smT2)*SP(3) - (3*MV2*shat + smT2)*SP(4)) &
!                 ))/ &
!             (shat*smT**2*(-(E**y*shat) + smT)**2* &
!               (-(MV2*shat) + smT2)**2) + &
!            (16*(-MV2 + E**y*smT)*LC(3)* &
!               (E**(4*y)*shat**2*smT*(2*MV2*shat - smT2)*smT2*SP(3) - &
!                 E**(3*y)*shat*smT2* &
!                  (4*MV2**2*shat**2 - 5*MV2*shat*smT2 + 3*smT2**2)*SP(3) &
!                  + smT**5*(-2*MV2*shat + smT2)*SP(4) + &
!                 E**(2*y)*smT*smT2* &
!                  ((MV2*(MV2 - shat)*shat**2 + shat*(-MV2 + shat)*smT2 + &
!                       smT2**2)*SP(3) + &
!                    shat**2*(-2*MV2*shat + smT2)*SP(4)) + &
!                 E**y*shat*smT2**2* &
!                  (-(smT2*(SP(3) + 2*SP(4))) + &
!                    MV2*shat*(SP(3) + 4*SP(4)))))/ &
!             (E**y*shat**2*smT**2*(-(E**y*shat) + smT)**2* &
!               (-(MV2*shat) + smT2)**2) - &
!            (4*(MV2 - E**y*smT)*LC(5)*SP(5))/(shat*(MV2*shat - smT2)) + &
!            (4*(MV2 - E**y*smT)*(-2*MV2*smT + E**y*(MV2*shat + smT2))* &
!               LC(4)*SP(6))/((E**y*shat - smT)*(-(MV2*shat) + smT2)**2) &
!             + (16*E**y*smT*(MV2 - E**y*smT)*LC(10)*SP(6))/ &
!             (-(MV2*shat) + smT2)**2 + &
!            (16*(MV2**2*shat - E**y*smT**3)*LC(14)*SP(6))/ &
!             (E**y*smT*(-(MV2*shat) + smT2)**2) + &
!            (4*smT*(-MV2 + E**y*smT)* &
!               (MV2*shat + smT*(-2*E**y*shat + smT))*LC(2)*SP(7))/ &
!             (shat*(-(E**y*shat) + smT)*(-(MV2*shat) + smT2)**2) + &
!            (16*E**y*smT*(-MV2 + E**y*smT)*LC(8)*SP(7))/ &
!             (-(MV2*shat) + smT2)**2 + &
!            (16*(-(MV2**2*shat) + E**y*smT**3)*LC(12)*SP(7))/ &
!             (E**y*smT*(-(MV2*shat) + smT2)**2)) + &
!         kap*SI(6)*((4*smT*(-(E**y*MV2) + smT)* &
!               (3*E**y*smT*(MV2*shat + smT2) - shat*(5*MV2*shat + smT2))* &
!               LC(1))/(shat*(shat - E**y*smT)*(-(MV2*shat) + smT2)**2) + &
!            ((6*MV2)/(MV2*shat - smT2) + &
!               (-30*MV2**2*shat**2*smT + 56*MV2*shat*smT**3 - 14*smT**5)/ &
!                (E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
!               (16*(-(MV2*shat) + smT2))/(shat*(shat - E**y*smT)**2) - &
!               (2*smT2*(MV2*shat + 5*smT2))/ &
!                (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
!               (2*(4*shat*(-6*MV2 + shat) + 7*smT2))/ &
!                (shat**2*(shat - E**y*smT)))*LC(11) + &
!            2*(-5/(shat - E**y*smT) + &
!               (2*E**y*MV2*smT*(MV2*shat - 6*smT2) + &
!                  E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
!                  smT2*(3*MV2*shat + 7*smT2))/ &
!                (E**y*smT*(-(MV2*shat) + smT2)**2))*LC(15) - &
!            (4*(E**(2*y)*MV2*shat*smT + smT**3 - &
!                 E**y*MV2*(shat**2 + smT2))*LC(7)*SP(1))/ &
!             (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) + &
!            (8*(E**(2*y)*MV2*shat*smT + smT**3 - &
!                 E**y*MV2*(shat**2 + smT2))*LC(6)*SP(2))/ &
!             (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) - &
!            (4*(E**y*MV2 - smT)*LC(9)* &
!               (-4*E**(2*y)*MV2*shat**2*SP(3) + &
!                 4*E**(3*y)*MV2*shat*smT*SP(3) + &
!                 shat*(-3*MV2*shat + smT2)*SP(4) + &
!                 E**y*smT*(MV2*shat + smT2)*SP(4)))/ &
!             (E**y*shat*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2) + &
!            (16*(E**y*MV2 - smT)*LC(3)* &
!               (E**(4*y)*smT**3*(2*MV2*shat - smT2)*SP(3) - &
!                 2*E**(3*y)*shat*(2*MV2*shat - smT2)*smT2*SP(3) - &
!                 MV2*shat**3*smT*SP(4) + &
!                 2*E**y*shat*(MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)* &
!                  SP(4) - E**(2*y)*smT* &
!                  (shat**2*(-2*MV2*shat + smT2)*SP(3) + &
!                    (MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4))))/ &
!             (E**(2*y)*shat**2*(shat - E**y*smT)**2* &
!               (-(MV2*shat) + smT2)**2) + &
!            4*LC(13)*((E**y*MV2*SP(3))/(MV2*shat*smT - smT**3) - &
!               (MV2*SP(3))/(MV2*shat**2 - shat*smT2) - &
!               (4*smT*(4*MV2*shat + shat**2 - 4*smT2)*SP(4))/ &
!                (E**y*shat**3*(MV2*shat - smT2)) + &
!               (12*(MV2*shat - smT2)*SP(4))/ &
!                (shat**2*(shat - E**y*smT)**2) + &
!               (4*smT**3*SP(4))/ &
!                (E**(3*y)*shat*(-(MV2*shat) + smT2)**2) - &
!               (4*smT2**2*SP(4))/ &
!                (E**(2*y)*shat**2*(-(MV2*shat) + smT2)**2) + &
!               (-16*smT2*SP(4) + shat**2*(SP(3) + 4*SP(4)))/ &
!                (shat**3*(shat - E**y*smT))) - &
!            (4*(E**(2*y)*MV2*shat*smT + smT**3 - &
!                 E**y*MV2*(shat**2 + smT2))*LC(5)*SP(5))/ &
!             (shat*smT*(shat - E**y*smT)*(MV2*shat - smT2)) + &
!            4*(1/(-shat**2 + E**y*shat*smT) + &
!               (-2*E**y*MV2**2*shat*smT + &
!                  E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
!                  (3*MV2*shat - smT2)*smT2)/ &
!                (E**y*shat*smT*(-(MV2*shat) + smT2)**2))*LC(4)*SP(6) + &
!            (16*smT*(-(E**y*MV2) + smT)*LC(10)*SP(6))/ &
!             (-(MV2*shat) + smT2)**2 + &
!            (16*(E**y*MV2 - smT)*(shat*smT + E**y*(MV2*shat - 2*smT2))* &
!               LC(14)*SP(6))/ &
!             (E**(2*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2) - &
!            (4*(2*E**(2*y)*MV2*shat*(MV2*shat - smT2) + &
!                 3*MV2*shat*smT2 - smT2**2 + &
!                 E**y*MV2*smT*(-3*MV2*shat + smT2))*LC(2)*SP(7))/ &
!             (E**y*shat*smT*(-(MV2*shat) + smT2)**2) - &
!            (16*smT*(-(E**y*MV2) + smT)*LC(8)*SP(7))/ &
!             (-(MV2*shat) + smT2)**2 + &
!            (16*(-(E**y*MV2) + smT)* &
!               (shat*smT + E**y*(MV2*shat - 2*smT2))*LC(12)*SP(7))/ &
!             (E**(2*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2)) + &
!         kap*SI(4)*((MV2*(-16*MV2**2*shat**2*smT**3 + &
!                 37*MV2*shat*smT**5 - 9*smT**7 - &
!                 6*E**(8*y)*smT**5*(MV2*shat + smT2) + &
!                 E**(7*y)*shat*smT2**2*(MV2*shat + 23*smT2) - &
!                 E**(5*y)*shat*smT2* &
!                  (20*MV2**2*shat**2 + 131*MV2*shat*smT2 - 79*smT2**2) - &
!                 E**(4*y)*smT*(MV2*shat - smT2)* &
!                  (16*MV2**2*shat**2 + 32*MV2*shat*smT2 - 45*smT2**2) + &
!                 E**(6*y)*smT**3* &
!                  (46*MV2**2*shat**2 + 5*MV2*shat*smT2 - 27*smT2**2) + &
!                 E**y*shat*smT2* &
!                  (56*MV2**2*shat**2 - 113*MV2*shat*smT2 + 33*smT2**2) - &
!                 E**(3*y)*shat* &
!                  (80*MV2**3*shat**3 - 348*MV2**2*shat**2*smT2 + &
!                    285*MV2*shat*smT2**2 - 89*smT2**3) - &
!                 E**(2*y)*smT* &
!                  (8*MV2**3*shat**3 + 86*MV2**2*shat**2*smT2 - &
!                    103*MV2*shat*smT2**2 + 33*smT2**3))*LC(1))/ &
!             (2q0*E**(3*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (MV2*(3*E**y*shat*smT**5*(MV2*shat - 5*smT2) - &
!                 3*E**(9*y)*shat*smT**5*(MV2*shat - smT2) + &
!                 8*MV2**2*shat**2*smT2**2 - 7*MV2*shat*smT2**3 + &
!                 5*smT2**4 + 2*E**(8*y)*smT2**2* &
!                  (12*MV2**2*shat**2 - 20*MV2*shat*smT2 + 5*smT2**2) - &
!                 2*E**(7*y)*shat*smT**3* &
!                  (16*MV2**2*shat**2 - 35*MV2*shat*smT2 + 13*smT2**2) - &
!                 2*E**(3*y)*shat*smT**3* &
!                  (10*MV2**2*shat**2 - 59*MV2*shat*smT2 + 31*smT2**2) + &
!                 4*E**(5*y)*shat*smT* &
!                  (8*MV2**3*shat**3 - 47*MV2**2*shat**2*smT2 + &
!                    49*MV2*shat*smT2**2 - 19*smT2**3) + &
!                 E**(2*y)*smT2* &
!                  (-64*MV2**3*shat**3 + 88*MV2**2*shat**2*smT2 - &
!                    61*MV2*shat*smT2**2 + 25*smT2**3) + &
!                 E**(6*y)*smT2* &
!                  (-56*MV2**3*shat**3 + 144*MV2**2*shat**2*smT2 - &
!                    111*MV2*shat*smT2**2 + 35*smT2**3) + &
!                 E**(4*y)*(-(MV2*shat) + smT2)* &
!                  (-128*MV2**3*shat**3 + 112*MV2**2*shat**2*smT2 - &
!                    80*MV2*shat*smT2**2 + 45*smT2**3))*LC(11))/ &
!             (2q0*E**(5*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (MV2*(-6*MV2*shat*smT**5 - 4*smT**7 + &
!                 E**(6*y)*(-6*MV2**2*shat**2*smT**3 - &
!                    46*MV2*shat*smT**5 + 32*smT**7) - &
!                 2*E**(8*y)*smT**5*(MV2*shat - 6*smT2) + &
!                 10*E**(2*y)*MV2*shat*smT**3*(7*MV2*shat - 5*smT2) - &
!                 2*E**y*shat*(3*MV2*shat - 13*smT2)*smT2**2 - &
!                 E**(9*y)*shat*(MV2*shat - smT2)*smT2**2 - &
!                 E**(3*y)*shat*smT2* &
!                  (32*MV2**2*shat**2 + 57*MV2*shat*smT2 - 29*smT2**2) + &
!                 E**(7*y)*shat*smT2* &
!                  (8*MV2**2*shat**2 - 7*MV2*shat*smT2 - 21*smT2**2) - &
!                 E**(5*y)*shat*(4*MV2*shat - 19*smT2)* &
!                  (4*MV2**2*shat**2 + MV2*shat*smT2 - smT2**2) + &
!                 8*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
!                  (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2))* &
!               LC(15))/ &
!             (2q0*E**(4*y)*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
!            (4*MV2*(-(E**y*shat) + smT)*LC(7)*SP(1))/ &
!             (shat*smT*(-(MV2*shat) + smT2)) + &
!            (LC(6)*(8*E**y*MV2*shat*SP(2) - 8*MV2*smT*SP(2)))/ &
!             (MV2*shat**2*smT - shat*smT**3) + &
!            (MV2*LC(9)*(-8*E**(9*y)*shat*smT**7*SP(3) + &
!                 4*E**(10*y)*MV2*shat*smT2**3*SP(3) - &
!                 2*MV2*shat*smT2**3*SP(4) + &
!                 2*E**(7*y)*shat*smT**5* &
!                  (6*(5*MV2*shat - 3*smT2)*SP(3) - &
!                    (MV2*shat + smT2)*SP(4)) - &
!                 2*E**y*shat*smT**5* &
!                  (2*(MV2*shat - smT2)*SP(3) - (MV2*shat + smT2)*SP(4)) &
!                  + E**(8*y)*smT2**2* &
!                  (4*(-8*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2)* &
!                     SP(3) + smT2*(MV2*shat + smT2)*SP(4)) - &
!                 2*E**(3*y)*shat*smT**3*(2*MV2*shat - smT2)* &
!                  (2*(MV2*shat - 3*smT2)*SP(3) + &
!                    (5*MV2*shat + smT2)*SP(4)) + &
!                 E**(4*y)*smT2* &
!                  (4*smT2*(4*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
!                       3*smT2**2)*SP(3) + &
!                    (-(MV2*shat) + smT2)* &
!                     (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(4)) - &
!                 2*E**(5*y)*shat*smT**3* &
!                  (2*MV2**2*shat**2*(22*SP(3) - 5*SP(4)) + &
!                    3*MV2*shat*smT2*(-18*SP(3) + SP(4)) + &
!                    smT2**2*(22*SP(3) + SP(4))) + &
!                 E**(6*y)*smT2* &
!                  (64*MV2**3*shat**3*SP(3) + &
!                    3*smT2**3*(4*SP(3) + SP(4)) + &
!                    MV2*shat*smT2**2*(16*SP(3) + SP(4)) - &
!                    4*MV2**2*shat**2*smT2*(23*SP(3) + 2*SP(4))) + &
!                 E**(2*y)*smT2**2* &
!                  (smT2**2*(4*SP(3) + SP(4)) + &
!                    4*MV2**2*shat**2*(3*SP(3) + 4*SP(4)) - &
!                    MV2*shat*smT2*(20*SP(3) + 13*SP(4)))))/ &
!             (E**(4*y)*shat*smT**2*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*MV2*smT2*LC(3)* &
!               (E**(9*y)*shat*smT**3*(3*MV2*shat - smT2)*smT2*SP(3) + &
!                 E**(10*y)*smT2**3*(-2*MV2*shat + smT2)*SP(3) + &
!                 smT2**3*(-2*MV2*shat + smT2)*SP(4) - &
!                 2*E**(3*y)*shat*smT**3* &
!                  ((4*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
!                     SP(3) + MV2*shat*(7*MV2*shat - 4*smT2)*SP(4)) + &
!                 E**(4*y)*smT2* &
!                  (2*(-5*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 - &
!                       7*MV2*shat*smT2**2 + 2*smT2**3)*SP(3) - &
!                    3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
!                     (2*MV2*shat - smT2)*SP(4)) + &
!                 2*E**(5*y)*shat*smT* &
!                  ((2*MV2*shat - smT2)* &
!                     (7*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(3) + MV2*shat* &
!                     (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
!                     SP(4)) + &
!                 E**(6*y)*smT2* &
!                  (-3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
!                     (2*MV2*shat - smT2)*SP(3) + &
!                    2*(-5*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 - &
!                       7*MV2*shat*smT2**2 + 2*smT2**3)*SP(4)) - &
!                 2*E**(7*y)*shat*smT*smT2* &
!                  (11*MV2**2*shat**2*SP(3) + 2*smT2**2*SP(3) + &
!                    MV2*shat*smT2*(-10*SP(3) + SP(4))) + &
!                 E**(8*y)*smT2**2* &
!                  (smT2**2*(4*SP(3) + SP(4)) + &
!                    MV2**2*shat**2*(14*SP(3) + SP(4)) - &
!                    MV2*shat*smT2*(16*SP(3) + SP(4))) + &
!                 E**y*shat*smT**5* &
!                  (-(smT2*SP(3)) + MV2*shat*(SP(3) + 2*SP(4))) + &
!                 E**(2*y)*smT2**2* &
!                  (smT2**2*(SP(3) + 4*SP(4)) + &
!                    MV2**2*shat**2*(SP(3) + 14*SP(4)) - &
!                    MV2*shat*smT2*(SP(3) + 16*SP(4)))))/ &
!             (E**(5*y)*shat**2*smT**3*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (MV2*LC(13)*(E**(10*y)*smT**5*(MV2*shat - smT2)*SP(3) - &
!                 E**(11*y)*shat*(MV2*shat - smT2)*smT2**2*SP(3) - &
!                 4*MV2*shat*smT**5*SP(4) - &
!                 12*E**(3*y)*shat*(5*MV2*shat - 3*smT2)*smT2**2*SP(4) + &
!                 8*E**y*shat*smT2**3*SP(4) + &
!                 4*E**(2*y)*smT**3* &
!                  (8*MV2**2*shat**2 - 5*MV2*shat*smT2 - smT2**2)*SP(4) + &
!                 E**(9*y)*shat*(MV2*shat - smT2)*smT2* &
!                  (8*MV2*shat*SP(3) - 5*smT2*SP(3) + 4*smT2*SP(4)) - &
!                 E**(7*y)*shat* &
!                  ((MV2*shat - smT2)* &
!                     (16*MV2**2*shat**2 - 20*MV2*shat*smT2 + 7*smT2**2)* &
!                     SP(3) - 4*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
!                     smT2*SP(4)) - &
!                 E**(8*y)*smT**3* &
!                  (3*(MV2*shat - smT2)*(2*MV2*shat - smT2)*SP(3) + &
!                    4*(3*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
!                     SP(4)) - &
!                 E**(4*y)*smT*(MV2*shat - smT2)* &
!                  (2*MV2*shat*smT2*(SP(3) - 14*SP(4)) + &
!                    64*MV2**2*shat**2*SP(4) - smT2**2*(SP(3) + 12*SP(4))) &
!                   + E**(6*y)*smT* &
!                  (8*MV2**3*shat**3*SP(3) - &
!                    16*MV2**2*shat**2*smT2*(SP(3) + SP(4)) - &
!                    3*smT2**3*(SP(3) + 4*SP(4)) + &
!                    MV2*shat*smT2**2*(11*SP(3) + 20*SP(4))) + &
!                 E**(5*y)*shat*smT2* &
!                  (4*MV2**2*shat**2*(SP(3) + 22*SP(4)) + &
!                    smT2**2*(3*SP(3) + 44*SP(4)) - &
!                    MV2*shat*smT2*(7*SP(3) + 108*SP(4)))))/ &
!             (E**(6*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
!            (4*MV2*(-(E**y*shat) + smT)*LC(5)*SP(5))/ &
!             (shat*smT*(-(MV2*shat) + smT2)) + &
!            (MV2*(-2*MV2*smT**5 + 2*E**(8*y)*MV2*smT**5 + &
!                 2*E**(2*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) - &
!                 2*E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) + &
!                 E**(9*y)*smT2**2*(-(MV2*shat) + smT2) + &
!                 2*E**y*smT2**2*(MV2*shat + smT2) + &
!                 E**(7*y)*smT2* &
!                  (8*MV2**2*shat**2 - 15*MV2*shat*smT2 + 3*smT2**2) + &
!                 E**(5*y)*(-4*MV2*shat + smT2)* &
!                  (4*MV2**2*shat**2 - 13*MV2*shat*smT2 + 5*smT2**2) + &
!                 E**(3*y)*smT2* &
!                  (-16*MV2**2*shat**2 - MV2*shat*smT2 + 5*smT2**2))* &
!               LC(4)*SP(6))/ &
!             (E**(4*y)*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-(MV2*smT**5) + E**(8*y)*MV2*smT**5 - &
!                 E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) + &
!                 6*E**(4*y)*MV2*smT*(-(MV2*shat) + smT2)**2 - &
!                 E**(7*y)*smT2**2*(MV2*shat + smT2) + &
!                 E**(5*y)*smT2* &
!                  (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
!                 E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
!                    smT2**2) + &
!                 E**(2*y)*MV2*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
!                 E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(10)*SP(6))/ &
!             (E**(3*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-2*E**(4*y)*MV2**2*shat*(8*MV2*shat - 5*smT2)* &
!                  (MV2*shat - smT2) + &
!                 E**(6*y)*MV2*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
!                  smT2 - MV2**2*shat*smT2**2 - &
!                 E**(8*y)*MV2*(MV2*shat - 2*smT2)*smT2**2 + &
!                 E**y*smT**5*(MV2*shat + smT2) - &
!                 E**(3*y)*smT**3* &
!                  (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) + &
!                 E**(2*y)*MV2*smT2* &
!                  (8*MV2**2*shat**2 - 5*MV2*shat*smT2 - smT2**2) + &
!                 E**(7*y)*smT**3* &
!                  (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(5*y)*smT* &
!                  (4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(14)*SP(6))/ &
!             (E**(5*y)*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (MV2*(MV2*shat*smT**5 + smT**7 + &
!                 E**(8*y)*(-3*MV2*shat*smT**5 + smT**7) + &
!                 E**(6*y)*(22*MV2**2*shat**2*smT**3 - &
!                    22*MV2*shat*smT**5 + 4*smT**7) + &
!                 2*E**(9*y)*shat*(MV2*shat - smT2)*smT2**2 + &
!                 2*E**(2*y)*smT**3*(MV2*shat + smT2)* &
!                  (-3*MV2*shat + 2*smT2) - &
!                 E**y*shat*smT2**2*(MV2*shat + 3*smT2) + &
!                 E**(3*y)*shat*smT2* &
!                  (8*MV2**2*shat**2 + 13*MV2*shat*smT2 - 9*smT2**2) + &
!                 2*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
!                  (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2) - &
!                 E**(7*y)*shat*smT2* &
!                  (16*MV2**2*shat**2 - 27*MV2*shat*smT2 + 7*smT2**2) + &
!                 E**(5*y)*shat* &
!                  (32*MV2**3*shat**3 - 88*MV2**2*shat**2*smT2 + &
!                    55*MV2*shat*smT2**2 - 11*smT2**3))*LC(2)*SP(7))/ &
!             (E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(MV2*smT**5 - E**(8*y)*MV2*smT**5 + &
!                 E**(6*y)*MV2*smT**3*(7*MV2*shat - 5*smT2) - &
!                 6*E**(4*y)*MV2*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**(7*y)*smT2**2*(MV2*shat + smT2) + &
!                 E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
!                    smT2**2) - &
!                 E**(2*y)*MV2*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(5*y)*smT2* &
!                  (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
!                 E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(8)*SP(7))/ &
!             (E**(3*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(2*E**(4*y)*MV2**2*shat*(8*MV2*shat - 5*smT2)* &
!                  (MV2*shat - smT2) - &
!                 E**(6*y)*MV2*(MV2*shat - 3*smT2)*(2*MV2*shat - smT2)* &
!                  smT2 + MV2**2*shat*smT2**2 + &
!                 E**(8*y)*MV2*(MV2*shat - 2*smT2)*smT2**2 - &
!                 E**y*smT**5*(MV2*shat + smT2) + &
!                 E**(3*y)*smT**3* &
!                  (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
!                 E**(7*y)*smT**3* &
!                  (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(2*y)*MV2*smT2* &
!                  (-8*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2) - &
!                 E**(5*y)*smT* &
!                  (4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(12)*SP(7))/ &
!             (E**(5*y)*smT*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
!         kap*SI(3)*(((-2*shat*smT**7*(5*MV2*shat + smT2) + &
!                 3*E**(12*y)*shat*smT**7*(MV2*shat + 3*smT2) + &
!                 E**(11*y)*smT2**3* &
!                  (MV2*(13*MV2 - 6*shat)*shat**2 - &
!                    2*shat*(23*MV2 + 9*shat)*smT2 + 9*smT2**2) + &
!                 E**(10*y)*smT**5* &
!                  (MV2*shat**2*(-16*MV2**2 - 41*MV2*shat + 3*shat**2) + &
!                    shat*(37*MV2**2 + 35*MV2*shat + 9*shat**2)*smT2 + &
!                    3*(-3*MV2 + 10*shat)*smT2**2) + &
!                 2*E**y*smT2**2* &
!                  (-8*MV2**2*shat**4 + &
!                    MV2*shat**2*(5*MV2 + 26*shat)*smT2 + &
!                    2*(2*MV2 - 3*shat)*shat*smT2**2 + 3*smT2**3) + &
!                 2*E**(9*y)*smT2**2* &
!                  (3*MV2**2*shat**3*(-15*MV2 + 7*shat) + &
!                    MV2*shat**2*(184*MV2 + 33*shat)*smT2 - &
!                    2*shat*(74*MV2 + 21*shat)*smT2**2 + 21*smT2**3) + &
!                 E**(8*y)*smT**3* &
!                  (2*MV2**2*shat**3* &
!                     (64*MV2**2 + 77*MV2*shat - 9*shat**2) - &
!                    MV2*shat**2*(358*MV2**2 + 595*MV2*shat + 47*shat**2)* &
!                     smT2 + shat* &
!                     (239*MV2**2 + 324*MV2*shat + 41*shat**2)*smT2**2 + &
!                    33*(-MV2 + shat)*smT2**3) + &
!                 E**(2*y)*smT**3* &
!                  (16*MV2**2*shat**4*(8*MV2 + shat) - &
!                    2*MV2*shat**3*(94*MV2 + 21*shat)*smT2 + &
!                    shat*(-6*MV2**2 + 45*MV2*shat + 14*shat**2)* &
!                     smT2**2 - 3*(2*MV2 + 3*shat)*smT2**3) + &
!                 E**(6*y)*smT*(-(MV2*shat) + smT2)* &
!                  (16*MV2**2*shat**3* &
!                     (16*MV2**2 + 32*MV2*shat + 3*shat**2) - &
!                    32*MV2*shat**2*(16*MV2**2 + 34*MV2*shat + 3*shat**2)* &
!                     smT2 + shat* &
!                     (304*MV2**2 + 607*MV2*shat + 69*shat**2)*smT2**2 + &
!                    (-45*MV2 + 7*shat)*smT2**3) + &
!                 E**(3*y)*smT2* &
!                  (-8*MV2**3*shat**4*(17*MV2 + 15*shat) + &
!                    6*MV2**2*shat**3*(11*MV2 + shat)*smT2 + &
!                    MV2*shat**2*(139*MV2 + 156*shat)*smT2**2 - &
!                    6*shat*(21*MV2 + 11*shat)*smT2**3 + 33*smT2**4) + &
!                 2*E**(7*y)*smT2* &
!                  (16*MV2**3*shat**4*(5*MV2 + shat) - &
!                    MV2**2*shat**3*(301*MV2 + 15*shat)*smT2 + &
!                    MV2*shat**2*(504*MV2 + 101*shat)*smT2**2 - &
!                    2*shat*(149*MV2 + 39*shat)*smT2**3 + 39*smT2**4) - &
!                 E**(4*y)*smT* &
!                  (8*MV2**3*(19*MV2 - shat)*shat**5 - &
!                    2*MV2**2*shat**3* &
!                     (68*MV2**2 + 459*MV2*shat + 61*shat**2)*smT2 + &
!                    MV2*shat**2* &
!                     (226*MV2**2 + 1033*MV2*shat + 157*shat**2)*smT2**2 &
!                     - 3*shat*(47*MV2**2 + 121*MV2*shat + 17*shat**2)* &
!                     smT2**3 + 3*(9*MV2 + 4*shat)*smT2**4) + &
!                 2*E**(5*y)*(16*MV2**4*shat**5*(8*MV2 + 5*shat) - &
!                    4*MV2**3*shat**4*(69*MV2 + 61*shat)*smT2 + &
!                    MV2**2*shat**3*(-71*MV2 + 95*shat)*smT2**2 + &
!                    3*MV2*shat**2*(133*MV2 + 39*shat)*smT2**3 - &
!                    24*shat*(10*MV2 + 3*shat)*smT2**4 + 36*smT2**5))* &
!               LC(1))/ &
!             (2q0*E**(4*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)* &
!               (-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            ((-2*E**(13*y)*shat*smT2**4*(MV2*shat + 2*smT2) - &
!                 shat**2*smT**7*(MV2*shat + 5*smT2) + &
!                 E**y*shat*smT2**3* &
!                  (MV2*shat**2*(-15*MV2 + 2*shat) + &
!                    10*shat*(3*MV2 + shat)*smT2 + 3*smT2**2) + &
!                 E**(12*y)*smT**5* &
!                  (8*MV2**2*(MV2 - 3*shat)*shat**3 + &
!                    3*MV2*shat**2*(-2*MV2 + 17*shat)*smT2 + &
!                    (5*MV2 - 9*shat)*shat*smT2**2 + 5*smT2**3) + &
!                 E**(2*y)*smT**5* &
!                  (MV2*shat**3*(16*MV2**2 + 27*MV2*shat - shat**2) - &
!                    shat**2*(MV2**2 + 15*MV2*shat + 5*shat**2)*smT2 - &
!                    shat*(43*MV2 + 36*shat)*smT2**2 + 10*smT2**3) + &
!                 E**(3*y)*smT2**2* &
!                  (MV2**2*(104*MV2 - 21*shat)*shat**4 - &
!                    MV2*shat**2*(24*MV2**2 + 308*MV2*shat + 43*shat**2)* &
!                     smT2 + shat*(40*MV2**2 + 199*MV2*shat + 58*shat**2)* &
!                     smT2**2 + 5*(-2*MV2 + shat)*smT2**3) + &
!                 E**(11*y)*smT2**2* &
!                  (8*MV2**2*shat**4*(16*MV2 + 5*shat) - &
!                    8*MV2*shat**2*(MV2**2 + 30*MV2*shat + 10*shat**2)* &
!                     smT2 + shat*(7*MV2**2 + 107*MV2*shat + 22*shat**2)* &
!                     smT2**2 - (5*MV2 + 19*shat)*smT2**3) + &
!                 E**(10*y)*smT**3* &
!                  (-8*MV2**2*shat**4* &
!                     (25*MV2**2 + 23*MV2*shat + 2*shat**2) + &
!                    MV2*shat**3*(224*MV2**2 + 201*MV2*shat + 31*shat**2)* &
!                     smT2 + shat**2* &
!                     (85*MV2**2 + 39*MV2*shat - 9*shat**2)*smT2**2 - &
!                    shat*(145*MV2 + 56*shat)*smT2**3 + 30*smT2**4) + &
!                 E**(4*y)*smT**3* &
!                  (-2*MV2**2*shat**4* &
!                     (64*MV2**2 + 82*MV2*shat - 5*shat**2) + &
!                    2*MV2*shat**3* &
!                     (48*MV2**2 + 205*MV2*shat + 13*shat**2)*smT2 + &
!                    shat**2*(298*MV2**2 - 85*MV2*shat - 24*shat**2)* &
!                     smT2**2 - shat*(287*MV2 + 107*shat)*smT2**3 + &
!                    45*smT2**4) - &
!                 E**(5*y)*smT2* &
!                  (4*MV2**3*(28*MV2 - 15*shat)*shat**5 + &
!                    MV2**2*shat**3* &
!                     (-192*MV2**2 - 720*MV2*shat + 83*shat**2)*smT2 + &
!                    MV2*shat**2* &
!                     (416*MV2**2 + 1326*MV2*shat + 149*shat**2)*smT2**2 &
!                     - shat*(247*MV2**2 + 662*MV2*shat + 136*shat**2)* &
!                     smT2**3 + 5*(7*MV2 + 2*shat)*smT2**4) - &
!                 E**(9*y)*smT2* &
!                  (16*MV2**3*(MV2 - 4*shat)*shat**5 - &
!                    MV2**2*shat**3* &
!                     (200*MV2**2 + 780*MV2*shat + 63*shat**2)*smT2 + &
!                    MV2*shat**2* &
!                     (360*MV2**2 + 1267*MV2*shat + 197*shat**2)*smT2**2 &
!                     - shat*(197*MV2**2 + 594*MV2*shat + 94*shat**2)* &
!                     smT2**3 + (25*MV2 + 37*shat)*smT2**4) + &
!                 2*E**(6*y)*smT* &
!                  (4*MV2**3*shat**5* &
!                     (32*MV2**2 + 36*MV2*shat - 3*shat**2) - &
!                    2*MV2**3*shat**4*(116*MV2 + 257*shat)*smT2 + &
!                    3*MV2*shat**3* &
!                     (-26*MV2**2 + 190*MV2*shat + 11*shat**2)*smT2**2 + &
!                    shat**2*(454*MV2**2 - 106*MV2*shat - 21*shat**2)* &
!                     smT2**3 - 2*shat*(150*MV2 + 41*shat)*smT2**4 + &
!                    40*smT2**5) + &
!                 E**(8*y)*smT* &
!                  (8*MV2**3*shat**5* &
!                     (48*MV2**2 + 23*MV2*shat - shat**2) - &
!                    2*MV2**2*shat**4* &
!                     (496*MV2**2 + 484*MV2*shat + 21*shat**2)*smT2 + &
!                    2*MV2*shat**3*(14*MV2 + shat)*(13*MV2 + 35*shat)* &
!                     smT2**2 + &
!                    shat**2*(644*MV2**2 - 153*MV2*shat - 32*shat**2)* &
!                     smT2**3 - shat*(506*MV2 + 135*shat)*smT2**4 + &
!                    70*smT2**5) - &
!                 E**(7*y)*(32*MV2**4*shat**6*(8*MV2 + shat) + &
!                    4*MV2**3*shat**4* &
!                     (96*MV2**2 - 14*MV2*shat - 53*shat**2)*smT2 + &
!                    MV2**2*shat**3* &
!                     (-1024*MV2**2 - 1556*MV2*shat + 95*shat**2)*smT2**2 &
!                     + MV2*shat**2* &
!                     (992*MV2**2 + 2276*MV2*shat + 221*shat**2)*smT2**3 &
!                     - shat*(397*MV2**2 + 978*MV2*shat + 160*shat**2)* &
!                     smT2**4 + (45*MV2 + 34*shat)*smT2**5))*LC(11))/ &
!             (2q0*E**(6*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)**2* &
!               (-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            ((-(shat*smT2**3*(3*MV2*shat + 7*smT2)) + &
!                 E**(12*y)*shat*smT2**3*(MV2*shat + 9*smT2) - &
!                 E**(11*y)*smT**5* &
!                  (MV2*shat**2*(5*MV2 + 3*shat) + &
!                    shat*(7*MV2 + 17*shat)*smT2 + 8*smT2**2) + &
!                 E**y*smT**5*(3*MV2*shat**2*(MV2 + 2*shat) + &
!                    shat*(5*MV2 + 14*shat)*smT2 + 12*smT2**2) + &
!                 E**(10*y)*smT2**2* &
!                  (MV2*shat**3*(-19*MV2 + 2*shat) + &
!                    shat*(6*MV2**2 - MV2*shat + 8*shat**2)*smT2 + &
!                    4*(MV2 + 10*shat)*smT2**2) - &
!                 E**(8*y)*shat*smT2* &
!                  (2*MV2**2*shat**3*(41*MV2 + 18*shat) + &
!                    MV2*shat*(70*MV2**2 + 68*MV2*shat + shat**2)*smT2 - &
!                    (50*MV2**2 + 32*MV2*shat + 17*shat**2)*smT2**2 - &
!                    48*smT2**3) + &
!                 E**(9*y)*smT**3* &
!                  (MV2**2*shat**3*(62*MV2 + 55*shat) + &
!                    MV2*shat**2*(-12*MV2 + 19*shat)*smT2 - &
!                    2*shat*(5*MV2 + 27*shat)*smT2**2 - 20*smT2**3) + &
!                 E**(3*y)*smT**3* &
!                  (-(MV2**2*shat**3*(46*MV2 + 79*shat)) + &
!                    2*MV2*shat**2*(9*MV2 + 10*shat)*smT2 + &
!                    shat*(-32*MV2 + 39*shat)*smT2**2 + 40*smT2**3) - &
!                 E**(6*y)*(-(MV2*shat) + smT2)* &
!                  (-8*MV2**2*shat**4*(4*MV2 + shat) + &
!                    8*MV2*shat**2*(8*MV2**2 + 11*MV2*shat + shat**2)* &
!                     smT2 - shat*(16*MV2 + shat)*(4*MV2 + 3*shat)* &
!                     smT2**2 + (24*MV2 + 13*shat)*smT2**3) + &
!                 E**(4*y)*smT2* &
!                  (2*MV2**2*shat**4*(41*MV2 + 22*shat) + &
!                    MV2*shat**2*(6*MV2**2 + 28*MV2*shat - 11*shat**2)* &
!                     smT2 + shat*(46*MV2**2 + 22*MV2*shat - 13*shat**2)* &
!                     smT2**2 - 2*(16*MV2 + 31*shat)*smT2**3) + &
!                 E**(7*y)*shat*smT* &
!                  (4*MV2**3*shat**3*(-12*MV2 + 17*shat) + &
!                    MV2**2*(170*MV2 - 37*shat)*shat**2*smT2 + &
!                    MV2*shat*(-41*MV2 + 55*shat)*smT2**2 - &
!                    (41*MV2 + 46*shat)*smT2**3) + &
!                 E**(5*y)*smT* &
!                  (4*MV2**3*(4*MV2 - 5*shat)*shat**4 - &
!                    MV2**2*shat**3*(154*MV2 + 83*shat)*smT2 + &
!                    MV2*shat**2*(133*MV2 + 47*shat)*smT2**2 + &
!                    shat*(-75*MV2 + 16*shat)*smT2**3 + 40*smT2**4) - &
!                 E**(2*y)*smT2**2* &
!                  (3*MV2*(shat**2 + 2*smT2)**2 + &
!                    shat*smT2*(7*shat**2 + 39*smT2) - &
!                    MV2**2*(31*shat**3 + 2*shat*smT2)))*LC(15))/ &
!             (2q0*E**(5*y)*(E**y*shat - smT)*(-shat + E**y*smT)* &
!               (-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(7)*SP(1))/ &
!             (shat*(-shat + E**y*smT)*(MV2*shat - smT2)) + &
!            (8*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(6)*SP(2))/ &
!             (shat*(shat - E**y*smT)*(MV2*shat - smT2)) + &
!            (LC(9)*(4*E**(15*y)*shat**2*smT2**4*SP(3) - &
!                 4*E**(14*y)*shat*smT**7*(2*shat**2 + 3*smT2)*SP(3) + &
!                 2*shat*smT**7*(2*MV2*shat - smT2)*SP(4) - &
!                 E**y*smT2**3* &
!                  (MV2*shat**2*(-MH2 + 4*MV2 + 11*shat) + &
!                    (MH2 - MV2 - 5*shat)*shat*smT2 + smT2**2)*SP(4) - &
!                 2*E**(13*y)*smT2**3* &
!                  (2*(shat**2*(4*MV2**2 + 7*MV2*shat - shat**2) - &
!                       shat*(8*MV2 + 11*shat)*smT2 + smT2**2)*SP(3) - &
!                    MV2*shat**3*SP(4)) + &
!                 E**(12*y)*smT**5* &
!                  (4*(MV2*shat**2*(3*MV2**2 + 5*MV2*shat + 14*shat**2) + &
!                       shat*(-5*MV2**2 + 9*MV2*shat - 13*shat**2)*smT2 + &
!                       (MV2 - 14*shat)*smT2**2)*SP(3) - &
!                    shat**2*(MV2*shat*(-MH2 + 2*MV2 + 3*shat) + &
!                       (MH2 + 4*MV2 + shat)*smT2)*SP(4)) + &
!                 E**(2*y)*smT**5* &
!                  (4*MV2*shat**2*smT2*SP(3) + &
!                    (-(MV2*(MH2 + 21*MV2 - 10*shat)*shat**3) + &
!                       shat*(MV2**2 + 36*MV2*shat + (MH2 - 4*shat)*shat)* &
!                        smT2 + (MV2 - 7*shat)*smT2**2)*SP(4)) + &
!                 E**(11*y)*smT2**2* &
!                  (4*(MV2*shat**3*(29*MV2**2 + 9*MV2*shat - 7*shat**2) + &
!                       shat**2*(-75*MV2**2 - 57*MV2*shat + 5*shat**2)* &
!                        smT2 + shat*(46*MV2 + 39*shat)*smT2**2 - &
!                       4*smT2**3)*SP(3) + &
!                    shat*(MV2*shat**3*(-14*MV2 + shat) + &
!                       shat*(-(MH2*MV2) + 4*MV2**2 + 21*MV2*shat + &
!                          shat**2)*smT2 + (MH2 + 2*MV2 + shat)*smT2**2)* &
!                     SP(4)) + &
!                 E**(3*y)*smT2**2* &
!                  (4*shat*(5*MV2**2*shat**3 - &
!                       MV2*shat*(MV2 + 12*shat)*smT2 - &
!                       (MV2 - 4*shat)*smT2**2)*SP(3) + &
!                    (MV2*shat**3* &
!                        (-8*(MH2 - 4*MV2)*MV2 + 76*MV2*shat - 3*shat**2) &
!                        + shat**2* &
!                        (11*MH2*MV2 - 34*MV2**2 - 97*MV2*shat + shat**2)* &
!                        smT2 + shat*(-3*MH2 + 8*MV2 + 21*shat)*smT2**2 - &
!                       4*smT2**3)*SP(4)) + &
!                 E**(10*y)*smT**3* &
!                  (-4*(MV2**2*shat**3* &
!                        (24*MV2**2 + 41*MV2*shat + 19*shat**2) - &
!                       2*MV2*shat**2* &
!                        (26*MV2**2 + 43*MV2*shat + 29*shat**2)*smT2 + &
!                       shat*(29*MV2**2 + 7*MV2*shat + 33*shat**2)* &
!                        smT2**2 - 3*(MV2 - 9*shat)*smT2**3)*SP(3) + &
!                    shat*(MV2**2*shat**3*(-8*MH2 + 16*MV2 + 19*shat) + &
!                       MV2*(11*MH2 + 15*MV2 - 9*shat)*shat**2*smT2 - &
!                       (2*MV2**2 + 32*MV2*shat + shat*(3*MH2 + 8*shat))* &
!                        smT2**2 + smT2**3)*SP(4)) - &
!                 E**(4*y)*smT**3* &
!                  (4*shat*(2*MV2**2*shat**3*(13*MV2 + 4*shat) - &
!                       MV2*shat**2*(41*MV2 + 17*shat)*smT2 + &
!                       (-MV2**2 + 8*MV2*shat + 6*shat**2)*smT2**2 + &
!                       3*smT2**3)*SP(3) + &
!                    (MV2**2*shat**4*(-8*MH2 + 36*MV2 + 73*shat) + &
!                       MV2*shat**2* &
!                        (8*MV2**2 + 30*MV2*shat + &
!                          (11*MH2 - 82*shat)*shat)*smT2 + &
!                       shat*(-MV2**2 - 56*MV2*shat + &
!                          shat*(-3*MH2 + 17*shat))*smT2**2 + &
!                       (-3*MV2 + 8*shat)*smT2**3)*SP(4)) + &
!                 E**(9*y)*smT2* &
!                  (-4*(MV2**2*shat**4* &
!                        (40*MV2**2 - 5*MV2*shat - 9*shat**2) - &
!                       2*MV2*shat**3* &
!                        (71*MV2**2 + 19*MV2*shat - 9*shat**2)*smT2 + &
!                       shat**2*(189*MV2**2 + 112*MV2*shat - 9*shat**2)* &
!                        smT2**2 - shat*(89*MV2 + 65*shat)*smT2**3 + &
!                       6*smT2**4)*SP(3) + &
!                    (2*MV2**2*(2*MV2 - 3*shat)*shat**5 - &
!                       2*MV2*shat**3* &
!                        (-4*(MH2 - 4*MV2)*MV2 + 39*MV2*shat + shat**2)* &
!                        smT2 + &
!                       shat**2* &
!                        (-11*MH2*MV2 + 12*MV2**2 + 47*MV2*shat + &
!                          4*shat**2)*smT2**2 + &
!                       shat*(3*MH2 + 13*MV2 + 9*shat)*smT2**3 - smT2**4)* &
!                     SP(4)) - &
!                 2*E**(5*y)*smT2* &
!                  (2*(-3*MV2**2*shat**4*(MV2 + shat)*(8*MV2 + shat) + &
!                       2*MV2*shat**3*(7*MV2**2 + 3*MV2*shat + 3*shat**2)* &
!                        smT2 + &
!                       shat**2*(27*MV2**2 + 45*MV2*shat - 2*shat**2)* &
!                        smT2**2 - shat*(19*MV2 + 24*shat)*smT2**3 + &
!                       smT2**4)*SP(3) + &
!                    (MV2**2*shat**4* &
!                        (-8*MH2*MV2 + 20*MV2**2 + 20*MV2*shat - &
!                          11*shat**2) + &
!                       MV2*shat**3* &
!                        (12*MH2*MV2 - 58*MV2**2 - 95*MV2*shat + &
!                          11*shat**2)*smT2 + &
!                       shat**2* &
!                        (-5*(MH2 - 8*MV2)*MV2 + 81*MV2*shat - 2*shat**2)* &
!                        smT2**2 + &
!                       (MH2 - 11*MV2 - 17*shat)*shat*smT2**3 + 3*smT2**4) &
!                      *SP(4)) + &
!                 E**(6*y)*smT* &
!                  (4*(4*MV2**3*(8*MV2 - shat)*shat**5 - &
!                       MV2**2*shat**3* &
!                        (24*MV2**2 + 175*MV2*shat + 37*shat**2)*smT2 + &
!                       2*MV2*shat**2* &
!                        (20*MV2**2 + 97*MV2*shat + 31*shat**2)*smT2**2 - &
!                       shat*(19*MV2**2 + 46*MV2*shat + 25*shat**2)* &
!                        smT2**3 + (MV2 - 14*shat)*smT2**4)*SP(3) + &
!                    (4*MV2**3*shat**5*(-4*MH2 + 20*MV2 + 15*shat) - &
!                       MV2**2*shat**3* &
!                        (8*MV2**2 + 188*MV2*shat + &
!                          3*shat*(-8*MH2 + 53*shat))*smT2 + &
!                       MV2*shat**2* &
!                        (16*MV2**2 + 118*MV2*shat + &
!                          shat*(-10*MH2 + 119*shat))*smT2**2 - &
!                       shat*(11*MV2**2 + 16*MV2*shat - &
!                          2*(MH2 - 14*shat)*shat)*smT2**3 + &
!                       (3*MV2 - 2*shat)*smT2**4)*SP(4)) + &
!                 E**(8*y)*smT* &
!                  (4*(4*MV2**3*shat**4* &
!                        (12*MV2**2 + 24*MV2*shat + shat**2) - &
!                       2*MV2**2*shat**3* &
!                        (64*MV2**2 + 147*MV2*shat + 28*shat**2)*smT2 + &
!                       MV2*shat**2* &
!                        (121*MV2**2 + 274*MV2*shat + 89*shat**2)*smT2**2 &
!                        - shat*(44*MV2**2 + 53*MV2*shat + 41*shat**2)* &
!                        smT2**3 + 3*(MV2 - 9*shat)*smT2**4)*SP(3) + &
!                    (4*MV2**3*shat**5*(4*MH2 - 2*MV2 + 5*shat) + &
!                       MV2**2*(-24*MH2 + 8*MV2 - 27*shat)*shat**4*smT2 + &
!                       MV2*shat**2* &
!                        (16*MV2**2 + 88*MV2*shat + &
!                          shat*(10*MH2 + 41*shat))*smT2**2 - &
!                       shat*(13*MV2**2 + 68*MV2*shat + &
!                          2*shat*(MH2 + 11*shat))*smT2**3 + &
!                       (MV2 + 2*shat)*smT2**4)*SP(4)) - &
!                 2*E**(7*y)*(2* &
!                     (16*MV2**4*shat**5*(3*MV2 + 2*shat) - &
!                       4*MV2**2*shat**4* &
!                        (24*MV2**2 + 26*MV2*shat + 3*shat**2)*smT2 + &
!                       MV2*shat**3* &
!                        (-27*MV2**2 + 22*MV2*shat + 17*shat**2)*smT2**2 &
!                        + shat**2*(136*MV2**2 + 95*MV2*shat - 7*shat**2)* &
!                        smT2**3 - shat*(71*MV2 + 56*shat)*smT2**4 + &
!                       4*smT2**5)*SP(3) + &
!                    (4*MV2**3*shat**6*(5*MV2 + 2*shat) - &
!                       2*MV2**2*shat**4* &
!                        (-4*MH2*MV2 + 4*MV2**2 + 15*MV2*shat + &
!                          8*shat**2)*smT2 + &
!                       MV2*shat**3* &
!                        (-12*MH2*MV2 + 14*MV2**2 + 3*MV2*shat + &
!                          11*shat**2)*smT2**2 + &
!                       shat**2* &
!                        (5*MH2*MV2 + 9*MV2**2 + 24*MV2*shat - 3*shat**2)* &
!                        smT2**3 - shat*(MH2 + 13*(MV2 + shat))*smT2**4 + &
!                       2*smT2**5)*SP(4))))/ &
!             (E**(5*y)*shat*(-(E**y*shat) + smT)**2*(-shat + E**y*smT)* &
!               (-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (LC(13)*(-(E**(13*y)*shat*(MV2*shat - smT2)*smT2**3*SP(3)) - &
!                 4*shat**2*smT**7*SP(4) + &
!                 4*E**y*shat*smT2**3*(shat**2 + 3*smT2)*SP(4) + &
!                 4*E**(2*y)*smT**5* &
!                  (MV2*shat**2*(4*MV2 + 7*shat) - &
!                    8*shat*(MV2 + shat)*smT2 + smT2**2)*SP(4) + &
!                 E**(12*y)*smT**5* &
!                  ((MV2*shat - smT2)*(3*MV2*shat + 2*shat**2 - smT2)* &
!                     SP(3) + 4*MV2*shat*smT2*SP(4)) - &
!                 E**(4*y)*smT**3* &
!                  ((MV2*shat - smT2)* &
!                     (smT2**2 + MV2*shat*(2*shat**2 + smT2))*SP(3) + &
!                    4*(2*MV2**2*shat**3*(16*MV2 + 5*shat) - &
!                       40*MV2*shat**2*(2*MV2 + shat)*smT2 + &
!                       shat*(47*MV2 + 24*shat)*smT2**2 - 4*smT2**3)*SP(4) &
!                    ) - 2*E**(10*y)*smT**3* &
!                  ((MV2*shat - smT2)* &
!                     (MV2*shat**2*(6*MV2 + 5*shat) - &
!                       3*shat*(2*MV2 + shat)*smT2 + 2*smT2**2)*SP(3) + &
!                    2*(2*MV2**2*shat**3*(13*MV2 + 4*shat) - &
!                       MV2*shat**2*(46*MV2 + 17*shat)*smT2 + &
!                       2*shat*(10*MV2 + 3*shat)*smT2**2 - smT2**3)*SP(4)) &
!                   + E**(9*y)*smT2* &
!                  ((MV2*shat - smT2)* &
!                     (2*MV2*shat**2*(MV2 + shat)*(4*MV2 + 3*shat) - &
!                       shat*(8*MV2**2 + 7*MV2*shat + 3*shat**2)*smT2 + &
!                       (3*MV2 + shat)*smT2**2)*SP(3) + &
!                    4*(3*MV2**2*shat**3*(MV2 + shat)*(8*MV2 + shat) - &
!                       2*MV2*shat**2* &
!                        (20*MV2**2 + 7*MV2*shat + 3*shat**2)*smT2 + &
!                       shat*(19*MV2**2 - 28*MV2*shat + 2*shat**2)* &
!                        smT2**2 - (MV2 - 18*shat)*smT2**3)*SP(4)) + &
!                 E**(5*y)*smT2* &
!                  ((MV2*shat - smT2)* &
!                     (2*MV2*shat**3*(MV2 + shat) - &
!                       shat**2*(5*MV2 + shat)*smT2 + &
!                       (MV2 + 4*shat)*smT2**2)*SP(3) + &
!                    4*(3*MV2**2*shat**3* &
!                        (8*MV2**2 + 3*MV2*shat + 3*shat**2) - &
!                       2*MV2*shat**2* &
!                        (26*MV2**2 + 3*MV2*shat + 9*shat**2)*smT2 + &
!                       shat*(29*MV2**2 - 40*MV2*shat + 9*shat**2)* &
!                        smT2**2 + (-3*MV2 + 31*shat)*smT2**3)*SP(4)) - &
!                 2*E**(8*y)*smT* &
!                  (-((shat - smT)*(shat + smT)*(4*MV2*shat - 3*smT2)* &
!                       (-(MV2*shat) + smT2)**2*SP(3)) + &
!                    2*(4*MV2**3*shat**4*(-8*MV2 + shat) + &
!                       2*MV2**2*shat**3*(74*MV2 + 17*shat)*smT2 - &
!                       4*MV2*shat**2*(45*MV2 + 14*shat)*smT2**2 + &
!                       shat*(74*MV2 + 23*shat)*smT2**3 - 4*smT2**4)*SP(4) &
!                    ) + 2*E**(6*y)*smT* &
!                  ((MV2*shat - smT2)* &
!                     (4*MV2**2*shat**4 + &
!                       MV2*(2*MV2 - 3*shat)*shat**2*smT2 + &
!                       shat*(2*MV2 + shat)*smT2**2 - 2*smT2**3)*SP(3) + &
!                    4*(2*MV2**3*shat**4*(16*MV2 + shat) - &
!                       MV2**2*shat**3*(97*MV2 + 22*shat)*smT2 + &
!                       MV2*shat**2*(109*MV2 + 36*shat)*smT2**2 - &
!                       shat*(46*MV2 + 17*shat)*smT2**3 + 3*smT2**4)*SP(4) &
!                    ) - E**(7*y)* &
!                  ((MV2*shat - smT2)* &
!                     (8*MV2**2*shat**4*(MV2 + shat) - &
!                       8*MV2*shat**3*(MV2 + shat)*smT2 + &
!                       shat*(6*MV2**2 + 13*MV2*shat + 3*shat**2)* &
!                        smT2**2 - (3*MV2 + 5*shat)*smT2**3)*SP(3) + &
!                    4*(16*MV2**4*shat**4*(3*MV2 + 2*shat) - &
!                       4*MV2**2*shat**3* &
!                        (32*MV2**2 + 25*MV2*shat + 3*shat**2)*smT2 + &
!                       MV2*shat**2* &
!                        (121*MV2**2 + 56*MV2*shat + 17*shat**2)*smT2**2 &
!                        - (11*MV2 - 7*shat)*(4*MV2 - shat)*shat* &
!                        smT2**3 + 3*(MV2 - 11*shat)*smT2**4)*SP(4)) - &
!                 E**(11*y)*smT2**2* &
!                  ((MV2*shat - smT2)* &
!                     (2*MV2**2*shat + shat**3 + 2*shat*smT2 - &
!                       MV2*(shat**2 + smT2))*SP(3) + &
!                    4*shat*(12*MV2*shat*smT2 - 4*smT2**2 + &
!                       MV2**2*(-5*shat**2 + smT2))*SP(4)) - &
!                 E**(3*y)*smT2**2* &
!                  (shat*smT2*(-(MV2*shat) + smT2)*SP(3) + &
!                    4*(3*MV2**3*shat**2 - &
!                       5*shat*smT2*(shat**2 + 3*smT2) + &
!                       MV2**2*(shat**3 - 5*shat*smT2) + &
!                       MV2*(7*shat**4 + 17*shat**2*smT2 + smT2**2))*SP(4) &
!                    )))/ &
!             (E**(7*y)*shat*(shat - E**y*smT)**2*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*LC(3)*(-(E**(16*y)*shat**2*(2*MV2*shat - smT2)*smT2**4* &
!                    SP(3)) + 3*E**(15*y)*shat*smT**7* &
!                  (2*MV2*shat**2*(MV2 + shat) - &
!                    shat*(2*MV2 + shat)*smT2 + smT2**2)*SP(3) - &
!                 MV2*shat**3*smT2**4*SP(4) + &
!                 E**y*shat*smT**7* &
!                  (3*MV2*shat**2*(MV2 + shat) - 2*MV2*shat*smT2 + &
!                    2*smT2**2)*SP(4) - &
!                 E**(14*y)*smT2**3* &
!                  ((2*MV2*shat**3*(2*MV2**2 + MV2*shat + 3*shat**2) + &
!                       shat**2*(-4*MV2**2 + 2*MV2*shat - 3*shat**2)* &
!                        smT2 + shat*(2*MV2 + 3*shat)*smT2**2 + smT2**3)* &
!                     SP(3) + shat**2*smT2*(-2*MV2*shat + smT2)*SP(4)) - &
!                 E**(2*y)*smT2**3* &
!                  (shat**2*smT2*(-2*MV2*shat + smT2)*SP(3) + &
!                    (MV2*shat**3*(2*MV2**2 + 2*MV2*shat + 3*shat**2) - &
!                       MV2**2*shat**2*smT2 + &
!                       shat*(MV2 + 5*shat)*smT2**2 + smT2**3)*SP(4)) + &
!                 E**(3*y)*smT**5* &
!                  (shat*(MV2**2*shat**4 - &
!                       MV2*shat**2*(2*MV2 + 7*shat)*smT2 + &
!                       3*shat*(-MV2 + shat)*smT2**2 + 2*smT2**3)*SP(3) + &
!                    (MV2*shat**4*(-18*MV2**2 - 11*MV2*shat + shat**2) + &
!                       MV2*shat**2*(MV2**2 + 34*MV2*shat + 10*shat**2)* &
!                        smT2 - &
!                       shat*(MV2**2 + 25*MV2*shat - 4*shat**2)*smT2**2 + &
!                       (MV2 + 12*shat)*smT2**3)*SP(4)) - &
!                 E**(4*y)*smT2**2* &
!                  ((2*MV2**3*shat**5 + &
!                       4*MV2*(2*MV2 - shat)*shat**4*smT2 + &
!                       shat**2*(-4*MV2**2 - 28*MV2*shat + shat**2)* &
!                        smT2**2 + 11*shat**2*smT2**3 + smT2**4)*SP(3) + &
!                    (-2*MV2**2*shat**4* &
!                        (8*MV2**2 + 27*MV2*shat + 8*shat**2) + &
!                       MV2*shat**3* &
!                        (22*MV2**2 + 75*MV2*shat + 12*shat**2)*smT2 + &
!                       shat**2*(-10*MV2**2 - 55*MV2*shat + shat**2)* &
!                        smT2**2 + shat*(-4*MV2 + 25*shat)*smT2**3 + &
!                       5*smT2**4)*SP(4)) + &
!                 E**(11*y)*smT**3* &
!                  ((2*MV2**2*shat**5* &
!                        (9*MV2**2 - 26*MV2*shat - 7*shat**2) + &
!                       MV2*shat**3* &
!                        (-8*MV2**3 - 150*MV2**2*shat + 35*MV2*shat**2 + &
!                          16*shat**3)*smT2 + &
!                       shat**2* &
!                        (14*MV2**3 + 212*MV2**2*shat - 3*MV2*shat**2 - &
!                          4*shat**3)*smT2**2 + &
!                       shat*(-12*MV2**2 - 135*MV2*shat + 5*shat**2)* &
!                        smT2**3 + 4*(MV2 + 10*shat)*smT2**4)*SP(3) + &
!                    (MV2**2*shat**5*(2*MV2**2 + 9*MV2*shat - shat**2) + &
!                       2*MV2**2*shat**4*(5*MV2 + 6*shat)*smT2 + &
!                       3*MV2*(MV2 - 12*shat)*shat**3*smT2**2 - &
!                       2*shat*(MV2**2 + 15*MV2*shat - 6*shat**2)* &
!                        smT2**3 + (MV2 + 12*shat)*smT2**4)*SP(4)) + &
!                 E**(5*y)*smT**3* &
!                  ((MV2**2*shat**5*(2*MV2**2 - 7*MV2*shat - shat**2) + &
!                       MV2*shat**4*(14*MV2**2 + 51*MV2*shat + shat**2)* &
!                        smT2 - &
!                       shat**3*(MV2**2 + 64*MV2*shat + shat**2)* &
!                        smT2**2 + &
!                       shat*(-2*MV2**2 - 31*MV2*shat + 17*shat**2)* &
!                        smT2**3 + (MV2 + 13*shat)*smT2**4)*SP(3) + &
!                    (2*MV2**2*shat**5* &
!                        (MV2**2 - 22*MV2*shat - 3*shat**2) - &
!                       MV2*(2*MV2 - shat)*shat**3* &
!                        (4*MV2**2 + 45*MV2*shat + 4*shat**2)*smT2 + &
!                       2*MV2*shat**2* &
!                        (7*MV2**2 + 61*MV2*shat - 13*shat**2)*smT2**2 + &
!                       shat*(-12*MV2**2 - 83*MV2*shat + 18*shat**2)* &
!                        smT2**3 + 2*(2*MV2 + 15*shat)*smT2**4)*SP(4)) - &
!                 E**(10*y)*smT2* &
!                  ((4*MV2**3*shat**5* &
!                        (16*MV2**2 + 47*MV2*shat + 6*shat**2) - &
!                       2*MV2**2*shat**4* &
!                        (86*MV2**2 + 285*MV2*shat + 41*shat**2)*smT2 + &
!                       2*MV2*shat**3* &
!                        (69*MV2**2 + 290*MV2*shat + 33*shat**2)*smT2**2 &
!                        - 2*shat**2* &
!                        (17*MV2**2 + 140*MV2*shat + 7*shat**2)*smT2**3 + &
!                       shat*(-12*MV2 + 67*shat)*smT2**4 + 10*smT2**5)* &
!                     SP(3) + (4*MV2**3*shat**6*(7*MV2 + 5*shat) + &
!                       2*MV2**2*shat**4* &
!                        (2*MV2**2 - 11*MV2*shat - 5*shat**2)*smT2 + &
!                       MV2*shat**3* &
!                        (20*MV2**2 + 70*MV2*shat - 11*shat**2)*smT2**2 + &
!                       shat**2*(-19*MV2**2 - 98*MV2*shat + 4*shat**2)* &
!                        smT2**3 + shat*(-7*MV2 + 31*shat)*smT2**4 + &
!                       5*smT2**5)*SP(4)) - &
!                 E**(6*y)*smT2* &
!                  ((4*MV2**3*shat**6*(3*MV2 + shat) + &
!                       2*MV2**2*shat**4* &
!                        (2*MV2**2 - 11*MV2*shat + 7*shat**2)*smT2 + &
!                       2*MV2*shat**3* &
!                        (11*MV2**2 + 61*MV2*shat - 8*shat**2)*smT2**2 + &
!                       shat**2*(-22*MV2**2 - 144*MV2*shat + shat**2)* &
!                        smT2**3 + shat*(-6*MV2 + 41*shat)*smT2**4 + &
!                       5*smT2**5)*SP(3) + &
!                    (4*MV2**3*shat**5* &
!                        (8*MV2**2 + 19*MV2*shat - 2*shat**2) - &
!                       2*MV2**2*shat**4* &
!                        (46*MV2**2 + 141*MV2*shat + 5*shat**2)*smT2 + &
!                       2*MV2*shat**3* &
!                        (31*MV2**2 + 152*MV2*shat + 4*shat**2)*smT2**2 + &
!                       4*shat**3*(-41*MV2 + shat)*smT2**3 + &
!                       3*shat*(-6*MV2 + 17*shat)*smT2**4 + 10*smT2**5)* &
!                     SP(4)) + &
!                 E**(13*y)*smT**5* &
!                  ((2*MV2*shat**4*(-19*MV2**2 - 14*MV2*shat + shat**2) + &
!                       shat**2* &
!                        (MV2**3 + 74*MV2**2*shat + 38*MV2*shat**2 - &
!                          shat**3)*smT2 - &
!                       shat*(MV2**2 + 50*MV2*shat + 7*shat**2)*smT2**2 + &
!                       (MV2 + 17*shat)*smT2**3)*SP(3) - &
!                    shat*(-2*smT2**2*(shat**2 + smT2) + &
!                       MV2*shat*smT2*(4*shat**2 + 3*smT2) + &
!                       MV2**2*(shat**4 + 2*shat**2*smT2))*SP(4)) + &
!                 E**(12*y)*smT2**2* &
!                  ((2*MV2**2*shat**4* &
!                        (16*MV2**2 + 50*MV2*shat + 19*shat**2) - &
!                       2*MV2*shat**3*(3*MV2 + shat)*(9*MV2 + 23*shat)* &
!                        smT2 + &
!                       shat**2*(30*MV2**2 + 94*MV2*shat + 11*shat**2)* &
!                        smT2**2 - 29*shat**2*smT2**3 - 5*smT2**4)*SP(3) &
!                     + (MV2*shat**3*smT2*(2*shat**2 + 23*smT2) - &
!                       smT2**2*(shat**4 + 9*shat**2*smT2 + smT2**2) + &
!                       MV2**2* &
!                        (2*shat**6 - 7*shat**4*smT2 + 4*shat**2*smT2**2)) &
!                      *SP(4)) + &
!                 E**(9*y)*smT* &
!                  (MV2*smT2**2* &
!                     ((27*shat**6 - 146*shat**4*smT2 - &
!                          169*shat**2*smT2**2 + 6*smT2**3)*SP(3) + &
!                       (2*shat**6 - 97*shat**4*smT2 - &
!                          90*shat**2*smT2**2 + 4*smT2**3)*SP(4)) + &
!                    16*MV2**5* &
!                     (shat**4*smT2*SP(3) + shat**6*(8*SP(3) + SP(4))) + &
!                    2*shat*smT2**3* &
!                     (-3*shat**4*SP(3) + 5*smT2**2*(5*SP(3) + 3*SP(4)) + &
!                       shat**2*smT2*(15*SP(3) + 14*SP(4))) + &
!                    2*MV2**4*shat**3* &
!                     (5*shat**2*smT2*(-23*SP(3) + SP(4)) + &
!                       smT2**2*(-23*SP(3) + SP(4)) + &
!                       shat**4*(97*SP(3) + 29*SP(4))) + &
!                    MV2**3*shat**2* &
!                     (2*shat**2*smT2**2*(SP(3) - 12*SP(4)) + &
!                       2*shat**6*(9*SP(3) + 5*SP(4)) + &
!                       smT2**3*(47*SP(3) + 10*SP(4)) - &
!                       shat**4*smT2*(469*SP(3) + 88*SP(4))) + &
!                    MV2**2*shat*smT2* &
!                     (73*shat**2*smT2**2*(3*SP(3) + SP(4)) - &
!                       shat**6*(39*SP(3) + 10*SP(4)) - &
!                       smT2**3*(23*SP(3) + 14*SP(4)) + &
!                       shat**4*smT2*(391*SP(3) + 114*SP(4)))) - &
!                 E**(8*y)*(32*MV2**5*shat**5* &
!                     (smT2*(SP(3) + SP(4)) + shat**2*(2*SP(3) + SP(4))) &
!                     + 8*MV2**4*shat**4*(shat - smT)*(shat + smT)* &
!                     (shat**2*(7*SP(3) + 3*SP(4)) + &
!                       2*smT2*(5*SP(3) + 4*SP(4))) + &
!                    smT2**4*(6*shat**4*(-SP(3) + SP(4)) + &
!                       10*smT2**2*(SP(3) + SP(4)) + &
!                       shat**2*smT2*(73*SP(3) + 54*SP(4))) - &
!                    2*MV2**3*shat**3*smT2* &
!                     (2*shat**4*(23*SP(3) + 7*SP(4)) - &
!                       smT2**2*(43*SP(3) + 27*SP(4)) + &
!                       shat**2*smT2*(175*SP(3) + 101*SP(4))) + &
!                    2*MV2**2*shat**2*smT2**2* &
!                     (-(smT2**2*(13*SP(3) + 3*SP(4))) + &
!                       shat**4*(17*SP(3) + 7*SP(4)) + &
!                       shat**2*smT2*(254*SP(3) + 139*SP(4))) + &
!                    MV2*shat*smT2**3* &
!                     (2*shat**4*(7*SP(3) - 5*SP(4)) - &
!                       4*smT2**2*(4*SP(3) + 5*SP(4)) - &
!                       shat**2*smT2*(304*SP(3) + 187*SP(4)))) + &
!                 E**(7*y)*smT* &
!                  (shat*smT2**3* &
!                     ((-4*shat**4 + 35*shat**2*smT2 + 35*smT2**2)* &
!                        SP(3) + 8*smT2*(4*shat**2 + 5*smT2)*SP(4)) + &
!                    MV2**2*shat*smT2* &
!                     (2*(-9*shat**6 + 141*shat**4*smT2 + &
!                          34*shat**2*smT2**2 - 7*smT2**3)*SP(3) + &
!                       (-7*shat**6 + 197*shat**4*smT2 + &
!                          167*shat**2*smT2**2 - 23*smT2**3)*SP(4)) + &
!                    MV2*smT2**2* &
!                     (2*(7*shat**6 - 84*shat**4*smT2 - &
!                          53*shat**2*smT2**2 + 2*smT2**3)*SP(3) + &
!                       (5*shat**6 - 98*shat**4*smT2 - &
!                          123*shat**2*smT2**2 + 6*smT2**3)*SP(4)) + &
!                    16*MV2**5* &
!                     (shat**4*smT2*SP(4) + shat**6*(SP(3) + 4*SP(4))) + &
!                    2*MV2**4*shat**3* &
!                     (smT2**2*(SP(3) - 23*SP(4)) + &
!                       shat**4*(45*SP(3) + 41*SP(4)) - &
!                       shat**2*smT2*(11*SP(3) + 51*SP(4))) + &
!                    MV2**3*shat**2* &
!                     (2*shat**2*smT2**2*(12*SP(3) - 23*SP(4)) + &
!                       2*shat**6*(5*SP(3) + SP(4)) + &
!                       smT2**3*(10*SP(3) + 47*SP(4)) - &
!                       shat**4*smT2*(224*SP(3) + 213*SP(4))))))/ &
!             (E**(6*y)*shat**2*(-(E**y*shat) + smT)**2* &
!               (shat - E**y*smT)**2*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*E**y*smT*(MV2 + shat - 2*smT*Cosh(y))*LC(5)*SP(5))/ &
!             (shat*(-shat + E**y*smT)*(MV2*shat - smT2)) + &
!            ((-3*MV2*shat*smT2**3 + smT2**4 + &
!                 E**(12*y)*smT2**3*(MV2*shat + smT2) + &
!                 E**y*smT**5*(3*MV2*shat*(MV2 + 2*shat) + &
!                    (MV2 - 2*shat)*smT2) - &
!                 E**(11*y)*smT**5* &
!                  (MV2*shat*(MV2 + 3*shat) + (3*MV2 + shat)*smT2) + &
!                 E**(10*y)*smT2**2* &
!                  (MV2*shat**2*(-3*MV2 + 2*shat) + &
!                    MV2*(2*MV2 + 3*shat)*smT2 + 4*smT2**2) + &
!                 E**(3*y)*smT**3* &
!                  (-(MV2**2*shat**2*(22*MV2 + 39*shat)) + &
!                    2*MV2*shat*(7*MV2 + 22*shat)*smT2 + &
!                    (4*MV2 - 9*shat)*smT2**2) + &
!                 E**(9*y)*smT**3* &
!                  (3*MV2**2*shat**2*(2*MV2 + 5*shat) + &
!                    MV2*(12*MV2 - 5*shat)*shat*smT2 - &
!                    2*(7*MV2 + 3*shat)*smT2**2) + &
!                 E**(8*y)*smT2* &
!                  (-2*MV2**2*shat**3*(13*MV2 + 6*shat) + &
!                    MV2*shat*(-14*MV2**2 - 4*MV2*shat + 7*shat**2)* &
!                     smT2 + (10*MV2**2 + 8*MV2*shat + shat**2)*smT2**2 + &
!                    8*smT2**3) + &
!                 E**(4*y)*smT2* &
!                  (2*MV2**2*shat**3*(13*MV2 + 10*shat) + &
!                    MV2*shat*(14*MV2**2 + 28*MV2*shat - 19*shat**2)* &
!                     smT2 + (-10*MV2**2 - 50*MV2*shat + 3*shat**2)* &
!                     smT2**2 + 10*smT2**3) + &
!                 E**(6*y)*(-(MV2*shat) + smT2)* &
!                  (8*MV2**2*shat**3*(4*MV2 + shat) - &
!                    8*MV2*shat**2*(3*MV2 + shat)*smT2 + &
!                    3*shat*(-4*MV2 + shat)*smT2**2 + 11*smT2**3) + &
!                 E**(7*y)*smT* &
!                  (4*MV2**3*shat**3*(4*MV2 + 9*shat) - &
!                    MV2**2*shat**2*(14*MV2 + 45*shat)*smT2 + &
!                    MV2*shat*(23*MV2 + 31*shat)*smT2**2 - &
!                    (17*MV2 + 14*shat)*smT2**3) + &
!                 E**(5*y)*smT* &
!                  (4*MV2**3*shat**3*(4*MV2 + 3*shat) - &
!                    3*MV2**2*shat**2*(22*MV2 + 25*shat)*smT2 + &
!                    MV2*shat*(45*MV2 + 71*shat)*smT2**2 - &
!                    (3*MV2 + 16*shat)*smT2**3) + &
!                 E**(2*y)*smT2**2* &
!                  (MV2**2*(15*shat**2 - 2*smT2) + &
!                    smT2*(shat**2 + 5*smT2) - &
!                    3*MV2*(shat**3 + 8*shat*smT2)))*LC(4)*SP(6))/ &
!             (E**(5*y)*(E**y*shat - smT)*(-shat + E**y*smT)* &
!               (-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(E**y*(MV2 + shat)*smT**5 - &
!                 E**(9*y)*(MV2 + shat)*smT**5 + &
!                 E**(7*y)*(MV2 + shat)*smT**3*(7*MV2*shat - 5*smT2) - &
!                 2*E**(6*y)*(5*MV2*shat - 4*smT2)*smT2**2 - smT2**3 + &
!                 E**(10*y)*smT2**3 - &
!                 E**(3*y)*(MV2 + shat)*smT*(-3*MV2*shat + smT2)* &
!                  (-2*MV2*shat + smT2) - &
!                 6*E**(5*y)*(MV2 + shat)*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**(8*y)*smT2**2*(-6*MV2*shat + 5*smT2) - &
!                 E**(2*y)*smT2* &
!                  (-4*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) + &
!                 E**(4*y)*(8*MV2**3*shat**3 - 4*MV2**2*shat**2*smT2 - &
!                    6*MV2*shat*smT2**2 + 4*smT2**3))*LC(10)*SP(6))/ &
!             (E**(4*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(E**(11*y)*smT**7 + shat*smT2**3 - &
!                 E**(10*y)*shat*smT2**2*(MV2**2 + 2*smT2) - &
!                 E**y*smT**5*(shat**2 + 2*smT2) + &
!                 E**(9*y)*shat*smT**3* &
!                  (-4*MV2**2*shat + 4*MV2*smT2 + shat*smT2) + &
!                 E**(2*y)*smT2**2* &
!                  (-(MV2*shat*(MV2 + 6*shat)) + 2*(MV2 + 3*shat)*smT2) + &
!                 E**(8*y)*MV2*smT2* &
!                  (2*MV2*shat**2*(7*MV2 + 5*shat) - &
!                    shat*(17*MV2 + 7*shat)*smT2 + 5*smT2**2) - &
!                 E**(3*y)*smT**3* &
!                  (MV2*(MV2 - 7*shat)*shat**2 + &
!                    shat*(-13*MV2 + 5*shat)*smT2 + 9*smT2**2) + &
!                 E**(4*y)*smT2* &
!                  (8*MV2**3*shat**2 - MV2*shat*(19*MV2 + 17*shat)*smT2 + &
!                    (9*MV2 + 13*shat)*smT2**2) - &
!                 2*E**(5*y)*smT* &
!                  (MV2**2*shat**3*(-5*MV2 + 3*shat) + &
!                    2*MV2*(5*MV2 - 3*shat)*shat**2*smT2 + &
!                    shat*(-13*MV2 + 3*shat)*smT2**2 + 7*smT2**3) - &
!                 E**(7*y)*smT* &
!                  (2*MV2**2*shat**3*(7*MV2 + 3*shat) - &
!                    MV2*shat**2*(MV2 + 5*shat)*smT2 + &
!                    shat*(-17*MV2 + shat)*smT2**2 + 8*smT2**3) + &
!                 2*E**(6*y)*(4*MV2**3*shat**3*(-2*MV2 + shat) + &
!                    MV2**2*shat**2*(19*MV2 + shat)*smT2 - &
!                    MV2*shat*(17*MV2 + 9*shat)*smT2**2 + &
!                    (6*MV2 + 5*shat)*smT2**3))*LC(14)*SP(6))/ &
!             (E**(6*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            ((-2*E**(11*y)*shat*smT**7 + &
!                 E**y*(MV2 + 2*shat)*smT**5*(3*MV2*shat - smT2) + &
!                 smT2**3*(-3*MV2*shat + smT2) + &
!                 E**(2*y)*smT2**2*(-3*MV2*shat + smT2)* &
!                  (-5*MV2*shat + shat**2 + 5*smT2) + &
!                 E**(10*y)*smT2**2* &
!                  (2*MV2*shat**3 + 3*MV2*shat*smT2 + smT2**2) - &
!                 E**(9*y)*smT**3* &
!                  (5*MV2**2*shat**3 + MV2*(MV2 - 14*shat)*shat*smT2 + &
!                    (MV2 + 9*shat)*smT2**2) - &
!                 E**(3*y)*smT**3* &
!                  (MV2**2*shat**2*(22*MV2 + 39*shat) - &
!                    2*MV2*shat*(11*MV2 + 21*shat)*smT2 + &
!                    (4*MV2 + 9*shat)*smT2**2) + &
!                 E**(8*y)*smT2* &
!                  (2*MV2**2*(MV2 - 6*shat)*shat**3 + &
!                    MV2*shat**2*(-17*MV2 + 7*shat)*smT2 + &
!                    shat*(4*MV2 + shat)*smT2**2 + 5*smT2**3) + &
!                 E**(4*y)*smT2* &
!                  (2*MV2**2*shat**3*(13*MV2 + 10*shat) + &
!                    MV2*(5*MV2 - 19*shat)*shat**2*smT2 + &
!                    3*shat*(-11*MV2 + shat)*smT2**2 + 10*smT2**3) + &
!                 E**(7*y)*smT* &
!                  (44*MV2**3*shat**4 + &
!                    MV2**2*(6*MV2 - 61*shat)*shat**2*smT2 + &
!                    2*MV2*shat*(MV2 + 21*shat)*smT2**2 - &
!                    (4*MV2 + 17*shat)*smT2**3) + &
!                 E**(5*y)*smT* &
!                  (4*MV2**3*shat**3*(4*MV2 + 3*shat) - &
!                    MV2**2*shat**2*(32*MV2 + 63*shat)*smT2 + &
!                    2*MV2*shat*(11*MV2 + 32*shat)*smT2**2 - &
!                    (6*MV2 + 17*shat)*smT2**3) + &
!                 E**(6*y)*(-8*MV2**3*shat**4*(4*MV2 + shat) + &
!                    4*MV2**2*shat**3*(9*MV2 + 4*shat)*smT2 - &
!                    MV2*shat**2*(3*MV2 + 11*shat)*smT2**2 + &
!                    3*shat*(-5*MV2 + shat)*smT2**3 + 10*smT2**4))*LC(2)* &
!               SP(7))/ &
!             (E**(5*y)*shat*(E**y*shat - smT)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-(E**y*(MV2 + shat)*smT**5) + &
!                 E**(9*y)*(MV2 + shat)*smT**5 + &
!                 E**(8*y)*(6*MV2*shat - 5*smT2)*smT2**2 + &
!                 2*E**(6*y)*(5*MV2*shat - 4*smT2)*smT2**2 + smT2**3 - &
!                 E**(10*y)*smT2**3 + &
!                 E**(3*y)*(MV2 + shat)*smT*(-3*MV2*shat + smT2)* &
!                  (-2*MV2*shat + smT2) + &
!                 6*E**(5*y)*(MV2 + shat)*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**(7*y)*(MV2 + shat)*smT**3*(-7*MV2*shat + 5*smT2) + &
!                 E**(2*y)*smT2* &
!                  (-4*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) + &
!                 E**(4*y)*(-8*MV2**3*shat**3 + 4*MV2**2*shat**2*smT2 + &
!                    6*MV2*shat*smT2**2 - 4*smT2**3))*LC(8)*SP(7))/ &
!             (E**(4*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-(E**(11*y)*smT**7) - shat*smT2**3 + &
!                 E**(10*y)*shat*smT2**2*(MV2**2 + 2*smT2) + &
!                 E**y*smT**5*(shat**2 + 2*smT2) + &
!                 E**(9*y)*shat*smT**3* &
!                  (4*MV2**2*shat - (4*MV2 + shat)*smT2) + &
!                 E**(2*y)*smT2**2* &
!                  (MV2*shat*(MV2 + 6*shat) - 2*(MV2 + 3*shat)*smT2) - &
!                 E**(8*y)*MV2*smT2* &
!                  (2*MV2*shat**2*(7*MV2 + 5*shat) - &
!                    shat*(17*MV2 + 7*shat)*smT2 + 5*smT2**2) + &
!                 E**(3*y)*smT**3* &
!                  (MV2*(MV2 - 7*shat)*shat**2 + &
!                    shat*(-13*MV2 + 5*shat)*smT2 + 9*smT2**2) - &
!                 E**(4*y)*smT2* &
!                  (8*MV2**3*shat**2 - MV2*shat*(19*MV2 + 17*shat)*smT2 + &
!                    (9*MV2 + 13*shat)*smT2**2) + &
!                 2*E**(5*y)*smT* &
!                  (MV2**2*shat**3*(-5*MV2 + 3*shat) + &
!                    2*MV2*(5*MV2 - 3*shat)*shat**2*smT2 + &
!                    shat*(-13*MV2 + 3*shat)*smT2**2 + 7*smT2**3) + &
!                 E**(7*y)*smT* &
!                  (2*MV2**2*shat**3*(7*MV2 + 3*shat) - &
!                    MV2*shat**2*(MV2 + 5*shat)*smT2 + &
!                    shat*(-17*MV2 + shat)*smT2**2 + 8*smT2**3) - &
!                 2*E**(6*y)*(4*MV2**3*shat**3*(-2*MV2 + shat) + &
!                    MV2**2*shat**2*(19*MV2 + shat)*smT2 - &
!                    MV2*shat*(17*MV2 + 9*shat)*smT2**2 + &
!                    (6*MV2 + 5*shat)*smT2**3))*LC(12)*SP(7))/ &
!             (E**(6*y)*(-shat + E**y*smT)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
!         kap*SI(0)*((12*MV2*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(1))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (6*MV2*(smT2 + E**(2*y)*(-4*MV2*shat + 3*smT2))*LC(11))/ &
!             (E**(2*y)*(MV2*shat - smT2)* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (20*MV2*shat*smT*LC(15)*Sinh(y))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) + &
!            (8*(-smT**3 - E**(3*y)*MV2*(4*MV2*shat - 3*smT2) + &
!                 E**(2*y)*smT*(3*MV2*shat - 2*smT2) + &
!                 E**(4*y)*smT*(MV2*shat - smT2) + E**y*MV2*smT2)*LC(13)* &
!               SP(4))/ &
!             (E**(3*y)*(-shat + E**y*smT)*(MV2*shat - smT2)* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (8*LC(3)*(E**(6*y)*smT**5*SP(3) + smT**5*SP(4) - &
!                 E**(3*y)*(-2*MV2*shat + smT2)* &
!                  (2*shat*smT2 + MV2*(-shat**2 + smT2))*(SP(3) + SP(4)) &
!                  - E**y*smT2* &
!                  (shat*(-(MV2*shat) + smT2)*SP(3) + &
!                    (MV2 + shat)*smT2*SP(4)) - &
!                 E**(5*y)*smT2* &
!                  ((MV2 + shat)*smT2*SP(3) + &
!                    shat*(-(MV2*shat) + smT2)*SP(4)) + &
!                 E**(2*y)*smT* &
!                  ((-2*MV2**2*shat**2 + smT2**2)*SP(3) + &
!                    2*(-(MV2*shat) + smT2)**2*SP(4)) + &
!                 E**(4*y)*smT* &
!                  (2*(-(MV2*shat) + smT2)**2*SP(3) + &
!                    (-2*MV2**2*shat**2 + smT2**2)*SP(4))))/ &
!             (E**(2*y)*shat*(E**y*shat - smT)*(shat - E**y*smT)* &
!               (MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
!             + (4*LC(9)*(-2*E**(5*y)*smT**3*SP(3) + &
!                 2*E**(4*y)*MV2*smT2*SP(3) - MV2*smT2*SP(4) - &
!                 E**(3*y)*smT* &
!                  (4*smT2*SP(3) + MV2*shat*(-6*SP(3) + SP(4))) + &
!                 E**y*smT*(-2*smT2*SP(3) + MV2*shat*(2*SP(3) + SP(4))) + &
!                 E**(2*y)*MV2* &
!                  (-8*MV2*shat*SP(3) + smT2*(6*SP(3) + SP(4)))))/ &
!             (E**y*(E**y*shat - smT)*(MV2*shat - smT2)* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (8*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(10)*SP(6))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (8*(2*MV2*shat - (2*smT2*Cosh(y))/E**y)*LC(14)*SP(6))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) + &
!            (8*MV2*smT*LC(4)*Sinh(y)*SP(6))/ &
!             ((-(MV2*shat) + smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
!              + (8*(-2*MV2*shat + smT2 + E**(2*y)*smT2)*LC(8)*SP(7))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (8*(-2*MV2*shat + smT2 + smT2/E**(2*y))*LC(12)*SP(7))/ &
!             ((MV2*shat - smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) - &
!            (8*MV2*smT*LC(2)*Sinh(y)*SP(7))/ &
!             ((-(MV2*shat) + smT2)*(-2*MV2*shat + smT2 + smT2*Cosh(2*y))) &
!            ) + kap*SI(5)*(((4*smT2**2*(-(MV2*shat) + smT2)**2 - &
!                 3*E**(9*y)*shat*smT**5*(MV2*shat + 3*smT2) + &
!                 E**(7*y)*shat*smT**3* &
!                  (18*MV2**2*shat**2 + 47*MV2*shat*smT2 - 41*smT2**2) + &
!                 E**(8*y)*smT2**2* &
!                  (3*MV2**2*shat**2 + 17*MV2*shat*smT2 + 4*smT2**2) - &
!                 2*E**y*shat*smT**3* &
!                  (8*MV2**2*shat**2 - 21*MV2*shat*smT2 + 7*smT2**2) + &
!                 3*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
!                  (16*MV2**2*shat**2 - 32*MV2*shat*smT2 + 23*smT2**2) + &
!                 E**(2*y)*smT2* &
!                  (24*MV2**3*shat**3 - 31*MV2**2*shat**2*smT2 - &
!                    33*MV2*shat*smT2**2 + 16*smT2**3) + &
!                 E**(6*y)*smT2* &
!                  (-44*MV2**3*shat**3 - 61*MV2**2*shat**2*smT2 + &
!                    17*MV2*shat*smT2**2 + 16*smT2**3) - &
!                 E**(3*y)*shat*smT* &
!                  (8*MV2**3*shat**3 + 122*MV2**2*shat**2*smT2 - &
!                    157*MV2*shat*smT2**2 + 51*smT2**3) + &
!                 E**(4*y)*(-16*MV2**4*shat**4 + &
!                    148*MV2**3*shat**3*smT2 - &
!                    59*MV2**2*shat**2*smT2**2 - 25*MV2*shat*smT2**3 + &
!                    24*smT2**4))*LC(1))/ &
!             (2q0*E**(4*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            ((-4*E**(9*y)*smT2**2*(-(MV2*shat) + smT2)**2 + &
!                 shat*smT**5*(MV2*shat + 5*smT2) - &
!                 2*E**(2*y)*shat*smT**3* &
!                  (5*MV2**2*shat**2 + 13*MV2*shat*smT2 - 12*smT2**2) + &
!                 E**y*smT2**2* &
!                  (MV2**2*shat**2 - 9*MV2*shat*smT2 - 4*smT2**2) + &
!                 E**(8*y)*shat*smT**3* &
!                  (16*MV2**2*shat**2 - 31*MV2*shat*smT2 + 9*smT2**2) + &
!                 6*E**(4*y)*shat*smT* &
!                  (4*MV2**3*shat**3 - 11*MV2*shat*smT2**2 + 7*smT2**3) + &
!                 2*E**(6*y)*shat*smT* &
!                  (4*MV2**3*shat**3 + 21*MV2**2*shat**2*smT2 - &
!                    35*MV2*shat*smT2**2 + 16*smT2**3) - &
!                 E**(7*y)*smT2* &
!                  (24*MV2**3*shat**3 - 27*MV2**2*shat**2*smT2 - &
!                    25*MV2*shat*smT2**2 + 16*smT2**3) - &
!                 E**(3*y)*smT2* &
!                  (8*MV2**3*shat**3 - 69*MV2**2*shat**2*smT2 + &
!                    9*MV2*shat*smT2**2 + 16*smT2**3) + &
!                 E**(5*y)*(16*MV2**4*shat**4 - 96*MV2**3*shat**3*smT2 + &
!                    51*MV2**2*shat**2*smT2**2 + 17*MV2*shat*smT2**3 - &
!                    24*smT2**4))*LC(11))/ &
!             (2q0*E**(5*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (shat*LC(15)*(-((MV2*shat - smT2)* &
!                    (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)) + &
!                 smT*(4*smT*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
!                     Cosh(2*y) + &
!                    (-(MV2*shat*smT**3) + smT**5)*Cosh(4*y) + &
!                    4*MV2*(18*MV2**2*shat**2 + 25*MV2*shat*smT2 - &
!                       13*smT2**2)*Sinh(y) + &
!                    10*smT*(-8*MV2**2*shat**2 + MV2*shat*smT2 + &
!                       3*smT2**2)*Sinh(2*y) + &
!                    4*MV2*(3*MV2*shat - 13*smT2)*smT2*Sinh(3*y) + &
!                    5*smT**3*(MV2*shat + 3*smT2)*Sinh(4*y))))/ &
!             (2q0*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*LC(3)*(-(E**(10*y)*smT**5*(2*MV2*shat - smT2)*SP(3)) + &
!                 2*E**(9*y)*MV2**2*shat*smT2**2*SP(3) - &
!                 MV2*shat*smT**5*SP(4) + &
!                 2*E**y*MV2**2*shat*smT2**2*SP(4) + &
!                 2*E**(5*y)*MV2**2*shat* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + 2*smT2**2)* &
!                  (SP(3) + SP(4)) - &
!                 2*E**(3*y)*MV2**2*shat*smT2* &
!                  (smT2*(SP(3) - 4*SP(4)) + 7*MV2*shat*SP(4)) - &
!                 E**(4*y)*smT* &
!                  (2*(5*MV2**3*shat**3 - 9*MV2**2*shat**2*smT2 + &
!                       7*MV2*shat*smT2**2 - 2*smT2**3)*SP(3) + &
!                    MV2*shat*(2*MV2*shat - 5*smT2)*(MV2*shat - smT2)* &
!                     SP(4)) - &
!                 E**(6*y)*smT* &
!                  (3*(3*MV2*shat - 2*smT2)*(MV2*shat - smT2)* &
!                     (2*MV2*shat - smT2)*SP(3) + &
!                    2*MV2*shat* &
!                     (5*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)*SP(4) &
!                    ) - 2*E**(7*y)*MV2**2*shat*smT2* &
!                  (7*MV2*shat*SP(3) + smT2*(-4*SP(3) + SP(4))) + &
!                 E**(8*y)*smT**3* &
!                  (-16*MV2*shat*smT2*SP(3) + 4*smT2**2*SP(3) + &
!                    MV2**2*shat**2*(14*SP(3) + SP(4))) + &
!                 E**(2*y)*smT**3* &
!                  (smT2**2*SP(3) - MV2*shat*smT2*(SP(3) + 4*SP(4)) + &
!                    MV2**2*shat**2*(SP(3) + 6*SP(4)))))/ &
!             (E**(5*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (LC(9)*(4*E**(10*y)*shat*smT2**3*SP(3) - &
!                 4*E**(9*y)*smT**5*(MV2*shat + smT2)*SP(3) + &
!                 shat*smT2**2*(-3*MV2*shat + smT2)*SP(4) - &
!                 E**(8*y)*shat*smT2**2* &
!                  (4*(7*MV2*shat - 5*smT2)*SP(3) - &
!                    (MV2*shat + smT2)*SP(4)) + &
!                 2*E**(6*y)*shat*smT2* &
!                  (18*(-(MV2*shat) + smT2)**2*SP(3) - &
!                    (3*MV2*shat - 2*smT2)*(MV2*shat + smT2)*SP(4)) - &
!                 E**y*smT**3*(4*(-(MV2*shat) + smT2)**2*SP(3) - &
!                    MV2*shat*(MV2*shat + 3*smT2)*SP(4)) + &
!                 E**(7*y)*smT**3* &
!                  (4*(7*MV2*shat - 4*smT2)*(MV2*shat + smT2)*SP(3) - &
!                    MV2*shat*(MV2*shat + 3*smT2)*SP(4)) - &
!                 2*E**(4*y)*shat* &
!                  (-2*smT2*(12*MV2**2*shat**2 - 17*MV2*shat*smT2 + &
!                       7*smT2**2)*SP(3) + &
!                    (MV2*shat - smT2)* &
!                     (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(4)) - &
!                 E**(5*y)*smT* &
!                  (4*(MV2*shat + 2*smT2)* &
!                     (6*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(3) + MV2*shat* &
!                     (-16*MV2**2*shat**2 + MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(4)) + &
!                 E**(3*y)*(-16*smT**7*SP(3) + &
!                    MV2**2*shat**2*smT**3*(-4*SP(3) + SP(4)) + &
!                    3*MV2*shat*smT**5*(12*SP(3) + SP(4)) - &
!                    8*MV2**3*shat**3*smT*(SP(3) + 2*SP(4))) + &
!                 2*E**(2*y)*shat*smT2* &
!                  (2*smT2**2*(2*SP(3) + SP(4)) + &
!                    MV2**2*shat**2*(6*SP(3) + 11*SP(4)) - &
!                    MV2*shat*smT2*(12*SP(3) + 11*SP(4)))))/ &
!             (E**(4*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (LC(13)*(E**(10*y)*shat*(MV2*shat - smT2)*smT2**2*SP(3) + &
!                 4*E**(2*y)*shat*(7*MV2*shat - 5*smT2)*smT2**2*SP(4) - &
!                 4*shat*smT2**3*SP(4) + &
!                 4*E**y*smT**5*(MV2*shat + smT2)*SP(4) - &
!                 E**(9*y)*smT**3*(MV2*shat - smT2)* &
!                  (MV2*shat*(SP(3) - 4*SP(4)) + 4*smT2*SP(4)) + &
!                 E**(3*y)*smT**3* &
!                  (MV2*shat*(MV2*shat - smT2)*SP(3) - &
!                    4*(7*MV2*shat - 4*smT2)*(MV2*shat + smT2)*SP(4)) + &
!                 E**(7*y)*smT* &
!                  (MV2*shat*(MV2*shat - smT2)*(4*MV2*shat - smT2)* &
!                     SP(3) + 4*(2*MV2*shat - smT2)* &
!                     (MV2**2*shat**2 + MV2*shat*smT2 - 4*smT2**2)*SP(4)) &
!                  - E**(5*y)*smT* &
!                  (MV2*shat*(MV2*shat - smT2)*(4*MV2*shat - smT2)* &
!                     SP(3) - 4*(MV2*shat + 2*smT2)* &
!                     (6*MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(4)) + &
!                 E**(6*y)*shat* &
!                  ((MV2*shat - smT2)* &
!                     (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
!                     SP(3) - 4*smT2* &
!                     (12*MV2**2*shat**2 - 17*MV2*shat*smT2 + 7*smT2**2)* &
!                     SP(4)) - &
!                 E**(8*y)*shat*smT2* &
!                  (6*MV2**2*shat**2*(SP(3) + 2*SP(4)) - &
!                    3*MV2*shat*smT2*(3*SP(3) + 8*SP(4)) + &
!                    smT2**2*(3*SP(3) + 8*SP(4))) - &
!                 E**(4*y)*shat*(MV2*shat - smT2)*smT2* &
!                  (2*MV2*shat*(SP(3) + 18*SP(4)) - &
!                    smT2*(SP(3) + 36*SP(4)))))/ &
!             (E**(6*y)*shat*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 - &
!                 E**(6*y)*shat*smT**3*(7*MV2*shat - 5*smT2) + &
!                 6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 - &
!                 E**(7*y)*smT2**2*(MV2*shat + smT2) + &
!                 E**(5*y)*smT2* &
!                  (8*MV2**2*shat**2 + MV2*shat*smT2 - 3*smT2**2) - &
!                 E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
!                    smT2**2) + &
!                 E**(2*y)*shat*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
!                 E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(10)*SP(6))/ &
!             (E**(3*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 + &
!                 E**(2*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
!                 6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**y*smT2**2*(MV2*shat + smT2) + &
!                 E**(7*y)*smT2* &
!                  (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
!                 E**(6*y)*shat*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(3*y)*smT2* &
!                  (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
!                 E**(5*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(14)*SP(6))/ &
!             (E**(5*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (2*LC(4)*(-(MV2*shat) + smT2 + &
!                 (12*MV2*shat*smT*(MV2*shat - smT2)*(MV2 - smT*Cosh(y))* &
!                    Sinh(y))/(-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2 - &
!                 (2*smT*(2*MV2*(MV2*shat + smT2) + &
!                      (-5*MV2*shat*smT + smT**3)*Cosh(y))*Sinh(y))/ &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y)))*SP(6))/ &
!             (-(MV2*shat) + smT2)**2 + &
!            (4*(shat*smT**5 - E**(8*y)*shat*smT**5 + &
!                 E**(6*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
!                 6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**(7*y)*smT2**2*(MV2*shat + smT2) + &
!                 E**y*smT2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + &
!                    smT2**2) - &
!                 E**(2*y)*shat*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(5*y)*smT2* &
!                  (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
!                 E**(3*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(8)*SP(7))/ &
!             (E**(3*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
!            (4*(-(shat*smT**5) + E**(8*y)*shat*smT**5 + &
!                 E**(2*y)*shat*smT**3*(7*MV2*shat - 5*smT2) - &
!                 6*E**(4*y)*shat*smT*(-(MV2*shat) + smT2)**2 + &
!                 E**y*smT2**2*(MV2*shat + smT2) + &
!                 E**(7*y)*smT2* &
!                  (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) - &
!                 E**(6*y)*shat*smT* &
!                  (6*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2) + &
!                 E**(3*y)*smT2* &
!                  (-8*MV2**2*shat**2 - MV2*shat*smT2 + 3*smT2**2) + &
!                 E**(5*y)*(4*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                    7*MV2*shat*smT2**2 + 3*smT2**3))*LC(12)*SP(7))/ &
!             (E**(5*y)*(-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!            (LC(2)*((MV2*shat - smT2)* &
!                  (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2) + &
!                 smT*(-4*smT*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
!                     Cosh(2*y) + smT**3*(MV2*shat - smT2)*Cosh(4*y) + &
!                    4*MV2*(-2*MV2*shat + smT2)*(5*MV2*shat + smT2)* &
!                     Sinh(y) + &
!                    2*smT*(16*MV2**2*shat**2 - 13*MV2*shat*smT2 + &
!                       smT2**2)*Sinh(2*y) + &
!                    4*MV2*smT2*(MV2*shat + smT2)*Sinh(3*y) + &
!                    smT**3*(-5*MV2*shat + smT2)*Sinh(4*y)))*SP(7))/ &
!             ((-(MV2*shat) + smT2)**2* &
!               (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
!         SI(13)*(kap*(((-12*MV2*shat**2*smT2**4 + &
!                    E**(14*y)*shat*smT2**4*(5*MV2*shat + 7*smT2) - &
!                    2*E**(13*y)*smT**7* &
!                     (2*MV2**2*shat**2 + 9*MV2*shat*smT2 + smT2**2) + &
!                    E**(12*y)*smT2**3* &
!                     (MV2**2*(-8*MT2 + 3*MV2 - 51*shat)*shat**2 + &
!                       MV2*(-8*MT2 + 17*MV2 - 47*shat)*shat*smT2 + &
!                       2*(8*MT2 + 2*MV2 + 19*shat)*smT2**2) + &
!                    E**(11*y)*smT**5* &
!                     (42*MV2**3*shat**3 + 143*MV2**2*shat**2*smT2 - &
!                       76*MV2*shat*smT2**2 - 13*smT2**3) + &
!                    E**y*smT**5* &
!                     (4*MV2**3*shat**3 - MV2**2*shat**2*smT2 + &
!                       24*MV2*shat*smT2**2 - 3*smT2**3) + &
!                    E**(2*y)*smT2**2* &
!                     (-8*MV2**3*shat**3*(4*MT2 + MV2 + shat) + &
!                       MV2**2*shat**2*(104*MT2 + 21*MV2 + 145*shat)* &
!                        smT2 - &
!                       MV2*shat*(88*MT2 + 41*(MV2 + 2*shat))*smT2**2 + &
!                       (16*MT2 + 4*MV2 + 5*shat)*smT2**3) + &
!                    E**(10*y)*smT2**2* &
!                     (2*MV2**3*shat**3*(40*MT2 - 15*MV2 + 79*shat) + &
!                       MV2**2*shat**2*(40*MT2 - 155*MV2 + 81*shat)* &
!                        smT2 + &
!                       MV2*(-200*MT2 + 45*MV2 - 194*shat)*shat*smT2**2 + &
!                       (80*MT2 + 20*MV2 + 87*shat)*smT2**3) + &
!                    2*E**(5*y)*smT**3* &
!                     (292*MV2**4*shat**4 - 356*MV2**3*shat**3*smT2 + &
!                       101*MV2**2*shat**2*smT2**2 + &
!                       43*MV2*shat*smT2**3 - 20*smT2**4) + &
!                    E**(3*y)*smT**3* &
!                     (56*MV2**4*shat**4 - 254*MV2**3*shat**3*smT2 + &
!                       35*MV2**2*shat**2*smT2**2 + 84*MV2*shat*smT2**3 - &
!                       17*smT2**4) - &
!                    2*E**(7*y)*smT*(-(MV2*shat) + smT2)* &
!                     (-48*MV2**4*shat**4 + 108*MV2**3*shat**3*smT2 - &
!                       162*MV2**2*shat**2*smT2**2 + &
!                       37*MV2*shat*smT2**3 + 25*smT2**4) - &
!                    E**(9*y)*smT**3* &
!                     (152*MV2**4*shat**4 + 204*MV2**3*shat**3*smT2 - &
!                       379*MV2**2*shat**2*smT2**2 + &
!                       108*MV2*shat*smT2**3 + 35*smT2**4) + &
!                    E**(4*y)*smT2* &
!                     (-16*MV2**4*shat**4*(-10*MT2 + 3*(MV2 + shat)) + &
!                       2*MV2**3*(-360*MT2 + 67*MV2 - 139*shat)*shat**3* &
!                        smT2 + &
!                       MV2**2*shat**2*(920*MT2 + 131*MV2 + 357*shat)* &
!                        smT2**2 - &
!                       MV2*shat*(440*MT2 + 117*MV2 + 193*shat)*smT2**3 + &
!                       10*(8*MT2 + 2*MV2 + 3*shat)*smT2**4) + &
!                    2*E**(8*y)*smT2* &
!                     (4*MV2**4*(-20*MT2 + 17*MV2 - 9*shat)*shat**4 + &
!                       MV2**3*(-120*MT2 + 141*MV2 - 101*shat)*shat**3* &
!                        smT2 + &
!                       MV2**2*shat**2*(440*MT2 - 111*MV2 + 139*shat)* &
!                        smT2**2 + &
!                       2*MV2*(-160*MT2 + MV2 - 73*shat)*shat*smT2**3 + &
!                       2*(40*MT2 + 10*MV2 + 27*shat)*smT2**4) + &
!                    E**(6*y)*(-32*MV2**11 + &
!                       32*MV2**5* &
!                        (MV2**6 + (-4*MT2 + MV2)*shat**5 + shat**6) + &
!                       8*MV2**4*shat**4*(128*MT2 - 47*MV2 + 7*shat)* &
!                        smT2 + &
!                       2*MV2**3*(-952*MT2 + 83*MV2 - 11*shat)*shat**3* &
!                        smT2**2 - &
!                       2*MV2**2*(-824*MT2 + MV2 - 155*shat)*shat**2* &
!                        smT2**3 - &
!                       MV2*shat*(100*(8*MT2 + MV2) + 273*shat)*smT2**4 + &
!                       (40*(4*MT2 + MV2) + 77*shat)*smT2**5))*LC(1))/ &
!                (4q0*E**(6*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((-6*shat*smT2**5 - &
!                    E**(14*y)*shat*(MV2*shat - 3*smT2)*smT2**3* &
!                     (2*MV2*shat + smT2) + &
!                    2*E**y*smT**9*(5*MV2*shat + smT2) + &
!                    E**(13*y)*smT**7* &
!                     (5*MV2**2*shat**2 - 19*MV2*shat*smT2 + 2*smT2**2) - &
!                    E**(2*y)*smT2**3* &
!                     (-(MV2**2*shat**2*(-4*MT2 + MV2 + shat)) + &
!                       MV2*(-20*MT2 + 9*MV2 - 59*shat)*shat*smT2 + &
!                       2*(8*MT2 + 2*MV2 + 15*shat)*smT2**2) + &
!                    E**(3*y)*smT**5* &
!                     (-4*MV2**3*shat**3 - 93*MV2**2*shat**2*smT2 + &
!                       37*MV2*shat*smT2**2 + 12*smT2**3) + &
!                    E**(12*y)*smT2**2* &
!                     (8*MV2**3*shat**3*(4*MT2 + MV2 + 3*shat) - &
!                       MV2**2*shat**2*(92*MT2 + 25*MV2 + 73*shat)*smT2 + &
!                       MV2*shat*(76*MT2 + 33*MV2 + 8*shat)*smT2**2 + &
!                       (-4*(4*MT2 + MV2) + 11*shat)*smT2**3) - &
!                    E**(4*y)*smT2**2* &
!                     (2*MV2**3*shat**3*(-20*MT2 + 5*MV2 + 6*shat) + &
!                       MV2**2*shat**2*(220*MT2 - 95*MV2 + 157*shat)* &
!                        smT2 + &
!                       MV2*(-260*MT2 + 5*MV2 - 164*shat)*shat*smT2**2 + &
!                       (80*MT2 + 20*MV2 + 61*shat)*smT2**3) + &
!                    E**(11*y)*smT**3* &
!                     (-88*MV2**4*shat**4 + 210*MV2**3*shat**3*smT2 - &
!                       37*MV2**2*shat**2*smT2**2 - 49*MV2*shat*smT2**3 + &
!                       12*smT2**4) + &
!                    2*E**(7*y)*smT*(-(MV2*shat) + smT2)* &
!                     (40*MV2**4*shat**4 + 28*MV2**3*shat**3*smT2 - &
!                       121*MV2**2*shat**2*smT2**2 + &
!                       38*MV2*shat*smT2**3 + 20*smT2**4) + &
!                    E**(5*y)*smT**3* &
!                     (36*MV2**4*shat**4 + 222*MV2**3*shat**3*smT2 - &
!                       287*MV2**2*shat**2*smT2**2 + &
!                       59*MV2*shat*smT2**3 + 30*smT2**4) + &
!                    E**(10*y)*smT2* &
!                     (16*MV2**4*shat**4*(-10*MT2 + 3*MV2 + shat) + &
!                       2*MV2**3*shat**3*(300*MT2 - 47*MV2 + 10*shat)* &
!                        smT2 + &
!                       MV2**2*shat**2*(-740*MT2 - 71*MV2 + 49*shat)* &
!                        smT2**2 + &
!                       MV2*(380*MT2 + 77*MV2 - 27*shat)*shat*smT2**3 - &
!                       4*(20*MT2 + 5*MV2 - 2*shat)*smT2**4) - &
!                    E**(6*y)*smT2* &
!                     (-8*MV2**4*shat**4*(-16*MT2 + 4*MV2 + 5*shat) + &
!                       2*MV2**3*shat**3*(-372*MT2 + 153*MV2 + 2*shat)* &
!                        smT2 + &
!                       2*MV2**2*shat**2*(572*MT2 - 83*MV2 + 4*shat)* &
!                        smT2**2 - &
!                       MV2*shat*(688*MT2 + 28*MV2 + 123*shat)*smT2**3 + &
!                       (40*(4*MT2 + MV2) + 61*shat)*smT2**4) + &
!                    2*E**(9*y)*smT* &
!                     (32*MV2**5*shat**5 - 178*MV2**4*shat**4*smT2 + &
!                       197*MV2**3*shat**3*smT2**2 - &
!                       83*MV2**2*shat**2*smT2**3 - 13*MV2*shat*smT2**4 + &
!                       15*smT2**5) - &
!                    2*E**(8*y)* &
!                     (16*MV2**5*shat**5*(-4*MT2 + MV2 + shat) - &
!                       4*MV2**4*shat**4*(-92*MT2 + 26*MV2 + 25*shat)* &
!                        smT2 + &
!                       MV2**3*shat**3*(-700*MT2 + 71*MV2 + 222*shat)* &
!                        smT2**2 + &
!                       MV2**2*(692*MT2 - 29*MV2 - 111*shat)*shat**2* &
!                        smT2**3 - &
!                       2*MV2*(188*MT2 + 17*MV2 - 3*shat)*shat*smT2**4 + &
!                       4*(20*MT2 + 5*MV2 + 3*shat)*smT2**5))*LC(11))/ &
!                (4q0*E**(8*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
!               (LC(15)*(2*(MV2*shat - smT2)* &
!                     (8*MV2**4*shat**4 - 32*MV2**3*shat**3*smT2 + &
!                       35*MV2**2*shat**2*smT2**2 - 21*MV2*shat*smT2**3 + &
!                       5*smT2**4) + &
!                    smT*(2*shat*(-(MV2*shat) + smT2)* &
!                        (-32*MV2**3*shat**3 + 64*MV2**2*shat**2*smT2 - &
!                          52*MV2*shat*smT2**2 + 15*smT2**3)*Cosh(y) - &
!                       smT*(-(MV2*shat) + smT2)* &
!                        (-16*MV2**3*shat**3 + 72*MV2**2*shat**2*smT2 - &
!                          56*MV2*shat*smT2**2 + 15*smT2**3)*Cosh(2*y) + &
!                       2*shat*smT2*(-(MV2*shat) + smT2)* &
!                        (32*MV2**2*shat**2 - 34*MV2*shat*smT2 + &
!                          11*smT2**2)*Cosh(3*y) + &
!                       2*smT**3*(MV2*shat - smT2)* &
!                        (MV2**2*shat**2 - 7*MV2*shat*smT2 + 3*smT2**2)* &
!                        Cosh(4*y) + &
!                       10*shat*smT2**2*(-2*MV2*shat + smT2)* &
!                        (-(MV2*shat) + smT2)*Cosh(5*y) + &
!                       smT**7*(MV2*shat - smT2)*Cosh(6*y) + &
!                       2*shat*smT2**3*(-(MV2*shat) + smT2)*Cosh(7*y) + &
!                       2*(8*MV2**4*(-28*MT2 + 13*MV2 - 2*shat)*shat**4 + &
!                          2*MV2**3*(100*MT2 + 86*MV2 - 67*shat)*shat**3* &
!                           smT2 + &
!                          4*MV2**2*shat**2*(26*MT2 - 25*MV2 + 12*shat)* &
!                           smT2**2 + &
!                          MV2*(-112*MT2 + 8*MV2 - 85*shat)*shat* &
!                           smT2**3 + (32*MT2 + 16*MV2 + 37*shat)*smT2**4) &
!                         *Sinh(y) - &
!                       smT*(328*MV2**4*shat**4 + &
!                          196*MV2**3*shat**3*smT2 - &
!                          546*MV2**2*shat**2*smT2**2 + &
!                          187*MV2*shat*smT2**3 + 35*smT2**4)*Sinh(2*y) + &
!                       4*smT2* &
!                        (MV2**3*shat**3*(36*MT2 - 2*MV2 + 53*shat) + &
!                          MV2**2*shat**2*(2*MT2 - 69*MV2 + 41*shat)* &
!                           smT2 + &
!                          MV2*(-62*MT2 + 9*(MV2 - 8*shat))*shat* &
!                           smT2**2 + 3*(8*MT2 + 4*MV2 + 11*shat)*smT2**3) &
!                         *Sinh(3*y) + &
!                       2*smT**3* &
!                        (23*MV2**3*shat**3 + 138*MV2**2*shat**2*smT2 - &
!                          67*MV2*shat*smT2**2 - 14*smT2**3)*Sinh(4*y) + &
!                       4*smT2**2* &
!                        (MV2**2*(-2*MT2 + MV2 - 15*shat)*shat**2 + &
!                          MV2*(-6*MT2 + 5*MV2 - 28*shat)*shat*smT2 + &
!                          2*(4*MT2 + 2*MV2 + 9*shat)*smT2**2)*Sinh(5*y) &
!                        - smT**5* &
!                        (6*MV2**2*shat**2 + 27*MV2*shat*smT2 + &
!                          7*smT2**2)*Sinh(6*y) + &
!                       2*shat*smT2**3*(3*MV2*shat + 7*smT2)*Sinh(7*y))))/ &
!                (4q0*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) - &
!               (2*(-2*E**(3*y)*MV2*shat**2*smT + E**y*shat*smT**3 + &
!                    E**(5*y)*shat*smT**3 - MV2*shat*smT2 + &
!                    E**(4*y)*smT2*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*(4*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
!                       smT2**2))*LC(7)*SP(1))/ &
!                (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
!               (8*smT*(shat + smT*Cosh(y))*LC(6)*Sinh(y)*SP(2))/ &
!                (shat*(MV2*shat - smT2)) + &
!               (LC(13)*(-(E**(15*y)*smT2**3*(-(MV2*shat) + smT2)* &
!                       (-(MV2*shat) + 2*smT2)*SP(3)) + &
!                    4*shat*smT**9*SP(4) - &
!                    2*E**y*smT2**4*(3*MV2*shat + smT2)*SP(4) - &
!                    4*E**(2*y)*smT**7* &
!                     (MV2*shat*(2*MT2 - MV2 + 11*shat) - &
!                       (2*MT2 + MV2 + 6*shat)*smT2)*SP(4) + &
!                    E**(14*y)*shat*smT**5* &
!                     (-((MV2*shat - smT2)*(MV2**2 - MV2*shat + 2*smT2)* &
!                          SP(3)) - 4*MV2*shat*smT2*SP(4)) - &
!                    E**(3*y)*smT2**3* &
!                     (MV2*shat*(MV2*shat - smT2)*SP(3) - &
!                       4*(5*MV2*shat - 3*smT2)*(3*MV2*shat + smT2)*SP(4)) &
!                      - E**(7*y)*smT2*(-(MV2*shat) + smT2)* &
!                     ((-56*MV2**3*shat**3 + 92*MV2**2*shat**2*smT2 - &
!                          51*MV2*shat*smT2**2 + 10*smT2**3)*SP(3) + &
!                       8*(-(MV2*shat) + smT2)**2*(12*MV2*shat + 5*smT2)* &
!                        SP(4)) - &
!                    E**(11*y)*smT2* &
!                     (5*(MV2*shat - smT2)* &
!                        (8*MV2**3*shat**3 - 20*MV2**2*shat**2*smT2 + &
!                          15*MV2*shat*smT2**2 - 4*smT2**3)*SP(3) - &
!                       4*(2*MV2*shat - smT2)*(MV2*shat + smT2)* &
!                        (4*MV2**2*shat**2 - 11*MV2*shat*smT2 + &
!                          3*smT2**2)*SP(4)) - &
!                    2*E**(5*y)*smT2**2* &
!                     (-((MV2*shat - smT2)* &
!                          (7*MV2**2*shat**2 - 6*MV2*shat*smT2 + smT2**2)* &
!                          SP(3)) + &
!                       5*(3*MV2*shat + smT2)* &
!                        (6*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
!                        SP(4)) + &
!                    E**(4*y)*smT**5* &
!                     (-(shat*(MV2*(-MV2 + shat) - 2*smT2)* &
!                          (MV2*shat - smT2)*SP(3)) + &
!                       4*(-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
!                          MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
!                          5*(2*MT2 + MV2 + 3*shat)*smT2**2)*SP(4)) + &
!                    E**(12*y)*smT**3* &
!                     (-(shat*(MV2*shat - smT2)* &
!                          (2*MV2**2*shat*(-5*MV2 + 3*shat) + &
!                            MV2*(5*MV2 - 17*shat)*smT2 + 6*smT2**2)*SP(3) &
!                          ) + &
!                       4*(-2*MV2**3*shat**3*(2*MT2 + MV2 + shat) + &
!                          2*MV2**2*shat**2*(6*MT2 + 3*MV2 + 8*shat)* &
!                           smT2 - &
!                          MV2*shat*(10*MT2 + 7*MV2 + 10*shat)*smT2**2 + &
!                          (2*MT2 + MV2 + shat)*smT2**3)*SP(4)) + &
!                    E**(6*y)*smT**3* &
!                     (shat*(MV2*shat - smT2)* &
!                        (2*MV2**2*shat*(-5*MV2 + 3*shat) + &
!                          MV2*(5*MV2 - 17*shat)*smT2 + 6*smT2**2)*SP(3) &
!                        - 20*(2*MV2**3*shat**3* &
!                           (6*MT2 - 3*MV2 + 5*shat) - &
!                          8*MV2**2*shat**2*(3*MT2 + 2*shat)*smT2 + &
!                          MV2*shat*(16*MT2 + 4*MV2 + 13*shat)*smT2**2 - &
!                          2*(2*MT2 + MV2 + 2*shat)*smT2**3)*SP(4)) + &
!                    2*E**(9*y)* &
!                     (2*(2*MV2*shat - smT2)*(-(MV2*shat) + smT2)**2* &
!                        (8*MV2**2*shat**2 - 8*MV2*shat*smT2 + 5*smT2**2)* &
!                        SP(3) + &
!                       smT2*(64*MV2**4*shat**4 - &
!                          72*MV2**3*shat**3*smT2 + &
!                          12*MV2**2*shat**2*smT2**2 + &
!                          31*MV2*shat*smT2**3 - 15*smT2**4)*SP(4)) + &
!                    2*E**(8*y)*smT* &
!                     (-(shat*(MV2*shat - smT2)* &
!                          (4*MV2**3*shat**2*(-3*MV2 + shat) + &
!                            3*MV2**2*(3*MV2 - 5*shat)*shat*smT2 - &
!                            2*MV2*(MV2 - 4*shat)*smT2**2 - 2*smT2**3)* &
!                          SP(3)) + &
!                       2*(4*MV2**4*shat**4*(12*MT2 - 3*MV2 + 5*shat) - &
!                          2*MV2**3*shat**3*(62*MT2 + 7*MV2 + 15*shat)* &
!                           smT2 + &
!                          12*MV2**2*shat**2*(12*MT2 + 2*MV2 + 5*shat)* &
!                           smT2**2 - &
!                          2*MV2*shat*(44*MT2 + 14*MV2 + 25*shat)* &
!                           smT2**3 + 5*(4*MT2 + 2*MV2 + 3*shat)*smT2**4)* &
!                        SP(4)) - &
!                    2*E**(10*y)*smT* &
!                     (-(shat*(MV2*shat - smT2)* &
!                          (4*MV2**3*shat**2*(-3*MV2 + shat) + &
!                            3*MV2**2*(3*MV2 - 5*shat)*shat*smT2 - &
!                            2*MV2*(MV2 - 4*shat)*smT2**2 - 2*smT2**3)* &
!                          SP(3)) - &
!                       2*(-4*MV2**4*shat**4*(-4*MT2 + MV2 + shat) + &
!                          2*MV2**3*(-34*MT2 + MV2 - 15*shat)*shat**3* &
!                           smT2 + &
!                          4*MV2**2*shat**2*(22*MT2 + 7*MV2 + 11*shat)* &
!                           smT2**2 - &
!                          MV2*shat*(46*MT2 + 21*MV2 + 27*shat)*smT2**3 + &
!                          (10*MT2 + 5*MV2 + 6*shat)*smT2**4)*SP(4)) - &
!                    2*E**(13*y)*smT2**2* &
!                     (smT2**3*(5*SP(3) + SP(4)) - &
!                       MV2**3*shat**3*(5*SP(3) + 2*SP(4)) + &
!                       MV2**2*shat**2*smT2*(17*SP(3) + 6*SP(4)) - &
!                       MV2*shat*smT2**2*(17*SP(3) + 9*SP(4)))))/ &
!                (2q0*E**(9*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               (LC(9)*(-4*E**(16*y)*shat*smT**9*SP(3) + &
!                    2*E**(15*y)*smT2**4*(3*MV2*shat + smT2)*SP(3) + &
!                    2*shat*smT**7*(2*MV2*shat - smT2)*SP(4) + &
!                    E**y*smT2**3* &
!                     (-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2)* &
!                     SP(4) + 2*E**(14*y)*smT**7* &
!                     (2*(MV2*shat*(2*MT2 - MV2 + 11*shat) - &
!                          (2*MT2 + MV2 + 6*shat)*smT2)*SP(3) - &
!                       MV2*shat**2*SP(4)) + &
!                    E**(2*y)*shat*smT**5* &
!                     (4*MV2*shat*smT2*SP(3) + &
!                       (MV2**2*(-4*MT2 + MV2 - 43*shat)*shat + &
!                          MV2*(4*MT2 + 3*MV2 + 44*shat)*smT2 - 11*smT2**2 &
!                          )*SP(4)) - &
!                    E**(12*y)*smT**5* &
!                     (4*(-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
!                          MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
!                          5*(2*MT2 + MV2 + 3*shat)*smT2**2)*SP(3) + &
!                       shat*(MV2**2*(-4*MT2 + MV2 - 23*shat)*shat + &
!                          MV2*(4*MT2 + 3*MV2 + 14*shat)*smT2 - smT2**2)* &
!                        SP(4)) + &
!                    2*E**(13*y)*smT2**3* &
!                     (-2*(5*MV2*shat - 3*smT2)*(3*MV2*shat + smT2)* &
!                        SP(3) + &
!                       (2*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) &
!                       ) + E**(4*y)*smT**3* &
!                     (4*(2*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
!                          2*MV2**2*shat**2*(6*MT2 + 3*MV2 + 8*shat)* &
!                           smT2 + &
!                          MV2*shat*(10*MT2 + 7*MV2 + 10*shat)*smT2**2 - &
!                          (2*MT2 + MV2 + shat)*smT2**3)*SP(3) + &
!                       shat*(2*MV2**3*shat**2* &
!                           (20*MT2 - 5*MV2 + 74*shat) - &
!                          MV2**2*shat*(60*MT2 + 25*MV2 + 229*shat)* &
!                           smT2 + &
!                          MV2*(20*MT2 + 15*MV2 + 126*shat)*smT2**2 - &
!                          23*smT2**3)*SP(4)) + &
!                    E**(10*y)*smT**3* &
!                     (20*(2*MV2**3*shat**3*(6*MT2 - 3*MV2 + 5*shat) - &
!                          8*MV2**2*shat**2*(3*MT2 + 2*shat)*smT2 + &
!                          MV2*shat*(16*MT2 + 4*MV2 + 13*shat)*smT2**2 - &
!                          2*(2*MT2 + MV2 + 2*shat)*smT2**3)*SP(3) + &
!                       shat*(2*MV2**3*(-20*MT2 + 5*MV2 - 42*shat)* &
!                           shat**2 + &
!                          MV2**2*shat*(60*MT2 + 25*MV2 + 97*shat)*smT2 - &
!                          MV2*(20*MT2 + 15*MV2 + 36*shat)*smT2**2 + &
!                          smT2**3)*SP(4)) + &
!                    E**(11*y)*smT2**2* &
!                     (10*(3*MV2*shat + smT2)* &
!                        (6*MV2**2*shat**2 - 8*MV2*shat*smT2 + 3*smT2**2)* &
!                        SP(3) + &
!                       (-42*MV2**3*shat**3 + 48*MV2**2*shat**2*smT2 - &
!                          33*MV2*shat*smT2**2 + 11*smT2**3)*SP(4)) - &
!                    2*E**(8*y)*smT* &
!                     (2*(4*MV2**4*shat**4*(12*MT2 - 3*MV2 + 5*shat) - &
!                          2*MV2**3*shat**3*(62*MT2 + 7*MV2 + 15*shat)* &
!                           smT2 + &
!                          12*MV2**2*shat**2*(12*MT2 + 2*MV2 + 5*shat)* &
!                           smT2**2 - &
!                          2*MV2*shat*(44*MT2 + 14*MV2 + 25*shat)* &
!                           smT2**3 + 5*(4*MT2 + 2*MV2 + 3*shat)*smT2**4)* &
!                        SP(3) + &
!                       shat*(12*MV2**4*(-4*MT2 + 2*MV2 - 3*shat)* &
!                           shat**3 + &
!                          MV2**3*shat**2*(84*MT2 + 3*MV2 + 34*shat)* &
!                           smT2 - &
!                          MV2**2*shat*(44*MT2 + 13*MV2 + 11*shat)* &
!                           smT2**2 + &
!                          2*MV2*(4*MT2 + 3*MV2 - 3*shat)*smT2**3 + &
!                          4*smT2**4)*SP(4)) - &
!                    2*E**(6*y)*smT* &
!                     (2*(-4*MV2**4*shat**4*(-4*MT2 + MV2 + shat) + &
!                          2*MV2**3*(-34*MT2 + MV2 - 15*shat)*shat**3* &
!                           smT2 + &
!                          4*MV2**2*shat**2*(22*MT2 + 7*MV2 + 11*shat)* &
!                           smT2**2 - &
!                          MV2*shat*(46*MT2 + 21*MV2 + 27*shat)*smT2**3 + &
!                          (10*MT2 + 5*MV2 + 6*shat)*smT2**4)*SP(3) + &
!                       shat*(4*MV2**4*shat**3* &
!                           (12*MT2 - 6*MV2 + 17*shat) - &
!                          MV2**3*shat**2*(84*MT2 + 3*MV2 + 130*shat)* &
!                           smT2 + &
!                          MV2**2*shat*(44*MT2 + 13*MV2 + 127*shat)* &
!                           smT2**2 - &
!                          MV2*(8*MT2 + 6*MV2 + 61*shat)*smT2**3 + &
!                          11*smT2**4)*SP(4)) + &
!                    E**(9*y)*smT2* &
!                     (-8*(MV2*shat - smT2)**3*(12*MV2*shat + 5*smT2)* &
!                        SP(3) + &
!                       (148*MV2**4*shat**4 - 294*MV2**3*shat**3*smT2 + &
!                          254*MV2**2*shat**2*smT2**2 - &
!                          113*MV2*shat*smT2**3 + 25*smT2**4)*SP(4)) + &
!                    E**(3*y)*smT2**2* &
!                     (-4*MV2**3*shat**3*(SP(3) - 4*SP(4)) + &
!                       smT2**3*(2*SP(3) + 7*SP(4)) + &
!                       4*MV2**2*shat**2*smT2*(3*SP(3) + 7*SP(4)) - &
!                       MV2*shat*smT2**2*(18*SP(3) + 35*SP(4))) + &
!                    2*E**(5*y)*smT2* &
!                     (MV2**3*shat**3*smT2*(36*SP(3) - 41*SP(4)) + &
!                       2*smT2**4*(3*SP(3) + 5*SP(4)) - &
!                       2*MV2**4*shat**4*(8*SP(3) + 9*SP(4)) + &
!                       6*MV2**2*shat**2*smT2**2*(3*SP(3) + 16*SP(4)) - &
!                       MV2*shat*smT2**3*(28*SP(3) + 57*SP(4))) - &
!                    2*E**(7*y)* &
!                     (6*MV2**2*shat**2*smT2**3*(2*SP(3) - 31*SP(4)) + &
!                       8*MV2**4*shat**4*smT2*(8*SP(3) - 19*SP(4)) + &
!                       40*MV2**5*shat**5*SP(4) - &
!                       15*smT2**5*(SP(3) + SP(4)) + &
!                       3*MV2**3*shat**3*smT2**2*(-24*SP(3) + 77*SP(4)) + &
!                       MV2*shat*smT2**4*(31*SP(3) + 82*SP(4)))))/ &
!                (2q0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               (LC(3)*(-2*E**(15*y)*MV2*shat*smT**7*(3*MV2*shat - smT2)* &
!                     SP(3) + E**(16*y)*shat*(3*MV2*shat - smT2)*smT2**4* &
!                     SP(3) - 4*E**y*MV2**2*shat**2*smT**7*SP(4) + &
!                    2*MV2*shat**2*smT2**4*SP(4) - &
!                    E**(14*y)*smT2**3* &
!                     ((MV2**2*shat**2*(8*MT2 - 4*MV2 + 31*shat) - &
!                          2*MV2*shat*(6*MT2 + 13*shat)*smT2 + &
!                          (4*MT2 + 5*shat)*smT2**2)*SP(3) + &
!                       MV2*shat**2*(3*MV2*shat - smT2)*SP(4)) + &
!                    2*E**(3*y)*MV2**2*shat**2*smT**5* &
!                     (2*smT2*SP(3) + 19*MV2*shat*SP(4) - 11*smT2*SP(4)) &
!                     + 2*E**(7*y)*MV2*shat*smT* &
!                     ((32*MV2**4*shat**4 - 46*MV2**3*shat**3*smT2 + &
!                          36*MV2**2*shat**2*smT2**2 - &
!                          17*MV2*shat*smT2**3 + 5*smT2**4)*SP(3) - &
!                       10*MV2*shat*(3*MV2*shat - 2*smT2)* &
!                        (MV2*shat - smT2)*smT2*SP(4)) + &
!                    E**(2*y)*smT2**3* &
!                     (-2*MV2**2*shat**3*SP(3) + &
!                       (4*MV2**2*(-2*MT2 + MV2 - 5*shat)*shat**2 + &
!                          3*MV2*shat*(4*MT2 + 3*shat)*smT2 + &
!                          (-4*MT2 + shat)*smT2**2)*SP(4)) - &
!                    E**(10*y)*smT2* &
!                     (2*(4*MV2**4*shat**4*(30*MT2 - 15*MV2 + 11*shat) + &
!                          12*MV2**3*(-25*MT2 + 5*MV2 - 6*shat)*shat**3* &
!                           smT2 + &
!                          4*MV2**2*shat**2*(70*MT2 - 5*MV2 + 18*shat)* &
!                           smT2**2 - &
!                          2*MV2*shat*(60*MT2 + 17*shat)*smT2**3 + &
!                          5*(4*MT2 + shat)*smT2**4)*SP(3) + &
!                       (80*MV2**4*shat**4*(MT2 + shat) - &
!                          40*MV2**3*shat**3*(5*MT2 + MV2 + 2*shat)* &
!                           smT2 + &
!                          10*MV2**2*shat**2*(20*MT2 + 2*MV2 + shat)* &
!                           smT2**2 + &
!                          MV2*shat*(-100*MT2 + 17*shat)*smT2**3 - &
!                          5*(-4*MT2 + shat)*smT2**4)*SP(4)) - &
!                    E**(6*y)*smT2* &
!                     ((16*MV2**4*shat**4*(5*MT2 + 4*shat) - &
!                          8*MV2**3*shat**3*(25*MT2 + 5*MV2 + 8*shat)* &
!                           smT2 + &
!                          MV2**2*shat**2*(200*MT2 + 20*MV2 + 23*shat)* &
!                           smT2**2 - &
!                          2*MV2*shat*(50*MT2 + shat)*smT2**3 + &
!                          (20*MT2 + shat)*smT2**4)*SP(3) + &
!                       5*(8*MV2**4*shat**4*(6*MT2 - 3*MV2 + shat) + &
!                          24*MV2**3*(-5*MT2 + MV2)*shat**3*smT2 + &
!                          MV2**2*(112*MT2 - 8*MV2 - 5*shat)*shat**2* &
!                           smT2**2 + &
!                          MV2*shat*(-48*MT2 + 5*shat)*smT2**3 - &
!                          2*(-4*MT2 + shat)*smT2**4)*SP(4)) - &
!                    2*E**(11*y)*MV2*shat*smT**3* &
!                     (-10*smT2**3*SP(3) + &
!                       7*MV2*shat*smT2**2*(10*SP(3) - SP(4)) + &
!                       5*MV2**2*shat**2*smT2*(-28*SP(3) + SP(4)) + &
!                       10*MV2**3*shat**3*(9*SP(3) + SP(4))) + &
!                    2*E**(13*y)*MV2*shat*smT**5* &
!                     (5*smT2**2*SP(3) + &
!                       MV2*shat*smT2*(-27*SP(3) + SP(4)) + &
!                       MV2**2*shat**2*(30*SP(3) + SP(4))) - &
!                    2*E**(5*y)*MV2*shat*smT**3* &
!                     (-(smT2**3*SP(3)) + &
!                       MV2**2*shat**2*smT2*(14*SP(3) - 65*SP(4)) + &
!                       MV2*shat*smT2**2*(-7*SP(3) + 25*SP(4)) + &
!                       2*MV2**3*shat**3*(SP(3) + 25*SP(4))) + &
!                    4*E**(9*y)*MV2*shat*smT* &
!                     (5*smT2**4*SP(3) + &
!                       MV2*shat*smT2**3*(-34*SP(3) + SP(4)) + &
!                       24*MV2**4*shat**4*(SP(3) + SP(4)) + &
!                       MV2**2*shat**2*smT2**2*(72*SP(3) + 7*SP(4)) - &
!                       MV2**3*shat**3*smT2*(67*SP(3) + 27*SP(4))) + &
!                    E**(12*y)*smT2**2* &
!                     (shat*(20*MV2**3*shat**2*(-2*MV2 + 5*shat) + &
!                          MV2**2*(20*MV2 - 133*shat)*shat*smT2 + &
!                          65*MV2*shat*smT2**2 - 10*smT2**3)*SP(3) + &
!                       shat*(28*MV2**3*shat**3 - &
!                          MV2**2*shat*(4*MV2 + 19*shat)*smT2 + smT2**3)* &
!                        SP(4) + &
!                       4*MT2*(MV2*shat - smT2)* &
!                        (5*(-2*MV2*shat + smT2)**2*SP(3) + &
!                          (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
!                           SP(4))) + &
!                    E**(4*y)*smT2**2* &
!                     (4*MT2*(MV2*shat - smT2)* &
!                        ((2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
!                           SP(3) + 5*(-2*MV2*shat + smT2)**2*SP(4)) + &
!                       shat*(MV2*shat* &
!                           (4*MV2**2*(5*shat**2 - smT2) - &
!                             9*MV2*shat*smT2 - smT2**2)*SP(3) + &
!                          (20*MV2**3*shat**2*(-2*MV2 + 3*shat) + &
!                             MV2**2*(20*MV2 - 47*shat)*shat*smT2 + &
!                             4*MV2*shat*smT2**2 + 5*smT2**3)*SP(4))) + &
!                    E**(8*y)*(8*MT2*(MV2*shat - smT2)* &
!                        (24*MV2**4*shat**4 - 48*MV2**3*shat**3*smT2 + &
!                          49*MV2**2*shat**2*smT2**2 - &
!                          25*MV2*shat*smT2**3 + 5*smT2**4)* &
!                        (SP(3) + SP(4)) + &
!                       shat*(MV2*shat*smT2**4*(29*SP(3) - 38*SP(4)) - &
!                          5*smT2**5*(SP(3) - 2*SP(4)) + &
!                          58*MV2**2*shat**2*smT2**3*(-SP(3) + SP(4)) - &
!                          48*MV2**6*shat**4*(SP(3) + SP(4)) + &
!                          8*MV2**5*shat**3*(2*shat**2 + 3*smT2)* &
!                           (SP(3) + SP(4)) + &
!                          8*MV2**3*shat*smT2**2* &
!                           (shat**2*(5*SP(3) - 7*SP(4)) + &
!                             smT2*(SP(3) + SP(4))) + &
!                          8*MV2**4*shat**2*smT2* &
!                           (-3*smT2*(SP(3) + SP(4)) + &
!                             shat**2*(SP(3) + 5*SP(4)))))))/ &
!                (E**(8*y)*shat**2*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               (2*(E**(5*y)*shat*smT**3 + &
!                    E**y*shat*smT*(2*MV2*shat - smT2) - &
!                    2*E**(3*y)*shat*smT*(2*MV2*shat - smT2) - &
!                    E**(4*y)*MV2*shat*smT2 + smT2*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*(4*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
!                       smT2**2))*LC(5)*SP(5))/ &
!                (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
!               ((-2*E**(14*y)*shat*smT**7*(3*MV2*shat - 2*smT2) + &
!                    2*shat*smT**7*(2*MV2*shat - smT2) - &
!                    4*E**y*MV2**2*shat**2*smT2**3 + &
!                    E**(13*y)*smT2**3* &
!                     (2*MV2**2*shat**2 + 3*MV2*shat*smT2 - smT2**2) - &
!                    E**(2*y)*shat*smT**5* &
!                     (-2*MV2**2*(-2*MT2 + MV2 - 22*shat)*shat - &
!                       MV2*(4*MT2 + 2*MV2 + 43*shat)*smT2 + 9*smT2**2) + &
!                    E**(12*y)*shat*smT**5* &
!                     (-2*MV2**2*(-2*MT2 + MV2 - 32*shat)*shat - &
!                       MV2*(4*MT2 + 2*MV2 + 73*shat)*smT2 + 19*smT2**2) &
!                     - E**(10*y)*shat*smT**3* &
!                     (2*MV2**3*shat**2*(20*MT2 - 10*MV2 + 109*shat) - &
!                       2*MV2**2*shat*(30*MT2 + 5*MV2 + 172*shat)*smT2 + &
!                       MV2*(20*MT2 + 10*MV2 + 183*shat)*smT2**2 - &
!                       35*smT2**3) + &
!                    E**(4*y)*shat*smT**3* &
!                     (2*MV2**3*shat**2*(20*MT2 - 10*MV2 + 77*shat) - &
!                       2*MV2**2*shat*(30*MT2 + 5*MV2 + 106*shat)*smT2 + &
!                       MV2*(20*MT2 + 10*MV2 + 93*shat)*smT2**2 - &
!                       13*smT2**3) - &
!                    E**(3*y)*smT2**2* &
!                     (-42*MV2**3*shat**3 + 28*MV2**2*shat**2*smT2 - &
!                       3*MV2*shat*smT2**2 + smT2**3) - &
!                    E**(11*y)*smT2**2* &
!                     (20*MV2**3*shat**3 + 16*MV2**2*shat**2*smT2 - &
!                       25*MV2*shat*smT2**2 + 5*smT2**3) - &
!                    2*E**(6*y)*shat*smT* &
!                     (12*MV2**4*shat**3*(4*MT2 - 3*MV2 + 6*shat) + &
!                       MV2**3*(-84*MT2 + 18*MV2 - 85*shat)*shat**2* &
!                        smT2 + &
!                       2*MV2**2*shat*(22*MT2 + MV2 + 19*shat)*smT2**2 - &
!                       MV2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3 - smT2**4) - &
!                    10*E**(9*y)*smT2* &
!                     (-6*MV2**4*shat**4 + MV2**3*shat**3*smT2 + &
!                       8*MV2**2*shat**2*smT2**2 - 6*MV2*shat*smT2**3 + &
!                       smT2**4) - &
!                    5*E**(5*y)*smT2* &
!                     (28*MV2**4*shat**4 - 42*MV2**3*shat**3*smT2 + &
!                       22*MV2**2*shat**2*smT2**2 - 5*MV2*shat*smT2**3 + &
!                       smT2**4) + &
!                    2*E**(7*y)*(MV2*shat - smT2)* &
!                     (24*MV2**4*shat**4 - 48*MV2**3*shat**3*smT2 + &
!                       49*MV2**2*shat**2*smT2**2 - 25*MV2*shat*smT2**3 + &
!                       5*smT2**4) + &
!                    2*E**(8*y)*shat*smT* &
!                     (4*MV2**4*shat**3*(12*MT2 - 9*MV2 + 26*shat) + &
!                       MV2**3*(-84*MT2 + 18*MV2 - 181*shat)*shat**2* &
!                        smT2 + &
!                       2*MV2**2*shat*(22*MT2 + MV2 + 77*shat)*smT2**2 - &
!                       4*MV2*(2*MT2 + MV2 + 19*shat)*smT2**3 + 14*smT2**4 &
!                       ))*LC(4)*SP(6))/ &
!                (2q0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((E**(5*y)*(-32*MV2**4*shat**4*smT - &
!                       12*MV2**3*shat**3*smT**3 + &
!                       48*MV2**2*shat**2*smT**5 - 27*MV2*shat*smT**7 + &
!                       3*smT**9) + 2*shat*smT2**4 - &
!                    2*E**(14*y)*shat*smT2**4 - &
!                    E**y*smT**7*(3*MV2*shat + smT2) + &
!                    E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
!                    2*E**(11*y)*smT**5*(3*MV2*shat + smT2)* &
!                     (-5*MV2*shat + 3*smT2) + &
!                    4*E**(7*y)*smT*(-(MV2*shat) + smT2)**3* &
!                     (8*MV2*shat + 3*smT2) + &
!                    2*E**(2*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
!                       (2*MT2 + MV2 + 5*shat)*smT2) - &
!                    2*E**(12*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
!                       (2*MT2 + MV2 + 5*shat)*smT2) - &
!                    2*E**(10*y)*smT2**2* &
!                     (-10*MV2**2*(-2*MT2 + MV2 - 3*shat)*shat**2 - &
!                       5*MV2*shat*(6*MT2 + MV2 + 6*shat)*smT2 + &
!                       (10*MT2 + 5*MV2 + 11*shat)*smT2**2) - &
!                    2*E**(3*y)*smT**3* &
!                     (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
!                       10*MV2*shat*smT2**2 + smT2**3) + &
!                    E**(9*y)*smT**3* &
!                     (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
!                       11*MV2*shat*smT2**2 + 13*smT2**3) + &
!                    2*E**(4*y)*smT2* &
!                     (4*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
!                       2*MV2**2*(2*MT2 + 11*MV2 - 9*shat)*shat**2*smT2 + &
!                       MV2*(-6*MT2 + 7*MV2 - 18*shat)*shat*smT2**2 + &
!                       (2*MT2 + MV2 + 7*shat)*smT2**3) + &
!                    2*E**(8*y)*smT2* &
!                     (8*MV2**3*shat**3*(7*MT2 - 4*MV2 + 2*shat) + &
!                       6*MV2**2*(-18*MT2 + MV2 - 3*shat)*shat**2*smT2 + &
!                       2*MV2*shat*(34*MT2 + 7*(MV2 + 2*shat))*smT2**2 - &
!                       (16*MT2 + 8*MV2 + 11*shat)*smT2**3) - &
!                    2*E**(6*y)* &
!                     (-8*MV2**4*shat**4*(-4*MT2 + MV2 + shat) - &
!                       4*MV2**3*(16*MT2 + 5*MV2 - 7*shat)*shat**3*smT2 + &
!                       6*MV2**2*(10*MT2 + 3*MV2 - shat)*shat**2* &
!                        smT2**2 - 2*MV2*(18*MT2 + 7*MV2)*shat*smT2**3 + &
!                       (8*MT2 + 4*MV2 + shat)*smT2**4))*LC(10)*SP(6))/ &
!                (E**(6*y)*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((-2*E**(14*y)*MV2*shat**2*smT2**3 + 2*shat*smT2**4 - &
!                    E**y*smT**7*(3*MV2*shat + smT2) + &
!                    E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
!                    2*E**(3*y)*smT**5*(5*MV2*shat - 3*smT2)* &
!                     (3*MV2*shat + smT2) + &
!                    4*E**(7*y)*smT*(MV2*shat - smT2)**3* &
!                     (8*MV2*shat + 3*smT2) + &
!                    2*E**(2*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 11*shat)*shat + &
!                       (2*MT2 + MV2 + 6*shat)*smT2) - &
!                    2*E**(12*y)*smT2**2* &
!                     (-10*MV2**2*shat**3 + &
!                       MV2*shat*(-2*MT2 + MV2 + 4*shat)*smT2 + &
!                       (2*MT2 + MV2 + shat)*smT2**2) + &
!                    2*E**(4*y)*smT2**2* &
!                     (-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
!                       MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
!                       5*(2*MT2 + MV2 + 3*shat)*smT2**2) + &
!                    2*E**(11*y)*smT**3* &
!                     (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
!                       10*MV2*shat*smT2**2 + smT2**3) - &
!                    E**(5*y)*smT**3* &
!                     (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
!                       11*MV2*shat*smT2**2 + 13*smT2**3) - &
!                    2*E**(10*y)*smT2* &
!                     (4*MV2**3*shat**3*(2*MT2 + MV2 + 9*shat) - &
!                       2*MV2**2*shat**2*(2*MT2 + 11*MV2 + 19*shat)* &
!                        smT2 + &
!                       MV2*shat*(-6*MT2 + 7*MV2 + 13*shat)*smT2**2 + &
!                       (2*MT2 + MV2)*smT2**3) + &
!                    2*E**(6*y)*smT2* &
!                     (8*MV2**3*(-7*MT2 + 4*MV2 - 6*shat)*shat**3 + &
!                       2*MV2**2*shat**2*(54*MT2 - 3*MV2 + 37*shat)* &
!                        smT2 - &
!                       MV2*shat*(68*MT2 + 14*MV2 + 59*shat)*smT2**2 + &
!                       2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3) + &
!                    E**(9*y)*smT* &
!                     (32*MV2**4*shat**4 + 12*MV2**3*shat**3*smT2 - &
!                       48*MV2**2*shat**2*smT2**2 + 27*MV2*shat*smT2**3 - &
!                       3*smT2**4) + &
!                    2*E**(8*y)* &
!                     (-8*MV2**4*(-4*MT2 + MV2 - 3*shat)*shat**4 - &
!                       4*MV2**3*shat**3*(16*MT2 + 5*MV2 + 9*shat)*smT2 + &
!                       6*MV2**2*shat**2*(10*MT2 + 3*MV2 + 9*shat)* &
!                        smT2**2 - &
!                       2*MV2*shat*(18*MT2 + 7*MV2 + 18*shat)*smT2**3 + &
!                       (8*MT2 + 4*MV2 + 9*shat)*smT2**4))*LC(14)*SP(6))/ &
!                (E**(8*y)*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((-2*shat*smT**9 + &
!                    E**(14*y)*shat*smT**7*(MV2*shat + smT2) + &
!                    E**(13*y)*smT2**3* &
!                     (-(MV2**2*shat**2) - 4*MV2*shat*smT2 + smT2**2) + &
!                    E**y*smT2**3* &
!                     (6*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2) - &
!                    2*E**(2*y)*shat*smT**5* &
!                     (MV2**2*shat*(-2*MT2 + MV2 + shat) + &
!                       MV2*(2*MT2 + MV2 - 10*shat)*smT2 + 4*smT2**2) + &
!                    E**(4*y)*shat*smT**3* &
!                     (2*MV2**3*shat**2*(-20*MT2 + 10*MV2 + 11*shat) + &
!                       2*MV2**2*(30*MT2 + 5*MV2 - 38*shat)*shat*smT2 - &
!                       5*MV2*(4*MT2 + 2*MV2 - 9*shat)*smT2**2 - &
!                       13*smT2**3) + &
!                    2*E**(10*y)*shat*smT**3* &
!                     (5*MV2**3*shat**2*(4*MT2 - 2*MV2 + shat) - &
!                       5*MV2**2*(6*MT2 + MV2 - shat)*shat*smT2 + &
!                       5*MV2*(2*MT2 + MV2)*smT2**2 + smT2**3) + &
!                    E**(11*y)*smT2**2* &
!                     (8*MV2**3*shat**3 + 35*MV2**2*shat**2*smT2 - &
!                       33*MV2*shat*smT2**2 + 6*smT2**3) + &
!                    E**(3*y)*smT2**2* &
!                     (-62*MV2**3*shat**3 + 67*MV2**2*shat**2*smT2 - &
!                       27*MV2*shat*smT2**2 + 6*smT2**3) + &
!                    E**(6*y)*shat*smT* &
!                     (-8*MV2**4*shat**3*(-12*MT2 + 9*MV2 + 10*shat) + &
!                       2*MV2**3*shat**2*(-84*MT2 + 18*MV2 + 83*shat)* &
!                        smT2 + &
!                       4*MV2**2*(22*MT2 + MV2 - 22*shat)*shat*smT2**2 + &
!                       MV2*(-8*(2*MT2 + MV2) + 43*shat)*smT2**3 - &
!                       11*smT2**4) + &
!                    2*E**(8*y)*shat*smT* &
!                     (12*MV2**4*shat**3*(-4*MT2 + 3*MV2 + 2*shat) + &
!                       MV2**3*(84*MT2 - 18*MV2 - 35*shat)*shat**2*smT2 - &
!                       2*MV2**2*shat*(22*MT2 + MV2 + 7*shat)*smT2**2 + &
!                       4*MV2*(2*MT2 + MV2 + 3*shat)*smT2**3 - 2*smT2**4) &
!                     + 2*E**(7*y)*(-(MV2*shat) + smT2)* &
!                     (56*MV2**4*shat**4 - 104*MV2**3*shat**3*smT2 + &
!                       105*MV2**2*shat**2*smT2**2 - &
!                       52*MV2*shat*smT2**3 + 10*smT2**4) + &
!                    E**(9*y)*smT2* &
!                     (-12*MV2**4*shat**4 - 94*MV2**3*shat**3*smT2 + &
!                       164*MV2**2*shat**2*smT2**2 - &
!                       93*MV2*shat*smT2**3 + 15*smT2**4) + &
!                    E**(5*y)*smT2* &
!                     (204*MV2**4*shat**4 - 362*MV2**3*shat**3*smT2 + &
!                       247*MV2**2*shat**2*smT2**2 - &
!                       84*MV2*shat*smT2**3 + 15*smT2**4) - &
!                    E**(12*y)*shat*smT**5* &
!                     (4*MT2*MV2*(MV2*shat - smT2) - &
!                       (MV2*shat + smT2)*(2*MV2*(MV2 - 4*shat) + 3*smT2)) &
!                    )*LC(2)*SP(7))/ &
!                (2q0*E**(7*y)*shat*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((-2*MV2*shat**2*smT2**3 + 2*E**(14*y)*shat*smT2**4 + &
!                    E**y*smT**7*(3*MV2*shat + smT2) - &
!                    E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
!                    2*E**(11*y)*smT**5*(5*MV2*shat - 3*smT2)* &
!                     (3*MV2*shat + smT2) + &
!                    4*E**(7*y)*smT*(MV2*shat - smT2)**3* &
!                     (8*MV2*shat + 3*smT2) + &
!                    2*E**(12*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 11*shat)*shat + &
!                       (2*MT2 + MV2 + 6*shat)*smT2) - &
!                    2*E**(2*y)*smT2**2* &
!                     (-10*MV2**2*shat**3 + &
!                       MV2*shat*(-2*MT2 + MV2 + 4*shat)*smT2 + &
!                       (2*MT2 + MV2 + shat)*smT2**2) + &
!                    2*E**(10*y)*smT2**2* &
!                     (-10*MV2**2*(-2*MT2 + MV2 - 4*shat)*shat**2 - &
!                       MV2*shat*(30*MT2 + 5*MV2 + 44*shat)*smT2 + &
!                       5*(2*MT2 + MV2 + 3*shat)*smT2**2) + &
!                    2*E**(3*y)*smT**3* &
!                     (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
!                       10*MV2*shat*smT2**2 + smT2**3) - &
!                    E**(9*y)*smT**3* &
!                     (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
!                       11*MV2*shat*smT2**2 + 13*smT2**3) - &
!                    2*E**(4*y)*smT2* &
!                     (4*MV2**3*shat**3*(2*MT2 + MV2 + 9*shat) - &
!                       2*MV2**2*shat**2*(2*MT2 + 11*MV2 + 19*shat)* &
!                        smT2 + &
!                       MV2*shat*(-6*MT2 + 7*MV2 + 13*shat)*smT2**2 + &
!                       (2*MT2 + MV2)*smT2**3) + &
!                    2*E**(8*y)*smT2* &
!                     (8*MV2**3*(-7*MT2 + 4*MV2 - 6*shat)*shat**3 + &
!                       2*MV2**2*shat**2*(54*MT2 - 3*MV2 + 37*shat)* &
!                        smT2 - &
!                       MV2*shat*(68*MT2 + 14*MV2 + 59*shat)*smT2**2 + &
!                       2*(8*MT2 + 4*MV2 + 9*shat)*smT2**3) + &
!                    E**(5*y)*smT* &
!                     (32*MV2**4*shat**4 + 12*MV2**3*shat**3*smT2 - &
!                       48*MV2**2*shat**2*smT2**2 + 27*MV2*shat*smT2**3 - &
!                       3*smT2**4) + &
!                    2*E**(6*y)* &
!                     (-8*MV2**4*(-4*MT2 + MV2 - 3*shat)*shat**4 - &
!                       4*MV2**3*shat**3*(16*MT2 + 5*MV2 + 9*shat)*smT2 + &
!                       6*MV2**2*shat**2*(10*MT2 + 3*MV2 + 9*shat)* &
!                        smT2**2 - &
!                       2*MV2*shat*(18*MT2 + 7*MV2 + 18*shat)*smT2**3 + &
!                       (8*MT2 + 4*MV2 + 9*shat)*smT2**4))*LC(8)*SP(7))/ &
!                (E**(6*y)*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2) + &
!               ((E**(9*y)*(-32*MV2**4*shat**4*smT - &
!                       12*MV2**3*shat**3*smT**3 + &
!                       48*MV2**2*shat**2*smT**5 - 27*MV2*shat*smT**7 + &
!                       3*smT**9) - 2*shat*smT2**4 + &
!                    2*E**(14*y)*shat*smT2**4 + &
!                    E**y*smT**7*(3*MV2*shat + smT2) - &
!                    E**(13*y)*smT**7*(3*MV2*shat + smT2) + &
!                    2*E**(3*y)*smT**5*(3*MV2*shat + smT2)* &
!                     (-5*MV2*shat + 3*smT2) + &
!                    4*E**(7*y)*smT*(-(MV2*shat) + smT2)**3* &
!                     (8*MV2*shat + 3*smT2) - &
!                    2*E**(2*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
!                       (2*MT2 + MV2 + 5*shat)*smT2) + &
!                    2*E**(12*y)*smT2**3* &
!                     (MV2*(-2*MT2 + MV2 - 10*shat)*shat + &
!                       (2*MT2 + MV2 + 5*shat)*smT2) - &
!                    2*E**(4*y)*smT2**2* &
!                     (-10*MV2**2*(-2*MT2 + MV2 - 3*shat)*shat**2 - &
!                       5*MV2*shat*(6*MT2 + MV2 + 6*shat)*smT2 + &
!                       (10*MT2 + 5*MV2 + 11*shat)*smT2**2) - &
!                    2*E**(11*y)*smT**3* &
!                     (2*MV2**3*shat**3 - 21*MV2**2*shat**2*smT2 + &
!                       10*MV2*shat*smT2**2 + smT2**3) + &
!                    E**(5*y)*smT**3* &
!                     (92*MV2**3*shat**3 - 96*MV2**2*shat**2*smT2 + &
!                       11*MV2*shat*smT2**2 + 13*smT2**3) + &
!                    2*E**(10*y)*smT2* &
!                     (4*MV2**3*shat**3*(2*MT2 + MV2 + shat) - &
!                       2*MV2**2*(2*MT2 + 11*MV2 - 9*shat)*shat**2*smT2 + &
!                       MV2*(-6*MT2 + 7*MV2 - 18*shat)*shat*smT2**2 + &
!                       (2*MT2 + MV2 + 7*shat)*smT2**3) + &
!                    2*E**(6*y)*smT2* &
!                     (8*MV2**3*shat**3*(7*MT2 - 4*MV2 + 2*shat) + &
!                       6*MV2**2*(-18*MT2 + MV2 - 3*shat)*shat**2*smT2 + &
!                       2*MV2*shat*(34*MT2 + 7*(MV2 + 2*shat))*smT2**2 - &
!                       (16*MT2 + 8*MV2 + 11*shat)*smT2**3) - &
!                    2*E**(8*y)* &
!                     (-8*MV2**4*shat**4*(-4*MT2 + MV2 + shat) - &
!                       4*MV2**3*(16*MT2 + 5*MV2 - 7*shat)*shat**3*smT2 + &
!                       6*MV2**2*(10*MT2 + 3*MV2 - shat)*shat**2* &
!                        smT2**2 - 2*MV2*(18*MT2 + 7*MV2)*shat*smT2**3 + &
!                       (8*MT2 + 4*MV2 + shat)*smT2**4))*LC(12)*SP(7))/ &
!                (E**(8*y)*(MV2*shat - smT2)**3* &
!                  (-2*MV2*shat + smT2 + smT2*Cosh(2*y))**2)) + &
!            (8*kapPr*(MV2*shat - smT2*Cosh(2*y))*(SP(3) + SP(4))* &
!               (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)))/ &
!             (-(MV2*shat) + smT2)**2) + &
!         SI(16)*(kap*(((-2*E**(2*y)*shat*smT2* &
!                     (2*MT2*shat - MV2*shat + smT2) + &
!                    2*E**y*smT*(-(MV2*shat) + smT2)* &
!                     (shat*(12*MT2 - 3*MV2 + shat) + 4*smT2) + &
!                    4*shat*(MT2*shat*(4*MV2*shat - 3*smT2) - &
!                       (-(MV2*shat) + smT2)**2))*LC(1))/ &
!                (shat**2*(MV2*shat - smT2)) + &
!               ((shat*smT2*(-((2*MT2 + MV2)*shat) + smT2) - &
!                    2*E**y*smT* &
!                     (MV2*(-4*MT2 + MV2)*shat**2 + &
!                       (12*MT2 - 5*MV2)*shat*smT2 + 4*smT2**2) - &
!                    2*E**(2*y)*shat* &
!                     (MT2*shat*(4*MV2*shat - 5*smT2) - &
!                       (-(MV2*shat) + smT2)**2))*LC(11))/ &
!                (E**(2*y)*shat**2*(MV2*shat - smT2)) + &
!               smT*LC(15)*(E**(-y) + &
!                  (8*MT2*shat*smT - 2*MV2*shat*smT + 2*smT**3 + &
!                     4*MT2*shat**2*Sinh(y))/(MV2*shat**2 - shat*smT2)) + &
!               (4*LC(7)*(2*MT2*shat*(-2*MV2*shat + smT2) + &
!                    (-(MV2*shat) + smT2)**2 + &
!                    shat*smT*(2*MT2*shat*Cosh(y) + &
!                       (2*MT2*shat - MV2*shat + smT2)*Sinh(y)))*SP(1))/ &
!                (shat**2*(MV2*shat - smT2)) - &
!               (4*LC(6)*(2*MT2*shat*(-2*MV2*shat + smT2) + &
!                    (-(MV2*shat) + smT2)**2 + &
!                    shat*smT*(2*MT2*shat*Cosh(y) + &
!                       (2*MT2*shat - MV2*shat + smT2)*Sinh(y)))*SP(2))/ &
!                (shat**2*(MV2*shat - smT2)) + &
!               (4*smT2*(4*MT2*shat - MV2*shat + smT2)*LC(13)*SP(4))/ &
!                (E**(2*y)*shat**2*(MV2*shat - smT2)) - &
!               (2*smT*LC(9)*(2*E**(3*y)*smT* &
!                     (4*MT2*shat - MV2*shat + smT2)*SP(3) - &
!                    2*E**(2*y)*MT2*shat**2*SP(4) + &
!                    shat*(2*MT2*shat - MV2*shat + smT2)*SP(4)))/ &
!                (E**y*shat**2*(MV2*shat - smT2)) + &
!               (4*LC(3)*(E**(3*y)*shat*smT2* &
!                     (2*MT2*shat - MV2*shat + smT2)*SP(3) - &
!                    2*E**(2*y)*smT* &
!                     ((-(MV2*shat) + smT2)**2 + &
!                       MT2*shat*(-4*MV2*shat + 3*smT2))*SP(3) + &
!                    smT*(-(MV2*shat**2*(-8*MT2 + 2*MV2 + shat)) + &
!                       shat*(-6*MT2 + 4*MV2 + shat)*smT2 - 2*smT2**2)* &
!                     SP(4) + E**y*shat* &
!                     (2*MT2*shat*(2*MV2*shat - smT2) - &
!                       (-(MV2*shat) + smT2)**2)*SP(4)))/ &
!                (E**y*shat**3*(MV2*shat - smT2)) + &
!               (2*smT*(-2*MT2*shat + &
!                    E**(2*y)*(10*MT2*shat - MV2*shat + smT2))*LC(4)*SP(6) &
!                  )/(E**y*shat*(MV2*shat - smT2)) + &
!               (4*(E**y*shat + smT)*(4*MT2*shat - MV2*shat + smT2)* &
!                  LC(14)*SP(6))/(E**y*shat*(MV2*shat - smT2)) + &
!               (2*smT*(-6*MT2*shat - &
!                    E**(2*y)*(2*MT2*shat - MV2*shat + smT2))*LC(2)*SP(7)) &
!                 /(E**y*shat*(MV2*shat - smT2)) + &
!               (4*(4*MT2*shat - MV2*shat + smT2)*LC(8)*SP(7))/ &
!                (MV2*shat - smT2) - &
!               (4*smT*(4*MT2*shat - MV2*shat + smT2)*LC(12)*SP(7))/ &
!                (E**y*shat*(MV2*shat - smT2))) + &
!            kapPr*(-2*E**y*smT*(SP(3)*SP(5) + SP(2)*SP(6)) - &
!               (2*smT*(SP(4)*SP(5) + SP(1)*SP(7)))/E**y + &
!               (4*(4*MT2*shat - MV2*shat + smT2)*(SP(3) + SP(4))* &
!                  (MV2*shat*SP(5) - smT2*SP(5) - shat*SP(6)*SP(7)))/ &
!                (shat*(MV2*shat - smT2)))) + &
!         SI(12)*(kap*(((-2*shat*smT*(2*MV2*shat - 3*smT2)* &
!                     (MV2*shat - smT2)*(2*MV2*shat - smT2) - &
!                    E**(4*y)*shat**3*smT**3*(5*MV2*shat + 7*smT2) + &
!                    E**(3*y)*shat**2*smT2* &
!                     (8*MV2**2*shat**2 + 17*MV2*shat*smT2 + 11*smT2**2) &
!                     + E**(2*y)*shat*smT* &
!                     (MV2**2*(8*MT2 - 19*MV2 - 4*shat)*shat**3 + &
!                       MV2*shat**2*(8*MT2 + 35*MV2 + 10*shat)*smT2 - &
!                       2*shat*(8*MT2 + 32*MV2 + 3*shat)*smT2**2 + &
!                       12*smT2**3) + &
!                    E**y*(8*MV2**3*shat**4*(4*MT2 + MV2 + shat) - &
!                       MV2**2*shat**3*(13*(8*MT2 + MV2) + 20*shat)* &
!                        smT2 + &
!                       MV2*shat**2*(88*MT2 + 5*MV2 + 12*shat)*smT2**2 + &
!                       4*(-4*MT2 + 5*MV2)*shat*smT2**3 - 8*smT2**4))* &
!                  LC(1))/(E**y*shat**2*(MV2*shat - smT2)**3) + &
!               ((E**y*shat - smT)* &
!                  (shat*(MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2 + &
!                    E**(3*y)*shat**2*smT*(MV2*shat - 3*smT2)* &
!                     (2*MV2*shat + smT2) + &
!                    E**(2*y)*shat* &
!                     (-4*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 + &
!                       5*MV2*shat*smT2**2 + 2*smT2**3) - &
!                    E**y*smT*(2*(2*MT2 - MV2)*MV2**2*shat**3 + &
!                       MV2*shat**2*(-20*MT2 + 16*MV2 + shat)*smT2 + &
!                       (16*MT2 - 16*MV2 - shat)*shat*smT2**2 + 8*smT2**3) &
!                    )*LC(11))/(E**(3*y)*shat**2*(MV2*shat - smT2)**3) + &
!               ((E**y*shat - smT)* &
!                  (shat*smT**3*(MV2*shat - smT2) + &
!                    2*E**(3*y)*shat**2*smT2*(MV2*shat + 4*smT2) - &
!                    E**(2*y)*shat*smT*(3*MV2*shat + smT2)* &
!                     (MV2*shat + 4*smT2) + &
!                    E**y*(MV2**2*shat**3*(-4*MT2 + MV2 + 2*shat) + &
!                       MV2*(-12*MT2 + 7*MV2 - 5*shat)*shat**2*smT2 + &
!                       shat*(16*MT2 + 3*shat)*smT2**2 + 2*smT2**3))* &
!                  LC(15))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(-3*E**(3*y)*MV2*shat**3*smT + E**(4*y)*shat**3*smT2 + &
!                    shat*(MV2*shat - smT2)*smT2 + &
!                    E**(2*y)*shat* &
!                     (MV2*shat**2*(MV2 + 2*shat) + &
!                       (3*MV2 - 2*shat)*shat*smT2 - smT2**2) + &
!                    E**y*smT*(MV2*(MV2 - 3*shat)*shat**2 + &
!                       shat*(-4*MV2 + 3*shat)*smT2 + 2*smT2**2))*LC(7)* &
!                  SP(1))/(E**(2*y)*shat**2*(-(MV2*shat) + smT2)**2) + &
!               (2*(E**(2*y)*shat**2*(-2*MV2 + shat) + &
!                    2*E**(3*y)*shat**2*smT - shat*smT2 + &
!                    2*E**y*smT*(-(MV2*shat) + smT2))*LC(6)*SP(2))/ &
!                (E**(2*y)*shat**2*(MV2*shat - smT2)) + &
!               (LC(13)*(2*E**(4*y)*shat**2*smT*(-(MV2*shat) + smT2)**2* &
!                     SP(3) - 4*smT*smT2* &
!                     (2*MV2*(-2*MT2 + MV2)*shat**2 + &
!                       (4*MT2 - MV2)*shat*smT2 + smT2**2)*SP(4) + &
!                    4*E**y*shat*smT2* &
!                     (MV2*(-4*MT2 + 3*MV2)*shat**2 + &
!                       2*(2*MT2 + MV2)*shat*smT2 + smT2**2)*SP(4) - &
!                    2*E**(3*y)*shat* &
!                     ((2*MV2*shat + shat**2 - smT2)* &
!                        (-(MV2*shat) + smT2)**2*SP(3) - &
!                       4*MV2*shat**3*smT2*SP(4)) + &
!                    2*E**(2*y)*shat**2*smT* &
!                     ((-(MV2*shat) + smT2)**2*SP(3) - &
!                       2*MV2*shat*(MV2*shat + 5*smT2)*SP(4))))/ &
!                (E**(3*y)*shat**2*(MV2*shat - smT2)**3) + &
!               (LC(9)*(8*E**(7*y)*shat**4*smT2**2*SP(3) - &
!                    4*E**(6*y)*shat**3*smT**3*(3*MV2*shat + 5*smT2)* &
!                     SP(3) + 2*shat*smT**5*(MV2*shat - smT2)*SP(4) - &
!                    2*E**y*shat**2*smT2* &
!                     (MV2**2*shat*(4*MT2 - MV2 + shat) + &
!                       MV2*(-4*MT2 - MV2 + shat)*smT2 - 2*smT2**2)*SP(4) &
!                     - 4*E**(5*y)*shat**2*smT2* &
!                     ((2*MV2*shat**2*(2*MT2 - MV2 + shat) - &
!                          shat*(4*MT2 + 5*MV2 + 2*shat)*smT2 - 5*smT2**2) &
!                         *SP(3) - MV2*shat**3*SP(4)) - &
!                    2*E**(4*y)*shat*smT* &
!                     (-2*(MV2**2*shat**3*(-6*MV2 + shat) + &
!                          2*MV2*shat**2*(4*MT2 + 7*MV2 + 2*shat)*smT2 - &
!                          shat*(8*MT2 + 19*MV2 + 5*shat)*smT2**2 + &
!                          3*smT2**3)*SP(3) + &
!                       shat**2* &
!                        (4*MV2**2*shat**2 + 3*MV2*shat*smT2 + smT2**2)* &
!                        SP(4)) + &
!                    2*E**(3*y)* &
!                     (2*(2*MV2**3*shat**4*(4*MT2 + 2*MV2 + shat) - &
!                          2*MV2**2*shat**3*(12*MT2 + 5*MV2 + 4*shat)* &
!                           smT2 + &
!                          4*MV2*shat**2*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
!                          shat*(-4*MT2 + MV2 + 2*shat)*smT2**3 - smT2**4) &
!                         *SP(3) + &
!                       MV2*shat**3* &
!                        (-(MV2*shat**2*(4*MT2 - 3*MV2 + shat)) + &
!                          shat*(4*MT2 + 3*MV2 + shat)*smT2 + 6*smT2**2)* &
!                        SP(4)) - &
!                    2*E**(2*y)*shat*smT* &
!                     (2*(MV2*shat - smT2)* &
!                        (2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
!                        SP(3) + &
!                       (-2*MV2**2*shat**3*(4*MT2 - 2*MV2 + shat) + &
!                          MV2*shat**2*(8*MT2 + shat)*smT2 + &
!                          shat*(5*MV2 + shat)*smT2**2 - smT2**3)*SP(4)))) &
!                 /(E**(2*y)*shat**2*(E**y*shat - smT)* &
!                  (MV2*shat - smT2)**3) + &
!               (LC(3)*(-4*E**(6*y)*shat**4*smT**3*(3*MV2*shat - smT2)* &
!                     SP(3) + 4*E**(5*y)*shat**3*smT2* &
!                     (4*MV2**2*shat**2 + 7*MV2*shat*smT2 - 3*smT2**2)* &
!                     SP(3) - 8*smT**3* &
!                     (-(MV2**2*shat**3*(-4*MT2 + 2*MV2 + shat)) + &
!                       MV2*shat**2*(-6*MT2 + 3*MV2 + 2*shat)*smT2 - &
!                       shat*(-2*MT2 + 3*MV2 + shat)*smT2**2 + smT2**3)* &
!                     SP(4) + 4*E**(4*y)*shat**2*smT* &
!                     ((MV2**2*shat**3*(8*MT2 + 6*MV2 + 3*shat) - &
!                          4*MV2*shat**2*(3*MT2 + 8*MV2 + shat)*smT2 + &
!                          shat*(4*MT2 + 19*MV2 + shat)*smT2**2 - &
!                          5*smT2**3)*SP(3) + &
!                       MV2*shat**3*(3*MV2*shat - smT2)*SP(4)) - &
!                    4*E**(2*y)*smT* &
!                     ((-2*MV2**3*shat**4*(4*MT2 + 3*shat) + &
!                          MV2**2*shat**3*(16*MT2 + 7*shat)*smT2 - &
!                          2*MV2*shat**2*(6*MT2 - 3*MV2 + shat)*smT2**2 + &
!                          shat*(4*MT2 - 6*MV2 + shat)*smT2**3 + 2*smT2**4 &
!                          )*SP(3) + &
!                       shat**3* &
!                        (-2*MV2**2*shat**2*(-4*MT2 + 4*MV2 + shat) + &
!                          MV2*shat*(-12*MT2 - 5*MV2 + 4*shat)*smT2 + &
!                          (4*MT2 + MV2 - 2*shat)*smT2**2)*SP(4)) - &
!                    4*E**y*shat*smT2* &
!                     (2*MV2**2*shat**2*(MV2*shat - smT2)*SP(3) + &
!                       (2*MV2**2*shat**3*(-8*MT2 + 5*MV2 + 2*shat) - &
!                          8*MV2*shat**2*(-3*MT2 + MV2 + shat)*smT2 + &
!                          shat*(-8*MT2 + 9*MV2 + 4*shat)*smT2**2 - &
!                          3*smT2**3)*SP(4)) - &
!                    4*E**(3*y)*shat* &
!                     ((4*MV2**3*shat**4*(4*MT2 + MV2 + shat) - &
!                          2*MV2**2*shat**3*(16*MT2 + 5*MV2 + shat)* &
!                           smT2 - &
!                          2*MV2*shat**2*(-12*MT2 + 5*MV2 + shat)* &
!                           smT2**2 + (-8*MT2 + 13*MV2)*shat*smT2**3 - &
!                          5*smT2**4)*SP(3) + &
!                       shat**2* &
!                        (2*MV2**3*shat**3 + 10*MV2**2*shat**2*smT2 - &
!                          5*MV2*shat*smT2**2 + smT2**3)*SP(4))))/ &
!                (E**(2*y)*shat**3*(E**y*shat - smT)*(MV2*shat - smT2)**3) &
!                 + (2*(-(E**(3*y)*shat**2*smT2) + &
!                    shat*smT*(MV2**2 - MV2*shat + smT2) + &
!                    E**(2*y)*shat*smT*(MV2*shat + 2*smT2) + &
!                    E**y*(MV2*shat**2*(MV2 + shat) - &
!                       shat*(5*MV2 + shat)*smT2 + smT2**2))*LC(5)*SP(5))/ &
!                (E**y*shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(E**y*shat - smT)* &
!                  (-(E**(2*y)*MV2*smT*(7*MV2*shat - 3*smT2)) + &
!                    MV2*smT*(MV2*shat - smT2) + &
!                    2*E**(3*y)*(3*MV2*shat - 2*smT2)*smT2 + &
!                    E**y*(MV2**2*shat*(-4*MT2 + MV2 + shat) + &
!                       MV2*(4*MT2 + MV2 - 3*shat)*smT2 + 2*smT2**2))* &
!                  LC(4)*SP(6))/(E**(2*y)*(MV2*shat - smT2)**3) - &
!               (4*(E**y*shat - smT)*smT* &
!                  ((MV2 - E**y*smT)* &
!                     (MV2*shat + smT*(-2*E**y*shat + smT)) + &
!                    4*MT2*(-(MV2*shat) + smT2))*LC(10)*SP(6))/ &
!                (-(MV2*shat) + smT2)**3 + &
!               (4*(E**y*shat - smT)*smT* &
!                  (2*MV2**2*shat**2 + smT**3*(E**y*shat + smT) + &
!                    MV2*shat*(E**y*shat*(2*E**y*shat - 5*smT) - smT2) + &
!                    4*MT2*shat*(-(MV2*shat) + smT2))*LC(14)*SP(6))/ &
!                (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(3*MV2*shat*smT2*(-(MV2*shat) + smT2) - &
!                    E**(4*y)*shat**2*smT2*(MV2*shat + smT2) + &
!                    E**y*shat*smT* &
!                     (MV2**2*shat*(-4*MT2 + 5*shat) + &
!                       2*MV2*(2*MT2 + MV2 - 2*shat)*smT2 - smT2**2) + &
!                    2*E**(3*y)*shat*smT* &
!                     (MV2**2*shat**2 + MV2*shat*smT2 + smT2**2) + &
!                    E**(2*y)*(-2*MV2**2*shat**3*(-2*MT2 + MV2 + shat) + &
!                       MV2*shat**2*(-4*MT2 + 2*MV2 + shat)*smT2 + &
!                       shat*(-7*MV2 + shat)*smT2**2 + smT2**3))*LC(2)* &
!                  SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) - &
!               (4*(E**y*shat - smT)*smT* &
!                  (2*MV2**2*shat**2 - MV2*shat*smT*(3*E**y*shat + smT) + &
!                    4*MT2*shat*(-(MV2*shat) + smT2) + &
!                    smT2*(E**y*shat*(2*E**y*shat - smT) + smT2))*LC(12)* &
!                  SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) - &
!               (4*(E**y*shat - smT)*LC(8)* &
!                  (2*MV2**2*shat**2 - &
!                    E**(2*y)*smT**3*(smT - 4*shat*Cosh(y)) + &
!                    E**y*smT*(MV2*(-4*MT2 + MV2 - 2*shat)*shat + &
!                       (4*MT2 + MV2)*smT2 - &
!                       MV2*shat*smT*(5*Cosh(y) + Sinh(y))))*SP(7))/ &
!                (E**y*(MV2*shat - smT2)**3)) + &
!            (4*kapPr*(E**y*shat - smT)* &
!               (-((MV2*shat - smT2)* &
!                    (smT2*(SP(3) + SP(4))*SP(5) + shat**2*SP(2)*SP(6))) &
!                  + 2*shat*(MV2*shat*SP(3) + smT2*SP(4))*SP(6)*SP(7) + &
!                 E**y*shat*smT* &
!                  ((MV2*shat - smT2)* &
!                     ((SP(3) + SP(4))*SP(5) + SP(2)*SP(6)) - &
!                    2*shat*(SP(3) + SP(4))*SP(6)*SP(7))))/ &
!             (E**y*shat*(-(MV2*shat) + smT2)**2)) + &
!         SI(9)*(kap*(((12*MV2*shat**3*smT2**2 + &
!                    2*E**(3*y)*shat*smT*(MV2*shat - smT2)* &
!                     (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2) + &
!                    E**y*shat*smT* &
!                     (4*MV2**3*shat**3 - 23*MV2**2*shat**2*smT2 - &
!                       2*MV2*shat*smT2**2 - 3*smT2**3) + &
!                    E**(2*y)*(-4*MV2**3*shat**5 + &
!                       MV2**2*shat**3*(-8*MV2 + 7*shat)*smT2 + &
!                       4*MV2*(9*MV2 - shat)*shat**2*smT2**2 + &
!                       shat*(-24*MV2 + shat)*smT2**3 + 8*smT2**4))*LC(1)) &
!                 /(E**(2*y)*shat**2*(MV2*shat - smT2)**3) + &
!               (smT2*(2*E**(3*y)*MV2*shat**2*smT*(MV2*shat - smT2) + &
!                    6*shat**2*smT2**2 - &
!                    E**y*shat*smT* &
!                     (-2*MV2**2*shat**2 + 13*MV2*shat*smT2 + smT2**2) + &
!                    E**(2*y)*(MV2**2*(MV2 - shat)*shat**3 - &
!                       7*MV2**2*shat**2*smT2 + &
!                       shat*(20*MV2 + shat)*smT2**2 - 8*smT2**3))*LC(11)) &
!                 /(E**(4*y)*shat**2*(MV2*shat - smT2)**3) + &
!               (smT*(-2*shat**2*smT2*(2*MV2*shat + 3*smT2) + &
!                    E**y*shat*smT* &
!                     (3*MV2**2*shat**2 + 15*MV2*shat*smT2 + 2*smT2**2) - &
!                    E**(2*y)*(MV2**2*(MV2 - 2*shat)*shat**3 + &
!                       MV2*shat**2*(3*MV2 + 4*shat)*smT2 + &
!                       2*(4*MV2 - shat)*shat*smT2**2 - 2*smT2**3))*LC(15) &
!                  )/(E**(3*y)*shat*(MV2*shat - smT2)**3) - &
!               (2*smT*LC(7)*(MV2*(MV2 - 2*shat)*shat**2 + &
!                    shat*(-4*MV2 + shat)*smT2 + 2*smT2**2 + &
!                    shat*smT*(3*MV2*shat - smT2)*Cosh(y) + &
!                    shat*smT*(-(MV2*shat) + smT2)*Sinh(y))*SP(1))/ &
!                (E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
!               (2*(E**y*shat*smT*(2*MV2*shat - smT2) - &
!                    smT2*(-2*MV2*shat + shat**2 + 2*smT2))*LC(6)*SP(2))/ &
!                (E**y*shat**2*(MV2*shat*smT - smT**3)) + &
!               (LC(13)*(-2*E**(4*y)*shat**2*(MV2*shat - smT2)*smT2* &
!                     (MV2**2 + smT2)*SP(3) - &
!                    2*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
!                     (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2)*SP(3) &
!                      - 8*shat**2*smT2**3*SP(4) + &
!                    4*E**y*shat*smT**3*smT2*(3*MV2*shat + smT2)*SP(4) + &
!                    4*E**(2*y)*smT2**2* &
!                     (2*MV2*shat**3 - shat*(3*MV2 + 2*shat)*smT2 + &
!                       smT2**2)*SP(4) + &
!                    2*E**(3*y)*shat*smT*(MV2*shat - smT2)*smT2* &
!                     (MV2*shat*SP(3) - 2*(MV2*shat + smT2)*SP(4))))/ &
!                (E**(5*y)*shat**2*smT*(MV2*shat - smT2)**3) + &
!               (LC(9)*(-4*E**(4*y)*smT2* &
!                     (-2*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                       3*MV2*shat*smT2**2 + smT2**3)*SP(3) + &
!                    4*shat**2*smT2**2*(-2*MV2*shat + smT2)*SP(4) + &
!                    4*E**y*shat*smT*smT2* &
!                     (MV2**2*shat**2 + 2*MV2*shat*smT2 - smT2**2)*SP(4) &
!                     + 2*E**(3*y)*shat*smT* &
!                     (2*MV2*shat*smT2*(MV2*shat + 3*smT2)*SP(3) + &
!                       (MV2*shat - smT2)* &
!                        (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2)* &
!                        SP(4)) + &
!                    2*E**(2*y)*shat**2*smT2* &
!                     (-4*MV2*shat*smT2*SP(3) + &
!                       (2*MV2**2*(shat**2 - smT2) - 3*MV2*shat*smT2 + &
!                          smT2**2)*SP(4))))/ &
!                (E**(3*y)*shat**2*smT*(MV2*shat - smT2)**3) + &
!               (4*LC(3)*(-(E**(5*y)*shat*(MV2*shat - smT2)*smT2* &
!                       (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)* &
!                       SP(3)) + &
!                    E**(4*y)*smT* &
!                     (2*MV2**3*shat**5 + &
!                       4*MV2**2*(MV2 - shat)*shat**3*smT2 + &
!                       3*MV2*shat**2*(-2*MV2 + shat)*smT2**2 + &
!                       (6*MV2 - shat)*shat*smT2**3 - 2*smT2**4)*SP(3) - &
!                    2*MV2*shat**3*smT**5*SP(4) + &
!                    4*E**y*MV2**2*shat**3*smT2**2*SP(4) - &
!                    E**(3*y)*shat*smT2* &
!                     (2*MV2**2*shat**2*(MV2*shat + smT2)*SP(3) + &
!                       (MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2*SP(4)) &
!                     + E**(2*y)*smT*smT2* &
!                     (2*MV2**2*shat**4*SP(3) + &
!                       (2*MV2**2*shat**4 - &
!                          3*MV2*shat**2*(2*MV2 + shat)*smT2 + &
!                          shat*(6*MV2 + shat)*smT2**2 - 2*smT2**3)*SP(4)) &
!                    ))/(E**(4*y)*shat**3*smT*(MV2*shat - smT2)**3) - &
!               (2*(E**y*MV2*shat*(MV2 + shat)*smT + &
!                    smT2*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*(2*MV2**2*shat**2 - 4*MV2*shat*smT2 + &
!                       smT2**2))*LC(5)*SP(5))/ &
!                (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(E**y*MV2*shat*(5*MV2*shat - smT2)*smT2 + &
!                    2*shat*smT**3*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*shat*smT* &
!                     (-(MV2**2*(MV2 - 5*shat)*shat) - &
!                       MV2*(MV2 + 6*shat)*smT2 + smT2**2) + &
!                    E**(3*y)*(-2*MV2**3*shat**3 + MV2**2*shat**2*smT2 + &
!                       smT2**3))*LC(4)*SP(6))/ &
!                (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
!               (4*(E**y*MV2 - smT)*smT2* &
!                  (-2*shat*smT + E**y*(MV2*shat + smT2))*LC(10)*SP(6))/ &
!                (E**(2*y)*(-(MV2*shat) + smT2)**3) + &
!               (4*(-shat + E**y*smT)*smT2*LC(14)* &
!                  (-3*MV2*shat*smT + smT**3 - &
!                    2*shat*((MV2*shat - 2*smT2)*Cosh(y) + &
!                       MV2*shat*Sinh(y)))*SP(6))/ &
!                (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(2*shat*smT**5 + &
!                    E**(2*y)*shat*smT* &
!                     (2*MV2**3*shat + MV2*shat*smT2 - smT2**2) - &
!                    E**y*smT2*(3*MV2**2*shat**2 + smT2**2) + &
!                    E**(3*y)*(MV2*shat - smT2)* &
!                     (4*MV2**2*shat**2 - 7*MV2*shat*smT2 + 2*smT2**2))* &
!                  LC(2)*SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3) - &
!               (4*(E**y*MV2 - smT)*smT* &
!                  (-2*MV2*shat**2 + E**y*smT*(MV2*shat + smT2))*LC(8)* &
!                  SP(7))/(E**(2*y)*(-(MV2*shat) + smT2)**3) + &
!               (4*smT**3*(shat - E**y*smT)* &
!                  (2*shat*smT + E**y*(-3*MV2*shat + smT2))*LC(12)*SP(7))/ &
!                (E**(4*y)*shat*(MV2*shat - smT2)**3)) + &
!            (4*kapPr*(E**(2*y)*shat*smT**3*(MV2*shat - smT2)*SP(2)* &
!                  SP(6) + shat*smT*smT2*(SP(3) + SP(4))* &
!                  (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)) + &
!                 E**y*smT2*(-((MV2*shat - smT2)* &
!                       (smT2*(SP(3) + SP(4))*SP(5) + shat**2*SP(2)*SP(6)) &
!                       ) + 2*shat*(MV2*shat*SP(3) + smT2*SP(4))*SP(6)* &
!                     SP(7))))/(E**(2*y)*shat*smT*(-(MV2*shat) + smT2)**2) &
!            ) + SI(8)*(kap*((smT* &
!                  (12*MV2*shat**2*smT2 - &
!                    E**(2*y)*shat*(-(MV2*shat) + smT2)* &
!                     (MV2**2 - 5*MV2*shat + 5*smT2) - &
!                    E**(4*y)*shat*smT2*(5*MV2*shat + 7*smT2) + &
!                    2*E**(3*y)*smT* &
!                     (4*MV2**2*shat**2 + 7*MV2*shat*smT2 + smT2**2) + &
!                    E**y*smT*(-11*MV2**2*shat**2 - 16*MV2*shat*smT2 + &
!                       3*smT2**2))*LC(1))/(E**y*(-(MV2*shat) + smT2)**3) &
!                + (smT*(6*shat*smT2**2 + &
!                    E**(4*y)*shat*(MV2*shat - 3*smT2)* &
!                     (2*MV2*shat + smT2) + &
!                    E**(2*y)*shat*(-(MV2*shat) + smT2)* &
!                     (MV2*(MV2 + shat) + 2*smT2) - &
!                    2*E**y*smT* &
!                     (-(MV2**2*shat**2) + 6*MV2*shat*smT2 + smT2**2) - &
!                    E**(3*y)*smT* &
!                     (MV2**2*shat**2 - 15*MV2*shat*smT2 + 2*smT2**2))* &
!                  LC(11))/(E**(3*y)*(-(MV2*shat) + smT2)**3) + &
!               (shat*smT*(E**(2*y)*shat*smT*(MV2*shat - smT2) + &
!                    2*shat*smT*(2*MV2*shat + 3*smT2) - &
!                    2*E**(4*y)*shat*smT*(MV2*shat + 4*smT2) + &
!                    E**(3*y)*(3*MV2*shat + smT2)*(MV2*shat + 4*smT2) - &
!                    E**y*(5*MV2**2*shat**2 + 12*MV2*shat*smT2 + &
!                       3*smT2**2))*LC(15))/ &
!                (E**(2*y)*(MV2*shat - smT2)**3) + &
!               (2*(MV2*shat*smT - E**(3*y)*shat*smT2 - &
!                    E**y*shat*(2*MV2**2 + smT2) + &
!                    E**(2*y)*smT*(2*MV2*shat + smT2))*LC(7)*SP(1))/ &
!                (E**y*(-(MV2*shat) + smT2)**2) + &
!               (4*smT*LC(6)*Sinh(y)*SP(2))/(-(MV2*shat) + smT2) + &
!               (LC(9)*(-8*E**(6*y)*shat*smT2**2*SP(3) + &
!                    4*E**(5*y)*smT**3*(3*MV2*shat + smT2)*SP(3) + &
!                    4*E**(4*y)*shat*smT2* &
!                     (-2*smT2*SP(3) + MV2*shat*(2*SP(3) - SP(4))) + &
!                    4*shat*(2*MV2*shat - smT2)*smT2*SP(4) + &
!                    2*E**y*smT* &
!                     (-4*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) &
!                     - 2*E**(2*y)*shat* &
!                     (-4*MV2*shat*smT2*SP(3) + &
!                       (MV2*shat - smT2)*(MV2*(MV2 + shat) - smT2)*SP(4)) &
!                      + 4*E**(3*y)*smT* &
!                     (-3*MV2*shat*smT2*SP(3) + smT2**2*SP(3) + &
!                       2*MV2**2*shat**2*(-SP(3) + SP(4)))))/ &
!                (E**(2*y)*(MV2*shat - smT2)**3) + &
!               (4*LC(3)*(E**(6*y)*smT**3*(3*MV2*shat - smT2)*SP(3) - &
!                    4*E**(5*y)*MV2**2*shat*smT2*SP(3) + &
!                    2*MV2*shat*smT**3*SP(4) - &
!                    4*E**y*MV2**2*shat*smT2*SP(4) + &
!                    4*E**(3*y)*MV2**3*shat**2*(SP(3) + SP(4)) - &
!                    E**(2*y)*smT* &
!                     (-3*MV2*shat*smT2*SP(4) + smT2**2*SP(4) + &
!                       2*MV2**2*shat**2*(SP(3) + SP(4))) + &
!                    E**(4*y)*(-2*smT**5*SP(3) + &
!                       MV2*shat*smT**3*(7*SP(3) + SP(4)) - &
!                       MV2**2*shat**2*smT*(5*SP(3) + 3*SP(4)))))/ &
!                (E**(3*y)*(MV2*shat - smT2)**3) + &
!               (LC(13)*(-2*E**(5*y)*MV2*shat*smT*(MV2*shat - smT2)* &
!                     SP(3) - 8*E**(2*y)*shat*(MV2*shat - smT2)*smT2* &
!                     SP(4) + 8*shat*smT2**2*SP(4) - &
!                    4*E**y*smT**3*(3*MV2*shat + smT2)*SP(4) + &
!                    2*E**(4*y)*MV2*shat* &
!                     ((MV2 + shat)*(MV2*shat - smT2)*SP(3) - &
!                       4*shat*smT2*SP(4)) - &
!                    2*E**(3*y)*smT* &
!                     (MV2**2*shat**2*(SP(3) - 4*SP(4)) + &
!                       2*smT2**2*SP(4) - MV2*shat*smT2*(SP(3) + 6*SP(4))) &
!                    ))/(E**(4*y)*(MV2*shat - smT2)**3) + &
!               (2*(-2*MV2*shat*smT + smT**3 + E**(3*y)*shat*smT2 + &
!                    E**y*shat*(2*MV2**2 + smT2) - &
!                    E**(2*y)*smT*(MV2*shat + 2*smT2))*LC(5)*SP(5))/ &
!                (E**y*(-(MV2*shat) + smT2)**2) - &
!               (2*(E**(3*y)*(-10*MV2**2*shat**2*smT + &
!                       5*MV2*shat*smT**3 + smT**5) + &
!                    2*E**y*MV2*shat*smT*(3*MV2*shat - smT2) + &
!                    2*E**(4*y)*shat*(3*MV2*shat - 2*smT2)*smT2 + &
!                    2*shat*smT2*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*shat*(MV2*shat - smT2)*(2*MV2*shat + smT2))* &
!                  LC(4)*SP(6))/(E**(2*y)*(MV2*shat - smT2)**3) - &
!               (4*shat*(2*E**(4*y)*MV2*shat**2*smT - 2*shat*smT**3 + &
!                    2*E**(2*y)*shat*smT*(MV2*shat - smT2) + &
!                    E**y*smT2*(3*MV2*shat + smT2) + &
!                    E**(3*y)*(-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + &
!                       smT2**2))*LC(14)*SP(6))/ &
!                (E**(3*y)*(MV2*shat - smT2)**3) + &
!               (8*E**y*shat*smT2* &
!                  (3*MV2*shat + smT2 - 4*shat*smT*Cosh(y))*LC(10)* &
!                  Sinh(y)*SP(6))/(MV2*shat - smT2)**3 + &
!               (2*(2*E**(2*y)*MV2*shat*(-MV2 + shat)*(MV2*shat - smT2) - &
!                    2*shat*smT2**2 + &
!                    E**(4*y)*shat*smT2*(MV2*shat + smT2) - &
!                    E**(3*y)*smT*(MV2*shat + smT2)**2 + &
!                    E**y*smT**3*(3*MV2*shat + smT2))*LC(2)*SP(7))/ &
!                (E**(2*y)*(MV2*shat - smT2)**3) - &
!               (4*shat*(2*MV2*shat**2*smT - 2*E**(4*y)*shat*smT**3 + &
!                    2*E**(2*y)*shat*smT*(MV2*shat - smT2) + &
!                    E**(3*y)*smT2*(3*MV2*shat + smT2) + &
!                    E**y*(-2*MV2**2*shat**2 - 3*MV2*shat*smT2 + smT2**2)) &
!                   *LC(8)*SP(7))/(E**y*(MV2*shat - smT2)**3) - &
!               (8*shat*smT2*(3*MV2*shat + smT2 - 4*shat*smT*Cosh(y))* &
!                  LC(12)*Sinh(y)*SP(7))/(E**y*(MV2*shat - smT2)**3)) + &
!            (4*kapPr*shat*(smT*(-(MV2*shat) + smT2)*(SP(3) + SP(4))* &
!                  SP(5) + smT* &
!                  ((-(MV2*shat) + smT2)*SP(1) + &
!                    2*shat*(SP(3) + SP(4))*SP(6))*SP(7) + &
!                 E**(2*y)*smT* &
!                  (-((MV2*shat - smT2)* &
!                       ((SP(3) + SP(4))*SP(5) + SP(2)*SP(6))) + &
!                    2*shat*(SP(3) + SP(4))*SP(6)*SP(7)) + &
!                 E**y*((MV2*shat - smT2)* &
!                     (2*MV2*(SP(3) + SP(4))*SP(5) + shat*SP(2)*SP(6)) + &
!                    (shat*(MV2*shat - smT2)*SP(1) - &
!                       2*(MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7))) &
!               )/(E**y*(-(MV2*shat) + smT2)**2)) + &
!         SI(15)*(kap*((smT*(-(E**(3*y)*shat*smT2*(11*MV2*shat + smT2)* &
!                       (MV2*shat + 2*smT2)) + &
!                    E**(4*y)*shat**2*smT**3*(5*MV2*shat + 7*smT2) + &
!                    E**(2*y)*shat*smT* &
!                     (MV2**2*shat**2*(-14*MT2 + 7*MV2 + 2*shat) + &
!                       MV2*(-8*MT2 + 25*MV2 - 6*shat)*shat*smT2 + &
!                       2*(11*MT2 + 2*(MV2 + shat))*smT2**2) - &
!                    E**y*(MV2**3*shat**3*(-4*MT2 + MV2 + 2*shat) + &
!                       MV2**2*(-26*MT2 + 9*MV2 - 6*shat)*shat**2*smT2 + &
!                       2*MV2*shat*(14*MT2 + MV2 + shat)*smT2**2 + &
!                       2*(MT2 + shat)*smT2**3) + &
!                    2*shat*smT*(-(MV2*shat) + smT2)* &
!                     (MV2*smT2 + 2*MT2*(-(MV2*shat) + smT2)))*LC(1))/ &
!                (shat*(MV2*shat - smT2)**3) + &
!               ((-(E**(4*y)*shat**2*(MV2*shat - 3*smT2)*smT2* &
!                       (2*MV2*shat + smT2)) + &
!                    shat*(MV2*shat - smT2)*smT2* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
!                    2*E**(3*y)*MV2*shat**2*smT* &
!                     (2*MV2**2*shat**2 - 4*MV2*shat*smT2 - 7*smT2**2) - &
!                    2*E**y*smT* &
!                     (MV2**3*shat**4 + &
!                       2*MV2**2*(-3*MT2 + MV2 - shat)*shat**2*smT2 + &
!                       MV2*shat*(3*MT2 + MV2 + 2*shat)*smT2**2 + &
!                       (3*MT2 - shat)*smT2**3) + &
!                    E**(2*y)*shat* &
!                     (2*(4*MT2 - MV2)*MV2**3*shat**3 + &
!                       3*MV2**2*(-8*MT2 + MV2)*shat**2*smT2 + &
!                       MV2*shat*(6*MT2 + 15*MV2 + shat)*smT2**2 + &
!                       (10*MT2 + 2*MV2 - shat)*smT2**3))*LC(11))/ &
!                (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               ((-2*E**(4*y)*shat**2*smT**3*(MV2*shat + 4*smT2) + &
!                    E**(3*y)*shat*smT2*(5*MV2*shat + smT2)* &
!                     (MV2*shat + 4*smT2) + &
!                    shat*smT*(MV2*shat - smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
!                    E**(2*y)*shat*smT* &
!                     (2*(3*MT2 - 2*MV2)*MV2**2*shat**2 + &
!                       MV2*shat*(18*MT2 - 18*MV2 + shat)*smT2 - &
!                       (24*MT2 + 8*MV2 + shat)*smT2**2) + &
!                    E**y*(MV2**3*(-4*MT2 + MV2)*shat**3 + &
!                       MV2**2*(-14*MT2 + 5*MV2 - 2*shat)*shat**2*smT2 + &
!                       2*MV2*shat*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
!                       8*MT2*smT2**3))*LC(15))/ &
!                (E**y*(MV2*shat - smT2)**3) + &
!               (2*(shat**2*(MV2 - E**y*smT)**3 + &
!                    2*MT2*(-2*MV2*shat + smT*(E**y*shat + smT))* &
!                     (MV2*shat - smT2))*LC(7)*SP(1))/ &
!                (shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(MV2 - E**y*smT)*(shat + E**y*smT)*LC(6)*SP(2))/ &
!                (MV2*shat - smT2) + &
!               (LC(13)*(-2*E**(4*y)*shat*(-(MV2*shat*smT) + smT**3)**2* &
!                     SP(3) - 4*E**y*shat*smT* &
!                     (MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                       2*MV2*(MT2 + 2*MV2)*shat*smT2 + &
!                       (2*MT2 + MV2)*smT2**2)*SP(4) + &
!                    4*smT2*(MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(4) &
!                      + 2*E**(3*y)*shat*smT* &
!                     ((MV2 + shat)*(-(MV2*shat) + smT2)**2*SP(3) - &
!                       4*MV2*shat**2*smT2*SP(4)) + &
!                    2*E**(2*y)*shat* &
!                     ((MV2*shat - 2*smT2)*(-(MV2*shat) + smT2)**2* &
!                        SP(3) + 6*MV2*shat*smT2*(MV2*shat + smT2)*SP(4))) &
!                  )/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*LC(9)*(-4*E**(6*y)*shat**2*smT**5*SP(3) + &
!                    2*E**(5*y)*shat*smT2**2*(5*MV2*shat + smT2)*SP(3) + &
!                    shat*smT*(MV2*shat - smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2)*SP(4) + &
!                    E**y*MV2*shat* &
!                     (shat*(MV2**2 - MV2*shat + smT2)* &
!                        (MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(4) &
!                      + 2*E**(4*y)*shat*smT**3* &
!                     (2*(MV2*shat*(3*MT2 - 2*MV2 + shat) - &
!                          (3*MT2 + MV2 + shat)*smT2)*SP(3) - &
!                       MV2*shat**2*SP(4)) + &
!                    E**(3*y)*smT2* &
!                     (2*(MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
!                          MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
!                          (2*MT2 + shat)*smT2**2)*SP(3) + &
!                       shat*(6*MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)* &
!                        SP(4)) - &
!                    E**(2*y)*shat*smT* &
!                     (-2*(MV2*shat - smT2)* &
!                        (4*MT2*(-(MV2*shat) + smT2) + &
!                          MV2*(MV2*shat + smT2))*SP(3) + &
!                       MV2*(-(MV2*shat**2*(6*MT2 - 5*MV2 + shat)) + &
!                          shat*(6*MT2 + shat)*smT2 + smT2**2)*SP(4))))/ &
!                (E**y*shat*(MV2*shat - smT2)**3) + &
!               (4*LC(3)*(-2*E**(4*y)*MV2*shat*smT**3*(4*MV2*shat - smT2)* &
!                     SP(3) + E**(5*y)*shat*(3*MV2*shat - smT2)*smT2**2* &
!                     SP(3) + MV2*smT* &
!                     (MV2**2*shat**2*(-8*MT2 + 2*MV2 + shat) - &
!                       2*MV2*shat*(-5*MT2 + shat)*smT2 + &
!                       (-2*MT2 + shat)*smT2**2)*SP(4) + &
!                    E**(2*y)*MV2*smT* &
!                     ((MV2**2*shat**2*(8*MT2 - 2*MV2 + 5*shat) - &
!                          2*MV2*shat*(5*MT2 + 3*shat)*smT2 + &
!                          (2*MT2 + shat)*smT2**2)*SP(3) + &
!                       6*MV2**2*shat**3*SP(4)) - &
!                    E**(3*y)*smT2* &
!                     ((MV2**2*shat**2*(10*MT2 - 7*MV2 + 3*shat) + &
!                          MV2*(-14*MT2 + MV2 - 4*shat)*shat*smT2 + &
!                          (4*MT2 + shat)*smT2**2)*SP(3) + &
!                       MV2*shat**2*(3*MV2*shat - smT2)*SP(4)) + &
!                    E**y*(-2*(MV2*shat - smT2)* &
!                        (MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                          6*MT2*MV2*shat*smT2 - 2*MT2*smT2**2)*SP(3) + &
!                       (MV2**3*(12*MT2 - 3*MV2 - 2*shat)*shat**3 + &
!                          MV2**2*shat**2*(-22*MT2 - 3*MV2 + 5*shat)* &
!                           smT2 + 2*MV2*(7*MT2 - 2*shat)*shat*smT2**2 + &
!                          (-4*MT2 + shat)*smT2**3)*SP(4))))/ &
!                (E**y*shat*(MV2*shat - smT2)**3) + &
!               (2*LC(5)*(-2*MT2*(-2*MV2*shat + smT*(E**y*shat + smT))* &
!                     (MV2*shat - smT2) + &
!                    shat*(-(MV2**3*shat) - &
!                       MV2**2*shat*(shat - 2*E**y*smT) + &
!                       smT**3* &
!                        (-2*smT + E**y*(shat + E**(2*y)*shat - E**y*smT)) &
!                         + E**y*MV2*smT* &
!                        (-shat**2 + smT2 + &
!                          shat*smT*(Cosh(y) - 5*Sinh(y)))))*SP(5))/ &
!                (shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(-2*E**(4*y)*shat*smT**3*(3*MV2*shat - 2*smT2) + &
!                    smT*(MV2*shat - smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) - &
!                    E**(3*y)*smT2* &
!                     (-12*MV2**2*shat**2 + 5*MV2*shat*smT2 + smT2**2) + &
!                    E**(2*y)*smT* &
!                     (MV2**2*shat**2*(14*MT2 - 7*MV2 + shat) - &
!                       MV2*shat*(22*MT2 + shat)*smT2 + &
!                       (8*MT2 + MV2)*smT2**2) + &
!                    E**y*MV2*(shat*(MV2**2 - MV2*shat + smT2)* &
!                        (MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
!                  LC(4)*SP(6))/(E**y*(MV2*shat - smT2)**3) + &
!               (4*E**y*smT*(-(shat*(MV2 - E**y*smT)**2* &
!                       (MV2*shat + smT*(-2*E**y*shat + smT))) + &
!                    2*MT2*(2*MV2*shat + smT*(-3*E**y*shat + smT))* &
!                     (MV2*shat - smT2))*LC(10)*SP(6))/ &
!                (-(MV2*shat) + smT2)**3 - &
!               (4*(E**y*shat - smT)* &
!                  (MV2*shat*(MV2 - E**y*smT)* &
!                     (MV2*shat + smT*(-2*E**y*shat + smT)) - &
!                    2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*LC(14)* &
!                  SP(6))/(E**y*(MV2*shat - smT2)**3) + &
!               (2*(E**(4*y)*shat**2*smT**3*(MV2*shat + smT2) - &
!                    2*E**(3*y)*MV2*shat**2*smT2*(MV2*shat + 2*smT2) + &
!                    3*shat*smT*(MV2*shat - smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
!                    E**(2*y)*shat*smT* &
!                     (MV2**2*shat**2*(-4*MT2 + MV2 + 2*shat) + &
!                       MV2*(2*MT2 + 5*MV2 - shat)*shat*smT2 + &
!                       (2*MT2 - shat)*smT2**2) - &
!                    2*E**y*(MV2**3*shat**4 + &
!                       MV2**2*shat**2*(-4*MT2 + MV2 + shat)*smT2 + &
!                       MV2*(5*MT2 - 2*shat)*shat*smT2**2 - MT2*smT2**3))* &
!                  LC(2)*SP(7))/(E**y*shat*(MV2*shat - smT2)**3) + &
!               (4*(shat*(-MV2 + E**y*smT)* &
!                     (MV2*shat + smT*(-2*E**y*shat + smT))* &
!                     (MV2*(shat + E**y*smT) - (1 + E**(2*y))*smT2) + &
!                    2*MT2*(MV2*shat - smT2)* &
!                     (2*MV2*shat*(shat + E**y*smT) + &
!                       (-((2 + 3*E**(2*y))*shat) + E**y*smT)*smT2))* &
!                  LC(8)*SP(7))/(MV2*shat - smT2)**3 + &
!               (4*smT*(shat*(MV2 - E**y*smT)**2* &
!                     (MV2*shat + smT*(-2*E**y*shat + smT)) - &
!                    2*MT2*(2*MV2*shat + smT*(-3*E**y*shat + smT))* &
!                     (MV2*shat - smT2))*LC(12)*SP(7))/ &
!                (E**y*(-(MV2*shat) + smT2)**3)) - &
!            (2*kapPr*(smT*(-(MV2*shat) + smT2)* &
!                  ((MV2*shat - smT2)*SP(4)*SP(5) + MV2*shat*SP(1)*SP(7)) &
!                  + E**(3*y)*shat*smT2* &
!                  ((MV2*shat - smT2)* &
!                     (2*(SP(3) + SP(4))*SP(5) + SP(2)*SP(6)) - &
!                    4*shat*(SP(3) + SP(4))*SP(6)*SP(7)) + &
!                 E**(2*y)*smT* &
!                  (-((MV2*shat - smT2)* &
!                       (((3*MV2*shat + smT2)*SP(3) + 4*MV2*shat*SP(4))* &
!                          SP(5) + shat*(MV2 + shat)*SP(2)*SP(6))) + &
!                    shat*(shat*(-(MV2*shat) + smT2)*SP(1) + &
!                       2*(3*MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7) &
!                    ) - E**y*((MV2*shat - smT2)* &
!                     (8*MT2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) - &
!                       2*shat* &
!                        (MV2**2*SP(3) + (MV2*(MV2 + shat) - smT2)*SP(4))* &
!                        SP(5) - shat*smT2*SP(2)*SP(6)) + &
!                    2*shat*(MV2*shat*(-(MV2*shat) + smT2)*SP(1) + &
!                       (4*MT2*(-(MV2*shat) + smT2) + &
!                          MV2*(MV2*shat + smT2))*(SP(3) + SP(4))*SP(6))* &
!                     SP(7))))/(E**y*(-(MV2*shat) + smT2)**2)) + &
!         SI(14)*(kap*((smT*(-12*MV2*shat**3*smT**3 + &
!                    2*E**(4*y)*shat*smT*(MV2*shat - smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) + &
!                    E**y*shat*smT2* &
!                     (23*MV2**2*shat**2 + 14*MV2*shat*smT2 - smT2**2) + &
!                    E**(2*y)*shat*smT* &
!                     (MV2**2*shat**2*(28*MT2 - 13*MV2 + 3*shat) - &
!                       2*MV2*shat*(10*MT2 + 11*MV2 + 2*shat)*smT2 + &
!                       (-8*MT2 - MV2 + shat)*smT2**2) + &
!                    E**(3*y)*(-(MV2**3*shat**3*(8*MT2 - 2*MV2 + shat)) - &
!                       4*MV2**2*shat**2*(4*MT2 - 2*MV2 + shat)*smT2 + &
!                       MV2*shat*(20*MT2 + 2*MV2 + 7*shat)*smT2**2 + &
!                       2*(2*MT2 - shat)*smT2**3))*LC(1))/ &
!                (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               ((-6*shat**2*smT2**3 - &
!                    E**y*shat*smT**3* &
!                     (MV2**2*shat**2 - 17*MV2*shat*smT2 - 2*smT2**2) + &
!                    2*E**(4*y)*shat*(MV2*shat - smT2)* &
!                     (MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                       (6*MT2 - MV2)*MV2*shat*smT2 + &
!                       (-2*MT2 + MV2)*smT2**2) + &
!                    E**(2*y)*shat*smT2* &
!                     (MV2**2*shat**2*(-2*MT2 + 2*MV2 + shat) + &
!                       2*(11*MT2 - 8*MV2)*MV2*shat*smT2 - &
!                       (20*MT2 + 4*MV2 + shat)*smT2**2) + &
!                    E**(3*y)*smT* &
!                     (-(MV2**3*shat**3*(-4*MT2 + MV2 + 3*shat)) + &
!                       MV2**2*shat**2*(-22*MT2 + 5*MV2 + 4*shat)*smT2 + &
!                       MV2*(14*MT2 + 2*MV2 - 3*shat)*shat*smT2**2 + &
!                       2*(2*MT2 + shat)*smT2**3))*LC(11))/ &
!                (E**(4*y)*shat*(MV2*shat - smT2)**3) + &
!               ((2*shat**2*smT**3*(2*MV2*shat + 3*smT2) - &
!                    E**y*shat*smT2* &
!                     (8*MV2**2*shat**2 + 19*MV2*shat*smT2 + 3*smT2**2) + &
!                    E**(2*y)*shat*smT* &
!                     (5*MV2**2*(-2*MT2 + MV2)*shat**2 + &
!                       2*MV2*(-5*MT2 + 9*MV2)*shat*smT2 + &
!                       (20*MT2 + 7*MV2)*smT2**2) - &
!                    E**(3*y)*(MV2*shat + 4*smT2)* &
!                     (MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
!                  LC(15))/(E**(3*y)*(MV2*shat - smT2)**3) + &
!               (2*(MV2*shat**2*smT2 - &
!                    E**y*shat**2*smT*(2*MV2**2 + smT2) + &
!                    E**(3*y)*shat*smT* &
!                     (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2) + &
!                    E**(2*y)*(MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                       MV2*shat*(6*MT2 + shat)*smT2 + &
!                       (-2*MT2 + shat)*smT2**2))*LC(7)*SP(1))/ &
!                (E**(2*y)*shat*(-(MV2*shat) + smT2)**2) - &
!               (2*(E**y*MV2 - smT)*(E**y*shat + smT)*LC(6)*SP(2))/ &
!                (E**(2*y)*(MV2*shat - smT2)) + &
!               (4*LC(3)*(E**(6*y)*smT2*(-(MV2*shat) + smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2)*SP(3) + &
!                    2*E**(5*y)*MV2*smT* &
!                     (MV2**2*shat**2*(-4*MT2 + MV2 + shat) + &
!                       MV2*(5*MT2 - shat)*shat*smT2 - MT2*smT2**2)*SP(3) &
!                     - 6*E**y*MV2**2*shat**2*smT**3*SP(4) + &
!                    2*MV2*shat**2*smT2**2*SP(4) + &
!                    E**(2*y)*smT2* &
!                     (-2*MV2**2*shat**3*SP(3) + &
!                       (-(MV2**2*shat**2*(8*MT2 - 6*MV2 + shat)) + &
!                          MV2*shat*(10*MT2 + shat)*smT2 - 2*MT2*smT2**2)* &
!                        SP(4)) + &
!                    2*E**(3*y)*MV2*smT* &
!                     (MV2*shat**2*(2*MV2*shat + smT2)*SP(3) + &
!                       (MV2**2*shat**2*(4*MT2 - MV2 + shat) - &
!                          MV2*shat*(5*MT2 + shat)*smT2 + MT2*smT2**2)* &
!                        SP(4)) + &
!                    E**(4*y)*((-2*MV2**3*shat**3*(-4*MT2 + MV2 + shat) - &
!                          4*MV2**2*(3*MT2 + MV2 - shat)*shat**2*smT2 - &
!                          3*MV2*shat*(-2*MT2 + shat)*smT2**2 + &
!                          (-2*MT2 + shat)*smT2**3)*SP(3) + &
!                       (MV2*shat - smT2)* &
!                        ((4*MT2 - MV2)*MV2**2*shat**2 - &
!                          6*MT2*MV2*shat*smT2 + 2*MT2*smT2**2)*SP(4))))/ &
!                (E**(4*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*smT*LC(9)*(-2*E**(5*y)*smT* &
!                     (MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2))*SP(3) &
!                      - 2*shat**2*smT2*(-2*MV2*shat + smT2)*SP(4) + &
!                    E**y*shat*smT*(-7*MV2**2*shat**2 + smT2**2)*SP(4) + &
!                    E**(4*y)*shat* &
!                     (2*(MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                          2*MV2*(MT2 + 2*MV2)*shat*smT2 + &
!                          (2*MT2 + MV2)*smT2**2)*SP(3) + &
!                       (MV2*shat - smT2)* &
!                        (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2)*SP(4)) + &
!                    E**(3*y)*smT* &
!                     (-6*MV2*shat**2*(MV2*shat + smT2)*SP(3) + &
!                       (MV2**2*shat**2*(8*MT2 - 2*MV2 + 3*shat) - &
!                          2*MV2*shat*(5*MT2 + 2*shat)*smT2 + &
!                          (2*MT2 + shat)*smT2**2)*SP(4)) - &
!                    E**(2*y)*shat* &
!                     (-4*MV2*shat**2*smT2*SP(3) + &
!                       (MV2**2*shat**2*(10*MT2 - 3*MV2 + 2*shat) - &
!                          MV2*shat*(14*MT2 + 4*MV2 + 3*shat)*smT2 + &
!                          (4*MT2 + MV2 + shat)*smT2**2)*SP(4))))/ &
!                (E**(3*y)*shat*(MV2*shat - smT2)**3) + &
!               (LC(13)*(-2*E**(6*y)*shat*smT*(MV2*shat - smT2)* &
!                     (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2)*SP(3) + &
!                    2*E**(5*y)*(-(MV2*shat) + smT2)* &
!                     (MV2**2*(-4*MT2 + MV2)*shat**2 + &
!                       MV2*shat*(6*MT2 + shat)*smT2 + &
!                       (-2*MT2 + shat)*smT2**2)*SP(3) + &
!                    8*shat**2*smT**5*SP(4) - &
!                    4*E**y*shat*smT2**2*(5*MV2*shat + smT2)*SP(4) - &
!                    8*E**(2*y)*shat*smT**3* &
!                     (MV2*shat*(3*MT2 - 2*MV2 + shat) - &
!                       (3*MT2 + MV2 + shat)*smT2)*SP(4) - &
!                    2*E**(3*y)*smT2* &
!                     (MV2*shat**2*(MV2*shat - smT2)*SP(3) + &
!                       2*(MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
!                          MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
!                          (2*MT2 + shat)*smT2**2)*SP(4)) + &
!                    2*E**(4*y)*shat*smT*(MV2*shat - smT2)* &
!                     (shat*(2*MV2**2 + smT2)*SP(3) - &
!                       2*(4*MT2*(-(MV2*shat) + smT2) + &
!                          MV2*(MV2*shat + smT2))*SP(4))))/ &
!                (E**(5*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(E**y*MV2*shat*smT*(3*MV2*shat + shat**2 - smT2) + &
!                    shat*smT2*(-2*MV2*shat + smT2) + &
!                    E**(2*y)*(MV2**2*shat**2*(4*MT2 - MV2 + shat) - &
!                       2*MV2*shat*(3*MT2 + 2*shat)*smT2 + &
!                       (2*MT2 + shat)*smT2**2) + &
!                    E**(3*y)*shat*smT* &
!                     (MV2*smT2 + 2*MT2*(-(MV2*shat) + smT2)))*LC(5)*SP(5) &
!                  )/(E**(2*y)*shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(E**(2*y)*(MV2**2*(-10*MT2 + 6*MV2 - 3*shat)*shat**2* &
!                        smT + 2*MV2*shat*(7*MT2 + shat)*smT**3 + &
!                       (-4*MT2 + shat)*smT**5) + &
!                    2*shat*smT**3*(2*MV2*shat - smT2) - &
!                    3*E**y*MV2*shat*(3*MV2*shat - smT2)*smT2 + &
!                    4*E**(4*y)*smT*(-(MV2*shat) + smT2)* &
!                     (MV2*(-2*MT2 + MV2)*shat + 2*MT2*smT2) - &
!                    E**(3*y)*MV2* &
!                     (MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
!                       MV2*(2*MT2 + MV2 - 2*shat)*shat*smT2 + &
!                       (2*MT2 + 5*shat)*smT2**2))*LC(4)*SP(6))/ &
!                (E**(3*y)*(MV2*shat - smT2)**3) + &
!               (4*smT*(-2*shat**2*smT**3 + &
!                    E**y*shat*smT2*(5*MV2*shat + smT2) - &
!                    2*E**(2*y)*shat*smT* &
!                     (3*MT2*(-(MV2*shat) + smT2) + &
!                       MV2*(2*MV2*shat + smT2)) + &
!                    E**(3*y)*(MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
!                  LC(10)*SP(6))/(E**(2*y)*(-(MV2*shat) + smT2)**3) + &
!               (4*(2*shat**2*smT2**2 - &
!                    E**y*shat*smT**3*(5*MV2*shat + smT2) - &
!                    2*E**(2*y)*shat*smT2* &
!                     (MV2*shat*(3*MT2 - 2*MV2 + shat) - &
!                       (3*MT2 + MV2 + shat)*smT2) - &
!                    E**(3*y)*smT* &
!                     (MV2**2*(-4*MT2 + MV2 - 3*shat)*shat**2 + &
!                       MV2*shat*(2*MT2 + MV2 + 2*shat)*smT2 + &
!                       (2*MT2 + shat)*smT2**2) - &
!                    E**(4*y)*shat*(MV2*shat - smT2)* &
!                     (4*MT2*(-(MV2*shat) + smT2) + MV2*(MV2*shat + smT2)) &
!                    )*LC(14)*SP(6))/(E**(4*y)*(MV2*shat - smT2)**3) + &
!               (2*(-2*shat**2*smT**5 + &
!                    E**(4*y)*shat*smT*(MV2*shat - smT2)* &
!                     (2*MT2*MV2*shat - (2*MT2 + MV2)*smT2) + &
!                    E**y*shat*smT2* &
!                     (3*MV2**2*shat**2 + 2*MV2*shat*smT2 + smT2**2) - &
!                    E**(2*y)*shat*smT* &
!                     (MV2**2*(-2*MT2 + 5*MV2)*shat**2 + &
!                       MV2*shat*(-2*MT2 + shat)*smT2 + &
!                       (4*MT2 + MV2 - shat)*smT2**2) + &
!                    E**(3*y)*(2*MV2**3*(-4*MT2 + MV2)*shat**3 + &
!                       MV2**2*shat**2*(12*MT2 + shat)*smT2 - &
!                       6*MT2*MV2*shat*smT2**2 + (2*MT2 - shat)*smT2**3))* &
!                  LC(2)*SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3) + &
!               (4*(-shat + E**y*smT)* &
!                  (2*MV2*shat**2*smT2 - &
!                    E**y*MV2*shat*smT*(3*MV2*shat + smT2) + &
!                    E**(2*y)*(MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
!                  LC(8)*SP(7))/(E**(2*y)*(MV2*shat - smT2)**3) + &
!               (4*smT*(-2*shat**2*smT**3 + &
!                    E**y*shat*smT2*(5*MV2*shat + smT2) - &
!                    2*E**(2*y)*shat*smT* &
!                     (3*MT2*(-(MV2*shat) + smT2) + &
!                       MV2*(2*MV2*shat + smT2)) + &
!                    E**(3*y)*(MV2**2*shat*(MV2*shat + smT2) - &
!                       2*MT2*(MV2*shat - smT2)*(2*MV2*shat + smT2)))* &
!                  LC(12)*SP(7))/(E**(4*y)*(MV2*shat - smT2)**3)) + &
!            (2*kapPr*(E**(3*y)*smT*(MV2*shat - smT2)* &
!                  ((MV2*shat - smT2)*SP(3)*SP(5) + MV2*shat*SP(2)*SP(6)) &
!                  + shat*smT2* &
!                  (-2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) + &
!                    ((-(MV2*shat) + smT2)*SP(1) + &
!                       4*shat*(SP(3) + SP(4))*SP(6))*SP(7)) + &
!                 E**y*smT*((MV2*shat - smT2)* &
!                     ((4*MV2*shat*SP(3) + (3*MV2*shat + smT2)*SP(4))* &
!                        SP(5) + shat**2*SP(2)*SP(6)) + &
!                    shat*((MV2 + shat)*(MV2*shat - smT2)*SP(1) - &
!                       2*(3*MV2*shat + smT2)*(SP(3) + SP(4))*SP(6))*SP(7) &
!                    ) + E**(2*y)* &
!                  (2*(MV2*shat - smT2)* &
!                     (4*MT2*(MV2*shat - smT2)*(SP(3) + SP(4))*SP(5) + &
!                       shat*(-(((MV2*(MV2 + shat) - smT2)*SP(3) + &
!                               MV2**2*SP(4))*SP(5)) - &
!                          MV2*shat*SP(2)*SP(6))) + &
!                    shat*(smT2*(-(MV2*shat) + smT2)*SP(1) + &
!                       2*(4*MT2*(-(MV2*shat) + smT2) + &
!                          MV2*(MV2*shat + smT2))*(SP(3) + SP(4))*SP(6))* &
!                     SP(7))))/(E**(2*y)*(-(MV2*shat) + smT2)**2)) + &
!         SI(11)*(kap*(((12*MV2*shat**4*smT**3 - &
!                    2*E**(4*y)*shat*smT**5*(MV2*shat - smT2) - &
!                    E**(2*y)*shat*smT* &
!                     (MV2**2*shat**3*(8*MT2 + 10*MV2 + shat) + &
!                       (8*MT2 - 57*MV2)*MV2*shat**2*smT2 - &
!                       shat*(16*MT2 - 24*MV2 + shat)*smT2**2 - 13*smT2**3 &
!                       ) + E**(3*y)*smT2* &
!                     (MV2**2*shat**3*(8*MT2 + 4*MV2 + shat) + &
!                       2*MV2*shat**2*(4*MT2 - 14*MV2 + shat)*smT2 + &
!                       (-16*MT2 + 20*MV2 - 3*shat)*shat*smT2**2 - &
!                       8*smT2**3) + &
!                    E**y*shat**2* &
!                     (4*MV2**3*shat**3 - 23*MV2**2*shat**2*smT2 - &
!                       14*MV2*shat*smT2**2 - 3*smT2**3))*LC(1))/ &
!                (E**y*shat**2*(MV2*shat - smT2)**3) + &
!               ((-6*shat**4*smT**5 + &
!                    2*E**(5*y)*shat*(MV2*shat - 2*smT2)* &
!                     (MV2*shat - smT2)*(2*MV2*shat - smT2)*smT2 + &
!                    3*E**y*shat**3*smT2**2*(3*MV2*shat + 5*smT2) - &
!                    E**(2*y)*shat**2*smT* &
!                     (MV2**2*(4*MT2 + 11*MV2 - shat)*shat**3 - &
!                       5*MV2*(4*MT2 + 5*MV2)*shat**2*smT2 + &
!                       shat*(16*MT2 + 44*MV2 + shat)*smT2**2 + 6*smT2**3) &
!                      + E**(3*y)*shat* &
!                     (4*MV2**3*(8*MT2 + 2*MV2 - shat)*shat**4 + &
!                       8*MV2**2*shat**3*(-11*MT2 + MV2 + shat)*smT2 - &
!                       2*MV2*shat**2*(-28*MT2 + 28*MV2 + 5*shat)* &
!                        smT2**2 + shat*(83*MV2 + 6*shat)*smT2**3 - &
!                       19*smT2**4) - &
!                    E**(4*y)*smT* &
!                     (16*MV2**3*(2*MT2 + MV2)*shat**4 - &
!                       MV2**2*shat**3*(92*MT2 + 5*(7*MV2 + shat))*smT2 + &
!                       MV2*shat**2*(76*MT2 + 17*MV2 + 4*shat)*smT2**2 + &
!                       shat*(-16*MT2 + 16*MV2 + shat)*smT2**3 - 8*smT2**4 &
!                       ))*LC(11))/ &
!                (E**(3*y)*shat**2*(-shat + E**y*smT)* &
!                  (MV2*shat - smT2)**3) + &
!               ((5*E**y*shat**2*smT*(MV2*shat + smT2)* &
!                     (MV2*shat + 2*smT2) - &
!                    2*shat**3*smT2*(2*MV2*shat + 3*smT2) + &
!                    E**(3*y)*smT* &
!                     (MV2**2*shat**3*(-4*MT2 + 5*MV2 + 2*shat) - &
!                       MV2*shat**2*(12*MT2 + 5*MV2 + 4*shat)*smT2 + &
!                       2*shat*(8*MT2 + 6*MV2 + shat)*smT2**2 - 2*smT2**3) &
!                      - E**(2*y)*shat* &
!                     (MV2**2*shat**3*(-4*MT2 + MV2 + 2*shat) - &
!                       4*MV2*shat**2*(3*MT2 - 3*MV2 + shat)*smT2 + &
!                       shat*(16*MT2 + 11*MV2 + 2*shat)*smT2**2 + &
!                       6*smT2**3))*LC(15))/ &
!                (E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(-(shat**2*smT**3) - &
!                    E**(3*y)*shat*(MV2*shat - 2*smT2)*smT2 + &
!                    E**(2*y)*smT* &
!                     (MV2*shat**2*(MV2 + 3*shat) - &
!                       shat*(4*MV2 + 5*shat)*smT2 + 2*smT2**2) - &
!                    E**y*shat* &
!                     (MV2*shat**2*(3*MV2 + 2*shat) - &
!                       shat*(8*MV2 + 3*shat)*smT2 + 3*smT2**2))*LC(7)* &
!                  SP(1))/(E**y*shat**2*(-(MV2*shat) + smT2)**2) + &
!               (2*(-2*shat**2*smT + E**(3*y)*shat*smT2 + &
!                    2*E**(2*y)*smT*(-(MV2*shat) + smT2) - &
!                    E**y*shat*(-6*MV2*shat + shat**2 + 4*smT2))*LC(6)* &
!                  SP(2))/(E**y*shat**2*(MV2*shat - smT2)) + &
!               (2*(shat - E**y*smT)*LC(9)* &
!                  (-2*E**(4*y)*smT2* &
!                     (2*MV2*(-2*MT2 + MV2)*shat**2 + &
!                       (4*MT2 - MV2)*shat*smT2 + smT2**2)*SP(3) + &
!                    4*E**y*MV2**2*shat**3*smT*SP(4) + &
!                    2*shat**2*smT2*(-2*MV2*shat + smT2)*SP(4) + &
!                    E**(3*y)*shat*smT* &
!                     (2*MV2*shat*(MV2*shat + 3*smT2)*SP(3) + &
!                       smT2*(-(MV2*shat) + smT2)*SP(4)) + &
!                    E**(2*y)*shat**2* &
!                     (-4*MV2*shat*smT2*SP(3) + &
!                       (2*MV2**2*shat*(2*MT2 + shat) - &
!                          MV2*(4*MT2 + 2*MV2 + 3*shat)*smT2 + smT2**2)* &
!                        SP(4))))/(E**(2*y)*shat**2*(MV2*shat - smT2)**3) &
!                + (LC(13)*(-2*E**(7*y)*shat*smT**5*(MV2*shat - smT2)* &
!                     SP(3) + 2*E**(6*y)*shat**2*(MV2*shat - smT2)*smT2* &
!                     (MV2**2 + 3*smT2)*SP(3) + 8*shat**4*smT2**2*SP(4) - &
!                    4*E**y*shat**3*smT**3*(3*MV2*shat + 5*smT2)*SP(4) - &
!                    4*E**(2*y)*shat**2*smT2* &
!                     (2*MV2*shat**2*(2*MT2 - MV2 + shat) - &
!                       shat*(4*MT2 + 5*MV2 + 2*shat)*smT2 - 5*smT2**2)* &
!                     SP(4) - 2*E**(5*y)*shat*smT*(MV2*shat - smT2)* &
!                     (shat*(2*MV2**2*shat + MV2*smT2 + 3*shat*smT2)* &
!                        SP(3) + &
!                       2*(2*MV2**2*shat**2 - 5*MV2*shat*smT2 + smT2**2)* &
!                        SP(4)) - &
!                    2*E**(3*y)*shat*smT* &
!                     (MV2*shat**3*(MV2*shat - smT2)*SP(3) - &
!                       2*(MV2**2*shat**3*(-6*MV2 + shat) + &
!                          2*MV2*shat**2*(4*MT2 + 7*MV2 + 2*shat)*smT2 - &
!                          shat*(8*MT2 + 19*MV2 + 5*shat)*smT2**2 + &
!                          3*smT2**3)*SP(4)) + &
!                    2*E**(4*y)* &
!                     (shat**3*(MV2*shat - smT2)* &
!                        (MV2**2*shat + 2*MV2*smT2 + shat*smT2)*SP(3) + &
!                       2*(2*MV2**3*shat**4*(4*MT2 + 2*MV2 + shat) - &
!                          2*MV2**2*shat**3*(12*MT2 + 5*MV2 + 4*shat)* &
!                           smT2 + &
!                          4*MV2*shat**2*(5*MT2 + 2*MV2 + shat)*smT2**2 + &
!                          shat*(-4*MT2 + MV2 + 2*shat)*smT2**3 - smT2**4) &
!                         *SP(4))))/ &
!                (E**(4*y)*shat**2*(-shat + E**y*smT)* &
!                  (MV2*shat - smT2)**3) + &
!               (4*LC(3)*(E**(7*y)*shat*(MV2*shat - smT2)*smT2**3*SP(3) + &
!                    E**(6*y)*smT**3* &
!                     (-2*MV2**2*shat**3*(-4*MT2 + 2*MV2 + shat) + &
!                       MV2*shat**2*(-12*MT2 + 6*MV2 + shat)*smT2 + &
!                       shat*(4*MT2 - 6*MV2 + shat)*smT2**2 + 2*smT2**3)* &
!                     SP(3) + 2*MV2*shat**5*smT**3*SP(4) - &
!                    4*E**y*MV2*shat**4*smT2*(MV2*shat + smT2)*SP(4) + &
!                    E**(5*y)*shat*smT2* &
!                     ((2*MV2**2*shat**3*(-8*MT2 + 5*MV2 + 2*shat) + &
!                          MV2*shat**2*(24*MT2 - 5*(2*MV2 + shat))*smT2 + &
!                          shat*(-8*MT2 + 12*MV2 + shat)*smT2**2 - &
!                          4*smT2**3)*SP(3) + &
!                       smT2*(-2*MV2*shat + smT2)*(-(MV2*shat) + smT2)* &
!                        SP(4)) - &
!                    E**(2*y)*shat**2*smT* &
!                     (2*MV2**2*shat**4*SP(3) + &
!                       (8*MT2*MV2**2*shat**3 + &
!                          MV2*shat**2*(-12*MT2 - 14*MV2 + shat)*smT2 + &
!                          (4*(MT2 + MV2) - shat)*shat*smT2**2 - 2*smT2**3 &
!                          )*SP(4)) + &
!                    E**(3*y)*shat* &
!                     (2*MV2**2*shat**4*(MV2*shat + 3*smT2)*SP(3) + &
!                       (2*MV2**3*(8*MT2 - shat)*shat**4 + &
!                          8*MV2**2*shat**3*(-4*MT2 + shat)*smT2 + &
!                          MV2*(24*MT2 - 16*MV2 - 7*shat)*shat**2* &
!                           smT2**2 + &
!                          shat*(-8*MT2 + 12*MV2 + shat)*smT2**3 - &
!                          4*smT2**4)*SP(4)) + &
!                    E**(4*y)*smT* &
!                     (-(shat**2* &
!                          (2*MV2**2*shat**3*(-4*MT2 + 4*MV2 + shat) - &
!                            3*MV2*shat**2*(-4*MT2 + shat)*smT2 + &
!                            shat*(-4*MT2 + 6*MV2 + shat)*smT2**2 - &
!                            2*smT2**3)*SP(3)) + &
!                       (2*MV2**3*shat**4*(-4*MT2 + shat) + &
!                          2*MV2**2*(8*MT2 - 5*shat)*shat**3*smT2 + &
!                          MV2*shat**2*(-12*MT2 + 6*MV2 + 11*shat)* &
!                           smT2**2 + &
!                          shat*(4*MT2 - 3*(2*MV2 + shat))*smT2**3 + &
!                          2*smT2**4)*SP(4))))/ &
!                (E**(3*y)*shat**3*(-shat + E**y*smT)* &
!                  (MV2*shat - smT2)**3) + &
!               (2*(shat*smT*(2*MV2*shat - smT2) - E**(3*y)*smT2**2 + &
!                    E**(2*y)*shat*smT*(MV2*(MV2 + shat) + smT2) - &
!                    E**y*(MV2*shat**2*(3*MV2 + shat) - 2*MV2*shat*smT2 + &
!                       smT2**2))*LC(5)*SP(5))/ &
!                (E**y*shat*(-(MV2*shat) + smT2)**2) + &
!               (2*(shat - E**y*smT)* &
!                  (E**y*MV2*shat*smT*(5*MV2*shat - smT2) + &
!                    2*shat*smT2*(-2*MV2*shat + smT2) + &
!                    E**(3*y)*smT*(-(MV2*shat) + smT2)* &
!                     (3*MV2*shat + smT2) + &
!                    E**(2*y)*shat* &
!                     (MV2**2*shat*(4*MT2 - MV2 + shat) - &
!                       MV2*(4*MT2 + MV2 - 2*shat)*smT2 - 3*smT2**2))* &
!                  LC(4)*SP(6))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (4*smT*(shat - E**y*smT)* &
!                  (2*shat*smT2 - E**y*smT*(3*MV2*shat + smT2) + &
!                    E**(2*y)*(4*MT2*(-(MV2*shat) + smT2) + &
!                       MV2*(MV2*shat + smT2)))*LC(10)*SP(6))/ &
!                (E**y*(-(MV2*shat) + smT2)**3) + &
!               (4*(-2*shat**3*smT**3 + &
!                    2*E**(4*y)*MV2*shat**2*smT*(MV2*shat - smT2) + &
!                    3*E**y*shat**2*smT2*(MV2*shat + smT2) - &
!                    2*E**(2*y)*shat*smT* &
!                     (MV2*(-2*MT2 + MV2 - shat)*shat**2 + &
!                       shat*(2*MT2 + MV2 + shat)*smT2 + smT2**2) + &
!                    E**(3*y)*(2*MV2**2*(MV2 - shat)*shat**3 - &
!                       4*MV2*(MT2 + MV2)*shat**2*smT2 + &
!                       shat*(4*MT2 + 5*MV2 + 2*shat)*smT2**2 - smT2**3))* &
!                  LC(14)*SP(6))/(E**(3*y)*shat*(MV2*shat - smT2)**3) + &
!               (2*(E**(4*y)*MV2*shat*(MV2*shat - smT2)*smT2 + &
!                    2*shat**2*smT2**2 - &
!                    3*E**y*shat*smT*(MV2**2*shat**2 + smT2**2) + &
!                    E**(3*y)*shat*smT* &
!                     (-(MV2**2*shat*(-4*MT2 + 2*MV2 + shat)) - &
!                       4*MT2*MV2*smT2 + smT2**2) - &
!                    E**(2*y)*(4*(MT2 - MV2)*MV2**2*shat**3 + &
!                       MV2*(-4*MT2 + 3*MV2 - shat)*shat**2*smT2 + &
!                       shat*(-6*MV2 + shat)*smT2**2 + smT2**3))*LC(2)* &
!                  SP(7))/(E**(2*y)*shat*(MV2*shat - smT2)**3) + &
!               (4*smT*(-shat + E**y*smT)* &
!                  (2*MV2*shat**2 + E**y*(-5*MV2*shat*smT + smT**3) + &
!                    E**(2*y)*(4*MT2*(-(MV2*shat) + smT2) + &
!                       MV2*(MV2*shat + smT2)))*LC(8)*SP(7))/ &
!                (E**y*(-(MV2*shat) + smT2)**3) + &
!               (4*(2*shat**3*smT**3 - &
!                    3*E**y*shat**2*smT2*(MV2*shat + smT2) + &
!                    2*E**(2*y)*shat*smT* &
!                     (MV2*(-2*MT2 + MV2)*shat**2 + &
!                       (2*MT2 + MV2)*shat*smT2 + smT2**2) + &
!                    E**(3*y)*(-2*MV2**3*shat**3 + &
!                       4*MV2*(MT2 + MV2)*shat**2*smT2 - &
!                       (4*MT2 + 5*MV2)*shat*smT2**2 + smT2**3))*LC(12)* &
!                  SP(7))/(E**(3*y)*shat*(MV2*shat - smT2)**3)) - &
!            (4*kapPr*(shat - E**y*smT)* &
!               (shat*smT*(-((MV2*shat - smT2)*(SP(3) + SP(4))*SP(5)) + &
!                    ((-(MV2*shat) + smT2)*SP(1) + &
!                       2*shat*(SP(3) + SP(4))*SP(6))*SP(7)) + &
!                 E**y*((MV2*shat - smT2)*smT2*(SP(3) + SP(4))*SP(5) + &
!                    shat*(shat*(MV2*shat - smT2)*SP(1) - &
!                       2*(smT2*SP(3) + MV2*shat*SP(4))*SP(6))*SP(7))))/ &
!             (E**y*shat*(-(MV2*shat) + smT2)**2)) + &
!         SI(10)*(kap*((E**y*(4*shat*smT*(MV2*shat - smT2)**3 + &
!                    2*shat*smT**3*smT2*(-(MV2*shat) + smT2) - &
!                    E**(3*y)*shat**2*smT2**2*(5*MV2*shat + 7*smT2) - &
!                    2*E**(2*y)*shat*smT* &
!                     ((-(MV2*shat*smT) + smT**3)**2 - &
!                       smT2*(3*MV2**2*shat**2 + 8*MV2*shat*smT2 + &
!                          smT2**2)) + &
!                    E**y*(2*(-(MV2*shat*smT) + smT**3)**2* &
!                        (-3*MV2*shat + shat**2 + 4*smT2) - &
!                       shat*smT2* &
!                        (MV2**2*shat**2*(MV2 + 2*shat) + &
!                          3*MV2*(3*MV2 - 2*shat)*shat*smT2 + &
!                          2*(MV2 + 2*shat)*smT2**2)))*LC(1))/ &
!                (shat**2*(MV2*shat - smT2)**3) + &
!               ((E**(3*y)*shat**2*(MV2*shat - 3*smT2)*smT2* &
!                     (2*MV2*shat + smT2) - &
!                    shat*smT*(MV2*shat - smT2)* &
!                     (2*MV2**2*shat**2 - 4*MV2*shat*smT2 + smT2**2) + &
!                    E**(2*y)*shat*smT* &
!                     (-4*MV2**3*shat**3 + 9*MV2**2*shat**2*smT2 + &
!                       5*MV2*shat*smT2**2 + 2*smT2**3) + &
!                    E**y*(8*MV2**4*shat**4 - 10*MV2**3*shat**3*smT2 - &
!                       MV2*shat**2*(16*MV2 + shat)*smT2**2 + &
!                       shat*(20*MV2 + shat)*smT2**3 - 8*smT2**4))*LC(11)) &
!                 /(E**y*shat**2*(MV2*shat - smT2)**3) + &
!               ((2*E**(3*y)*shat**2*smT**3*(MV2*shat + 4*smT2) - &
!                    E**(2*y)*shat*smT2*(3*MV2*shat + smT2)* &
!                     (MV2*shat + 4*smT2) + &
!                    shat*(MV2*shat - smT2)* &
!                     (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2) + &
!                    E**y*smT*(-(MV2**2*shat**3*(3*MV2 + 2*shat)) + &
!                       3*MV2*shat**2*(5*MV2 + shat)*smT2 - &
!                       shat*(4*MV2 + shat)*smT2**2 + 2*smT2**3))*LC(15))/ &
!                (shat*(MV2*shat - smT2)**3) + &
!               (2*(E**(3*y)*shat**2*smT**3 + &
!                    shat*(MV2*shat - smT2)*(2*MV2*shat - smT2) - &
!                    E**(2*y)*shat*smT2*(MV2*shat + smT2) - &
!                    E**y*smT*(MV2*shat**2*(MV2 + 2*shat) - &
!                       2*shat*(2*MV2 + shat)*smT2 + 2*smT2**2))*LC(7)* &
!                  SP(1))/(shat**2*(-(MV2*shat) + smT2)**2) + &
!               (2*(-2*MV2*shat**2 + &
!                    E**y*smT*(2*MV2*shat + shat**2 - 2*smT2) + shat*smT2) &
!                   *LC(6)*SP(2))/(shat**2*(MV2*shat - smT2)) + &
!               (LC(13)*(2*E**(3*y)*shat*smT*smT2*(-(MV2*shat) + smT2)**2* &
!                     SP(3) - 4*E**y*MV2*shat**2*smT**3* &
!                     (MV2*shat + 3*smT2)*SP(4) + &
!                    4*smT2*(-2*MV2**3*shat**3 + 6*MV2**2*shat**2*smT2 - &
!                       3*MV2*shat*smT2**2 + smT2**3)*SP(4) - &
!                    2*E**(2*y)*shat**2*smT2* &
!                     ((-(MV2*shat) + smT2)**2*SP(3) - &
!                       4*MV2*shat*smT2*SP(4))))/ &
!                (E**y*shat**2*smT*(MV2*shat - smT2)**3) - &
!               (4*LC(3)*(E**(4*y)*shat**2*smT**3*(3*MV2*shat - smT2)* &
!                     smT2*SP(3) - &
!                    E**(3*y)*shat*smT2* &
!                     ((-(MV2*shat*smT) + smT**3)**2 + &
!                       MV2*shat*(5*MV2*shat - smT2)*smT2)*SP(3) + &
!                    2*smT**3*(-2*MV2*shat + smT2)* &
!                     (MV2**2*shat**2 - MV2*shat*smT2 + smT2**2)*SP(4) + &
!                    E**(2*y)*smT* &
!                     (smT2*(-3*MV2**2*shat**4 + &
!                          2*MV2*shat**2*(3*MV2 + 2*shat)*smT2 - &
!                          shat*(6*MV2 + shat)*smT2**2 + 2*smT2**3)*SP(3) &
!                        + MV2*shat**3*smT2*(-3*MV2*shat + smT2)*SP(4)) + &
!                    E**y*shat*smT2* &
!                     (-2*MV2**2*shat**2*smT2*(SP(3) - 2*SP(4)) - &
!                       3*MV2*shat*smT2**2*SP(4) + smT2**3*SP(4) + &
!                       2*MV2**3*shat**3*(SP(3) + SP(4)))))/ &
!                (shat**3*smT*(MV2*shat - smT2)**3) + &
!               (2*LC(9)*(4*E**(5*y)*shat**2*smT2**3*SP(3) - &
!                    2*E**(4*y)*shat*smT**3*smT2*(3*MV2*shat + smT2)* &
!                     SP(3) + E**y*MV2*shat**2*smT2* &
!                     (MV2*(MV2 - shat)*shat + (MV2 + shat)*smT2)*SP(4) + &
!                    shat*smT*(MV2*shat - smT2)* &
!                     (2*MV2**2*shat**2 - 2*MV2*shat*smT2 + smT2**2)*SP(4) &
!                      - 2*E**(3*y)*smT2**2* &
!                     ((2*MV2*shat**3 - shat*(3*MV2 + 2*shat)*smT2 + &
!                          smT2**2)*SP(3) - MV2*shat**3*SP(4)) + &
!                    E**(2*y)*shat*smT*smT2* &
!                     (2*MV2**2*shat**2*(SP(3) - 2*SP(4)) + &
!                       MV2*shat*smT2*SP(4) - smT2**2*(2*SP(3) + SP(4))))) &
!                 /(shat**2*smT*(MV2*shat - smT2)**3) + &
!               (2*E**y*smT*(MV2*shat*(-MV2 + shat) + &
!                    E**y*smT*(MV2*shat + smT2 - 2*shat*smT*Cosh(y)))* &
!                  LC(5)*SP(5))/(shat*(-(MV2*shat) + smT2)**2) + &
!               (2*smT*(MV2*shat*smT*(MV2*shat - smT2) + &
!                    2*E**(3*y)*shat*(3*MV2*shat - 2*smT2)*smT2 - &
!                    E**(2*y)*smT* &
!                     (5*MV2**2*shat**2 + MV2*shat*smT2 - 2*smT2**2) + &
!                    E**y*shat* &
!                     (MV2**2*(MV2 - 3*shat)*shat + &
!                       MV2*(MV2 + 5*shat)*smT2 - 2*smT2**2))*LC(4)*SP(6)) &
!                 /(shat*(MV2*shat - smT2)**3) + &
!               (4*E**(2*y)*(-MV2 + E**y*smT)* &
!                  (MV2*shat + smT*(-2*E**y*shat + smT))*smT2*LC(10)*SP(6) &
!                  )/(-(MV2*shat) + smT2)**3 + &
!               (4*(-2*MV2**3*shat**3 - &
!                    2*MV2**2*shat**2*(E**y*shat - 3*smT)*smT + &
!                    smT**5*(-(E**y*shat) + smT) + &
!                    MV2*shat*(2*E**y*shat - 3*smT)*(E**y*shat + smT)*smT2 &
!                    )*LC(14)*SP(6))/(shat*(MV2*shat - smT2)**3) - &
!               (2*smT*(E**(2*y)*(-5*MV2*shat*smT**3 + smT**5) + &
!                    3*MV2*shat*smT*(-(MV2*shat) + smT2) + &
!                    E**(3*y)*shat*smT2*(MV2*shat + smT2) + &
!                    E**y*shat* &
!                     (-(MV2*shat*smT2) - smT2**2 + &
!                       2*MV2**2*(shat**2 + smT2)))*LC(2)*SP(7))/ &
!                (shat*(MV2*shat - smT2)**3) + &
!               (4*E**y*smT2*(2*MV2*shat*smT - 2*smT**3 + &
!                    2*E**(3*y)*shat*smT2 - &
!                    E**(2*y)*smT*(3*MV2*shat + smT2) + &
!                    E**y*(MV2*(MV2 - 2*shat)*shat + (MV2 + 2*shat)*smT2)) &
!                   *LC(8)*SP(7))/(-(MV2*shat) + smT2)**3 + &
!               (4*(2*MV2**3*shat**3 + &
!                    3*MV2*shat*smT**3*(E**y*shat + smT) - &
!                    6*MV2**2*shat**2*smT2 - &
!                    smT2**2*(E**y*shat*(2*E**y*shat - smT) + smT2))* &
!                  LC(12)*SP(7))/(shat*(MV2*shat - smT2)**3)) + &
!            (4*kapPr*(shat*smT**3*(MV2*shat - smT2)*SP(1)*SP(7) + &
!                 E**(2*y)*shat*smT*smT2*(SP(3) + SP(4))* &
!                  (MV2*shat*SP(5) - smT2*SP(5) - 2*shat*SP(6)*SP(7)) + &
!                 E**y*smT2*(smT2*(-(MV2*shat) + smT2)*(SP(3) + SP(4))* &
!                     SP(5) + shat* &
!                     (shat*(-(MV2*shat) + smT2)*SP(1) + &
!                       2*(smT2*SP(3) + MV2*shat*SP(4))*SP(6))*SP(7))))/ &
!             (shat*smT*(-(MV2*shat) + smT2)**2))
! 
!      
! 
!      calc_MassiveBox_QP= dcmplx(calc_MassiveBox)
!      
! return              
! end function               

