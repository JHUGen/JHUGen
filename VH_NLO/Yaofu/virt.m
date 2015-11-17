
TIReduction = Get["/afs/cern.ch/user/m/maschulz/lib/MathLib/TIReduction.dat"];


LSP[VecID1_,VecID2_]:=Dot[VecID1,VecID2];
Unprotect[Dot];
Dot[VecID1_,VecID2_] := Dot[VecID2,VecID1]   /; Not[ OrderedQ[ {VecID1,VecID2} ] ];
Dot[-VecID1_, VecID2_]:=-Dot[VecID1,VecID2];
Dot[ VecID1_,-VecID2_]:=-Dot[VecID1,VecID2];
Dot[VecID1_-VecID2_,VecID3_] := Dot[VecID1,VecID3]-Dot[VecID2,VecID3];
Dot[VecID1_,VecID2_-VecID3_] := Dot[VecID1,VecID2]-Dot[VecID1,VecID3];
Dot[VecID1_+VecID2_,VecID3_] := Dot[VecID1,VecID3]+Dot[VecID2,VecID3];
Dot[VecID1_,VecID2_+VecID3_] := Dot[VecID1,VecID2]+Dot[VecID1,VecID3];




insertMasses = { p1.p1 -> m1^2,
                 p2.p2 -> m2^2,
                 p3.p3 -> m3^2,
                 p0.p0 -> m0^2
};

insertMassless = {m0 -> 0,
                  m1 -> 0,
                  m2 -> 0,
                  m3 -> 0
};

insertDimD = {DimST -> 4-2\[CurlyEpsilon]};


(* B1 *)
TCInput=TC[2,1,p1,m0,m1];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 



(* B00 *)
TCInput=TC[2,0,0,p1,m0,m1];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 



(* B11 *)
TCInput=TC[2,1,1,p1,m0,m1];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 


(* C1 *)
TCInput=TC[3,1,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify]


(* C2 *)
TCInput=TC[3,2,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify]


(* C00 *)
TCInput=TC[3,0,0,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 




(* C11 *)
TCInput=TC[3,1,1,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 




(* C22 *)
TCInput=TC[3,2,2,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
(*TheCoeff = TheCoeff    //.insirtDimD;*)
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 



(* ::InheritFromParent:: *)
(*((-1+DimST) (p1.p1)^2 (p1.p2-p2.p2) SI[2,p1,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)-(((-2+DimST) (p1.p2)^3+(p1.p1)^2 p2.p2+p1.p1 p1.p2 ((-2+DimST) p1.p2+(3-2 DimST) p2.p2)) SI[2,p2,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+(((-2+DimST) (p1.p2)^3+p1.p1 p1.p2 (2 p1.p2+(3-2 DimST) p2.p2)+(p1.p1)^2 (-(-1+DimST) p1.p2+2 (-2+DimST) p2.p2)) SI[2,-p1+p2,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+((p1.p1)^2 ((-2+DimST) (p1.p2)^2-2 (-1+DimST) p1.p2 p2.p2+p2.p2 (p1.p1+(-1+DimST) p2.p2)) SI[3,p1,p2,0,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)*)


(* C12 *)
TCInput=TC[3,1,2,p1,p2,m0,m1,m2];
TheCoeff = TCInput //.TIReduction;
(*TheCoeff = TheCoeff    //.insertMasses;*)
TheCoeff = TheCoeff    //.insertMassless;
TheCoeff = TheCoeff    // Expand;
TheCoeff = Collect[ TheCoeff, SI[___], Simplify] 




B1[m0_,m1_,p1_]:=SI[1,m0]/(2 p1.p1)-SI[1,m1]/(2 p1.p1)-1/2 SI[2,p1,m0,m1];
B11[m0_,m1_,p1_,DimST_]:=SI[1,m0]/(4 (-1+DimST))+SI[1,m1]/(4 (-1+DimST))+(p1.p1 SI[2,p1,m0,m1])/(4-4 DimST);
C1[m0_,m1_,m2_,p1_,p2_]:=(p1.p2 SI[2,p1,m0,m1])/(2 (p1.p2)^2-2 p1.p1 p2.p2)+(p2.p2 SI[2,p2,m0,m2])/(-2 (p1.p2)^2+2 p1.p1 p2.p2)+((-p1.p2+p2.p2) SI[2,-p1+p2,m1,m2])/(2 ((p1.p2)^2-p1.p1 p2.p2))+((p1.p1-p1.p2) p2.p2 SI[3,p1,p2,m0,m1,m2])/(2 ((p1.p2)^2-p1.p1 p2.p2));
C2[m0_,m1_,m2_,p1_,p2_]:=(p1.p1 SI[2,p1,m0,m1])/(-2 (p1.p2)^2+2 p1.p1 p2.p2)+(p1.p2 SI[2,p2,m0,m2])/(2 (p1.p2)^2-2 p1.p1 p2.p2)+((p1.p1-p1.p2) SI[2,-p1+p2,m1,m2])/(2 (p1.p2)^2-2 p1.p1 p2.p2)+(p1.p1 (p1.p2-p2.p2) SI[3,p1,p2,m0,m1,m2])/(-2 (p1.p2)^2+2 p1.p1 p2.p2);
C00[m0_,m1_,m2_,p1_,p2_,DimST_]:=(p1.p1 (p1.p2-p2.p2) SI[2,p1,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2))+((-p1.p1+p1.p2) p2.p2 SI[2,p2,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2))+(p1.p2 (-p1.p1+2 p1.p2-p2.p2) SI[2,-p1+p2,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2))+(p1.p1 p2.p2 (p1.p1-2 p1.p2+p2.p2) SI[3,p1,p2,0,0,0])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2))
C11[m0_,m1_,m2_,p1_,p2_,DimST_]:=(p1.p2 SI[1,m0])/(4 p1.p1 (p1.p2)^2-4 (p1.p1)^2 p2.p2)+((-2 (p1.p2)^2+p1.p1 p2.p2+p1.p2 p2.p2) SI[1,m1])/(4 p1.p1 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))+((p1.p2-p2.p2) SI[1,m2])/(4 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))-(((-2+DimST) (p1.p2)^3+(3-2 DimST) p1.p1 p1.p2 p2.p2+(-2+DimST) (p1.p2)^2 p2.p2+p1.p1 (p2.p2)^2) SI[2,p1,m0,m1])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+((-1+DimST) (-p1.p1+p1.p2) (p2.p2)^2 SI[2,p2,m0,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+(((-2+DimST) (p1.p2)^3+2 (p1.p2)^2 p2.p2+2 (-2+DimST) p1.p1 (p2.p2)^2+p1.p2 p2.p2 ((3-2 DimST) p1.p1-(-1+DimST) p2.p2)) SI[2,-p1+p2,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+((p2.p2)^2 ((-1+DimST) (p1.p1)^2+(-2+DimST) (p1.p2)^2+p1.p1 (-2 (-1+DimST) p1.p2+p2.p2)) SI[3,p1,p2,m0,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2);
C22[m0_,m1_,m2_,p1_,p2_,DimST_]:=(p1.p2 SI[1,m0])/(4 (p1.p2)^2 p2.p2-4 p1.p1 (p2.p2)^2)+((-p1.p1+p1.p2) SI[1,m1])/(4 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))+((-2 (p1.p2)^2+p1.p1 (p1.p2+p2.p2)) SI[1,m2])/(4 p2.p2 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))+((-1+DimST) (p1.p1)^2 (p1.p2-p2.p2) SI[2,p1,m0,m1])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)-(((-2+DimST) (p1.p2)^3+(p1.p1)^2 p2.p2+p1.p1 p1.p2 ((-2+DimST) p1.p2+(3-2 DimST) p2.p2)) SI[2,p2,m0,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+(((-2+DimST) (p1.p2)^3+p1.p1 p1.p2 (2 p1.p2+(3-2 DimST) p2.p2)+(p1.p1)^2 (-(-1+DimST) p1.p2+2 (-2+DimST) p2.p2)) SI[2,-p1+p2,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+((p1.p1)^2 ((-2+DimST) (p1.p2)^2-2 (-1+DimST) p1.p2 p2.p2+p2.p2 (p1.p1+(-1+DimST) p2.p2)) SI[3,p1,p2,m0,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2);
C12[m0_,m1_,m2_,p1_,p2_,DimST_]:=-(SI[1,m0]/(4 (p1.p2)^2-4 p1.p1 p2.p2))+((p1.p2-p2.p2) SI[1,m1])/(4 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))+((-p1.p1+p1.p2) SI[1,m2])/(4 (p1.p1-2 p1.p2+p2.p2) (-(p1.p2)^2+p1.p1 p2.p2))-(p1.p1 ((p1.p2)^2+(-2+DimST) p1.p1 p2.p2-(-1+DimST) p1.p2 p2.p2) SI[2,p1,m0,m1])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)-(p2.p2 ((p1.p2)^2+p1.p1 (-(-1+DimST) p1.p2+(-2+DimST) p2.p2)) SI[2,p2,m0,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+(((-2+DimST) (p1.p1)^2 p2.p2+(p1.p2)^2 ((-4+DimST) p1.p2+p2.p2)+p1.p1 ((p1.p2)^2-3 (-2+DimST) p1.p2 p2.p2+(-2+DimST) (p2.p2)^2)) SI[2,-p1+p2,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)+(p1.p1 p2.p2 (p1.p1 (-(-1+DimST) p1.p2+(-2+DimST) p2.p2)+p1.p2 (DimST p1.p2-(-1+DimST) p2.p2)) SI[3,p1,p2,m0,m1,m2])/(4 (-2+DimST) ((p1.p2)^2-p1.p1 p2.p2)^2)


MV[p1_,p2_,DimST_]:=J.eps*(4p1.p2(SI[3,p1,p2,0,0,0]+2C12[0,0,0,p1,p2,DimST]+C1[0,0,0,p1,p2]-C2[0,0,0,p1,p2])+(3DimST-4)C00[0,0,0,p1,p2,DimST])


MVNLO = MV[p1,p2,DimST] //.insertDimD;


MVNLO = MVNLO //.insertMasses;


MVNLO=Limit[MVNLO, m1->0];


MVNLO=Limit[MVNLO, m2->0] //Simplify


MVNLOprefact = (4 (-1+\[CurlyEpsilon]) SI[2,p1,0,0]-4 (-1+\[CurlyEpsilon]) SI[2,p2,0,0]-4 SI[2,-p1+p2,0,0]+7 \[CurlyEpsilon] SI[2,-p1+p2,0,0]-8 p1.p2 SI[3,p1,p2,0,0,0]+8 \[CurlyEpsilon] p1.p2 SI[3,p1,p2,0,0,0])/(2 (-1+\[CurlyEpsilon]))


MVNLOprefact = Series[MVNLOprefact,{\[CurlyEpsilon],Infinity,2}]
