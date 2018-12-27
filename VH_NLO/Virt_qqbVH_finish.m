InitPath = $HomeDirectory<>"/lib/MathLib/";
FeynArtsPath = $HomeDirectory<>"/lib/FeynArts-3.5/";
FeynArtsToFormPath = $HomeDirectory<>"/lib/FeynArtsToForm/";
ProjectPath = $HomeDirectory<>"/1git_JHUGen/VH_NLO/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[ ProjectPath<>"Virt_qqbVH_output.dat" ];

Clear[TCReplList,SIList,TheRedAmpList];

TIReduction = Get[ InitPath<>"TIReduction.dat" ];


insertMasses = { p1.p1 -> 0,
                 p2.p2 -> 0,
                 p3.p3 -> MZ^2,
                 p4.p4 -> MH^2,
                 p1.p2 -> shat/2,
                 p1.p3 ->(that-MZ^2)/(-2),
                 p1.p4 ->(uhat-MH^2)/(-2),
                 p2.p3 ->(uhat-MZ^2)/(-2) ,
                 p2.p4 ->(that-MH^2)/(-2) ,
                 p3.p4 ->(shat-MZ^2-MH^2)/2
};

TheAmpList    = StringReplace[AmpList,"["->"M["] // ToExpression;
TheRedAmpList = StringReplace[AmpList,"["->"redM["] // ToExpression;


   (* generate list of all tensor integral coefficients *)
   TCList={};
   For[i=1,i<=Length[TheAmpList],i++,
         TCList = Union[TCList, Cases[ TheAmpList[[i]] ,TC[___],100] ];
   ];

   (* simplify this list *)
   TCReplList={};
   For[i=1, i<=Length[TCList],i++,
         Print[i,"/",Length[TCList],": simplifying ",TCList[[i]]];
         
         TheCoeff = TCList[[i]] //.TIReduction;
         TheCoeff = TheCoeff    //.insertMasses;
         TheCoeff = TheCoeff    // Expand;
         TheCoeff = SimplifyTC[ TheCoeff ];
         TheCoeff = Collect[ TheCoeff, SI[___], Simplify];
         TCReplList = Append[TCReplList,TCList[[i]]-> TheCoeff ];
   ];
 

  SIList={};
  For[i=1,i<=Length[TheAmpList],i++,
         Print["reducing diagram ",i,"/",Length[TheAmpList]];
         dummy = TheAmpList[[i]] //. TCReplList; 
         Print["expanding in D-4"];
         dummy = Expand[ dummy ] //. InsertFinitePart   /. DSTm4->0;
         Print["converting to invariants"];
         dummy = SIInvariants[ dummy,insertMasses];         
         dummy = SIOrder[ dummy ]  //. insertMasses;
         dummy = dummy //. { SI[1,0]->0 };
         Print["simplifying amplitude"];
         dummy = dummy //. PropDenom[x_]->1/x    // SimplAmp;
         Evaluate[ TheRedAmpList[[i]] ] = dummy//Simplify;

         (* generate list of integrals in terms of invariants *)
         SIList = Union[SIList, Cases[ TheRedAmpList[[i]] ,SI[___],100] ] //SIOrder //Simplify;
   ];

SIList // TableForm



AnalyticSI={
     SI[2,0,0,0] -> DeltaUV - DeltaIR,
     SI[2,shat_,0,0] -> DeltaUV + 2 + Log[MURen^2/shat] + I*Pi,
     SI[3,0,shat_,0,0,0,0] -> 1/shat * ( DeltaIR2 + DeltaIR*Log[-MURen^2/shat] + 1/2*Log[-MURen^2/shat]^2 -Pi^2/6)
};




AnalyticRes = 1/(16 Pi^4) TheRedAmpList[[3]] //.{GS^2->4*Pi*alphaS,EL^2->4*Pi*alpha, ICol[1]->1} // FullSimplify[#,shat+that+uhat==MZ^2+MH^2]&

Tree = EL^2*CW^(-1)*SW^(-1)*MZ*cI*( IZ[3,-1]*SME[1,-1] + IZ[3,1]*SME[1,1])/(shat-MZ^2)  //.{GS^2->4*Pi*alphaS,EL^2->4*Pi*alpha}



AnalyticResFF = AnalyticRes / Tree  // FullSimplify;

AnalyticResFF = AnalyticResFF //.  AnalyticSI // FullSimplify



WFRC = -(GS^2/4/Pi)*CF/4/Pi*(DeltaUV - DeltaIR)//.{GS^2->4*Pi*alphaS,EL^2->4*Pi*alpha, ICol[1]->1, CF-> (NCol^2-1)/2/NCol}


AnalyticResFFUVren = AnalyticResFF + WFRC // FullSimplify
