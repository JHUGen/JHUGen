(* ::Package:: *)

InitPath = "/home/schulze/lib/MathLib/";
FeynArtsPath = "/home/schulze/lib/FeynArts-3.5/";
FeynArtsToFormPath = "/home/schulze/lib/FeynArtsToForm/";
ProjectPath = "/home/schulze/projects/JHUGenerator/JHUGen-origin/VH_NLO/";

Get[ InitPath<>"StandardInit.m" ];
Get[ InitPath<>"SIOrder.m" ];
Get[ ProjectPath<>"Tree_ggVH_output.dat" ];



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
         
         TheCoeff = TCList[[i]] //.TIReduction //Expand;
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




TheRedAmpList // ByteCount  ;
% /1024^2 // N 





TotalAmp = TheRedAmpList[[3]]+TheRedAmpList[[4]]+TheRedAmpList[[5]]+TheRedAmpList[[6]]+TheRedAmpList[[7]]+TheRedAmpList[[8]]+TheRedAmpList[[9]]+TheRedAmpList[[10]]+TheRedAmpList[[11]]+TheRedAmpList[[12]];
TotalAmp // ByteCount  ;
% /1024^2 // N





TotalAmp = TotalAmp  // SimplAmp;
TotalAmp // ByteCount  ;
% /1024^2 // N



TotalAmp  // FortranForm >> ProjectPath<>"test.f"



