InitPath =  $HomeDirectory<>"/lib/MathLib/";
FeynArtsPath =  $HomeDirectory<>"/lib/FeynArts-3.9/";
FeynArtsToFormPath =  $HomeDirectory<>"/lib/FeynArtsToForm/";
ProjectPath =  $HomeDirectory<>"/1git_JHUGen/VH_NLO/";
Unprotect[Diagonal];
Clear[Diagonal];

Get[ FeynArtsPath<>"FeynArts.m" ];
Get[ FeynArtsToFormPath<>"ExtractAmp.m" ];


TheProcess = {F[3, {1}],-F[3, {1}]} -> {V[2],S[1]};

TheLevel = {Classes};
TheGenericModel = "Lorentz";
TheModel = SMQCD;
PaintIt  = True;
WriteAmp = True;

SetOptions[InsertFields,
    	InsertionLevel -> TheLevel,
    	GenericModel -> TheGenericModel ,
    	Model -> TheModel,
        ExcludeParticles -> { SV[2],SV[3] },
        LastSelections-> {},
        Restrictions -> { NoLightFHCoupling }
];

SetOptions[Paint, Numbering->Simple, ColumnsXRows->{5,7} ];
CKM = IndexDelta;



(* generate tree level diagrams *)
TheTreeTop  = CreateTopologies[ 0, 2 -> 2 ];
TheTreeDiag = InsertFields[ TheTreeTop, TheProcess ];

If[PaintIt,
   Paint[TheTreeDiag, DisplayFunction->(Display[ProjectPath<>"Tree_qqbVH.ps",#]&) ];
];

TheTreeAmp   = CreateFeynAmp[ TheTreeDiag,   PreFactor-> 1 ];
Do[  TheTreeAmp[[i,3]]  *=DID[i] ,{i,1,Length[ TheTreeAmp]} ];

MyTreeAmp   = ExtractAmp[ TheTreeAmp ];
If[WriteAmp,
   AmpToFORM[ Join[MyTreeAmp] ,ProjectPath<>"Tree_qqbVH_input.frm"];
];

