(* ::Package:: *)

(*
   Copyright 2019 UniversalQCompiler (https://github.com/Q-Compiler/UniversalQCompiler)
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*)

(*Switch off precision warnings*)
Off[N::meprec]
Off[Floor::meprec]
Off[Eigensystem::eivec0]
(*Helper methods for tests for tests*)
isZero[x_]:=Chop[x]==0;
isZeroMatrix[mat_]:=Tr[Abs[Flatten[Chop[mat]]]]==0;
isDiagonal[m_]:=Norm[Chop[DiagonalMatrix[Diagonal[m]]-m]]==0;
isIdentity[m_]:=Norm[Chop[m-IdentityMatrix[Length[m]]]]==0;
isIdentityUpToPhase[m_]:=If[Chop[FullSimplify[Norm[m[[1,1]]]]]==0,False,
If[Chop[FullSimplify[Norm[m[[1,1]]]-1]]==0,
Chop[FullSimplify[Norm[m/(m[[1,1]])-IdentityMatrix[Length[m]]]]]==0,
False
]
];
IsEqualUptoPhaseForUnitaries[a_,b_]:=Module[{},
If[Dimensions[a]==Dimensions[b], ,Return[False]];
Return[isIdentityUpToPhase[a.ConjugateTranspose[b]]]
]
(*Create a random n\[Times]m permutation matrix*)
randPermutMat[m_,n_]:=Module[{perm},(
perm=IdentityMatrix[n][[All,RandomSample[Range[m]]]];
perm=perm[[;;,1;;m]];
perm
)]
(*Input dimension instead of the number of qubits*)
CreateIsometryFromListDim[st_,M_]:=CreateOperationFromGateList[st,Log2[M]]
NCreateIsometryFromListDim[st_,M_]:=NCreateOperationFromGateList[st,Log2[M]]
(*Multiply the gates in the list representation st to get the matrix representaiton of the operator represented by st. 
Note that n denotes the total number of qubits.
*)
MultiplyGates[st_,n_]:=CreateOperationFromGateList[st,n];
(*Numerically multiply the gates in the list representation st to get the matrix 
representaiton of the operator represented by st.  New version is faster than previous
*)
NMultiplyGates[st_,n_]:=NCreateOperationFromGateList[st,n]

(*Some functions for random numbers and matrices.
RR\[Rule]Random Real Number;
RC\[Rule]Random Complex Number;
RG[n]\[Rule]Random Complex Matrix of size n x n;
RU[n]\[Rule]Random Unitary Matrix of size n x n;
*)
RR:=RandomReal[NormalDistribution[0,1]];
RC:=RR+I*RR;
RG[n_]:=Table[RC,{n},{n}];
RU[n_]:=Orthogonalize[RG[n]];


(*Helpers for checking mathematical decompositions*)


checkIsoToUnitary[iso_]:=Module[{u},( 
u=IsoToUnitary[iso];
isIdentity[u.ConjugateTranspose[u]]&&isZero[Norm[iso-u[[All,1;;Dimensions[iso][[2]]]]]]
)]


(*Tests for single methods for mathematical decompositions*)


testAppendColums:=Module[{},( 
checkIsoToUnitary[PickRandomIsometry[2,4]]&&
checkIsoToUnitary[randPermutMat[2,4]]
)]


(*Unit tests for some basic methods*)


(*Helpers for checking basic methods*)


checkSimplifyGateLists[st_]:=Module[{u,stNew},( 
stNew=SimplifyGateList[st];
isIdentityUpToPhase[ConjugateTranspose[CreateOperationFromGateList[st]].CreateOperationFromGateList[stNew]]
)]


checkRzAngle[angle_]:=isIdentityUpToPhase[RzM[RzAngle[RzM[angle,1]],1].ConjugateTranspose[RzM[angle,1]]]


checkRxAngle[angle_]:=isIdentityUpToPhase[RxM[RxAngle[RxM[angle,1]],1].ConjugateTranspose[RxM[angle,1]]]


checkRyAngle[angle_]:=isIdentityUpToPhase[RyM[RyAngle[RyM[angle,1]],1].ConjugateTranspose[RyM[angle,1]]]


(*Tests for checking basic methods*)


testSimplifyGateLists:=Module[{},( 
If[
checkSimplifyGateLists[DecIsometry[PickRandomIsometry[2,4]]]&&
checkSimplifyGateLists[ColumnByColumnDec[N[randPermutMat[4,16]]]],
True
,
Print["Error in SimplifyGateLists"];
]
)]


testRzAngle:=Module[{},( 
If[
checkRzAngle[0]&&checkRzAngle[0.]&&checkRzAngle[Pi/2]&&checkRzAngle[Pi/2.]&&
checkRzAngle[3Pi]&&checkRzAngle[3.Pi]&&checkRzAngle[0.456]&&checkRzAngle[4]
,
True
,
Print["Error in RzAngle"];
]
)]


testRxAngle:=Module[{},( 
If[
checkRxAngle[0]&&checkRxAngle[0.]&&checkRxAngle[Pi/2]&&checkRxAngle[Pi/2.]&&
checkRxAngle[3Pi]&&checkRxAngle[3.Pi]&&checkRxAngle[0.456]&&checkRxAngle[4]
,
True
,
Print["Error in RxAngle"];
]
)]


testRyAngle:=Module[{},( 
If[
checkRyAngle[0]&&checkRyAngle[0.]&&checkRyAngle[Pi/2]&&checkRyAngle[Pi/2.]&&
checkRyAngle[3Pi]&&checkRyAngle[3.Pi]&&checkRyAngle[0.456]&&checkRyAngle[4]
,
True
,
Print["Error in RyAngle"];
]
)]


(*All tests for basic methods*)


testAllBasicMethods:=Module[{},( If[
testSimplifyGateLists&&
testRzAngle&&
testRxAngle&&
testRyAngle
,Print["All tests for the basic methods pass."],,
Print["testAllBasicMethods did neither return True nor False"]
]
)
]


(*Tests for (uniformly) multi controlled gates*)


(*Unit tests for MCGs*)


(*Helpers for checking the Decomposition of MCGs*)


(*Note: the matrix mat is required to be a unitary for this test.*)
checkApplyMCG[MCG_,mat_]:=Module[{mcg,invertMat,i,gates},
mcg=CreateMCG[MCG[[1]],Join[MCG[[2]],MCG[[3]]],MCG[[4]],MCG[[5]]];
(*Account for zero controls*)
gates={};
For[i=1,i<=MCG[[5]],i++,
If[MemberQ[MCG[[3]],i],
AppendTo[gates,RxM[Pi]],
AppendTo[gates,IdentityMatrix[2]];
];
];
invertMat=KroneckerProduct@@gates;
mcg=invertMat.mcg.invertMat;
isIdentityUpToPhase[mcg.mat.ConjugateTranspose[ApplyMCG[MCG,mat]]]
]


checkDecToffoliMultiControlUpToDiagonal[control_,target_,n_]:=Module[{i,st,mat},( 
st=DecToffoliMultiControlUpToDiagonal[control,target,n];
mat=NMultiplyGates[st,n];
isDiagonal[mat.ConjugateTranspose[CreateMCToffoli[control,target,n]]]
)]


checkDecToffoliMultiControl[control_,target_,n_]:=Module[{i,st,mat},( 
st=DecToffoliMultiControl[control,target,n];
mat=NMultiplyGates[st,n];
isIdentityUpToPhase[mat.ConjugateTranspose[CreateMCToffoli[control,target,n]]
]
)]


checkDecMCSpecialUnitary[su_,control_,target_,n_]:=Module[{i,st,mat},( 
st=DecMCSpecialUnitary[su,control,target,n];
mat=NMultiplyGates[st,n];
isIdentityUpToPhase[mat.ConjugateTranspose[CreateMCG[su,control,target,n]]]
)]


checkDecMCG[u_,control_,target_,n_]:=Module[{i,st,mat},( 
st=DecMCG[u,control,target,n];
mat=NMultiplyGates[st,n]; 
isIdentityUpToPhase[mat.ConjugateTranspose[CreateMCG[u,control,target,n]]]
)]


checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[u_,control_,target_,n_]:=Module[{i,st,mat},( 
st=DecMCSpecialUnitaryUpToDiagonal[u,control,target,n,ReturnDiagonal->True];
mat=NMultiplyGates[st,n];
FullSimplify[isIdentityUpToPhase[mat.ConjugateTranspose[CreateMCG[u,control,target,n]]]]
)]


checkDecMCSpecialUnitaryUpToDiagonal[u_,control_,target_,n_]:=Module[{i,st,mat},( 
st=DecMCSpecialUnitaryUpToDiagonal[u,control,target,n];
mat=NMultiplyGates[st,n];
FullSimplify[isDiagonal[mat.ConjugateTranspose[CreateMCG[u,control,target,n]]]]
)]


(*Tests for single methods for MGGs*)


testApplyMCG:=Module[{error},( 
error=0;
If[Quiet[Check[checkApplyMCG[{RU[2],{1,2,4},{},3,5},RU[2^5]],error=1;False]] &&
Quiet[Check[checkApplyMCG[{RU[2],{2},{3,5},1,5},RU[2^5]],error=2;False]],True,
Print["Error in ApplyMCG[] with error message code ",error];False
]
)]


testDecToffoliMultiControlUpToDiagonal:=Module[{error},( 
error=0;
If[Quiet[Check[checkDecToffoliMultiControlUpToDiagonal[{1,2,4},3,5],error=1;False]] &&
Quiet[Check[checkDecToffoliMultiControlUpToDiagonal[{1,2,4,5},3,7],error=2;False]],True,
Print["Error in DecToffoliMultiControlUpToDiagonal[] with error message code ",error];False
]
)]


testDecToffoliMultiControl:=Module[{error},( 
error=0;
If[Quiet[Check[checkDecToffoliMultiControl[{1,3},2,3],error=1;False]] &&
Quiet[Check[checkDecToffoliMultiControl[{1,2,5},3,5],error=2;False]]&&
Quiet[Check[checkDecToffoliMultiControl[{1,2,3,4},5,6],error=3;False]],True,
Print["Error in DecToffoliMultiControl[] with error message code ",error];False
]
)]


testDecMCSpecialUnitary:=Module[{M,error},( 
error=0;
M=RU[2];
M=M/Det[M]^(1/2);(*M is a random special unitary matrix*)
If[(*Quiet[Check[checkDecMCSpecialUnitary[N[-IdentityMatrix[2]],{1,2,3,4,5,6,8},7,8],error=1;False]] &&*)
Quiet[Check[checkDecMCSpecialUnitary[M,{1,2,3,4,5,6,7},8,8],error=2;False]],True,
Print["Error in DecMCSpecialUnitary[] with error message code ",error];False
]
)]


testDecMCG:=Module[{M,error},(
error=0; 
M=RU[2];
If[Quiet[Check[checkDecMCG[-IdentityMatrix[2],{1,2,3,4},6,6],error=1;False]] &&
Quiet[Check[checkDecMCG[M,{1,2,5},4,5],error=2;False]] &&
Quiet[Check[checkDecMCG[M,{3},1,3],error=3;False]]&&
Quiet[Check[checkDecMCG[M,{3,1},2,3],error=4;False]],True,
Print["Error in DecMCG[] with error message code ",error];False
]
)]


(*Here it is enough to test for less than 8 controls, since otherwise we use the method DecMCSpecialUnitary[], 
which is already tested above.*)
testDecMCSpecialUnitaryUpToDiagonalWithDiagonal:=Module[{M,error},( 
error=0; 
M=RU[2];
M=M/Det[M]^(1/2);(*M is a random special unitary matrix*)
If[Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[N[-IdentityMatrix[2]],{1,2,3,4},5,5],error=1;False]] &&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[M,{1,2,5},4,5],error=2;False]] &&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[M,{3},1,3],error=3;False]]&&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[M,{},1,3],error=4;False]]&&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[M,{3,1},2,3],error=5;False]]&&
(*Also check the method for analytic calculations:*)
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonalWithDiagonal[IdentityMatrix[2],{1,2,3},4,4],error=6;False]],
True,
Print["Error in DecMCGspecialUnitaryUpToDiagonal[] (with option returning the diagonal gate) with error message code ",error];False
]
)]


(*Here it is enough to test for less than 8 controls, since otherwise we use the method DecMCSpecialUnitary[], 
which is already tested above.*)
testDecMCGspecialUnitaryUpToDiagonal:=Module[{M,error},( 
error=0; 
M=RU[2];
M=M/Det[M]^(1/2);(*M is a random special unitary matrix*)
If[Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[N[-IdentityMatrix[2]],{1,2,3,4},5,5],error=1;False]] &&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[M,{1,2,5},4,5],error=2;False]] &&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[M,{3},1,3],error=3;False]]&&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[M,{},1,3],error=4;False]]&&
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[M,{3,1},2,3],error=5;False]]&&
(*Also check the method for analytic calculations:*)
Quiet[Check[checkDecMCSpecialUnitaryUpToDiagonal[N[IdentityMatrix[2]],{1,2,3},4,4],error=6;False]],
True,
Print["Error in DecMCGspecialUnitaryUpToDiagonal[] with error message code ",error];False
]
)]


(*Test all Methods for MCG*)


testAllMCGMethods:=Module[{},( If[
testApplyMCG&&
testDecToffoliMultiControlUpToDiagonal&&
testDecToffoliMultiControl&&
testDecMCG &&
testDecMCSpecialUnitary&&
testDecMCSpecialUnitaryUpToDiagonalWithDiagonal&&
testDecMCGspecialUnitaryUpToDiagonal
,Print["All tests for MCGs pass."],,
Print["testAllMCGMethods did neither return True nor False"]
]
)
]


(*Unit tests for UCGs*)


checkApplyUCG[UCG_,mat_]:=Module[{ucg},( 
ucg=CreateUCG[UCG[[1]],UCG[[2]],UCG[[3]],UCG[[4]]];
Chop[Norm[ucg.mat-ApplyUCG[UCG,mat]]]==0
)]


(*Helpers for checking the Decomposition of UCGs*)


checkDecUCGUpToDiagonalWithDiagonal[u_,control_,target_,n_]:=Module[{i,st,mat,dia},( 
st=DecUCGUpToDiagonal[u,control,target,n,ReturnDiagonal->True];
dia=st[[-1]];
mat=ApplyDiag[{dia[[2]],dia[[3]],n},MultiplyGates[Drop[st,-1],n]];
FullSimplify[isIdentityUpToPhase[mat.ConjugateTranspose[CreateUCG[u,control,target,n]]]]
)]


checkDecUCGUpToDiagonal[u_,control_,target_,n_]:=Module[{i,st,mat},( 
st=DecUCGUpToDiagonal[u,control,target,n];
mat=MultiplyGates[st,n];
FullSimplify[isDiagonal[mat.ConjugateTranspose[CreateUCG[u,control,target,n]]]]
)]


(*checks UC Z-rotations (set rotAxis=3) and Y-rotations (set rotAxis=2)*)


RotGate[\[Alpha]_,i_]:=
MatrixExp[I \[Alpha]/2 PauliMatrix[i]]
Options[checkDecUCZY] = {DecUCYWithCZ ->  False};
checkDecUCZY[angles_,control_,target_,n_,rotAxis_,OptionsPattern[]]:=Module[{i,st,mat,op,str,m},( 
If[rotAxis==3,
st=DecUCZ[angles,control,target,n],
If[OptionValue[DecUCYWithCZ],
st=DecUCY[angles,control,target,n,TwoQubitGate->"Cz"],
st=DecUCY[angles,control,target,n]
];
];
mat=MultiplyGates[st,n];
m={};
For[i=1,i<= Length[angles],i++,AppendTo[m,RotGate[angles[[i]],rotAxis]]];
isIdentityUpToPhase[mat.ConjugateTranspose[CreateUCG[m,control,target,n]]]
)]


(*Tests for methods for UGGs*)


testApplyUCG:=Module[{i,m0,m1,m2,m3,error},( 
error=0;
(*generate set with unitaries*)
m0={};
For[i=1,i<= 2^3,i++,AppendTo[m0,IdentityMatrix[2]]];
m1={};
For[i=1,i<= 2^4,i++,AppendTo[m1,N[IdentityMatrix[2]]]];
m2={{{0,1},{1,0}}};
For[i=1,i<= 2^3-1,i++,AppendTo[m2,RU[2]]];
If[Quiet[Check[checkApplyUCG[{m0,{1,2,3},4,4},RU[2^4]],error=1;False]] &&
Quiet[Check[checkApplyUCG[{m1,{1,2,3,5},4,5},RU[2^5]],error=2;False]] &&
Quiet[Check[checkApplyUCG[{m1,{1,2,5},4,5},RU[2^5]],error=3;False]]&&
Quiet[Check[checkApplyUCG[{m1,{1,5,2},4,5},RU[2^5]],error=4;False]]&&
Quiet[Check[checkApplyUCG[{m1,{5,1,2},4,6},RU[2^6]],error=5;False]]&&
Quiet[Check[checkApplyUCG[{m2,{2,3,4},1,6},RU[2^6]],error=6;False]],
True,
Print["Error in ApplyUCG[] with error message code ",error];False
]
)]


testDecUCG:=Module[{i,m0,m1,m2,m3,error},( 
error=0;
(*generate set with unitaries*)
m0={};
For[i=1,i<= 2^2,i++,AppendTo[m0,IdentityMatrix[2]]];
m1={};
For[i=1,i<= 2^4,i++,AppendTo[m1,N[IdentityMatrix[2]]]];
m2={{{0,1},{1,0}}};
For[i=1,i<= 2^3-1,i++,AppendTo[m2,RU[2]]];
(*
m3={};
For[i=1,i<= 2^7,i++,AppendTo[m3,RU[2]]];
*)
If[Quiet[Check[checkDecUCGUpToDiagonalWithDiagonal[m0,{1,3},2,3],error=1;False]] &&
Quiet[Check[checkDecUCGUpToDiagonalWithDiagonal[m1,{1,2,3,5},4,5],error=2;False]] &&
Quiet[Check[checkDecUCGUpToDiagonalWithDiagonal[m2,{2,4,3},1,6],error=3;False]] (*&&
Quiet[Check[checkDecUCGUpToDiagonalWithDiagonal[m3,{1,2,3,5,6,7,8},4,8],error=4;False]]
*),
True,
Print["Error in DecUCG[] with error message code ",error];False
]
)]


testDecUCGUpToDiagonal:=Module[{i,m0,m1,m2,m3,error},( 
error=0;
(*generate set with unitaries*)
m0={};
For[i=1,i<= 2^2,i++,AppendTo[m0,IdentityMatrix[2]]];
m1={};
For[i=1,i<= 2^4,i++,AppendTo[m1,N[IdentityMatrix[2]]]];
m2={{{0,1},{1,0}}};
For[i=1,i<= 2^3-1,i++,AppendTo[m2,RU[2]]];
(*
m3={};
For[i=1,i<= 2^7,i++,AppendTo[m3,RU[2]]];
*)
If[Quiet[Check[checkDecUCGUpToDiagonal[m0,{1,2},3,3],error=1;False]] &&
Quiet[Check[checkDecUCGUpToDiagonal[m1,{1,2,3,5},4,5],error=2;False]] &&
Quiet[Check[checkDecUCGUpToDiagonal[m2,{4,2,3},1,6],error=3;False]] (*&&
Quiet[Check[checkDecUCGUpToDiagonal[m3,{1,2,3,5,6,7,8},4,8],error=4;False]]
*),
True,
Print["Error in DecUCGUpToDiagonal[] with error message code ",error];False
]
)]


testDecUCZY:=Module[{i,m1,m2,error},( 
error=0;
(*generate angles for UCZ*)
m1={};
For[i=1,i<= 2^4,i++,AppendTo[m1,0]];
m2={};
For[i=1,i<= 2^4,i++,AppendTo[m2,RandomReal[]]];
If[Quiet[Check[checkDecUCZY[m1,{1,3,5,6},2,6,3],error=1;False]] && 
Quiet[Check[checkDecUCZY[m1,{1,3,5,6},2,6,2],error=2;False]]&& 
Quiet[Check[checkDecUCZY[m1,{6,3,5,1},2,6,2,DecUCYWithCZ ->  True],error=3;False]]&& 
Quiet[Check[checkDecUCZY[m2,{1,2,4,5},3,6,3],error=4;False]]&& 
Quiet[Check[checkDecUCZY[m2,{1,2,4,5},3,6,2],error=5;False]]&& 
Quiet[Check[checkDecUCZY[m2,{1,4,2,5},3,6,2,DecUCYWithCZ ->  True],error=6;False]],
True,
Print["Error in DecUCGZ[] with error message code ",error];False
]
)]


(*Test UGGs*)


testUCGs:=If[testApplyUCG&&testDecUCG &&testDecUCGUpToDiagonal&& testDecUCZY,Print["All tests for UCGs pass."],,Print["testUCGs did neither return True nor False"]]


(*Unit tests for diagonal gates*)


checkApplyDiag[diagGate_,mat_]:=Module[{diag},( 
diag=DiagMat[diagGate[[1]],diagGate[[2]],diagGate[[3]]];
Chop[Norm[diag.mat-ApplyDiag[diagGate,mat]]]==0
)
]


(*Helpers for checking decomposition of diagonal gates*)


checkDecDiagGate[dia_]:=Module[{st,m},( 
st=DecDiagGate[dia];
m=MultiplyGates[st,Log2[Length[dia]]];
isIdentityUpToPhase[ConjugateTranspose[DiagonalMatrix[dia]].m]
)
]


(*Test for Diagonal gates*)


testApplyDiag:=Module[{error},( 
error=0;
If[Quiet[Check[checkApplyDiag[{Exp[I RandomReal[{0,2Pi},2^1]],{2},2},RU[2^2]],error=1;False]]&&
Quiet[Check[checkApplyDiag[{Exp[I RandomReal[{0,2Pi},2^2]],{1,3},3},RU[2^3]],error=2;False]]&&
Quiet[Check[checkApplyDiag[{Exp[I RandomReal[{0,2Pi},2^3]],{3,1,4},5},IdentityMatrix[2^5]],error=3;False]]
,
True,
Print["Error in ApplyDiag[] with error message code ",error];False
]
)]


testDecDiagGate:=Module[{error},( 
error=0;
If[Quiet[Check[checkDecDiagGate[Exp[I RandomReal[{0,2Pi},2^1]]],error=1;False]]&&
Quiet[Check[checkDecDiagGate[Exp[I RandomReal[{0,2Pi},2^2]]],error=2;False]]&& 
Quiet[Check[checkDecDiagGate[Exp[I RandomReal[{0,2Pi},2^3]]],error=3;False]]&& 
Quiet[Check[checkDecDiagGate[ConstantArray[1,2^3]],error=4;False]],
True,
Print["Error in DecDiagGate[] with error message code ",error];False
]
)]


(*Test all methods for diagonal gates*)


testAllDiagGateMethods:=Module[{},( 
If[testDecDiagGate &&  testApplyDiag,Print["All tests for diagonal gates pass."],,
Print["testAllDiagGateMethods did neither return True nor False"]
]
)
]


(*Unit tests for column-by-column decomposition*)


(*Helpers for checking the CC-decomposition*)


checkColumnByColumnDec[v_]:=Module[{i,st,iso,st2,iso2},( 
st=ColumnByColumnDec[v];
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
st2=ColumnByColumnDec[v,FirstColumn->"StatePreparation"];
iso2=NCreateIsometryFromListDim[st2,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].iso]&&isIdentityUpToPhase[ConjugateTranspose[v].iso2]
)
]


checkColumnByColumnDecExact[v_]:=Module[{i,st,iso},( 
st=ColumnByColumnDec[v,Simp->False];
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].iso]
)
]


(*Test methods for CC-decomposition*)


testColumnByColumnDec:=Module[{error,v1,v2,v3,vExact1},( 
error=0;
v1=CreateOperationFromGateList[{Ancilla[0,1],Rx[1,2],CNOT[1,2],Rz[1,1],Rz[1,1]}];
v2=CreateOperationFromGateList[NGateList[{Ancilla[0,1],Rx[1,3],CNOT[3,1],Rz[1,3],Rz[1,1],CNOT[2,3],Ry[1,1]}]];
v3=CreateOperationFromGateList[NGateList[{Ancilla[0,1],Ancilla[0,4],Rx[1,3],CNOT[3,1],Rz[2.5,3],Rz[0.7,1],CNOT[2,3],Ry[0.4,1],Rz[0.4,1],Ry[\[Pi]/2,2],Rz[\[Pi],2],CNOT[2,4]}]];
vExact1=CreateOperationFromGateList[{Ancilla[0,1],Rx[Pi/4,2],Rx[Pi/4,1],CNOT[1,2],CNOT[3,2],CNOT[2,3],CNOT[3,1]}];
If[Quiet[Check[checkColumnByColumnDec[PickRandomIsometry[2^2,2^4]],error=1;False]]&&
Quiet[Check[checkColumnByColumnDec[PickRandomIsometry[2^3,2^3]],error=2;False]]&&
Quiet[Check[checkColumnByColumnDec[randPermutMat[2^1,2^2]],error=3;False]]&&
Quiet[Check[checkColumnByColumnDec[N[randPermutMat[2^2,2^3]]],error=4;False]]&&
Quiet[Check[checkColumnByColumnDec[PickRandomIsometry[1,2^2]],error=5;False]]&&
Quiet[Check[checkColumnByColumnDec[v1],error=6;False]]&&
Quiet[Check[checkColumnByColumnDec[v2],error=7;False]]&&
Quiet[Check[checkColumnByColumnDec[v3],error=8;False]]&&
Quiet[Check[checkColumnByColumnDecExact[vExact1],error=9;False]]&&
Quiet[Check[checkColumnByColumnDecExact[randPermutMat[8,8]],error=10;False]]
(*&&
Quiet[Check[checkColumnByColumnDec[PickRandomIsometry[2^3,2^3]],error=6;False]]&&
Quiet[Check[checkColumnByColumnDec[PickRandomIsometry[2^4,2^4]],error=7;False]]*),
True,
Print["Error in ColumnByColumnDec[] with error message code ",error];False
]
)]


(*Test CC-decomposition*)


testCCDec:=If[testDecUCGUpToDiagonal&&testColumnByColumnDec,Print["All tests for the column-by-column decomposition pass."],,
Print["ttestCCDec did neither return True nor False"]]


(*Unit tests for DecIsometry*)


(*Helpers for testing DecIsometry*)


checkDecIsometry[v_]:=Module[{i,st,iso,st2,iso2},( 
st=DecIsometry[v];
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].iso]
)
]


(*Tests for DecIsometry*)


testDecIsometry:=Module[{error,v1,v2,v3},( 
error=0;
v1=CreateOperationFromGateList[{Ancilla[0,1],Rx[Pi,2],CNOT[1,2],Rz[Pi,1],Rz[Pi,1]}];
v2=CreateOperationFromGateList[NGateList[{Ancilla[0,1],Rx[1,3],CNOT[3,1],Rz[1,3],Rz[1,1],CNOT[2,3],Ry[1,1]}]];
v3=CreateOperationFromGateList[NGateList[{Ancilla[0,1],Ancilla[0,4],Rx[1,3],CNOT[3,1],Rz[2.5,3],Rz[0.7,1],CNOT[2,3],Ry[0.4,1],Rz[0.4,1],Ry[\[Pi]/2,2],Rz[\[Pi],2],CNOT[2,4]}]];
If[Quiet[Check[checkDecIsometry[PickRandomIsometry[2^2,2^4]],error=1;False]]&&
Quiet[Check[checkDecIsometry[PickRandomIsometry[2^1,2^2]],error=2;False]]&&
Quiet[Check[checkDecIsometry[randPermutMat[2^1,2^2]],error=3;False]]&&
Quiet[Check[checkDecIsometry[N[randPermutMat[2^2,2^3]]],error=4;False]]&&
Quiet[Check[checkDecIsometry[PickRandomIsometry[1,2^2]],error=5;False]]&&
Quiet[Check[checkDecIsometry[PickRandomIsometry[2^4,2^4]],error=6;False]]&&
Quiet[Check[checkDecIsometry[v1],error=7;False]]&&
Quiet[Check[checkDecIsometry[v2],error=8;False]]&&
Quiet[Check[checkDecIsometry[v3],error=9;False]]
,
True,
Print["Error in DecIsometry[] with error message code ",error];False
]
)]


testIsometryDecompositions:=If[testDecIsometry,Print["All tests for DecIsometry pass."],,
Print["testIsometryDecompositions did neither return True nor False"]]


(*Unit tests for Isometries on a small number of qubits*)


(*Helpers for checking the methods for isoemtries on a small number of qubits*)


checkStatePrep1Qubit[s_]:=Module[{i,st,iso},( 
st=StatePrep1Qubit[s];
iso=NCreateIsometryFromListDim[st,Dimensions[s][[1]]];
isIdentityUpToPhase[ConjugateTranspose[s].iso]
)
]


checkStatePrep2Qubits[s_]:=Module[{i,st,iso},( 
st=StatePrep2Qubits[s];
iso=NCreateIsometryFromListDim[st,Dimensions[s][[1]]];
isIdentityUpToPhase[ConjugateTranspose[s].iso]
)
]


checkStatePrep3Qubits[s_,class_]:=Module[{i,st,iso},( 
st=StatePrep3Qubits[s];
iso=NCreateIsometryFromListDim[st,Dimensions[s][[1]]];
isIdentityUpToPhase[ConjugateTranspose[s].iso] &&
(CNOTCount[st]==class)
)
]


checkDecIso12[v_,analytic_:False]:=Module[{i,st,iso},( 
st=DecIso12[v];
If[analytic,
iso=CreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].N[iso]]
,
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].iso]
]
)
]


(*Test methods for isometries on a small number of qubits*)


testStatePrep1Qubit:=Module[{error},( 
error=0;
If[Quiet[Check[checkStatePrep1Qubit[PickRandomPsi[2]],error=1;False]]&&
Quiet[Check[checkStatePrep1Qubit[FPickRandomPsi[2,4]],error=2;False]]&&
Quiet[Check[checkStatePrep1Qubit[I*{{0},{1}}],error=3;False]]&&
Quiet[Check[checkStatePrep1Qubit[{{1},{0}}],error=4;False]]&&
Quiet[Check[checkStatePrep1Qubit[CreateOperationFromGateList[{Ancilla[0,1],Rz[1/3,1],Rz[1/3,1]}]],error=5;False]]
,
True,
Print["Error in statePrep1Qubit[] with error message code ",error];False
]
)]


testStatePrep2Qubits:=Module[{error,psi},( 
error=0;
If[Quiet[Check[checkStatePrep2Qubits[PickRandomPsi[4]],error=1;False]]&&
Quiet[Check[checkStatePrep2Qubits[FPickRandomPsi[4,4]],error=2;False]]&&
(psi = {\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
FractionBox["1", 
SqrtBox["2"]]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{"0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{"0"}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
FractionBox["1", 
SqrtBox["2"]]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)};
Quiet[Check[checkStatePrep2Qubits[psi],error=3;False]])&&
Quiet[Check[checkStatePrep2Qubits[CreateOperationFromGateList[{Ancilla[0,1],Ancilla[0,2],Rz[1/3,1],Rz[1/3,1],CNOT[2,1],Rz[1/3,1],Rz[1/3,2]}]],error=4;False]]
,
True,
Print["Error in statePrep2Qubits[] with error message code ",error];False
]
)]


testStatePrep3Qubits := Module[{state,TestBoolean, theta, ru},
TestBoolean = True;

(*class 0*)
state = Transpose[{Flatten[KroneckerProduct[PickRandomPsi[2],PickRandomPsi[2], PickRandomPsi[2]]]}];
If[checkStatePrep3Qubits[state,0], ,Print["3 qubit state preparation failed on a Class 0 state"]; TestBoolean = False];

state = NCreateIsometryFromListDim[{Ancilla[0,3],Ancilla[0,2],Ancilla[0,1],Rz[0.8,1],Rx[0.1,2],Rx[0.8,3],Rx[0.4,1],Ry[0.4, 3]}, 8];
If[checkStatePrep3Qubits[state,0], ,Print["3 qubit state preparation failed on a Class 0 state"]; TestBoolean = False];

(*class 1*)
state = Transpose[{Flatten[KroneckerProduct[PickRandomPsi[2],PickRandomPsi[4]]]}];
If[checkStatePrep3Qubits[state,1], ,Print["3 qubit state preparation failed on a Class 1 state"]; TestBoolean = False];

state = NCreateIsometryFromListDim[{Ancilla[0,3],Ancilla[0,2],Ancilla[0,1],Rz[0.8,1],Rx[0.1,2],Rx[0.8,3],Ry[0.8,1],Rz[0.1,2],CNOT[1, 2],Rx[0.4,1],Ry[0.4,3],Ry[0.8,1],Rx[0.1,2],Rx[0.8,3]},8];
If[checkStatePrep3Qubits[state,1], ,Print["3 qubit state preparation failed on a Class 1 state"]; TestBoolean = False];

(*class 2*)
 theta = RandomVariate[UniformDistribution[{0, \[Pi]/2}]];ru = RU[2]; state = Transpose[{Cos[theta] *Flatten[KroneckerProduct[ru[[All,1]],PickRandomPsi[2], PickRandomPsi[2]]]+  Sin[theta]*Flatten[KroneckerProduct[ru[[All,2]],PickRandomPsi[2], PickRandomPsi[2]]]}];
If[checkStatePrep3Qubits[state,2], ,Print["3 qubit state preparation failed on a Class 2 state, with the orthogonal factors on the first qubit"]; TestBoolean = False];

theta = RandomVariate[UniformDistribution[{0, \[Pi]/2}]];ru = RU[2]; state = Transpose[{Cos[theta] *Flatten[KroneckerProduct[PickRandomPsi[2],ru[[All,1]], PickRandomPsi[2]]]+  Sin[theta]*Flatten[KroneckerProduct[PickRandomPsi[2],ru[[All,2]], PickRandomPsi[2]]]}];
If[checkStatePrep3Qubits[state,2], ,Print["3 qubit state preparation failed on a Class 2 state, with the orthogonal factors on the second qubit"]; TestBoolean = False];

 theta = RandomVariate[UniformDistribution[{0, \[Pi]/2}]];ru = RU[2]; state = Transpose[{Cos[theta] *Flatten[KroneckerProduct[PickRandomPsi[2], PickRandomPsi[2],ru[[All,1]]]]+  Sin[theta]*Flatten[KroneckerProduct[PickRandomPsi[2], PickRandomPsi[2],ru[[All,2]]]]}];
If[checkStatePrep3Qubits[state,2], ,Print["3 qubit state preparation failed on a Class 2 state, with the orthogonal factors on the third qubit"]; TestBoolean = False];

state = NCreateIsometryFromListDim[{Ancilla[0,3],Ancilla[0,2],Ancilla[0,1],Rz[0.8,1],Rx[0.1,2],Rx[0.8,3],Ry[0.8,1],Rz[0.1,2],CNOT[1,2],Rx[0.4,2],Ry[0.4,3],CNOT[2,3],Ry[0.9,3],Ry[0.8,1],Rx[0.1,2],Rx[0.8,3]},8];
If[checkStatePrep3Qubits[state,2], ,Print["3 qubit state preparation failed on a Class 2 state"]; TestBoolean = False];

state = NCreateIsometryFromListDim[{Ancilla[0,3],Ancilla[0,2],Ancilla[0,1],Rz[0.8,1],Rx[0.1,2],Rx[0.8,3],Ry[0.8,1],Rz[0.1,2],CNOT[1, 2],Rx[0.4,1],Ry[0.4,3],CNOT[3,1],Ry[0.8,1],Rx[0.1,2],Rx[0.8,3]},8];
If[checkStatePrep3Qubits[state,2], ,Print["3 qubit state preparation failed on a Class 2 state"]; TestBoolean = False];

(*class 3*)
state = PickRandomPsi[8];
If[checkStatePrep3Qubits[state,3], ,Print["3 qubit state preparation failed on a random 3 qubit state (assumed to be class 3)"]; TestBoolean = False];

TestBoolean
]


testDecIso12:=Module[{error,v1},
error=0;
v1=CreateOperationFromGateList[{Ancilla[0,1],Rx[1,2],CNOT[1,2],Rz[1,1],Rz[1,1]}];
If[Quiet[Check[checkDecIso12[PickRandomIsometry[2,4]],error=1;False]]&&
Quiet[Check[checkDecIso12[v1,True],error=2;False]]
,
True,
Print["Error in DecIso12[] with error message code ",error];False
]
]


(*Test Isometries on a small number of qubits*)


testIsoSmall:=If[testStatePrep1Qubit && testStatePrep2Qubits && testStatePrep3Qubits && testDecIso12 ,Print["All tests for the isometries on a small number of qubits pass."],,
Print["testIsoSmall did neither return True nor False"]]


(*Unit tests for decomposition of single-qubit unitaries*)


(*Unit tests for decomposition of two-qubit unitaries*)


(*Helpers for testing ZYZDec and XYXDec*)


checkZYZDec[u_] := Module[{st},
st = ZYZDec[u,1];
If[isIdentityUpToPhase[MultiplyGates[st,1].ConjugateTranspose[u]], 
Return[True], 
Print["ZYZDec failed on input\n", MatrixFormOp[u]];
Return[False]
];
]


checkXYXDec[u_] := Module[{st},
st = XYXDec[u,1];
If[isIdentityUpToPhase[MultiplyGates[st,1].ConjugateTranspose[u]], 
Return[True], 
Print["ZYZDec failed on input\n", MatrixFormOp[u]]; 
Return[False]
];
]


(*Tests for ZYZDec and XYXDec*)


testZYZDec := Module[{isWorking},
isWorking=checkZYZDec[PickRandomIsometry[2,2]]&&checkZYZDec[{{0,1},{1,0}}]&&checkZYZDec[N[{{0,1},{1,0}}]]&&checkZYZDec[CreateOperationFromGateList[{Rz[Pi/2,1]}]]&&
checkZYZDec[CreateOperationFromGateList[{Ry[Pi/3,1]}]]&&checkZYZDec[CreateOperationFromGateList[{Ry[Pi/3,1]}]];
If[Not[isWorking],
Print["ZYZDec failed"]];
isWorking
];


testXYXDec := Module[{isWorking},
isWorking=checkXYXDec[PickRandomIsometry[2,2]]&&checkXYXDec[{{0,1},{1,0}}]&&checkXYXDec[N[{{0,1},{1,0}}]]
&&checkXYXDec[CreateOperationFromGateList[{Rz[Pi/3,1]}]]&&
checkXYXDec[CreateOperationFromGateList[{Ry[Pi/3,1]}]]&&checkXYXDec[CreateOperationFromGateList[{Ry[Pi/3,1]}]];
If[Not[isWorking],
Print["ZYZDec failed"]];
isWorking
];


(*All tests for single-qubit unitaries*)


testDecSingleQubit := Module[{out},
If[testZYZDec && testXYXDec
,
Print["All tests for ZYZDec and XYXDec pass"];
out=True;
,
Print["Error(s) in ZYZDec or XYXDec"];
out=False;
,
Print["testDecSingleQubit did neither return True nor False"]
];
out
]



(*Helpers for testing Dec2Qubit *)


checkDec2QubitGateUpToDiagonal[u_] := Module[{diag, st},
st = DecUnitary2Qubits[u,{1,2}, UpToDiagonal-> True];
If[isIdentityUpToPhase[NMultiplyGates[st,2].ConjugateTranspose[u]], 
Return[True], 
Print["Dec2Qubit with option UpToDiagonal->True failed on input\n", MatrixFormOp[u]]; 
Return[False]
];
]

(*Helpers for testing Dec2Qubit *not* up to a diagonal *)
checkDecUnitary2Qubits[u_,numCNOT_:Null] := Module[{diag, st,out,isCNOTCOUNTCorrect},
st = DecUnitary2Qubits[u,{1,2}, UpToDiagonal->False];
isCNOTCOUNTCorrect=If[numCNOT===Null,True,CNOTCount[st]==numCNOT];
If[isIdentityUpToPhase[NMultiplyGates[st,2].ConjugateTranspose[u]]&&isCNOTCOUNTCorrect, 
out=True, 
If[isIdentityUpToPhase[NMultiplyGates[st,2].ConjugateTranspose[u]],,Print["Dec2Qubit failed on input\n", MatrixForm[u]]];
 If[isCNOTCOUNTCorrect,,Print["Dec2Qubit outputs wrong number of C-NOT gates on input\n", MatrixForm[u]]];
out=False
];
out
]


(*Tests for the Method Dec2Qubit *)


testDec2QubitGateUpToDiagonalNoCnots := Module[{productGate1,productGate2,testResult1,testResult2},
productGate1 = KroneckerProduct[RU[2], RU[2]];
productGate2=KroneckerProduct[randPermutMat[2,2],randPermutMat[2,2]]; 
testResult1 = checkDec2QubitGateUpToDiagonal[productGate1];
testResult2 = checkDec2QubitGateUpToDiagonal[productGate2];
If[Not[testResult1&&testResult2],
Print["Dec2Qubit with option UpToDiagonal->True failed on an input requiring no cnot gates"]];
testResult1&&testResult2
];


testDec2QubitGateUpToDiagonalOneCnot := Module[{oneCnotGate1,oneCnotGate2,testResult1,testResult2},
oneCnotGate1 = KroneckerProduct[RU[2], RU[2]].CNOTM[1,2,2].KroneckerProduct[RU[2], RU[2]];
oneCnotGate2= KroneckerProduct[randPermutMat[2,2], randPermutMat[2,2]].CNOTM[1,2,2].KroneckerProduct[randPermutMat[2,2], randPermutMat[2,2]];
testResult1 = checkDec2QubitGateUpToDiagonal[oneCnotGate1];
testResult2= checkDec2QubitGateUpToDiagonal[oneCnotGate2];
If[Not[testResult1&&testResult2],
Print["Dec2Qubit with option UpToDiagonal->True failed on an input requiring one cnot gate"]
];
testResult1&&testResult2
];


testDec2QubitGateUpToDiagonalGeneric := Module[{genericGate,permGate,testResult1,testResult2},
genericGate = RU[4];
permGate=randPermutMat[4,4];
testResult1 = checkDec2QubitGateUpToDiagonal[genericGate];
testResult2=checkDec2QubitGateUpToDiagonal[permGate];
If[Not[testResult1&&testResult2],
Print["Dec2Qubit with option UpToDiagonal->True failed on a generic input gate"];
];
testResult1&&testResult2
];


testDec2QubitNoCnots := Module[{productGate,testResult,v1,v2},
productGate = KroneckerProduct[RU[2], RU[2]];
v1=CreateOperationFromGateList[NGateList[{Rx[Pi,2],Rz[Pi/4,1],Rz[Pi/3,2],Ry[Pi/6,2]}]];
v2=CreateOperationFromGateList[{Rx[Pi,2],Rz[Pi/4,1],Rz[Pi/3,2],Ry[Pi/6,2]}];
testResult = checkDecUnitary2Qubits[productGate,0]&&checkDecUnitary2Qubits[v1,0]&&checkDecUnitary2Qubits[v2,0];
If[Not[testResult],
Print["Dec2Qubit failed on an input requiring no cnot gates"];
];
testResult
];


testDec2QubitOneCnot := Module[{oneCnotGate,testResult,v1,v2},
oneCnotGate = KroneckerProduct[RU[2], RU[2]].CNOTM[1,2,2].KroneckerProduct[RU[2], RU[2]];
v1=CreateOperationFromGateList[NGateList[{Rx[Pi,2],CNOT[1,2],Rz[Pi/4,1],Rz[Pi/3,2],Ry[Pi/6,2]}]];
v2=CreateOperationFromGateList[{Ry[Pi/3,1],Ry[Pi,2],CNOT[2,1]}];
testResult = checkDecUnitary2Qubits[oneCnotGate,1]&&checkDecUnitary2Qubits[v1,1]&&checkDecUnitary2Qubits[v2,1];
If[Not[testResult],
Print["Dec2Qubit failed on an input requiring one cnot gate"];
];
testResult
];


testDec2QubitTwoCnot := Module[{oneCnotGate,testResult,v1,v2,v3,v3Prime},
v1=CreateOperationFromGateList[NGateList[{Rx[Pi,2],CNOT[1,2],Rz[Pi/4,1],Rz[Pi/3,2],CNOT[2,1],Ry[Pi/6,2]}]];
v2=CreateOperationFromGateList[NGateList[{Rx[Pi,2],CNOT[1,2],Rz[Pi/4,1],Rz[Pi/3,2],CNOT[1,2],Ry[Pi/6,2],Rz[Pi/2,2],Ry[Pi/3,1]}]];
v3Prime=CreateOperationFromGateList[{Ancilla[0,1],Ry[Pi/3,1],Ry[Pi,2],CNOT[2,1]}];
v3=AppendCols[v3Prime];
testResult = checkDecUnitary2Qubits[v1,2]&&checkDecUnitary2Qubits[v2,2]&&checkDecUnitary2Qubits[v3,2];
If[Not[testResult],
Print["Dec2Qubit failed on an input requiring two cnot gates"];
];
testResult
];


testDec2QubitGeneric := Module[{genericGate,testResult,v1,v2},
genericGate = RU[4];
v1=CreateOperationFromGateList[NGateList[{Rx[Pi,2],CNOT[1,2],Rz[Pi/4,1],Rz[Pi/3,2],CNOT[2,1],Ry[Pi/6,2],Ry[Pi/5,1],CNOT[2,1]}]];
v2=CreateOperationFromGateList[{CNOT[1,2],CNOT[2,1],CNOT[1,2]}];
testResult = checkDecUnitary2Qubits[genericGate,3]&&checkDecUnitary2Qubits[v1,3]&&checkDecUnitary2Qubits[v2,3];
If[Not[testResult],
Print["Dec2Qubit failed on a generic input gate"];
];
testResult
];


(*All tests for Dec2Qubit*)


testDec2Qubit := Module[{out},
If[testDec2QubitGateUpToDiagonalNoCnots && testDec2QubitGateUpToDiagonalOneCnot && testDec2QubitGateUpToDiagonalGeneric &&
testDec2QubitNoCnots && testDec2QubitOneCnot &&testDec2QubitTwoCnot && testDec2QubitGeneric
,
Print["All tests for Dec2Qubit pass"];
out=True;
,
Print["Error(s) in Dec2Qubit"];
out=False,
Print["testDec2Qubit did neither return True nor False"]
];
out
]



(*Unit tests for the Quantum Shannon Decompostion*)


(*Helpers for testing the Quantum Shannon Decompostion*)


checkQSD[testUnitary_,simp_:True]:= Module[{n,op, str,st,opV,opVCorrect,stV,stVCorrect},
st = QSD[testUnitary,Simp->simp];
stV =ConjugateTranspose[testUnitary].NCreateIsometryFromListDim[st,Dimensions[testUnitary][[1]]];
stVCorrect = isIdentityUpToPhase[stV];
If[Not[stVCorrect], 
Print["QSD (st) failed on input\n", MatrixForm[testUnitary]]; 
];
Return[stVCorrect];
]


(*Tests for Methods for the Quantum Shannon Decomposition*)


(*Test QSC*)
testQSD := Module[{vExact1,n,passedAllTestsSoFar,stV,stVCorrect,opV,opVCorrect,st,m},
passedAllTestsSoFar = True;
For[n=1, n < 4, n++,
If[Not[checkQSD[RU[2^n]]],
passedAllTestsSoFar = False;
Print["QSD failed on ", n, " qubit input"];
];
];
(*In the past the 8x8 identity matrix caused problems, so we check that it works*)
If[Not[checkQSD[N[IdentityMatrix[8]]]], Print["QSD (st) failed on input: N[IdentityMatrix[8]]"]; passedAllTestsSoFar = False;];
(*CheckQSD for exact input*)
vExact1=CreateOperationFromGateList[{Ancilla[0,1],Rx[Pi/4,2],Rx[Pi/2,1],CNOT[1,2],CNOT[3,2],CNOT[2,3],CNOT[3,1]}];
If[checkQSD[vExact1,False],,Print["QSD (st) failed on exact input v= ",vExact1]; passedAllTestsSoFar = False;];
(*In the past the 8x8 identity matrix caused problems, so we check that it works*)
If[FullSimplify[Not[checkQSD[IdentityMatrix[8]]]], Print["QSD (st) failed on input: IdentityMatrix[8]"]; passedAllTestsSoFar = False;];
(*Tests for isometries*)
For[n=2, n < 4, n++,
For[m=1, m <n , m++,
If[Not[checkQSD[PickRandomIsometry[2^m,2^n]]],
passedAllTestsSoFar = False;
Print["QSD failed on isometry from ", m, " to ",n, " qubits" ];
];
];
];
Return[passedAllTestsSoFar];
];


(*Test the Quantum Shannon Decomposition*)


testQSDAll:=If[testQSD,Print["All tests for the QSD pass."],,
Print["testQSDAll did neither return True nor False"]]


(*Unit tests for Stinespring*)


(*Helper methods for checks for Stinespring*)


Options[checkStinespring]={Exact->False};
checkStinespring[chan_,OptionsPattern[]] := 
 Module[{u, anc,uc,ancc,i, out = True,ch},
 {u, anc} = StinespringQubit[chan,TryToCompress->False];
 For[i = 1, i <= Dimensions[chan][[1]], i++, 
   out = And[out, 
     isZeroMatrix[KroneckerProduct[BraV[i - 1, 2^anc],IdentityMatrix[
           Dimensions[chan][[2]]]].u - chan[[i]]]]];
           out=And[out,isIdentity[ConjugateTranspose[u].u]];
  If[Not[OptionValue[Exact]],ch=MinimizeKrausRank[chan];
  {uc, ancc} = StinespringQubit[chan]; 
  For[i = 1, i <= Dimensions[ch][[1]], i++, 
   out = And[out, 
     isZeroMatrix[KroneckerProduct[BraV[i - 1, 2^ancc],IdentityMatrix[
           Dimensions[ch][[2]]]].uc - ch[[i]]]]];out=
           And[out,isIdentity[ConjugateTranspose[uc].uc]]];out]


(*Test of Stinespring*)


(* Note that for a channel from m to n qubits, the number of Kraus
   operators must be at least 2^(m-n) *)
testStinespring := 
 Module[{n, m, p, passedAllTestsSoFar, stV, stVCorrect, opV, 
   opVCorrect,chan},
  passedAllTestsSoFar = True;
  For[m=2,m<6,m++, 
   For[n=2,n<6,n++,For[p=Ceiling[2^(m-n)],p<5+Ceiling[2^(m-n)],p++,
     If[Not[checkStinespring[PickRandomChannel[2^m, 2^n, p]]],
      passedAllTestsSoFar = False;
      Print["StinespringQubit failed on a generic channel from ", m, 
       " qubits to ", n, " qubits with ", p, " Kraus operators"]]]]]; If[passedAllTestsSoFar, 
   Print["All tests for StinespringQubit on generic inputs pass"]];
   For[p=2,p<6,p++,
   chan=FPickRandomChannel[2,2,p,2];If[Not[checkStinespring[chan,Exact->True]],
      passedAllTestsSoFar = False;Print["StinespringQubit failed on an exact
channel from 1 qubit to 1 qubits with ", p, " Kraus operators"]]];
    If[passedAllTestsSoFar, 
   Print["All tests for StinespringQubit pass"],,
   Print["testStinespring  did neither return True nor False"]]]


(*Unit tests for POVMs*)


(*Helper method for checking POVMs*)


checkPOVM[POVM_] := 
 Module[{r,i,u,anc, out = True}, {u, anc} = POVMToIsometry[POVM]; 
 r=PickRandomRho[Dimensions[POVM][[2]]];
  For[i = 1, i <= Dimensions[POVM][[1]], i++, 
   out = And[out, 
     isZeroMatrix[MatrixPower[KroneckerProduct[BraV[i - 1, 2^anc],IdentityMatrix[
           Dimensions[POVM][[2]]]].u ,2]- POVM[[i]]]];
     out=And[out,isZero[Tr[KroneckerProduct[BraV[i - 1, 2^anc],IdentityMatrix[
           Dimensions[POVM][[2]]]].u.r.CT[u].KroneckerProduct[KetV[i - 1, 2^anc],IdentityMatrix[
           Dimensions[POVM][[2]]]]]-Tr[POVM[[i]].r]]]]; 
           And[out,isIdentity[ConjugateTranspose[u].u]]]


(*Test POVMs*)


testPOVM := 
 Module[{n, m, p, passedAllTestsSoFar, stV, stVCorrect, opV, 
   opVCorrect,chan},
  passedAllTestsSoFar = True;
  For[m=2,m<6,m++, 
   For[p=Ceiling[2^(m-n)],p<5+Ceiling[2^(m-n)],p++,
     If[Not[checkPOVM[PickRandomPOVM[2^m, p]]],
      passedAllTestsSoFar = False;
      Print["POVMToIsometry failed on a generic POVM on ", m, 
       " qubits with ", p, " elements"]]]]; If[passedAllTestsSoFar, 
   Print["All tests for POVMToIsometry on generic inputs pass"]];
   For[p=2,p<10,p++,
   If[Not[checkPOVM[Table[IdentityMatrix[4]/p,{i,1,p}]]],
      passedAllTestsSoFar = False;
      Print["POVMToIsometry failed on an exact POVM with ", p, " elements"]]];
    If[passedAllTestsSoFar, 
   Print["All tests for POVMToIsometry pass"],,
    Print["testPOVM did neither return True nor False"]]]


(*Unit tests for Knill decomposition*)


(*Check Methods used in Knill's decomposition*)


CheckXToY[x_,y_]:=Module[{q},
q = XToYTransform[x,y];
isIdentity[N[CT[q].q]]&&isZeroMatrix[N[q.x-y]]
]


CheckIsoToUnitarySpecial[v_]:=Module[{u,n,m},
{n,m}=Dimensions[v];
u = IsoToUnitarySpecial[v];
Return[isIdentity[CT[u].u]&&isZeroMatrix[u[[All,1;;m]]-v]&&(Count[Chop[Eigenvalues[u]-1],0]>=n-m)]
]


CheckUED[u_]:=Module[{angles,anglemats,vecs,length,q,p,n,i},
(*check vector output*)
n = Dimensions[u][[1]];
{angles,vecs} = UnitaryEigenvalueDecomp[u];
length = Dimensions[angles][[1]];
p= IdentityMatrix[n];
For[i=1,i<=length,i++,p=p.IsoToUnitary[vecs[[i]]].(
SparseArray[{1,1}->(Exp[I*angles[[i]]]-1),{n,n}]+IdentityMatrix[n]
).CT[IsoToUnitary[vecs[[i]]]]];
Return[isZeroMatrix[p-u]]
]


Options[checkKnillDec] = {UseDec -> "QSD"};
checkKnillDec[v_,OptionsPattern[]]:=Module[{i,st,iso},( 
st=KnillDec[v,UseDec->OptionValue[UseDec]];
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[N[ConjugateTranspose[v].iso]]
)
]


(*Test Methods used in Knill's decomposition*)


TestXToY := Module[{x,y,q,error,u},
(*Define fractional x and y:*)
x = {\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
RowBox[{
RowBox[{"2", " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
SqrtBox[
FractionBox["2", "13"]]}]},
{
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
SqrtBox[
FractionBox["2", "13"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
RowBox[{
RowBox[{"2", " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
SqrtBox[
FractionBox["2", "13"]]}]},
{
RowBox[{
RowBox[{
RowBox[{"-", "2"}], " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
SqrtBox[
FractionBox["2", "13"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
RowBox[{
RowBox[{
RowBox[{"-", "25"}], " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
FractionBox["1", 
SqrtBox["26"]]}]},
{
RowBox[{
RowBox[{"25", " ", 
SqrtBox[
FractionBox["2", "2821"]]}], "+", 
FractionBox["1", 
SqrtBox["26"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
RowBox[{
SqrtBox[
FractionBox["2", "13"]], "+", 
FractionBox["17", 
SqrtBox["5642"]]}]},
{
RowBox[{
SqrtBox[
FractionBox["2", "13"]], "-", 
FractionBox["17", 
SqrtBox["5642"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)};
y = {\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
RowBox[{
RowBox[{"7", " ", 
SqrtBox[
FractionBox["5", "598"]]}], "-", 
FractionBox["47", 
SqrtBox["648830"]]}]},
{
RowBox[{
RowBox[{"7", " ", 
SqrtBox[
FractionBox["5", "598"]]}], "+", 
FractionBox["47", 
SqrtBox["648830"]]}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
FractionBox[
RowBox[{"4960", "-", 
RowBox[{"647", " ", 
SqrtBox["217"]}]}], 
RowBox[{"31", " ", 
SqrtBox["1019590"]}]]},
{
FractionBox[
RowBox[{"4960", "+", 
RowBox[{"647", " ", 
SqrtBox["217"]}]}], 
RowBox[{"31", " ", 
SqrtBox["1019590"]}]]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
FractionBox[
RowBox[{"808", "+", 
RowBox[{"131", " ", 
SqrtBox["217"]}]}], 
RowBox[{"217", " ", 
SqrtBox["1430"]}]]},
{
FractionBox[
RowBox[{
RowBox[{"-", "808"}], "+", 
RowBox[{"131", " ", 
SqrtBox["217"]}]}], 
RowBox[{"217", " ", 
SqrtBox["1430"]}]]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\),\!\(\*
TagBox[
RowBox[{"(", 
TagBox[GridBox[{
{
FractionBox[
RowBox[{"3", " ", 
RowBox[{"(", 
RowBox[{"7", "+", 
RowBox[{"3", " ", 
SqrtBox["217"]}]}], ")"}]}], 
RowBox[{"7", " ", 
SqrtBox["910"]}]]},
{
FractionBox[
RowBox[{"3", " ", 
RowBox[{"(", 
RowBox[{"7", "-", 
RowBox[{"3", " ", 
SqrtBox["217"]}]}], ")"}]}], 
RowBox[{"7", " ", 
SqrtBox["910"]}]]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.5599999999999999]}, Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}, "RowsIndexed" -> {}, "Items" -> {}, "ItemsIndexed" -> {}}],
Column], ")"}],
Function[BoxForm`e$, MatrixForm[BoxForm`e$]]]\)};
If[Quiet[Check[CheckXToY[x,y],error=1;False]]&&
(x=PickRandomIsometry[16,64];u = PickRandomUnitary[Dimensions[x][[1]]];y=u.x;
Quiet[Check[CheckXToY[x,y],error=2;False]])&&
(x=PickRandomIsometry[25,50];u = PickRandomUnitary[Dimensions[x][[1]]];y=u.x;
Quiet[Check[CheckXToY[x,y],error=3;False]])
,
True,
Print["Error in XToY with error message code ",error];False
]
]


TestIsoToUnitarySpecial:=Module[{error},( 
error=0;
If[Quiet[Check[CheckIsoToUnitarySpecial[PickRandomIsometry[16,64]],error=1;False]]&&
Quiet[Check[CheckIsoToUnitarySpecial[PickRandomIsometry[25,50]],error=1;False]]
,
True,
Print["Error in IsoToUnitarySpecial with error message code ",error];False
]
)]


TestUnitaryEigenvalueDecomp:=Module[{error},( 
error=0;
If[Quiet[Check[CheckUED[IsoToUnitarySpecial[PickRandomIsometry[16,64]]],error=1;False]]&&
Quiet[Check[CheckUED[IsoToUnitarySpecial[PickRandomIsometry[25,50]]],error=1;False]]
,
True,
Print["Error in UnitaryEigenvalueDecomp with error message code ",error];False
]
)]


TestKnillDec1:=Module[{error},( 
error=0;
If[Quiet[Check[checkKnillDec[PickRandomIsometry[2,2^4]],error=1;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[1,2^3]],error=2;False]]&&
Quiet[Check[checkKnillDec[randPermutMat[2^1,2^2]],error=3;False]]&&
Quiet[Check[checkKnillDec[N[randPermutMat[2^2,2^4]]],error=4;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[1,2^2]],error=5;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[2^3,2^3]],error=6;False]],
True,
Print["Error in KnillDec[] (using the QSD for decomposing the unitaries) with error message code ",error];False
]
)]


TestKnillDec2:=Module[{error,wrongDecomposition},( 
error=0;
If[Quiet[Check[checkKnillDec[PickRandomIsometry[2,8],UseDec->"KnillDec"],error=1;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[8,8],UseDec->"KnillDec"],error=2;False]]&&
Quiet[Check[checkKnillDec[randPermutMat[2^1,2^2],UseDec->"KnillDec"],error=3;False]]&&
Quiet[Check[checkKnillDec[N[randPermutMat[2^2,2^3]],UseDec->"KnillDec"],error=4;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[1,2^2],UseDec->"KnillDec"],error=5;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[2^3,2^3],UseDec->"KnillDec"],error=6;False]],
True,
Print["Error in KnillDec[] (using the Knill decomposition for decomposing the unitaries) with error message code ",error];False
]
)]


TestKnillDec3:=Module[{error,wrongDecomposition},( 
error=0;
If[Quiet[Check[checkKnillDec[PickRandomIsometry[2,8],UseDec->"ColumnByColumnDec"],error=1;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[8,8],UseDec->"ColumnByColumnDec"],error=2;False]]&&
Quiet[Check[checkKnillDec[randPermutMat[2^1,2^2],UseDec->"ColumnByColumnDec"],error=3;False]]&&
Quiet[Check[checkKnillDec[N[randPermutMat[2^2,2^3]],UseDec->"ColumnByColumnDec"],error=4;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[1,2^2],UseDec->"ColumnByColumnDec"],error=5;False]]&&
Quiet[Check[checkKnillDec[PickRandomIsometry[2^3,2^3],UseDec->"ColumnByColumnDec"],error=6;False]],
True,
Print["Error in KnillDec[] (using the ColumnByColumnDec for decomposing the unitaries) with error message code ",error];False
]
)]


FTestKnillDec := Module[{diag,exactDiagonalMat,bellIso,bellUnit,directUnit,error},
(*tests with fractional (exact) inputs to Knill decomp. based on diagonal/ block diagonal inputs(based on Bell states). knill decomp does NOT work with arbitrary exact inputs*)
(*Construct matrices:*)
diag = {Exp[I*RandomInteger[15]*Pi/15]};
Do[AppendTo[diag,Exp[I*RandomInteger[15]*Pi/15]],{i,3}];
exactDiagonalMat = DiagonalMatrix[diag];
bellIso = Transpose[{{1,0,0,1},{0,1,1,0}}/2^(1/2)];
bellUnit= {{1,0,0,-1},{0,1,1,0},{0,-1,1,0},{1,0,0,1}}/2^(1/2);
For[i = 1, i <= 10, i++, directUnit =  
    FullSimplify[
      DirectSum[FPickRandomUnitary[2, 5], FPickRandomUnitary[2, 5]]]; 
 If[directUnit.CT[directUnit] == IdentityMatrix[4], Break[]]]; If[
 i == 11, Print[
  "Warning: Fractional Knill check on non-unitary input"]];
error=0;
If[Quiet[Check[checkKnillDec[exactDiagonalMat],error=1;False]]&&
Quiet[Check[checkKnillDec[bellIso],error=2;False]]&&
Quiet[Check[checkKnillDec[bellUnit],error=3;False]]&&
Quiet[Check[checkKnillDec[directUnit],error=4;False]],
True,
Print["Error in KnillDec[] for exact inputs with error message code ",error];False
]
]



(*All tests for Knill's decomposition*)


testKnill := Module[{},(
If[TestXToY &&TestIsoToUnitarySpecial&&TestUnitaryEigenvalueDecomp&&TestKnillDec1
&&TestKnillDec2&&TestKnillDec3&&FTestKnillDec,
Print["All tests for Knill's decomposition pass"],,
  Print["testKnill did neither return True nor False"]
]
)]


(*Unit tests for state preparaion (Plesch-Brukner decomposition)*)


(*Checks for state preparation*)


CheckStatePreparation[v_,recLevel_]:=Module[{st,iso},
st=StatePreparation[v,Range[Length[v]],recLevel];
iso=NCreateIsometryFromListDim[st,Dimensions[v][[1]]];
isIdentityUpToPhase[ConjugateTranspose[v].iso]
]


(*Tests for state preparation*)


TestStatePreparation:=Module[{error},( 
error=0;
If[Quiet[Check[CheckStatePreparation[PickRandomPsi[2^4],level->1],error=1;False]]&&
Quiet[Check[CheckStatePreparation[PickRandomPsi[2^4],level->2],error=2;False]]&&
Quiet[Check[CheckStatePreparation[PickRandomPsi[2^4],level->3],error=3;False]]&&
Quiet[Check[CheckStatePreparation[N[randPermutMat[1,8]],level->1],error=1;False]]&&
Quiet[Check[CheckStatePreparation[N[randPermutMat[1,8]],level->2],error=2;False]]&&
Quiet[Check[CheckStatePreparation[N[randPermutMat[1,8]],level->3],error=3;False]]
,
True,
Print["Error in StatePrepRecursive with error message code ",error];False
]
)]


(*All tests for state preparation*)


testStatePreparationAll := Module[{},(
If[TestStatePreparation,
Print["All tests for State preparation pass"],,
 Print["testStatePreparationAll did neither return True nor False"]
]
)]


CheckCNOTtoXX[v_]:=Module[{st,st2,iso,ch},
st=DecIsometryGeneric[v];st2=CNOTRotationsToXXRGates[st];
iso=CreateOperationFromGateList[st2,Log[2,Dimensions[v][[1]]]];ch=CT[v].iso;
Chop[ch/ch[[1,1]]-IdentityMatrix[Dimensions[iso][[2]]],10^-6]==0*IdentityMatrix[Dimensions[iso][[2]]]]

CheckXXtoCNOT[v_]:=Module[{st,st2,st3,iso,ch},
st=DecIsometryGeneric[v];st2=CNOTRotationsToXXRGates[st];
st3=XXRGatesToCNOTRotations[st2];
iso=CreateOperationFromGateList[st3,Log[2,Dimensions[v][[1]]]];ch=CT[v].iso;
Chop[ch/ch[[1,1]]-IdentityMatrix[Dimensions[iso][[2]]],10^-6]==0*IdentityMatrix[Dimensions[iso][[2]]]]

TestCNOTtoXX:=Module[{error=0},If[Quiet[Check[CheckCNOTtoXX[PickRandomIsometry[2^2,2^2]],error=1;False]]&&
Quiet[Check[CheckCNOTtoXX[PickRandomIsometry[2^2,2^3]],error=2;False]]&&
Quiet[Check[CheckCNOTtoXX[PickRandomIsometry[2^2,2^4]],error=3;False]]&&
Quiet[Check[CheckCNOTtoXX[PickRandomIsometry[2^3,2^3]],error=4;False]]&&
Quiet[Check[CheckCNOTtoXX[PickRandomIsometry[2^3,2^4]],error=5;False]]&&
Quiet[Check[CheckCNOTtoXX[randPermutMat[4,4]],error=6;False]]&&
Quiet[Check[CheckCNOTtoXX[N[randPermutMat[8,8]]],error=7;False]]
,
True,
Print["Error in CheckCNOTtoXX with error message code ",error];False
]]

TestXXtoCNOT:=Module[{error=0},If[Quiet[Check[CheckXXtoCNOT[PickRandomIsometry[2^2,2^2]],error=1;False]]&&
Quiet[Check[CheckXXtoCNOT[PickRandomIsometry[2^2,2^3]],error=2;False]]&&
Quiet[Check[CheckXXtoCNOT[PickRandomIsometry[2^2,2^4]],error=3;False]]&&
Quiet[Check[CheckXXtoCNOT[PickRandomIsometry[2^3,2^3]],error=4;False]]&&
Quiet[Check[CheckXXtoCNOT[PickRandomIsometry[2^3,2^4]],error=5;False]]&&
Quiet[Check[CheckXXtoCNOT[randPermutMat[4,4]],error=6;False]]&&
Quiet[Check[CheckXXtoCNOT[N[randPermutMat[8,8]]],error=7;False]]
,
True,
Print["Error in CheckXXtoCNOT with error message code ",error];False
]]

testXXCNOTAll := Module[{},(
If[TestCNOTtoXX&&TestXXtoCNOT,
Print["All tests for converting to and from XX gates pass"],,
 Print["testXXCNOTAll did not return True or False"]
]
)]


(*Run all tests*)
runAllTests:=Module[{},(testAllBasicMethods;testUCGs;testAllDiagGateMethods;testIsoSmall;testCCDec;testDec2Qubit;testDecSingleQubit;testQSDAll;testQSD;testStatePreparationAll;testAllMCGMethods;testKnill;testIsometryDecompositions;testStinespring;testPOVM;testXXCNOTAll)]


Timing[runAllTests]


(*Switch precision warnings on again*)
On[N::meprec]
On[Floor::meprec]
On[Eigensystem::eivec0]
