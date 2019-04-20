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

BeginPackage["UniversalQCompiler`",{"QI`"}];

(*ToDos:
1. Add "do not simplify" option to methods, which allows to work with fixed circuit topologies.
2. Clear up unncessesary applications of Reverse[] (restructure the code to make it more redable)
3. Improve efficiency by removing For-loops and if statements with more efficient alternativs provided by Mathematica (and do not copy (long)lists)
4. Add ZXZ and XZX decompositions for single-qubit unitaries and use them for simplifying gate sequences 
(for instance if one of the Zs in ZXZ is actually identity, the X might commute out of a neighbouring target gate)
5. Create code for transforming gate sequences consisting of C-NOT and single-qubit gates into ones consisting of MS gates and single qubit rotations adapted for trapped ions.
6. Implementation of multi-controlled-Toffoli gates using ancillas to lower the C-NOT count.
*)

(*Methods to handle and simplify gate sequence*)
(*CreateIsometryFromList::usage="CreateIsometryFromList[st,(n)] creates the operator corresponding to the list (optionally, the total number of qubits n can be determined)."
NCreateIsometryFromList::usage="CreateIsometryFromList[st,(n)] creates the operator (with numerical numbers) corresponding to the list (optionally, the total number of qubits n can be determined)."
CreateChannelFromList::usage="CreateChannelFromList[st, (n)] multiplies the gates in the list st on n qubits and outputs the channel represented."
NCreateChannelFromList::usage="NCreateChannelFromList[st, (n)] multiplies the gates in the list st on n qubits numericallyand outputs the channel represented."
CreateInstrumentFromList::usage="CreateInstrumentFromList[st, (n)] multiplies the gates in the list st on n qubits and outputs the instrument represented."
NCreateInstrumentFromList::usage="NCreateInstrumentFromList[st, (n)] multiplies the gates in the list st on n qubits numericallyand outputs the instrument represented."
*)
CreatePOVMFromGateList::usage="CreatePOVMFromGateList[st, (n)] multiplies the gates in the list st on n qubits and outputs the POVM represented (optionally, the total number of qubits n can be determined)."
NCreatePOVMFromGateList::usage="NCreatePOVMFromGateList[st, (n)] multiplies the gates in the list st on n qubits numerically and outputs the POVM represented (optionally, the total number of qubits n can be determined)."
CreateOperationFromGateList::usage="CreateOperationFromGateList[st, (n)] generates the isometry/channel/instrument represented by st (optionally, the total number of qubits n can be determined)."
NCreateOperationFromGateList::usage="CreateOperationFromGateList[st, (n)] generates the isometry/channel/instrument represented by st numerically (optionally, the total number of qubits n can be determined)."
NGateList::usage="NGateList[st] makes the gate parameters numerical. [Helpful to speed up further processing of the gate list.]"
RelabelQubits::usage="RelabelQubits[st,numIn,numOut] relabels qubits in the gate sequence st (in list form) with qubit number in the list numIn with the qubit numbers given in the list numOut "
CNOTCount::usage="CNOTCount[st] returns the number of C-NOT gates in the gate sequence st."
InverseGateList::usage="InverseGateList[st] takes the inverse of a gate sequence st in list format."
AdjustAngle::usage="AdjustAngle[st] adjusts the angles in a gate sequence st to lie in [0,2 \[Pi]). Note that this may change the global phase by a factor -1."
ListFormToStr::usage="ListFormToStr[st] transforms a list in list format to a list containing the corresponding string representations of the gates."
SimplifyGateList::usage="SimplifyGateList[st] simplifies a gate sequence on in list format st by merging and commuting single-qubit and C-NOT gates (see the full documentation for more details)."
NSimplifyGateList::usage"NSimplifyGateList[st] simplifies a gate sequence (numerically) in list format st by merging and commuting single-qubit and C-NOT gates (see the full documentation for more details)."
NumberOfQubits::usage"NumberOfQubits[st] finds the largest qubit number featured in the list format st ."
(*Create gates in list form*)
CNOT::usage="CNOT[i,j] creates a C-NOT gate in list form with control qubit i and target qubit j."
CZ::usage="CZ[i,j] creates a controlled-Z gate in list form with control qubit i and target qubit j."
XX::usage="XX[phi,i,j] creates an XX gate in list form with parameter phi acting on qubits i,j."
RGate::usage="RGate[theta,phi,i] creates an R gate in list form with parameters theta,phi acting on qubit i."
Diag::usage="Diag[entr,act] creates a  is a diagonal gate with diagonal entries given as a list entr and acting on the qubits listed in act."
Rx::usage="Rx[angle,act] creates an Rx gate in list form on the qubit act and with angle angle."
Ry::usage="Ry[diag,act] creates an Ry gate in list form on the qubit act and with angle angle."
Rz::usage="Rz[diag,act] creates an Rz gate in list form on the qubit act and with angle angle."
Mmt::usage="Mmt[act] creates a measurement in list form on the qubit act."
TrOut::usage="TrOut[act] creates a tracing out operation in list form on the qubit act."
Ancilla::usage="Ancilla[i_,act_] indicates that the qubit with number act is an ancilla starting in state |i> with i=0,1."
PostSelect::usage="PostSelect[i_,act_] indicates that the qubit with number act is postselected to state |i> with i=0,1."
(*Transform gates to matrices*)
MatrixFormOp::usage="MatrixFormOp[op] prints all the matrices in the list op in matrix form."
GateTypes::usage="GateTypes[] prints the convention for gate types used in this package."
ListFormToOp::usage="ListFormToOp[st,(n)] transforms a sequence of gates given in list from to a list containing the matrix repressentations of the gates (optionally, the total number of qubits n can be determined)."
RxM::usage="RxM[\[Alpha],i,n] creates a rotation gate around the x axis with angle \[Alpha] on qubit with number i, where the total number of qubits is n."
RyM::usage="RyM[\[Alpha],i,n] creates a rotation gate around the y axis with angle \[Alpha] on qubit with number i, where the total number of qubits is n."
RzM::usage="RzM[\[Alpha],i,n] creates a rotation gate around the z axis with angle \[Alpha] on qubit with number i, where the total number of qubits is n."
CNOTM::usage="CNOTM[i,j,n] creates a C-NOT gate on n qubits with control on the ith qubit and target on the jth qubit." 
CZM::usage="CZM[i,j,n] creates a controlled-Z gate on n qubits with control on the ith qubit and target on the jth qubit." 
XXM::usage="XXM[phi,{i,j},n] creates an XX gate with parameter phi acting on the i,j qubits with n qubits in total." 
RGateM::usage="RGateM[theta,phi,i,n] creates an R-gate with parameters theta,phi acting on the ith qubit with n qubits in total." 
DiagMat::usage="DiagMat[diag,act,n] creates a diagonal matrix with diagonal entries in the list diag representing a diagonal gate on n qubits acting on the qubits in the list act." 
RzAngle::usage="RzAngle[Rz(\[Theta])] extracts the rotation angle \[Theta] from a rotation gate Rz (\[Theta])."
RyAngle::usage="RyAngle[Ry(\[Theta])] extracts the rotation angle \[Theta] from a rotation gate Ry (\[Theta])."
RxAngle::usage="RxAngle[Rx(\[Theta])] extracts the rotation angle \[Theta] from a rotation gate Rx (\[Theta])."
(*Apply gates efficiently*)
ApplyDiag::usage="TBA."
ApplyMCG::usage="TBA."
ApplyUCG::usage="TBA."
(*Visualization of ciruits*)
PrintCircuit::usage="PrintCircuit[st, (n)] prints a circuit st given in list form. Optionally, one can determine the total number of qubits n."
LatexQCircuit::usage="LatexQCircuit[st, (n)] provides an circuit in the form of the Latex package QCircuit. Optionally, one can determine the total number of qubits n."
(*Matrix decompositions*)
IsoToUnitary::usage="IsoToUnitary[v] expands an isometry v to an unitary by appending additional columns."
(*KronFactor::usage="KronFactor[m=a\[CircleTimes]b] for a 4\[Cross]4 matrix m (and two 2\[Cross]2 matrices a,b) returns c*a and 1/c*b for some constant c."*)
KronFactorUnitaryDim4::usage="KronFactorUnitaryDim4[u=a\[CircleTimes]b] for a 4\[Cross]4 unitary u (and two 2\[Cross]2 unitaries a,b) returns c*a and 1/c*b for some constant c."
KronFactorVector::usage="KronFactorVector[v=VecProd[a,b]] for a four-dimensional vector v, returns two two-dimensional vectors c*a and 1/c*b for some constant c."
VecProd::usage="VecProd[a,b] generates the tensor product of a,b."
XToYTransform::usage="TBA."
UnitaryEigenvalueDecomp::usage="TBA."
IsoToUnitarySpecial::usage="IsoToUnitarySpecial[v] extends the isometry v to a unitary such that the unitary has as many eigenvalues equal to 1 as possible."
(*Matrix decomposition*)
ZYZDecomposition::usage="ZYZDecomposition[u] for an single-qubit unitary u, returns {a,b,c,d} such that u=e^{i d}.R_z(c).R_y(b).R_z(a)."
XYXDecomposition::usage="XYXDecomposition[u] for an single-qubit unitary u, returns {a,b,c,d} such that u=e^{i d}.R_x(c).R_y(b).R_x(a)."
RQDecomposition::usage="RQDecomposition[m] for a square matrix m returns {r,q}, where r is an upper triangular matrix and q is an orthogonal/unitary matrix such that r.ConjugateTranspose[q] = m."
QLDecomposition::usage="QLDecomposition[m] for a square matrix m returns {q,l}, where l is a lower triangular matrix and q is an orthogonal/unitary matrix such that ConjugateTranspose[q].l = m."
CSD::usage="CSD[u] returns {m1,m2,m3}, the Cosine-Sine decomposition of u, i.e., u=m1.m2.m3."
(*Decompositions for single-qubit gates*)
ZYZDec::usage="ZYZDec[u,(action)] for an single-qubit unitary u acting on qubit action (default: action=1), returns the ZYZ decompostition in list form."
XYXDec::usage="XYXDec[u,(action)] for an single-qubit unitary u acting on qubit action (default: action=1), returns the XYX decompostition in list form."
(*Uniformly controlled rotations*)
DecUCY::usage="See the documentation."
DecUCZ::usage="See the documentation."
(*Quantum Shannon Decomposition*)
DemultiplexUC::usage="See the documentation."
QSD::usage="QSD[v,(action)] decomposes an isometry v from m qubits to n qubits into single-qubit gates and C-NOTs. Optionally, one can determine the qubits the isometry is acting on by providing their numers in teh list action (default: action=Range[n])."
(*Decompose diagonal gate*)
DecDiagGate::usage="DecDiagGate[{\!\(\*SubscriptBox[\(d\), \(1\)]\),...,\!\(\*SubscriptBox[\(d\), \(2^n\)]\)},action] decomposes a diagonal gate on n qubits listed in action (default:action=Range[n]) with diagonal entries \!\(\*SubscriptBox[\(d\), \(1\)]\),...,\!\(\*SubscriptBox[\(d\), \(2^n\)]\)."
(*Decompose two-qubit gates*)
DecUnitary2Qubits::usage="DecUnitary2Qubits[u,(action={n1,n2})] decomposes a two-qubit gate u acting on the two qubits n1 and n2 (default: action={1,2}) (see the full documentation for further options)."
(*Decompositions multi controlled Toffoli gates*)
DecToffoliMultiControlUpToDiagonal::usage="TBA."
DecToffoliUpToDiagonal::usage="TBA."
DecToffoli::usage="TBA."
DecToffoliMultiControl::usage="DecToffoliMultiControl[control,target,n] decomposes a controlled NOT gate on n>=3 qubits controlling on the qubits in the list control and acting on the target qubit target."
(*Decompositions multi controlled single-qubit gates*)
DecMCG::usage="DecToffoliMultiControl[u,control,target,n] decomposes a controlled single-unitary gate u on n qubits controlled on qubits in the list control and acting on the target qubit target."
DecMCSpecialUnitary::usage="TBA."
DecMCSpecialUnitaryUpToDiagonal::usage="TBA."
(*Methods to create and handle multi controlled gates*)
CreateMCToffoli::usage="CreateMCToffoli[control,target,n] creates a matrix corresponding to a multi controlled Toffoli gate on n qubits with controls given in the list control and target qubit nubmer target."
CreateMCG::usage="TBA."
CreateUCG::usage="TBA."
ExtractMCG::usage="TBA."
ExtractUCG::usage="TBA."
(*Methods to decompose uniformly controlled gates*)
DecUCGUpToDiagonal::usage="TBA."
(*Column-by-column decomposition*)
ColumnByColumnDec::usage="ColumnByColumnDec[v,(action)] decomposes the isometry v from m to n qubits on the qubits listed in action (default: action=Range[n]) into single-qubit rotations and C-NOT gates (see the full documentation for additional options)."
(*Decomposition for isometries on a small number of qubits*)
StatePrep1Qubit::usage="TBA."
StatePrep2Qubits::usage="TBA."
StatePrep3Qubits::usage="TBA."
DecIso12::usage="DecIso12[v,(action)] decomposes an isometry v from 1 to 2 qubits acting on the qubits whose numbers are given in the list action (default: action={1,2}) into a sequence of single-qubit and C-NOT gates using a bespoke optimization for this case."
(*Optimal decomposition of an isometry*)
DecIsometry::usage="DecIsometry[v,(action)] decomposes an isometry v from m to n qubits acting on the n qubits whose numbers are given in the list action (default: action=Range[n]) into a sequence of single-qubit and C-NOT gates using the decomposition scheme that achieves the lowest known C-NOT count."
DecIsometryGeneric::usage="DecIsometryGeneric[v,(action)] decomposes an isometry v from m to n qubits acting on the n qubits whose numbers are given in the list action (default: action=Range[n]) into a sequence of single-qubit and C-NOT gates using the decomposition scheme that achieves the lowest known C-NOT count for a generic isometry of the given dimensions."
(*State preparation (Plesch and Brukner)*)
StatePreparation::usage="StatePreparation[v,(action)] prepares an n qubit state (i.e. a \!\(\*SuperscriptBox[\(2\), \(n\)]\)-dimensional normalized column vector) v on the qubits listed in action (default: action=Range[n]) into single-qubit and C-NOT gates (see the full documentation for more details)."
(*Knill's decomposition for isometries*)
KnillDec::usage="KnillDec[v, (action)] decomposes an isometry v from m to n qubits acting on the qubits whose numbers are listed in action (default: action=Range[n]) into a sequence of single-qubit and C-NOT gates based on Knill's decomposition (see the full documentation for more details)."
(*Methods for channels and POVMs*)
KrausToChoi::usage="KrausToChoi[chan] gives the Choi state corresponding to the channel chan."
ChoiToKraus::usage="ChoiToKraus[st] gives the channel corresponding to Choi state st."
MinimizeKrausRank::usage="MinimizeKrausRank[chan] gives another representation of the channel chan potentially with fewer Kraus operators."
StinespringQubit::usage="TBA."
POVMToIsometry::usage="TBA."
(*Decompose channels and POVMs in the quantum ciruit model*)
DecChannelInQCM::usage="TBA."
DecPOVMInQCM::usage="TBA."
PrepareForQASM::usage="PrepareForQASM[gatelist] takes a list of gates and prepares it into a form suitable for use with the python script that converts to QASM."
PickRandomCircuitIsometry::usage="PickRandomCircuitIsometry[m,n,t] creates a random circuit from m qubits to n qubits with t CNOTs, where an arbitrary single qubit unitary is included at the start and after each CNOT. With Option TotGates->True, the value of t is the total number of gates and these are placed randomly."
RxRGateDecomp::usage="RxRGateDecomp[u] takes a single qubit unitary u and outputs (a,b,c,d) such that u is equal to Rx[a] followed by R[b,c] up to the phase E^(I*d)."
ReplaceCNOTWithXX::usage="ReplaceCNOTWithXX[st] takes a gate list and replaces all CNOT gates with XX gates and additional single qubit rotations."
ReplaceXXWithCNOT::usage="ReplaceXXWithCNOT[st] takes a gate list and replaces all XX gates with CNOTs and additional single qubit rotations."
CNOTRotationsToXXRGates::usage="CNOTRotationsToXXRGates[st] takes a gate list and replaces all CNOT and single-qubit rotations by XX and R gates."
XXRGatesToCNOTRotations::usage="XXRGatesToCNOTRotations[st] takes a gate list and replaces all XX and R gates with CNOTs and single-qubit rotations."

Begin["`Private`"];
(*Set debug to True to run additional tests during running the methods.*)
debug=False;
(*Set analyzeAnalyticDecUnitary2Qubits to True to get out intermediat results for DecUnitary2Qubits.*)
analyzeAnalyticDecUnitary2Qubits=False;
(*Set analyzeAnalyticCCDec to True to get out intermediat results for ColumnByColumnDec.*)
analyzeAnalyticCCDec=False;
(*Set analyzeAnalyticCCDec to True to get out intermediat results for SimplifyGateList.*)
analyzeAnalyticSimplifyGateList=False;
(*Set analyzeAnalyticCCDec to True to get out intermediat results for QSD.*)
analyzeAnalyticQSD=False;
(*Switch off warnings that are generated by Null arguments (e.g., If[True,,a=5] would return a warning)*)
Off[Syntax::com]

(*----------- Methods to handle and simplify gate sequences in list representation (public)------------*)

(*If one adds support for additional gates to one of the following methods,
 one should update all of them.*)
 
(*Displays each element in a list in matrix form*)
MatrixFormOp[m_]:=If[Length[Dimensions[m]]<=2,If[m=={},{},MatrixForm[m]],MatrixForm[#]&/@m]


Options[CreatePOVMFromGateList]={DropZero->"Last"};
CreatePOVMFromGateList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{mat,i},
IsListForm[st,"CreatePOVMFromGateList"];
If[OptionValue[DropZero]==="None",mat=CreateChannelFromList[DeleteCases[st,{4,0,_}],n,{POVM->True,DropZero->False}];Map[CT[#].#&,mat],
If[OptionValue[DropZero]==="All",mat=CreateChannelFromList[DeleteCases[st,{4,0,_}],n,{POVM->True,DropZero->True}];Map[CT[#].#&,mat],
If[OptionValue[DropZero]==="Last",mat=CreatePOVMFromGateList[st,n,DropZero->"None"];
For[i=Length[mat],i>=1,i--,If[Chop[mat[[i]]]==0*mat[[i]],mat=Drop[mat,-1],Break[]]];mat,Throw[StringForm["Unknown option for DropZero in CreatePOVMFromGateList. Valid options are None, All or Last"]]]]
]]

Options[NCreatePOVMFromGateList]={DropZero->"Last"};
NCreatePOVMFromGateList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{},
IsListForm[st,"NCreatePOVMFromGateList"];
CreatePOVMFromGateList[NGateList[st],n,{DropZero->OptionValue[DropZero]}]
]

(*Relabel qubits in list numIn in the list format input st with numbers in numOut. Note that numIn must be of the same length as numOut.*)
RelabelQubits[st_,numIn_,numOut_]:=Module[{pos,st2,act},
IsListForm[st,"RelabelQubits"];
st2=st;
Do[
Which[
st2[[i]][[1]]==-2,
act=Flatten[If[MemberQ[numIn,#],numOut[[Position[numIn,#][[1]]]],#]&/@st2[[i]][[3]]];
st2=ReplacePart[st2,i->{-2,st2[[i]][[2]],act}]
,
MemberQ[{-1,0},st2[[i]][[1]]],st2=ReplacePart[st2,i->{st2[[i]][[1]],pos=Flatten[Position[numIn,st2[[i]][[2]]]];If[Length[pos]!=0,numOut[[pos]][[1]],st2[[i]][[2]]],pos=Flatten[Position[numIn,st2[[i]][[3]]]];If[Length[pos]!=0,numOut[[pos]][[1]],st2[[i]][[3]]]}],
MemberQ[{1,2,3,4,5,6},st2[[i]][[1]]],st2=ReplacePart[st2,i->{st2[[i]][[1]],st2[[i]][[2]],pos=Flatten[Position[numIn,st2[[i]][[3]]]];If[Length[pos]!=0,numOut[[pos]][[1]],st2[[i]][[3]]]}],
st2[[i]][[1]]==100,act=Flatten[If[MemberQ[numIn,#],numOut[[Position[numIn,#][[1]]]],#]&/@st2[[i]][[3]]];st2=ReplacePart[st2,i->{100,st2[[i]][[2]],act}],
st2[[i]][[1]]==101,st2=ReplacePart[st2,i->{101,st2[[i]][[2]],pos=Flatten[Position[numIn,st2[[i]][[3]]]];If[Length[pos]!=0,numOut[[pos]][[1]],st2[[i]][[3]]]}],
True,Throw[StringForm["Gate type is not supported by RelabelQubits."]]
],
{i,Length[st]}];
st2]

(*Counts the number of cnot gates in the list st*)
CNOTCount[st_]:=
Module[{i,count},(
IsListForm[st,"CNOTCount"];
count=0;
For[i=1,i<=Length[st],i++,
If[st[[i]][[1]]==0,count=count+1]
];
count
)]

(*
Reverses the elements in the list representation along with changing the signs of the angles in the rotation matrices. 
This gives the list representation for the decomposition of a matrix that is the conjugate transpose of the one represented by the input 
*)
InverseGateList[st_]:=
Module[{numstr,i,st2},(
IsListForm[st,"InverseGateList"];
st2=Reverse[st];
Do[
Which[
MemberQ[{-2},st2[[i]][[1]]],st2=ReplacePart[st2,i->{st2[[i]][[1]],Conjugate[st2[[i]][[2]]],st2[[i]][[3]]}],
MemberQ[{-1,0},st2[[i]][[1]]],,
MemberQ[{1,2,3},st2[[i]][[1]]],st2=ReplacePart[st2,i->{st2[[i]][[1]],-st2[[i]][[2]],st2[[i]][[3]]}],
st2[[i]][[1]]==4,Throw[StringForm["Measurements/Traces in InverseGateList[] are not reversable."]],
st2[[i]][[1]]==5||st2[[i]][[1]]==6,st2=ReplacePart[st2,i->{If[st2[[i]][[2]]==5,6,5],st2[[i]][[2]],st2[[i]][[3]]}],
st2[[i]][[1]]==100,st2=ReplacePart[st2,i->{100,-st2[[i]][[2]],st2[[i]][[3]]}],
st2[[i]][[1]]==101,st2=ReplacePart[st2,i->{101,{-st2[[i]][[2]][[1]],st2[[i]][[2]][[2]]},st2[[i]][[3]]}],
True,Throw[StringForm["Unknown gate `` in InverseGateList[] cannot be reversed.",st2[[i]]]]
],
{i,Length[st]}
];
st2
)]

(* Adjusts the angles of the gates in the list representation and adjust them to be between 0 and 2\[Pi] *)
AdjustAngle[st_]:=
Module[{numstr,i,st2},(
IsListForm[st,"AdjustAngle"];
If[Length[Dimensions[st]]==1,
If[st=={},
Return[{}],
Return[AdjustAngle[{st}][[1]]]
]
,
st2=st;
];
Do[
If[MemberQ[{1,2,3},st2[[i]][[1]]],st2=ReplacePart[st2,i->{st2[[i]][[1]],AdjustAngleHelp[st2[[i]][[2]]],st2[[i]][[3]]}]],
{i,Length[st]}
];
st2
)]

GateTypes[]:=Module[{},Print["Gate types for UniversalQCompiler"];
Print["{-2,diag,act}: diagonal gate with entries diag on the qubits listed in act"];
Print["{-1,n,m}: CZ with qubit n the control and m the target"];
Print["{0,n,m}: CNOT with qubit n the control and m the target"];
Print["{1,t,n}: x-rotation by angle t for qubit n"];
Print["{2,t,n}: y-rotation by angle t for qubit n"];
Print["{3,t,n}: z-rotation by angle t for qubit n"];
Print["{4,0,n}: trace out qubit n"];
Print["{4,1,n}: measure qubit n in the computational basis"];
Print["{5,0,n}: qubit n starts in state |0>"];
Print["{5,1,n}: qubit n starts in state |1>"];
Print["{6,0,n}: qubit n is postselected on |0>"];
Print["{6,1,n}: qubit n is postselected on |1>"];
Print["{100,t,{n,m}}: XX-gate with angle t on qubits n,m"];
Print["{101,{t,p},n}: R-gate with angles t,p on qubit n"]]

(*Transforms a list in list format to a list containing the corresponding matrices*)
ListFormToOp[list_,n_:Null]:=Module[{numQubits=n},
IsListForm[list,"ListFormToOp"];
If[n==Null,numQubits=If[Dimensions[list]=={3},NumberOfQubits[{list}],NumberOfQubits[list]]];If[Dimensions[list]=={3},
Which[
MemberQ[{1,2,3},list[[1]]],RotGateM[list[[2]],list[[1]],list[[3]],numQubits],
list[[1]]==0,CNOTM[list[[2]],list[[3]],numQubits],
list[[1]]==-1,CZM[list[[2]],list[[3]],numQubits],
list[[1]]==100,XXM[list[[2]],list[[3]][[1]],list[[3]][[2]],numQubits],
list[[1]]==101,RGateM[list[[2]][[1]],list[[2]][[2]],list[[3]],numQubits],
list[[1]]==-2,If[Dimensions[list[[2]]]!= {},DiagMat[list[[2]],list[[3]],numQubits],Throw[StringForm["Unspecified diagonal gate cannot be represented as a matrix using ListFormToOp[]"]]],
list[[1]]==4,Throw[StringForm["Measurements/Trace cannot be represented by matrices using ListFormToOp[]"]],
list[[1]]==5,Throw[StringForm["Starting in |0> or |1> cannot be represented by matrices using ListFormToOp[]"]],
list[[1]]==6,Throw[StringForm["Postselection cannot be represented by matrices using ListFormToOp[]"]]
],ListFormToOp[#,numQubits]&/@list]]

(*Transforms a list in list format to a list containing the corresponding string representations of the gates*)
ListFormToStr[list_]:=
Module[{},
IsListForm[list,"ListFormToStr"];
If[Length[Dimensions[list]]==1,If[list=={},{},ListFormToStrSingleGate[list]],ListFormToStrSingleGate/@list]
]

Options[SimplifyGateList]={FullSimp->True};
SimplifyGateList[st_,OptionsPattern[]]:=Module[{traceout,ancillain,ancillaout,out},
IsListForm[st,"SimplifyGateList"];If[Cases[st,{100,_,_}]==={}&&Cases[st,{101,_,_}]==={},
traceout=Cases[st,{4,_,_}];ancillain=Cases[st,{5,_,_}];ancillaout=Cases[st,{6,_,_}];
out=Join[ancillain,Reverse[SimplifyGateListReverseGateOrder[Reverse[DeleteCases[st,x_/;x[[1]]>=4]],OptionValue[FullSimp]]],ancillaout,traceout],
Throw[StringForm["Gatetypes 100 and 101 cannot be used in SimplifyGateList[]"]]];out]

NSimplifyGateList[st_]:=Module[{},
IsListForm[st,"NSimplifyGateList"];
SimplifyGateList[NGateList[st]]
]

(*----------- Methods to handle and simplify gate sequences in list representation (private)------------*)

(*Create an n qubit isometry from list form. Multiplies the unitaries described in the list (in reversed order!)  and outputs the first m columns*)
(* Use FullSimp\[Rule]False to avoid attempts to use FullSimplify *)
Options[CreateIsometryFromList]={FullSimp->True};
CreateIsometryFromList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{mat,mat2,i,k,ancillain,ancillainnums,ancillainvals,ancillaout,ancillaoutnums,ancillaoutvals,st2,id,rest,n1=n},
IsListForm[st];
If[n===Null,n1=NumberOfQubits[st]];ancillain=SortBy[Cases[st,{5,_,_}],Last];
If[ancillain==={},ancillainnums={},ancillainnums=Transpose[ancillain][[3]];ancillainvals=Transpose[ancillain][[2]]];ancillaout=SortBy[Cases[st,{6,_,_}],Last];
If[ancillaout==={},ancillaoutnums={},ancillaoutnums=Transpose[ancillaout][[3]];ancillaoutvals=Transpose[ancillaout][[2]]];
st2=DeleteCases[st,{x_/;x==5||x==6,_,_}];mat={{1}};k=0;
For[i=1,i<=n1,i++,mat=KroneckerProduct[mat,If[MemberQ[ancillainnums,i],k++;KetV[ancillainvals[[k]],2],IdentityMatrix[2]]]];
isAnalytic=True;
For[i=1,i<=Length[st2],i++,
If[isAnalyticGate[st2[[i]]],,isAnalytic=False];
If[isAnalytic,
If[OptionValue[FullSimp],mat=FullSimplifyNoRoots[ListFormToOp[st2[[i]],n1].mat],mat=Simplify[ListFormToOp[st2[[i]],n1].mat],Print["CreateIsometryFromList: Error"]],
mat=ListFormToOp[st2[[i]],n1].mat
]
];
If[ancillaoutnums=={},mat,mat2={{1}};k=0;For[i=1,i<=n1,i++,(* usually ancillaoutvals will be all 1s, so create ket 0, sometimes (for instrument generation) we want to create ket 1 *)mat2=KroneckerProduct[mat2,If[MemberQ[ancillaoutnums,i],k++;KetV[ancillaoutvals[[k]],2],IdentityMatrix[2]]]];
CT[mat2].mat]]

(*Create an n qubit isometry from list form. Multiplies the unitaries described in the list (numerically) and outputs the first m columns*)
NCreateIsometryFromList[st_,n_:Null]:=CreateIsometryFromList[NGateList[st],n]

Options[CreateChannelFromList]={POVM->False,DropZero->True,FullSimp->True};
CreateChannelFromList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{mat,i,traces,tracesnums,postsel,postselnums,posn,chanout,st2,dims,n1=n},If[Not[OptionValue[POVM]]&&MemberQ[st,{4,1,_}],Print["CreateChannelFromList: measurement gate type found"]];
If[n===Null,n1=NumberOfQubits[st]];traces=Cases[st,{4,_,_}];If[traces==={},tracesnums={},tracesnums=Transpose[traces][[3]]];postsel=Cases[st,{6,_,_}];If[postsel==={},postselnums={},postselnums=Transpose[postsel][[3]]];If[Dimensions[Intersection[tracesnums,postselnums]]=={0},,Print["CreateChannelFromList: Cannot postselect on zero and measure/trace on the same qubit."]];st2=DeleteCases[st,{4,_,_}];
mat=CreateIsometryFromList[st2,n1,FullSimp->OptionValue[FullSimp]];
chanout={};posn={};dims={};For[i=1,i<=n1,i++,If[MemberQ[postselnums,i],,If[MemberQ[tracesnums,i],posn=Insert[posn,1,-1];dims=Insert[dims,{1,2},-1],posn=Insert[posn,2,-1];dims=Insert[dims,{2,2},-1]]]];For[i=0,i<=2^(Length[tracesnums])-1,i++,chanout=Insert[chanout,Tensor[BraV[i,2^(Length[tracesnums])],IdentityMatrix[2^(n1-Length[postselnums]-Length[tracesnums])],posn,dims].mat,-1]];
If[OptionValue[DropZero],For[i=Length[chanout],i>=1,i--,If[Chop[chanout[[i]]]==0*chanout[[i]],chanout=Drop[chanout,{i}]]]];chanout]

Options[NCreateChannelFromList]={POVM->False,DropZero->True};
NCreateChannelFromList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=CreateIsometryFromList[NGateList[st],n,{POVM->OptionValue[POVM],DropZero->OptionValue[DropZero]}]

Options[CreateInstrumentFromList]={DropZero->True,FullSimp->True};(* using DropZero here prevents identification using Length[Dimensions[out]], where out is the output of CreateOperationFromGateList *)
CreateInstrumentFromList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{mat,i,j,traces,tracesnums,postsel,postselnums,posn,mmt,mmtnums,inst,chanout,st2,st3,digs,dims,n1=n},
If[n===Null,n1=NumberOfQubits[st]];traces=Cases[st,{4,0,_}];If[traces==={},tracesnums={},tracesnums=Transpose[traces][[3]]];postsel=Cases[st,{6,_,_}];If[postsel==={},postselnums={},postselnums=Transpose[postsel][[3]]];mmt=Cases[st,{4,1,_}];If[mmt==={},mmtnums={},mmtnums=Transpose[mmt][[3]]];If[Dimensions[Intersection[tracesnums,postselnums,mmtnums]]=={0},,Print["CreateInstrumentFromList: Cannot have combinations of postselect on zero/measure/trace on the same qubit."]];inst={};st2=DeleteCases[st,{4,1,_}];
For[j=1,j<=2^(Length[mmtnums]),j++,digs=IntegerDigits[j-1,2,Length[mmtnums]];st3=st2;For[i=1,i<=Length[mmtnums],i++,st3=Insert[st3,{6,digs[[i]],mmtnums[[i]]},-1]];
inst=Insert[inst,CreateChannelFromList[st3,n1,{DropZero->OptionValue[DropZero],FullSimp->OptionValue[FullSimp]}],-1]];inst]

Options[NCreateInstrumentFromList]={DropZero->True};(* using DropZero here prevents identification using Length[Dimensions[out]], where out is the output of CreateOperationFromGateList *)
NCreateInstrumentFromList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=CreateInstrumentFromList[NGateList[st],n,{DropZero->OptionValue[DropZero]}]

(*ToDo: Improve efficiency by implementing the application of C-NOTs and single-qubit rotations efficiently*)
Options[CreateOperationFromGateList]={FullSimp->True};
CreateOperationFromGateList[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{four0s,four1s},four1s=MemberQ[st,{4,1,_}];
If[four1s,CreateInstrumentFromList[st,n,{DropZero->False,FullSimp->OptionValue[FullSimp]}],four0s=MemberQ[st,{4,0,_}];
If[four0s,CreateChannelFromList[st,n,FullSimp->OptionValue[FullSimp]],CreateIsometryFromList[st,n,FullSimp->OptionValue[FullSimp]]]]]

NCreateOperationFromGateList[st_,n_: Null]:=Module[{four0s,four1s},four1s=MemberQ[st,{4,1,_}];
If[four1s,NCreateInstrumentFromList[st,n,DropZero->False],four0s=MemberQ[st,{4,0,_}];
If[four0s,NCreateChannelFromList[st,n],NCreateIsometryFromList[st,n]]]]

ListFormToStrSingleGate[list_]:=Switch[list[[1]],
-2,ToString[StringForm["D(`1`)(`2`)",list[[2]],list[[3]] ]],
-1,CZNotationStr[list[[2]],list[[3]]],
0,CNOTNotationStr[list[[2]],list[[3]]],
1,ToString[StringForm["Rx(`2`)(`1`)",list[[3]],list[[2]] ]],
2,ToString[StringForm["Ry(`2`)(`1`)",list[[3]],list[[2]] ]],
3,ToString[StringForm["Rz(`2`)(`1`)",list[[3]],list[[2]] ]],
4,If[list[[2]]==0,ToString[StringForm["Tr(`1`)",list[[3]] ]],ToString[StringForm["M(`1`)",list[[3]]]]],
5,ToString[StringForm["|`1`>(`1`)",list[[2]],list[[3]]]],
6,ToString[StringForm["<`1`|(`1`)",list[[2]],list[[3]]]],
100,ToString[StringForm["XX(`1`)(`2`,`3`))",list[[2]],list[[3,1]],list[[3,2]]]],
101,ToString[StringForm["RGate(`1`,`2`)(`3`))",list[[2,1]],list[[2,2]],list[[3]]]],
_,Throw[StringForm["Unknown gate type found in ListFormToStr[]."]]
]

(*Input: st: gate sequence containing only single-qubit and C-NOT gates, n: total number of qubits.
Merges the single-qubit gates if this helps to reduce the gate count.
Commutes Rz rotation with the control and Rx rotations with the target of C-NOT gates. Moreover, two following C-NOT gates
are canceled out.
Output: Simplified gate sequence.*)
SimplifyGateListReverseGateOrder[st_,FullSimp_:True] := Module[{numQubits, maxNum,A,stNew,controlQubit,targetQubit,mergedContr,commutedContr,mergedTar,commutedTarg},( 
numQubits=NumberOfQubits[st];
A=ConstantArray[{},numQubits]; (*List to save the collected single-qubit gates on the n qubits. The list of a qubit is updated after merging.*)
stNew={};
Do[
If[st[[i]][[1]]==0,
(*case: C-NOT gate*)
controlQubit=st[[i]][[2]];
targetQubit=st[[i]][[3]];
If[analyzeAnalyticSimplifyGateList,Print["CNOT gate with control qubit ",controlQubit," and target qubit ", targetQubit," found."]];
If[analyzeAnalyticSimplifyGateList,Print["Call MergeAndCommuteSQG with input ",{A[[controlQubit]],False,True}]];
{mergedContr,commutedContr}=MergeAndCommuteSQG[A[[controlQubit]],False,True,FullSimp];
If[analyzeAnalyticSimplifyGateList,Print["Merged gates on control qubit ",controlQubit," are given by ",mergedContr]];
If[analyzeAnalyticSimplifyGateList,Print["Commuted gates on control qubit ",controlQubit," are given by ",commutedContr]];
If[analyzeAnalyticSimplifyGateList,Print["Call MergeAndCommuteSQG with input ",{A[[targetQubit]],True,False}]];
{mergedTar,commutedTarg}=MergeAndCommuteSQG[A[[targetQubit]],True,False,FullSimp];
If[analyzeAnalyticSimplifyGateList,Print["Merged gates on target qubit ",targetQubit," are given by ",mergedTar]];
If[analyzeAnalyticSimplifyGateList,Print["Commuted gates on target qubit ",targetQubit," are given by ",commutedTarg]];
A[[controlQubit]]={};
A[[targetQubit]]={};
(*Add merged gates:*)
stNew=Join[stNew,mergedContr];
stNew=Join[stNew,mergedTar];
(*Append C-NOT:*)
stNew=AppendCNOT[stNew,st[[i]]];
(*If two C-NOTs were cancelled out, we have to add the last single-qubit gates again to the array A*)
If[Length[stNew]>=1&&stNew[[-1]][[1]]!= 0,
If[analyzeAnalyticSimplifyGateList,Print["Call ExtractSQG with input ",{stNew,controlQubit}]];
{stNew,A[[controlQubit]]}=ExtractSQG[stNew,controlQubit];
If[analyzeAnalyticSimplifyGateList,Print["Call ExtractSQG with input ",{stNew,targetQubit}]];
{stNew,A[[targetQubit]]}=ExtractSQG[stNew,targetQubit];
];
(*Upade gate array:*)
A[[controlQubit]]=Join[A[[controlQubit]],commutedContr];
A[[targetQubit]]=Join[A[[targetQubit]],commutedTarg];
,
(*case: single-qubit gate acting on st[[i]][[3]]*)
AppendTo[A[[st[[i]][[3]]]],st[[i]]];
] 
,{i,1,Length[st]}
];
(*Add the single-qubit gates at the beginning of the circuit to stNew*)
Do[
If[analyzeAnalyticSimplifyGateList,Print["Call MergeAndCommuteSQG with input ", {A[[i]],False,False}]];
stNew=Join[stNew,MergeAndCommuteSQG[A[[i]],False,False,FullSimp][[1]]]
,{i,1,numQubits}
];
If[analyzeAnalyticSimplifyGateList,Print["Call AdjustAngle with input ", stNew]];
AdjustAngle[stNew]
)]

(*Gets a sequence of single-qubit gates as an input and the information about the previous action of the two qubit gate: control or target of 
the C-NOT gate. If "isTarget" and "isControl" are wrong, then we assume that we can not commute anything to the left in the circuit.
Outputs the merged single-qubit gates as well as the gate that can be commuted trough the previous two-qubit gate.*)
MergeAndCommuteSQG[stSQGsInp_,isTarget_,isControl_,FullSimp_:True] := Module[{angles,stSQGs,mergedGates,commutedGate,gate},(
If[isTarget&&isControl,Throw[StringForm["isTarget and isControl are set to True in the method MergeAndCommuteSQG[]. This setting is impossible and not allowed as an input."]];Return[{}]];
(*Merge gates if we have two single-qubit rotations of the same kind:*)
If[analyzeAnalyticSimplifyGateList,Print["Call MergeSameRot with input ",stSQGsInp]];
stSQGs=MergeSameRot[stSQGsInp,FullSimp];
mergedGates={};
commutedGate={}; 
If[Length[stSQGs]>= 1 && Length[stSQGs]<= 2,
(*Commute single-qubit gate if possible*)
If[isTarget&&stSQGs[[-1]][[1]]==1,
commutedGate={stSQGs[[-1]]};
stSQGs=Delete[stSQGs,-1]
,
If[isControl&&stSQGs[[-1]][[1]]==3,
commutedGate={stSQGs[[-1]]};
stSQGs=Delete[stSQGs,-1];
];
];
(*Define merged gates as the ones that were not commuted*)
mergedGates=stSQGs;
];
If[Length[stSQGs]>= 3,
If[analyzeAnalyticSimplifyGateList,Print["Call MultiplySQGate with input ",stSQGs]];
gate=MultiplySQGates[stSQGs,FullSimp];
If[isTarget,
If[analyzeAnalyticSimplifyGateList,Print["Call XYXDecomposition with input ",gate]];
angles=XYXDecomposition[gate];
If[analyzeAnalyticSimplifyGateList,Print["Call of XYXDecomposition finished. "]];
angles=Reverse[angles];
If[Chop@angles[[3]]==0,
mergedGates={};
commutedGate=MergeSameRot[{{1,If[FullSimp,FullSimplifyNoRoots[Chop@angles[[4]]+Chop@angles[[2]]],Simplify[Chop@angles[[4]]+Chop@angles[[2]]]],stSQGs[[-1]][[3]]}},FullSimp];
,
mergedGates=MergeSameRot[{{1,Chop@angles[[2]],stSQGs[[-1]][[3]]},{2,Chop@angles[[3]],stSQGs[[-1]][[3]]}},FullSimp];
commutedGate=MergeSameRot[{{1,Chop@angles[[4]],stSQGs[[-1]][[3]]}},FullSimp];
]
,
If[analyzeAnalyticSimplifyGateList,Print["Call ZYZDecomposition with input ",gate]];
angles=ZYZDecomposition[gate];
If[analyzeAnalyticSimplifyGateList,Print["Call of XYXDecomposition finished. "]];
angles=Reverse[angles];
If[isControl,
If[Chop@angles[[3]]==0,
mergedGates={};
commutedGate=MergeSameRot[{{3,If[FullSimp,FullSimplifyNoRoots[Chop@angles[[4]]+Chop@angles[[2]]],Simplify[Chop@angles[[4]]+Chop@angles[[2]]]],stSQGs[[-1]][[3]]}},FullSimp]
,
mergedGates=MergeSameRot[{{3,Chop@angles[[2]],stSQGs[[-1]][[3]]},{2,Chop@angles[[3]],stSQGs[[-1]][[3]]}},FullSimp];
commutedGate=MergeSameRot[{{3,Chop@angles[[4]],stSQGs[[-1]][[3]]}},FullSimp]
]
,
(*The following line is executed in the case where there is no control or targeg previously (e.g., at the start of the circuit)*)
mergedGates=MergeSameRot[{{3,Chop@angles[[2]],stSQGs[[-1]][[3]]},{2,Chop@angles[[3]],stSQGs[[-1]][[3]]},{3,Chop@angles[[4]],stSQGs[[-1]][[3]]}},FullSimp];
];
];
];
{mergedGates,commutedGate}
)]

(*Merge the single-qubit rotations of the same kind in a list stSQGs of single-qubit rotations. Remove
rotations that rotate with angle zero.*)
MergeSameRot[stSQGs_,FullSimp_:True] := Module[{stSQGsSimplified,counter},(
If[Length[stSQGs]==0,
{}
,
If[Length[stSQGs]==1,
If[rotIsZero[stSQGs[[1]]],{},stSQGs]
,
stSQGsSimplified={};
Do[
If[Length[stSQGsSimplified]>=1,
If[stSQGsSimplified[[-1]][[1]]==stSQGs[[i]][[1]],
stSQGsSimplified[[-1]]={stSQGsSimplified[[-1]][[1]],If[FullSimp,FullSimplifyNoRoots[stSQGsSimplified[[-1]][[2]]+stSQGs[[i]][[2]]],Simplify[stSQGsSimplified[[-1]][[2]]+stSQGs[[i]][[2]]]],stSQGsSimplified[[-1]][[3]]};
If[rotIsZero[stSQGsSimplified[[-1]]],stSQGsSimplified=Delete[stSQGsSimplified,-1]]
,
If[rotIsZero[stSQGs[[i]]],,AppendTo[stSQGsSimplified,stSQGs[[i]]]]
]
,
If[rotIsZero[stSQGs[[i]]],,AppendTo[stSQGsSimplified,stSQGs[[i]]]]
]
,
{i,1,Length[stSQGs]}
];
stSQGsSimplified
]
]
)]

(*Append a C-NOT gate and cancel it if two are following each other*)
AppendCNOT[stIn_,cnot_]:=Module[{ctr,tar,st,sameContr,sameTarget},(
st=stIn;
ctr=cnot[[2]];
tar=cnot[[3]];
Do[ 
If[st[[-j]][[1]]!= 0&&(st[[-j]][[3]]==tar||st[[-j]][[3]]==ctr),AppendTo[st,cnot];Goto[End]];
If[st[[-j]][[1]]==0&&(st[[-j]][[2]]==tar||st[[-j]][[2]]==ctr||st[[-j]][[3]]==tar||st[[-j]][[3]]==ctr),
(*If we have another C-NOT with the same target qubit, but another control qubit, or another C-NOT 
with the same control qubit and anothertarget qubit  they commute 
(and we should search further for C-NOTs that may cancel.)*)
If[(st[[-j]][[2]]==ctr&&st[[-j]][[3]]!= tar)||(st[[-j]][[2]]!= ctr &&st[[-j]][[3]]== tar)
,
,
If[st[[-j]][[2]]==ctr&&st[[-j]][[3]]==tar,
st=Delete[st,{-j}];Goto[End],
(*Case wher the considered C-NOT cannot be cancelled and does not commute*)
AppendTo[st,cnot];
Goto[End]
]
]
]
,{j,1,Length[st]}];
AppendTo[st,cnot];
Label[End];
st
)];

(*Extracts the single-qubit gates at the qubit with number quNum until the next C-NOT gate*)
ExtractSQG[stIn_,quNum_]:=Module[{st,SQGs,toDelete},(
st=stIn;
SQGs={};
toDelete={};
Do[ 
If[(st[[-j]][[1]]==0&&(st[[-j]][[2]]==quNum||st[[-j]][[3]]==quNum))
||(st[[-j]][[3]]==quNum)
,
If[st[[-j]][[1]]==0,
Goto[End],
AppendTo[SQGs,st[[-j]]];
AppendTo[toDelete,{-j}]
]
]
,{j,1,Length[st]}];
Label[End];
st=Delete[st,toDelete];
{st,Reverse[SQGs]}
)];

(*---------------------------------Some methods to generate Matrices (public)---------------------------------*)

(*Rotation matrices*)
RxM[\[Alpha]_,i_:1,n_:1]:=RotGateM[\[Alpha],1,i,If[n<=i,i,n]]
RyM[\[Alpha]_,i_:1,n_:1]:=RotGateM[\[Alpha],2,i,If[n<=i,i,n]]
RzM[\[Alpha]_,i_:1,n_:1]:=RotGateM[\[Alpha],3,i,If[n<=i,i,n]]

(*Generates a CNOT gate with control on i, target on j and total number of qubits n.
Note that for inputs i and j must be labelled according to the notation -least significant qubit is labelled n and most significant qubit is labelled 1. *)
CNOTM[i_,j_,n1_:Null]:=
Module[{n,iden,gate,sup,sdown},(
n=If[n1===Null,Max[i,j],n1];
sup={1,0};
sdown={0,1};
iden=Table[IdentityMatrix[2],{n}];
gate=KroneckerProduct@@ReplacePart[iden,i->Outer[Times,sup,sup]]+KroneckerProduct@@ReplacePart[iden,{i->Outer[Times,sdown,sdown],j->PauliMatrix[1]}];
gate
)]

(*Generates a controlled z gate with control on i and  target on j, wher n is the total number of qubits (for controlled z, control and target don't matter, i.e, are interchangeable).
Note that for the input to this gate, i and j must be labelled according to the notation -least significant qubit is labelled n and most significant qubit is labelled 1. *)
CZM[i_,j_,n1_:Null]:=
Module[{n,iden,gate,sup,sdown},(
n=If[n1===Null,Max[i,j],n1];
sup={1,0};
sdown={0,1};
iden=Table[IdentityMatrix[2],{n}];
gate=KroneckerProduct@@ReplacePart[iden,i->Outer[Times,sup,sup]]+KroneckerProduct@@ReplacePart[iden,{i->Outer[Times,sdown,sdown],j->PauliMatrix[3]}];
gate
)]

XXM[phi_,i_:Null,j_:Null,n_:Null]:=Module[{out,sys,n1,i1=i,j1=j},If[i1===Null,If[j1===Null,i1=1;j1=2,If[j1==1,i1=2,i1=1]]];
If[j1===Null,If[i1==2,j1=1,j1=2]];
If[n===Null,n1=Max[i1,j1],n1=n];
out={{Cos[phi],0,0,-I*Sin[phi]},{0,Cos[phi],-I*Sin[phi],0},{0,-I*Sin[phi],Cos[phi],0},{-I*Sin[phi],0,0,Cos[phi]}};
sys=Insert[Insert[Table[1,n1-2],2,Min[i1,j1]],2,Max[i1,j1]];
Tensor[IdentityMatrix[2^n1/4],out,sys,Table[2,n1]]]

RGateM[th_,phi_,i_:1,n_:Null]:=Module[{out,sys,n1},
If[n===Null,n1=i,n1=n];
out={{Cos[th/2],-I*E^(-I*phi)*Sin[th/2]},{-I*E^(I*phi)*Sin[th/2],Cos[th/2]}};
sys=Insert[Table[1,n1-1],2,i];Tensor[IdentityMatrix[2^n1/2],out,sys,Table[2,n1]]]

(*Created the diagonal matrix with diagonal entries in diag and acting on the qubits listed in 
the list act, where we consider n qubits in total.*)
DiagMat[diag_,act1_:Null,n1_:Null]:=
Module[{n,act,diagReordered,newPos,dimDiag,dimId,pos,i,actSorted},(
n=If[n1===Null,If[act1===Null,Log2[Length[diag]],Max[act1]],n1];
act=If[act1===Null,Range[n],act1];
actSorted=Sort[act];
If[Log2[Length[diag]]>1,
newPos=Table[Position[actSorted,i][[1,1]],{i,act}];
diagReordered=Flatten[ExchangeSystems[Transpose[{diag}],newPos,ConstantArray[2,Log2[Length[diag]]]]],
diagReordered=diag
];
dimDiag=2^(Length[actSorted]);
dimId=2^(n-Length[actSorted]);
pos={};
For[i=1,i<=n,i++,
If[MemberQ[actSorted,i],
AppendTo[pos,1],
AppendTo[pos,2]
]
];
Tensor[DiagonalMatrix[diagReordered], IdentityMatrix[dimId], pos, ConstantArray[2,n]]
)]

(*For an isometry from m to n qubits (with m<n), this method adds extra columns to the matrix to make it an unitary on n qubits.*)
IsoToUnitary[m_]:=
Module[{vec,mat,dim,norm},(
dim=Dimensions[m];
mat=m\[Transpose];
mat=Join[mat,Table[ConstantArray[0,dim[[1]]],{dim[[1]]-dim[[2]]}]];
vec=Orthogonalize[NullSpace[mat]];
(*norm=Norm/@vec;
vec=vec;*)
mat[[-Length[vec];;-1]]=ConjSimplify[vec];
mat\[Transpose]
)]

(*Take a rectangular matrix with more columns than rows and adds zero rows until it is square *)
FillZero[r_] := 
 Module[{dr, dc}, {dr, dc} = Dimensions[r]; 
  If[dc - dr > 0, Join[r, ConstantArray[0, {dc - dr, dc}]], r]]
  
  
(*Extract angles from rotation matrices*)

  (*For a Rz(\[Theta]) matrix, returns the value of \[Theta]*)
RzAngle[u_]:=
Module[{th},(
th=2*Arg[u[[1,1]]];
th//Chop
)]

  (*For a Ry(\[Theta]) matrix, returns the value of \[Theta]*)
RyAngle[u_]:=
Module[{th},(
th=2*ArcTan2[u[[1,1]],u[[1,2]]];
th//Chop
)]

  (*For a Rx(\[Theta]) matrix, returns the value of \[Theta]*)
RxAngle[u_]:=
Module[{th},(
th=2*ArcTan2[u[[1,1]],Im@u[[1,2]]];
th//Chop
)]  

(*----------------------------------------Apply gates efficiently (public)---------------------------------*)

(*Apply a diagonal gate diagGate={diagonalEntries,actionQubits,bits} efficiently to matrix.*)
(*ToDo: Improve efficiency by parallelizing the computation.*)
ApplyDiag[diagGate_,mat_]:=Module[{n,diag,i,rowNum,act,diagIndex},(
rowNum=Dimensions[mat][[1]];
act=diagGate[[2]];
n=diagGate[[3]];
(*Construct diagonal*)
diag={};
For[i=1,i<= rowNum,i++,
diagIndex=FromDigits[IntegerDigits[i-1,2,n][[act]],2]+1;
AppendTo[diag,diagGate[[1]][[diagIndex]]]
];
diag*mat
)]

(*Apply a multi controlled gate MCG={gate,oneControls,zeroControls,target,bits} efficiently to a matrix mat 
(controls on zeros are also supported).*)
ApplyMCG[MCG_,mat_]:=Module[{inPos,j,con,gate,bits,target,control, ind1,ind2,mat2,index,diag,i,colNum,act,diagIndex,freeQubits,binaryInd,controlZero,controlOne},(
colNum=Dimensions[mat][[2]];
gate=MCG[[1]];
controlOne=MCG[[2]];
controlZero=MCG[[3]];
bits=MCG[[5]];
target=MCG[[4]];
control=Sort[Join[controlOne,controlZero]];
freeQubits=bits-Length[control]-1;
(*Construct a list of indexes  where the controled qubits acts on (ignoring the bit of the target qubit)*)
binaryInd=Tuples[{1, 0}, freeQubits];
For[i=1,i<=Length[control],i++,
con=control[[i]];
If[con> target,inPos=con-1,inPos=con];
If[MemberQ[controlOne,con],
binaryInd=Insert[#,1,inPos]&/@binaryInd,
binaryInd=Insert[#,0,inPos]&/@binaryInd
];
];
mat2=mat;
Do[
ind1=FromDigits[Insert[index,0,target],2]+1;
ind2=FromDigits[Insert[index,1,target],2]+1;
For[j=1,j<=colNum,j++,
{{mat2[[ind1]][[j]]},{mat2[[ind2]][[j]]}}=gate.{{mat[[ind1]][[j]]},{mat[[ind2]][[j]]}};
];
,{index,binaryInd}];
mat2
)]

(*Efficient way to apply a uniformly controlled gate UCG={gates,control,target,bits} to a matrix mat.*)
(*ToDo: Improve efficiency further by parallelizing the loop in parallel (unfortunately, ParallelDo seems not to work properly).*)
ApplyUCG[UCG_,mat_]:=Module[{newPos,orderForSQGs,controlSorted,gatesOrdered,mat2,colNum,rowNum,bits,rep,gateIndex,gates,target, control,spacing,i},(
{gates,control,target,bits}=UCG;
controlSorted=Sort[control];
If[Length[control]>1,
newPos=Table[Position[controlSorted,i][[1,1]],{i,control}];
orderForSQGs=Flatten[ExchangeSystems[Transpose[{Range[2^Length[control]]}],newPos,ConstantArray[2,Length[control]]]];
gatesOrdered=gates[[Flatten[orderForSQGs]]],
gatesOrdered=gates
];
spacing=2^(bits-target);
rep=2^(bits-Length[control]-1);
mat2=ConstantArray[0,Dimensions[mat]];(*ToDo: Improve efficiency by not copying matrix*)
Do[
i=Floor[j/spacing]*spacing+j;
gateIndex=FromDigits[IntegerDigits[i,2,bits][[control]],2]+1;
{mat2[[i+1]],mat2[[i+1+spacing]]}=gates[[gateIndex]].{mat[[i+1]],mat[[i+1+spacing]]};
,{j,0,2^(bits-1)-1}];
mat2
)]

(*----------------------------------------Visualization of cirucits (pucblic)---------------------------------*)

(*Export a gates sequence st (given in list form) to Latex using the package QCircuit.*)
Options[LatexQCircuit]={AnglePrecision->0};
LatexQCircuit[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:= Module[{draw=False},If[OptionValue[AnglePrecision]>0,draw=True];IsListForm[st,"LatexQCircuit"];qCircuitFromGridForm[computeGridForm[st, n],{DrawRotationAngles->draw,Digits->OptionValue[AnglePrecision]}]]

(*Print a gate sequence st (given in list form)*)
Options[PrintCircuit]={AnglePrecision->0};
PrintCircuit[st_,n:Except[_?OptionQ]:Null,OptionsPattern[]]:= Module[{draw=False},If[OptionValue[AnglePrecision]>0,draw=True];IsListForm[st,"PrintCircuit"];drawGraphFromGridForm[computeGridForm[st, n],{DrawRotationAngles->draw,Digits->OptionValue[AnglePrecision]}]]

(*----------------------------------------Visualization of cirucits (private)---------------------------------*)

computeGridForm[st_, n_:Null(*,maxGatesPerLine_:Infinity*)]:= Module[{xxType,postSelSetectingOps,initializingOps,gate,maxNum,numQubits,gatesOnEachWireSoFar,type, controlOrParameter, target,czType,cnotType,controlWireLength, targetWireLength,controlGate,shortestWire,longestWireLength,blankMustBeLeftBlankGate,blankWhichCanBeOverwrittenGate,minIndex, maxIndex,index,leftMostPositionThatCanBeFilled,overWriteOrAppend,
positionForNewGate,positionsOfblackMustBeLeftBlankGateOnTarget,positionsOfblackMustBeLeftBlankGateOnControl,diagType,setActingQubits,positionsOfblackMustBeLeftBlankGate,wire},
(*set some relevant constants*)
initializingOps={};
postSelSetectingOps={};
diagType=-2;
czType = -1;cnotType = 0;xxType = 100; 
controlGate = "control";
blankWhichCanBeOverwrittenGate = "blankWhichCanBeOverwritten";
blankMustBeLeftBlankGate = "mustBeLeftBlank";

leftMostPositionThatCanBeFilled[wire_] :=Module[{fp},fp=FirstPosition[wire,blankWhichCanBeOverwrittenGate,Length[wire] +1];(*We want to give back a number, but FirstPosition gives back a list in genearal, and the number one if the input is the empty set*)If[Length[fp]==0,,fp=fp[[1]]];fp];
(*Set the number of qubits if not already set as an input*)
maxNum=NumberOfQubits[st];
numQubits = Switch[n, 
Null, maxNum, 
_, n];
If[numQubits<maxNum,Print["Error in computeGridForm[]. The total number of qubits must be larger than the largest qubit numer the input circuit is acting on. Null was returned."];Return[Null]]; 
(*Initialize grid*)
gatesOnEachWireSoFar = ConstantArray[{}, numQubits];
(*Add gates:*)
Do[
{type, controlOrParameter, target} = gate;

If[type == czType || type == cnotType||type == xxType || type == diagType,
(*Code for multiqubit gates*)
If[type==diagType ||type == xxType,
{minIndex, maxIndex}={Min[target],Max[target]}
,
{minIndex, maxIndex} = If[target < controlOrParameter, {target,controlOrParameter} , {controlOrParameter, target}];
];
positionForNewGate = Max[Map[leftMostPositionThatCanBeFilled, gatesOnEachWireSoFar[[minIndex;;maxIndex]]]];
(*Fill in the grid for the entries between the control and the target (or the first and last qubit the diagonal gate is acting on)*)
For[index = minIndex+1, index < maxIndex, index++,
gatesOnEachWireSoFar[[index]] = PadRightHelp[gatesOnEachWireSoFar[[index]] , positionForNewGate, blankWhichCanBeOverwrittenGate];
gatesOnEachWireSoFar[[index, positionForNewGate]]  = blankMustBeLeftBlankGate; 
];
(*Fill up all entries where the multi-qubit gate is acting on with "blankMustBeLeftBlankGate" until the current gate in two steps:*)
If[type==diagType ||type == xxType,
setActingQubits=target,
setActingQubits={controlOrParameter,target}
];
Do[
positionsOfblackMustBeLeftBlankGate=
(*Step 1: Replace "blankWhichCanBeOverwrittenGate" entries*)Flatten[Position[gatesOnEachWireSoFar[[actQubit]],blankWhichCanBeOverwrittenGate]];
gatesOnEachWireSoFar[[actQubit]][[positionsOfblackMustBeLeftBlankGate]]= blankMustBeLeftBlankGate;
(*Step 2: Fill up with " blankMustBeLeftBlankGate" entries*)
gatesOnEachWireSoFar[[actQubit]] = PadRightHelp[gatesOnEachWireSoFar[[actQubit]] , positionForNewGate, blankMustBeLeftBlankGate];
gatesOnEachWireSoFar[[actQubit,positionForNewGate]] = {type, controlOrParameter,target};
,
{actQubit,setActingQubits}];
,
(*Initializing and post selection operations*)
If[type==5,
AppendTo[initializingOps,gate]
,
If[type==6,
AppendTo[postSelSetectingOps,gate]
,
(*Code for single qubit gates*)
positionForNewGate = leftMostPositionThatCanBeFilled[gatesOnEachWireSoFar[[target]]];
If[positionForNewGate == Length[gatesOnEachWireSoFar[[target]]] +1,
AppendTo[gatesOnEachWireSoFar[[target]], {type, controlOrParameter,target}],
gatesOnEachWireSoFar[[target]][[positionForNewGate]]= {type, controlOrParameter,target}
]
]
]
]
,{gate, st}];

longestWireLength = Max[Map[Length, gatesOnEachWireSoFar]];
For[wire = 1, wire <=  Length[gatesOnEachWireSoFar], wire ++,
gatesOnEachWireSoFar[[wire]] = PadRightHelp[gatesOnEachWireSoFar[[wire]],longestWireLength+2,"mustBeLeftBlank" ,1];
];
(*ToDo: Finish:
(*Restructure gridForm (use several rows for the representation)*)
gridLength=Length[gatesOnEachWireSoFar[[1]]];
numOfRows=Ceiling[gridLength/maxGatesPerLine];
Do[
If[i\[Equal]0,
grid=gridForm[[All,i*24+1;;(i+1)*24]],
grid=Join[grid,ConstantArray["mustBeLeftBlank",{1,24}]];
grid=Join[grid,gridForm[[All,i*24+1;;(i+1)*24]]]
],
{i,{0,1}}
]
*)
Return[{gatesOnEachWireSoFar,initializingOps,postSelSetectingOps}];
]

PadRightHelp[list_,position_,entry_,margin_:0]:=Module[{out},
If[position>=Length[list],
out=PadRight[list,position,entry,margin],
out=list;
];
out
]

Options[drawGraphFromGridForm]={DrawRotationAngles->False,Digits->2};
drawGraphFromGridForm[gridForm1_,OptionsPattern[]]:= Module[{post,init,out,postSelSet,initSet,newTraceOuts,traceOutWires,positions,gridForm,edgeRenderingFunction,actQubits,text, posTemp,longestWireLength,vertexCoordRules,edges,wire, gateIndex,sortTarget,box,vertexRenderingFunction,vertexName},
longestWireLength = Max[Map[Length, gridForm1[[1]]]];
(*Move tracing out operations to the end (Remark: Removing this leads to a "cutting" of the white lines drawn after the tracing out operations
with the control line of C-NOT gates. ToDo: Do not move tracing out to the end and fix the mentioned "bug".*)
initSet=gridForm1[[2]];
postSelSet=gridForm1[[3]];
positions=Position[gridForm1[[1]],{4,0,_}];
gridForm=ReplacePart[gridForm1[[1]],positions->"mustBeLeftBlank"];
Do[
gridForm[[pos[[1]],longestWireLength-1]]=gridForm1[[1]][[Delete[pos,0]]]
,{pos,positions}];
(*Create vertexCoordRules*)
vertexCoordRules = {};
vertexName[wire_, gateIndex_]:= ((wire -1)*longestWireLength+gateIndex-1);
For[wire =1, wire <= Length[gridForm], wire ++,
For[gateIndex = 1, gateIndex <= Length[gridForm[[wire]]], gateIndex++,
AppendTo[vertexCoordRules, vertexName[wire,gateIndex]-> {gateIndex,-wire}];
];
];
edges = {};
(*edges between adjacent gates on the same wire*)
For[wire =1, wire <= Length[gridForm], wire++,
For[gateIndex = 1, gateIndex < Length[gridForm[[wire]]], gateIndex++,
AppendTo[edges, vertexName[wire, gateIndex]-> vertexName[wire, gateIndex+1]];
];
(*edges for controls*)
For[gateIndex = 1, gateIndex <= Length[gridForm[[wire]]], gateIndex++,
If[Head[gridForm[[wire, gateIndex]]]=== List &&gridForm[[wire, gateIndex]][[1]] <= 0&&
gridForm[[wire, gateIndex]][[3]]==wire,
AppendTo[edges, vertexName[wire, gateIndex]->  vertexName[gridForm[[wire, gateIndex]][[2]],gateIndex]];
];
];
];
vertexRenderingFunction[pos_, name_] := Module[{wireIndex, positionAlongWire,type, controlOrParameter,gate,target}, 
{wireIndex,positionAlongWire} = QuotientRemainder[name, longestWireLength] + {1,1};
gate = gridForm[[wireIndex, positionAlongWire]];
(*If a qubits starts in a fixed state, draw \ket{0}*)
init={};
If[positionAlongWire==1&&MemberQ[initSet,{_,_,wireIndex}],
If[MemberQ[initSet,{_,0,wireIndex}],
init={Black, Text["|0>", pos-{0.2,0}]},
init={Black, Text["|1>", pos-{0.2,0}]}
]
];
(*If a qubits are postSelSetected in a fixed state, draw \bra{0}*)
post={};
If[positionAlongWire==longestWireLength&&MemberQ[postSelSet,{_,_,wireIndex}],
If[MemberQ[postSelSet,{_,0,wireIndex}],
post={Black, Text["<0|", pos+{0.2,0}]},
post={Black, Text["<1|", pos+{0.2,0}]}
]
];
out=Which[TrueQ[gate == "mustBeLeftBlank"],{},TrueQ[gate == "blankWhichCanBeOverwritten"],{},
{type, controlOrParameter,target} = gate;TrueQ[(type==0 ||type==-1)&& controlOrParameter==-pos[[2]]],
(*Control is recognized (from a gate type 0 or -1)*)
{Black, Disk[pos, 0.1]},
type == 1, 
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text[If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["Rx(`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"Rx"], pos]},
type == 2,
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text[If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["Ry(`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"Ry"], pos]},
type == 3,
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text[If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["Rz(`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"Rz"], pos]},
type == 101,
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text[If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["R(`1`,`2`)",ScientificForm[controlOrParameter[[1]],OptionValue[Digits]],ScientificForm[controlOrParameter[[2]],OptionValue[Digits]]]],"R"], pos]},
type == 4,
If[controlOrParameter==0,
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text["Tr", pos]},
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, Text["Mmt", pos]}
],
type == 0,
{White, EdgeForm[Black], Disk[pos , 0.2], Black, Line[{pos - {0,0.2}, pos + {0,0.2}}] , Line[{pos - {0.2,0}, pos + {0.2,0}}] },
type == -1,
{Black, Disk[pos, 0.05]},
type==-2 || type== 100,
If[Min[target]==-pos[[2]],
(*Draw the whole diagonal or XX gate*)
If[Max[target]==Min[target],
{White, EdgeForm[Black], Rectangle[pos - {0.2,0.2}, pos + {0.2,0.2}], Black, If[type == -2,Text["\[CapitalDelta]", pos],Text["XX", pos]]},
sortTarget=Sort[target];
(*Draw box for diagonal gate*)
box={White, EdgeForm[Black], Rectangle[{pos[[1]] -0.2,-sortTarget[[-1]]-0.2},{pos[[1]] +0.2,pos[[2]]+0.2}]};
(*Text for diag or XX gate*)
If[EvenQ[sortTarget[[-1]]-sortTarget[[1]]],
text={Black, Text[If[type == -2,
If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["\[CapitalDelta](`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"\[CapitalDelta]"],
If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["XX(`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"XX"]],
pos+{0,0.5-(sortTarget[[-1]]-sortTarget[[1]])/2.}]},
text={Black, Text[If[type == -2,
If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["\[CapitalDelta](`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"\[CapitalDelta]"],
If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["XX(`1`)",ScientificForm[controlOrParameter,OptionValue[Digits]]]],"XX"]],
pos+{0,-(sortTarget[[-1]]-sortTarget[[1]])/2.}]}
];
(*Do print the wire for the qubits where the diagonal gate is anot cting on*)
actQubits={};
Do[
posTemp=pos;posTemp[[2]]=-act;
AppendTo[actQubits,{Directive[Thick, Black],Line[{posTemp - {0.5,0}, posTemp + {0.5,0}}]}];
,{act,Complement[Range[Min[target],Max[target]],target]}];
Join[box,actQubits,text]
],
{}
]
,
True,
{}
];
Join[init,out,post]
];

edgeRenderingFunction[pos_,name_,a_] :=Module[{wireIndex,positionAlongWire,positionMst,out},
{wireIndex,positionAlongWire} = QuotientRemainder[name[[1]], longestWireLength] + {1,1};
(*Find position of last measurement gate*)
positionMst=Position[gridForm[[wireIndex,1;;positionAlongWire]],_?(If[Length[#]==0,False,MemberQ[{4},#[[1]]]]&),{1},Heads->False];
If[Length[positionMst]==0,
out={Directive[Thick, Black], Line[pos]},(*No measurement on the wire before the current position*)
positionMst=positionMst[[1]][[1]];
If[gridForm[[wireIndex,positionMst]][[2]]==0,
out={Directive[Thick, White], Line[pos]},(*qubit was traced out earlier*)
out={ Black, Line[pos+{{0,0.05},{0,0.05}}],Line[pos-{{0,0.05},{0,0.05}}]}(*qubit was measured out earlier*)
];
];
out
];

Return[GraphPlot[edges,VertexCoordinateRules->vertexCoordRules ,VertexRenderingFunction->vertexRenderingFunction, EdgeRenderingFunction->edgeRenderingFunction]];
];

Options[qCircuitFromGridForm]={DrawRotationAngles->False,Digits->2};
qCircuitFromGridForm[gridForm1_,OptionsPattern[]]:= Module[{gridForm,initSet,postSelSet,string,wire,gateIndex,token,type,controlOrParameter,gate,target},
gridForm=gridForm1[[1]];
initSet=gridForm1[[2]];
postSelSet=gridForm1[[3]];
string = "{\n";
For[wire = 1, wire <= Length[gridForm], wire++,
For[gateIndex = 1, gateIndex <= Length[gridForm[[wire]]], gateIndex++,
(*If a qubits starts in a fixed state, add <0/1|*)
If[gateIndex==1&&MemberQ[initSet,{_,_,wire}],
If[MemberQ[initSet,{_,0,wire}],
string=StringJoin[string, "\\lstick{\\left| 0 \\right>}"],
string=StringJoin[string, "\\lstick{\\left| 1 \\right>}"]
]
];
string = StringJoin[string, "& "];
(*If a qubits are postSelSetected in a fixed state, add <0/1|*)
If[gateIndex==Length[gridForm[[wire]]]&&MemberQ[postSelSet,{_,_,wire}],
If[MemberQ[postSelSet,{_,0,wire}],
string=StringJoin[string, "\\rstick{\\left< 0 \\right|}"],
string=StringJoin[string, "\\rstick{\\left< 1 \\right|}"]
]
];
gate = gridForm[[wire, gateIndex]] ;
token = Which[
TrueQ[gate  == "mustBeLeftBlank"], wireType[wire,gateIndex,gridForm],
TrueQ[gate  == "blankWhichCanBeOverwritten"], wireType[wire,gateIndex,gridForm],
{type, controlOrParameter,target} = gate;
type==0, If[wire==target,"\\targ",ToString[StringForm["\\ctrl{``}", target-controlOrParameter]]],
type == 1, If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["\\gate{R_x(\\textnormal{$`1`$})}",ParamToLatex[controlOrParameter,OptionValue[Digits]]]],"\\gate{R_x}"],
type ==2, If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["\\gate{R_y(\\textnormal{$`1`$})}",ParamToLatex[controlOrParameter,OptionValue[Digits]]]],"\\gate{R_y}"],
type ==3, If[OptionValue[DrawRotationAngles],ToStringStandard[StringForm["\\gate{R_z(\\textnormal{$`1`$})}",ParamToLatex[controlOrParameter,OptionValue[Digits]]]],"\\gate{R_z}"],
type == 101, If[OptionValue[DrawRotationAngles],
ToStringStandard[StringForm["\\gate{R(\\textnormal{$`1`$,$`2`$})}",ParamToLatex[controlOrParameter[[1]],OptionValue[Digits]],ParamToLatex[controlOrParameter[[2]],OptionValue[Digits]]]],
"\\gate{R}"
],
type == 4, "\\meter",
type == -1, If[wire==target,"\\control \\qw",ToString[StringForm["\\ctrl{``}", target-controlOrParameter]]],
type == -2, If[wire==Min[target],ToString[StringForm["\\multigate{``}{\\Delta}",Max[target]-Min[target]]], "\\ghost{\\Delta}"],
type == 100, If[wire==Min[target],
If[OptionValue[DrawRotationAngles],
ToStringStandard[StringForm["\\multigate{`1`}{XX(\\textnormal{$`2`$})}",Max[target]-Min[target],ParamToLatex[controlOrParameter,OptionValue[Digits]]]],
ToStringStandard[StringForm["\\multigate{``}{XX}",Max[target]-Min[target]]]
],
If[OptionValue[DrawRotationAngles],
ToStringStandard[StringForm["\\ghost{`1`}",ToStringStandard[StringForm["XX(\\textnormal{$`2`$})",Max[target]-Min[target],ParamToLatex[controlOrParameter,OptionValue[Digits]]]]]],
"\\ghost{XX}"
]
]
];
string = StringJoin[string, token, " "];
];
string = StringJoin[string, "\\\\\n"];
];
string = StringDrop[string, -3];
string = StringJoin[string, "\n}"];
Return[string];
];

(*Helper methods*)
(*Create Strings*)
ParamToLatex[param_,digits_]:=
ToString[TeXForm[ScientificForm[param,digits]]];

ToStringStandard[str_]:=ToString[str,FormatType->StandardForm];

wireType[wireIndex_,positionAlongWire_,gridForm_]:=Module[{out,positionMst},
positionMst=Position[gridForm[[wireIndex,1;;positionAlongWire]],_?(If[Length[#]==0,False,MemberQ[{4},#[[1]]]]&)];
If[Length[positionMst]==0,
out="\\qw ",
positionMst=positionMst[[1]][[1]];
If[gridForm[[wireIndex,positionMst]][[2]]==0,
out=" ",
out="\\cw "
];
];
out
]

(*----------------------------------------Matrix decompositions (public)---------------------------------*)

(*Performs the cosine-sine decomposition for an arbitrary even dimensional unitary matrix and outputs the resultant matrices.*)
(*ToDo: Improve momory requirements by giving the outputs in an efficient representation.*)
Options[CSD]={FullSimp->True};
CSD[u_,OptionsPattern[]]:=Module[{dim, u11, u12, u21, u22, v1, c, w1, v2, w2, r, l, m, dimpos, 
  posns}, 
  dim = Dimensions[u][[1]]; {{u11, u12}, {u21, u22}} = 
  Partition[u, {dim/2, dim/2}]; 
  If[analyzeAnalyticQSD,Print["Call ReorderedSVD in CSD with input u11=",u11]];
  {v1, c, w1} = ReorderedSVD[u11,OptionValue[FullSimp]];
  If[analyzeAnalyticQSD,Print["Finished running ReorderedSVD in CSD with output: ",{v1, c, w1}]];
 If[OptionValue[FullSimp],w1=CTSimplify[w1],w1=CT[w1]];
 If[analyzeAnalyticQSD,Print["Call QLDecompositionPos in CSD with input FullSimplifyNoRoots[u21.CTSimplify[w1]]=",FullSimplifyNoRoots[u21.CTSimplify[w1]]]];
 {v2,l}=QLDecompositionPos[If[OptionValue[FullSimp],FullSimplifyNoRoots[u21.CTSimplify[w1]],Simplify[u21.CT[w1]]],OptionValue[FullSimp]]; 
  If[analyzeAnalyticQSD,Print["Finished running LDecompositionPos in CSD with output: ",{v2, l}]];
 v2 = CTSimplify[v2]; 
 If[analyzeAnalyticQSD,Print["Call  RQDecompositionNeg in CSD with input FullSimplifyNoRoots[CTSimplify[v1].u12=",FullSimplifyNoRoots[CTSimplify[v1].u12]]];
 {r, w2} = 
  RQDecompositionNeg[If[OptionValue[FullSimp],FullSimplifyNoRoots[CTSimplify[v1].u12],Simplify[CT[v1].u12]],OptionValue[FullSimp]]; 
   If[analyzeAnalyticQSD,Print["Finished running RQDecompositionNeg in CSD with output: ",{r, w2}]];
 If[OptionValue[FullSimp],w2=CTSimplify[w2],w2=CT[w2]];
(* The next part makes necessary adjustments if any of the singular values are equal to1 *)
 posns = Flatten[Position[Chop[N[Diagonal[c] - 1]], 0, {1}]]; 
 dimpos = Dimensions[posns][[1]]; 
 If[dimpos > 0, If[OptionValue[FullSimp],m=Take[CTSimplify[v2].u22.CTSimplify[w2],dimpos,dimpos],m=Take[CT[v2].u22.CT[w2],dimpos,dimpos]]; 
  If[m != IdentityMatrix[dimpos], 
   If[dim/2 - dimpos != 0, 
    v2 = v2.DirectSum[m, IdentityMatrix[dim/2 - dimpos]], 
    v2 = v2.m]]]; {ArrayFlatten[{{v1, 0*v1}, {0*v1, v2}}], 
  ArrayFlatten[{{c, -l}, {l, c}}], 
  ArrayFlatten[{{w1, 0*w1}, {0*w1, w2}}]}]
  
  (* Finds the QL and RQ decompositions of square matrices, i.e., M = CT[Q].L or M = R.CT[Q] with Q unitary, L left triangular, R right triangular *)
RQDecomposition[M_]:=Module[{Q,R},{Q,R}=QRDecomp[CT[Reverse[M]]];{M.CT[Reverse[Q]],CT[Reverse[Q]]}]
QLDecomposition[M_]:=Module[{Q,R},{R,Q}=RQDecomposition[CT[M]];{CT[Q],CT[R]}]

(*
(*If ta is a tensor and the kronecker product of 2 matrices, this function finds the two matrices such that a\[CircleTimes]b=ta.
Note that this function first rearranges the matrix such that if it is a kronecker product of matrices a,b, the rearranged 
matrix is the outer product between a and b rearranged to vectors. The Singular Value Decomposition (SVD) actually gives a 
decomposition of a matrix in terms of outer product of vectors or rank 1 matrices. If our original matrix was actually a 
tensor product of 2 matrices then the rearranged matrix should be a rank 1 matrix since it would be the outer product of 
2 columns. As a result, the SVD would only have 1 term or a single non zero term in the diagonal. Since the SVD built from 
Mathematica sorts the diagonal in descending order, this corresponds to the first entry and the corresponding vectors become 
the first column of u and first column of v\[Conjugate] 
*)
KronFactor[ta_]:=
Module[{a,b,u,w,v,blockMatrix, R,vecM1,vecM2,vectorize,scaleFactor,aTrace,m1,m2,n1,n2},
(*ToDo: Should also work for general matrices sizes. Seems not to work\[Rule]fix it.*)
m1=2;
m2=2;
n1=2;
n2=2;
(*see Approximation with Kroneker products - C.F. van Loan and N. Pitsianis*)
blockMatrix = Table[ta[[(i-1)m2 +1;;i m2,(j-1)n2 +1;; j n2]], {j, 1, n1},{i,1,m1}];
vectorize = Flatten [Transpose[#]] &;
R = Map[vectorize, Flatten[blockMatrix,1]];
{u,w,v} = SingularValueDecomposition[R];
vecM1 = Sqrt[w[[1,1]]] * CT[u][[1]];
vecM2 =Sqrt[w[[1,1]]] * Transpose[v][[1]]; (*That this is the transpose, and not the conjugate transpose is deliberate, two conjugates cancel out*)
a = CT[ArrayReshape[vecM1 ,{m1,n1}]];
b =  CT[ArrayReshape[vecM2,{m2,n2}]];
(*If the input was unitary we should return unitary matrices*)
If[UnitaryMatrixQ[ta] ,
scaleFactor =Sqrt[a[[1]].Conjugate[a[[All]][[1]]] ];
a = a/scaleFactor;
b = scaleFactor * b;
];
(*If the input was Hermitian we should return Hermitian matrices*)
If[HermitianMatrixQ[ta],
aTrace = Tr[a];
If[ Chop[Im[aTrace]] != 0,
a = Abs[aTrace]/aTrace * a;
b= aTrace/Abs[aTrace] * b;
];
];
{a,b}
]
*)

(* This is a faster method than KronFactor for analytic 4x4 inputs u. 
	u must be such that u=u1\[CircleTimes]u2  *)
Options[KronFactorUnitaryDim4]={FullSimp->True};
KronFactorUnitaryDim4[u_,OptionsPattern[]]:=Module[{a,b,c,p,i,max,imax},
IsQubitIsometry[u,"KronFactorUnitaryDim4"];
p[0]=IdentityMatrix[2];p[1]={{0,1},{1,0}};p[2]={{0,I},{-I,0}};p[3]={{1,0},{0,-1}};
max=0;imax=0;
For[i=0,i<=3,i++,
a=Tr[Abs[Flatten[PartialTrace[KroneckerProduct[IdentityMatrix[2],p[i]].N[u],2,2,2]]]];
If[a>max,max=a;imax=i]];
If[max==0,Print["Error in KronFactorUnitaryDim4: the variable 'max' cannot be equal to zero."]];
(* There is a special case if u=u1\[CircleTimes]u2 with Tr[u2]=0.  In this case we rotate u2 using
a Pauli matrix. Choosing the Pauli matrix can safely be done numerically *)
If[OptionValue[FullSimp],a=FullSimplifyNoRoots[PartialTrace[KroneckerProduct[IdentityMatrix[2],p[imax]].u,2,2,2]];
c=FullSimplifyNoRoots[((a.CTSimplify[a])[[1,1]])^(1/2)];a=FullSimplifyNoRoots[a/c];
b=p[imax].SimplifyTrigo[PartialTrace[KroneckerProduct[CTSimplify[a],p[imax]].u,2,2,1]/2];,
a=Simplify[PartialTrace[KroneckerProduct[IdentityMatrix[2],p[imax]].u,2,2,2]];
c=Simplify[((a.CT[a])[[1,1]])^(1/2)];a=Simplify[a/c];
b=p[imax].Simplify[PartialTrace[KroneckerProduct[CT[a],p[imax]].u,2,2,1]/2]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix a in KronFactorUnitaryDim4: ",a];
Print["Matrix b in KronFactorUnitaryDim4: ",b]
];
{a,b}]

(*If ta is a vector which is a tensor product of 2 vectors, this function finds two vectors such that VecProd[a,b]=ta
*)
KronFactorVector[ta_]:=
Module[{a,b,c,u,w,v},(
a=ArrayReshape[ta,{2,2}];
{u,w,v}=SingularValueDecomposition[a];
b=w[[1,1]]*u[[All,1]];
c=Conjugate[v[[All,1]]];
{b,c}
)]

(*Generates the tensor product of the vectors a and b*)
VecProd[a_,b_]:=Delete[KroneckerProduct[{a},{b}],0]

(*----------------------------------------Decomposition of single-qubit gates (public)----------------------------------------*)

(*For a 2 x 2 unitary matrix, returns the angles (a,b,c,d) such that 
u=Exp[id]*Rz[c].Ry[b].Rz[a].
*)
(*While evaluating the angles, this function uses the inbuilt Arg and ArcCos functions of Mathematica. 
Arg function returns the answer between 0 and $2\pi$ and ArcCos between 0 and $\pi$. 
Therefore to take into account all possibilities this function makes a list of all possible values 
for a, b, c, d and checks which combination satisfies the requirement. Although there is no unique answer, 
the function returns the first list of values which satisfies the equation.*)
ZYZDecomposition[u0_] :=
Module[{theta, psi1, psi2,t,u}, (
u=Simplify[u0];
If[Dimensions[u][[1]] == 2,
If[Chop[N[u[[1,1]]]] == 0,
psi1 = ArgAndSimplify[-(u[[1,2]]/u[[2,1]])];
Return[SimplifyTrigo[{0,\[Pi],psi1,ArgAndSimplify[u[[1,2]]Exp[-I psi1/2]]}]];
];
If[Chop[N[u[[1,2]]]] == 0,
psi1 = ArgAndSimplify[u[[1,1]]/u[[2,2]]];
Return[SimplifyTrigo[{0,0,psi1,ArgAndSimplify[u[[1,1]]Exp[-I psi1/2]]}]];
];
t =  2 ArcCos[NormSimplify[u[[1,1]]]];
{psi1, psi2} = If[isBiggerThanZero[Sin[t/2]-Cos[t/2]],
{ArgAndSimplify[-(u[[1,1]]/u[[2,1]])],ArgAndSimplify[u[[1,1]]/u[[1,2]]]},
{-ArgAndSimplify[-(u[[2,1]]/u[[1,1]])],-ArgAndSimplify[u[[1,2]]/u[[1,1]]]}
];
Return[SimplifyTrigo[{psi2,t,psi1,ArgAndSimplify[u[[1,1]]Exp[-I (psi1 +psi2)/2]]}]]
,
u
]
)]

(*For a 2 x 2 unitary matrix, returns the angles a,b,c,d such that 
u=Exp[ia]*Rz[b].Ry[c].Rz[d]
[The method is based on ZYZDecomposition[]]*)
XYXDecomposition[u_] :=
Module[{basisTrans,uPrime}, (
basisTrans=RotGate[-Pi/2,2];
uPrime=Simplify[ConjugateTranspose[basisTrans].u.basisTrans];
ZYZDecomposition[uPrime]
)]

(*For a 2 x 2 unitary matrix acting on qubit action, returns the ZYZ decomposition in list form*)
Options[ZYZDec]={Simp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
ZYZDec[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] :=Module[{st,angles,actionQ}, (
IsQubitIsometry[u,"ZYZDec"];
angles=ZYZDecomposition[u];
Switch[action,
Null,actionQ=1,
_,actionQ=action
];
st={{3,angles[[3]],actionQ},{2,angles[[2]],actionQ},{3,angles[[1]],actionQ}};
If[OptionValue[Simp],
st=SimplifyGateListReverseGateOrder[st]
];
Reverse[st]
)
]

(*For a 2 x 2 unitary matrix acting on qubit action, returns the XYX decomposition in list form*)
Options[XYXDec]={Simp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
XYXDec[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] :=Module[{st,angles,actionQ}, (
IsQubitIsometry[u,"XYXDec"];
angles=XYXDecomposition[u];
Switch[action,
Null,actionQ=1,
_,actionQ=action
];
st={{1,angles[[3]],actionQ},{2,angles[[2]],actionQ},{1,angles[[1]],actionQ}};
If[OptionValue[Simp],
st=SimplifyGateListReverseGateOrder[st]
];
Reverse[st]
)
]

(*----------------------------------------Decomposition of uniformly controlled rotation gates (public)----------------------------------------*)

(*
Decomposes a uniformly controlled Z-rotation. The Input is given by:angles,controls,target,n. Thereby, angles is a list containing the angles of the controlled Z-rotation gates, controls is a list containing the controls, target is the qubit we are acting on and n is the total number of qubits. 
*)
DecUCZ[angles_,controls_,target_,n_]:=
Module[{st,UCZ},(
UCZ={angles,controls,target,n};
st=DecUCZYHelp[UCZ,3];
Reverse[st]
)]

(*
Decomposes a uniformly controlled Y-rotation. If TwoQubitGate\[Rule]"Cz", controlled-Z rotations (instead of C-NOT gates) are used in the decomposition. The Input is given by:angles,controls,target,n. Thereby, angles is a list containing the angles of the controlled Y-rotation gates, controls is a list containing the controls, target is the qubit we are acting on and n is the total number of qubits. 

Issue- answer may have discrepancies on the order of 10^-7 
(Dot@@op may differ from a by this order in some of the matrix entries. Chop the output with appropriate threshold to reduce to zero)
*)
Options[DecUCY] = {TwoQubitGate->  "CNOT"};
DecUCY[angles_,controls_,target_,n_,OptionsPattern[]]:=
Module[{st,UCY},(
UCY={angles,controls,target,n};
If[OptionValue[TwoQubitGate] =="CNOT",
st=DecUCZYHelp[UCY,2],
st=DecUCZYHelp[UCY,2,DecUCYWithCZ ->True]
];
Reverse[st]
)]

(*----------------------------------------Decomposition of uniformly controlled rotation gates (private)----------------------------------------*)

(*Helps to decompose a uniformly controled Z-rotation or Y-rotation gates. Only removes the control on the least significant qubit.
The input and output uses the notation UCZ={angles,controls,target,n} for a uniformly controlled Z-rotation or Y-rotation [note that dependent on the case, we interpret the list angles as rotations around the Z-axes or Y-axes on the Bloch sphere respectively] with controls in the list controls, and acting on target qubit target. Note that n denotes the total number of qubits.*)
(*Note: The gate sequence given in the output is given in reversed order (the first gate is the one which is applied at the end of the circuit)*)
Options[DecUCZYHelp1] = {DecUCYWithCZ ->  False};
DecUCZYHelp1[UCZ_,OptionsPattern[]]:=
Module[{half,d1,d2,u1,u2,rots,controls,target,n,i,out},(
{rots,controls,target,n}=UCZ;
If[Length[controls]==0,
	out=UCZ;,
half=2^(Length[controls]-1);d1={};d2={};
For[i=1,i<=half,i++,d1=AppendTo[d1,1/2*(rots[[2*i-1]]+rots[[2*i]])]];
For[i=1,i<=half,i++,d2=AppendTo[d2,1/2*(rots[[2*i-1]]-rots[[2*i]])]];
If[OptionValue[DecUCYWithCZ],
out={{d1,Drop[controls,-1],target,n},CZ[controls[[-1]],target],{d2,Drop[controls,-1],target,n},CZ[controls[[-1]],target]};
,
	out={{d1,Drop[controls,-1],target,n},CNOT[controls[[-1]],target],{d2,Drop[controls,-1],target,n},CNOT[controls[[-1]],target]};
]
];
out
)]

(*Recursively decomposes a uniformly controlled Z-rotation (set rotGate=3) or uniformly controlle Y-rotation (set rotGate=2) gate acting on the qubit with number target and control on qubits in the list controls, where the input is UCZ={angles,controls,target,n}.
This decomposition takes care of cancelling the extra CNOTs. 

[Setting the option "DecUCYWithCZ" to True, uses controlled Z-rotations instead of C-NOTs to decompose uniformly controlled Y-rotation gates.]
 *)

(*Note: The gate sequence given in the output is given in reversed order (the first gate is the one which is applied at the end of the circuit)*)

Options[DecUCZYHelp] = {DecUCYWithCZ ->  False};
DecUCZYHelp[UCZ_,rotGate_,OptionsPattern[]]:=
Module[{d1,d2,st,st1,st2,angles,controls,target,n,stTemp,CZOption},(
{angles,controls,target,n}=UCZ;
CZOption=OptionValue[DecUCYWithCZ];
If[Length[controls]==0,
st={{rotGate,angles[[1]],target}};
st,
stTemp=DecUCZYHelp1[UCZ,DecUCYWithCZ->CZOption];
If[Length[controls]==1,
st1=DecUCZYHelp[stTemp[[1]],rotGate,DecUCYWithCZ->CZOption];
st2=DecUCZYHelp[stTemp[[3]],rotGate,DecUCYWithCZ->CZOption];
st=Join[st1,{stTemp[[2]]},st2,{stTemp[[4]]}]
,
st1=DecUCZYHelp[stTemp[[1]],rotGate,DecUCYWithCZ->CZOption];
st1=Drop[st1,-1];(*Drop last C-NOT (cancels with first one of InverseGateList[st2] below*)
st2=DecUCZYHelp[stTemp[[3]],rotGate,DecUCYWithCZ->CZOption];
st2=Reverse[st2];
st2=Drop[st2,1];(*Reverse necessary since we are cancelling the two cnots from different decompositions and so the second decomposition is in reverse order to get the cnots adjacent*)
st=Join[st1,{stTemp[[2]]},st2,{stTemp[[4]]}];
]
];
st
)]

(*
Decomposes a uniformly controlled Y-rotations with C-NOTs and one controlled-Z rotation at the end (this structure is useful to have for the QSD). The Input is given by:angles,controls,target,n. 
Thereby, angles is a list containing the angles of the controlled Y-rotation gates, controls is a list containing the controls, target is the qubit we are acting on and n is the total number of qubits. 
*)
DecUCYHelpCNOTAndRz[angles_,controls_,target_,n_]:=
Module[{st,UCY,k,i},(
UCY={angles,controls,target,n};
st=DecUCZYHelp[UCY,2,DecUCYWithCZ -> True];
For[k=1,k<=Length[st]/2-1,k++,
(*Replace a Cz gate with a C-NOT*)
i=2*k;
st[[i-1]][[2]]=st[[i-1]][[2]]+\[Pi]/2;
st[[i]][[1]]=0; 
st[[i+1]][[2]]=st[[i+1]][[2]]-\[Pi]/2;
];
Reverse[st]
)]

(*-------------------------------------------Quantum Shannon Decomposition (public) --------------------------------------------*)

(*Demultiplexes a multiplexed gate with a single control on the most significant qubit. The input is given as u0,u1, where u0 
is apllied to all the qubits apart from the most significant one, if the most 
significant one is in the zero-basis state. And u1 is applied if the most significant qubit is in the one-basis state. 
Outputs the three matrices of the decomposition. 
Note that the two unitaries are gates on one qubit less than the input matrix.
Refer - Synthesis of Quantum-Logic Circuits by Vivek V. Shende et al.*)
Options[DemultiplexUC]={FullSimp->True};
DemultiplexUC[u0_,u1_,OptionsPattern[]]:=
Module[{v,d,w,block,dim},(
IsQubitIsometry[u0,"DemultiplexUC"];
IsQubitIsometry[u1,"DemultiplexUC"];
If[analyzeAnalyticQSD,Print["Call EigensystemExact in DemultiplexUC with input: ",FullSimplifyNoRoots[u0.CTSimplify[u1]]]];
{d,v}=EigensystemExact[If[OptionValue[FullSimp],FullSimplifyNoRoots[u0.CTSimplify[u1]],Simplify[u0.CT[u1]]]];
If[analyzeAnalyticQSD,Print["Finished running EigensystemExact in DemultiplexUC with output: ",{d,v}]];
v=Transpose[v]; (*output of eigensystem gives the eigenvectors as rows*)
d=Sqrt[d];
w=d*CT[v].u1;
{w,d,v}
)]

(*v is an isometry  of size 2^n x 2^m for some m,n\[Element]\[DoubleStruckCapitalN] and m\[LessEqual]n;
Decomposes am arbitrary isometry into single-qubit operators and CNOT gates using the Quantum Shannon decomposition [Refer - "Synthesis of Quantum-Logic Circuits" by Vivek V. Shende et al. and "Quantum Circuits for Isometries"].*)
Options[QSD] = {Simp->True,IgnoreAncilla->False,FullSimp->True};
QSD[v_,action_:Null,OptionsPattern[]]:=
Module[{actionQ,res,strres,stres,bits,gate,op,str,st,i,diag,diagOp,diagStr,diagSt,leastSignificantQubitsGatesIndices,twoQubitGate,twoQubitSt,stLSQI,free,dim,u},
IsQubitIsometry[v,"QSD"];
(*Extend the isometry to a unitary*)
dim=Dimensions[v];
bits=Log2[dim[[1]]];
free=Log2[dim[[1]]/dim[[2]]];
If[dim[[1]]!=dim[[2]],
If[analyzeAnalyticQSD,Print["Call IsoToUnitary on v= ",v]];
u=IsoToUnitary[v],
If[analyzeAnalyticQSD,Print["Finished running IsoToUnitary"]];
u=v;
];
If[bits==1,st=Reverse[ZYZDec[u]],
If[bits==2,
If[analyzeAnalyticQSD,Print["Call DecUnitary2Qubits on u= ",u]];
st=Reverse[DecUnitary2Qubits[u,UpToDiagonal-> False]];
If[analyzeAnalyticQSD,Print["Finished DecUnitary2Qubits with output st= ",st]];
If[debug,
If[isListEqualOpUpToPhase[Reverse[st],u],,Print["Debugging error in QSD: DecUnitary2Qubits failed on input u= ",u]];
];
,
If[analyzeAnalyticQSD,Print["Call CSDecRec with input: ",{u,Range[bits],bits,free}]];
st=CSDecRec[u,Range[bits],bits,free,OptionValue[FullSimp]]; (*The first element of st is a diagonal gate on the least significant two qubits*)
If[analyzeAnalyticQSD,Print["Finished running CSDecRec."]];
If[debug,
If[isListEqualOpUpToPhase[Reverse[st],u,free],,Print["Debugging error in QSD: CSDecRec failed on input ",{u,Range[bits],bits,free}]];
];
diag=DiagonalMatrix[st[[1]][[2]]];
st=Drop[st,1];
(*Take elements until we run into a CNOT or CZ with target or control on a gate that is not one of the two least significant qubits*)
leastSignificantQubitsGatesIndices=ExtractTwoQubitGatesIndices[st,bits-1,bits];
stLSQI=If[OptionValue[FullSimp],FullSimplifyNoRoots[RelabelQubits[st[[leastSignificantQubitsGatesIndices]],{bits-1,bits},{1,2}]],Simplify[RelabelQubits[st[[leastSignificantQubitsGatesIndices]],{bits-1,bits},{1,2}]]]; (*Relabel qubits to create matrix representation on a susystem (of all the considered qubits)*)
twoQubitGate=diag.CreateIsometryFromList[Reverse[stLSQI],2]; 
If[analyzeAnalyticQSD,Print["Call DecUnitary2Qubits in QSD with input twoQubitGate= ",twoQubitGate]];
twoQubitSt=Reverse[DecUnitary2Qubits[twoQubitGate,{1,2},{UpToDiagonal->False,Simp->OptionValue[Simp]}]];
If[analyzeAnalyticQSD,Print["Finished running DecUnitary2Qubits in QSD."]];
twoQubitSt=RelabelQubits[twoQubitSt,{1,2},{bits-1,bits}];(*Restore correct labeling of qubit w.r.t. all the considered qubits.*)
leastSignificantQubitsGatesIndices = Map[{#}&,leastSignificantQubitsGatesIndices]; (*The syntax for passing multiple indices to "Delete" requires this *);
st = Delete[st,leastSignificantQubitsGatesIndices];
st = Join[twoQubitSt, st]]];
(* Add indicator of qubits that start in 0, note that gates will be reversed later, so add to end *)
If[OptionValue[IgnoreAncilla],,For[i=1,i<=free,i++,st=Insert[st,{5,0,i},-1]]];
actionQ=
Switch[action, 
Null,Range[bits], 
_, action
];
If[actionQ==Range[bits],,st=RelabelQubits[st,Range[bits],actionQ]];
If[OptionValue[Simp],st=SimplifyGateList[Reverse[st]],st=Reverse[st]];
st];

(*-------------------------------------------Quantum Shannon Decomposition (private) --------------------------------------------*)

(*Transform matrix representation of an uniformly controled Y-rotation acting on the most significant qubit 
and controlled on all other qubits to UCY notation, i.e., UCY={angles,controls,target,n}*)
(*ToDo: Improve efficiency*)
MatrixRepToUCY[mat_]:=Module[{gateList,n,angles,i},(
n=Log2[Dimensions[mat][[1]]];
gateList=ExtractUCG[mat,Range[2,n],1,n];
angles={};
For[i=1,i<=  Length[gateList],i++,
AppendTo[angles,RyAngle[gateList[[i]]]];
];
{angles,Range[2,n],1,n}
)]

(*Transform a diagonal matrix representation (note that the input dia is a list containing the diagonal entries) of a uniformly controled Z-rotation acting on the most significant qubit 
and controlled on all other qubits to UCZ notation, i.e., UCZ={angles,controls,target,n}*)
DiagRepToUCZ[dia_]:=Module[{gateList,n,angles,len,arg,rot},(
len=Length[dia];
arg=Arg[dia];
rot=Table[(arg[[i]]-arg[[len/2+i]])/2,{i,1,len/2}];
{2*rot,Range[2,Log2[len]],1,Log2[len]}
)]

(*This function recursively applies CSD and QSDec to decompose an arbitrary unitary into single-qubit rotations and CNOT gates. The input a is the unitary matrix of size 2^m, where m=Length[act]. The input act is a list containing the numbers of the qubits we are acting on (during the recursion this list is reduced in each step). 
The recursion stops at two-qubit gates which are decomposed by DecUnitary2Qubits upto diagonals. The diagonals are merged with the gate in front of it before its decomposition. 

The first gate in the list ouput is a diagonal gate, represented as {-2,diag,act}, where diag is a list containing the diagonal entries and act is a list containing the qubit numbers the diagonal gate is acting on. *)
CSDecRec[a_,act_,n_,free_,FullSim_:True]:=
Module[{UCYDebug,st,st1,st2,st3,dim,partitioned,a1,a2,st3Debug,b1,b2,m1,m2,m3,m4,diag,UCY,block,u1,u2,m3Merged,dim2,stLast,opLast,op3First,nM3,actSorted,m1Merged},(
dim=Dimensions[a][[1]];
actSorted=Sort[act];(*Order the action qubits for simplicity*)
If[dim==4,
If[analyzeAnalyticQSD,Print["Call DecUnitary2Qubits in CSDecRec with input a= ",a]];
st=DecUnitary2Qubits[a,{1,2}, UpToDiagonal-> True];
If[analyzeAnalyticQSD,Print["Finished running DecUnitary2Qubits in CSDecRec."]];
diag=st[[-1]][[2]];
st=Drop[st,-1];
st=Reverse[RelabelQubits[st,{1,2},{n-1,n}]];
Prepend[st,{-2,diag,{actSorted[[-2]],actSorted[[-1]]}}]
,
If[analyzeAnalyticQSD,Print["Call CSD in CSDecRec with input a= ",a]];
{m1,m2,m3}=CSD[a,FullSimp->FullSim];
If[analyzeAnalyticQSD,Print["Finished running CSD in CSDecRec ."]];
If[debug,
If[isZeroMatrix[m1.m2.m3-a],,Print["Debugging error in CSDecRec: CSD failed on input u= ",a]];
If[isIdentity[m1.ConjugateTranspose[m1]]&&isIdentity[m2.ConjugateTranspose[m2]]&&isIdentity[m3.ConjugateTranspose[m3]],,Print["Error in CSDecRec: CSD outputs non unitary matrices for the input u= ",a]]
];
dim2=Dimensions[m3][[1]];
If[analyzeAnalyticQSD,Print["Call MatrixRepToUCY in CSDecRec with input m2= ",m2]];
UCY=MatrixRepToUCY[m2];
If[analyzeAnalyticQSD,Print["Finished running MatrixRepToUCY in CSDecRec."]];
If[debug,
UCYDebug=ReplacePart[UCY,1->(RyM[#1,1,1]&/@UCY[[1]])];
If[isZeroMatrix[CreateUCG@@UCYDebug-m2],,Print["Debugging error in CSDecRec: MatrixRepToUCY failed on input ",m2]];
];
If[analyzeAnalyticQSD,Print["Call DecUCYHelpCNOTAndRz in CSDecRec with input: ",{UCY[[1]],actSorted[[2;;]],actSorted[[1]],n}]];
st2=DecUCYHelpCNOTAndRz[UCY[[1]],actSorted[[2;;]],actSorted[[1]],n];
If[analyzeAnalyticQSD,Print["Finished running DecUCYHelpCNOTAndRz in CSDecRec."]];
st2=Reverse[st2];
(*ToDo: Improve efficiency of the following code. We just need to get u1,u2.*)
stLast=CZ[actSorted[[-1]],actSorted[[1]]];(*Operator corresponding to st2[[-1]]*)
stLast=RelabelQubits[{stLast},{actSorted[[-1]],actSorted[[1]]},{actSorted[[-1]]-actSorted[[1]]+1,1}][[1]];(*Relabiling qubits to prepare stLast for converison to a matrix*)
opLast=ListFormToOp[stLast,Log2[dim]];
m3Merged=opLast.m3; (*Merge the RZ gate into m3*)
st2=Drop[st2,-1];(*Ignore last entry in st2, since it was merged into m3*)
block=Partition[m3Merged,{dim2/2,dim2/2}];
u1=block[[1,1]];
u2=block[[2,2]];
If[free==0,
If[analyzeAnalyticQSD,Print["Call QSDecRec in CSDecRec with input: ",{{u1,u2},actSorted,n,free}]];
st3=QSDecRec[{u1,u2},actSorted,n,free,FullSim];
If[analyzeAnalyticQSD,Print["Finished running QSDecRec in CSDecRec."]];
If[debug,
st3Debug=Reverse[RelabelQubits[st3,actSorted,Range[Length[actSorted]]]];
If[isListEqualOpUpToPhase[st3Debug,DirectSum[u1,u2],free],,Print["Error in CSDecRec: QSDecRec failed on input ",{u1,u2}]];
],
If[FullSim,u1=FullSimplifyNoRoots[u1],u1=Simplify[u1]];
If[analyzeAnalyticQSD,Print["Call CSDecRec in CSDecRec with input: ",{u1,Drop[actSorted,1],n,If[free==0,0,free-1]}]];
st3=CSDecRec[u1,Drop[actSorted,1],n,If[free==0,0,free-1],FullSim];
If[analyzeAnalyticQSD,Print["Finished running SDecRec in CSDecRec."]];
];
nM3=Log2[Dimensions[m3][[1]]];
op3First=DiagMat[st3[[1]][[2]],{nM3-1,nM3},nM3]; (*extract diagonal gate from m3*)
m1Merged=m1.op3First;(*ToDo: Improve efficiency of multiplication with diagonal matrix*)
st3=Drop[st3,1];
(*As above, improve efficiency in the following preparation to get u1,u2.*)
block=Partition[m1Merged,{dim2/2,dim2/2}];
If[FullSim,u1=FullSimplifyNoRoots[block[[1,1]]];
u2=FullSimplifyNoRoots[block[[2,2]]],u1=Simplify[block[[1,1]]];
u2=Simplify[block[[2,2]]]];
If[analyzeAnalyticQSD,Print["Call QSDecRec in CSDecRec with input: ",{{u1,u2},actSorted,n,0}]];
st1=QSDecRec[{u1,u2},actSorted,n,0,FullSim];
If[analyzeAnalyticQSD,Print["Finished running QSDecRec in CSDecRec."]];
st=Join[st1,st2,st3];
st
]
)]

(*The input gates={u0,u1} containes the uniformly controlled unitaries of a single-controlled gate, where we control on the most significant qubit in the list act. The input n is the total number of qubits.
Note that after each CSDecRec invocation, it merges the diagonal gate at front to the gate before it before sending it to CSDecRec 
Note that the input unitaries u0,u1 must have dimesion\[GreaterEqual]4*)
QSDecRec[gates_,act_,n_,free_,FullSim_:True]:=
Module[{v,d,w,dim,st,st1,st2,st3,UCZ,nV,op3First,fr},(
If[analyzeAnalyticQSD,Print["Call DemultiplexUC in QSDecRec with input: ",{gates[[1]],gates[[2]]}]];
{w,d,v}=DemultiplexUC[gates[[1]],gates[[2]],FullSimp->FullSim];
If[analyzeAnalyticQSD,Print["Finished running DemultiplexUC in QSDecRec."]];
If[analyzeAnalyticQSD,Print["Call CSDecRec in QSDecRec with input: ",{w,Drop[act,1],n,If[free==0,0,free-1]}]];
st3=CSDecRec[w,Drop[act,1],n,If[free==0,0,free-1],FullSim];
If[analyzeAnalyticQSD,Print["Finished running CSDecRec in QSDecRec."]];
nV=Log2[Dimensions[v][[1]]];
op3First=DiagMat[st3[[1]][[2]],{nV-1,nV},nV];
v=v.op3First;(*merge diagonal to operator v*)
st3=Drop[st3,1];
If[analyzeAnalyticQSD,Print["Call CSDecRec in QSDecRec with input: ",{v,Drop[act,1],n,0}]];
st1=CSDecRec[v,Drop[act,1],n,0,FullSim];
If[analyzeAnalyticQSD,Print["Finished running CSDecRec in QSDecRec."]];
(*ToDo: Improve efficiency. Do not construct the full diagonal matrix.*)
d=Join[d,Conjugate[d]];
UCZ=DiagRepToUCZ[d];
st2=DecUCZ[UCZ[[1]],act[[2;;]],act[[1]],n];
st=Join[st1,Reverse[st2],st3];
st
)]

(*Helper method: gives all the indices of the gates in st on the two qubits with number i and j, until the qubits interact with another qubit.*)
ExtractTwoQubitGatesIndices[st_,i_,j_]:=
Module[{notInteractedi,notInteractedj,n,ind,controlAndActionQ},(
notInteractedi=True;
notInteractedj=True;
n=1;
ind={};
While[n<=Length[st]&&(notInteractedi||notInteractedj),
If[st[[n]][[1]]==0||st[[n]][[1]]==-1,
controlAndActionQ={st[[n]][[2]],st[[n]][[3]]};
If[MemberQ[controlAndActionQ,i],
If[MemberQ[controlAndActionQ,j],(*C-NOT is acting on qubit i and j*)
AppendTo[ind,n]
,
(*C-NOT is acting on qubit i but not on j*)
notInteractedi=False
]
,
(*C-NOT is not acting on qubit i*)
If[MemberQ[controlAndActionQ,j],
notInteractedj=False;(*C-NOT is not acting on qubit i, but on j*)
,
(*C-NOT is not acting on qubit i and not on j*)
]
],
If[(st[[n]][[3]]==i &&notInteractedi)||(st[[n]][[3]]==j &&notInteractedj),AppendTo[ind,n]]
];
n++;
];
ind
)]

(*-------------------------------------------Decomposition of diagonal gates (public) --------------------------------------------*)

(* 
d\[Rule]list containing 2^n diagonal entries of a unitary on n qubit, for some natural number n.
Decomposes a diagonal gate into single qubit gates, Rz and CNOT gates.
*)
DecDiagGate[dia_,action_:Null]:=
Module[{n,actionQ,st},(
n=Log2[Length[dia]];
actionQ=
Switch[action, 
Null,Range[n], 
_, action
];
st=DecDiagGateRec[dia,Range[n],n];
If[actionQ==Range[n],,st=RelabelQubits[st,Range[n],actionQ]];
Reverse[st]
)]

(*-------------------------------------------Decomposition of diagonal gates (private) --------------------------------------------*)

(* 
d\[Rule]list containing 2^m diagonal entries (with m\[GreaterEqual] 1) of a unitary on n qubit, where m=Length[act] and act is a list containing the numbers of the m qubits the diagonal gate is acting on.
Decomposes a diagonal gate (up to gloabal phase) into single qubit gates, Rz and CNOT gates.
*)
DecDiagGateRec[dia_,act_,n_]:=
Module[{len,arg,newdia,rot,st,st1,st2,Zangle},(
If[Length@dia==2,
Zangle=ZYZDecomposition[DiagonalMatrix[ dia]][[3]];
st={{3,Zangle,act[[-1]]}};
,
len=Length[dia]/2;
arg=Arg[dia];
newdia=Table[(arg[[i]]+arg[[len+i]])/2,{i,1,len}];
newdia=Exp[I newdia];
rot=Table[(arg[[i]]-arg[[len+i]])/2,{i,1,len}];
st1=DecUCZ[2*rot,Range[act[[2]],act[[-1]]],act[[1]],n];
st2=DecDiagGateRec[newdia,Drop[act,1],n];
st=Join[Reverse[st1],st2];
];
st
)]

(*-------------------------------------------Decomposition of two qubit gates (public) --------------------------------------------*)

(*
Given a 4 x 4 unitary matrix acting on the qubits in the list action, this function decomposes it into CNOTs and single-qubit rotations up to a diagonal matrix.
The method returns {st,d}, where st is a list of gates in the list form and d is a list containing the diagonal entries of the diagonal gate d.
Refer- http://arxiv.org/abs/quant-ph/0308045v3 (Small circuit structure in two qubit operators)
*)
Options[DecUnitary2Qubits]={UpToDiagonal->False,precision->10^-10,Simp->True,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
DecUnitary2Qubits[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]]:= Module[{actionQ,diag,st,qBits},
IsQubitIsometry[u,"DecUnitary2Qubits"];
qBits=2;
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[OptionValue[UpToDiagonal]==True,
{diag,st}=DecUnitary2QubitsHelp[SimplifyTrigo[u],2,{UpToDiagonal->OptionValue[UpToDiagonal],precision->OptionValue[precision],FullSimp->OptionValue[FullSimp]}];
st=RelabelQubits[st,{1,2},{2,1}];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
Append[If[OptionValue[Simp],SimplifyGateList[Reverse[st],FullSimp->OptionValue[FullSimp]],Reverse[st]],{-2,diag,actionQ}]
,
st=DecUnitary2QubitsHelp[SimplifyTrigo[u],2,{UpToDiagonal->OptionValue[UpToDiagonal],precision->OptionValue[precision],FullSimp->OptionValue[FullSimp]}];
st=RelabelQubits[st,{1,2},{2,1}];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
If[OptionValue[Simp],SimplifyGateList[Reverse[st],FullSimp->OptionValue[FullSimp]],Reverse[st]]
]]

(*-------------------------------------------Decomposition of two qubit gates (private) --------------------------------------------*)
(*ToDo: Make some of the following methods public and tidy up the code (e.g., never use the notation that the first 
qubit is the least significant one.*)

(*
Given two commuting, diagonalisable matrices u and v, SimultaneouslyDiagonalize returns the unitary that simultaneously diagonalizes them.
{valA, valB, vecs} = SimultaneouslyDiagonalize[A,B] => ,
A = CT[vecs].DiagonalMatrix[valA].vecs,
B = CT[vecs].DiagonalMatrix[valB].vecs

JointEigensystem returns the same in a slightly different form:
{valA, valB, vecs} = JointEigensystem[A,B] => ,
A = Transpose[vecs].DiagonalMatrix[valA].Conjugate[vecs],
B = Transpose[vecs].DiagonalMatrix[valB].Conjugate[vecs]
The precision option is used to decide how close two eigenvalues need to be to be considered equal for the purposes of the algorithm so that machine precision doesn't affect anything
*)
Options[SimultaneouslyDiagonalize]={FullSimp->True};
SimultaneouslyDiagonalize[A_,B_,precision:Except[_?OptionQ]:Null,OptionsPattern[]]:=
(*Work with high precision*)
Block[{$MaxExtraPrecision=500},
Module[{vals, vecs, unitary, blocksizes, startrow, out, i, Bblock,
  vals2, vec2, unitary2, A1, B1, mA, nA, mB, nB, valA, valB,
  vec,ord},{mA, nA} = Dimensions[A]; {mB, nB} = Dimensions[B];
 If[A.B - B.A != ConstantArray[0, {mA, nA}],A1=Chop[A];B1=Chop[B],A1=A;B1=B,A1=Simplify[A];B1=Simplify[B];]; {vals,unitary}=EigensystemExact[A1,precision];
 blocksizes =
  Transpose[
    Tally[vals,
     If[NumericQ[precision], Chop[N[#1 - #2], precision] == 0 &,
      Chop[N[#1 - #2]] == 0 &]]][[2]]; startrow = 1; out = {};
 For[i = 1, i <= Dimensions[blocksizes][[1]], i++,
  Bblock = Chop[Take[unitary.B1.If[OptionValue[FullSimp],CTSimplify[unitary],CT[unitary]], {startrow,
     startrow - 1 + blocksizes[[i]]}, {startrow,
     startrow - 1 + blocksizes[[i]]}]];
  startrow = startrow + blocksizes[[i]];
  If[Bblock == DiagonalMatrix[Diagonal[Bblock]],
 unitary2 = IdentityMatrix[Dimensions[Bblock][[1]]],
 {vals2,unitary2}=If[OptionValue[FullSimp],FullSimplifyNoRoots[EigensystemExact[Bblock,precision]],Simplify[EigensystemExact[Bblock,precision]]]];
  If[i == 1, out = unitary2, 
  out = DirectSum[out, unitary2]]];out = out.unitary;
  If[OptionValue[FullSimp],
 {valA,valB,vec}={FullSimplifyNoRoots[Diagonal[out.A1.CTSimplify[out]]],FullSimplifyNoRoots[Diagonal[out.B1.CTSimplify[out]]],FullSimplifyNoRoots[out]},
 {valA,valB,vec}={Simplify[Diagonal[out.A1.CT[out]]],Simplify[Diagonal[out.B1.CT[out]]],Simplify[out]}];(* The next If is not needed,
 but can be useful to flag problems *)
(* If[Chop[N[CT[vec].DiagonalMatrix[valA].vec - A1]] !=
    ConstantArray[0, {mA, nA}] ||
   Chop[N[CT[vec].DiagonalMatrix[valB].vec - B1]] !=
      ConstantArray[0, {mA, nA}],
  Print["Error in SimultaneouslyDiagonalize with inputs ", A, ", ",
   B, ", ", precision]]*); {valA, valB, vec}]
   ]

JointEigensystem[A_,B_,precision_:Null,FullSmp_:True] := 
 Module[{out1, out2, out3},{out1, out2, out3} = 
   SimultaneouslyDiagonalize[A, B, precision,FullSimp->FullSmp]; {out1, out2, 
   ConjSimplify[out3]}]  
   
     
   RealTrace3[u_]:=Module[{a,b,aAbs,bAbs,thetaPrime,phiPrime,theta,phi,sigma,delta,alpha,beta},
a =  u[[1,4]]u[[4,1]] - u[[1,3]]u[[4,2]] - u[[1,2]]u[[4,3]] + u[[1,1]]u[[4,4]];
b = -u[[2,4]]u[[3,1]] + u[[2,3]]u[[3,2]] + u[[2,2]]u[[3,3]] - u[[2,1]]u[[3,4]];
aAbs = Abs[a];
bAbs = Abs[b];
alpha = Arg[a];
beta = Arg[b];
delta = Arg[Det[u]];
theta = -alpha -ArcTan2[aAbs - bAbs Cos[delta - beta-alpha], -bAbs Sin[delta - beta-alpha]];
phi = -delta - theta;
DiagonalMatrix[{Exp[I theta/2],Exp[I phi/2],Exp[I phi/2],Exp[I theta/2]}]
];

(*
The Mathematica documentation does not gurantee any ordering for complex numbers which have equal real parts and imaginary parts equal in magnitude (e.g. -I and I)
This function first compares by real part then by imaginary part. So 
Re(x) > Re(y) => x \[Succeeds] y
Re(x) \[Equal] Re(y) and Im(x) \[Succeeds] Im(y) => x > y
*)
(*ComplexOrderingFunction[x_,y_]:= Which[Re[x] > Re[y], -1, Re[x] < Re[y], 1, Im[x] > Im[y], -1, Im[x] < Im[y], 1,True,0];*)
ComplexOrderingFunction[x_,y_]:= Which[N[Arg[x]] > N[Arg[y]], -1, N[Arg[x]] < N[Arg[y]], 1, N[Abs[x]] > N[Abs[y]], -1, N[Abs[x]] <N[Abs[y]], 1,True,0];

(*Ref-Minimal Universal Two-Qubit CNOT-based Circuits
Gamma function as defined definition IV.1 in the reference for 4 x 4 matrices
*)
Gam[a_]:=
Module[{i,t},(
t=KroneckerProduct[PauliMatrix[2],PauliMatrix[2]];
(*For[i=2,i<n,i++,t=KroneckerProduct[t,PauliMatrix[2]]];*)
a.t.a\[Transpose].t
)]

(*
Given u,v 4x4 unitary matrices such that there exist 2x2 unitary matrices a,b,c,d and a complex number \[Alpha] with \[LeftBracketingBar]\[Alpha]\[LeftBracketingBar]=1 such that
u = \[Alpha] a\[CircleTimes]b v c\[CircleTimes]d
returns {a,b,c,d} as 2x2 unitary matrices
*)
TwoQubitFindMatchingProductMatrices[u_,v_,precision_:10^-10,FullSim_:True] := Module[{su,sv,sv2,m,counter,u2,u3,v2,v3,uRealVal,uImagVal,uvec,uval,vRealVal,vImagVal,vvec,vval,vvec2,uvec2,uval2,vval2,
i,j,mat1,mat2,uOrder,vOrder,a,b,c,d},
If[FullSim,su=FullSimplifyNoRoots[u/(Det[u]^(1/4))];sv=FullSimplifyNoRoots[v/(Det[v]^(1/4))],su=Simplify[u/(Det[u]^(1/4))];sv=Simplify[v/(Det[v]^(1/4))]];
counter=0;
While[
(*Norm[Sort[Chop[Eigenvalues[Simplify[Gam[su]]]],ComplexOrderingFunction] - Sort[Chop[Eigenvalues[Gam[sv]]],ComplexOrderingFunction]]  > precision && counter<4,*)
Norm[Chop[EigenvaluesExact[Simplify[Gam[su]]]-EigenvaluesExact[Gam[sv]]]]>precision && counter<4,
sv = sv * I;
counter = counter+1;
];
m={{1,0,0,I},{0,I,1,0},{0,I,-1,0},{1,0,0,-I}}/Sqrt[2];
u2=Simplify[CT[m].su.m];
If[FullSim,u3=FullSimplifyNoRoots[u2.u2\[Transpose]],u3=Simplify[u2.u2\[Transpose]]];
v2=Simplify[CT[m].sv.m];
If[FullSim,v3=FullSimplifyNoRoots[v2.v2\[Transpose]],v3=Simplify[v2.v2\[Transpose]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix u2 in TwoQubitFindMatchingProductMatrices: ",u2];
Print["Matrix u3 in TwoQubitFindMatchingProductMatrices: ",u3];
Print["Matrix v2 in TwoQubitFindMatchingProductMatrices: ",v2];
Print["Matrix v3 in TwoQubitFindMatchingProductMatrices: ",v3];
];
(*
For gates we can simulate with one cnot it is the case that Re[u3] and Re[v3] end up being zero
Ignoring Re[u3] and Re[v3] entirely in this case wastes less precision and is more efficient
*)
If[AllTrue[Flatten[u3], Abs[Re[#]]< precision &],
uRealVal = {0,0,0,0};
{uImagVal, uvec} = If[FullSim,FullSimplifyNoRoots[EigensystemExact[Simplify[Im[u3]]]],Simplify[EigensystemExact[Simplify[Im[u3]]]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Eigensystem output {uImagVal, uvec}  in TwoQubitFindMatchingProductMatrices: ",{uImagVal, uvec} ];
];
,
{uRealVal,uImagVal,uvec}=If[FullSim,FullSimplifyNoRoots[JointEigensystem[Simplify[Re[u3]],Simplify[Im[u3]],precision,FullSim]],Simplify[JointEigensystem[Simplify[Re[u3]], Simplify[Im[u3]],precision,FullSim]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Eigensystem output {uRealVal,uImagVal,uvec}  in TwoQubitFindMatchingProductMatrices: ",{uRealVal,uImagVal,uvec} ];
];
];
uval = uRealVal + I uImagVal;
If[FullSim,uvec=FullSimplifyNoRoots[Map[Normalize, uvec]],uvec=Simplify[Map[Normalize, uvec]]];
If[AllTrue[Flatten[v3], Abs[Re[#]]< precision &],
vRealVal = {0,0,0,0};
{vImagVal, vvec} = If[FullSim,FullSimplifyNoRoots[EigensystemExact[Simplify[Im[v3]]]],Simplify[EigensystemExact[Simplify[Im[v3]]]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Eigensystem output {vImagVal, vvec}  in TwoQubitFindMatchingProductMatrices: ",{vImagVal, vvec} ];
];
,
{vRealVal,vImagVal,vvec}=If[FullSim,FullSimplifyNoRoots[JointEigensystem[Simplify[Re[v3]],Simplify[Im[v3]],precision,FullSim]],Simplify[JointEigensystem[Simplify[Re[v3]],Simplify[Im[v3]],precision,FullSim]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Eigensystem output {vRealVal,vImagVal,vvec}  in TwoQubitFindMatchingProductMatrices: ",{vRealVal,vImagVal,vvec} ];
];
];
vval = vRealVal + I vImagVal;
vvec = If[FullSim,FullSimplifyNoRoots[Map[Normalize, vvec]],Simplify[Map[Normalize, vvec]]];

uOrder = Ordering[Chop[uval,precision],All,ComplexOrderingFunction];
vOrder = Ordering[Chop[vval,precision],All,ComplexOrderingFunction];
uvec = uvec[[uOrder]];
vvec = vvec[[vOrder]];

If[Det[vvec] < 0,
vvec[[1]]=-vvec[[1]];
];
If[Det[uvec] < 0 ,
uvec[[1]]=-uvec[[1]];
];
mat1=If[FullSim,FullSimplifyNoRoots[m.Transpose[uvec].vvec.CTSimplify[m]],Simplify[m.Transpose[uvec].vvec.CT[m]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix mat1 in TwoQubitFindMatchingProductMatrices: ",mat1];
];
mat2=If[FullSim,FullSimplifyNoRoots[m.CTSimplify[v2].Transpose[vvec].Simplify[uvec.u2.CTSimplify[m]]],Simplify[m.CT[v2].Transpose[vvec].Simplify[uvec.u2.CT[m]]]];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix mat2 in TwoQubitFindMatchingProductMatrices: ",mat2];
];
{a,b} = KronFactorUnitaryDim4[mat1,FullSimp->FullSim];
{c,d} = KronFactorUnitaryDim4[mat2,FullSimp->FullSim];
{a,b,c,d}
];

(*
Given a 4 x 4 unitary matrix, this function decomposes it into CNOTs and single-qubit rotations up to a diagonal matrix
u\[Rule] 4x4 unitary matrix 
n\[Rule] total number of qubits (the output is given on n qubits)
returns {d, ops}, where ops is a list of operators in the "list form" format and d is a list such that DiagonalMatrix[d].(Dot @@ Map[ListFormToOp[#,2]&, ops]) \[Equal] \[Alpha] u for some phase \[Alpha]
Refer- http://arxiv.org/abs/quant-ph/0308045v3 (Small circuit structure in two qubit operators)
[Notation for DecUnitary2QubitsHelp is such that the most significant qubit in the two qubit gate is labelled n while the least significant is labelled n-1]*)
Options[DecUnitary2QubitsHelp] = { UpToDiagonal ->  False, precision ->  10^-10,FullSimp->True};
DecUnitary2QubitsHelp[u_, n_, OptionsPattern[]]:= Module[{m,det,su,diag,gu,realTraceU,traceFixingMatrix,eigenVals,st,threeCnotsInOutput},
m={{1,0,0,I},{0,I,1,0},{0,I,-1,0},{1,0,0,-I}}/Sqrt[2];
det=Simplify[Det[u]];
su=Simplify[u*ConjSimplify[det^(1/4)]];
gu=Simplify[Gam[su]];
eigenVals =  If[OptionValue[FullSimp],FullSimplifyNoRoots[Eigenvalues[gu]],Simplify[Eigenvalues[gu]]];
(*
by construction su is an su(4) matrix 
we now check if su requires 0, 1 or 2 cnots to decompose
*)
(*Prop III.1 - no cnots needed if \[Chi][\[Gamma][u]](x) = (x+1)^4 or (x-1)^4*)
If[ Fold[And, Map[Abs[#-1] < OptionValue[precision] &,N[eigenVals]]] ||  Fold[And,Map[Abs[#+1] < OptionValue[precision]&,N[eigenVals]]], 
Module[{m1, m2,a1,b1,c1,d1,a2,b2,c2,d2}, 
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix su in DecUnitary2QubitsHelp (case: zero C-NOTs): ",su]
];
{m1, m2} = KronFactorUnitaryDim4[su,FullSimp->OptionValue[FullSimp]];
{d1,c1,b1,a1} = ZYZDecomposition[m1];
{d2,c2,b2,a2}= ZYZDecomposition[m2];
st = {
{3,b1,n},
{2,c1,n},
{3,d1,n},
{3,b2,n-1},
{2,c2,n-1},
{3,d2,n-1}
};

If[OptionValue[UpToDiagonal],
Return[{{1,1,1,1}, st}],
Return[st]
];
];
];

(*Prop III.2 - one cnot needed if \[Chi][\[Gamma][u]](x) = (x+I)^2(x-I)^2*)
If[Length[Select[N[eigenVals], Abs[# - I] < OptionValue[precision] &]] == 2 && Length[Select[N[eigenVals], Abs[# + I] < OptionValue[precision] &]] == 2,
Return[Module[{a,b,c,d,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a4,b4,c4,d4},
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix su in DecUnitary2QubitsHelp (case: one C-NOTs): ",su]
];
{a,b,c,d} = TwoQubitFindMatchingProductMatrices[su,CNOTM[1,2,2],OptionValue[precision],OptionValue[FullSimp]];
{d1,c1,b1,a1} = ZYZDecomposition[a];
{d2,c2,b2,a2} = ZYZDecomposition[b];
{d3,c3,b3,a3} = ZYZDecomposition[c];
{d4,c4,b4,a4} = ZYZDecomposition[d];
st = {
{3,b1,n},
{2,c1,n},
{3,d1,n},
{3,b2,n-1},
{2,c2,n-1},
{3,d2,n-1},
{0,n,n-1},
{3,b3,n},
{2,c3,n},
{3,d3,n},
{3,b4,n-1},
{2,c4,n-1},
{3,d4,n-1}
};

If[OptionValue[UpToDiagonal],
{{1,1,1,1}, st},
st
]
]
]
];

(*We now require two C-NOTS*)
(*
It might be that tr[Gam[su]] is not real
we can fix this, at the cost of being out by a diagonal gate
*)

If[Abs[Im[Tr[N[Gam[su]]]]] > OptionValue[precision] && OptionValue[UpToDiagonal],
traceFixingMatrix = Simplify[RealTrace3[su]];
su = Simplify[traceFixingMatrix.su];
,
traceFixingMatrix = IdentityMatrix[4];
];

If[Abs[Im[Tr[N[Gam[su]]]]] < OptionValue[precision],
(*Prop III.3 - two cnot needed if Tr[\[Gamma][u]] is real*)
Module[{x1,x2,del,phi,eVals,v,a,b,c,d,a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,a4,b4,c4,d4},
eVals = Sort[If[OptionValue[FullSimp],FullSimplifyNoRoots[Eigenvalues[Gam[su]]],Simplify[Eigenvalues[Gam[su]]]],ComplexOrderingFunction];
x1=Arg[eVals[[1]]];
x2=Arg[eVals[[3]]];
If[OptionValue[FullSimp],
del=Arg[FullSimplifyNoRoots[eVals[[1]]*eVals[[3]]]]/2;
phi=Arg[FullSimplifyNoRoots[eVals[[1]]/eVals[[3]]]]/2,del=Arg[Simplify[eVals[[1]]*eVals[[3]]]]/2;
phi=Arg[Simplify[eVals[[1]]/eVals[[3]]]]/2];
v=CNOTM[2,1,2].KroneckerProduct[SimplifyTrigo[RotGate[del,3]],SimplifyTrigo[RotGate[phi,1]]].CNOTM[2,1,2];
If[analyzeAnalyticDecUnitary2Qubits,
Print["Matrix su in DecUnitary2QubitsHelp (case: two C-NOTs): ",su];
Print["Matrix v in DecUnitary2QubitsHelp (case: two C-NOTs): ",v]
];
{a,b,c,d} = TwoQubitFindMatchingProductMatrices[su,v,OptionValue[precision],OptionValue[FullSimp]];
{d1,c1,b1,a1} = ZYZDecomposition[a];
{d2,c2,b2,a2} = ZYZDecomposition[b];
{d3,c3,b3,a3} = ZYZDecomposition[c];
{d4,c4,b4,a4} = ZYZDecomposition[d];
st = {
{3,b1,n},
{2,c1,n},
{3,d1,n},
{3,b2,n-1},
{2,c2,n-1},
{3,d2,n-1},
{0,n-1,n},
{3,del,n},
{1,phi,n-1},
{0,n-1,n},
{3,b3,n},
{2,c3,n},
{3,d3,n},
{3,b4,n-1},
{2,c4,n-1},
{3,d4,n-1}
};
If[OptionValue[UpToDiagonal],
Return[{Conjugate[Diagonal[traceFixingMatrix]], st}],
Return[st]
];
];
];

(*
Otherwise we need three cnots
see "http://web.eecs.umich.edu/~imarkov/pubs/jour/pra04-univ.pdf" - Minimal Universal Two-Qubit CNOT -based Circuits Theorem VI.3
*)
su = CNOTM[1,2,2].su;
traceFixingMatrix = Simplify[RealTrace3[su]];
su = Simplify[traceFixingMatrix.su];
st = DecUnitary2QubitsHelp[su, n, {precision-> OptionValue[precision], UpToDiagonal->False,FullSimp->OptionValue[FullSimp]}];
Return[Join[{{3,-Arg[traceFixingMatrix[[1,1]]/traceFixingMatrix[[2,2]]],n-1},{0,n,n-1}}, st]];
];

(*-------------------------------------------Basic helper methods (private) --------------------------------------------*)

(*Note that there is a bug in ArcTan[]: One gets 
ArcTan[1, 0. + 0 I]=0, but
ArcTan[1, 0. + 0. I]=-1.5708+0.I
This is fixed in ArcTan2 by ignoring the imaginary part (hence, ArcTan2 should only be used with real inputs).
*)
ArcTan2[x_, y_] := If[Chop[x] == 0 && Chop[y] == 0, 0, ArcTan[Re[x],Re[y]]]


isZeroMatrix[mat_]:=Tr[Abs[Flatten[Chop[mat]]]]==0;
isIdentity[m_]:=Norm[Chop[m-IdentityMatrix[Length[m]]]]==0;
isIdentityUpToPhase[m_]:=If[Chop[FullSimplifyNoRoots[Norm[m[[1,1]]]]]==0,False,
If[Chop[FullSimplifyNoRoots[Norm[m[[1,1]]]-1]]==0,
Chop[FullSimplifyNoRoots[Norm[m/(m[[1,1]])-IdentityMatrix[Length[m]]]]]==0,
False
]
];

isAnalyticGate[gate_]:=Module[{},
If[Length[gate[[2]]]==0,out=Not[MachineNumberQ[gate[[2]]]];,out=True;Do[If[MachineNumberQ[gate[[2]][[i]]],out=False;Goto[LabelEnd]],{i,1,Length[gate[[2]]]}];
Label[LabelEnd];
];
out
]

isListEqualOpUpToPhase[st_,mat_,free_:0]:=If[st!={},isIdentityUpToPhase[ConjugateTranspose[mat[[All,1;;2^(Log2[Dimensions[mat][[1]]]-free)]]].(CreateIsometryFromList[st][[All,1;;2^(Log2[Dimensions[mat][[1]]]-free)]])],isIdentityUpToPhase[mat]]

(*Generates the rotation matrix corresponding to the input angle \[Alpha] and pauli matrix i*)
RotGate[\[Alpha]_,i_]:=If[Dimensions[N[\[Alpha]]]==={},
MatrixExp[I \[Alpha]/2 PauliMatrix[i]],Throw[StringForm["RotGate called with input `1` rather than a single number",\[Alpha]]]]

RotGateM[\[Alpha]_,rotAxis_,i_,n_]:=
KroneckerProduct[IdentityMatrix[2^(i-1)],RotGate[\[Alpha],rotAxis],IdentityMatrix[2^(n-i)]]

(* Equivalent to the Sign function applied to (real parts of) elements of a list \
except SignZeroPos[0]=1 and SignZeroNeg[0]=-1 *)
SignZeroPos[x_] := 
 Module[{out = {}, i}, 
  For[i = 1, i <= Dimensions[x][[1]], i++, 
   out = Insert[out, If[Re[x[[i]]] < 0, -1, 1], -1]]; out]

SignZeroNeg[x_] := 
 Module[{out = {}, i}, 
  For[i = 1, i <= Dimensions[x][[1]], i++, 
   out = Insert[out, If[Re[x[[i]]] > 0, 1, -1], -1]]; out]
   
   (*Helpful for making the gate list numerical without making integers numerical *)
NGateList[st_]:=Module[{i,out={},a,b,c},For[i=1,i<=Length[st],i++,{a,b,c}=st[[i]];
out=Insert[out,If[MemberQ[{-2,1,2,3,100,101},a],{a,N[b],c},{a,b,c}],-1]];out]

(* Takes an angle and outputs angle' where 0\[LessEqual]angle'<2\[Pi]*)
AdjustAngleHelp[angl_]:=Mod[angl,2Pi]
   
(*Generates a string denoting the control and target of a cnot gate.
i denotes the control while j denotes the target.
*)
CNOTNotationStr[i_,j_]:=ToString[StringForm["C(`1`)(`2`)",i,j]]

CNOT[i_,j_]:={0,i,j}
Diag[diag_,act_]:={-2,diag,act}
Rx[angle_,act_]:={1,angle,act}
Ry[diag_,act_]:={2,diag,act}
Rz[diag_,act_]:={3,diag,act}
Mmt[act_]:={4,1,act}
TrOut[act_]:={4,0,act}
Ancilla[i_,act_]:={5,i,act}
PostSelect[i_,act_]:={6,i,act}

(*Generates a string denoting the control and target of a controlled z gate. 
i denotes the control while j denotes the target.
*)
CZNotationStr[i_,j_]:=ToString[StringForm["Cz(`1`)(`2`)",i,j]]

CZ[i_,j_]:={-1,i,j}

XX[phi_,i_,j_]:={100,phi,{i,j}}

RGate[th_,phi_,i_]:={101,{th,phi},i}

(*Given a,b,c,d, returns the matrix specified by e^ia*Rz(b).Ry(c).Rz(d) *)
ApplyZYZ[a_,b_,c_,d_]:=Exp[I*d]*RotGate[c,3].RotGate[b,2].RotGate[a,3]

(*Given a,b,c,d, returns the matrix specified by e^ia*Rx(b).Ry(c).Rx(d) *)
ApplyXYX[a_,b_,c_,d_]:=Exp[I*d]*RotGate[c,1].RotGate[b,2].RotGate[a,1]

(*Helper: Check if expression is bigger than zero*)
isBiggerThanZero[expr_]:=If[Chop[N[expr]]>0,True,False,Throw[StringForm["Error in isBiggerThanZero[]. It could not be determined
if the expression is bigger than zero or not."]]]

(*Helper: Calculate Arg after Simplifying*)
ArgAndSimplify[x_,timeConstraint_:300]:=TimeConstrained[FullSimplifyNoRoots[Quiet[Arg[ComplexExpand[FullSimplifyNoRoots[x]]],{N::meprec}]],timeConstraint,Arg[x]];

(*Helper: FullSimplify without using implicit "root" equations as outputs*)
FullSimplifyNoRoots[x_,timeConstraint_:300]:=
If[StringMatchQ[ToString[x],"*Re[*"]||StringMatchQ[ToString[x],"*Im[*"]||StringMatchQ[ToString[x],"*Arg[*"],
(*If we take Imaginary or Real parts or Arguments in the expression, we remove them by using ComplexExpand*)
TimeConstrained[FullSimplify[ComplexExpand[TrigToExp[x]],ComplexityFunction->(1000 Count[#,Root[__],All]+LeafCount[#]&)],timeConstraint,x],
TimeConstrained[FullSimplify[x,ComplexityFunction->(1000 Count[#,Root[__],All]+LeafCount[#]&)],timeConstraint,x]
]

(*
SimplifyIfLong[expr_, length_] :=
Module[{dim}, dim = Dimensions[N[expr]];
 If[dim === {},
  If[StringLength[ToString[expr, FormatType -> InputForm]] > length,
   FullSimplifyNoRoots[expr], expr],
  Map[SimplifyIfLong[#, length] &, expr, Length[dim]]]]
*)

(*Helper: FullSimplify trigo expressions*)
SimplifyTrigo[x_]:=FullSimplifyNoRoots[TrigToExp[x]]

NormalForm[z_]:=Re[z]+I*Im[z]
PolarForm[z_]:=Abs[z]*E^(I*Arg[z])

(*Takes the conjugates transpose. To handle analytic expressions, it computes the polar form, takes the transpose, and then transforms back to the normal form*)
CTSimplify[m_]:=Module[{},
out=Simplify[CT[m]];
If[StringMatchQ[ToString[out],"*Conjugate*"],
out=SimplifyTrigo[NormalForm[CT[ComplexExpand[m]]]]
];
out
]

ConjSimplify[m_]:=Module[{},
out=Simplify[Conjugate[m]];
If[StringMatchQ[ToString[out],"*Conjugate*"],
out=SimplifyTrigo[NormalForm[Conjugate[PolarForm[m]]]]
];
out
]

NormSimplify[mat_]:=Module[{},
out=Simplify[Norm[mat]];
If[StringMatchQ[ToString[mat],"*Norm[*"],
out=SimplifyTrigo[Norm[ComplexExpand[mat]]]
];
out
]

(*Takes a two dimensional vector vec and a basis state n (0 or 1) as an input. 
Outputs the SU[2] u, such that u.vec=r (1,0) if n=0 or u.vec=r (0,1) if n is 1, 
where r is the norm of the vector vec. *)
QubitInvert[vec_,n_]:=
Module[{psi,phi,mat,r,gate},(
If[MachineNumberQ[vec[[1,1]]]==True ||MachineNumberQ[vec[[2,1]]]==True ,
mat=N[IdentityMatrix[2]],
mat=IdentityMatrix[2]
];
r=NormSimplify[vec];
If[Chop[N[r]] == 0 || (n==0 && Chop[N[vec[[2,1]]]]==0) || (n==1 && Chop[N[vec[[1,1]]]]==0),
mat,
psi=Simplify[vec/r];
phi=FullSimplifyNoRoots[(-CTSimplify[psi].mat[[2]])[[1]]*mat[[1]]+(CTSimplify[psi].mat[[1]])[[1]]*mat[[2]]];
gate=KroneckerProduct[mat[[n+1]],CTSimplify[psi]]+KroneckerProduct[mat[[-(n+1)]],ConjSimplify[phi]*(-1)^n],
Throw[StringForm["Error: If condition in QubitInvert could did neither evaluate to True nor to False"]]
]
)]

(*Helper function: Reorders diagonal elements when the least significant qubit is swiched inbetween the others. Input is given by the diagonal
entries dia and the insert position pos*)
RelabelTwoQubitsDiag[dia_,pos_]:=Module[{i,bits,dia2,digits},(
bits=Log[2,Length[dia]];
dia2=dia;
For[i=0,i<Length[dia],i++,
digits=IntegerDigits[i,2,bits];
digits=Insert[digits,digits[[bits]],pos];
digits=Delete[digits,-1];
dia2[[FromDigits[digits,2]+1]]=dia[[i+1]];
];
dia2
)]

(*Helper function: Transform Rz rotation on most significant qubit to a diagonal gate*)
UCZToDiagonal[gates_]:=Module[{upperBlock,lowerBlock},(
upperBlock=Table[gates[[i]][[1,1]],{i,Length[gates]}];
lowerBlock=Table[gates[[i]][[2,2]],{i,Length[gates]}];
Join[upperBlock,lowerBlock]
)]  

(*Helper function: Transform 2\[Cross]2 block diagonal matrix (with zero non diagonal entries) 
into diagonal matrix*)
BlockToDiagonal[gates_]:=Module[{firstEntries,secondEntries},(
firstEntries=Table[gates[[i]][[1,1]],{i,Length[gates]}];
secondEntries=Table[gates[[i]][[2,2]],{i,Length[gates]}];
Riffle[firstEntries,secondEntries]
)]  

(*Helper function: checks if all gates are equal to Identity*)
IsUCGEqualToIdendity[gates_]:=Module[{},(
AllTrue[isIdentity/@gates,#==True&]
)]
isIdentity[m_]:=Chop[Norm[m-IdentityMatrix[Length[m]]]]==0;

(*Normalises the rows of matrix a - useful to form orthonormal bases from orthogonal bases*)
NormaliseRows[a_] := Module[{i,b},
b = a;
For[i=1, i<= Dimensions[a][[1]], i++, b[[i]] = Normalize[a[[i]]]];
Return[b]
]

(*Normalises the columns of matrix a - useful to form orthonormal bases from orthogonal bases*)
NormaliseCols[a_] := Module[{b},
Return[Transpose[NormaliseRows[Transpose[a]]]]
]

(*Multiply the single-qubit gates in the list representation st to get the matrix representaiton of the operator represented by st. (Ignores the qubit number and only
considers one qubit) *)
MultiplySQGates[st_,FullSimp_:True]:=Module[{mat,i,stTemp},(mat=IdentityMatrix[2^1];
For[i=1,i<=Length[st],i++,
stTemp={st[[i]][[1]],st[[i]][[2]],1};
mat=mat.SimplifyTrigo[ListFormToOp[stTemp,1]]];
mat
)];

(*Allows to change the global phase*)
rotIsZero[st_] := Module[{},(
If[Chop[Mod[N[st[[2]]],2Pi]]==0,True,False,Throw[StringForm["Error: rotIsZero did neither return True nor False"]]]
)]

(*Finds the highest qubit number in st (on which is acted non trivially)*)
NumberOfQubits[st_]:=Module[{},(
If[st=={},0,
Max[Map[ If[#[[1]] == -1 || #[[1]] == 0, 
Max[#[[2]], #[[3]]],
If[#[[1]]==-2||#[[1]]==100,Max[#[[3]]],
 #[[3]]]] &  ,st]]
 ]
 )]
 
 (*------------------------------------------- Adapted matrix decompositions (private)--------------------------------------------*)

(* Uses Mathematica's SingularValueDecomposition and then reorders the matrices such that 
in the diagonal matrix dp any singular values that is equal to 1 comes first *)
ReorderedSVD[m_,FullSimp_:True] := 
  Module[{dim, dimpos, up, dp, vp, posns, ord, permmat, pos, i}, 
    dim = Dimensions[m][[1]]; {up, dp, vp} = 
      SingularValueDecomposition[m]; 
    posns = Flatten[Position[Chop[N[Diagonal[dp] - 1]], 0, {1}]]; 
    If[Dimensions[posns][[1]] != 0, 
   ord = Join[posns, 
     Complement[Table[i, {i, 1, dim}], posns]]; {up[[ord, ord]], 
    dp[[ord, ord]], vp[[ord, ord]]}, If[FullSimp,FullSimplifyNoRoots[{up, dp, vp}],Simplify[{up, dp, vp}]] ]]  
    
    (* Does the QL decomposition but making the diagonal entries of the left triangular matrix L positive *)
QLDecompositionPos[m_,FullSimp_:True] := 
 Module[{q, l, s}, {q, l} = QLDecomposition[m]; 
  s = DiagonalMatrix[SignZeroPos[Diagonal[l]]];If[FullSimp, FullSimplifyNoRoots[{s.q, s.l}],Simplify[{s.q, s.l}]]]

(* Does the RQ decomposition making the diagonal entries of the right trianglar matrix R negative *)
RQDecompositionNeg[m_,FullSimp_:True] := 
 Module[{q, r, s}, {r, q} = RQDecomposition[m]; 
  s = DiagonalMatrix[-SignZeroNeg[Diagonal[r]]];If[FullSimp,FullSimplifyNoRoots[{r.s, q.s}],Simplify[{r.s, q.s}]]]

(*-------------------------------------------Decomposition of multi-controlled gates (public) --------------------------------------------*)

(*Decomposes a k-control Toffoli gate on n\[GreaterEqual]5 qubits with 3\[LessEqual]k\[LessEqual]\[LeftCeiling]n/2\[RightCeiling] up to a diagonal gate. It outputs the answer in terms of single-qubit gates and CNOTs up to a diagonal gate. The output is given in list representation and the (not specified) diagonal gate is ignored*)
DecToffoliMultiControlUpToDiagonal[tofcontrol_,toftarget_,tofbits_]:=
Module[{op1,st,st1,st2,st3,gate,left,centre,temp,xx,yy},(
op1=DecToffoliHalfCircuit[tofcontrol,toftarget,tofbits];
gate=op1[[1]];
op1=Drop[op1,1];
left=op1[[1;;Floor[Length[op1]/2]]];
centre=op1[[Floor[Length[op1]/2+1]]];
left=DecToffoliUpToDiagonal2/@left;
left=Reverse/@left;
centre=DecToffoliUpToDiagonal2@centre;
centre=Reverse@centre;
st1=left[[All,1;;4]];
temp=Delete[#,0]&/@st1;
st1=Reverse[left[[All,4;;7]]];
st1=Delete[#,0]&/@st1;
st1=Join[temp,centre,st1];
(*op2={ToffoliMultiC@gate};*)
st2=DecToffoliUpToDiagonal2@gate;
st2=Reverse@st2;
(*diag={Inverse[Dot@@op2].ToffoliMultiC[gate]};*)
st3=st2[[4;;7]];
st2=st2[[1;;4]];
st=Join[st1,st2,st1,st3];
Reverse[st]
)]

(*Decomposes a Toffoli gate into single-qubit rotations and CNOT gates up to a diagonal gate. Input: tofcontrols is a set containing the numbers of the two control qubits of the Toffoli gate. Output: the gate sequence in list form. Note that the diagonal is not given with the output*)
DecToffoliUpToDiagonal[tofcontrol_,toftarget_,tofbits_]:=
Module[{op,a,str,st},(
st={{2,-Pi/4,toftarget},CNOT[tofcontrol[[1]],toftarget],
{2,-Pi/4,toftarget},CNOT[tofcontrol[[2]],toftarget],
{2,Pi/4,toftarget},CNOT[tofcontrol[[1]],toftarget],
{2,Pi/4,toftarget}};
st
)]

(*Decomposes a Toffoli gate (up to global phase) into single-qubit rotations and CNOTs. Input: tofcontrols is a set containing the numbers of the two control qubits of the Toffoli gate. Output: gives the output in list form*)
DecToffoli[tofcontrol_,toftarget_,tofbits_]:=
Module[{st},(
st={{3,-Pi/4,tofcontrol[[2]]}(*E*),{3,Pi/2,toftarget}(*C*),CNOT[tofcontrol[[2]],toftarget],
{2,-Pi/4,toftarget}(*B*),CNOT[tofcontrol[[1]],toftarget],CNOT[tofcontrol[[1]],tofcontrol[[2]]],  
{3,Pi/4,tofcontrol[[2]]}(*E^{-1}*),{2,Pi/4,toftarget}(*B^{-1}*),CNOT[tofcontrol[[1]],tofcontrol[[2]]],
CNOT[tofcontrol[[2]],toftarget],{2,-Pi/4,toftarget}(*B*),CNOT[tofcontrol[[1]],toftarget],{3,-Pi/4,tofcontrol[[1]]}(*E*),
{2,Pi/4,toftarget}(*A part1*),{3,-Pi/2,toftarget}(*A part2*)};
st
)]

(*Decomposes k-controlled Toffoli gates (up to global phase) on n\[GreaterEqual]3 qubits for k\[GreaterEqual] 2. For k\[GreaterEqual] 3 we require that n\[GreaterEqual]k+2, i.e., we require one ancilla qubit 
(that can be in an arbitrary initial state and is given back in the same state). Outputs the decomposition into single-qubit rotations and C-NOTs in list form. 
Remark: Note that more than one qubit may be used as an ancilla, so you can not run the gate sequence in parallel with other gates on the qubits we neither control nor act on.*)
DecToffoliMultiControl[tofcontrol_,toftarget_,tofbits_]:=
Module[{k,st},(
k=Length[tofcontrol];
If[k==2,
st=Reverse[DecToffoli[tofcontrol,toftarget,tofbits]],
If[k<=  Ceiling[tofbits/2],
st=DecToffoliMultiControl1[tofcontrol,toftarget,tofbits],
st=DecToffoliMultiControl2[tofcontrol,toftarget,tofbits]
]
];
Reverse[st]
)
]

(*-------------------------------------------Decomposition of multi-controlled gates (private) --------------------------------------------*)

(*ToDo: Clear up gate order\[Rule]remove "Reverse[]".*)

(*Helper method: Takes as an input a multi-control Toffoli gate (described by tofcontrol,toftarget,n) 
with n\[GreaterEqual]5, where n denotes the total number of qubits, and where the number of controls k satisfies 3\[LessEqual]k\[LessEqual]\[LeftCeiling]n/2\[RightCeiling]. 
It outputs the first halfe of the decomposition in terms of Toffoli gates in notation form, i.e., given as tof={controls,target,n}*)
DecToffoliHalfCircuit[tofcontrol_,toftarget_,tofbits_]:=
Module[{free,op,i,j,gate,temp},(
free=Table[i,{i,1,tofbits}];
free=Complement[free,tofcontrol,{toftarget}];
op={{tofcontrol[[1;;2]],free[[1]],tofbits}};
For[i=3;j=1,i<Length[tofcontrol],i++;j++,
gate={{tofcontrol[[i]],free[[j]]},free[[j+1]],tofbits};
PrependTo[op,gate];
AppendTo[op,gate];
];
temp=op;
gate={{tofcontrol[[i]],free[[j]]},toftarget,tofbits};
PrependTo[op,gate];
(*AppendTo[op,gate];
op=Join[op,temp];*)
op
)]

(*Decomposes a k-control Toffoli gate on n\[GreaterEqual]5 qubits with 3\[LessEqual]k\[LessEqual]\[LeftCeiling]n/2\[RightCeiling]. 
It outputs the answer in terms of single-qubit gate rotations and CNOTs. 
Answer is given in list representation. Note that the output is split into action and 
reverse part to simlify the use of the Method in the decomposition for the MCGs
*)
DecToffoliMultiControlHelp[tofcontrol_,toftarget_,tofbits_]:=
Module[{op1,st,st1,st2,st3,gate,diag,res,left,centre,temp3,xx,yy},(
op1=DecToffoliHalfCircuit[tofcontrol,toftarget,tofbits];
gate=op1[[1]];
op1=Drop[op1,1];
left=op1[[1;;Floor[Length[op1]/2]]];
centre=op1[[Floor[Length[op1]/2+1]]];
left=DecToffoliUpToDiagonal2/@left;
left=Reverse/@left;
centre=DecToffoliUpToDiagonal2@centre;
centre=Reverse@centre;
st1=left[[All,1;;4]];
temp3=Delete[#,0]&/@st1;
st1=Reverse[left[[All,4;;7]]];
st1=Delete[#,0]&/@st1;
st1=Join[temp3,centre,st1];
{st2}=DecToffoliHelp@gate;
st2=Reverse@st2;
(*diag={Inverse[Dot@@op2].ToffoliMultiC[gate]};*)
st2=Drop[st2,5];
st3=st2;
st2=InverseGateList[st2];
st=Join[st2,st1,st3];
{st,st1}
)]

(*Decomposes a k-control Toffoli gate on n\[GreaterEqual]5 qubits with 3\[LessEqual]k\[LessEqual]\[LeftCeiling]n/2\[RightCeiling]. 
It outputs the answer in terms of single-qubit gate rotations and CNOTs. 
Answer is given in list representation.*)
DecToffoliMultiControl1[tofcontrol_,toftarget_,tofbits_]:=Module[{op,str,st,op1,str1,st1},(
{st,st1}=DecToffoliMultiControlHelp[tofcontrol,toftarget,tofbits];
st=Join[st1,st];
st
)]

(*Same as DecToffoli, but with input given as a list tof={controls,target,n}*)
DecToffoliHelp[{tofcontrol_,toftarget_,tofbits_}]:={DecToffoli[tofcontrol,toftarget,tofbits]}

(*Decomposes k-controlled Toffoli gate on n\[GreaterEqual]5 qubits for k\[LessEqual] n-2, i.e., requires one ancilla qubit 
(that can be in an arbitrary state and is given back in the same state). Outputs the Decomposition splitted in two parts. 
The two parts are required for DecMCSpecialUnitary[]. *)
DecToffoliMultiControl2Help[tofcontrol_,toftarget_,tofbits_]:=
Module[{control2,control1,st1,st2,st2Reverse,stPart1,stReversePart,controlAndAncilla,numberControl2,control2temp,ancilla},(
numberControl2=Ceiling[(Length[tofcontrol]+2)/2];
control2temp=tofcontrol[[-(numberControl2-1);;-1]];
ancilla=Complement[Range[tofbits],Join[tofcontrol,{toftarget}]][[1]];(*Find free ancilla*)
control1=Complement[tofcontrol,control2temp];
control2=Sort[Join[control2temp,{ancilla}]];
If[Length[control1]==2,
st1=Reverse[DecToffoli[control1,ancilla,tofbits]],
st1=Reverse[DecToffoliMultiControlUpToDiagonal[control1,ancilla,tofbits]];
];
{st2,st2Reverse}=DecToffoliMultiControlHelp[control2,toftarget,tofbits];
stPart1=Join[st1,Reverse[st2],Reverse[st2Reverse],st1,Reverse[st2]];
stReversePart=Reverse[st2Reverse];
{stPart1,stReversePart}
)
]

(*Decomposes k-controlled Toffoli gate on n\[GreaterEqual]5 qubits for \[LeftCeiling]n/2\[RightCeiling]<k\[LessEqual] n-2, i.e., requires one ancilla qubit 
(that can be in an arbitrary state and is given back in the same state). Outputs the Decomposition into single-qubit rotations and C-NOTs in list form. *)
DecToffoliMultiControl2[tofcontrol_,toftarget_,tofbits_]:=
Module[{},(
{stPart1,stReversePart}=DecToffoliMultiControl2Help[tofcontrol,toftarget,tofbits];
Reverse[Join[stPart1,stReversePart]]
)
]

(*Helper method: Same as DecToffoliUpToDiagonal, but the input is given as a list tof={controls, target,n}*)

DecToffoliUpToDiagonal2[{tofcontrol_,toftarget_,tofbits_}]:=DecToffoliUpToDiagonal[tofcontrol,toftarget,tofbits]

(*---------------------------------------Methods to decompose multi controlled single-qubit gates (public)------------------------------------*)

(*Decomposes a one-controlled single-qubit unitary into CNOTs and single-qubit gates (see Lemma 6 in "Quantum circuits for isometries"). Takes as input the single-qubit gate u to be applied and the qubit number of the control and target qubit as well as the total number of qbits. Note that since this is single control, the value of control that is input should be a number(not a list like in the cases of multiple controls). Outputs a decomposition in matrix, string form and list form.*)
DecSingleMCG[u_, control_, target_, bits_] :=
 Module[{\[Alpha], \[Beta], \[Gamma], \[Delta], a, b, c,  mat, tt, v, st}, (
{\[Delta],\[Gamma],\[Beta],\[Alpha]} =ZYZDecomposition[u];
a = RotGate[\[Beta], 3].RotGate[\[Gamma]/2, 2];
b = RotGate[-\[Gamma]/2, 2].RotGate[-(\[Delta] + \[Beta])/2, 3];
c = RotGate[(\[Delta] - \[Beta])/2, 3];
st={{3,Chop@(\[Delta] - \[Beta])/2,target},CNOT[control,target],{3,Chop@-(\[Delta] + \[Beta])/2,target},{2,Chop@-\[Gamma]/2,target},CNOT[control,target],{2,Chop@\[Gamma]/2,target},{3,Chop@\[Beta],target}};
If[Chop[N[\[Alpha]]]!=0,
PrependTo[st,{3,Chop@-\[Alpha],control}];
];
st
)]

(*Decomposes a multi-controlled gate into single-qubit rotations and CNOTs.*)
DecMCG[u_,control_,target_,bits_]:=
Module[{v,st,st1,st2,st3,st2Reverse,toffoliGateMC},(
If[Length[control]==1,
st=DecSingleMCG[u,control[[1]],target,bits],
v=MatrixPower[u,1/2];
st1=Reverse[DecSingleMCG[v,control[[-1]],target,bits]];
st2=Reverse[DecSingleMCG[CTSimplify[v],control[[-1]],target,bits]];
If[Length[control]==2,
st3=Reverse[DecSingleMCG[v,Complement[control,{control[[-1]]}][[1]],target,bits]];
toffoliGateMC={CNOT[Drop[control,-1][[1]],control[[-1]]]};
,
	st3=Reverse[DecMCG[v,Complement[control,{control[[-1]]}],target,bits]];
If[Length[control]==3,
toffoliGateMC=Reverse[DecToffoliUpToDiagonal[Drop[control,-1],control[[-1]],bits]],
toffoliGateMC=Reverse[DecToffoliMultiControl[Drop[control,-1],control[[-1]],bits]];
]
];
st=Join[Reverse[st1],Reverse[toffoliGateMC],Reverse[st2],Reverse[toffoliGateMC],Reverse[st3]];
st=st;
];
st
)]

(*Decomposes a k-controlled special 2\[Cross]2 unitary u\[Element]SU(2) into single-qubit rotations and CNOTs. Works for the case k\[GreaterEqual]7. Answer is given in both matrix and string form*)
DecMCSpecialUnitary[u_,control_,target_,bits_]:=
Module[{gate,temp,angles,a,b,cAngle,sta,stb,st2,st2Reverse,st1,st,toffoliPart1,reversePart,stc},(
angles=ZYZDecomposition[u];
(*Since u is a special unitary, we know that the global phase is +1 or -1. We can remove the -1 with the z-rotation.*)
If[Sign[Re[N[ApplyZYZ[angles[[1]],angles[[2]],angles[[3]],0]]]]==Sign[Re[N[u]]],,angles[[3]]=angles[[3]]+2Pi,Throw[StringForm["Error: If condition in DecMCSpecialUnitary did neither return True nor False"]]];
(*Notation of the gates corresponding to Lemma 7.9 here: https://arxiv.org/pdf/quant-ph/9503016.pdf. Note that the gates are a,b,c are in reversed order in Lemma 4.3. *)
a = ApplyZYZ[0,angles[[2]]/2,angles[[3]],0];
b = ApplyZYZ[-(angles[[1]] + angles[[3]])/2,-angles[[2]]/2,0,0];
cAngle =(angles[[1]] - angles[[3]])/2;
sta=Reverse[DecSingleMCG[a,control[[-1]],target,bits]];
stb=Reverse[DecSingleMCG[b,control[[-1]],target,bits]];
(*The gate c is a controlled z gate for which there exists a simpler decomposition*)
stc={CNOT[control[[-1]],target],{3,Chop@(-cAngle/2),target},CNOT[control[[-1]],target],{3,Chop@(cAngle/2),target}};
(*Decompose the Toffoli gates*)
{toffoliPart1,st2Reverse}=DecToffoliMultiControl2Help[Drop[control,-1],target,bits];
st=Join[Reverse[stc],toffoliPart1,Reverse[stb],InverseGateList[toffoliPart1],Reverse[sta]];
st
)];

(*Decompose a  k-controlled special unitary u\[Element]SU(2). Note that we don't need any ancillas for the decomposition. 
There is an option to return also the diagonal gate in the list form as 
{-2,list of diagonal entries, action qubits}. Note that if no diagonal
gate is returned (however ReturnDiagonal is set to True), the diagonal gate is equal to the identity operation.*)
Options[DecMCSpecialUnitaryUpToDiagonal] = {ReturnDiagonal->False,FullSimp->True};
DecMCSpecialUnitaryUpToDiagonal[u_,control_,target_,bits_,OptionsPattern[]]:=
Module[{len,st,m,i,diag,op},(
len=Length@control;
If[len==0,
If[OptionValue[ReturnDiagonal]==False,
Drop[ZYZDec[u,target],-1],
st=ZYZDec[u,target];
If[st[[-1]][[1]]==3,
op=ListFormToOp[{st[[-1]][[1]],st[[-1]][[2]],1},1];
st=Drop[st,-1];
diag={-2,Diagonal[op],{target}};
AppendTo[st,diag];
];
st,
Throw[StringForm["Error: If condition checking 'OptionValue[ReturnDiagonal]==False' in DecMCSpecialUnitaryUpToDiagonal did neither return True nor False"]]
],
If[len<=7,
m={};
(*Check if the input was numerical or analytical (analytic input slows down the computation)*)
If[MachineNumberQ[u[[1,1]]]==True  || MachineNumberQ[u[[1,2]]]|| MachineNumberQ[u[[2,1]]]|| MachineNumberQ[u[[2,2]]],
For[i=1,i<2^Length[control],i++,AppendTo[m,N[IdentityMatrix[2]]]],
For[i=1,i<2^Length[control],i++,AppendTo[m,IdentityMatrix[2]]];
];
AppendTo[m,u];
If[OptionValue[ReturnDiagonal]==False,
st=DecUCGUpToDiagonal[m,control,target,bits,FullSimp->OptionValue[FullSimp]];,
st=DecUCGUpToDiagonalHelp[m,control,target,bits,OptionValue[FullSimp]],
Throw[StringForm["Error: If condition checking 'OptionValue[ReturnDiagonal]==False' in DecMCSpecialUnitaryUpToDiagonal did neither return True nor False"]]
];
st
,
DecMCSpecialUnitary[u,control,target,bits]
]
]
)]

(*---------------------------------------Methods to decompose multi controlled single-qubit gates (private)------------------------------------*)
(*ToDo: Fix this improved version of DecTwoControlGate[]. The list output is not correct.*)

(*
(*Decomposes a two-controlled gate into single-qubit rotations and CNOTs. CNOT count=6 (see Lemma 8 in "Quantum circuits for isometries"). The decomposition is a generalized version of the code for DecToffoli (note if this code is used to decompose a Toffoli gate it will also contain some identity matrices)*)
DecTwoControlGate[u_,control_,target_,bits_]:=
Module[{v,op,a,b,c,str,st,e1,e2,\[Alpha], \[Beta], \[Gamma], \[Delta]},(
v=MatrixPower[u,0.5];
{\[Alpha], \[Beta], \[Gamma], \[Delta]} =ZYZDecomposition[v];

a=KroneckerProduct@@ReplacePart[Table[IdentityMatrix[2],{bits}],target->#]&/@{RotGate[\[Beta],3],RotGate[\[Gamma]/2,2]};

b=KroneckerProduct@@ReplacePart[Table[IdentityMatrix[2],{bits}],target->#]&/@{RotGate[-\[Gamma]/2, 2],RotGate[-(\[Delta] + \[Beta])/2, 3]};

c=KroneckerProduct@@ReplacePart[Table[IdentityMatrix[2],{bits}],target->RotGate[(\[Delta] - \[Beta])/2, 3]];

e1=KroneckerProduct@@ReplacePart[Table[IdentityMatrix[2],{bits}],control[[2]]->RotGate[-\[Alpha], 3]*Exp[I \[Alpha]/2]];
e2=KroneckerProduct@@ReplacePart[Table[IdentityMatrix[2],{bits}],control[[1]]->RotGate[-\[Alpha], 3]*Exp[I \[Alpha]/2]];

op={a[[1]],a[[2]],e2,CNOT[control[[1]],target,bits],b[[1]],b[[2]],CNOT[control[[2]],target,bits],CNOT[control[[1]],control[[2]],bits],b[[2]]\[ConjugateTranspose],b[[1]]\[ConjugateTranspose],e1\[ConjugateTranspose],CNOT[control[[1]],control[[2]],bits],CNOT[control[[1]],target,bits],b[[1]],b[[2]],CNOT[control[[2]],target,bits],c,e1};
str={ToString[StringForm["Rz(`2`)(`1`)",Chop@\[Beta],target]],ToString[StringForm["Ry(`2`)(`1`)",Chop@\[Gamma]/2,target]],ToString[StringForm["Rz(`2`)(`1`)",Chop@-\[Alpha],control[[1]]]],CNOTNotationStr[control[[1]],target],ToString[StringForm["Ry(`2`)(`1`)",Chop@-\[Gamma]/2,target]],ToString[StringForm["Rz(`2`)(`1`)",Chop@-(\[Delta] + \[Beta])/2,target]],CNOTNotationStr[control[[2]],target],CNOTNotationStr[control[[1]],control[[2]]],ToString[StringForm["Rz(`2`)(`1`)",Chop@(\[Delta] + \[Beta])/2,target]],ToString[StringForm["Ry(`2`)(`1`)",Chop@\[Gamma]/2,target]],ToString[StringForm["Rz(`2`)(`1`)",Chop@\[Alpha],control[[1]]]],CNOTNotationStr[control[[1]],control[[2]]],CNOTNotationStr[control[[1]],target],ToString[StringForm["Ry(`2`)(`1`)",Chop@-\[Gamma]/2,target]],ToString[StringForm["Rz(`2`)(`1`)",Chop@-(\[Delta] + \[Beta])/2,target]],CNOTNotationStr[control[[2]],target],ToString[StringForm["Rz(`2`)(`1`)",Chop@(\[Delta] - \[Beta])/2,target]],ToString[StringForm["Rz(`2`)(`1`)",Chop@-\[Alpha],control[[2]]]]};st={{3,Chop@\[Beta],target},{2,Chop@\[Gamma]/2,target},{3,Chop@-\[Alpha],control[[1]]},CNOT[control[[1]],target],{2,Chop@-\[Gamma]/2,target},{3,Chop@-(\[Delta] + \[Beta])/2,target},CNOT[control[[2]],target],CNOT[control[[1]],control[[2]]],{3,Chop@(\[Delta] + \[Beta])/2,target},{2,Chop@\[Gamma]/2,target},{3,Chop@\[Alpha],control[[1]]},CNOT[control[[1]],control[[2]]],CNOT[control[[1]],target],{2,Chop@-\[Gamma]/2,target},{3,Chop@-(\[Delta] + \[Beta])/2,target},CNOT[control[[2]],target],{3,Chop@(\[Delta] - \[Beta])/2,target},{3,Chop@-\[Alpha],control[[2]]}};

If[Chop[\[Alpha]]==0,
op=Delete[op,{{3},{11},{18}}];
str=Delete[str,{{3},{11},{18}}];st=Delete[st,{{3},{11},{18}}];
];
{op,str,st}
)]
*)

(*---------------------------------------Methods to create and handle multi controlled gates (public)------------------------------------*)

(*Generates a matrix representing a multi-control Toffoli gate described by
control,target,n where control is a list of the control qubits,target is an integer denoting the target qubit and n is the total number of qubits*)
CreateMCToffoli[control_,target_,n_]:=
Module[{part,vec,gate,u0,u1,spinUp,spinDown},(
spinUp={1,0};
spinDown={0,1};
u0=Outer[Times,spinUp,spinUp];
u1=Outer[Times,spinDown,spinDown];
vec=Table[IdentityMatrix[2],{n}];
If[Length[vec]>1,
gate=KroneckerProduct@@vec,
vec[[1]]
];
part=Fold[ReplacePart[#1,#2->u1]&,vec,control];
If[Length[vec]>1,
gate=gate-(KroneckerProduct@@part),
gate-part[[1]]
];
part=ReplacePart[part,target->PauliMatrix[1]];
If[Length[vec]>1,
gate=gate+(KroneckerProduct@@part),
gate=gate+part[[1]]
];
gate
)]

(*Generates the matrix representation of a multi controlled 2\[Cross]2 unitary u with control qubits in the list control and a target qubit target, and acting on a total number of qubits n*)
CreateMCG[u_,control_,target_,n_]:=
Module[{part,vec,gate,u1,spinDown},(
spinDown={0,1};
u1=Outer[Times,spinDown,spinDown];
vec=Table[IdentityMatrix[2],{n}];
If[Length[vec]>1,
gate=KroneckerProduct@@vec,
vec[[1]]
];
part=Fold[ReplacePart[#1,#2->u1]&,vec,control];
If[Length[vec]>1,
gate=gate-(KroneckerProduct@@part),
gate=gate-part[[1]]
];
part=ReplacePart[part,target->u];
If[Length[vec]>1,
gate=gate+(KroneckerProduct@@part),
gate=gate+part[[1]]
];
gate
)]

(*Extracts the single-qubit gate from a multi controlled gate (mcg). Note that u is the matrix representation of the whole mcg.*)
ExtractMCG[u_,control_,target_,bits_]:=
Module[{gate,vec,up,down,row,column},(
row=Total[2^(bits-control)]+1;
gate=ConstantArray[0,{2,2}];
gate=Table[u[[i,j]],{i,{row,row+2^(bits-target)}},{j,{row,row+2^(bits-target)}}];
gate
)]

(*Creates the matrix representation of a uniformly controlled gate*)
CreateUCG[gates_,control_,target_,bits_]:=
Module[{ucg,ex,bin,i,j,k,m,n,row,temp,mat,spinUp,spinDown},(
spinUp={1,0};
spinDown={0,1};
mat={KroneckerProduct[spinUp,spinUp],KroneckerProduct[spinDown,spinDown]};
ex=2^(bits-control);
ucg=ConstantArray[0,{2^bits,2^bits}];
For[k=0,k<2^Length[control],k++,
bin=IntegerDigits[k,2,Length[control]];
temp=Table[IdentityMatrix[2],{bits}];
temp=ReplacePart[temp,Table[control[[i]]->mat[[bin[[i]]+1]],{i,1,Length[control]}]];
temp=ReplacePart[temp,target->gates[[k+1]]];
If[Length[temp]>1,
temp=KroneckerProduct@@temp,
temp=temp[[1]]
];
ucg=ucg+temp;
];
ucg
)]

(*Extracts the single-qubit gates from a uniformly controlled gate (ucg) matrix representation. Note that u is the the matrix representation of the whole ucg.*)
ExtractUCG[u_,control_,target_,bits_]:=
Module[{gates,ex,bin,i,j,k,row},(
ex=2^(bits-control);
gates={};
For[k=0,k<2^Length[control],k++,
bin=IntegerDigits[k,2,Length[control]];
row=Total[bin*ex]+1;
AppendTo[gates,Table[u[[i,j]],{i,{row,row+2^(bits-target)}},{j,{row,row+2^(bits-target)}}]];
];
gates
)]

(*---------------------------------------Decomposition of uniformly controlled gates (public)------------------------------------*)
(*Decomposes a UCG up to a diagonal gate into CNOTs and single-qubit rotations. 
 Answer is given in list form. The diagonal gate is ignored in the output by default. If ReturnDiagonal\[Rule]"True",
 the output is {diagGate,st} with diagGate = {\[Minus]2,diag,actionQubits}*)
Options[DecUCGUpToDiagonal] = {ReturnDiagonal->False,FullSimp->True};
DecUCGUpToDiagonal[gates_,control_,target_,bits_,OptionsPattern[]]:=Module[{st},(
st=DecUCGUpToDiagonalHelp[gates,control,target,bits,OptionValue[FullSimp]];
st=RelabelQubits[st,Sort[control],control];
If[OptionValue[ReturnDiagonal]==False,
Drop[st,-1]
,
st,
Throw[StringForm["Error: If condition checking 'OptionValue[ReturnDiagonal]==False' in ecUCGUpToDiagonal did neither return True nor False"]]
]
)
]

(*---------------------------------------Decomposition of uniformly controlled gates (private)------------------------------------*)

(*Invokes the DecUCGHelp function to decompose the multi-control ucg gate and decomposes the single-qubit gates into rotations. Diagonal gate remains as it is. Answer is given in list form. The diagonal gate is indicated in the output.*)
(*ToDo: Recognize if some of the controls are not necessary, e.g., for gates={a,b,a,b}*)
DecUCGUpToDiagonalHelp[gates_,controlIn_,target_,bits_,FullSim_:True]:=
Module[{dim,op,op1,op2,str,str1,str2,st,st1,st2,a,u,i,j,max,targetTemp,isSwitched,diagRelabled,diag,actionDiag,control},(
(*Switch target to an additional ancilla qubit if it is not already least significant (compared to the control qubits) *)
control=Sort[controlIn];
max=Max[control];
isSwitched=False;
If[max> target,isSwitched=True;targetTemp=bits+1,targetTemp=target];
(*Decompose the UCG with switched target*)
st1=DecUCGUpToDiagonalHelp2[gates,control,targetTemp,bits,FullSim];
(*Decompose single qubit gates*)
st={};
For[i=2,i<=Length[st1],i=i+2,
AppendTo[st,st1[[i-1]]];
If[analyzeAnalyticCCDec,Print["Calculating ZYZDecomposition in DecUCGUpToDiagonalHelp with input u= ",st1[[i]][[1]]]];
a=ZYZDecomposition[st1[[i]][[1]]];
If[analyzeAnalyticCCDec,Print["Finished calculating ZYZDecomposition in DecUCGUpToDiagonalHelp"]];
st2={{3,Chop@a[[3]],targetTemp},{2,Chop@a[[2]],targetTemp},{3,Chop@a[[1]],targetTemp}};
st=Join[st,st2];
];
If[isSwitched==True,
(*Move the target back*)
diag=st[[1]];
diag[[3]]=Sort[diag[[3]]];
(*Relabel the output list*)
st=RelabelQubits[Drop[st,1],{bits+1,target},{target,bits+1}];
(*Reorder the elements of the diagonal gate*)
diagRelabled=RelabelTwoQubitsDiag[diag[[2]],FirstPosition[diag[[3]],x_/;x>target]];
actionDiag=Delete[diag[[3]],Position[diag[[3]],bits+1]];
actionDiag=Sort[Append[actionDiag,target]];
PrependTo[st,{-2,diagRelabled,actionDiag}];
];
Reverse[st]
)]

(*Decomposes a uniformly controlled gate into single-qubit gates(not decomposed into rotations) and CNOTs upto a diagonal. 
Works for any number of controls. The number of the target quit must be higher than all the control qubit numbers (i.e., the target should be least significant compared to the control qubits).*)
DecUCGUpToDiagonalHelp2[gates_,control_,target_,bits_,FullSim_:True]:=
Module[{r1temp,tempR,a,b,u,v,r,r1,r2,d,x,x1,det,h,op,op1,op2,str,str1,str2,st,st1,st2,sys,matu,matv,start,temp,temp2,vals,vecs,rnew,bini,binj,i,j,matr,diag,RzhList,hList,uprime,vprime,RzList,diagRz,singleRZ},(
If[analyzeAnalyticCCDec,Print["Calling DecUCGUpToDiagonalHelp2 with input {gates,control,target,bits}= ",{gates,control,target,bits}]];
a=gates[[1;;Length[gates]/2]];
b=gates[[Length[gates]/2+1;;-1]];
If[FullSim,x=FullSimplifyNoRoots[Table[a[[i]].CTSimplify[b[[i]]],{i,1,Length[a]}]];
det=FullSimplifyNoRoots[Det/@x],x=Simplify[Table[a[[i]].CT[b[[i]]],{i,1,Length[a]}]];
det=Simplify[Det/@x]];
If[analyzeAnalyticCCDec,Print["Calculating x1 in DecUCGUpToDiagonalHelp2"]];
x1=x/Sqrt[det];
If[FullSim,x1=FullSimplifyNoRoots[#[[1,1]]&/@x1],x1=Simplify[#[[1,1]]&/@x1]];
If[analyzeAnalyticCCDec,Print["Finished caclulating x1 in DecUCGUpToDiagonalHelp2  with output=",x1]];
If[analyzeAnalyticCCDec,Print["Caclulating r1,r2 in DecUCGUpToDiagonalHelp2 with {det,x1}= ",{det,x1}]];
If[FullSim,r1=FullSimplifyNoRoots[Table[Exp[I (Pi/2-Arg[det[[i]]]/2-Arg[x1[[i]]])/2],{i,1,Length[x]}]];
r2=FullSimplifyNoRoots[Table[Exp[I (Pi/2-Arg[det[[i]]]/2+Arg[x1[[i]]]+Pi)/2],{i,1,Length[x]}]],
r1=Simplify[Table[Exp[I (Pi/2-Arg[det[[i]]]/2-Arg[x1[[i]]])/2],{i,1,Length[x]}]];
r2=Simplify[Table[Exp[I (Pi/2-Arg[det[[i]]]/2+Arg[x1[[i]]]+Pi)/2],{i,1,Length[x]}]];];
(*ToDo: Replace with version that is easier to handle analytically...e.g., 
r1=FullSimplifyNoRoots[Table[tempR=Exp[I (Pi/4)]*(det[[i]])^(-1/4)*(x1[[i]])^(-1/2);If[Chop[N[NormSimplify[tempR]]]\[Equal]0,1,tempR/NormSimplify[tempR]],{i,1,Length[x]}]];
But: One has to handle the cases where we devide by zero!*)

If[analyzeAnalyticCCDec,Print["Finished caclulating r1,r2 in DecUCGUpToDiagonalHelp2"]];
r=Table[{{r1[[i]],0},{0,r2[[i]]}},{i,1,Length[r1]}];
If[analyzeAnalyticCCDec,Print["Calling EigensystemExact for inputs r[[i]].x[[i]].r[[i]]] with r= ",r," and x= ", x]];
If[FullSim,sys=FullSimplifyNoRoots[Table[EigensystemExact[Simplify[r[[i]].FullSimplifyNoRoots[x[[i]].r[[i]]]]],{i,1,Length[x]}]],
sys=Simplify[Table[EigensystemExact[Simplify[r[[i]].Simplify[x[[i]].r[[i]]]]],{i,1,Length[x]}]]];
If[analyzeAnalyticCCDec,Print["Finished call of EigensystemExact. Found sys= ",sys]];
vals=sys[[All,1]];
vecs=sys[[All,2]];
If[FullSim,vecs=FullSimplifyNoRoots[Map[Map[Normalize,#]&,vecs]],vecs=Simplify[Map[Map[Normalize,#]&,vecs]]];

For[i=1,i<=Length[sys[[All,1]]],i++,
If[Chop[N[vals[[i]]+{I,-I}]]=={0,0},
vecs[[i]]=Reverse[vecs[[i]]];
vals[[i]]=Reverse[vals[[i]]];
,,
Throw[StringForm["Error: If condition checking 'Chop[N[vals[[i]]+{I,-I}]]=={0,0}' in DecUCGUpToDiagonalHelp2 did neither return True nor False"]]
];
];
If[FullSim,d=FullSimplifyNoRoots[DiagonalMatrix/@Sqrt[vals]];
u=FullSimplifyNoRoots[Transpose/@vecs],d=Simplify[DiagonalMatrix/@Sqrt[vals]];
u=Simplify[Transpose/@vecs]];
If[analyzeAnalyticCCDec,Print["Calculating table v for decomposing UCG..."]];
If[FullSim,v=Table[SimplifyTrigo[d[[i]].CTSimplify[u[[i]]].CTSimplify[r[[i]]].b[[i]]],{i,1,Length[d]}],
v=Table[Simplify[d[[i]].CT[u[[i]]].CT[r[[i]]].b[[i]]],{i,1,Length[d]}]];
If[analyzeAnalyticCCDec,Print["Found table v: ",v]];

If[Length[control]==1,
st2={{v[[1]],{},target,bits}};
st1={{u[[1]],{},target,bits}};
,
st2=DecUCGUpToDiagonalHelp2[SimplifyTrigo[v],Drop[control,1],target,bits,FullSim];
If[analyzeAnalyticCCDec,Print["Finished running DecUCGUpToDiagonalHelp2"]];
(*Absorb diagonal gate diagRz into uniformly controlled gate u*)
diagRz=st2[[1]][[2]];
u=Table[u[[i]].DiagonalMatrix[{diagRz[[2*i-1]],diagRz[[2*i]]}],{i,Length[u]}];
st2=Drop[st2,1];
st1=DecUCGUpToDiagonalHelp2[SimplifyTrigo[u],Drop[control,1],target,bits,FullSim];
If[analyzeAnalyticCCDec,Print["Finished running DecUCGUpToDiagonalHelp2"]];
];
(*Absorb Hadamard gate into v*)
If[FullSim,st1=FullSimplifyNoRoots[Expand[st1]];
st2=FullSimplifyNoRoots[Expand[st2]],st1=Simplify[Expand[st1]];
st2=Simplify[Expand[st2]]];
h={{1,1},{1,-1}}/Sqrt[2];
st2[[1]][[1]]=h.st2[[1]][[1]];
(*Absorg Z-rotation and Hadamard gate into u*)
st1[[-1]][[1]]=st1[[-1]][[1]].RotGate[Pi/2,3].h;
(*Merge Rz gate on first qubit into the uniformly controlled Rz rotaiton gate r*)
If[FullSim,r=Join[CTSimplify/@r,r],r=Join[CT/@r,r]];(*Create full gate r*)
r=BlockToDiagonal[r];(*write as a diagonal gate*)
singleRZ=Table[RotGate[Pi/2,3],{i,Length[r]/2}];
diagRz=UCZToDiagonal[singleRZ];
r=Simplify[r*diagRz];(*Absorb Rz rotation into r gate*)
If[Length[control]!=1,
(*Merge two diagonal gates together at the end of the circuit*)
If[analyzeAnalyticCCDec,Print["Simplifying output in DecUCGUpToDiagonalHelp2..."]];
If[FullSim,r=FullSimplifyNoRoots[r*Join[st1[[1]][[2]],st1[[1]][[2]]]],r=Simplify[r*Join[st1[[1]][[2]],st1[[1]][[2]]]]];(*Enlarge diagonal st1[[1]][[2]] to one qubit more by using Join[..]*)
If[analyzeAnalyticCCDec,Print["Finished simplifying output in DecUCGUpToDiagonalHelp2"]];
st1=Drop[st1,1];
];
st=Join[{{-2,r,Join[control,{target}]}},st1,{CNOT[control[[1]],target]},st2];
st
)]

(*---------------------------------------Column-by-column decomposition (public)------------------------------------*)

(*Takes an 2^n \[Times] 2^m isometry with n\[GreaterEqual] 2 as an input and and decomposes it column by column into C-NOTs and single-qubit rotations. See "Quantum circuits for isometries" for the details of the decomposition.*)
Options[ColumnByColumnDec] = {FirstColumn->"UCG",Simp->True,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
Block[{$MaxExtraPrecision=500},
ColumnByColumnDec[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]]:=
Module[{actionQ,i,u2,st,st1,numCols,bits,free,dia,st2,j,isAnalytic,dim},(
IsQubitIsometry[u,"ColumnByColumnDec"];
If[StringQ[OptionValue[FirstColumn]],,Throw[StringForm["Error in ColumnByColumnDec: Value for option FirstColumn must be a String (UCG or StatePreparation)."]]];
If[OptionValue[FirstColumn]=="StatePreparation"||OptionValue[FirstColumn]=="UCG",,Throw[StringForm["Error in ColumnByColumnDec: String for option FirstColumn must be UCG or StatePreparation."]]];
isAnalytic=True;
If[OptionValue[FullSimp],u2=FullSimplifyNoRoots[u],u2=Simplify[u]];dim=Dimensions[u2];
For[i=1,i<=dim[[1]],i++,
For[j=1,i<=dim[[2]],i++,
If[MachineNumberQ[u2[[i,j]]]==True,isAnalytic=False;Goto["endOfLoop"]];
];
];
Label["endOfLoop"];
numCols=dim[[2]];
bits=Log2[dim[[1]]];free=Log2[dim[[1]]/dim[[2]]];
If[bits==1,
If[numCols==1,st=InverseGateList[Reverse[StatePrep1Qubit[u2,,{IgnoreAncilla->True}]]],
st=InverseGateList[Reverse[QSD[u2,,{Simp->OptionValue[Simp],IgnoreAncilla->True,FullSimp->OptionValue[FullSimp]}]]]
]
];
If[bits>=2,st={};
st1={};
st2={};
(*ToDo: Improve efficency by only saving the phases of the already implemented columns*)
For[i=0,i<numCols,i++,
If[i==0&& (OptionValue[FirstColumn]=="StatePreparation"),
st1=InverseGateList[Reverse[StatePreparation[Transpose[{u2[[All,1]]}],,{Simp->OptionValue[Simp],IgnoreAncilla->True,FullSimp->OptionValue[FullSimp]}]]];
u2=CreateIsometryFromList[Reverse[st1],bits,FullSimp->OptionValue[FullSimp]].u2;
, 
If[i==0&& (OptionValue[FirstColumn]!= "UCG"),Throw[StringForm["Error: Option valus 'FirstColumn' in ColumnByColumnDec is unknown."]]];
{st1, u2} = DecSingleColumn[u2,i, bits,isAnalytic,OptionValue[FullSimp]];
];
(*Simplify the expressions in the matrix in the case of an analytic calculation*)
If[isAnalytic,If[OptionValue[FullSimp],u2=FullSimplifyNoRoots[u2];
st=Join[FullSimplifyNoRoots[st1],st],u2=Simplify[u2];
st=Join[Simplify[st1],st]],st=Join[Simplify[st1],st]]
];
If[analyzeAnalyticCCDec,Print["Extract diagonal in ColumnByColumnDec..."]];
dia=ConjSimplify[Diagonal[u2]];
If[analyzeAnalyticCCDec,Print["Finished extracting diagonal in ColumnByColumnDec."]];
(*Check if it is necessary to decompose the diagonal gate, i.e., if it has length bigger than one and if it is not equal to identity up to a global phase shift.*)
If[Length[dia]>1 && dia/dia[[1]]!= ConstantArray[1,Length[dia]],
If[analyzeAnalyticCCDec,Print["Decompose diagonal in ColumnByColumnDec..."]];
st2=Reverse[DecDiagGate[dia]];
If[analyzeAnalyticCCDec,Print["Finished decomposing diagonal in ColumnByColumnDec with output st2=",st2]];
(*Applying the diagonal gate at the last Log(numCols) qubits*)
st2=RelabelQubits[st2,Range[Log[2,Length[dia]]],Range[Log[2,Length[dia]]]+bits-Log[2,Length[dia]]];
If[OptionValue[FullSimp],st=Join[FullSimplifyNoRoots[st2],st],st=Join[Simplify[st2],st]];
]];
(* Add indicator of qubits that start in 0: this will be run through reverse st then reverse befor output *)
For[i=1,i<=free,i++,st=Insert[st,{6,0,i},1]];
actionQ=
Switch[action,Null,Range[bits],_,action];
If[actionQ==Range[bits],,st=RelabelQubits[st,Range[bits],actionQ]];
If[OptionValue[Simp],
If[analyzeAnalyticCCDec,Print["Running SimplifyGateList on st=",FullSimplifyNoRoots[Reverse[InverseGateList[st]]]]];
If[OptionValue[FullSimp],SimplifyGateList[FullSimplifyNoRoots[Reverse[InverseGateList[st]]]],SimplifyGateList[Simplify[Reverse[InverseGateList[st]]],FullSimp->OptionValue[FullSimp]]]
,Reverse[InverseGateList[st]],
Throw[StringForm["Error: If condition checking 'OptionValue[Simp]' in ColumnByColumnDec did neither return True nor False"]]]
)]
]

(*---------------------------------------Column-by-column decomposition (private)------------------------------------*)

(*Decomposes a single column of an isometry from m to n\[GreaterEqual]2 qubits into single-qubit rotations and C-NOT gates. The input is the isometry isoTemp with already decomposed colums 1,...,colindex. The column colindex+1 is decomposed (note that colindex starts at zero). The number of bits is equal to n. If the computation should be done analytically, isAnalytic should be set to True.*)
DecSingleColumn[isoTemp_,colindex_, bits_,isAnalytic_,FullSim_:True] := Module[{ind1,ind2,iso,numCols,k,aSPlus1,bSPlus1,kBinary,bina,binb,j,l,s,k0,k1,vec,gate,st,st2,st1,fix,u2,control,target,angles,dia,mat,mcg, ucg,col,diag,start,bS},(
If[analyzeAnalyticCCDec,Print["Running DecSingleColumn with input {isoTemp,colindex, bits,isAnalytic}",{isoTemp,colindex, bits,isAnalytic}]];
iso=isoTemp;
col = iso[[All,colindex+1]];
st = {};
kBinary=IntegerDigits[colindex,2,bits];(*Binary representation of colindex*)
For[s=0,s<= bits-1,s++,
st1={};st2={};
Subscript[k, s]=kBinary[[bits-s]];
target=bits-s;
(*Notation: in "Quantum circuits for isometries" the variables aSPlus1,bSPlus1 and bS are denoted by Subsuperscript[a, s+1, k],Subsuperscript[b, s+1, k] and Subsuperscript[b, s, k], respectively.*)
If[s<bits-1,
aSPlus1=FromDigits[kBinary[[1;;bits-s-1]],2],
aSPlus1=0
];
bSPlus1=FromDigits[kBinary[[bits-s;;bits]],2];
If[s==0,
bS=0,
bS=FromDigits[kBinary[[bits-s+1;;bits]],2]
];
(*Check if a MCG is necessary and decompose it in the case that it is required*)
If[Subscript[k, s]==0&&bSPlus1!=0 &&Chop[col[[(2*aSPlus1+1)*2^s+bSPlus1+1]]]!=0,
If[analyzeAnalyticCCDec,Print["Decomposing MCG in DecSingleColumn..."]];
control=Flatten[Position[kBinary,1]];
(*Remove target qubit if it is included in control*)
control=Delete[control,Position[control,target]];
ind1=2*aSPlus1*2^s+bSPlus1;
ind2=(2*aSPlus1+1)*2^s+bSPlus1;
vec={{col[[ind1+1]]},{col[[ind2+1]]}};
If[analyzeAnalyticCCDec,Print["Calling Qubit Invert with input vec= ",vec]];
gate=QubitInvert[vec,0];
If[isAnalytic,If[FullSim,gate=FullSimplifyNoRoots[gate],gate=Simplify[gate]]];
If[analyzeAnalyticCCDec,Print["Calculated (and simplifid) the gate for the MCG: gate= ",gate]];
(*Decompose the MCG*)
st1=Reverse[DecMCSpecialUnitaryUpToDiagonal[gate,control,target,bits,{ReturnDiagonal->True,FullSimp->FullSim}]];
If[analyzeAnalyticCCDec,Print["Calculated decomposition of the MCG (using DecMCSpecialUnitaryUpToDiagonal). Outcome is: st1= ",st1]];
(*Update iso*)
iso=ApplyMCG[{gate,control,{},target,bits},iso];
If[st1[[1]][[1]]==-2,
diag=st1[[1]];st1=Drop[st1,1];
iso=ApplyDiag[{Conjugate[diag[[2]]],diag[[3]],bits},iso];
];
If[isAnalytic,iso=If[FullSim,FullSimplifyNoRoots[iso],Simplify[iso]]];
(*Update column vector*)
col=iso[[All,colindex+1]];
If[analyzeAnalyticCCDec,Print["Finished decomposing MCG in DecSingleColumn."]];
];
(*Create gates for UCG*)
control=Table[i,{i,1,target-1}];
If[isAnalytic,
gate=Table[IdentityMatrix[2],{2^(bits-s-1)}],
gate=Table[N[IdentityMatrix[2]],{2^(bits-s-1)}];
];
If[bSPlus1==0,start=aSPlus1,start=aSPlus1+1];
For[l=start,l<= 2^(bits-s-1)-1,l++, 
ind1=2*l*2^s+bS;
ind2=(2*l+1)*2^s+bS;
If[analyzeAnalyticCCDec,"Calculating gates for UCG"];
gate[[l+1]]=QubitInvert[{{col[[ind1+1]]},{col[[ind2+1]]}},Subscript[k, s]];
];
If[isAnalytic,gate=If[FullSim,FullSimplifyNoRoots[gate],Simplify[gate]]];
If[analyzeAnalyticCCDec,"Finished calculating gates for UCG with outcome gate= ",gate];
(*Decompose the uniformly controlled gate if it is not equal to the identity operator*)
If[IsUCGEqualToIdendity[gate],,
If[Length@control>0,
st2=Reverse[DecUCGUpToDiagonalHelp[gate,control,target,bits,FullSim]];
,
angles=ZYZDecomposition[gate[[1]]];
st2={{3,Chop@angles[[3]],target},{2,Chop@angles[[2]],target},{3,Chop@angles[[1]],target}};
];
(*Apply the UCG to the isometry*)
If[analyzeAnalyticCCDec,Print["Applying UCG in DecSingleColumn..."]];
iso=ApplyUCG[{gate,control,target,bits},iso];
If[analyzeAnalyticCCDec,Print["Finished applying UCG in DecSingleColumn."]];
(*Apply the diagonal gate to the isometry*)
If[st2[[1]][[1]]==-2,
diag=st2[[1]];
	   st2=Drop[st2,1];
	   If[analyzeAnalyticCCDec,Print["Applying diagona lgate in DecSingleColumn..."]];
iso=ApplyDiag[{Conjugate[diag[[2]]],diag[[3]],bits},iso];
If[analyzeAnalyticCCDec,Print["Finished applying diagonal gate in DecSingleColumn."]];
];
If[isAnalytic,iso=If[FullSim,FullSimplifyNoRoots[iso],Simplify[iso]]];
(*Updated col*)
col=iso[[All,colindex+1]];
];
st=Join[st2,st1,st];
];
If[analyzeAnalyticCCDec,Print["Finished running DecSingleColumn for column with index: ",colindex]];
{st,iso}
)]

(*---------------------------------------Decompositions on a small number of qubits (public)------------------------------------*)

(*State preparation on a single qubit*)
Options[StatePrep1Qubit]={IgnoreAncilla->False};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
StatePrep1Qubit[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{st,actionQ,qBits},
st=ZYZDec[IsoToUnitary[u],1];
If[Length[st]>0&&st[[1]][[1]]==3,st=Drop[st,1]];If[OptionValue[IgnoreAncilla],,st=Insert[st,{5,0,1},1]];
qBits=1;
actionQ=
Switch[action, 
Null,qBits, 
_, action
];
If[actionQ==qBits,,st=RelabelQubits[st,{qBits},{actionQ}]];
st
]

(*State preparation on two qubits based on [O. Giraud, M. Znidaric, Marko and B. Georgeot, Phys. Rev. A 80, 042309 (2009)].*)
Options[StatePrep2Qubits]={IgnoreAncilla->False};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
StatePrep2Qubits[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] :=  Module[{qBits,actionQ,a,\[CapitalSigma],b,schmidtParameter, localQubitUnitary1,localQubitUnitary2, st,a1,b1,c1,d1,a2,b2,c2,d2},
{a,\[CapitalSigma],b} = SingularValueDecomposition[ArrayReshape[u,{2,2}]];
schmidtParameter = ArcTan2[\[CapitalSigma][[1,1]], \[CapitalSigma][[2,2]]];
localQubitUnitary1 =a;(*qubitinvert[a[[All,1]],0];*)
localQubitUnitary2 =Conjugate[b];(*qubitinvert[b[[All,1]],0];*)
{a1,b1,c1,d1} = ZYZDecomposition[localQubitUnitary1];
{a2,b2,c2,d2} = ZYZDecomposition[localQubitUnitary2];
st = 
{
{3,Chop[c2],2},
{2,Chop[b2],2},
{3,Chop[a2],2},
{3,Chop[c1],1},
{2,Chop[b1],1},
{3,Chop[a1],1}
(*{0,1,2},
{2, Chop[-2*schmidtParameter],1}*)
};
If[Chop[N[\[CapitalSigma][[2,2]]]] != 0 ,
st =Join[st,
{{0,1,2},
{2, Chop[-2*schmidtParameter],1}}
]
];
;If[OptionValue[IgnoreAncilla],,st=Insert[Insert[st,{5,0,1},-1],{5,0,2},-1],
Throw[StringForm["Error: If condition checking 'OptionValue[IgnoreAncilla]' in DecMCSpecialUnitaryUpToDiagonal did neither return True nor False"]]
];
qBits=2;
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
Reverse[st]
];

(*Given an isometry from 1 to 2 qubits, returns a decomposition of a it into single-qubit and C-NOT gates. (cf. Appendix B1 in "Quantum circuits for isometries"*)
Options[DecIso12]={Simp->True,IgnoreAncilla->False,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
DecIso12[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] := Module[{actionQ,qBits,diag ,st,d, rotationAngle,v},
v=IsoToUnitary[u];
If[OptionValue[FullSimp],st=DecUnitary2Qubits[CTSimplify[v],{1,2},UpToDiagonal->True],st=DecUnitary2Qubits[CT[v],{1,2},UpToDiagonal->True]]; 
(*Decompose the inverse of v, so the diagonal ends up in the right place (this trick lets us save cnots by ignoring the action of the diagonal on the second qubit)*)
d=st[[-1]][[2]];
st=Drop[st,-1];
diag = DiagonalMatrix[d];
diag[[3,3]]  = 1;
diag[[4,4]] = 1;
rotationAngle =  Arg[diag[[1,1]] /diag[[2,2]]];
st = Append[st, {3,rotationAngle,2}];
st = InverseGateList[st];
qBits=2; If[OptionValue[IgnoreAncilla],,st=Insert[st,{5,0,1},1]];
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
If[OptionValue[Simp],SimplifyGateList[st,FullSimp->OptionValue[FullSimp]],st,
Throw[StringForm["Error: If condition checking 'OptionValue[IgnoreAncilla]' in DecIso12 neither returned True nor False"]]
]
]

(*--------------------- State preparation on three qubits (public)------------------------------------*)
(*
  See Optimal number of controlled-NOT gates to generate a three-qubit state 
PHYSICAL REVIEW A 77, 032320, 2008
https://journals.aps.org/pra/pdf/10.1103/PhysRevA.77.032320 
*)
Options[StatePrep3Qubits]={IgnoreAncilla->False,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
StatePrep3Qubits[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] :=Module[{st,actionQ,qBits},
st=StatePrep3QubitsReversed[u,OptionValue[FullSimp]];
qBits=3;If[OptionValue[IgnoreAncilla],,st=Join[st,{{5,0,1},{5,0,2},{5,0,3}}],
Throw[StringForm["Error: If condition checking 'OptionValue[IgnoreAncilla]' in StatePrep3Qubits did neither return True nor False"]]
];
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
Reverse[st]
]

(*--------------------- State preparation on three qubits (private)------------------------------------*)
(*
Following Generalized Schmidt decomposition and classification of three-quantum-bit states
A.Acin, A.Andrianov, L.Costa,E.Jane, J.I.Latorre and R.Tarrach
https://arxiv.org/pdf/quant-ph/0003050.pdf
*)
(*
Returns {U,a1,s1,b1,a2,s2,b2} such that 
specialQubitIndex \[Equal] 1 => u \[Equal] Flatten[KroneckerProduct[U[[1]], s1[[1,1]]a1[[1]],b1[[1]]]]+
Flatten[KroneckerProduct[U[[2], s2[[1,1]]KroneckerProduct[a2[[1]],b2[[1]]]+ s2[[2,2]]KroneckerProduct[a2[[2]],b2[[2]]]]]; 
 
specialQubitIndex \[Equal] 2 => u \[Equal] Flatten[KroneckerProduct[s1[[1,1]]a1[[1]],U[[1]] ,b1[[1]]]]+
Flatten[KroneckerProduct[s2[[1,1]]a2[[1]],U[[2]],b2[[1]]]] + Flatten[KroneckerProduct[s2[[2,2]]a2[[2]],U[[2]],b2[[2]]]];

specialQubitIndex \[Equal] 3 => u \[Equal] Flatten[KroneckerProduct[s1[[1,1]]a1[[1]],b1[[1]]],U[[1]]]]+
Flatten[KroneckerProduct[Flatten[s2[[1,1]]KroneckerProduct[a2[[1]],b2[[1]]] + s2[[2,2]]KroneckerProduct[a2[[2]],b2[[2]]] ],U[[2]]]];


*)
ThreeQubitSchmidtDecomposition[u_, specialQubitIndex_] := Module[{T1,T2,T1Prime, eVals,u11,u12,U,a1,s1,b1,a2,s2,b2,lambda,T2Prime,ans,T1Indices,T2Indices},
T1Indices = Select[Range[0,7], Function[n, IntegerDigits[n,2,3][[specialQubitIndex]] == 0]] +1;
T2Indices = Select[Range[0,7], Function[n, IntegerDigits[n,2,3][[specialQubitIndex]] == 1]] +1;
T1 = ArrayReshape[u[[T1Indices]],{2,2}];
T2 =  ArrayReshape[u[[T2Indices]],{2,2}];

eVals = Eigenvalues[{T2,T1}];

lambda = If[Chop[N[eVals[[1]]]] != 0, eVals[[1]], If[Chop[N[eVals[[2]]]] != 0, eVals[[2]], Throw["In threeQubitSchmidtDecomposition, neither generalised eigenValue is non-zero"]]];

If[N[eVals[[1]]]!=\[Infinity],u11 = Abs[eVals[[1]]]Sqrt[1/(1+Abs[eVals[[1]]]^2)];
u12 = -(u11/eVals[[1]]),u11 = Limit[Abs[x]Sqrt[1/(1+Abs[x]^2)],x->\[Infinity]];
u12 = Limit[-(u11/x),x->\[Infinity]]];If[Chop[N[u12]]!=0,
U = {{u11, u12}, {(1-u11^2)/u12,-u11}},U = {{u11, u12}, {Limit[(1-u11^2)/x,x->0],-u11}}];
T1Prime = U[[1,1 ]]T1 + U[[1,2]] T2;
{a1,s1,b1} = SingularValueDecomposition[T1Prime];
T2Prime = U[[2,1]] T1 + U[[2,2]] T2;
{a2,s2,b2} = SingularValueDecomposition[T2Prime];
a1 = Transpose[a1];
a2 = Transpose[a2];
b1 = CT[b1];
b2 = CT[b2];
U = Conjugate[U];

{U,a1,s1,b1,a2,s2,b2}
];

StatePrep3QubitsReversed[u_,FullSim_:True] := 
 Module[{a, sigma, b, qubit, unentangledQubitState, 
   entangledQubitsState, unentangledQubitSt, entangledQubitSt, 
   entangledQubitIndices, relableFn, secondSystemStateOne, 
   secondSystemStateTwo, a1, s1, b1, a2, s2, b2, U, u1, u2, u3, 
   highestNonSpecialQubit, otherQubit, gammaPrime, g1, g2, m1, m2, 
   relabelingList, w1, x1, y1, z1, w2, x2, y2, z2, w3, x3, y3, z3, 
   state, rho23, vals, vecs, i, theta, counter, a23, s23, b23, 
   schmidtParameter, localQubitUnitary1, localQubitUnitary2, vPrime1, 
   vPrime2},
  (*First see if we can separate out the each qubit - 
  this covers classes 0 & 1 from the paper*)
  (* see if the first qubit is in tensor product with the other two *)
{sigma, a, b} = SchmidtDecomposition[u, {2, 4}]; 
If[Chop[N[sigma[[1]]]] == 0,Print["StatePrep3Qubits: possible ordering issue in case 1 for input, ",u]];
  If[Chop[N[sigma[[2]]]] == 0,
   (* there is only one non-
   zero Schmidt coefficient so this state is a product state*)
   
   unentangledQubitState = a[[1]];
   entangledQubitsState = b[[1]]; 
   unentangledQubitSt = Reverse[StatePrep1Qubit[unentangledQubitState,,IgnoreAncilla->True]]; 
   entangledQubitIndices  = {2, 3}; 
   relableFn[listForm_] := If[listForm[[1]] > 0, 
     {listForm[[1]], listForm[[2]], 
      entangledQubitIndices[[listForm[[3]]]]}, 
     {listForm[[1]], entangledQubitIndices[[listForm[[2]]]], 
      entangledQubitIndices[[listForm[[3]]]]}]; 
   entangledQubitSt = 
    Map[relableFn, Reverse[StatePrep2Qubits[entangledQubitsState,,IgnoreAncilla->True]]]; 
   Return[Join[unentangledQubitSt, entangledQubitSt]]];
  
  (* see if the second qubit is in tensor product with the other two *)
  {sigma, a, b} = 
   SchmidtDecomposition[
    ExchangeSystems[u, {2, 1, 3}, {2, 2, 2}], {2, 4}]; 
  If[Chop[N[sigma[[1]]]] == 0,Print["StatePrep3Qubits: possible ordering issue in case 2 for input, ",u]];
If[Chop[N[sigma[[2]]]] == 0,
   (* there is only one non-
   zero Schmidt coefficient so this state is a product state*)
   
   unentangledQubitState = a[[1]];
   entangledQubitsState = b[[1]]; 
   unentangledQubitSt = 
    Map[{#[[1]], #[[2]], 2} &, 
     Reverse[StatePrep1Qubit[unentangledQubitState,,IgnoreAncilla->True]]]; 
   entangledQubitIndices = {1, 3}; 
   relableFn[listForm_] := 
    If[listForm[[1]] > 0, {listForm[[1]], listForm[[2]], 
      entangledQubitIndices[[listForm[[3]]]]}, {listForm[[1]], 
      entangledQubitIndices[[listForm[[2]]]], 
      entangledQubitIndices[[listForm[[3]]]]}]; 
   entangledQubitSt = 
    Map[relableFn, Reverse[StatePrep2Qubits[entangledQubitsState,,IgnoreAncilla->True]]]; 
   Return[Join[unentangledQubitSt, entangledQubitSt]]];
  
  (* see if the third qubit is in tensor product with the other two *)
  {sigma, a, b} = SchmidtDecomposition[u, {4, 2}]; 
  If[Chop[N[sigma[[1]]]] == 0,Print["StatePrep3Qubits: possible ordering issue in case 3 for input, ",u]];
If[Chop[N[sigma[[2]]]] == 0,
   (* there is only one non-
   zero Schmidt coefficient so this state is a product state*)
   
   unentangledQubitState = b[[1]]; entangledQubitsState = a[[1]]; 
   unentangledQubitSt = 
    Map[{#[[1]], #[[2]], 3} &, 
     Reverse[StatePrep1Qubit[unentangledQubitState,,IgnoreAncilla->True]]]; 
   entangledQubitSt = Reverse[StatePrep2Qubits[entangledQubitsState,,IgnoreAncilla->True]]; 
   Return[Join[unentangledQubitSt, entangledQubitSt]]];
  

(*If we get here the state is not of class 0 or 1*)
(*try class 2*)
For[qubit = 1, qubit <= 3, qubit++ ,
{U, a1,s1,b1,a2,s2,b2} = ThreeQubitSchmidtDecomposition[u, qubit];
If[Chop[N[s2[[2,2]]]] == 0,(*state now looks like |\[CapitalPsi]> = Cos[\[Theta]]|\[Alpha]\[Beta]\[Gamma]> + Sin[\[Theta]]|Subscript[\[Alpha], \[Perpendicular]]\[Beta]'\[Gamma]'>*)
(*First we find single qubit unitaries so  Subscript[u, 1]\[CircleTimes]Subscript[u, 2]\[CircleTimes]Subscript[u, 3]|\[CapitalPsi]> = Cos[\[Theta]]|000> + Sin[\[Theta]]|0> \[CircleTimes] Subscript[u, 2]|\[Beta]'> \[CircleTimes] Subscript[u, 3]|\[Gamma]'>*)
u1 = Conjugate[U];
u2 = QubitInvert[ArrayReshape[a1[[1]], {2, 1}],0];
u3= QubitInvert[ArrayReshape[b1[[1]], {2, 1}],0];

state = s1[[1,1]]Flatten[KroneckerProduct[U[[1]], a1[[1]],b1[[1]]]]+ s2[[1,1]]Flatten[KroneckerProduct[U[[2]],KroneckerProduct[a2[[1]],b2[[1]]]]]; 
(*Print[MatrixForm[s1[[1,1]]{1,0,0,0,0,0,0,0}  + s2[[1,1]]Flatten[ KroneckerProduct[{0,1}, u2.a2[[1]], u3.b2[[1]]]]]];*)
(*So now we prepareCos[\[Theta]]|000> + Sin[\[Theta]]|0> \[CircleTimes] Subscript[u, 2]|\[Beta]'> \[CircleTimes]Subscript[u, 3]|\[Gamma]'>, before applying the inverse of Subscript[u, 1]\[CircleTimes]Subscript[u, 2]\[CircleTimes]Subscript[u, 3]*)
highestNonSpecialQubit = If[qubit == 3, 2, 3];
otherQubit = If[qubit == 1, 2, 1];
gammaPrime = u3.(b2[[1]]);
{g1, g2} =Arg[gammaPrime];
{m1, m2} =Abs[gammaPrime];
theta = ArcTan2[m1,m2];
relableFn[listForm_] :=If[listForm[[1]] >0, 
{listForm[[1]], listForm[[2]], {qubit, otherQubit}[[listForm[[3]]]]}, 
{listForm[[1]], {qubit, otherQubit}[[listForm[[2]]]], {qubit, otherQubit}[[listForm[[3]]]]}];
{w1,x1,y1,z1} = ZYZDecomposition[CT[u1]];
{w2,x2,y2,z2} = ZYZDecomposition[CT[u2]];
{w3,x3,y3,z3} = ZYZDecomposition[CT[u3]];
Return[Catenate[{
{
{3, y1, qubit},{2, x1, qubit},{3, w1, qubit},
{3, y2, otherQubit},{2, x2, otherQubit},{3, w2, otherQubit},
{3, y3, highestNonSpecialQubit},{2, x3, highestNonSpecialQubit},{3, w3, highestNonSpecialQubit}
},
{
{3,g1 - g2,highestNonSpecialQubit},
{2, \[Pi]/2- theta,highestNonSpecialQubit},
{0,qubit, highestNonSpecialQubit},
{2, -(\[Pi]/2)+ theta,highestNonSpecialQubit}
},
Map[relableFn,Reverse[StatePrep2Qubits[{s1[[1,1]]Exp[-I (g1-g2)/2],0,0,0} + s2[[1,1]]Flatten[KroneckerProduct[{0,1},Exp[I (g1+g2)/2]u2.(a2[[1]])]],,IgnoreAncilla->True]]]
}
]
];
];
];
(*If we get here the state is not of class 0, 1 or 2*)
Return[StatePrep3QubitsGeneral[u,,FullSim]]
];

(*This will always require 3 CNOTs when doing state prep on 3 qubits -- not optimal in all cases  (see StatePrep3Qubits[] for an optinal method). 
The decomposition uses  Plesch and Brukner's decomposition scheme for state prep and then decomposes the last isometry using only two CNOTs.*)
StatePrep3QubitsGeneral[v_,action_:Null,FullSim_:True]:=Module[{qBits,actionQ,gates,isoGates,\[Alpha],a,b,mat1,mat2,anum,bnum,cnum,dnum},
      gates = {};
{\[Alpha],a,b} = SchmidtDecomposition[v, {2,4}];

(*find mats for last stages of Plesch and Bruckner state prep*)
     {mat1,mat2} = LastStageMat[a,b];
isoGates = Reverse[DecIso12[mat2,,{IgnoreAncilla->True,FullSimp->FullSim}]];
gates = Join[gates,RelabelQubits[isoGates,{1,2},{2,3}]];
{anum,bnum,cnum,dnum} = ZYZDecomposition[mat1];
      gates = Join[gates,{{3,cnum,1},{2,bnum,1},{3,anum,1}}];
AppendTo[gates,{0,1,3}];

{anum,bnum,cnum,dnum} = ZYZDecomposition[IsoToUnitary[Transpose[{\[Alpha]}]]];
gates = Join[gates,{{3,cnum,1},{2,bnum,1},{3,anum,1}}];
qBits=3;
actionQ=
Switch[action, 
Null,Range[qBits], 
_, Sort[action]
];
If[actionQ==Range[qBits],,gates=RelabelQubits[gates,Range[qBits],actionQ]];
gates
]

(*--------------------- Decomposition for isometries [best decomposition scheme is chosen] (public)------------------------------------*)

(*Takes an arbitrary isometry u as an input and chooses the best method (in terms of the required number of C-NOT gates) to decompose it into single-qubit rotations and C-NOT gates*)
Options[DecIsometry]={FullSimp->True,SpeedUp->False};
DecIsometry[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] := Module[{UseDecString,actionQ,qBits,dim,m,n,out,out1,free,i},
IsQubitIsometry[u,"DecIsometry"];
If[OptionValue[SpeedUp],UseDecString="QSD",UseDecString="DecIsometry",Throw[StringForm["Error: 'OptionValue[SpeedUp]' in DecIsometry is neither True nor False"]]];
If[Length[Dimensions[u]]==1,If[u=={},Return[{}],dim={1,Dimensions[u][[1]]}],dim=Reverse[Dimensions[u]]];
If[Chop[N[u]]==IdentityMatrix[dim[[2]]][[All,1;;dim[[1]]]],
out={};free=Log2[dim[[2]]/dim[[1]]];
 For[i=1,i<=free,i++,out=Insert[out,{5,0,i},1]],
Switch[
Map[Log2,dim],
{1,2},out=DecIso12[u,action,FullSimp->OptionValue[FullSimp]];out1=QSD[u,action,FullSimp->OptionValue[FullSimp]];If[CNOTCount[out]>=CNOTCount[out1],out=out1];
out1=KnillDec[u,action,{FullSimp->OptionValue[FullSimp],UseDec->UseDecString}];If[CNOTCount[out]>=CNOTCount[out1],out=out1];
out1=ColumnByColumnDec[u,action,{FirstColumn->"StatePreparation",FullSimp->OptionValue[FullSimp]}];If[CNOTCount[out]>=CNOTCount[out1],out=out1];,
{0,_},out=StatePreparation[u,action,FullSimp->OptionValue[FullSimp]];
out1=ColumnByColumnDec[u,action,{FullSimp->OptionValue[FullSimp]}];If[CNOTCount[out]>=CNOTCount[out1],out=out1]
,
{_,_},{m,n}=Map[Log2,dim];
out=QSD[u,action,FullSimp->OptionValue[FullSimp]];
out1=KnillDec[u,action,{FullSimp->OptionValue[FullSimp],UseDec->UseDecString}];If[CNOTCount[out]>=CNOTCount[out1],out=out1];out1=ColumnByColumnDec[u,action,{FirstColumn->"StatePreparation",FullSimp->OptionValue[FullSimp]}];If[CNOTCount[out]>=CNOTCount[out1],out=out1];
out1=ColumnByColumnDec[u,action,FullSimp->OptionValue[FullSimp]];If[CNOTCount[out]>=CNOTCount[out1],out=out1];
];
];
qBits=Log2[Dimensions[u][[1]]];
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,out=RelabelQubits[out,Range[qBits],actionQ]];
out
]

(*Takes an arbitrary isometry u as an input and chooses the best method (in terms of the required number of C-NOT gates for a generic gate of the given dimensions) to decompose it into single-qubit rotations and C-NOT gates*)
Options[DecIsometryGeneric] = {Simp->False,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
DecIsometryGeneric[u_,action:Except[_?OptionQ]:Null,OptionsPattern[]] := Module[{actionQ,qBits,dim,m,n,out},
IsQubitIsometry[u,"DecIsometryGeneric"];
If[Length[Dimensions[u]]==1, If[u=={},Return[{}],dim={1,Dimensions[u][[1]]}],dim=Reverse[Dimensions[u]]];
If[Chop[N[u]]==IdentityMatrix[dim[[2]]][[All,1;;dim[[1]]]],
out={};,
Switch[
Map[Log2 , dim],
{1,2},out=DecIso12[u,action,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
{0,_},out=StatePreparation[u,action,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
{_,_},{m,n}=Map[Log2 , dim] ;
If[m==n-1 ||m==n,
out=QSD[u,action,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}];,
out=ColumnByColumnDec[u,action,{FirstColumn->"StatePreparation",Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}];
];
];
];
qBits=Log2[Dimensions[u][[1]]];
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,out=RelabelQubits[out,Range[qBits],actionQ]];
out
]

(*----------------------- Methods for state preparation(public)------------------*)
Options[StatePreparation]={level->3,Simp->True,IgnoreAncilla->False,FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
StatePreparation[v_,action:Except[_?OptionQ]:Null,OptionsPattern[]]:=Module[{actionQ,qBits,st,i},
IsQubitIsometry[v,"StatePreparation"];
qBits=Log2[Length[v]];
st=StatePrepRecursive[v,OptionValue[level],FullSimp->OptionValue[FullSimp]];
If[OptionValue[IgnoreAncilla],,st=Join[st,Table[{5,0,i},{i,1,qBits}]],
Throw[StringForm["Error: If condition checking 'OptionValue[IgnoreAncilla]' in StatePreparation did neither return True nor False"]]
];
actionQ=
Switch[action, 
Null,Range[qBits], 
_, action
];
If[actionQ==Range[qBits],,st=RelabelQubits[st,Range[qBits],actionQ]];
If[OptionValue[Simp],SimplifyGateList[Reverse[st],FullSimp->OptionValue[FullSimp]],Reverse[st],
Throw[StringForm["Error: If condition checking 'OptionValue[Simp]' in StatePreparation did neither return True nor False"]]
]
]

(*----------------------- Methods for state preparation(private)------------------*)(*helper function - gives the list representation of CNOTs for stage 2 of state prep (as described by Plesch and Brukner) - i.e. a chain of CNOTs with targets running from firstQubit to the lastQubit ***USES NOTATION THAT 1 IS MOST SIGNIFICANT (TOP) AND N IS LEAST (LOWEST)*** *)
ChainCNOT[firstQubit_,lastQubit_,firstTargetQubit_] := Module[{n,st},
n = firstTargetQubit - firstQubit; (*this is the 'distance' of each qubit*)
Table[CNOT[i,i+n],{i,firstQubit,lastQubit}]
]

(*This method gives the gates/matrices required in stage 3 and 4,given that the schmidt decomposition vectors have been inputed -- need two functions: one for outputing the matrix representation of the isometries and one giving the gate output coming from the csd. *)
LastStageMat[a_,b_] := Module[{u1mat, u2mat,i},
(*create matrices and extend to unitary if neccesary*)
u1mat = a[[1]]; u2mat = b[[1]];
For[i=2,i<=Dimensions[a][[1]],i++,u1mat = Join[u1mat,a[[i]],2]];
For[i=2,i<=Dimensions[b][[1]],i++,u2mat = Join[u2mat,b[[i]],2]];
Return[{u1mat,u2mat}]
]

(*This is an implementation of the Plesch and Bruckner state preparation scheme (optimized in "Quantum Circuits for Isometries"). Note: due to the ordering of the first columns, e.g. |011> = |0>|11> rather than |110> = |11>|0>, the ancilla qubit for state preparation on an odd number of qubits is at the top of the second "section".*)
Options[StatePrepRecursive] = {FullUnitary-> False,DoDecomposeUnitaries->True,FullSimp->True}; (*If the Option FullUnitary is set to True, the isometry from (n-1)/2 to (n+1)/2 arising fo state preparation on an odd number n of qubits is extended to a full unitary. This option is required for the Knill decomposition.*)
(*ToDo: add option for chosing the decomposition scheme to decompose the isometries*)
(*ToDo: A diagonal gate can be extracted from the unitary appearing in step 3 and can be merged into state preparation in step 1 (this will save one C-NOT).*)
StatePrepRecursive[v_,level_,OptionsPattern[]]:=Module[{st2Qubit,numQubits,currentState,\[Alpha],a,b,mat1,mat2,gates,gates1,gates2,u,qubitNumA,qubitNumB,twoQubitGateIndices,diagEntries,diagEntriesInv,twoQubitGate,stTwoQubitUpToDiag,stTwoQubitUpToDiagInv},
numQubits = Log[2,Dimensions[v][[1]]]; 
If[IntegerQ[numQubits],
currentState = v;
gates = {};
While[numQubits>= level+1,
qubitNumA=IntegerPart[numQubits/2];
If[EvenQ[numQubits],
qubitNumB=(numQubits/2),
qubitNumB=IntegerPart[numQubits/2]+1
];
{\[Alpha],a,b} = SchmidtDecomposition[currentState,{2^qubitNumA,2^qubitNumB}];
{mat1,mat2} = LastStageMat[a,b];
If[OptionValue[DoDecomposeUnitaries],,Return[{{mat1,IsoToUnitary[mat2]}, 
Join[Reverse[StatePrepRecursive[Transpose[{\[Alpha]}],level,FullSimp->OptionValue[FullSimp]]],ChainCNOT[1,qubitNumA,qubitNumB+1]]}],
Throw[StringForm["Error: If condition checking 'OptionValue[DoDecomposeUnitaries]' in StatePrepRecursive did neither return True nor False"]]
];

If[OptionValue[FullUnitary]&&OddQ[numQubits],mat2=IsoToUnitary[mat2],,
Throw[StringForm["Error: If condition checking 'OptionValue[FullUnitary]' in StatePrepRecursive did neither return True nor False"]]
];
(*ANCILLA QUBIT IS AT THE TOP OF THE SECOND SECTION!!!*)

(*gates are reversed in sections so the overall structure is preserved but multiply gates works with gatelist *)
currentState = Transpose[{\[Alpha]}];
If[qubitNumA<  2 || OptionValue[FullUnitary],
gates1 = QSD[mat1,,IgnoreAncilla->True],
gates1 = InverseGateList[Reverse[QSD[CT[mat1],FullSimp->OptionValue[FullSimp]]]];
gates1=Reverse[gates1];
(*Decopose the inverse such that the circuit starts with a two-qubit gate that is decomposed into three C-NOTS (and not into 2 (up to diagonal))*)
twoQubitGateIndices=ExtractTwoQubitGatesIndices[gates1,qubitNumA-1,qubitNumA](*Extract the Indices of the two-qubit gate first applied in stage 3 (and coming from the QSD). The goal is to extract a two qubit diagonal and merge it into the state preparation in stage 1 (this allows to save one C-NOT gate)*);
twoQubitGate=CreateIsometryFromList[RelabelQubits[gates1[[twoQubitGateIndices]],{qubitNumA-1,qubitNumA},{1,2}],2,FullSimp->OptionValue[FullSimp]];
(*Take inverse of two qubit gate, such that the diagonal appears on the total left of stage 3 in the circuit*)
st2Qubit=DecUnitary2Qubits[ConjugateTranspose[twoQubitGate],{qubitNumA-1,qubitNumA},{UpToDiagonal->True,FullSimp->OptionValue[FullSimp]}];
{stTwoQubitUpToDiagInv,diagEntriesInv}={Drop[st2Qubit,-1],st2Qubit[[-1]][[2]]};
stTwoQubitUpToDiag=InverseGateList[stTwoQubitUpToDiagInv];
twoQubitGateIndices = Map[{#}&,twoQubitGateIndices]; (*The syntax for passing multiple indices to "Delete" requires this *);
gates1=Delete[gates1,twoQubitGateIndices];
gates1=Join[stTwoQubitUpToDiag,gates1];
(*Merge diagonal gate into currentState*)
diagEntries=Conjugate[diagEntriesInv];
currentState=ApplyDiag[{diagEntries,{qubitNumA-1,qubitNumA},qubitNumA},currentState];
];
gates1=DeleteCases[gates1,{x_/;x>=5,_,_}];
gates2 = RelabelQubits[QSD[mat2,,{IgnoreAncilla->True,FullSimp->OptionValue[FullSimp]}],Range[qubitNumB],Range[qubitNumA+1,numQubits]];
(*find the required gates*)
gates = Join[gates2,gates];
gates = Join[gates1,gates];
gates = Join[ChainCNOT[1,qubitNumA,qubitNumB+1],gates];
(*ancilla qubit at top leads to an additional shift in the odd case*)
numQubits=IntegerPart[numQubits/2];
];
(*gates for last state prep -- this depends on the level of recursion*)
If[numQubits == 1,gates=Join[StatePrep1Qubit[currentState,,IgnoreAncilla->True],gates]];
If[numQubits == 2,gates=Join[StatePrep2Qubits[currentState,,IgnoreAncilla->True],gates]];
If[numQubits == 3,gates=Join[StatePrep3Qubits[currentState,,{IgnoreAncilla->True,FullSimp->OptionValue[FullSimp]}],gates]];
Return[Reverse[gates]], 
Print["StatePrepRecursive: Input given does not have the right dimensions to be a state"]
]
];

(*----------------------- Methods for Knill's decomposition (public)------------------*)

(* For inputs x and y (n-by-m matrices) which satsify x^\[Dagger]x = y^\[Dagger]y, outputs u (n-by-n) such that ux = y - based on Lemma A.1 from Knill 1995: http://arxiv.org/pdf/quant-ph/9508006. Should output error message if x^\[Dagger]x \[NotEqual] y^\[Dagger]y or if the matrices have different dimension*)
XToYTransform[x_,y_]  := Module[{u,w,\[CapitalSigma],v,svd,a,b,c,n,m,\[Sigma],\[Sigma]inv,i},

If[Chop[N[CT[x].x  -CT[y].y]] == 0*IdentityMatrix[Dimensions[y][[2]]] && Dimensions[x] == Dimensions[y],
(*code given that the condition x^\[Dagger]x = y^\[Dagger]y holds*)

(*find dimensions*)
n = Dimensions[x][[1]]; m=Dimensions[x][[2]];

{w,\[CapitalSigma],v} = SingularValueDecomposition[y]; v = CT[v];

a = x.CT[v]; (*n-by-m*)
(*find diagonal elts - should be in vector form*)
\[Sigma] = Diagonal[\[CapitalSigma]]; \[Sigma]inv = {}; For[i = 1, i <= Dimensions[\[Sigma]][[1]], i++,
 If[Chop[N[\[Sigma][[i]]]] == 0, \[Sigma]inv =
   Insert[\[Sigma]inv, 0, -1], \[Sigma]inv =
   Insert[\[Sigma]inv, \[Sigma][[i]]^(-1), -1]]]; 

(*forming b*)
b = CT[a]; (*m-by-n*)
For[i = 1, i <= m, i++,b[[i]] = \[Sigma]inv[[i]]*b[[i]]];
c = CT[IsoToUnitary[CT[b]]];
u = w.c;
Return[u],
Return[Print["XToYTransform: The condition \!\(\*SuperscriptBox[\(x\), \(\[Dagger]\)]\)x = \!\(\*SuperscriptBox[\(y\), \(\[Dagger]\)]\)y does not hold or the matrices are of different dimension"]],
Return[Print["XToYTransform: The condition \!\(\*SuperscriptBox[\(x\), \(\[Dagger]\)]\)x = \!\(\*SuperscriptBox[\(y\), \(\[Dagger]\)]\)y does not hold or the matrices are of different dimension - are they even matrices?"]]
]
]

(*input v must be a n-by-m isometry with m\[LessEqual]n. The isometry v is then extended to a unitary maximizing the numbers of eigenvalues that are equal to one*)
IsoToUnitarySpecial[v_] := Module[{n,m,vTran,w,u,x,x1,x2, p,q,z},
(*find dimensions of v*)
{n,m}=Dimensions[v];
(*do nothing and return input if input is already unitary -- this is to stop nullspace being {} casuing errors*)
If[n==m,Return[v],(*continue*)]; 

(*Define W which is another isometry and extends v to a unitary*)
vTran = Transpose[v]; 
w = CT[Orthogonalize[NullSpace[vTran]]][[All,1;;n-m]];
w = NormaliseCols[w];  

(*Define u = [v^\[Dagger]/w^\[Dagger]] (with w^\[Dagger] underneath) *)
u = Join[CT[v],CT[w]];

(*Find an x such that ConjugateTranspose[v].x has the first m cols of x*)
x = GiveX[v];  x1 = x[[1;;m]]; x2 = x[[m+1;; n]];
p = CT[w].x;
q = XToYTransform[p,x2];
z = CT[Join[CT[v],q.CT[w]]];
Return[z]
]

(*UnitaryEigenvalueDecomp decomposes a unitary matrix as in the last part of Theorem 3.1 from: http://arxiv.org/pdf/quant-ph/9508006 - i.e. into state preparation matrices and matrices with {e^Subscript[i\[Theta], k],1,1,...,1} as diagonal elements.*)
UnitaryEigenvalueDecomp[u_] := Module[{angles,statePrepMats,nonIdEvals,nonIdEvecs,esys,i,n},

n = Dimensions[u][[1]];

(*take esys as eigensystem of input*)
esys = EigensystemExact[u]; 

(*initialise matrices*)
nonIdEvals ={};nonIdEvecs = {};angles = {}; statePrepMats = {};

(*find non identity eigenvalues and corresponding eigenvectors*)
For[i = 1, i<=n, i++, If[Chop[N[esys[[1]][[i]]-1]]==0,(*if identity eigenvalue, do nothing*) ,
AppendTo[nonIdEvals,esys[[1]][[i]]]]];
For[i = 1, i<=n, i++, If[Chop[N[esys[[1]][[i]]-1]]==0,(*if identity eigenvalue, do nothing*) ,
AppendTo[nonIdEvecs,Transpose[{esys[[2]][[i]]}]]]];

(*find angles given by non-identity eigenvalues*)
angles = Arg[nonIdEvals];

Return[{angles,nonIdEvecs}]
]

(*
Takes in an 2^n-by-2^m (m\[LessEqual]n) isometry and outputs a list of gates. See Quantum Circuits for Isometries Section IV B.*)
(*ToDo [Optimization]: If we work on a odd number of quits, we could implement the first and the last unitary 
(that comes from Plesch and Brukner's state preparation scheme) as an isometry*)
Options[KnillDec] = {UseDec->"DecIsometry",Simp->True,FirstColumn->"UCG",FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
KnillDec[v_, action:Except[_?OptionQ]:Null,OptionsPattern[]] := 
 Module[{angles, u, vectors, gates1, mats1, mats2, gates2, merge1, 
   merge2, len, free,gates, f, mat, i, num, n,actionQ},
  IsQubitIsometry[v,"KnillDec"];
  Switch[OptionValue[UseDec],
  "QSD", f[mat_] := QSD[mat,,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
  "KnillDec",f[mat_] := KnillDec[mat,,{UseDec->OptionValue[UseDec],Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
  "ColumnByColumnDec",f[mat_] := ColumnByColumnDec[mat,,{Simp->OptionValue[Simp],FirstColumn->OptionValue[FirstColumn],FullSimp->OptionValue[FullSimp]}],
  "DecIsometry", f[mat_] := DecIsometry[mat,FullSimp->OptionValue[FullSimp]],
  "DecIsometryGeneric", f[mat_] := DecIsometryGeneric[mat,FullSimp->OptionValue[FullSimp]],
  _,Throw[StringForm["Error: The option value for 'UseDec' in KnillDec was not recognized."]]
  ];
  
  u = Simplify[IsoToUnitarySpecial[v]];
  n = Log2[Dimensions[u][[1]]];
  
  (*Note: KnillDecomp does not work on a 1 qubit isometry, 
  so pass to QSD*)
  
  If[n == 1, gates=QSD[u,action,FullSimp->OptionValue[FullSimp]],
   {angles, vectors} = UnitaryEigenvalueDecomp[u];
   len = Dimensions[angles][[1]]; 
   gates = {};
   (*neccesary since there is a possibility that there are more than \
n-m identity eigenvalues*)
   
   If[len >= 2,
    (*first section (unmerged) *)
    
    gates = 
     Join[gates, 
      Reverse[InverseGateList[StatePrepRecursive[vectors[[len]], 1,{FullUnitary->True,FullSimp->OptionValue[FullSimp]}]]]];
    gates = Join[gates, MCGAngleDecomp[angles[[len]], n]];
    
    (*add all the remaining gates except for the very last state \
preparation
    ---START FROM LEN AND WORK BACKWARDS!!! (matrices representation \
reverse of gate representation!!!)*)
    {mats1, gates1} = 
     StatePrepRecursive[vectors[[len]], 1,{DoDecomposeUnitaries->False,FullSimp->OptionValue[FullSimp]}]; {mats2, gates2} = 
     StatePrepRecursive[vectors[[len - 1]], 1,{DoDecomposeUnitaries->False,FullSimp->OptionValue[FullSimp]}];
    
    For[i = len, i >= 2, i--,
     
     (*append gates from v_(len-i+1) before merged mats*)
     
     gates = Join[gates, gates1];
     
     (*find merged mats*)
     merge1 = CT[mats2[[1]]].mats1[[1]];
     merge2 = CT[mats2[[2]]].mats1[[2]];
     
     (*add gates from first merged matrix*)
     
     gates = Join[gates, f[merge1]];
     (*add relabeled gates from second merged matrix*)
     
     If[EvenQ[n],
      gates = 
       Join[gates, 
        RelabelQubits[f[merge2], Range[n/2] , Range[(n/2) + 1, n]]],
      gates = 
       Join[gates, 
        RelabelQubits[f[merge2], Range[(n + 1)/2] , 
         Range[(n + 1)/2, n]]]
      ];
     
     (*append gates from (conjugate transpose of) v_(len-
     i) after merged mats*)
     
     gates = Join[gates, InverseGateList[gates2]];
     
     (*angles*)
     
     gates = Join[gates, MCGAngleDecomp[angles[[i - 1]], n]];
     If[i >= 3, {mats1, gates1} = {mats2, gates2};
      {mats2, gates2} = StatePrepRecursive[vectors[[i - 2]], 1,{DoDecomposeUnitaries->False,FullSimp->OptionValue[FullSimp]}],
      (*do nothing*)];
     
     ];
    
    (*Append gates for V_1 SP*)
    
    gates = Join[gates, 
      Reverse[StatePrepRecursive[vectors[[1]], 1,{FullUnitary->True,FullSimp->OptionValue[FullSimp]}]]],
    
    If[len == 1, 
     gates = Join[gates, 
       Reverse[InverseGateList[StatePrepRecursive[vectors[[len]], 1,{FullUnitary->True,FullSimp->OptionValue[FullSimp]}]]]];
     gates = Join[gates, MCGAngleDecomp[angles[[len]], n]]; 
     gates = Join[gates, Reverse[StatePrepRecursive[vectors[[1]], 1,{FullUnitary->True,FullSimp->OptionValue[FullSimp]}]]]]
    ]];
(* Add indicator of qubits that start in 0 *)
 free=Log2[Dimensions[v][[1]]/Dimensions[v][[2]]];
 For[i=1,i<=free,i++,gates=Insert[gates,{5,0,i},1]];
   actionQ=
Switch[action, 
Null,Range[n], 
_, action
];
If[actionQ==Range[n],,gates=RelabelQubits[gates,Range[n],actionQ]];
   
If[OptionValue[Simp],SimplifyGateList[gates,FullSimp->OptionValue[FullSimp]],gates,
Throw[StringForm["Error: If condition checking 'OptionValue[Simp]' in KnillDec did neither return True nor False"]]
]
   ]
  
  
  (*----------------------- Methods for Knill's decomposition (private)------------------*)

(*For an n-by-m matrix v, where n\[LessEqual]m, gives x such that v^\[Dagger].x has first m rows equal to x*)
GiveX[vMat_] := Module[{nDim,mDim,idMbyN,v2,xMat,x1,x2},
nDim = Dimensions[vMat][[1]];mDim = Dimensions[vMat][[2]];
idMbyN = SparseArray[{i_,i_}->1,{mDim,nDim}]; 
v2 = CT[vMat] - idMbyN; 
xMat = Transpose[NullSpace[v2]][[All,1;;nDim-mDim]]; 
xMat = NormaliseCols[xMat];
Return[xMat]
]

(*MCGAngleDecomp is used to decompose the controlled "angle" gates used in the Knill decomposition*)
(*ToDo: Adapt notation to the "standard" one for controlled gates and make MCGAngleDecomp[] public*)
MCGAngleDecomp[angle_,numQubits_]:=Module[{gatelist,diag},
(*diagonal gate decomp beats multi-controlled, for numQubits \[LessEqual] 10*)
If[numQubits <= 10, 
diag = {Exp[I*angle]};
Do[AppendTo[diag,1],{i,2^(numQubits)-1}];
gatelist = DecDiagGate[diag],

gatelist = DecMCG[DiagonalMatrix[{1,Exp[I*angle]}],Range[numQubits-1],numQubits,numQubits];
For[i=numQubits, i>=1,i--,PrependTo[gatelist,{1,Pi,i}]];  
For[i=1,i<=numQubits,i++,AppendTo[gatelist,{1,Pi,i}]]; 
];
gatelist
]

(*----------------------- Methods for Channels and POVMs (public)------------------*)

(* Constructs the Choi state from a list of Kraus operators *)
KrausToChoi[set_] := 
 Module[{ent}, 
  ent = Sum[
     KroneckerProduct[KetV[i - 1, Dimensions[set[[1]]][[2]]], 
      KetV[i - 1, Dimensions[set[[1]]][[2]]]], {i, 1, 
      Dimensions[set[[1]]][[2]]}];ent=ent.CT[ent]; 
  Sum[KroneckerProduct[IdentityMatrix[Dimensions[set[[1]]][[2]]], 
     set[[k]]].ent.KroneckerProduct[
     IdentityMatrix[Dimensions[set[[1]]][[2]]], CT[set[[k]]]], {k, 1, 
    Dimensions[set][[1]]}]]

(* Constructs a set of Kraus operators from a Choi state. Note that this 
works more generally for arbitrary linear operators; if given a Choi state the
outputs set1 and set2 should be equal (may be useful for checking) *)
ChoiToKraus[state_, dA_, dB_] := 
 Module[{set1 = {}, set2 = {}, i, j, w1, w2, u, d, v}, {u, d, v} = 
   SingularValueDecomposition[state]; 
  For[i = 1, i <= dA*dB, i++, If[Chop[N[d[[i,i]]]]!=0,
   w1 = (d[[i, i]])^(1/2)*Transpose[{Transpose[u][[i]]}]; 
   w2 = (d[[i, i]])^(1/2)*Transpose[{Transpose[v][[i]]}]; 
   set1 = Insert[set1, Transpose[Partition[Flatten[w1], dB]], -1]; 
   set2 = Insert[set2, Transpose[Partition[Flatten[w2], dB]], -1]]]; 
  {set1, set2}]

(* Takes a channel and finds another representation of it potentially with fewer Kraus ops *)
MinimizeKrausRank[chan_] := 
 Module[{n, dA, dB}, {n, dB, dA} = Dimensions[chan]; 
  ChoiToKraus[KrausToChoi[chan], dA, dB][[1]]]

(* Takes in a channel (specified in terms of a list of Kraus \
operators) and returns (u,a) where u is an isometry and a is the \
number of qubit ancilla that need to be traced out to recover the \
original channel *)
Options[StinespringQubit] = {TryToCompress->True};
StinespringQubit[chan_, OptionsPattern[]] := 
 Module[{out, n, anc, chan2}, 
  If[OptionValue[TryToCompress], chan2 = MinimizeKrausRank[chan], 
   chan2 = chan,
   Throw[StringForm["Error: If condition checking 'OptionValue[TryToCompress]' in StinespringQubit did neither return True nor False"]]
   ]; n = Dimensions[chan2][[1]]; 
  anc = Ceiling[Log[2, n]]; {Sum[
    KroneckerProduct[Transpose[{UnitVector[2^anc, i]}], 
     chan2[[i]]], {i, 1, n}], anc}]

(* Takes in a POVM (specified in terms of a list of POVM \
elements) and returns (u,a) where u is an isometry and a is the \
number of qubit ancilla that need to be measured to recover the \
POVM *)
POVMToIsometry[POVM_] := 
 Module[{out, n, anc}, n = Dimensions[POVM][[1]]; 
  anc = Ceiling[Log[2, n]]; {Sum[
    KroneckerProduct[Transpose[{UnitVector[2^anc, i]}], 
     MatrixPower[POVM[[i]], 1/2]], {i, 1, n}], anc}]
     
(*----------------------- Decomposition of channels and POVMs (public)------------------*)

(* Takes in a channel (specified in terms of a list of Kraus \
operators) and returns a gate sequence (including tracing out operations \
 at the end of the circuit) in list format that implements the channel.*)
Options[DecChannelInQCM] = {TryToCompress->True,DecomposeIso->"DecIsometry",Simp->False,FirstColumn->"UCG",FullSimp->True};
(*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
DecChannelInQCM[chan_,actionAndAncilla:Except[_?OptionQ]:Null,OptionsPattern[]] := 
 Module[{anc,st,st1,v,a,n,ancillaQ,actionQ,aQ}, 
  {v,a}=StinespringQubit[chan, TryToCompress->OptionValue[TryToCompress]];
  Switch[OptionValue[DecomposeIso],
  "DecIsometry",st=DecIsometry[v,FullSimp->OptionValue[FullSimp]],
  "DecIsometryGeneric",st=DecIsometryGeneric[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
  "QSD",st=QSD[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
  "ColumnByColumnDec",st=ColumnByColumnDec[v,Null,{Simp->OptionValue[Simp],FirstColumn->OptionValue[FirstColumn],FullSimp->OptionValue[FullSimp]}],
  "KnillDec",st=KnillDec[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
    _,Throw[StringForm["An unknown decomposition was provided as as option value for DecomposeIso."]]
    ];
  st1=Table[{4,0,anc},{anc,Range[a]}];
  st=Join[st,st1];
  n=Log2[Dimensions[v][[1]]];
  ancillaQ=Switch[actionAndAncilla, 
Null,Range[a], 
_, aQ=actionAndAncilla[[1]];If[Length[aQ]>a,aQ=Delete[aQ,Transpose[{Range[a+1,Length[aQ]]}]]];aQ
];
  actionQ=Switch[actionAndAncilla, 
Null,a+Range[n-a], 
_, actionAndAncilla[[2]]
];
Switch[actionAndAncilla,
Null,(*no relabeling required*),
_,st=RelabelQubits[st,Range[n],Join[ancillaQ,actionQ]]
];
st
];

  (* Takes in a POVM (specified in terms of a list of POVM \
elements) and returns a gate sequence (including measurments) \
 in list format that implements the channel. *)
 Options[DecPOVMInQCM] = {DecomposeIso->"DecIsometry",Simp->False,FirstColumn->"UCG",FullSimp->True,PostMmt->False};
 (*Except[_?OptionQ] is a trick to allow for optional arguments (together with options). Without this trick, having something like f[x_,y:Null,OptionsPattern[]]:=...
would give an error calling f[x,option\[Rule]optionValue]*)
DecPOVMInQCM[POVM_,actionAndAncilla:Except[_?OptionQ]:Null,OptionsPattern[]] := 
 Module[{anc,st,st1,v,a,n,ancillaQ,actionQ,aQ,st2}, 
  {v,a}=POVMToIsometry[POVM];
  Switch[OptionValue[DecomposeIso],
  "DecIsometry",st=DecIsometry[v,FullSimp->OptionValue[FullSimp]],
  "DecIsometryGeneric",st=DecIsometryGeneric[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
  "QSD",st=QSD[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
   "ColumnByColumnDec",st=ColumnByColumnDec[v,Null,{Simp->OptionValue[Simp],FirstColumn->OptionValue[FirstColumn],FullSimp->OptionValue[FullSimp]}],
    "KnillDec",st=KnillDec[v,Null,{Simp->OptionValue[Simp],FullSimp->OptionValue[FullSimp]}],
    _,Throw[StringForm["An unknown decomposition was provided as as option value for DecomposeIso."]]
    ];
  st1=Table[{4,1,anc},{anc,Range[a]}];
  If[OptionValue[PostMmt],
  st2={},
  st2=Table[{4,0,anc},{anc,Range[a+1,Log2[Dimensions[v][[1]]]]}]
  ];
  st=Join[st,st1,st2];
  n=Log2[Dimensions[v][[1]]];
  ancillaQ=Switch[actionAndAncilla, 
Null,Range[a], 
_, aQ=actionAndAncilla[[1]];If[Length[aQ]>a,aQ=Delete[aQ,Transpose[{Range[a+1,Length[aQ]]}]]];aQ
];
  actionQ=Switch[actionAndAncilla, 
Null,a+Range[n-a], 
_, actionAndAncilla[[2]]
];
Switch[actionAndAncilla,
Null,(*no relabeling required*),
_,st=RelabelQubits[st,Range[n],Join[ancillaQ,actionQ]]
];
  st
  ]
  
  (*----------------------- Methods for checking the inputs of methods (private)------------------*)
  
IsQubitIsometry[v_,methodName_:"UNKNOWN"]:=Module[{numRow,numCol},
  numRow=Dimensions[v][[1]];
  numCol=Dimensions[v][[2]];
  If[IntegerQ[Log2[numRow]],,Throw[StringJoin["The number of rows of the input isometry in the method ",methodName ," is not a power of two."]]];
  If[IntegerQ[Log2[numCol]],,Throw[StringJoin["The number of columns of the input isometry in the method ",methodName ," is not a power of two."]]];
   If[numCol>numRow,Throw[StringJoin["The number of columns of the input matrix in the method ",methodName ," is bigger than the number of rows (i.e., it is not an isometry)."]]];
  If[isIdentity[ConjugateTranspose[N[v]].N[v]],,Throw[StringJoin["The input matrix in the method ",methodName ," is not an isometry since ConjugateTranspose[v].v]\[NotEqual]Id."]]]; 
  ]
  
IsListFormHelp[gate_,methodName_]:=Module[{},
  If[MemberQ[{-2,-1,0,1,2,3,4,5,6,100,101},gate[[1]]],,Throw[StringJoin["The gate ",ToString[gate]," appearing as an input in method ",methodName ," is of unknown type."]]];
  Which[
   MemberQ[{-2},gate[[1]]],If[gate[[2]]=={}&&gate[[3]]=={},Goto[LabelEnd];];If[Length[Dimensions[gate[[2]]]]==1&&Length[Dimensions[gate[[3]]]]==1,,Throw[StringJoin["The dimensions of the parameter lists of the diagonal gate ",ToString[gate],"  appearing in the input of method ",methodName," are not supported."]]];
   If[Log2[Length[gate[[2]]]]==Length[gate[[3]]],,Throw[ToString[StringForm["The diagonal gate `1` appearing as an input in method `2` does not contain 2^`3` entries.",gate,methodName,Length[gate[[3]]]]]]],
    MemberQ[{-1,0},gate[[1]]],If[IntegerQ[gate[[2]]]&&IntegerQ[gate[[3]]],,Throw[StringJoin["There is a control gate ",ToString[gate]," appearing as an input in method ",methodName ," that has a control or a target qubit number that is not an integer."]]];
    If[gate[[2]]==gate[[3]],Throw[StringJoin["There is a controlled gate ",ToString[gate]," appearing as an input in method ",methodName ," that has a control qubit number equal to the target qubit number (which is not supported)."]]],
  MemberQ[{1,2,3},gate[[1]]], If[NumericQ[gate[[2]]]&&IntegerQ[gate[[3]]],,Throw[StringJoin["There is a rotation gate ",ToString[gate]," appering as an input in method ",methodName ," that has incorrect parameters."]]],
  MemberQ[{4,5,6},gate[[1]]],If[MemberQ[{0,1},gate[[2]]]&&IntegerQ[gate[[3]]],,Throw[ToString[StringForm["The gate `1` appearing as an input in method `2` has unknown parameter types.",gate,methodName]]]],
 MemberQ[{100},gate[[1]]],If[Length[gate[[3]]]==2&&IntegerQ[gate[[3]][[1]]]&&IntegerQ[gate[[3]][[2]]]&&NumericQ[gate[[2]]],,Throw[StringJoin["There is an XX gate ",ToString[gate]," appering as an input in method ",methodName ," that has incorrect parameters."]]],
MemberQ[{101},gate[[1]]],If[Length[gate[[2]]]==2&&NumericQ[gate[[2]][[1]]]&&NumericQ[gate[[2]][[2]]]&&IntegerQ[gate[[3]]],,Throw[StringJoin["There is an R gate ",ToString[gate]," appering as an input in method ",methodName ," that has incorrect parameters."]]]
];
Label[LabelEnd];
  ];
  
IsListForm[st_,methodName_:"UNKNOWN"]:=Module[{postsel,postselnums,traceout,traceoutnums,measure,measurenums,int1,int2,int3},If[st=={},,If[Length[Dimensions[st]]==1, IsListFormHelp[st,methodName],postsel=Cases[st,{6,_,_}];If[postsel==={},postselnums={},postselnums=Transpose[postsel][[3]]];traceout=Cases[st,{4,0,_}];If[traceout==={},traceoutnums={},traceoutnums=Transpose[traceout][[3]]];measure=Cases[st,{4,1,_}];If[measure==={},measurenums={},measurenums=Transpose[measure][[3]]];int1=Intersection[postselnums,traceoutnums];int2=Intersection[postselnums,measurenums];int3=Intersection[traceoutnums,measurenums];If[int1==={},,Throw[StringJoin["For qubits ",ToString[int1]," appearing as an input in method ",methodName ," there is both a postselection and a trace."]]];If[int2==={},,Throw[StringJoin["For qubits ",ToString[int2]," appearing as an input in method ",methodName ," there is both a postselection and a measurement."]]];If[int3==={},,Throw[StringJoin["For qubits ",ToString[int3]," appearing as an input in method ",methodName ," there is both a trace and a measurement."]]];IsListFormHelp[#,methodName]&/@st]]]
  
PrepareForQASM[st_]:=Module[{out,traceouts,traceoutnums,i},out=NGateList[st];If[Cases[out,{5,_,_}]==={},,Print["Notice: ancillas found in the input have been removed.  Output gate sequence corresponds to the same operation as the input provided the ancillas in the input sequence start in the correct states."];out=DeleteCases[out,{5,_,_}]];If[Cases[out,{-2,_,_}]==={},,Throw["Error in PrepareForQASM: diagonal gates found in the input, which is not supported by QASM. You may want to decompose them using DecDiagGate[]."];out=DeleteCases[out,{-2,_,_}]];If[Cases[st,{6,_,_}]==={},,Print["Warning: postselection gate found in the input has been removed.  Output gate sequence may not be as intended."];out=DeleteCases[out,{6,_,_}]];traceouts=Cases[st,{4,0,_}];If[traceouts==={},,traceoutnums=Transpose[traceouts][[3]];Print["Notice: trace out gate found in the input has been replaced by measurement.  Forgetting the outcome will recover the same operation as the input."];out=DeleteCases[out,{4,0,_}];For[i=1,i<=Length[traceoutnums],i++,out=Insert[out,{4,1,traceoutnums[[i]]},-1]]];out] 

Options[PickRandomCircuitIsometry]={TotGates->False}
PickRandomCircuitIsometry[qubitsin_,qubitsout_,totcnots_,OptionsPattern[]]:=Module[{i,out,type,cnots,ctrl,targ}, 
  If[qubitsin>qubitsout,Throw[StringForm["PickRandomCircuitIsometry Error: number of input qubits must be smaller than number of output qubits."]]];
  out=Table[{5,0,j},{j,1,qubitsout-qubitsin}];
  If[OptionValue[TotGates],
  For[i=1,i<=totcnots,i++,ctrl=0;targ=0;type=RandomInteger[{0, 3}];
  If[type == 0,While[ctrl==targ,ctrl=RandomInteger[{1,qubitsout}];targ=RandomInteger[{1,qubitsout}]];
     out=Insert[out,{0,ctrl,targ},-1], out=Insert[out,{type,2*\[Pi]*RandomReal[],RandomInteger[{1,qubitsout}]},-1]]]
     , 
  out=Join[out,Table[{3,2*\[Pi]*RandomReal[],j},{j,qubitsout-qubitsin+1,qubitsout}]];
  out=Join[out,Table[{2,2*\[Pi]*RandomReal[],j},{j,1,qubitsout}]]; 
  out=Join[out,Table[{3,2*\[Pi]*RandomReal[],j},{j,1,qubitsout}]]; 
  For[i=1,i<=totcnots,i++,ctrl=0;targ=0; 
    While[ctrl==targ,ctrl=RandomInteger[{1,qubitsout}];targ=RandomInteger[{1,qubitsout}]];
    out=Insert[out,{0,ctrl,targ},-1]; 
    out=Insert[out,{2,2*\[Pi]*RandomReal[],ctrl},-1]; 
    out=Insert[out,{3,2*\[Pi]*RandomReal[],ctrl},-1]; 
    out=Insert[out,{2,2*\[Pi]*RandomReal[],targ},-1]; 
    out=Insert[out,{1,2*\[Pi]*RandomReal[],targ},-1]]];
out]
   
(* Commands for conversion to XX and R *)

(* outputs (a,b,c,d) such that U is equal to Rx[a] followed by R[b,c] up to the phase E^(I*d) *)
 RxRGateDecomp[U_]:=Module[{a,b,c,d,th,phi},
 {a,b,c,d}=XYXDecomposition[U];
 th=2*ArcCos[Cos[c]*Cos[b/2]];
 If[Chop[(Cos[c]*Cos[b/2])^2-1]==0,phi=-\[Pi]/2,
 If[Chop[c]>0,phi=ArcSin[Sin[b/2]/(1-(Cos[c]*Cos[b/2])^2)^(1/2)]-\[Pi],
 If[Chop[c]<0,phi=-ArcSin[Sin[b/2]/(1-(Cos[c]*Cos[b/2])^2)^(1/2)],
 If[Chop[c]==0,phi=-\[Pi]/2]]]];
 {a-c,th,phi,d}]

 (* Replaces all CNOT gates with XX gates and additional single qubit rotations *)
 ReplaceCNOTWithXX[st_]:=Module[{out=st,i,pos,cnotposns,ctrl,targ},
 cnotposns=Flatten[Position[st,{0,_,_}]];
 For[i=Length[cnotposns],i>=1,i--,pos=cnotposns[[i]];ctrl=out[[pos]][[2]];
 targ=out[[pos]][[3]];out=Delete[out,pos];out=Insert[out,Ry[\[Pi]/2,ctrl],pos];
 out=Insert[out,Rx[\[Pi]/2,ctrl],pos];out=Insert[out,Rx[\[Pi]/2,targ],pos];
 out=Insert[out,XX[\[Pi]/4,ctrl,targ],pos];out=Insert[out,Ry[-\[Pi]/2,ctrl],pos]];
 out]
 
 (* Replaces all XX gates with CNOTs and additional single qubit rotations *)
ReplaceXXWithCNOT[st_]:=Module[{out=st,i,j,pos,xxposns,gates},
xxposns=Flatten[Position[st,{100,_,_}]];
For[i=Length[xxposns],i>=1,i--,pos=xxposns[[i]];
gates=DecUnitary2Qubits[XXM[out[[pos]][[2]]],out[[pos]][[3]]];
out=Delete[out,pos];For[j=Length[gates],j>=1,j--,out=Insert[out,gates[[j]],pos]]];
out]  
      
  (* Replaces single qubit rotations after XX gates with X rotations before the XX and an R-gate after *)
ReplaceRotationsWithRGatesAfterXX[st_]:=Module[{out=st,i,j,pos,xxposns,ctrl,targ,toremove,ctrlst,targst,ctrlfin,targfin,a,b,c,d,U},
xxposns=Flatten[Position[st,{100,_,_}]];
For[i=Length[xxposns],i>=1,i--,pos=xxposns[[i]];ctrl=out[[pos]][[3]][[1]];
targ=out[[pos]][[3]][[2]];toremove={};ctrlst={};targst={};
ctrlfin=False;targfin=False;
For[j=pos+1,j<=Length[out],j++,If[ctrlfin&&targfin,Break[]];
If[(CheckGateForQubit[out[[j]],ctrl]&&Not[ctrlfin])||(CheckGateForQubit[out[[j]],targ]&&Not[targfin]),
If[MemberQ[{1,2,3},out[[j]][[1]]],toremove=Insert[toremove,{j},-1];
If[out[[j]][[3]]==ctrl,ctrlst=Insert[ctrlst,out[[j]],-1],If[out[[j]][[3]]==targ,targst=Insert[targst,out[[j]],-1],Print["ReplaceWithRAfterXX: Error"]]],If[CheckGateForQubit[out[[j]],ctrl],ctrlfin=True];
If[CheckGateForQubit[out[[j]],targ],targfin=True]]]];
(* The next line could presumably be uncommented *)
(* ctrlst=MergeSameRot[ctrlst];targst=MergeSameRot[targst];*)
out=Delete[out,toremove];If[Length[ctrlst]==1&&ctrlst[[1,1]]==1,out=Insert[out,ctrlst[[1]],pos];pos++,
ctrlst=Transpose[Insert[Drop[Transpose[ctrlst],-1],Table[1,{i,1,Length[ctrlst]}],-1]];
U=CreateOperationFromGateList[ctrlst];{a,b,c,d}=RxRGateDecomp[U];
If[Chop[b]!=0,out=Insert[out,RGate[b,c,ctrl],pos+1]];
out=Insert[out,Rx[a,ctrl],pos];pos++];If[Length[targst]==1&&targst[[1,1]]==1,out=Insert[out,targst[[1]],pos];pos++,
targst=Transpose[Insert[Drop[Transpose[targst],-1],Table[1,{i,1,Length[targst]}],-1]];
U=CreateOperationFromGateList[targst];{a,b,c,d}=RxRGateDecomp[U];
If[Chop[b]!=0,out=Insert[out,RGate[b,c,targ],pos+1]];out=Insert[out,Rx[a,targ],pos];pos++]];
out]
   
(* Replaces all RGates with single qubit rotations using XYXDecomposition *)
ReplaceRGatesWithRotations[st_]:=Module[{out=st,i,gate,pos,Rposns,a,b,c,d,U},
Rposns=Flatten[Position[out,{101,_,_}]];
For[i=Length[Rposns],i>=1,i--,pos=Rposns[[i]];gate=out[[pos]];
U=CreateOperationFromGateList[{Insert[Drop[gate,-1],1,-1]}];
{a,b,c,d}=XYXDecomposition[U];out=Delete[out,pos];
If[Chop[c]!=0,out=Insert[out,{1,c,gate[[3]]},pos]];
If[Chop[b]!=0,out=Insert[out,{2,b,gate[[3]]},pos]];
If[Chop[a]!=0,out=Insert[out,{1,a,gate[[3]]},pos]]];
out]     
      
(* checks whether the given gate acts on qubit n *)
CheckGateForQubit[gate_,n_]:=If[(gate[[1]]==0&&(gate[[2]]==n||gate[[3]]==n))||(gate[[1]]==1&&gate[[3]]==n)||(gate[[1]]==2&&gate[[3]]==n)||(gate[[1]]==3&&gate[[3]]==n)||(gate[[1]]==-1&&(gate[[2]]==n||gate[[3]]==n))||(gate[[1]]==4&&gate[[3]]==n)||(gate[[1]]==6&&gate[[3]]==n)||(gate[[1]]==-2&&MemberQ[gate[[3]],n])||(gate[[1]]==100&&MemberQ[gate[[3]],n])||(gate[[1]]==101&&gate[[3]]==n),True,False]      
    
(* Replaces initial rotations by RGates where initial means before any non-rotation gate NB: this passes through initial ancilla *)
ReplaceInitialRotationsByRGates[st_]:=Module[{out=st,i,j,pos,xxposns,ctrl,targ,toremove,ctrlst,targst,ctrlfin,targfin,a,b,c,d,U,ancillain},
ancillain=SortBy[Cases[st,{5,_,_}],Last];out=DeleteCases[out,{5,_,_}];pos=1;
For[ctrl=NumberOfQubits[st],ctrl>=1,ctrl--,toremove={};ctrlst={};ctrlfin=False;
For[j=pos,j<=Length[out],j++,If[ctrlfin,Break[]];If[CheckGateForQubit[out[[j]],ctrl]&&Not[ctrlfin],
If[MemberQ[{1,2,3},out[[j]][[1]]],
toremove=Insert[toremove,{j},-1];ctrlst=Insert[ctrlst,out[[j]],-1],ctrlfin=True]]];
(* ctrlst=MergeSameRot[ctrlst];targst=MergeSameRot[targst];*)
out=Delete[out,toremove];If[ctrlst==={},,If[Length[ctrlst]==1&&ctrlst[[1,1]]==1,If[Chop[ctrlst[[1]][[2]]]!=0,out=Insert[out,RGate[-ctrlst[[1]][[2]],0,ctrl],pos]],ctrlst=Transpose[Insert[Drop[Transpose[ctrlst],-1],Table[1,{i,1,Length[ctrlst]}],-1]];
U=CreateOperationFromGateList[ctrlst];{a,b,c,d}=RxRGateDecomp[U];
If[Chop[b]!=0,out=Insert[out,RGate[b,c,ctrl],pos]];
If[Chop[a]!=0,out=Insert[out,RGate[-a,0,ctrl],pos]]]]];
Join[ancillain,out]]      
          
CNOTRotationsToXXRGates[st_]:=Module[{out},out=ReplaceCNOTWithXX[st];out=ReplaceRotationsWithRGatesAfterXX[out];ReplaceInitialRotationsByRGates[out]]

XXRGatesToCNOTRotations[st_]:=Module[{out},out=ReplaceXXWithCNOT[st];out=ReplaceRGatesWithRotations[out];SimplifyGateList[out]]          
          
End[];

EndPackage[]
