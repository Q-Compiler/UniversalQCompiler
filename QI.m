(* ::Package:: *)

(*
   Copyright 2019 (https://github.com/rogercolbeck/QI)

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


BeginPackage["QI`"]

CT::usage="CT[] is short for ConjugateTranspose[]."

KetV::usage="KetV[i,d] gives |i> in dimension d."

BraV::usage="BraV[i,d] gives <i| in dimension d."

DM::usage="DM[] converts a vector (state) into a (density) matrix."

CircleTimes::usage="CircleTimes[] performs the tensor product."

Tensor::usage="Tensor[A, B, pos, desc] takes the tensor product of A and B placing the systems in the positions pos and with the sizes desc."

TensorPower::usage="TensorPower[M,power] takes the tensor product of M with itself power times."

ExchangeSystems::usage="ExchangeSystems[A,newpos,desc] takes an operator representing systems in one order and converts it to another, as specified by newpos."

DirectSum::usage="DirectSum[{a, b, ...}] takes the direct sum of square matrices a, b, ..."

QubitPartialTrace::usage="QubitPartialTrace[M, {i, j, ...}] traces out the i, j, ... qubits of M assuming that M comprises multiple qubit spaces."

PartialTrace::usage="PartialTrace[M, dim1, dim2, tr] traces out the tr=1 or tr=2 system of M, where M represents a bipartite system of dimensions dim1 and dim2."

PT::usage="PT[M, keep, desc] traces out the systems specified by keep of M which is composed of subsystems of sizes desc. Components of keep equal to 1 are kept and those equal to 0 are traced out."

BasisForm::usage="BasisForm[vec, desc] writes out a vector vec based on dimensions desc."

BasisFormS::usage="BasisFormS[vec, desc] writes out a vector vec based on dimensions desc."

Purify::usage="Purify[rho] gives a state that purifies rho (in vector form). The chosen purification takes a copy of the original eigenvector."

SchmidtDecomposition::usage="SchmidtDecomposition[v, desc] gives the Schmidt decomposition of v where desc is a two component list of the dimensions across which the decomposition is taken."

DiagonalizingUnitary::usage="DiagonalizingUnitary[M] outputs {U, Diag} where U is the unitary that diagonalizes normal matrix M to Diag."

EigenvaluesExact::usage="EigenvaluesExact[M,(prec)] outputs as Eigenvalues, except that it orders the outputs first by real and then imaginary part treating these as equal if within prec of each other (default that of Chop[])."

EigensystemExact::usage="EigensystemExact[M,(prec)] outputs as Eigensystem, except that it orders the outputs first by real and then imaginary part treating these as equal if within prec of each other (default that of Chop[]) and it ensures orthonormality of the eigenvectors even for exact inputs."

SimultaneouslyDiagonalize::usage="SimultaneouslyDiagonalize[A, B] give a unitary U that diagonalizes both commuting normal matrices A and B."

BlochSphere::usage="BlochSphere[rho] gives the Bloch sphere representation of the qubit density operator rho."

FromBlochSphere::usage="FromBlochSphere[{r1, r2, r3}] gives the density matrix for the corresponding point on the Bloch sphere."

AppendCols::usage="AppendCols[] takes an isometry and adds columns to form a unitary."

FillZero::usage="FillZero[] takes a matrix with more columns than rows and returns a square matrix filling in zeros."

MeasureBasis::usage="MeasureBasis[rho, basis, sys, desc] performs the quantum channel on rho corresponding to measuring in the basis specified by basis on the sys part of rho and forgetting the outcome (keeping the system), partitioned as desc. Without sys and desc arguments, the first part of rho is measured."

MeasurePOVM::usage="MeasurePOVM[rho, POVM, sys, desc] performs the quantum channel on rho equivalent to measuring the POVM on the sys part of rho and then tracing out the measured system. sys is a list of 2s with one 1 corresponding to the system to measure."

POVMIsometry::usage="POVMIsometry[rho, POVM, sys, desc] performs the isometry on the sys part of rho equivalent to the Stinespring dilation of the measurement channel."

ChoiState::usage="ChoiState[channel] computes the Choi state for the channel (given as a list of Kraus operators). ChoiState[set1, set2] does the same for two sets of operators defining a linear map."

ChoiChannel::usage="ChoiChannel[state, dA, dB] computes the channel corresponding to the given Choi state, outputting {set1, set2}. If the state given is positive, set1=set2."

ChannelCompress::usage="ChannelCompress[channel] takes a Kraus representation of a channel and returns another Kraus representation of the channel that may have fewer Kraus operators."

ExtremeChannelQ::usage="ExtremeChannelQ[channel] checks whether the channel is extremal or not."

PickRandomPsi::usage="PickRandomPsi[n] chooses a random pure state in dimension n."

RPickRandomPsi::usage="RPickRandomPsi[n] chooses a random real pure state in dimension n."

FPickRandomPsi::usage="FPickRandomPsi[n, prec] chooses a random real pure state with exact values in dimension n."

PickRandomPsi2::usage="PickRandomPsi2[n] chooses a random pure state in dimension n."

RPickRandomPsi2::usage="RPickRandomPsi2[n] chooses a random real pure state in dimension n."

FPickRandomPsi2::usage="FPickRandomPsi2[n, prec] chooses a random real pure state with exact values in dimension n."

PickRandomUnitary::usage="PickRandomUnitary[n] chooses a random unitary in dimension n."

RPickRandomUnitary::usage="RPickRandomUnitary[n] chooses a random real unitary in dimension n."

FPickRandomUnitary::usage="FPickRandomUnitary[n, prec] chooses a random real unitary with exact values in dimension n."

PickRandomIsometry::usage="PickRandomIsometry[dim1, dim2] chooses a random isometry from dimension dim1 to dim2."

RPickRandomIsometry::usage="RPickRandomIsometry[dim1, dim2] chooses a random real isometry from dimension dim1 to dim2."

FPickRandomIsometry::usage="FPickRandomIsometry[dim1, dim2, prec] chooses a random real isometry with exact values from dimension dim1 to dim2."

PickRandomRho::usage="PickRandomRho[n] chooses a random density operator in dimension n."

RPickRandomRho::usage="RPickRandomRho[n] chooses a random real density operator in dimension n."

FPickRandomRho::usage="FPickRandomRho[n, prec] chooses a random real density operator with exact values in dimension n."

PickRandomPOVM::usage="PickRandomPOVM[dim, num] chooses a random POVM in dimension dim with num elements."

RPickRandomPOVM::usage="RPickRandomPOVM[dim, num] chooses a random real POVM in dimension dim with num elements."

PickRandomMeasurement::usage="PickRandomMeasurement[dim, num] chooses a random measurement in dimension dim with num elements."

RPickRandomMeasurement::usage="RPickRandomMeasurement[dim, num] chooses a random real measurement in dimension dim with num elements."

PickRandomChannel::usage="PickRandomChannel[dim1, dim2, num] chooses a random channel from dimension dim1 to dim2 with num Kraus operators."

RPickRandomChannel::usage="RPickRandomChannel[dim1, dim2, num] chooses a random real channel from dimension dim1 to dim2 with num Kraus operators."

FPickRandomChannel::usage="RPickRandomChannel[dim1, dim2, num, prec] chooses a random real channel with exact values from dimension dim1 to dim2 with num Kraus operators."

Dist::usage="Dist[rho, sigma] computes the trace distance between the two density matrices."

QuickDist::usage="QuickDist[rho, sigma] computes the trace distance between the two density matrices."

Fidelity::usage="Fidelity[rho, sigma] computes the fidelity between the two density matrices."

OptimumPOVM::usage="OptimumPOVM[rho, sigma (,p1)] returns the POVM that optimally distinguishes the states, where p1 is the a priori probability of rho (with only two arguments, p1 is taken as 1/2"

Matrixlog::usage="Matrixlog[(b,) M] takes the logarithm of Hermitian matrix M base b.  With only one argument this gives the natural logarithm."

Matrixxlogx::usage="Matrixxlogx[(b,) M] takes the function x*log(x) to Hermitian matrix M, taking 0*log(0) to be 0.  The optional first argument selects the base."

ShanEntropy::usage="ShanEntropy[{p1, p2, ...}] gives the Shannon entropy of the given distribution. ShanEntropy[p] gives the binary entropy of p."

vNEntropy::usage="vNEntropy[rho, keep, desc] computes the conditional von Neumann entropy of rho, where keep is a list of 0s, 1s and 2s representing systems to ignore (0), take the entropy of (1) or condition on (2), and desc give the dimensions of the subsystems."

vNInfo::usage="vNInfo[rho, keep, desc] computes the conditional mutual information of rho, where keep is a list of 0s, 1s, 2s and 3s representing systems to ignore (0), take the information between (1 and 2) or condition on (3), and desc give the dimensions of the subsystems."

RelEnt::usage="RelEnt[(alpha,) A, B] computes the relative entropy of A and B.  Using alpha gives the Renyi relative entropy."

RenyiEnt::usage="RenyiEnt[alpha, rho, keep, desc] computes the conditional Renyi entropy of rho, where keep is a list of 0s, 1s and 2s representing systems to ignore (0), take the entropy of (1) or condition on (2), and desc give the dimensions of the subsystems."

RemoveIneqConstraints::usage="RemoveIneqConstraints[M, b] removes redundancy from a set of equations of the form Mx >= b."

NRemoveIneqConstraints::usage="NRemoveIneqConstraints[M, b] removes redundancy from a set of equations of the form Mx >= b, performing the calculation numerically."

RemoveDuplicateConstraints::usage=="RemoveDuplicateConstraints[M, b] removes duplicate constraints from a set of equations of the form Mx >= b."

Prep::usage="Prep[M, b] converts inequality constraints to all be greater than form."

FourierMotzkin::usage="FourierMotzkin[M, b, elim] performs Fourier-Motzkin elimination to a set of equations in the form Mx >= b removing the variables specified in the list elim."

Elim::usage="Elim[M, b, i] removes the ith variable from the set of equations Mx >= b."

IntDigs::usage="IntDigs[num, bases] writes the number num in terms of digits of the specified bases, with output length the same as the number of bases."

Progress::usage="Progress[i, imin, imax, num] gives num progress indications where i is in a loop from imin to imax."

ProgressTemporary::usage="ProgressTemporart[i, imin, imax, num] gives num progress indications where i is in a loop from imin to imax, removing after the evaluation is completed."

ThreadSolve::usage="ThreadSolve[eqns, soln, var] solves the set of equations eqns[[i]]=soln for variables var."

NThreadSolve::usage="NThreadSolve[eqns, soln, var] numerically solves the set of equations eqns[[i]]=soln for variables var."

\[Sigma]::usage="\[Sigma][i] gives the ith Pauli matrix."

Support::usage="Support[M] returns the projector onto the support of matrix M."

PosPart::usage="PosPart[M (,assum)] returns the strictly positive part of matrix M with optional assumptions assum (useful for symbolic matrices."

NegPart::usage="NegPart[M (,assum)] returns the strictly negative part of matrix M with optional assumptions assum (useful for symbolic matrices."

AntiDiagonalMatrix::usage="AntiDiagonalMatrix[M] turns a matrix into a vector of its diagonal values."

CreateMatrix::usage="CreateMatrix[W, r, c] generates a rxc matrix with elements W[i, j]."

CreateHermitianMatrix::usage="CreateHermitianMatrix[W, r, c] generates a rxc Hermitian matrix with elements W[i, j] or Conjugate[W[i, j]]."

CreateHermitianMatrixR::usage="CreateHermitianMatrixR[W, r, c] generates a rxc Hermitian matrix with elements WR[i, j] +/- I WC[i, j]."

CreateSymmetricMatrix::usage="CreateSymmetricMatrix[W, r, c] generates a rxc symmetric matrix with elements W[i, j]."

RemoveRepeatedInstances::usage="RemoveRepeatedInstances[list] removes any repeated elements of list, preserving the order."

SetMinus::usage="SetMinus[S1, S2] removes elements of S2 from S1, maintaining the order."

MajorizeQ::usage="MajorizeQ[A, B] checks whether A majorizes B, where A and B are lists of numbers."

ManualDeriv::usage="ManualDeriv[fn, vars, point, steps] estimates the derivative of a function fn of the variables vars at point point, with stepsizes steps."

QRDecomp::usage="QRDecomp[M] gives {Q, R} such that Q is unitary, R is right (upper) triangular and CT[Q].R = M."

RQDecomp::usage="RQDecomp[M] gives {R, Q} such that R is right (upper) triangular, Q is unitary and R.CT[Q] = M."

QLDecomp::usage="QLDecomp[M] gives {Q, L} such that Q is unitary, L is left (lower) triangular and CT[Q].L = M."

LQDecomp::usage="LQDecomp[M] gives {L, Q} such that L is left (lower) triangular, Q is unitary and L.CT[Q] = M."

QRDecompPos::usage="QRDecompPos[M] outputs as QRDecomp but ensures R has positive diagonal entries."

QRDecompNeg::usage="QRDecompNeg[M] outputs as QRDecomp but ensures R has negative diagonal entries."

RQDecompPos::usage="RQDecompPos[M] outputs as RQDecomp but ensures R has positive diagonal entries."

RQDecompNeg::usage="RQDecompNeg[M] outputs as RQDecomp but ensures R has negative diagonal entries."

QLDecompPos::usage="QLDecompPos[M] outputs as QLDecomp but ensures L has positive diagonal entries."

QLDecompNeg::usage="QLDecompNeg[M] outputs as QLDecomp but ensures L has negative diagonal entries."

LQDecompPos::usage="LQDecompPos[M] outputs as RQDecomp but ensures L has positive diagonal entries."

LQDecompNeg::usage="LQDecompNeg[M] outputs as RQDecomp but ensures L has negative diagonal entries."






Begin["`Private`"]

CT[M_]:=ConjugateTranspose[M]

KetV[i_,d_]:=Module[{},If[i>=d||i<0,Print["KetV: input should be between 0 and d-1"]];Transpose[{UnitVector[d,i+1]}]]

BraV[i_,d_]:=Module[{},If[i>=d||i<0,Print["BraV: input should be between 0 and d-1"]];Conjugate[{UnitVector[d,i+1]}]]

DM[vec_]:=vec.CT[vec]

CircleTimes[a__]:=KroneckerProduct[a]

Tensor[rho_,sigma_,sys_,desc_]:=Module[{desc1,desc2},If[Length[Dimensions[desc]]==1,If[Dimensions[rho][[1]]==Dimensions[rho][[2]]&&Dimensions[sigma][[1]]==Dimensions[sigma][[2]],desc1=desc;desc2=desc1,If[Dimensions[rho][[2]]==1&&Dimensions[sigma][[2]]==1,desc1=desc;desc2=ConstantArray[1,Length[desc]],Print["Tensor: Fourth argument needs to be pairs of dimensions for given input."]]],{desc1,desc2}=Transpose[desc]];If[Dimensions[rho][[1]]!=Tr[DeleteCases[(2-sys)*desc1,0],Times]||Dimensions[sigma][[1]]!=Tr[DeleteCases[(sys-1)*desc1,0],Times]||Dimensions[rho][[2]]!=Tr[DeleteCases[(2-sys)*desc2,0],Times]||Dimensions[sigma][[2]]!=Tr[DeleteCases[(sys-1)*desc2,0],Times],Print["Tensor: inputs have incorrect dimensions"]];ExchangeSystems[rho\[CircleTimes]sigma,Flatten[Join[Position[sys,1],Position[sys,2]]],Transpose[{desc1[[Flatten[Join[Position[sys,1],Position[sys,2]]]]],desc2[[Flatten[Join[Position[sys,1],Position[sys,2]]]]]}]]]

TensorPower[M_,power_]:=Module[{i,R},For[i=1;R={{1}},i<=power,i++,R=CircleTimes[R,M]];R]

ExchangeSystems[vec_,newpos_,desc_]:=Module[{list1,list2,i,desc1,desc2},If[Length[Dimensions[desc]]==1,If[Dimensions[vec][[1]]==Dimensions[vec][[2]],desc1=desc;desc2=desc1,If[Dimensions[vec][[2]]==1,desc1=desc;desc2=ConstantArray[1,Length[desc]],Print["ExchangeSystems: Third argument needs to be pairs of dimensions for given input."]]],{desc1,desc2}=Transpose[desc]];If[Tr[desc1,Times]!=Dimensions[vec][[1]]||Tr[desc2,Times]!=Dimensions[vec][[2]],Print["ExchangeSystems: Wrong dimensions"]];list1={};list2={};For[i=0,i<=Dimensions[vec][[1]]-1,i++,list1=Insert[list1,1+FromDigits[Permute[IntDigs[i,desc1],newpos],MixedRadix[Permute[desc1,newpos]]],-1]];For[i=0,i<=Dimensions[vec][[2]]-1,i++,list2=Insert[list2,1+FromDigits[Permute[IntDigs[i,desc2],newpos],MixedRadix[Permute[desc2,newpos]]],-1]];Transpose[Permute[Transpose[Permute[vec,list1]],list2]]]

DirectSum[a_,b_]:=Module[{dim1,dim2,out,m},dim1=Dimensions[a][[1]];dim2=Dimensions[b][[1]];out=Join[a,ConstantArray[0,{dim1,dim2}],2];out=Join[out,Join[ConstantArray[0,{dim2,dim1}],b,2]];out]

DirectSum[list_]:=Module[{dim1,dim2,out,m,number,i},number=Dimensions[list][[1]];out=list[[1]];For[i=2,i<=number,i++,out=DirectSum[out,list[[i]]]];out]

QubitPartialTrace[M_,sys_]:=PT[M,Table[If[MemberQ[sys,i],0,1],{i,1,Log[2,Dimensions[M][[1]]]}],ConstantArray[2,Log[2,Dimensions[M][[1]]]]]

PartialTrace[mat_,dim1_,dim2_,s_]:=Module[{},If[(s==1||s==2)&&Dimensions[mat][[1]]==Dimensions[mat][[2]]==dim1*dim2,If[s==1,PT[mat,{0,1},{dim1,dim2}],PT[mat,{1,0},{dim1,dim2}]],Print["PartialTrace: input error"]]]

PT[mat_,keep_,desc_]:=Module[{dim,parts,mat2,i},dim=Tr[DeleteCases[keep*desc,0],Times];parts=Join[Flatten[Position[keep,0]],Flatten[Position[keep,1]]];mat2=Partition[ExchangeSystems[mat,Permute[Range[Length[desc]],parts],desc],{dim,dim}];Sum[mat2[[i,i]],{i,1,Tr[desc,Times]/dim}]]

BasisForm[vec_,desc_]:=Module[{i,dim,v=Flatten[vec]},dim=Tr[desc,Times];For[i=1,i<=dim,i++,If[v[[i]]!=0,Print[v[[i]],"|"<>StringDrop[StringDrop[ToString[IntDigs[i-1,desc]],1],-1]<>">"]]]]

BasisFormS[vec_,desc_]:=Module[{i,dim,string,v=Flatten[vec]},string="";dim=Tr[desc,Times];For[i=1,i<=dim,i++,If[v[[i]]>0,string=string<>"+"<>ToString[v[[i]],FormatType->StandardForm]<>"|"<>StringDrop[StringDrop[ToString[IntDigs[i-1,desc]],1],-1]<>"> ",If[v[[i]]<0,string=string<>ToString[v[[i]],FormatType->StandardForm]<>"|"<>StringDrop[StringDrop[ToString[IntDigs[i-1,desc]],1],-1]<>"> "]]];string]

Purify[rho_]:=Module[{dim,vals,vecs},dim=Dimensions[rho][[1]];{vals,vecs}=Eigensystem[rho];Sum[(vals[[i]])^(1/2)*Transpose[{vecs[[i]]}]\[CircleTimes]Transpose[{vecs[[i]]}],{i,1,dim}]]

SchmidtDecomposition[vec_,sys_]:=Module[{coeffs,vecsa,vecsb,u,w,v,vecmat,i},If[Tr[sys,Times]!=Dimensions[vec][[1]],Print["SchmidtDecomposition: Wrong dimensions"];Break[]];vecmat=Partition[Flatten[vec],sys[[2]]];{u,w,v}=SingularValueDecomposition[vecmat];coeffs={};vecsa={};vecsb={};For[i=1,i<=Min[sys],i++,coeffs=Insert[coeffs,w[[i,i]],-1];vecsa=Insert[vecsa,Transpose[{Transpose[u][[i]]}],-1];vecsb=Insert[vecsb,Transpose[{CT[v][[i]]}],-1]];{coeffs,vecsa,vecsb}]

DiagonalizingUnitary[M_]:=Module[{vals,vecs,i,U},{vals,vecs}=Eigensystem[M];U={};For[i=1,i<=Dimensions[M][[1]],i++,U=Insert[U,Normalize[vecs[[i]]],-1]];{Transpose[U],DiagonalMatrix[vals]}]      

OrderingF[vals_,prec_:Null]:=Module[{blocksizes,ord,ord2,v,i,start,out,block,blockv},v=vals;ord=Ordering[Map[Re,v]];blocksizes=Transpose[Tally[Map[Re,v[[ord]]],If[NumericQ[prec],Chop[N[#1-#2],prec]==0&,Chop[N[#1-#2]]==0&]]][[2]];(* note that numerical parts are taken before Tally *) start=1;out={};For[i=1,i<=Dimensions[blocksizes][[1]],i++,block=Take[ord,{start,start-1+blocksizes[[i]]}];blockv=Take[v[[ord]],{start,start-1+blocksizes[[i]]}];ord2=Ordering[Map[Im,blockv]];block=block[[ord2]];start = start + blocksizes[[i]];out=Join[out,block]];out]

EigenvaluesExact[m_,prec_:Null]:=Module[{vals,ord},vals=Eigenvalues[m];Off[Tally::smtst];ord=OrderingF[N[vals]];On[Tally::smtst];vals[[ord]]]

EigensystemExact[m_,prec_:Null]:=Module[{vals,vecs,i,blocksizes,ord,out,block,block2,startrow},{vals,vecs}=Eigensystem[m];Off[Tally::smtst];ord=OrderingF[N[vals]];vals=vals[[ord]];vecs=vecs[[ord]];blocksizes=Transpose[Tally[vals,If[NumericQ[prec], Chop[N[#1 - #2], prec] == 0 &,Chop[N[#1 - #2]] == 0 &]]][[2]];On[Tally::smtst];startrow = 1;out={};For[i=1,i<=Dimensions[blocksizes][[1]],i++,block=Take[vecs,{startrow,startrow-1+blocksizes[[i]]}];block2=Simplify[Orthogonalize[block]];block2=Simplify[Map[Normalize[#]&,block2]];startrow = startrow + blocksizes[[i]];out=Join[out,block2]];{vals,out}]

SimultaneouslyDiagonalize[A_,B_,precision_:Null]:=Module[{vals,vecs,unitary,blocksizes,startrow,out,i,Bblock,vals2,vec2,unitary2,A1,B1,mA,nA,mB,nB,valA,valB,vec,ord},{mA,nA}=Dimensions[A];{mB,nB}=Dimensions[B];If[A.B-B.A!=ConstantArray[0,{mA,nA}],A1=Chop[A];B1=Chop[B],A1=A;B1=B,A1=Simplify[A];B1=Simplify[B]];{vals,unitary}=EigensystemExact[A1,precision];blocksizes=Transpose[Tally[vals,If[NumericQ[precision],Chop[N[#1-#2],precision]==0&,Chop[N[#1-#2]]==0&]]][[2]];startrow=1;out={};For[i=1,i<=Dimensions[blocksizes][[1]],i++,Bblock=Chop[Take[unitary.B1.CT[unitary],{startrow,startrow-1+blocksizes[[i]]},{startrow,startrow-1+blocksizes[[i]]}]];startrow=startrow+blocksizes[[i]];If[Bblock==DiagonalMatrix[Diagonal[Bblock]],unitary2=IdentityMatrix[Dimensions[Bblock][[1]]],{vals2,unitary2}=EigensystemExact[Bblock,precision]];If[i==1,out=unitary2,out=DirectSum[out,unitary2]]];out=out.unitary;{valA,valB,vec}={Diagonal[out.A1.CT[out]],Diagonal[out.B1.CT[out]],out};(* The next If is not needed, but can be useful to flag problems *)If[Chop[N[CT[vec].DiagonalMatrix[valA].vec-A1]]!=ConstantArray[0,{mA,nA}]||Chop[N[CT[vec].DiagonalMatrix[valB].vec - B1]]!=ConstantArray[0, {mA, nA}],Print["Error in SimultaneouslyDiagonalize with inputs ",A,", ",B,", ",precision]];{valA,valB,vec}]

BlochSphere[rho_]:=Module[{x,y,z},x=Tr[rho.\[Sigma][1]];y=Tr[rho.\[Sigma][2]];z=Tr[rho.\[Sigma][3]];{{x},{y},{z}}]     

FromBlochSphere[{{x_},{y_},{z_}}]:=Simplify[(1/2)*(IdentityMatrix[2]+x*\[Sigma][1]+y*\[Sigma][2]+z*\[Sigma][3])]

AppendCols[m_]:=Module[{vec,mat,dim},(dim=Dimensions[m];mat=Transpose[m];mat=Join[mat,Table[ConstantArray[0,dim[[1]]],{dim[[1]]-dim[[2]]}]];vec=NullSpace[mat];mat[[-Length[vec];;-1]]=Conjugate[vec];mat=Map[Normalize,mat];Transpose[mat])]

FillZero[r_]:=Module[{dr,dc},{dr,dc}=Dimensions[r];If[dc-dr>0,Join[r,ConstantArray[0,{dc-dr,dc}]],r]]

MeasureBasis[rho_,basis_]:=Module[{out,i,j,dimrho,dimU,proj},dimrho=Dimensions[rho][[1]];dimU=Dimensions[basis][[1]];out=rho*0;For[i=1,i<=dimU,i++,proj=Transpose[{basis[[i]]}].Conjugate[{basis[[i]]}];out=out+(proj\[CircleTimes]IdentityMatrix[dimrho/dimU]).rho.(proj\[CircleTimes]IdentityMatrix[dimrho/dimU])];out]

MeasureBasis[rho_,basis_,sys_,desc_]:=Module[{out,i,j,dimrho,dimU,proj},dimrho=Dimensions[rho][[1]];If[Tr[desc,Times]!=dimrho,Print["MeasureBasis: Dimensions of desc should match those of matrix"]];dimU=Dimensions[basis][[1]];out=rho*0;For[i=1,i<=dimU,i++,proj=Transpose[{basis[[i]]}].Conjugate[{basis[[i]]}];out=out+Tensor[proj,IdentityMatrix[dimrho/dimU],sys,desc].rho.Tensor[proj,IdentityMatrix[dimrho/dimU],sys,desc]];out]

MeasurePOVM[rho_,POVM_,sys_,desc_]:=Module[{i,j,dimU,pos1},pos1=Position[sys,1];If[Dimensions[pos1][[1]]!=1,Print["MeasurePOVM: Can only measure on one subsystem"]];pos1=Tr[pos1];If[Tr[desc,Times]!=Dimensions[rho][[1]],Print["MeasurePOVM: Dimensions of desc should match those of matrix"]];dimU=Dimensions[POVM][[1]];Sum[DM[KetV[i-1,dimU]]\[CircleTimes]PT[Tensor[POVM[[i]],IdentityMatrix[Tr[desc,Times]/Dimensions[POVM[[1]]][[1]]],sys,desc].rho,sys-1,desc],{i,1,dimU}]]

POVMIsometry[rho_,POVM_,sys_,desc_]:=Module[{out,i,j,dimrho,dimU,proj,newdesc,pos1,kraus},dimrho=Dimensions[rho][[1]];pos1=Position[sys,1];If[Dimensions[pos1][[1]]!=1,Print["POVMIsometry: Can only measure on one subsystem"]];pos1=Tr[pos1];If[Tr[desc,Times]!=dimrho,Print["POVMIsometry: Dimensions of desc should match those of matrix"]];dimU=Dimensions[POVM][[1]];newdesc=ReplacePart[desc,pos1->dimU];kraus=Sum[Transpose[{UnitVector[dimU,i]}]\[CircleTimes]Transpose[{UnitVector[dimU,i]}]\[CircleTimes]Tensor[MatrixPower[POVM[[i]],1/2],IdentityMatrix[Tr[DeleteCases[desc*(sys-1),0],Times]],sys,desc],{i,1,dimU}];kraus.rho.CT[kraus]]

ChoiState[set_]:=Module[{ent},ent=DM[Sum[(KetV[i-1,Dimensions[set[[1]]][[2]]]\[CircleTimes]KetV[i-1,Dimensions[set[[1]]][[2]]]),{i,1,Dimensions[set[[1]]][[2]]}]];Sum[(IdentityMatrix[Dimensions[set[[1]]][[2]]]\[CircleTimes]set[[k]]).ent.(IdentityMatrix[Dimensions[set[[1]]][[2]]]\[CircleTimes]CT[set[[k]]]),{k,1,Dimensions[set][[1]]}]]

ChoiState[set1_,set2_]:=Module[{ent},If[Dimensions[set1][[1]]!=Dimensions[set2][[1]],Print["ChoiState Error: sets have different numbers of elements"]];ent=DM[Sum[(KetV[i-1,Dimensions[set1[[1]]][[2]]]\[CircleTimes]KetV[i-1,Dimensions[set1[[1]]][[2]]]),{i,1,Dimensions[set1[[1]]][[2]]}]];Sum[(IdentityMatrix[Dimensions[set1[[1]]][[2]]]\[CircleTimes]set1[[k]]).ent.(IdentityMatrix[Dimensions[set1[[1]]][[2]]]\[CircleTimes]CT[set2[[k]]]),{k,1,Dimensions[set1][[1]]}]]

ChoiChannel[state_,dA_,dB_]:=Module[{set1={},set2={},i,j,w1,w2,u,d,v},{u,d,v}=SingularValueDecomposition[state];For[i=1,i<=dA*dB,i++,If[Chop[d[[i,i]]]!=0,w1=(d[[i,i]])^(1/2)*Transpose[{Transpose[u][[i]]}];w2=(d[[i,i]])^(1/2)*Transpose[{Transpose[v][[i]]}];set1=Insert[set1,Transpose[Partition[Flatten[w1],dB]],-1];set2=Insert[set2,Transpose[Partition[Flatten[w2],dB]],-1]]];Chop[{set1,set2}]]

ChannelCompress[chan_]:=Module[{n,dA,dB},{n,dB,dA}=Dimensions[chan];ChoiChannel[ChoiState[chan], dA, dB][[1]]]

ExtremeChannelQ[list_]:=Module[{dim=Dimensions[list][[1]],i,j,newlist},newlist={};For[i=1,i<=dim,i++,For[j=1,j<=dim,j++,newlist=Insert[newlist,Flatten[CT[list[[i]]].list[[j]]],-1]]];If[MatrixRank[Chop[newlist]]==dim^2,True,False]]

PickRandomPsi[n_]:=Module[{psi,list1,i,R,phi},list1={};For[i=1,i<=n,i++,phi=Random[]*2*\[Pi];list1=Insert[list1,R[i,1]->Random[]*(Cos[phi]+I*Sin[phi]),-1]];psi=CreateMatrix[R,n,1]/.list1;psi/(Tr[Conjugate[Transpose[psi]].psi])^(1/2)]

RPickRandomPsi[n_]:=Module[{psi,list1,i,R},list1={};For[i=1,i<=n,i++,list1=Insert[list1,R[i,1]->Random[],-1]];psi=CreateMatrix[R,n,1]/.list1;psi/(Tr[Transpose[psi].psi])^(1/2)]

RPickRandomPsip[n_,rank_]:=Module[{psi,list1,i,R},list1={};For[i=1,i<=n,i++,If[i<=rank,list1=Insert[list1,R[i,1]->Random[],-1],list1=Insert[list1,R[i,1]->0,-1]]];psi=CreateMatrix[R,n,1]/.list1;psi/(Tr[CT[psi].psi])^(1/2)]

FPickRandomPsi[n_,prec_]:=Module[{psi,list1,i,R},list1={};For[i=1,i<=n,i++,list1=Insert[list1,R[i,1]->1+Random[Integer,prec-1]/(1+Random[Integer,prec-1]),-1]];psi=CreateMatrix[R,n,1]/.list1;psi/(Tr[Transpose[psi].psi])^(1/2)]

FPickRandomPsi[n_,prec_,rank_]:=Module[{psi,list1,i,R},list1={};For[i=1,i<=n,i++,If[i<=rank,list1=Insert[list1,R[i,1]->1+Random[Integer,prec-1]/(1+Random[Integer,prec-1]),-1],list1=Insert[list1,R[i,1]->0,-1]]];psi=CreateMatrix[R,n,1]/.list1;psi/(Tr[Transpose[psi].psi])^(1/2)]

PickRandomPsi2[n_]:=Module[{i,tot=1,out={},re,phi},For[i=1,i<=n-1,i++,re=RandomReal[tot];tot=tot-re;phi=Random[]*2*\[Pi];out=Insert[out,re^(1/2)*(Cos[phi]+I*Sin[phi]),-1]];phi=Random[]*2*\[Pi];out=Insert[out,tot^(1/2)*(Cos[phi]+I*Sin[phi]),-1];Transpose[{out}]]

RPickRandomPsi2[n_]:=Module[{i,tot=1,out={},re},For[i=1,i<=n-1,i++,re=RandomReal[tot];tot=tot-re;out=Insert[out,re^(1/2),-1]];out=Insert[out,tot^(1/2),-1];Transpose[{out}]]

FPickRandomPsi2[n_,prec_]:=Module[{i,tot=1,out={},re},For[i=1,i<=n-1,i++,re=1+Random[Integer,prec-1];re=tot*re/(re+Random[Integer,prec-1]);tot=tot-re;out=Insert[out,re^(1/2),-1]];out=Insert[out,tot^(1/2),-1];Transpose[{out}]]

PickRandomUnitary[n_]:=Module[{i,v,M},M={};For[i=1,i<=n,i++,v[i]=Flatten[PickRandomPsi[n]];M=Insert[M,v[i],-1]];Orthogonalize[M]]

RPickRandomUnitary[n_]:=Module[{i,v,M},M={};For[i=1,i<=n,i++,v[i]=Flatten[RPickRandomPsi[n]];M=Insert[M,v[i],-1]];Orthogonalize[M]]

FPickRandomUnitary[n_,prec_]:=Module[{i,v,M},M={};For[i=1,i<=n,i++,v[i]=Flatten[FPickRandomPsi[n,prec]];M=Insert[M,v[i],-1]];Orthogonalize[M]]

PickRandomIsometry[dim1_,dim2_]:=Module[{i,M1,M2},M1={};M2={};For[i=1,i<=dim1,i++,M1=Insert[M1,Flatten[PickRandomPsi[dim1]],-1];M2=Insert[M2,Flatten[PickRandomPsi[dim2]],-1]];M1=Orthogonalize[M1];M2=Orthogonalize[M2];Sum[CT[{M2[[i]]}].{M1[[i]]},{i,1,dim1}]]

RPickRandomIsometry[dim1_,dim2_]:=Module[{i,M1,M2},M1={};M2={};For[i=1,i<=dim1,i++,M1=Insert[M1,Flatten[RPickRandomPsi[dim1]],-1];M2=Insert[M2,Flatten[RPickRandomPsi[dim2]],-1]];M1=Orthogonalize[M1];M2=Orthogonalize[M2];Sum[CT[{M2[[i]]}].{M1[[i]]},{i,1,dim1}]]

FPickRandomIsometry[dim1_,dim2_,prec_]:=  Module[{i,M1,M2,ins1,ins2,recurse=0},M1={};M2={};For[i=1,i<=dim1,i++;recurse++,If[recurse==500,Print["FPickRandomIsometry:: failed to find an isometry after 500 tries, so output may not be an isometry, try increasing prec"];M1=Insert[M1,ins1,-1];M2=Insert[M2,ins2,-1];Break[]];ins1=Flatten[FPickRandomPsi[dim1,prec]];ins2=Flatten[FPickRandomPsi[dim2,prec]];If[MemberQ[M1,ins1]||MemberQ[M2,ins2],i--,M1=Insert[M1,ins1,-1];M2=Insert[M2,ins2,-1]]];M1=Orthogonalize[M1];M2=Orthogonalize[M2];Sum[CT[{M2[[i]]}].{M1[[i]]},{i,1,dim1}]]

PickRandomRho[n_]:=Module[{U},U=PickRandomUnitary[n];Chop[U.DiagonalMatrix[Flatten[(RPickRandomPsi[n])^2]].CT[U]]]

PickRandomRho[n_,rank_]:=Module[{U},U=PickRandomUnitary[n];Chop[U.DiagonalMatrix[Flatten[(RPickRandomPsip[n,rank])^2]].CT[U]]]   
          
RPickRandomRho[n_]:=Module[{U},U=RPickRandomUnitary[n];Chop[U.DiagonalMatrix[Flatten[(RPickRandomPsi[n])^2]].CT[U]]]          

RPickRandomRho[n_,rank_]:=Module[{U},U=RPickRandomUnitary[n];Chop[U.DiagonalMatrix[Flatten[(RPickRandomPsip[n,rank])^2]].CT[U]]]

FPickRandomRho[n_,prec_]:=Module[{U},U=FPickRandomUnitary[n,prec];Chop[U.DiagonalMatrix[Flatten[(FPickRandomPsi[n,prec])^2]].CT[U]]]                

FPickRandomRho[n_,prec_,rank_]:=Module[{U},U=FPickRandomUnitary[n,prec];Chop[U.DiagonalMatrix[Flatten[(FPickRandomPsi[n,prec,rank])^2]].CT[U]]]   

PickRandomPOVM[dim_,numels_]:=Module[{out,i,tr,cand,left},out={};tr=1;left=IdentityMatrix[dim];For[i=1,i<=numels-1,i++,cand=Random[]*tr*PickRandomRho[dim];While[Min[Chop[Eigenvalues[left-cand]]]<=0,cand=Random[]*tr*PickRandomRho[dim]];out=Insert[out,cand,-1];left=left-cand;tr=tr-Tr[cand]];Insert[out,left,-1]]
    
RPickRandomPOVM[dim_,numels_]:=Module[{out,i,tr,cand,left},out={};tr=1;left=IdentityMatrix[dim];For[i=1,i<=numels-1,i++,cand=Random[]*tr*RPickRandomRho[dim];While[Min[Chop[Eigenvalues[left-cand]]]<=0,cand=Random[]*tr*RPickRandomRho[dim]];out=Insert[out,cand,-1];left=left-cand;tr=tr-Tr[cand]];Insert[out,left,-1]]

PickRandomMeasurement[dim_,numels_]:=Module[{out,i,tr,cand,left},out={};tr=1;left=IdentityMatrix[dim];For[i=1,i<=numels-1,i++,cand=Random[]*tr*PickRandomUnitary[dim].PickRandomRho[dim];While[Min[Chop[Eigenvalues[left-CT[cand].cand]]]<=0,cand=Random[]*tr*PickRandomUnitary[dim].PickRandomRho[dim]];out=Insert[out,cand,-1];left=left-CT[cand].cand;tr=tr-Tr[CT[cand].cand]];Insert[out,PickRandomUnitary[dim].MatrixPower[left,1/2],-1]]

RPickRandomMeasurement[dim_,numels_]:=Module[{out,i,tr,cand,left},out={};tr=1;left=IdentityMatrix[dim];For[i=1,i<=numels-1,i++,cand=Random[]*tr*RPickRandomUnitary[dim].RPickRandomRho[dim];While[Min[Chop[Eigenvalues[left-CT[cand].cand]]]<=0,cand=Random[]*tr*RPickRandomUnitary[dim].RPickRandomRho[dim]];out=Insert[out,cand,-1];left=left-CT[cand].cand;tr=tr-Tr[CT[cand].cand]];Insert[out,RPickRandomUnitary[dim].MatrixPower[left,1/2],-1]]

PickRandomChannel[dimA_,dimB_,n_]:=Module[{iso=PickRandomIsometry[dimA,dimB*n],list1={},i},For[i=1,i<=n,i++,list1=Insert[list1,(IdentityMatrix[dimB]\[CircleTimes]BraV[i-1,n]).iso,-1]];list1]

RPickRandomChannel[dimA_,dimB_,n_]:=Module[{iso=RPickRandomIsometry[dimA,dimB*n],list1={},i},For[i=1,i<=n,i++,list1=Insert[list1,(IdentityMatrix[dimB]\[CircleTimes]BraV[i-1,n]).iso,-1]];list1]

FPickRandomChannel[dimA_,dimB_,n_,prec_]:=Module[{iso=FPickRandomIsometry[dimA,dimB*n,prec],list1={},i},For[i=1,i<=n,i++,list1=Insert[list1,(IdentityMatrix[dimB]\[CircleTimes]BraV[i-1,n]).iso,-1]];list1]

Dist[A_,B_]:=(1/2)*Tr[MatrixPower[CT[(A-B)].(A-B),1/2]]

QuickDist[A_,B_]:=(1/2)*Tr[Abs[Eigenvalues[A-B]]] 

Fidelity[A_,B_]:=Tr[MatrixPower[MatrixPower[A,1/2].B.MatrixPower[A,1/2],1/2]]

OptimumPOVM[A_,B_]:=Module[{dim,D1,U,F,Q,S,Q1,S1,i},dim=Extract[Dimensions[A],1];{U,D1}=DiagonalizingUnitary[A-B];F=AntiDiagonalMatrix[D1];Q1={};S1={};For[i=1,i<=dim,i++,If[Chop[Extract[Extract[F,i],1]]>0,Q1=Insert[Q1,1,-1];S1=Insert[S1,0,-1],S1=Insert[S1,1,-1];Q1=Insert[Q1,0,-1]]];Q=DiagonalMatrix[Q1];S=DiagonalMatrix[S1];Q=U.Q.CT[U];S=U.S.CT[U];{Q,S}]
     
OptimumPOVM[A_,B_,pA_]:=Module[{dim,D1,U,F,Q,S,Q1,S1,i},dim=Extract[Dimensions[A],1];{U,D1}=DiagonalizingUnitary[pA*A-(1-pA)*B];F=AntiDiagonalMatrix[D1];Q1={};S1={};For[i=1,i<=dim,i++,If[Chop[Extract[Extract[F,i],1]]>0,Q1=Insert[Q1,1,-1];S1=Insert[S1,0,-1],S1=Insert[S1,1,-1];Q1=Insert[Q1,0,-1]]];Q=DiagonalMatrix[Q1];S=DiagonalMatrix[S1];Q=U.Q.CT[U];S=U.S.CT[U];{Q,S}]

Matrixlog[A_]:=Module[{D1,U},{U,D1}=DiagonalizingUnitary[A];U.DiagonalMatrix[Log[Flatten[AntiDiagonalMatrix[D1]]]].Conjugate[Transpose[U]]]
	
Matrixlog[b_,A_]:=Module[{D1,U},{U,D1}=DiagonalizingUnitary[A];U.DiagonalMatrix[Log[b,Flatten[AntiDiagonalMatrix[D1]]]].Conjugate[Transpose[U]]]	

Matrixlog0[A_]:=Module[{D1,U,f},f[x_]:=If[x==0,0,Log[x]];{U,D1}=DiagonalizingUnitary[A];D1=Chop[D1];U.DiagonalMatrix[Map[f,Flatten[Diagonal[D1]]]].CT[U]]

Matrixlog0[b_,A_]:=Module[{D1,U,f},f[x_]:=If[x==0,0,Log[b,x]];{U,D1}=DiagonalizingUnitary[A];D1=Chop[D1];U.DiagonalMatrix[Map[f,Flatten[Diagonal[D1]]]].CT[U]]

Matrixxlogx[A_]:=Module[{D1,U,f},f[x_]:=If[x==0,0,x*Log[x]];{U,D1}=DiagonalizingUnitary[A];D1=Chop[D1];U.DiagonalMatrix[Map[f,Flatten[Diagonal[D1]]]].CT[U]]

Matrixxlogx[b_,A_]:=Module[{D1,U,f},f[x_]:=If[x==0,0,x*Log[b,x]];{U,D1}=DiagonalizingUnitary[A];D1=Chop[D1];U.DiagonalMatrix[Map[f,Flatten[Diagonal[D1]]]].CT[U]]

MatrixPower0[A_,a_]:=Module[{U,D1,f},f[x_]:=If[x==0,0,x^a];{U,D1}=DiagonalizingUnitary[A];D1=Chop[D1];U.DiagonalMatrix[Map[f,Flatten[Diagonal[D1]]]].CT[U]]

ShanEntropy[q_]:=Module[{i,p,ent},ent=0;If[!VectorQ[q],If[q<0||q>1,Print["ShanEntropy: Invalid input"]];ent=Limit[-(p*Log[p]+(1-p)*Log[1-p])/Log[2],p->q],For[i=1,i<=Dimensions[q][[1]],i++,If[q[[i]]<0||q[[i]]>1,Print["ShanEntropy: Invalid input"]];ent=ent-Limit[p*Log[p]/Log[2],p->q[[i]]]]];ent]

vNEntropy[rho_,keep_,desc_]:=Module[{i,keep2={},keep3={},rhoAB,rhoB},For[i=1,i<=Dimensions[desc][[1]],i++,If[keep[[i]]==2,keep2=Insert[keep2,1,-1];keep3=Insert[keep3,1,-1],keep2=Insert[keep2,keep[[i]],-1];keep3=Insert[keep3,0,-1]]];rhoAB=PT[rho,keep2,desc];rhoB=PT[rho,keep3,desc];ShanEntropy[Chop[Eigenvalues[rhoAB]]]-ShanEntropy[Chop[Eigenvalues[rhoB]]]]

vNInfo[rho_,keep_,desc_]:=Module[{i,keep1={},keep2={},keep3={}},For[i=1,i<=Dimensions[desc][[1]],i++,If[keep[[i]]==1,keep1=Insert[keep1,1,-1],keep1=Insert[keep1,0,-1]];If[keep[[i]]==2,keep2=Insert[keep2,1,-1],keep2=Insert[keep2,0,-1]];If[keep[[i]]==3,keep3=Insert[keep3,1,-1],keep3=Insert[keep3,0,-1]]];-vNEntropy[rho,keep1+keep2+keep3,desc]-vNEntropy[rho,keep3,desc]+vNEntropy[rho,keep1+keep3,desc]+vNEntropy[rho,keep2+keep3,desc]]

RelEnt[A_,B_]:=Module[{suppB},suppB=Support[B];If[Chop[suppB.A.CT[suppB]-A]==0*A,(Tr[Matrixxlogx[2,A]]-Tr[A.Matrixlog0[2,B]])/Tr[A],\[Infinity]]]

RelEnt[alpha_,A_,B_]:=Module[{supp,out},If[alpha==1,out=RelEnt[A,B]];If[alpha>1,supp=Support[B];If[Chop[supp.A.CT[supp]-A]==0*A,out=(1/(alpha-1))*Log[2,Tr[MatrixPower[A,alpha].MatrixPower0[B,1-alpha]]],out=\[Infinity]]];If[0<alpha<1,out=(1/(alpha-1))*Log[2,Tr[MatrixPower[A,alpha].MatrixPower[B,1-alpha]]]];If[alpha==0,out=-Log[2,Tr[Support[A].B]]];If[alpha<0,supp=Support[A];If[Chop[supp.B.CT[supp]-B]==0*B,out=(1/(alpha-1))*Log[2,Tr[MatrixPower0[A,alpha].MatrixPower[B,1-alpha]]],out=\[Infinity]]];out]

RenyiEnt[alpha_,rho_,keep_,desc_]:=Module[{i,keep2={},keep3={},rhoAB,rhoB,out},If[Dimensions[rho][[1]]!=Tr[desc,Times],Print["RenyiEnt: Wrong dimensions"]];If[Dimensions[desc][[1]]==1&&keep=={1},If[alpha==1,out=ShanEntropy[Chop[Eigenvalues[rho]]],out=1/(1-alpha)*Log[2,Tr[MatrixPower[rho,alpha]]]],For[i=1,i<=Dimensions[desc][[1]],i++,If[keep[[i]]==2,keep2=Insert[keep2,1,-1];keep3=Insert[keep3,1,-1],keep2=Insert[keep2,keep[[i]],-1];keep3=Insert[keep3,0,-1]]];rhoAB=PT[rho,keep2,desc];rhoB=PT[rho,keep3,desc];out=-RelEnt[alpha,rhoAB,Tensor[IdentityMatrix[Dimensions[rhoAB][[1]]/Dimensions[rhoB][[1]]],rhoB,DeleteCases[keep,0],Delete[desc,Position[keep,0]]]]];out]

RemoveIneqConstraints[M_,b_]:=Module[{outM,outb,i,j,obj,Mp,bp,ans},Off[LinearProgramming::lpsub];outM=M;outb=b;If[Dimensions[Dimensions[b]]=={2},For[i=Dimensions[M][[1]],i>=1,i--,If[b[[i]][[2]]==-1,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[-obj,Mp,bp];If[ans.obj>b[[i]][[1]]||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp,-1]];If[b[[i]][[2]]==1,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[obj,Mp,bp];If[ans.obj<b[[i]][[1]]||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp]]]];If[Dimensions[Dimensions[b]]=={1},For[i=Dimensions[M][[1]],i>=1,i--,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[obj,Mp,bp];If[ans.obj<b[[i]]||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp]]];On[LinearProgramming::lpsub];{outM,outb}]

NRemoveIneqConstraints[M_,b_]:=Module[{outM,outb,i,j,obj,Mp,bp,ans},Off[LinearProgramming::lpsub];outM=M;outb=b;If[Dimensions[Dimensions[b]]=={2},For[i=Dimensions[M][[1]],i>=1,i--,If[b[[i]][[2]]==-1,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[-obj*1.0,Mp,bp];If[Chop[ans.obj-b[[i]][[1]]]>0||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp,-1]];If[b[[i]][[2]]==1,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[obj*1.0,Mp,bp];If[Chop[ans.obj-b[[i]][[1]]]<0||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp]]]];If[Dimensions[Dimensions[b]]=={1},For[i=Dimensions[M][[1]],i>=1,i--,obj=outM[[i]];Mp=Drop[outM,{i}];bp=Drop[outb,{i}];ans=LinearProgramming[obj*1.0,Mp,bp];If[Chop[ans.obj-b[[i]]]<0||ans.obj===Indeterminate,j=1,outM=Mp;outb=bp]]];On[LinearProgramming::lpsub];{outM,outb}]

RemoveDuplicateConstraints[M_,b_]:=Module[{i,j,Mp,Mg,Ml,bp,gathered,max,maxpos,min,minpos,pos,drop,posg,posl},drop={};If[Dimensions[Dimensions[b]]=={2},Mp=Join[M,Transpose[{b[[All,2]]}],2];gathered=Gather[Mp];Mg=Join[M,Sign[Transpose[{b[[All,2]]}]+1/2],2];Ml=Join[M,Sign[Transpose[{b[[All,2]]}]-1/2],2];For[i=1,i<=Dimensions[gathered][[1]],i++,If[Dimensions[gathered[[i]]][[1]]>=2&&gathered[[i]][[1]][[-1]]==1,pos=Position[Mp,gathered[[i]][[1]]];posg=Position[Mg,gathered[[i]][[1]]];If[Dimensions[posg]!=Dimensions[pos],drop=Join[drop,pos],max=-\[Infinity];For[j=1,j<=Dimensions[pos][[1]],j++,If[b[[pos[[j]][[1]]]][[1]]>max,max=b[[pos[[j]][[1]]]][[1]];maxpos=pos[[j]][[1]]]];drop=Join[drop,DeleteCases[pos,{maxpos}]]]];If[Dimensions[gathered[[i]]][[1]]>=2&&gathered[[i]][[1]][[-1]]==-1,pos=Position[Mp,gathered[[i]][[1]]];posl=Position[Ml,gathered[[i]][[1]]];If[Dimensions[posl]!=Dimensions[pos],drop=Join[drop,pos],min=\[Infinity];For[j=1,j<=Dimensions[pos][[1]],j++,If[b[[pos[[j]][[1]]]][[1]]<min,min=b[[pos[[j]][[1]]]][[1]];minpos=pos[[j]][[1]]]];drop=Join[drop,DeleteCases[pos,{minpos}]]]];If[Dimensions[gathered[[i]]][[1]]>=2&&gathered[[i]][[1]][[-1]]==0,pos=Position[Mp,gathered[[i]][[1]]];drop=Join[drop,DeleteCases[pos,{i}]]]]];If[Dimensions[Dimensions[b]]=={1},gathered=Gather[M];For[i=1,i<=Dimensions[gathered][[1]],i++,If[Dimensions[gathered[[i]][[1]]][[1]]>=2,pos=Position[M,gathered[[i]][[1]]];max=-\[Infinity];For[j=1,j<=Dimensions[pos][[1]],j++,If[b[[pos[[j]][[1]]]]>max,max=b[[pos[[j]][[1]]]];maxpos=pos[[j]][[1]]]];drop=Join[drop,DeleteCases[pos,{maxpos}]]]]];{Delete[M,drop],Delete[b,drop]}]

Prep[M_,b_]:=Module[{out,bout,i,j},If[Dimensions[Dimensions[b]]!={2},Print["Prep: Error, expecting b to have the form {{b1,s1},{b2,s2},...}"]];out={};bout={};For[i=1,i<=Dimensions[M][[1]],i++,If[b[[i]][[2]]==0,out=Insert[out,M[[i]],-1];out=Insert[out,-M[[i]],-1];bout=Insert[bout,b[[i]][[1]],-1];bout=Insert[bout,-b[[i]][[1]],-1],If[b[[i]][[2]]==1,out=Insert[out,M[[i]],-1];bout=Insert[bout,b[[i]][[1]],-1],If[b[[i]][[2]]==-1,out=Insert[out,-M[[i]],-1];bout=Insert[bout,-b[[i]][[1]],-1]]]]];{out,bout}]

FourierMotzkin[M_, b_, elim_]:=Module[{out, bout, outp, i, j, min, mini, na, nb, el},If[Dimensions[Dimensions[b]] == {1}, out = M; bout = b,Print["FourierMotzkin: Error, expecting b to have the form {b1,b2,...}.  Use Prep[M,b] first if in the form {{b1,s1},{b2,s2},...}"]];el = Sort[elim];For[j = 1, j <= Dimensions[elim][[1]], j++, min = \[Infinity];mini = 1; outp = Join[out, IdentityMatrix[Dimensions[out][[2]]]];For[i = 1, i <= Dimensions[el][[1]], i++,na = Dimensions[PosInstances[outp, el[[i]]]][[1]];nb = Dimensions[NegInstances[outp, el[[i]]]][[1]];If[Dimensions[outp][[1]] - na - nb + na*nb < min,min = Dimensions[outp][[1]] - na - nb + na*nb; mini = i]];{out,bout} = Elim[out, bout, el[[mini]]]; {out, bout}=RemoveIneqConstraints[out, bout]; el = Drop[el, {mini}];el=el+PadRight[PadLeft[{}, mini - 1], Dimensions[el], -1]]; {out,bout}]

PosInstances[M_,i_]:=Module[{j,out={}},For[j=1,j<=Dimensions[M][[1]],j++,If[M[[j]][[i]]>0,out=Insert[out,j,-1]]];out]

NegInstances[M_,i_]:=Module[{j,out={}},For[j=1,j<=Dimensions[M][[1]],j++,If[M[[j]][[i]]<0,out=Insert[out,j,-1]]];out]

Elim[M_, b_, i_]:=Module[{out, bout, Mp, bp, j, k, insta, instb, c1, c2, d}, out = {};bout = {}; Mp = Insert[M, UnitVector[Dimensions[M][[2]], i], -1];bp = Insert[b, 0, -1];For[j = 1, j <= Dimensions[Mp][[1]], j++,If[Mp[[j]][[i]] == 0, out = Insert[out, Drop[Mp[[j]], {i}], -1];bout=Insert[bout,bp[[j]],-1]]];insta=PosInstances[Mp,i];instb=NegInstances[Mp,i];For[j=1,j<=Dimensions[insta][[1]],j++,For[k=1,k<=Dimensions[instb][[1]],k++,c1=Mp[[insta[[j]]]][[i]];c2=-Mp[[instb[[k]]]][[i]];d=LCM[c1,c2];out=Insert[out,Drop[d/c1*Mp[[insta[[j]]]]+d/c2*Mp[[instb[[k]]]], {i}], -1];bout=Insert[bout,d/c1*bp[[insta[[j]]]] + d/c2*bp[[instb[[k]]]], -1]]];{out,bout}]

IntDigs[num_,bases_]:=Module[{out,i,prod,prevprod,nm=num},out=Table[0,{Dimensions[bases][[1]]}];prod=1;For[i=1,i<=Dimensions[bases][[1]],i++,prevprod=prod;prod=bases[[-i]];out=ReplacePart[out,-i->Mod[nm,prod]];nm=(nm-Mod[nm,prod])/bases[[-i]]];out]

Progress[i_,imin_,imax_,num_]:=If[Floor[num*(i-imin)/(imax-imin)]-Floor[num*(i-imin-1)/(imax-imin)]==1,Print[100*Floor[num*(i-imin)/(imax-imin)]/num,"%"]]

ProgressTemporary[i_,imin_,imax_,num_]:=If[Floor[num*(i-imin)/(imax-imin)]-Floor[num*(i-imin-1)/(imax-imin)]==1,PrintTemporary[100*Floor[num*(i-imin)/(imax-imin)]/num,"%"]]

ThreadSolve[eqns_,soln_,var_]:=Module[{i,out={}},For[i=1,i<=Dimensions[eqns][[1]],i++,out=Insert[out,var/.Solve[eqns[[i]]==soln,var],-1]];Flatten[out,1]]

NThreadSolve[eqns_,soln_,var_]:=Module[{i,out={}},For[i=1,i<=Dimensions[eqns][[1]],i++,out=Insert[out,var/.NSolve[eqns[[i]]==soln,var],-1]];Flatten[out,1]]

\[Sigma][i_]:=If[i==0,{{1,0},{0,1}},If[i==1,{{0,1},{1,0}},If[i==2,{{0,-I},{I,0}},If[i==3,\[Sigma][3]={{1,0},{0,-1}}]]]];

Support[rho_]:=Module[{dim,i,U,Di,vals,a},dim=Dimensions[rho][[1]];{U,Di}=DiagonalizingUnitary[rho];vals=Flatten[AntiDiagonalMatrix[Chop[Di]]];vals=Limit[vals^a,a->0,Direction->-1];U.DiagonalMatrix[vals].CT[U]]

PosPart[mat_]:=Module[{i,vals,vecs,proj,warn=0},{vals,vecs}=Chop[Eigensystem[mat]];proj=0*IdentityMatrix[Dimensions[mat][[1]]];For[i=1,i<=Dimensions[vals][[1]],i++,If[vals[[i]]>0,proj=proj+Transpose[{Normalize[vecs[[i]]]}].Conjugate[{Normalize[vecs[[i]]]}],Null,If[warn==0,warn=1;Print["PosPart: unable to determine a sign"]]]];proj.mat.proj]

PosPart[mat_,assum_]:=Module[{i,vals,vecs,proj,warn=0},{vals,vecs}=Chop[Eigensystem[mat]];proj=0*IdentityMatrix[Dimensions[mat][[1]]];For[i=1,i<=Dimensions[vals][[1]],i++,If[Simplify[vals[[i]]>0,assum],proj=proj+Transpose[{Normalize[vecs[[i]]]}].Conjugate[{Normalize[vecs[[i]]]}],Null,If[warn==0,warn=1;Print["PosPart: unable to determine a sign"]]]];proj.mat.proj]

NegPart[mat_]:=Module[{i,vals,vecs,proj,warn=0},{vals,vecs}=Chop[Eigensystem[mat]];proj=0*IdentityMatrix[Dimensions[mat][[1]]];For[i=1,i<=Dimensions[vals][[1]],i++,If[vals[[i]]<0,proj=proj+Transpose[{Normalize[vecs[[i]]]}].Conjugate[{Normalize[vecs[[i]]]}],Null,If[warn==0,warn=1;Print["NegPart: unable to determine a sign"]]]];proj.mat.proj]

NegPart[mat_,assum_]:=Module[{i,vals,vecs,proj,warn=0},{vals,vecs}=Chop[Eigensystem[mat]];proj=0*IdentityMatrix[Dimensions[mat][[1]]];For[i=1,i<=Dimensions[vals][[1]],i++,If[Simplify[vals[[i]]<0,assum],proj=proj+Transpose[{Normalize[vecs[[i]]]}].Conjugate[{Normalize[vecs[[i]]]}],Null,If[warn==0,warn=1;Print["NegPart: unable to determine a sign"]]]];proj.mat.proj]

AntiDiagonalMatrix[M_]:=Module[{k,n,v},n=Extract[Dimensions[M],1];v={Extract[Extract[M,1],1]};For[k=2,k<=n,k++,v=Insert[v,Extract[Extract[M,k],k],k]];Transpose[{v}]]

CreateMatrix[W_,r_,c_]:=Module[{i,j,out,t},out={};For[i=1,i<=r,i++,For[j=1;t={},j<=c,j++,t=Insert[t,W[i,j],-1]];out=Insert[out,t,-1]];out]
     
CreateHermitianMatrix[W_,r_,c_]:=Module[{i,j,out,t},out={};For[i=1,i<=r,i++,For[j=1;t={},j<=c,j++,If[i<=j,t=Insert[t,W[i,j],-1],t=Insert[t,Conjugate[W[j,i]],-1]]];out=Insert[out,t,-1]];out]   

CreateHermitianMatrixR[W_,r_,c_]:=Module[{i,j,out,t,WR,WC},WR=ToExpression[ToString[W]<>"R"];WC=ToExpression[ToString[W]<>"C"];out={};For[i=1,i<=r,i++,For[j=1;t={},j<=c,j++,If[i==j,t=Insert[t,WR[i,j],-1],If[i<j,t=Insert[t,WR[j,i]+I*WC[j,i],-1],t=Insert[t,WR[i,j]-I*WC[i,j],-1]]]];out=Insert[out,t,-1]];out] 
 
CreateSymmetricMatrix[W_,r_,c_]:=Module[{i,j,out,t},out={};For[i=1,i<=r,i++,For[j=1;t={},j<=c,j++,If[i<=j,t=Insert[t,W[i,j],-1],t=Insert[t,W[j,i],-1]]];out=Insert[out,t,-1]];out]   
     
RemoveRepeatedInstances[list_]:=Module[{dim,i,j,k,list2,listg,RRI},RRI[listg_,k_]:=Module[{j1,list1,exit},list1=listg;exit=0;For[i=k,i>=1&&exit==0,i--,For[j1=i-1,j1>=1&&exit==0,j1--,If[Extract[listg,i]==Extract[listg,j1],exit=1;list1=Delete[list1,i]]]];list1];list2=list;dim=Extract[Dimensions[list2],1];list2=RRI[list,dim];While[i>=1,list2=RRI[list2,i]];list2]

SetMinus[s1_,s2_]:=Module[{out,i},out=s1;For[i=1,i<=Dimensions[s2][[1]],i++,out=DeleteCases[out,s2[[i]]]];out]

MajorizeQ[A_,B_]:=Module[{dim,i,out,Ap,Bp},dim=Dimensions[A][[1]];Ap=Sort[A,Greater];Bp=Sort[B,Greater];out=True;For[i=1,i<=dim,i++,If[Sum[Ap[[j]],{j,1,i}]<Sum[Bp[[j]],{j,1,i}],out=False;Break]];out]

ManualDeriv[fn_,vars_,point_,steps_]:=Module[{i,n=Dimensions[vars][[1]],out={},zero},zero=Table[0,n-1];For[i=1,i<=n,i++,out=Insert[out,ReleaseHold[((fn/.Thread[vars->point+Insert[zero,steps[[i]],i]])-(fn/.Thread[vars->point]))/steps[[i]]],-1]];out]

QRDecomp[M_]:=Module[{q,r,dimq,dimm=Dimensions[M]},If[dimm[[1]]!=dimm[[2]],Print["QRDecomp warning: input not a square matrix, so output may not be of expected form."]];{q,r}=QRDecomposition[M];dimq=Dimensions[q];If[dimq[[1]]!=dimq[[2]],q=CT[AppendCols[CT[q]]];r=FillZero[r]];{q,r}]

RQDecomp[M_]:=Module[{Q,R},{Q,R}=QRDecomp[CT[Reverse[M]]];{M.CT[Reverse[Q]],CT[Reverse[Q]]}]

LQDecomp[M_]:=Module[{Q,R},{Q,R}=QRDecomp[CT[M]];{CT[R],CT[Q]}]

QLDecomp[M_]:=Module[{Q,R},{R,Q}=RQDecomp[CT[M]];{CT[Q],CT[R]}]

SignZeroPos[x_]:=Module[{out={},i},For[i=1,i<=Dimensions[x][[1]],i++,out=Insert[out,If[Re[x[[i]]]<0,-1,1],-1]];out]

SignZeroNeg[x_]:=Module[{out={},i},For[i=1,i<=Dimensions[x][[1]],i++,out=Insert[out,If[Re[x[[i]]]>0,1,-1],-1]];out]

QRDecompPos[m_]:=Module[{q,r,s},{q,r}=QRDecomp[m];s=DiagonalMatrix[SignZeroPos[Diagonal[r]]];{s.q,s.r}]

QRDecompNeg[m_]:=Module[{q,r,s},{q,r}=QRDecomp[m];s=DiagonalMatrix[-SignZeroNeg[Diagonal[r]]];{s.q,s.r}]

RQDecompPos[m_]:=Module[{q,r,s},{r,q}=RQDecomp[m];s=DiagonalMatrix[SignZeroPos[Diagonal[r]]];{r.s,q.s}]

RQDecompNeg[m_]:=Module[{q,r,s},{r,q}=RQDecomp[m];s=DiagonalMatrix[-SignZeroNeg[Diagonal[r]]];{r.s,q.s}]

LQDecompPos[m_]:=Module[{q,l,s},{l,q}=LQDecomp[m];s=DiagonalMatrix[SignZeroPos[Diagonal[l]]];{l.s,q.s}]

LQDecompNeg[m_]:=Module[{q,l,s},{l,q}=LQDecomp[m];s=DiagonalMatrix[-SignZeroNeg[Diagonal[l]]];{l.s,q.s}]

QLDecompPos[m_]:=Module[{q,l,s},{q,l}=QLDecomp[m];s=DiagonalMatrix[SignZeroPos[Diagonal[l]]];{s.q,s.l}]

QLDecompNeg[m_]:=Module[{q,l,s},{q,l}=QLDecomp[m];s=DiagonalMatrix[-SignZeroNeg[Diagonal[l]]];{s.q,s.l}]



End[ ]

EndPackage[ ]
