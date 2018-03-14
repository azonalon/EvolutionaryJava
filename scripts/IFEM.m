(* ::Package:: *)

(* ::Input::Initialization:: *)
$Assumptions=_\[Element]Reals;
\[Lambda]L[k_,\[Nu]_]=(k \[Nu])/((1+\[Nu])(1-2\[Nu]));
\[Mu]L[k_,\[Nu]_]=k /(2(1+\[Nu]));;
LameCoefficients[\[Kappa]_,\[Nu]_]={\[Lambda][k,\[Nu]], \[Mu][k,\[Nu]]};
Id= IdentityMatrix[2];

Ds[tr_]:=Transpose[tr[[1;;2]]-{tr[[3]],tr[[3]]}];
Bm[tr_]:=Inverse[Ds[tr]]

ElasticForce[P_,vA_]:=Module[{H},
H=-Abs[Det[Ds[vA]]]/2P.(Inverse[Ds[vA]])\[Transpose];
{H\[Transpose][[1]],H\[Transpose][[2]],-H\[Transpose][[1]]-H\[Transpose][[2]]}
]

F[vA_,vB_]:=Ds[vB].Bm[vA];
ScanLine[\[Alpha]_,v_,vA_,k_,\[Nu]_]:=Module[{P,f},
f=ElasticForce[VenantKirchhoffStress[F[vA,v],k,\[Nu]],vA];
P=VenantKirchhoffStress[F[vA,v+\[Alpha] f],k,\[Nu]];
{-Flatten[ElasticForce[P,vA]].Flatten[f],
VenantKirchhoffPotential[F[vA,v+\[Alpha] f],k,\[Nu]],f}
]
W[x_]:=Area[Simplex[x]];

StressHessian[potential_,XX_,xx_,\[Mu]_,\[Lambda]_]:=Module[{xT,XT},
xT=Table[x[2i+j],{i,0,2},{j,1,2}];
XT=Table[X[2i+j],{i,0,2},{j,1,2}];
-W[XT]Table[D[potential[F[XT,xT],\[Mu],\[Lambda]],x[i],x[j]],{i,1,6},{j,1,6}]/.{x[i_]-> Flatten[xx][[i]],X[i_]->  Flatten[XX][[i]]}
];
StressGradient[potential_,XX_,xx_,\[Mu]_,\[Lambda]_]:=Module[{xT,XT},
xT=Table[x[2i+j],{i,0,2},{j,1,2}];
XT=Table[X[2i+j],{i,0,2},{j,1,2}];
-W[XT]Table[D[potential[F[XT,xT],\[Mu],\[Lambda]],x[i]],{i,1,6}]/.{x[i_]-> Flatten[xx][[i]],X[i_]->  Flatten[XX][[i]]}
];
MatrixComponents[A_]:=Table[UnitVector[6,i].A[UnitVector[6,j]],{i,1,6},{j,1,6}];
ElasticModelTest[vA_,vInit_,dv_,\[Psi]_,f_,df_,\[Mu]_,\[Lambda]_]:=Module[{m1,m2,g1,g2},
m1=StressHessian[\[Psi],vA,vInit,\[Mu],\[Lambda]];
m2=MatrixComponents[
Flatten[ElasticForce[
df[F[vA,vInit],F[vA,Partition[#,2]],\[Mu],\[Lambda]],vA
]]&
];
Print[AllTrue[Chop[m1-m2],#==0&,2]];
g1=StressGradient[\[Psi],vA,vInit,\[Mu],\[Lambda]];
g2=Flatten[ElasticForce[f[F[vA,vInit],\[Mu],\[Lambda]],vA]];
Print[AllTrue[Chop[g1-g2],#==0&]];
{W[vA]\[Psi][F[vA,vInit],\[Mu],\[Lambda]],Partition[g2,1],Partition[m1.Flatten[dv],1]}
]
ElasticTriangleSimulation[vertices_, externalForce_,P_,dP_,t1_,kDamp_:0.3]:=Module[{state, rhs, lhs},
state[t_]=Table[i[j][t],{j,1,3},{i,{xc,yc}}];
rhs=ElasticForce[P[F[vertices,state[t]]],vertices]+
     kDamp  ElasticForce[dP[F[vertices,state[t]],F[vertices,state'[t]]],vertices]+
externalForce[t];
lhs=state''[t];
state[t]/.NDSolve[{lhs==rhs,
state[0]==vertices,state'[0]==0vertices},state[t],{t,0,t1}][[1]]
]



(* ::Input:: *)
(*Ds3D=({*)
(* {x1-x4, x2-x4, x3-x4},*)
(* {y1-y4, y2-y4, y3-y4},*)
(* {z1-z4, z2-z4, z3-z4}*)
(*})/.{z4->-z1-z2-z3,y4->-y1-y2-y3,x4->-x1-x2-x3};*)
(*Ds2D=({*)
(* {x1-x3, x2-x3},*)
(* {y1-y3, y2-y3}*)
(*});*)
(*a=Sqrt[(x3-x1)^2+(y3-y1)^2];*)
(*b=Sqrt[(x3-x2)^2+(y3-y2)^2];*)
(*c=Sqrt[(x1-x2)^2+(y1-y2)^2];*)
(*s=(a+b+c)/2;*)
(*Vh=Sqrt[s(s-a)(s-b)(s-c)];*)
(*FullSimplify[{Vh,Det[Ds2D]},Reals]*)



