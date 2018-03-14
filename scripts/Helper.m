(* ::Package:: *)

(* ::Input::Initialization:: *)
CoefficientValueTuples[x_,coeffs_]:=Module[{l,n},
l=CoefficientList[x,coeffs]//Simplify;
n=Length[coeffs];
Cases[Flatten[MapIndexed[
{Product[coeffs[[i]]^(#2[[i]]-1),{i,1,n}],#1}&,l
,{n}],n-1],Except[{_,0}]]
]
ApproximateHessian[f_,x_,y_,h_]:=1/h^2 {{f[x+2h ,y]-2f[x+h ,y]+f[x,y],f[x+h ,y+h ]-f[x+h ,y]-f[x,y+h ]+f[x,y]},
{f[x+h ,y+h ]-f[x+h ,y]-f[x,y+h ]+f[x,y],f[x,y+2h ]-2f[x ,y+h]+f[x,y]}
};
Atan2[x_,y_]:=Arg[I x+y];
Atan2[y_,x_]:=2 ArcTan[y/(Sqrt[x^2+y^2]+x)];
svdsubst={f^2+g^2->R^2,e^2+h^2->Q^2};
svdbacksubst={R-> Sqrt[f^2+g^2],Q-> Sqrt[e^2+h^2]};
hypot[x_,y_]:=Sqrt[x^2+y^2];

SVD2D[m_]:=Module[{a1,a2,sx,sy,\[Theta],\[Phi],e,f,g,h,Q,R,u1,u2,v1,v2},
e=1/2 (m[[1,1]]+m[[2,2]]);
f=1/2 (m[[1,1]]-m[[2,2]]);
g=1/2 (m[[2,1]]+m[[1,2]]);
h=1/2 (m[[2,1]]-m[[1,2]]);
Q=hypot[e,h];R=hypot[f,g];
sx=Q+R;sy=Q-R;
a1=Atan2[g,f];a2=Atan2[h,e];
(*a1=Arg[I g+ f];a2=Arg[I h+e];*)
\[Theta]=(a2-a1)/2;\[Phi]=(a2+a1)/2;
{Cos[\[Phi]],Sin[\[Phi]],sx,sy,Cos[\[Theta]],Sin[\[Theta]]}
]



(* ::Input:: *)
(*CoefficientValueTuples[x+3y+z+4+y z^7 8 +a y^2 z+u^7,{x,y,z}]*)
(*D[Sin[x y],{{x,y},2}]/.{x->0,y->0}*)
(*ApproximateHessian[Sin[#1 #2]&,0,0,0.00001]*)


(* ::Input:: *)
(**)
(*Block[{m},*)
(*m=RandomReal[{-1,1},{2,2}];*)
(*(*m={{1,0},{0,1}};*)*)
(*Print@SVD2D[m];*)
(*SingularValueDecomposition[m]*)
(*]*)
