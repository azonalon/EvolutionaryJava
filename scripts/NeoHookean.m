(* ::Package:: *)

(* ::Input::Initialization:: *)
NeoHookeanPotential[F_,\[Lambda]_,\[Mu]_]:=(
\[Mu] /2(Tr[F\[Transpose].F]-2)-\[Mu] Log[Det[F]]+\[Lambda]/2 Log[Det[F]]^2
)
NeoHookeanPotential[\[Sigma]x_,\[Sigma]y_,\[Lambda]_,\[Mu]_]:=1/2 (\[Mu] (-2+\[Sigma]x^2+\[Sigma]y^2)-2 \[Mu] Log[\[Sigma]x \[Sigma]y]+\[Lambda] Log[\[Sigma]x \[Sigma]y]^2)
NeoHookeanStress[F_,\[Lambda]_,\[Mu]_]:=(
\[Mu](F-Inverse[F]\[Transpose])+\[Lambda] Log[Det[F]] Inverse[F]\[Transpose]
)
NeoHookeanStressDifferential[F_,dF_,\[Lambda]_,\[Mu]_]:=(
\[Mu] dF+(\[Mu]-\[Lambda] Log[Det[F]])Inverse[F]\[Transpose].dF\[Transpose].Inverse[F]\[Transpose]+\[Lambda] Tr[Inverse[F]\[Transpose].dF\[Transpose]] Inverse[F]\[Transpose]
)


(* ::Input:: *)
(**)
(*vA={{0,0},{0,1},{1,0}};*)
(*vInit = vA+0.7{{0,1},{0,0},{0,0}};*)
(*dv = {{-.05,.1},{-.2,.3},{.4,.5}};*)
(*t1=10;*)
(*fExt[t_]={{1,1},{-1,0},{0,-1}}HeavisideTheta[3-t];*)
(*vertices[t_]=ElasticTriangleSimulation[*)
(*vA, fExt,*)
(*NeoHookeanStress[#,1,1]&,*)
(*NeoHookeanStressDifferential[#1,#2,1,1]&,t1,*)
(*.1];*)
(*ParametricPlot[vertices[t],{t,0,t1}]*)
(*Animate[*)
(*Graphics[Triangle[vertices[t]],Axes->True,PlotRange->{{-2,2},{-2,2}}]*)
(*,{t,0,t1}]*)
(**)
(**)
(**)
(**)


(* ::Input:: *)
(*F[vA, vertices[0]]*)
(*Plot[Det[F[vA, vertices[t]]],{t,0,t1},PlotRange->Full]*)


(* ::Input:: *)
(*FSVD = ({*)
(* {Cos[\[Phi]], -Sin[\[Phi]]},*)
(* {Sin[\[Phi]], Cos[\[Phi]]}*)
(*}).({*)
(* {sx, 0},*)
(* {0, sy}*)
(*}).({*)
(* {Cos[\[Theta]], -Sin[\[Theta]]},*)
(* {Sin[\[Theta]], Cos[\[Theta]]}*)
(*});*)
(*\[Psi]v[sx_,sy_]=FullSimplify[NeoHookeanPotential[FSVD,\[Lambda],\[Mu]]];*)
(*vf=-D[\[Psi]v[sx,sy],{{sx,sy}}]/.{\[Mu]->1,\[Lambda]->1};*)
(*vp=VectorPlot[Re[vf],{sx,-2,2},{sy,-2,2},VectorScale->{Automatic,Automatic,None},VectorStyle->Tiny,*)
(*VectorPoints->Fine];*)
(*grad=D[\[Psi]v[sx,sy],{{sx, sy}}];*)
(*Hess=FullSimplify[D[grad,{{sx,sy}}]];*)
(*temp=FullSimplify[Eigensystem[Hess]];*)
(*Simplify[temp[[1,1]]-temp[[1,2]]];*)
(*eigenv=temp[[2,2]];*)
(*cont = Re[(grad.eigenv)]/.{\[Mu]->1,\[Lambda]->1};*)
(*Show[{vp,ContourPlot[cont==0,{sx,-2,2},{sy,-2,2}]}]*)
(**)
