BeginPackage["SpecialRelativity`"]

mkFourVector::usage = "mkFourVector[t,x,y,z] - create a four-vector"
FourVectorQ::usage = "check if argument is a four-vector"
LightlikeQ::usage = "check if argument is a lightlike four-vector"
TimelikeQ::usage = "check if argument is a timelike four-vector"
SpacelikeQ::usage = "check if argument is a spacelike four-vector"

Dtau::usage = "Dtau[wordline][lambda]"
ProperTime::usage = "ProperTime[wordline][lambda]"
FourVelocity::usage = "FourVelocity[wordline][lambda]"
ProperVelocity::usage = "ProperVelocity[wordline][lambda]"
FourAcceleration::usage = "FourAcceleration[wordline][lambda]"
ProperAcceleration::usage = "ProperAcceleration[wordline][lambda]"

FramePosition::usage = "FramePosition[fourVector]"
FrameTime::usage = "FrameTime[fourVector]"
FrameVelocity::usage = "LocalVelocity[wordline][lambda]"
FrameAcceleration::usage = "LocalAcceleration[wordline][lambda]"


FrameTimeParametrization::usage = "FrameTimeParametrization[wordline][lambda]"
ProperTimeParametrization::usage = "ProperTimeParametrization[wordline][lambda]"

ProperTimeToFrameTime::usage = "ProperTimeToFrameTime[wordline][tau]"
FrameTimeToProperTime::usage = "FrameTimeToProperTime[wordline][t]"

(* Conversion Between SI and Natural Units *)
ToNaturalUnits::usage = "ToNaturalUnits[quantity]"
FromNaturalUnits::usage = "FromNaturalUnits[quantity, siUnits]"

LorenzBoost::usage = "LorenzBoost[{x_, y_, z_}, v_] - lorenz boost matrix in a given direction"

Begin["Private`"]

mkFourVector[{t_,x_,y_,z_}] := FourVector[t,x,y,z]
mkFourVector[t_,x_,y_,z_] := FourVector[t,x,y,z]
mkFourVector[t_,x_,y_] := FourVector[t,x,y,0]
mkFourVector[t_,x_] := FourVector[t,x,0,0]
mkFourVector[t_] := FourVector[t,0,0,0]

FourVectorQ[FourVector[__]] := True
FourVectorQ[x_] := False

MakeBoxes[FourVector[t_,x_,y_,z_],form_]:= Module[{
rowBox=RowBox[{"FourVector","[",
	MakeBoxes[t,form], ",", 
	MakeBoxes[x,form], ",", 
	MakeBoxes[y,form], ",",
	MakeBoxes[z,form],
  "]"}],
  expr=FourVector[t,x,y,z]},
  InterpretationBox[#1,#2]&@@{rowBox,expr}
]

(* Dot Products *)

FourVector /: Dot[FourVector[t1_,x1_,y1_,z1_],FourVector[t2_,x2_,y2_,z2_]] := t1 t2 - (x1 x2 + y1 y2 + z1 z2)
FourVector /: Dot[
  {{a11_,a12_,a13_,a14_},{a21_,a22_,a23_,a24_},{a31_,a32_,a33_,a34_},{a41_,a42_,a43_,a44_}},
  FourVector[t_,x_,y_,z_]
  ] := mkFourVector[
    {{a11,a12,a13,a14},{a21,a22,a23,a24},{a31,a32,a33,a34},{a41,a42,a43,a44}} . {t,x,y,z}
   ]

FourVector /: Norm[FourVector[v__]] := Sqrt[FourVector[v].FourVector[v]]

(* Algebra *)

FourVector /: Divide[FourVector[t_,x_,y_,z_],s_] := FourVector[t/s,x/s,y/s,z/s]
FourVector /: Times[FourVector[t_,x_,y_,z_],s_] := FourVector[t*s,x*s,y*s,z*s]
FourVector /: Plus[FourVector[t1_,x1_,y1_,z1_],FourVector[t2_,x2_,y2_,z2_]] := FourVector[t1+t2,x1+x2,y1+y2,z1+z2]


LightlikeQ[v_FourVector] := v.v == 0
TimelikeQ[v_FourVector] := v.v > 0
SpacelikeQ[v_FourVector] := v.v < 0

(* Working with word lines *)

FourVector /: Dt[FourVector[t_,x_,y_,z_]] := FourVector[Dt[t],Dt[x],Dt[y],Dt[z]]
FourVector /: Dt[FourVector[t_,x_,y_,z_], vars_] := FourVector[Dt[t, vars],Dt[x, vars],Dt[y, vars],Dt[z, vars]]

(* Proper Time and Time Transformations *)

Dtau[wordLine_] := Function[lambda,
Simplify[Norm[Dt[wordLine[lambda]]], Assumptions->{Dt[lambda]>0}]]

ProperTime[wordLine_] := Integrate[
Dtau[wordLine][lambda]/Dt[lambda],{lambda,0,#}, Assumptions->{Dt[lambda]>0}] &

ProperTimeToFrameTime[wordLine_]:=Function[tau, Module[
{lambda},
(FrameTime[wordLine[lambda]] /.#)& /@Solve[tau ==ProperTime[wordLine][lambda],lambda,Reals]]]

FrameTimeToProperTime[wordLine_]:=Function[t, Module[
{lambda},
(ProperTime[wordLine][lambda] /.#)& /@Solve[t == FrameTime[wordLine[lambda]],lambda,Reals]]]

ProperTimeDerivative[expr_, wordLine_]:= Function[lambda,
Simplify[Dt[expr[lambda]]/Dtau[wordLine][lambda], Assumptions->Dt[lambda]>0]]

(* Kinematics *)

FourVelocity[wordLine_] := ProperTimeDerivative[wordLine,wordLine]
ProperVelocity[wordLine_] := Norm[FourVelocity[wordLine][#]]&

FourAcceleration[wordLine_] := ProperTimeDerivative[FourVelocity[wordLine],wordLine]
ProperAcceleration[wordLine_]:=Norm[FourAcceleration[wordLine][#]]&

(* Parametrizations *)

FrameTimeParametrization[wordLine_] := Function[t, Module[{lambda},
(wordLine[lambda] /. #&) /@ Solve[t == FrameTime[wordLine[lambda]], lambda, Reals][[1]]]]

ProperTimeParametrization[wordLine_] := Function[t, Module[{lambda},
(wordLine[lambda] /. #&) /@ Solve[t == ProperTime[wordLine][lambda], lambda, Reals][[1]]]]

(* Local Frame Computations *)

FramePosition[FourVector[t_, x_, y_, z_]] := {x, y, z}
FrameTime[FourVector[t_, x_, y_, z_]] := t

FrameVelocity[wordLine_] := Function[lambda,
Module[
	{d=Dt[wordLine[lambda], lambda]}, 
	FramePosition[d/d[[1]]]
]]

FrameAcceleration[wordLine_] := Function[lambda,
Module[
	{
		d=Dt[wordLine[lambda], lambda],
		dv=Dt[FrameVelocity[wordLine][lambda], lambda]
	}, 
	dv/d[[1]]
]]


(* Unit Conversion *)

TimeUnitPower[dims_] := Module[
  {timeUnits = Cases[dims, {"TimeUnit", t_} -> t]},
  If[Length[timeUnits] == 0, 0,timeUnits[[1]] ]
]

ToNaturalUnits[q_] := Module[
  {si = UnitConvert[q], dims},
  dims = UnitDimensions[si];
  si*Quantity[3*^8, "Meters/Seconds"]^TimeUnitPower[dims]
]

FromNaturalUnits[q_, siUnits_] := Module[
  {dims},
  dims = UnitDimensions[Quantity[1, siUnits]];
  q/Quantity[3*^8, "Meters/Seconds"]^TimeUnitPower[dims]
]

(* Lorenz Boost *)

LorenzBoost[{x_, y_, z_}, v_] := Module[
  { 
    gamma = 1/Sqrt[1 - v^2], 
    nx, ny, nz
  },
  { nx, ny, nz } = Normalize[{x, y, z}];
  {
   {       gamma,          -gamma  v nx,         - gamma  v ny,          -gamma v nz  },
   { -gamma v nx, 1 + (gamma - 1)  nx^2,     (gamma - 1) nx ny,     (gamma - 1) nx nz },
   { -gamma v ny,     (gamma - 1) ny nx, 1 + (gamma - 1)  ny^2,     (gamma - 1) ny nz },
   { -gamma v nz,     (gamma - 1) nz nx,     (gamma - 1) nz ny, 1 + (gamma - 1) nz^2  }
  }
]

End[]

EndPackage[]
  