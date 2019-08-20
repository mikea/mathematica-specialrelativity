(* ::Package:: *)

BeginPackage["SpecialRelativity`"]

mkFourVector::usage = "mkFourVector[t,x,y,z] - create a four-vector"
FourVectorQ::usage = "check if argument is a four-vector"
LightlikeQ::usage = "check if argument is a lightlike four-vector"
TimelikeQ::usage = "check if argument is a timelike four-vector"
SpacelikeQ::usage = "check if argument is a spacelike four-vector"

Dtau::usage = "Dtau[worldLine][lambda]"
ProperTimeRate::usage = "ProperTimeRate[worldLine][lambda]"
ProperTime::usage = "ProperTime[worldLine][lambda]"

FourVelocity::usage = "FourVelocity[worldLine][lambda]"
ProperVelocity::usage = "ProperVelocity[worldLine][lambda]"
FourAcceleration::usage = "FourAcceleration[worldLine][lambda]"
ProperAcceleration::usage = "ProperAcceleration[worldLine][lambda]"

FramePosition::usage = "FramePosition[fourVector]"
FrameTime::usage = "FrameTime[fourVector]"
FrameVelocity::usage = "LocalVelocity[worldLine][lambda]"
FrameAcceleration::usage = "LocalAcceleration[worldLine][lambda]"


FrameTimeParametrization::usage = "FrameTimeParametrization[worldLine][lambda]"
ProperTimeParametrization::usage = "ProperTimeParametrization[worldLine][lambda]"

ProperTimeToFrameTime::usage = "ProperTimeToFrameTime[worldLine][tau]"
FrameTimeToProperTime::usage = "FrameTimeToProperTime[worldLine][t]"

(* Conversion Between SI and Natural Units *)
ToNaturalUnits::usage = "ToNaturalUnits[quantity]"
FromNaturalUnits::usage = "FromNaturalUnits[quantity, siUnits]"

LorenzBoost::usage = "LorenzBoost[{vx_, vy_, vz_}] - lorenz boost matrix in a given direction"
MCRFBoost::usage = "MCRFBoost[worldLine][tau] - lorenz boost of MCRF"

Begin["Private`"]

(* TODO: vector element access *)

(* Four Vectors *)

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

Dtau[worldLine_] := Function[lambda,
  Simplify[Norm[Dt[worldLine[lambda]]], Assumptions->{Dt[lambda]>0}]]
ProperTimeRate[worldLine_] := Dtau[worldLine]

ProperTime[worldLine_] := Integrate[
  Dtau[worldLine][lambda]/Dt[lambda],{lambda,0,#}, Assumptions->{Dt[lambda]>0}] &

ProperTimeToFrameTime[worldLine_]:=Function[tau, Module[
  {lambda},
  (FrameTime[worldLine[lambda]] /.#)& /@Solve[tau ==ProperTime[worldLine][lambda],lambda,Reals]]]

FrameTimeToProperTime[worldLine_]:=Function[t, Module[
  {lambda},
  (ProperTime[worldLine][lambda] /.#)& /@Solve[t == FrameTime[worldLine[lambda]],lambda,Reals]]]

ProperTimeDerivative[expr_, worldLine_]:= Function[lambda,
  Simplify[Dt[expr[lambda]]/Dtau[worldLine][lambda], Assumptions->Dt[lambda]>0]]

(* Kinematics *)

FourVelocity[worldLine_] := ProperTimeDerivative[worldLine,worldLine]
ProperVelocity[worldLine_] := Norm[FourVelocity[worldLine][#]]&

FourAcceleration[worldLine_] := ProperTimeDerivative[FourVelocity[worldLine],worldLine]
ProperAcceleration[worldLine_]:=Norm[FourAcceleration[worldLine][#]]&

(* Parametrizations *)

FrameTimeParametrization[worldLine_] := Function[t, Module[{lambda},
  (worldLine[lambda] /. #&) /@ Solve[t == FrameTime[worldLine[lambda]], lambda, Reals][[1]]]]

ProperTimeParametrization[worldLine_] := Function[t, Module[{lambda},
  (worldLine[lambda] /. #&) /@ Solve[t == ProperTime[worldLine][lambda], lambda, Reals][[1]]]]

(* Local Frame Computations *)

FramePosition[worldLine_] := Function[lambda, FramePosition[worldLine[lambda]]]
FramePosition[FourVector[t_, x_, y_, z_]] := {x, y, z}

FrameTime[worldLine_] := Function[lambda, FrameTime[worldLine[lambda]]]
FrameTime[FourVector[t_, x_, y_, z_]] := t

FrameVelocity[worldLine_] := Function[lambda,
  Module[
    {d=Dt[worldLine[lambda], lambda]}, 
    FramePosition[d/d[[1]]]
  ]
]

FrameAcceleration[worldLine_] := Function[lambda,
  Module[
    {
      d=Dt[worldLine[lambda], lambda],
      dv=Dt[FrameVelocity[worldLine][lambda], lambda]
    }, 
    dv/d[[1]]
  ]
]


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

LorenzBoost[{vx_, vy_, vz_}] := Module[
  { 
    v, gamma, nx, ny, nz
  },
  v = Sqrt[vx^2+vy^2+vz^2];
  gamma = 1/Sqrt[1 - v^2];
  { nx, ny, nz } = {vx, vy, vz}/v;
  {
   {       gamma,          -gamma  v nx,         - gamma  v ny,          -gamma v nz  },
   { -gamma v nx, 1 + (gamma - 1)  nx^2,     (gamma - 1) nx ny,     (gamma - 1) nx nz },
   { -gamma v ny,     (gamma - 1) ny nx, 1 + (gamma - 1)  ny^2,     (gamma - 1) ny nz },
   { -gamma v nz,     (gamma - 1) nz nx,     (gamma - 1) nz ny, 1 + (gamma - 1) nz^2  }
  }
]

MCRFBoost[worldLine_] := Function[tau, LorenzBoost[FrameVelocity[worldLine][tau]]]

End[]

EndPackage[]
  
