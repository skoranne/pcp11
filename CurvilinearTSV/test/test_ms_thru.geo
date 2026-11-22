SetFactory("OpenCASCADE");
// Microstrip through, 2-port. Units: um, L0=1e-6
L=60.0; W=4.0; T=2.0; H=10.0; GndT=2.0; box=40.0; air=30.0; gap=2.0;
Box(1) = {-box,-box,-GndT, 2*box,2*box, H+air+GndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,H}; // substrate eps=4.4
Box(3) = {-box,-box,-GndT, 2*box,2*box,GndT}; // bottom ground
// Signal trace split into 3: left feed, center, right feed with gaps for ports
Box(4) = {-L/2,-W/2,H, gap, W, T};       // port 1 gap
Box(5) = {-L/2+gap,-W/2,H, L-2*gap, W, T}; // signal trace
Box(6) = {L/2-gap,-W/2,H, gap, W, T};     // port 2 gap

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6}; Delete; };
airv()={}; subv()={}; gndv()={}; metal()={}; p1v()={}; p2v()={};
For i In {0:#v()-1}
  bb()=BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin; cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  If (zmax <= 0.1)
    gndv() += {v(i)};
  EndIf
  If ((zmin >= -0.1) * (zmax <= H+0.1) * (dz > 1))
    subv() += {v(i)};
  EndIf
  If (zmin >= H-0.1)
    If (Abs(cx - (-L/2+gap/2)) < 1.5) p1v() += {v(i)}; EndIf
    If (Abs(cx - (L/2-gap/2)) < 1.5) p2v() += {v(i)}; EndIf
    If ((Abs(cx - (-L/2+gap/2)) >= 1.5) * (Abs(cx - (L/2-gap/2)) >= 1.5)) metal() += {v(i)}; EndIf
  EndIf
  If ((zmax > 0.1) * (!((zmin >= -0.1) * (zmax <= H+0.1) * (dz > 1))) * (!(zmin >= H-0.1)))
    airv() += {v(i)};
  EndIf
EndFor
If(#airv()>0) Physical Volume(1)={airv()}; EndIf
If(#subv()>0) Physical Volume(2)={subv()}; EndIf
If(#gndv()>0) Physical Volume(3)={gndv()}; EndIf
If(#metal()>0) Physical Volume(4)={metal()}; EndIf
If(#p1v()>0) Physical Volume(5)={p1v()}; EndIf
If(#p2v()>0) Physical Volume(6)={p2v()}; EndIf

allv()={}; allv()+=airv(); allv()+=subv(); allv()+=gndv(); allv()+=metal(); allv()+=p1v(); allv()+=p2v();
bnd()={}; For i In {0:#allv()-1} bnd()+=Boundary{Volume{allv(i)};}; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.3;
For i In {0:#bnd()-1}
  s=Abs(bnd(i)); bb()=BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter = (xmin<-box+eps)+(xmax>box-eps)+(ymin<-box+eps)+(ymax>box-eps)+(zmin<-GndT-eps)+(zmax>H+air-eps); onOuter=(onOuter>0);
  isZFace = (dz < 0.5);
  onPortY = (Abs(cy) <= W/2+eps);
  isP1 = isZFace * (Abs(cz - H) < 0.6) * (xmin >= -L/2-eps) * (xmax <= -L/2+gap+eps) * onPortY;
  isP2 = isZFace * (Abs(cz - H) < 0.6) * (xmin >= L/2-gap-eps) * (xmax <= L/2+eps) * onPortY;
  If (isP1) p1() += {s}; EndIf
  If (isP2) p2() += {s}; EndIf
  If (onOuter * (!isP1) * (!isP2)) absb() += {s}; EndIf
EndFor
For i In {0:#gndv()-1}
  gb()=Boundary{Volume{gndv(i)};};
  For j In {0:#gb()-1}
    gs=Abs(gb(j)); bb()=BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    onOuter=(xmin<-box+eps)+(xmax>box-eps)+(ymin<-box+eps)+(ymax>box-eps)+(zmin<-GndT-eps)+(zmax>H+air-eps); onOuter=(onOuter>0);
    If(!onOuter) pec()+={gs}; EndIf
  EndFor
EndFor
If(#absb()>0) Physical Surface(10)={absb()}; EndIf
If(#pec()>0) Physical Surface(11)={pec()}; EndIf
If(#p1()>0) Physical Surface(12)={p1()}; EndIf
If(#p2()>0) Physical Surface(13)={p2()}; EndIf
Mesh.CharacteristicLengthMin=1.5; Mesh.CharacteristicLengthMax=8; Mesh.Algorithm3D=10;
