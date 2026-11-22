SetFactory("OpenCASCADE");
W=4.0; T=2.0; H=5.0; GndT=2.0; L=60.0; boxY=30.0; airTop=20.0;
x0=-L/2; x1=L/2;
Box(1) = {x0,-boxY,-GndT, L,2*boxY, H+airTop+GndT}; // air
Box(2) = {x0,-boxY,0, L,2*boxY,H}; // substrate
Box(3) = {x0,-boxY,-GndT, L,2*boxY,GndT}; // ground
Box(4) = {x0,-W/2,H, L,W,T}; // trace
// Port marker boxes (thin, on boundary) to force surface split
Box(5) = {x0,-W/2,-GndT, 0.5,W,H+T+GndT}; // port1 region
Box(6) = {x1-0.5,-W/2,-GndT, 0.5,W,H+T+GndT}; // port2 region
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6}; Delete; };
airv()={}; subv()={}; gndv()={}; metal()={};
For i In {0:#v()-1}
  bb()=BoundingBox Volume{v(i)};
  xmin=bb(0); xmax=bb(3); ymin=bb(1); ymax=bb(4); zmin=bb(2); zmax=bb(5); dz=zmax-zmin;
  inSub=(zmin>=-0.1)*(zmax<=H+0.1)*(dz>1)*(Abs(xmax-xmin)>1);
  inGnd=(zmax <= 0.1)*(dz<=GndT+0.5)*(Abs(xmax-xmin)>1);
  inTr=(zmin>=H-0.1)*(zmax<=H+T+0.1)*(Abs((ymin+ymax)/2)<=W/2+1)*(Abs(xmax-xmin)>1);
  If (inGnd) gndv()+={v(i)}; EndIf
  If (inSub) subv()+={v(i)}; EndIf
  If (inTr) metal()+={v(i)}; EndIf
  If ((!inGnd)*(!inSub)*(!inTr)) airv()+={v(i)}; EndIf
EndFor
If(#airv()>0) Physical Volume(1)={airv()}; EndIf
If(#subv()>0) Physical Volume(2)={subv()}; EndIf
If(#gndv()>0) Physical Volume(3)={gndv()}; EndIf
If(#metal()>0) Physical Volume(4)={metal()}; EndIf
allv()={}; allv()+=airv(); allv()+=subv(); allv()+=gndv(); allv()+=metal();
bnd()={}; For i In {0:#allv()-1} bnd()+=Boundary{Volume{allv(i)};}; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.3;
For i In {0:#bnd()-1}
  s=Abs(bnd(i)); bb()=BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onXmin=(xmin<x0+eps); onXmax=(xmax>x1-eps);
  onY=(ymin<-boxY+eps)+(ymax>boxY-eps); onZ=(zmin<-GndT-eps)+(zmax>H+airTop-eps);
  onOuter=(onXmin+onXmax+(onY>0)+(onZ>0))>0;
  isXFace=(dx<0.5); inPortY=(Abs(cy)<=W/2+eps); inPortZ=(zmin>=-GndT-eps)*(zmax<=H+T+eps);
  isP1=isXFace*onXmin*inPortY*inPortZ; isP2=isXFace*onXmax*inPortY*inPortZ;
  If(isP1) p1()+={s}; EndIf
  If(isP2) p2()+={s}; EndIf
  If(onOuter*(!isP1)*(!isP2)) absb()+={s}; EndIf
EndFor
For i In {0:#gndv()-1}
  gb()=Boundary{Volume{gndv(i)};};
  For j In {0:#gb()-1}
    gs=Abs(gb(j)); bb()=BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
    onXmin=(xmin<x0+eps); onXmax=(xmax>x1-eps);
    onY=(ymin<-boxY+eps)+(ymax>boxY-eps); onZ=(zmin<-GndT-eps)+(zmax>H+airTop-eps);
    onOuter=(onXmin+onXmax+(onY>0)+(onZ>0))>0;
    isXFace=(dx<0.5); inPortY=(Abs(cy)<=W/2+eps); inPortZ=(zmin>=-GndT-eps)*(zmax<=H+T+eps);
    isP1=isXFace*onXmin*inPortY*inPortZ; isP2=isXFace*onXmax*inPortY*inPortZ;
    If((!onOuter)*(!isP1)*(!isP2)) pec()+={gs}; EndIf
  EndFor
EndFor
If(#absb()>0) Physical Surface(10)={absb()}; EndIf
If(#pec()>0) Physical Surface(11)={pec()}; EndIf
If(#p1()>0) Physical Surface(12)={p1()}; EndIf
If(#p2()>0) Physical Surface(13)={p2()}; EndIf
Mesh.CharacteristicLengthMin=1.0; Mesh.CharacteristicLengthMax=6.0; Mesh.Algorithm3D=10;
