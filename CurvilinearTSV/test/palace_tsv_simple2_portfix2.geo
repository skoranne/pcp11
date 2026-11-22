SetFactory("OpenCASCADE");
// Simple 2-port TSV, trace-to-ground lumped ports. Units: um, L0=1e-6
tsv=10.0; tox=0.5; h=50.0; box=60.0; airTop=30.0; airBot=15.0; gndT=2.0;
traceW=6.0; traceT=2.0; padW=12.0; portL=2.0;

Box(1) = {-box,-box,-airBot-gndT, 2*box,2*box, h+airTop+airBot+gndT}; // air
Box(2) = {-box,-box,0, 2*box,2*box,h}; // Si
Box(3) = {-tsv/2-tox,-tsv/2-tox,0, tsv+2*tox,tsv+2*tox,h}; // oxide
Box(4) = {-tsv/2,-tsv/2,0, tsv,tsv,h}; // Cu via
Box(5) = {-box,-box,-gndT, 2*box,2*box,gndT}; // bottom ground
// Top trace: stop short of boundaries for port bridges
Box(6) = {-box+portL,-traceW/2,h, 2*box-2*portL,traceW,traceT}; // top trace
Box(7) = {-padW/2,-padW/2,-gndT-traceT, padW,padW,traceT}; // bottom pad
Box(8) = {-padW/2,-padW/2,-gndT, padW,padW,gndT}; // via-bottom connection
// Port bridge volumes (air/field region, will be in air) - force boundary split
Box(9)  = {box-portL,-traceW/2,0, portL,traceW,h+traceT};  // port1 bridge at +X
Box(10) = {-box,-traceW/2,-gndT-traceT, portL,traceW,gndT+traceT}; // port2 bridge at -X

v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,6,7,8,9,10}; Delete; };

air()={}; si()={}; ox()={}; cu()={}; gnd()={}; metal()={};
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  inSiZ=(zmin>=-0.1)*(zmax<=h+0.1)*(dz>10);
  inBotGndZ=(zmin>=-gndT-0.1)*(zmax<=0.1)*(dz<=gndT+0.5);
  inTopTrZ=(zmin>=h-0.1)*(zmax<=h+traceT+0.1)*(dz>0.5);
  inBotTrZ=(zmin>=-gndT-traceT-0.1)*(zmax<=-gndT+0.1)*(dz>0.5);
  isOuterXY=(dx>20)+(dy>20);
  inPadArea=(Abs(cx)<=padW/2+1)*(Abs(cy)<=padW/2+1)*(dx<=padW+2);
  inPortBridge = (dx<=portL+0.5) * (Abs(cy)<=traceW/2+1) * ((xmax>box-portL-0.5) + (xmin<-box+portL+0.5));
  If(inBotGndZ)
    If(inPadArea) metal()+={v(i)}; EndIf
    If(!inPadArea) gnd()+={v(i)}; EndIf
  EndIf
  If(inSiZ * (!isOuterXY))
    If(dx<=tsv+0.2 && dy<=tsv+0.2) cu()+={v(i)}; EndIf
    If(!(dx<=tsv+0.2 && dy<=tsv+0.2)) ox()+={v(i)}; EndIf
  EndIf
  If(inSiZ * isOuterXY) si()+={v(i)}; EndIf
  If(inTopTrZ + inBotTrZ) metal()+={v(i)}; EndIf
  If(!(inBotGndZ)*!(inSiZ)*!(inTopTrZ+inBotTrZ)) air()+={v(i)}; EndIf
EndFor

If(#air()>0) Physical Volume(1)={air()}; EndIf
If(#si()>0)  Physical Volume(2)={si()}; EndIf
If(#ox()>0)  Physical Volume(3)={ox()}; EndIf
If(#cu()>0)  Physical Volume(4)={cu()}; EndIf
If(#gnd()>0) Physical Volume(5)={gnd()}; EndIf
If(#metal()>0) Physical Volume(6)={metal()}; EndIf

allv()={}; allv()+=air(); allv()+=si(); allv()+=ox(); allv()+=cu(); allv()+=gnd(); allv()+=metal();
bnd()={}; For i In {0:#allv()-1} bnd()+=Boundary{Volume{allv(i)};}; EndFor
pec()={}; p1()={}; p2()={}; absb()={}; eps=0.4;
For i In {0:#bnd()-1}
  s=Abs(bnd(i)); bb()=BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  onOuter=(xmin<-box+eps)+(xmax>box-eps)+(ymin<-box+eps)+(ymax>box-eps)+(zmin<-airBot-gndT-traceT+eps)+(zmax>h+airTop+eps); onOuter=(onOuter>0);
  isXFace=(dx<0.5); inPortY=(Abs(cy)<=traceW/2+eps);
  inP1Z=(zmin>=-eps)*(zmax<=h+traceT+eps);
  inP2Z=(zmin>=-gndT-traceT-eps)*(zmax<=0+eps);
  isP1=isXFace*(xmax>box-eps)*inPortY*inP1Z;
  isP2=isXFace*(xmin<-box+eps)*inPortY*inP2Z;
  If(isP1) p1()+={s}; EndIf
  If(isP2) p2()+={s}; EndIf
  If(onOuter*(!isP1)*(!isP2)) absb()+={s}; EndIf
EndFor
For i In {0:#gnd()-1}
  gb()=Boundary{Volume{gnd(i)};};
  For j In {0:#gb()-1}
    gs=Abs(gb(j)); bb()=BoundingBox Surface{gs};
    xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
    cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2; dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
    onOuter=(xmin<-box+eps)+(xmax>box-eps)+(ymin<-box+eps)+(ymax>box-eps)+(zmin<-airBot-gndT-traceT+eps)+(zmax>h+airTop+eps); onOuter=(onOuter>0);
    isXFace=(dx<0.5); inPortY=(Abs(cy)<=traceW/2+eps);
    inP1Zs=(zmin>=-eps)*(zmax<=h+traceT+eps);
    inP2Zs=(zmin>=-gndT-traceT-eps)*(zmax<=0+eps);
    isP1s=isXFace*(xmax>box-eps)*inPortY*inP1Zs;
    isP2s=isXFace*(xmin<-box+eps)*inPortY*inP2Zs;
    If((!onOuter)*(!isP1s)*(!isP2s)) pec()+={gs}; EndIf
  EndFor
EndFor
If(#absb()>0) Physical Surface(10)={absb()}; EndIf
If(#pec()>0) Physical Surface(11)={pec()}; EndIf
If(#p1()>0) Physical Surface(12)={p1()}; EndIf
If(#p2()>0) Physical Surface(13)={p2()}; EndIf
Mesh.CharacteristicLengthMin=1.5; Mesh.CharacteristicLengthMax=10.0; Mesh.Algorithm3D=10;
Field[1]=Box; Field[1].VIn=1.0; Field[1].VOut=10; Field[1].XMin=-box-1; Field[1].XMax=-box+5; Field[1].YMin=-traceW-2; Field[1].YMax=traceW+2; Field[1].ZMin=-gndT-traceT-2; Field[1].ZMax=h+traceT+2;
Field[2]=Box; Field[2].VIn=1.0; Field[2].VOut=10; Field[2].XMin=box-5; Field[2].XMax=box+1; Field[2].YMin=-traceW-2; Field[2].YMax=traceW+2; Field[2].ZMin=-gndT-2; Field[2].ZMax=h+traceT+2;
Field[3]=Min; Field[3].FieldsList={1,2};
Background Field=3;
