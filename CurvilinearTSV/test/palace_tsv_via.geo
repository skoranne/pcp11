SetFactory("OpenCASCADE");

// Units: microns (Palace L0 = 1e-6)
r_cu = 5.0;     // Cu TSV radius
t_ox = 0.2;     // SiO2 liner thickness
r_ox = r_cu + t_ox;
h_si = 50.0;    // substrate height
box = 80.0;     // lateral box half-size
air_top = 40.0; // air above
air_bot = 20.0; // air below
port_h = 2.0;   // port annulus height
port_r_out = 12.0; // port outer radius (ground ring)

// Volumes
// Air + substrate container
Box(1) = {-box, -box, -air_bot, 2*box, 2*box, h_si + air_top + air_bot};

// Silicon substrate
Box(2) = {-box, -box, 0, 2*box, 2*box, h_si};

// SiO2 liner
Cylinder(3) = {0, 0, 0, 0, 0, h_si, r_ox};

// Cu TSV
Cylinder(4) = {0, 0, 0, 0, 0, h_si, r_cu};

// Top port annulus (air gap ring between TSV and ground ring) - as a thin cylinder shell
// We'll model port as a surface; create a small annular volume for tagging
Cylinder(5) = {0, 0, h_si, 0, 0, port_h, port_r_out};
Cylinder(6) = {0, 0, h_si, 0, 0, port_h, r_cu}; // subtract TSV

// Bottom port annulus
Cylinder(7) = {0, 0, -port_h, 0, 0, port_h, port_r_out};
Cylinder(8) = {0, 0, -port_h, 0, 0, port_h, r_cu};

BooleanDifference{ Volume{5}; Delete; }{ Volume{6}; Delete; }
BooleanDifference{ Volume{7}; Delete; }{ Volume{8}; Delete; }

// Fragment all
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2,3,4,5,7}; Delete; };

// Classify volumes
airV() = {}; siV() = {}; oxV() = {}; cuV() = {}; portTopV() = {}; portBotV() = {};
For i In {0 : #v()-1}
  bb() = BoundingBox Volume{v(i)};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;

  isAirBox = (dx > 150 && dy > 150);
  isSi = (!isAirBox && dz > 45 && dz < 55 && dx > 150);
  isCu = (!isAirBox && !isSi && Max(dx,dy) < 12 && cz > -1 && cz < h_si+1 && dz > 45);
  isOx = (!isAirBox && !isSi && !isCu && Max(dx,dy) < 14 && dz > 45 && cz > -1 && cz < h_si+1);
  isPortTop = (!isAirBox && cz > h_si-1 && dz <= port_h+1 && Max(dx,dy) < 30);
  isPortBot = (!isAirBox && cz < 1 && zmin < 0 && dz <= port_h+1 && Max(dx,dy) < 30);

  If (isAirBox)
    airV() += {v(i)};
  ElseIf (isSi)
    siV() += {v(i)};
  ElseIf (isCu)
    cuV() += {v(i)};
  ElseIf (isOx)
    oxV() += {v(i)};
  ElseIf (isPortTop)
    portTopV() += {v(i)};
  ElseIf (isPortBot)
    portBotV() += {v(i)};
  Else
    If (cz > h_si/2) airV() += {v(i)};
    Else siV() += {v(i)};
    EndIf
  EndIf
EndFor

If (#airV() > 0) Physical Volume(1) = {airV()}; EndIf
If (#siV() > 0)  Physical Volume(2) = {siV()}; EndIf
If (#oxV() > 0)  Physical Volume(3) = {oxV()}; EndIf
If (#cuV() > 0)  Physical Volume(4) = {cuV()}; EndIf
// port volumes merged into air for EM; we'll tag their bottom/top surfaces

// Boundary surfaces
allV() = {}; allV() += airV(); allV() += siV(); allV() += oxV(); allV() += cuV(); allV() += portTopV(); allV() += portBotV();
bnd() = {};
For i In {0 : #allV()-1}
  bnd() += Boundary{ Volume{allV(i)}; };
EndFor

pec() = {}; port1() = {}; port2() = {}; absb() = {};
eps = 0.3;
For i In {0 : #bnd()-1}
  s = Abs(bnd(i));
  bb() = BoundingBox Surface{s};
  xmin=bb(0); ymin=bb(1); zmin=bb(2); xmax=bb(3); ymax=bb(4); zmax=bb(5);
  cx=(xmin+xmax)/2; cy=(ymin+ymax)/2; cz=(zmin+zmax)/2;
  dx=xmax-xmin; dy=ymax-ymin; dz=zmax-zmin;

  onOuter = (xmin < -box+eps || xmax > box-eps || ymin < -box+eps || ymax > box-eps || zmin < -air_bot+eps || zmax > h_si+air_top-eps);
  // Port 1: top annular surface at z = h_si, r between r_cu and port_r_out
  rmax = Sqrt(Max(xmax*xmax+ymax*ymax, xmin*xmin+ymin*ymin));
  rmin_pt = Sqrt(Min(xmin*xmin+ymin*ymin, xmax*xmax+ymax*ymax));
  isPort1 = (Abs(zmin-h_si) < eps && Abs(zmax-h_si) < eps && rmax <= port_r_out+eps && rmax > r_cu);
  isPort2 = (Abs(zmin+port_h) < eps && Abs(zmax+port_h) < eps && rmax <= port_r_out+eps && rmax > r_cu);

  If (onOuter)
    absb() += {s};
  ElseIf (isPort1)
    port1() += {s};
  ElseIf (isPort2)
    port2() += {s};
  Else
    pec() += {s};
  EndIf
EndFor

Physical Surface(10) = {absb()};
Physical Surface(11) = {pec()};
If (#port1() > 0) Physical Surface(12) = {port1()}; EndIf
If (#port2() > 0) Physical Surface(13) = {port2()}; EndIf

Mesh.CharacteristicLengthMin = 1.0;
Mesh.CharacteristicLengthMax = 10.0;
Mesh.Algorithm3D = 10;
