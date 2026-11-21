SetFactory("OpenCASCADE");

// Units: microns, Palace L0=1e-6
// Simple 2-port twin TSV (square approximation) through Si with SiO2 liner
tsv_w = 10.0;   // TSV square width
t_ox = 0.5;     // oxide thickness
h_si = 50.0;    // substrate height
pitch = 40.0;   // center-to-center
box = 60.0;     // air margin
air_top = 30.0;
air_bot = 10.0;
trace_w = 4.0;
trace_t = 2.0;
trace_l = 20.0;
port_gap = 2.0;
gnd_t = 2.0;

// Air domain
Box(1) = {-pitch/2-box, -box, -air_bot-gnd_t, pitch+2*box, 2*box, h_si+air_top+air_bot+gnd_t};
// Si substrate
Box(2) = {-pitch/2-box, -box, 0, pitch+2*box, 2*box, h_si};
// Oxide TSV1
Box(3) = {-pitch/2-tsv_w/2-t_ox, -tsv_w/2-t_ox, 0, tsv_w+2*t_ox, tsv_w+2*t_ox, h_si};
// Oxide TSV2
Box(4) = { pitch/2-tsv_w/2-t_ox, -tsv_w/2-t_ox, 0, tsv_w+2*t_ox, tsv_w+2*t_ox, h_si};
// Cu TSV1
Box(5) = {-pitch/2-tsv_w/2, -tsv_w/2, 0, tsv_w, tsv_w, h_si};
// Cu TSV2
Box(6) = { pitch/2-tsv_w/2, -tsv_w/2, 0, tsv_w, tsv_w, h_si};
// Bottom ground
Box(7) = {-pitch/2-box, -box, -gnd_t, pitch+2*box, 2*box, gnd_t};
// Top trace 1 (outward to -x)
Box(8) = {-pitch/2-tsv_w/2-trace_l, -trace_w/2, h_si, trace_l, trace_w, trace_t};
// Top trace 2 (outward to +x)
Box(9) = { pitch/2+tsv_w/2, -trace_w/2, h_si, trace_l, trace_w, trace_t};
// Ground return traces on top (connect to outer ground walls conceptually) - extend to box edge
Box(10) = {-pitch/2-box, -trace_w/2, h_si, box - trace_l - tsv_w/2, trace_w, trace_t};
Box(11) = { pitch/2+tsv_w/2+trace_l, -trace_w/2, h_si, box - trace_l - tsv_w/2, trace_w, trace_t};
// Port gap 1
Box(12) = {-pitch/2-tsv_w/2-trace_l-port_gap, -trace_w/2, h_si, port_gap, trace_w, trace_t};
// Port gap 2
Box(13) = { pitch/2+tsv_w/2+trace_l, -trace_w/2, h_si, port_gap, trace_w, trace_t};

// Fragment
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{2:13}; Delete; };

Physical Volume(1) = {v()};

Mesh.CharacteristicLengthMin = 2.0;
Mesh.CharacteristicLengthMax = 12.0;
Mesh.Algorithm3D = 10;
