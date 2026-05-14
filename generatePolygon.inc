// -------------------------------------------------------------------
// MACRO: BuildBox
// Generates a box dynamically from points, lines, loops, and surfaces
// Inputs required before calling:
// - x0, y0, dx, dy
// - z0: Base height
// - dz: Extrusion thickness
// - lc: Mesh characteristic length
// Return sl_id
// -------------------------------------------------------------------
Macro BuildBox
  // Base points
  p = newp;
  Point(p)   = {x0, y0, z0, lc};
  Point(p+1) = {x0+dx, y0, z0, lc};
  Point(p+2) = {x0+dx, y0+dy, z0, lc};
  Point(p+3) = {x0, y0+dy, z0, lc};
  Point(p+4) = {x0, y0, z0+dz, lc};
  Point(p+5) = {x0+dx, y0, z0+dz, lc};
  Point(p+6) = {x0+dx, y0+dy, z0+dz, lc};
  Point(p+7) = {x0, y0+dy, z0+dz, lc};

  // Connecting lines
  l = newl;
  Line(l)   = {p, p+1};     // Bottom edges
  Line(l+1) = {p+1, p+2};
  Line(l+2) = {p+2, p+3};
  Line(l+3) = {p+3, p};
  Line(l+4) = {p+4, p+5};   // Top edges
  Line(l+5) = {p+5, p+6};
  Line(l+6) = {p+6, p+7};
  Line(l+7) = {p+7, p+4};
  Line(l+8) = {p, p+4};     // Vertical edges
  Line(l+9) = {p+1, p+5};
  Line(l+10)= {p+2, p+6};
  Line(l+11)= {p+3, p+7};

  // Curve loops defining the 6 faces
  cl = newcl;
  Curve Loop(cl)   = {l, l+1, l+2, l+3};          // Bottom face
  Curve Loop(cl+1) = {l+4, l+5, l+6, l+7};        // Top face
  Curve Loop(cl+2) = {l, l+9, -(l+4), -(l+8)};    // Front face
  Curve Loop(cl+3) = {l+1, l+10, -(l+5), -(l+9)}; // Right face
  Curve Loop(cl+4) = {l+2, l+11, -(l+6), -(l+10)};// Back face
  Curve Loop(cl+5) = {l+3, l+8, -(l+7), -(l+11)}; // Left face

  // Plane surfaces mapping to the curve loops
  s = news;
  Plane Surface(s)   = {cl};
  Plane Surface(s+1) = {cl+1};
  Plane Surface(s+2) = {cl+2};
  Plane Surface(s+3) = {cl+3};
  Plane Surface(s+4) = {cl+4};
  Plane Surface(s+5) = {cl+5};

  // Final Surface Loop enclosing the geometry
  sl_id = newsl;
  Surface Loop(sl_id) = {s, s+1, s+2, s+3, s+4, s+5};

  // Save specific surface IDs so we can assign them to Physical Groups later
  s_bot = s; s_top = s+1; s_fnt = s+2; s_rgt = s+3; s_bck = s+4; s_lft = s+5;
  all_surfs[] = { s_bot, s_top, s_fnt, s_rgt, s_bck, s_lft };
  surfs_side = { s_fnt, s_rgt, s_bck, s_lft };
  For i In {0 : #surfs_side[] - 1}
     // Force a structured grid on the side walls
     Transfinite Surface {surfs_side[i]};
     // Recombine triangles into quads (optional, but often better for thin layers)
     Recombine Surface {surfs_side[i]}; 
  EndFor
Return

// -------------------------------------------------------------------
// MACRO: BuildPolygon
// Generates an N-sided 3D polygon from coordinate arrays.
// Inputs required before calling:
// - px[]: Array of X coordinates for the vertices (ordered sequentially)
// - py[]: Array of Y coordinates for the vertices (ordered sequentially)
// - z0: Base height
// - dz: Extrusion thickness
// - lc: Mesh characteristic length
// return value sl_id and all_surfs, both are VERY important. sl_id is jus
// topological.
// -------------------------------------------------------------------
Macro BuildPolygon
  N = #px[]; // Automatically get the number of vertices

  // Initialize empty arrays to store entity IDs
  pts_bot[] = {}; pts_top[] = {};
  lines_bot[] = {}; lines_top[] = {}; lines_vert[] = {};
  surfs_side[] = {};

  // 1. Generate Points (Bottom and Top)
  For i In {0 : N-1}
    pts_bot[i] = newp; Point(pts_bot[i]) = {px[i], py[i], z0, lc};
    pts_top[i] = newp; Point(pts_top[i]) = {px[i], py[i], z0+dz, lc};
  EndFor

  // 2. Generate Lines
  For i In {0 : N-1}
    // Determine the next index, wrapping back to 0 at the end of the array
    next_i = (i == N-1) ? 0 : i+1;

    lines_bot[i]  = newl; Line(lines_bot[i])  = {pts_bot[i], pts_bot[next_i]};
    lines_top[i]  = newl; Line(lines_top[i])  = {pts_top[i], pts_top[next_i]};
    lines_vert[i] = newl; Line(lines_vert[i]) = {pts_bot[i], pts_top[i]};
  EndFor

  // 3. Generate Top and Bottom Surfaces
  cl_bot = newcl; Curve Loop(cl_bot) = lines_bot[];
  s_bot = news;   Plane Surface(s_bot) = {cl_bot};

  cl_top = newcl; Curve Loop(cl_top) = lines_top[];
  s_top = news;   Plane Surface(s_top) = {cl_top};

  // 4. Generate Side Surfaces dynamically
  For i In {0 : N-1}
    next_i = (i == N-1) ? 0 : i+1;
    cl_side = newcl;
    // Loop orientation: Bottom line (forward), Right vertical (forward), 
    // Top line (backward), Left vertical (backward)
    Curve Loop(cl_side) = {lines_bot[i], lines_vert[next_i], -lines_top[i], -lines_vert[i]};
    surfs_side[i] = news; Plane Surface(surfs_side[i]) = {cl_side};
  EndFor

  // 5. Create Final Surface Loop
  sl_id = newsl;
  Surface Loop(sl_id) = {s_bot, s_top, surfs_side[]};
  For i In {0 : #surfs_side[] - 1}
      // Force a structured grid on the side walls
      Transfinite Surface {surfs_side[i]};
      // Recombine triangles into quads (optional, but often better for thin layers)
      Recombine Surface {surfs_side[i]}; 
  EndFor
  // Collect all generated surfaces into a single array so we can assign them easily
  all_surfs[] = {s_bot, s_top, surfs_side[]};
Return
