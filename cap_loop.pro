// File   : cap.pro
// Author : Sandeep Koranne (C) 2026. All rights reserved.
// Purpose: getdp control file to calculate V and E field for rectangular
//        : conductors immersed in dielectric, using OpenCASCADE kernel
//        : cap.pro -msh cs.msh -solve EleSta_v -pos Map -v 4 
current_pass = 1 ;
Group {
  // Matching the Physical IDs from the .geo file
  Rect1     = Region[ {2:800} ]   ;
  Rect2     = Region[ 801 ];
  Rect3     = Region[ {802:1000} ];
  Diel      = Region[ 1 ];
  Vol_Ele = Region[ {Diel} ];  
}

Function {
  eps0 = 8.854187e-12;
  epsilon[] = ( Z[] < 1376.1 ) ? 3.9 : 4.5;
  v_rect1 = 0.0;
  v_rect2 = 1.0;
  v_rect3 = 0.0;
  //epsilon[AllRect] = 1.0 * eps0;
}

Constraint {
  { Name Dirichlet_Ele; Type Assign;
    Case {
       { Region Rect1           ; Value v_rect1; }
       { Region Rect2           ; Value v_rect2; }
       { Region Rect3           ; Value v_rect3; }
    }
  }
}

FunctionSpace {
  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Vol_Ele; Entity NodesOf[ Vol_Ele ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_Ele; }
    }
  }
}

Jacobian {
	 { Name Vol                      ;
	 Case {
//              { Region AllRect; Jacobian Sur; } 
	      { Region Vol_Ele                    ; Jacobian Vol; } }
	 }
}

Integration {
  { Name Int; 
    Case { 
      { Type Gauss; 
        Case { 
          // 3D Elements
          { GeoElement Tetrahedron ; NumberOfPoints 4 ; }
          { GeoElement Hexahedron  ; NumberOfPoints 8 ; }
          { GeoElement Prism       ; NumberOfPoints 6 ; }
          { GeoElement Pyramid     ; NumberOfPoints 8 ; } // <-- The fix

          // 2D & 1D Elements (for your boundary conditions)
          { GeoElement Triangle    ; NumberOfPoints 3 ; }
          { GeoElement Quadrangle  ; NumberOfPoints 4 ; }
          { GeoElement Line        ; NumberOfPoints 2 ; }
        }
      }
    }
  }
}

Formulation {
  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
//      { Name q; Type Global; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      // Standard Laplace: (epsilon * grad v, grad v')
        Integral { [ epsilon[] * Dof{d v} , {d v} ];
        In Vol_Ele; Jacobian Vol; Integration Int; }
	//GlobalTerm { [ Dof{q}, {v} ]; In All_Rects; }
    }
  }
}

Resolution {
  { Name EleSta_v;
    System {
      { Name Sys; NameOfFormulation Electrostatics_v; }
    }
    Operation {
          Generate[Sys]; 
          Solve[Sys]; 
          PostOperation[Map];

    }
  }
}

PostProcessing {
  { Name EleSta_v; NameOfFormulation Electrostatics_v;
    Quantity {
      { Name v; Value { Term { [ {v} ]; In Vol_Ele; Jacobian Vol; } } }
      { Name e; Value { Term { [ -{d v} ]; In Vol_Ele; Jacobian Vol; } } }
      // Capacitance via Energy: W = 1/2 * integral( epsilon * |E|^2 )
      { Name Energy; Value { 
          Integral { [ 0.5 * epsilon[] * Norm[-{d v}]^2 ]; 
          In Vol_Ele; Jacobian Vol; Integration Int; } 
        } 
      }
    }
  }
}


PostOperation {
  { Name Map; NameOfPostProcessing EleSta_v;
    Operation {
      Print[ v, OnElementsOf Vol_Ele, File "v.pos" ];
      Print[ e, OnElementsOf Vol_Ele, File "e.pos" ];
      // C = 2 * Energy (since V = 1)
      Print[ Energy[Vol_Ele], OnGlobal, Format Table, File "energy.txt" ];
    }
  }
}

PostProcessing {
  { Name DebugCoords; NameOfFormulation Electrostatics_v;
    Quantity {
      // 1. Wrap the Z[] function into a named Quantity
      { Name My_Z_Coord; Value { Term { [ Z[] ]; In Vol_Ele; } } }
    }
  }
}
  
PostOperation {
  { Name Print_Z; NameOfPostProcessing DebugCoords;
    Operation {
      Print[ My_Z_Coord, OnElementsOf Vol_Ele, File "z.txt" ];
    }
  }
}
