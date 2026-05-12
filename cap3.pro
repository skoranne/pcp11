// File   : cap3.pro
// Author : Sandeep Koranne (C) 2026. All rights reserved.
// Purpose: getdp control file to calculate V and E field for rectangular
//        : conductors immersed in dielectric, using OpenCASCADE kernel
//        : cap.pro -msh cs.msh -solve EleSta_v -pos Map -v 4 
current_pass = 1 ;
Group {
  // Matching the Physical IDs from the .geo file
  AllDiel      = Region[503:507]  ;
  SubDiel      = Region[503];
  FOXKDiel     = Region[504]     ;
  PSGDiel      = Region[505]     ;
  NILD2Diel    = Region[506];        
  NILD3Diel    = Region[507]     ;
  NILD4Diel    = Region[508]     ;
  NILD5Diel    = Region[509];    
  NILD6Diel    = Region[510]     ;
  PI1KDiel     = Region[511];    
  TargetRegion = Region[ {4:6,10:12,NILD2Diel} ];  
  OneV      = Region[618]       ;
  GND1       = Region[602:617]  ;  
  GND2       = Region[617:668]  ;
  All_Rects = Region [ 133:142 ];  
  Vol_Ele = Region[ {AllDiel, All_Rects} ];

}

Function {
  eps0 = 8.854187e-12           ;
  epsilon[SubDiel] = 3.9 * eps0 ;
  epsilon[FOXKDiel] = 3.9 * eps0 ;
  epsilon[PSGDiel]  = 3.9 * eps0 ;
  epsilon[NILD2Diel]= 4.05 * eps0 ;
  epsilon[NILD3Diel]= 4.5 * eps0 ;
  epsilon[NILD4Diel]= 4.2 * eps0 ;
  epsilon[NILD5Diel]= 4.1 * eps0 ;
  epsilon[NILD6Diel]= 4.0 * eps0 ;
  epsilon[PI1KDiel] = 2.94 * eps0 ;
  epsilon[All_Rects] = 1.0*eps0;
}

Constraint {
  { Name Dirichlet_Ele; Type Assign;
    Case {
       { Region OneV            ; Value 1.0; }
       { Region GND1             ; Value 0.0; }
       { Region GND2           ; Value 0.0; }              
    }
  }
}

FunctionSpace {
  { Name Hgrad_v_Ele; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Vol_Ele; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_Ele; }
    }
  }
}

Jacobian {
  { Name Vol; Case { { Region All; Jacobian Vol; } } }
}

Integration {
  { Name Int;
    Case { 
      { Type Gauss;
        Case { 
          { GeoElement Triangle;   NumberOfPoints 4; }
          { GeoElement Quadrangle ; NumberOfPoints 4; }
          { GeoElement Tetrahedron     ; NumberOfPoints 4; }
          { GeoElement Hexahedron      ;  NumberOfPoints 8; } 
        }
      }
    }
  }
}

Formulation {
  { Name Electrostatics_v; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_v_Ele; }
    }
    Equation {
      // Standard Laplace: (epsilon * grad v, grad v')
        Integral { [ epsilon[] * Dof{d v} , {d v} ];
        In Vol_Ele; Jacobian Vol; Integration Int; }
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

PostOperation {
  { Name CrossSectionMap; NameOfPostProcessing EleSta_v;
    Operation {
    Print[ v, OnElementsOf TargetRegion, File "vtarget.pos" ]; 
    Print[ e, OnElementsOf Vol_Ele, File "e.pos" ]; 
    }
  }
}

  
