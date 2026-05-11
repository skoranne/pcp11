// File   : cap.pro
// Author : Sandeep Koranne (C) 2026. All rights reserved.
// Purpose: getdp control file to calculate V and E field for rectangular
//        : conductors immersed in dielectric, using OpenCASCADE kernel
//        : cap.pro -msh cs.msh -solve EleSta_v -pos Map -v 4 
current_pass = 1 ;
Group {
  // Matching the Physical IDs from the .geo file
  Rect1     = Region[10]        ;
  Rect2     = Region[11]        ;
  Rect3     = Region[12];   
  Diel      = Region[13]; 
  BndRect1 = Region[110]; 
  BndRect2 = Region[111]; 
  BndRect3 = Region[113];
  Vol_Ele = Region[ {Diel, Rect1, Rect2, Rect3} ];
  All_Rects = Region [ {Rect1, Rect2, Rect3} ];
}

Function {
  eps0 = 8.854187e-12;
  epsilon[Diel] = 1.0 * eps0;
  epsilon[Rect1] = 1.0 * eps0;
  epsilon[Rect2] = 1.0 * eps0;
  epsilon[Rect3] = 1.0 * eps0;  
  
  v_rect1 = 0.0;
  v_rect2 = 1.0;
  v_rect3 = 0.0;
  
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
          { GeoElement Quadrangle; NumberOfPoints 4; }
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

  
