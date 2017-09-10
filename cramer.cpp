////////////////////////////////////////////////////////////////////////////////
// File    : cramer.cpp
// Author  : Sandeep Koranne (C) 2017
// Purpose : Teaching example for C++ for Advay
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <iomanip>

// The 2 equations are given to us as
// (ax + by = c) and
// (dx + ey = f)
// We solve them using Determinants

void PrettyPrintEqn( int cx, int cy, int rhs)
{
  std::cout << std::setw(3) << cx << "x + " 
	    << std::setw(3) << cy << "y = " 
	    << std::setw(3) << rhs << "\n";
}

void ReadEquation(int eqno, int &cx, int& cy, int &rhs)
{
  // reads the coefficient for X,Y and RHS
  std::cout << "From Eqn " << eqno << " enter coefficient for X:\n";
  std::cin >> cx;
  std::cout << "From Eqn " << eqno << " enter coefficient for Y:\n";
  std::cin >> cy;
  std::cout << "Enter RHS:\n";
  std::cin >> rhs;
}

bool SolveEquation( int cx1, int cy1, int rhs1,
		    int cx2, int cy2, int rhs2,
		    int& X, int& Y )
{
  // calculate the determinant as
  std::cout << "Calculating the determinant as:\n";
  std::cout << "| cx1   cy1 |\n"
	    << "|           |    = ( cx1*cy2 - cx2*cy1 )\n"
	    << "| cx2   cy2 |\n"
	    << "             \n";
  const int D = cx1*cy2 - cx2*cy1;

  std::cout << "|" << std::setw(3) << cx1 << "    " << std::setw(3) << cy1 << "|\n"
	    << "|          |   = " << D << "\n"
	    << "|" << std::setw(3) << cx2 << "    " << std::setw(3) << cy2 << "|\n\n";

  if( D == 0 ) {
    std::cerr << "Determinant is 0. No solution.\n";
    X = Y = 0;
    return false;
  }
  std::cout << "In the above equation the determinant is " << D << std::endl;
  std::cout << "Next we calculate the determinant after column substitution.\n";
  std::cout << "First replace the X column with RHS\n";
  std::cout << "| rhs1   cy1 |\n"
	    << "|           |    = ( rhs1*cy2 - rhs2*cy1 )\n"
	    << "| rhs2   cy2 |\n"
	    << "             \n";
  float fX = (rhs1*cy2 - rhs2*cy1);
  std::cout << "|" << std::setw(3) << rhs1 << "    " << std::setw(3) << cy1 << "|\n"
	    << "|          |   = " << fX << "\n"
	    << "|" << std::setw(3) << rhs2 << "    " << std::setw(3) << cy2 << "|\n\n";
  std::cout << "Then divide this determinant by the coefficient determinant " << D << "\n";

  std::cout << "Therefore X = " << fX << " / " << D << " = " << fX / D << "\n";
  fX = fX / D; 

  std::cout << "Next we calculate the determinant after column substitution.\n";
  std::cout << "First replace the Y column with RHS\n";
  std::cout << "| cx1   rhs1 |\n"
	    << "|           |    = ( cx1*rhs2 - cx2*rhs1 )\n"
	    << "| cx2   rhs2 |\n"
	    << "             \n";
  float fY = (cx1*rhs2 - cx2*rhs1);
  std::cout << "|" << std::setw(3) << cx1 << "    " << std::setw(3) << rhs1 << "|\n"
	    << "|          |   = " << fY << "\n"
	    << "|" << std::setw(3) << cx2 << "    " << std::setw(3) << rhs2 << "|\n\n";
  std::cout << "Then divide this determinant by the coefficient determinant " << D << "\n";
  std::cout << "Therefore Y = " << fY << " / " << D << " = " << fY / D << "\n";
  fY = fY / D; 
  X = ( fX );
  Y = ( fY );
  return true;
}

int main()
{
  std::cout << "Solving 2x2 linear equation using Cramer's rule.\n";
  int cx1, cy1, rhs1, cx2, cy2, rhs2;
  ReadEquation( 1, cx1, cy1, rhs1);
  ReadEquation( 2, cx2, cy2, rhs2);
  std::cout << "\nStarting to solve the system of equation...\n";
  PrettyPrintEqn( cx1, cy1, rhs1);
  PrettyPrintEqn( cx2, cy2, rhs2);
  int X, Y;
  bool status = SolveEquation( cx1, cy1, rhs1, cx2, cy2, rhs2, X, Y );
  if( status ) {
    std::cout << "Solution is X = " << X << " Y = " << Y << std::endl;
  }
  return 0;
}

  
