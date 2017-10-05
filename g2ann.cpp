////////////////////////////////////////////////////////////////////////////////
// File    : g2ann.cpp
// Author  : Sandeep Koranne (C) 2017
// Purpose : Convert a graph (Petri-Net) or ANN
//         :
//         : The training part is not described here, this is just a way
//         : to convert a DAG into computable C++ representation
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <string>
#include <sstream>

int main( int argc, char* argv[])
{
  std::string line;
  while( std::getline( std::cin, line ) ) {
    if( line == "" ) break;
    //std::cout << "Processing: " << str << std::endl;
    std::string token;
    std::stringstream stream( line );
    stream >> token;
    if( token == "" ) break;
    if( token == "#" ) continue;
    if( token == "INORDER") continue;
    if( token == "OUTORDER") continue;
    //if( ( token[0] != 'n' ) && ( token[0] != 'O' ) ) continue; // INORDER
    std::string LHS, RHS1, RHS2;
    bool neg1 = false, neg2 = false;
    LHS = token;
    if( LHS == "" ) break;
    stream >> token;
    if( token != "=" ) continue; // not an assignment line
    stream >> RHS1;
    if( RHS1[0] == '!' ) { neg1 = true; RHS1 = RHS1.substr(1); }
    stream >> token;
    assert( token == "*" );
    stream >> RHS2;
    // usually RHS2 has the trailing ;, get rid of it
    RHS2 = RHS2.substr( 0, RHS2.size()-1);
    if( RHS2[0] == '!' ) { neg2 = true; RHS2 = RHS2.substr(1); }
    
    std::cout << "Node processing. " 
	      << LHS << " = " << ( neg1 ? "!" : "" ) << RHS1
	      << " * " << ( neg2 ? "!" : "" ) << RHS2 << "\n";
  }
  return 0;
}

