////////////////////////////////////////////////////////////////////////////////
// File   : mesh_parser.cpp
// Author : Sandeep Koranne (C) 2018. All rights reserved.
// Purpose: GMSH mesh file parser for conversion.
//
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cstdlib>

#if 0
cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface{6};
Recombine Surface {6};
Mesh.SecondOrderIncomplete=0;

$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
9
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0.499999999998694 0 0
6 1 0.499999999998694 0
7 0.5000000000020591 1 0
8 0 0.5000000000020591 0
9 0.5000000000003766 0.5000000000003767 0
$EndNodes
$Elements
16
1 15 2 0 1 1
2 15 2 0 2 2
3 15 2 0 3 3
4 15 2 0 4 4
5 1 2 0 1 1 5
6 1 2 0 1 5 2
7 1 2 0 2 2 6
8 1 2 0 2 6 3
9 1 2 0 3 3 7
10 1 2 0 3 7 4
11 1 2 0 4 4 8
12 1 2 0 4 8 1
13 3 2 0 6 4 8 9 7
14 3 2 0 6 7 9 6 3
15 3 2 0 6 8 1 5 9
16 3 2 0 6 9 5 2 6
$EndElements

#endif

namespace MESH {
  struct Point
  {
    double x,y,z;
    Point( double ix=0, double iy=0, double iz=0 ): x( ix ), y( iy ), z( iz ) {}
  };
  typedef std::vector< Point > PointVector;

  struct Element
  {
    int id, physical_id, geometry_id;
    std::vector< int > pvec;
    void PrintElement( const PointVector& point_list );
  };
  typedef std::vector< Element > ElementVector;
  struct Mesh
  {
    Mesh() {}
    void ReadPoints( int N, std::istream& );
    void ReadElements( int N, std::istream& );
    static int GetNumberPoints(int id);
    void PrintMesh(std::ostream& COORD, std::ostream& E3, std::ostream& E4 );
    PointVector pvec;
    ElementVector evec;
  };
}

int MESH::Mesh::GetNumberPoints(int id)
{
  switch( id ) {
  case 15: return 1; // Point
  case 1:  return 2; // Line
  case 2:  return 3; // Triangle
  case 3:  return 4; // Quadrilateral
  case 4:  return 4; // Tetrahedron
  case 5:  return 8; // Hexahedron
  case 6:  return 6; // Prism
  case 7:  return 5; // Pyramid
  case 8:  return 3; // Second order line
  case 9:  return 6; // Seconrd order triangle
  case 11: return 10;// Second order tetrahedron
  default: return -1;
  }
  return -1;
}

void MESH::Element::PrintElement( const std::vector<Point>& point_list )
{
  return;
  std::cout << "Id = " << id << " Point = " << pvec.size() << std::endl;
  for( size_t i=0; i < pvec.size(); ++i ) {
    
  }
}

void MESH::Mesh::PrintMesh( std::ostream& COORD, std::ostream& E3, std::ostream& E4)
{
  // first print the coordinates
  COORD << "% Input-file for vertices generated from MESH\n";
  COORD << "% Node-number X Y\n";
  for( size_t i=1; i < pvec.size(); ++i ) {
    COORD << i << "\t" << pvec[i].x << "\t" << pvec[i].y << "\n";
  }
  std::cout << "Elements = " << evec.size() << std::endl;
  E3 << "% Input-file of triangles generated from MESH file.\n";
  E3 << "% Element-number / 1-node / 2-node/ 3-node\n";
  E4 << "% Input-file of parallelograms generated from MESH file.\n";
  E4 << "% Element-number / 1-node / 2-node/ 3-node / 4-node\n";
  
  size_t TRIANGLE_COUNT = 1, QUAD_COUNT = 1;
  for( size_t i=0; i < evec.size(); ++i ) {
    Element &E( evec[i] );
    if( ( E.id != 2 ) && ( E.id != 3 ) ) continue;
    E.PrintElement( pvec );
    if( E.id == 2 ) { // TRIANGLE
      assert( E.pvec.size() == 3 );
      E3 << TRIANGLE_COUNT++ << "\t" << E.pvec[0] << "\t" << E.pvec[1] << "\t" << E.pvec[2] << "\n";
    }
    else if( E.id == 3 ) { // QUADRILATERAL
      assert( E.pvec.size() == 4 );
      E4 << QUAD_COUNT++ << "\t" << E.pvec[0] << "\t" << E.pvec[1] << "\t" << E.pvec[2] << "\t" << E.pvec[3] << "\n";
    }
  }

  // Generate Neumann condition on the left side
  std::ofstream NEUMANN("neumann.dat");
  size_t NEUMANN_COUNT = 1;
  std::ofstream DIRICHLET("dirichlet.dat");
  size_t DIRICHLET_COUNT = 1;
  NEUMANN << "% Input-file for Neumann BC generated from MESH.\n";
  NEUMANN << "% Neumann-Edge-count / 1-Node / 2-Node\n";
  DIRICHLET << "% Input-file for Dirichlet BC generated from MESH.\n";
  DIRICHLET << "% Dirichlet-Edge-count / 1-Node / 2-Node\n";
  
  for( size_t i=1; i < evec.size(); ++i ) {
    Element &E( evec[i] );
    if( E.id != 1 ) continue;
    assert( E.pvec.size() == 2 ); // since this is a line
    #if 1
    if( pvec[E.pvec[0]].x == -1 || pvec[E.pvec[0]].x == 1 ||
	pvec[E.pvec[0]].y == -1 || pvec[E.pvec[0]].y == 1 ) {
    #endif
    #if 0
    if( pvec[E.pvec[0]].x == 0 || pvec[E.pvec[0]].x == 1 ||
	pvec[E.pvec[0]].y == 0 || pvec[E.pvec[0]].y == 1 ) {
    #endif
      DIRICHLET << DIRICHLET_COUNT++ << "\t" << E.pvec[0] << "\t" << E.pvec[1] << "\n";
    } else { //
      NEUMANN << NEUMANN_COUNT++ << "\t" << E.pvec[0] << "\t" << E.pvec[1] << "\n";      
    }      
  }
      
}


using namespace MESH;

void Mesh::ReadPoints( int N, std::istream& is )
{
  pvec.resize( N+1 );
  std::string line;
  for( int i=0; i < N; ++i ) {
    std::getline( is, line );
    std::stringstream sstr( line );
    int id;
    sstr >> id >> pvec[id].x >> pvec[id].y >> pvec[id].z;
    assert( ( id == i+1 ) && "Point id mismatch." );
  }
  std::getline( is, line );
  assert( line == "$EndNodes");
}

void Mesh::ReadElements( int N, std::istream& is )
{
  evec.resize( N+1 );
  std::string line;
  for( int i=0; i < N; ++i ) {
    std::getline( is, line );
    std::stringstream sstr( line );
    int id;
    int num_tags = 0;
    sstr >> id >> evec[id].id >> num_tags >> evec[id].geometry_id >> evec[id].physical_id;
    assert( num_tags == 2 );
    // now depending on the type of element we have to read the point list
    int num_points = Mesh::GetNumberPoints( evec[id].id );
    evec[id].pvec.resize( num_points );
    for( int j=0; j < num_points; ++j ) sstr >> evec[id].pvec[j];
    assert( ( id == i+1 ) && "Element id mismatch." );
  }
  std::getline( is, line );
  assert( line == "$EndElements");
}


static void ParseMesh( std::istream& is, Mesh& msh )
{
  std::string line;
  unsigned long line_number = 0;
  while( is ) {
    std::getline( is, line );
    std::cout << "L:" << line_number << "\t" << line << "\n";
    if( line_number == 0 ) assert( line == "$MeshFormat" );
    std::stringstream sstr( line );
    std::string token;
    sstr >> token;
    if( token == "$MeshFormat" ) {
      std::getline( is, line );
      std::getline( is, line );      
      line_number+= 3;
      continue;
    }
    if( token == "$Nodes" ) {
      std::getline( is, line );
      line_number += 2;
      int number_points = atoi( line.c_str() );
      line_number += number_points+1;
      std::cout << "Will read " << number_points << " points.\n";
      msh.ReadPoints( number_points, is );
    }
    if( token == "$Elements") {
      std::getline( is, line );
      line_number += 2;
      int number_elements = atoi( line.c_str() );
      line_number += number_elements+1;
      std::cout << "Will read " << number_elements << " elements.\n";
      msh.ReadElements( number_elements, is );
    }      
  }
}


static void Usage()
{
  std::cout << "./mesh_parser <msh-file> <output-file>\n";
}


int main( int argc, char* argv[] )
{
  if( argc != 3 ) {
    Usage();
    exit(-1);
  }
  std::ifstream ifs( argv[1] );
  if( !ifs ) {
    std::cout << "Cannot open file: " << argv[1] << " for reading.\n";
    exit(-1);
  }
  Mesh msh;  
  {
    std::string coordinate_file( argv[2] ), element3( argv[2] ), element4( argv[2] );
    coordinate_file += "_coordinates.dat";
    element3        += "_element3.dat";
    element4        += "_element4.dat";    
    std::ofstream coord_ofs( coordinate_file.c_str() );
    std::ofstream element3_ofs( element3.c_str() );
    std::ofstream element4_ofs( element4.c_str() );    
    ParseMesh( ifs, msh );
    msh.PrintMesh( coord_ofs, element3_ofs, element4_ofs);    
  }
  std::cout << "Processed mesh written to " << argv[2] << std::endl;
  return 0;
}

