////////////////////////////////////////////////////////////////////////////////
// File    : function_arrangements.cpp
// Author  : Sandeep Koranne (C) 2019. All rights reserved.
//
// Purpose : Understanding code logic for functions
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <vector>
#include <deque>
#include <functional>
#include <algorithm>
#include <iostream>
#include <string>
#include <typeinfo>
#include <type_traits>
#include <cmath>
#include <cstddef> // for std::byte in C++-17
#include <memory>
//#include <ranges> // We are using C++-17 so cannot use ranges yet

namespace CodeLogic
{
  class RelatedConcepts
  {
  public:
    RelatedConcepts() = delete;
    template <typename T>
    static void InitializeDataStore(T&);
  };
} // end of namespace CodeLogic

namespace DataStore
{
  enum class TypeStorage { INTEGER, CHAR, FLOAT };
  class IntegerStore
  {
  public:
    static constexpr bool is_a_data_store = true;
    static constexpr TypeStorage type = TypeStorage::INTEGER;
    IntegerStore() {}
    void Init() { _data.clear(); }
    IntegerStore( IntegerStore&& rhs ) : _data( rhs._data ) {
      std::cout << __PRETTY_FUNCTION__ << std::endl;
    }
    IntegerStore& operator=( IntegerStore&& rhs ) {
      if( this == &rhs ) return *this;
      std::cout << __PRETTY_FUNCTION__ << ":" << __LINE__ << std::endl;
      _data = std::move( rhs._data );
      return (*this);
    }
    IntegerStore( const IntegerStore& ) = delete;
    IntegerStore& operator=( const IntegerStore& ) = delete;
  private:
    std::deque<int> _data;
  };

  class FloatStore
  {
  public:
    static constexpr bool is_a_data_store = true;
    static constexpr TypeStorage type = TypeStorage::FLOAT;
    FloatStore() {}
    void Init() {}
    FloatStore( FloatStore&& rhs ) {}
    FloatStore( const FloatStore& ) = delete;
    FloatStore& operator=( const FloatStore& ) = delete;
  };

} // end of DataStore namespace

namespace CommonDataTypes
{
  using Indices = std::vector<unsigned long>;
  using OneDMesh  = std::vector<float>;
  using TwoDMesh  = std::vector<OneDMesh>;
  enum class NORM { L1, L2, LINF };
  static void CheckCommonDataTypes();
  double ComputeNorm( const OneDMesh&, NORM n );
  double ComputeNorm( const TwoDMesh&, NORM n );
}

[[nodiscard]] double CommonDataTypes::ComputeNorm( const OneDMesh& mesh, NORM n )
{
  switch( n ) {
  case NORM::L2: {
    double retval = 0.0;
    for( unsigned int i=0; i < mesh.size(); ++i ) {
      retval += mesh[i]*mesh[i];
    }
    retval = std::sqrt( retval );
    return retval;
  }
  default: return 0.0;
  }
  return 0.0;
}

double CommonDataTypes::ComputeNorm( const TwoDMesh& mesh, NORM n )
{
  int X [[gnu::unused]] = 10;
  switch( n ) {
  [[likely]] case NORM::L2: {
    double retval = 0.0;
    for( unsigned int i=0; i < mesh.size(); ++i ) {
      for( unsigned int j=0; j < mesh[i].size(); ++j ) {
	retval += mesh[i][j] * mesh[i][j];
      }
    }
    retval = sqrt( retval );
    return retval;
  }
  default: return 0.0;
  }
  return 0.0;
}

namespace FEMCommon
{
  using namespace CommonDataTypes;

  class UnitSquare
  {
  private:
    static constexpr unsigned int N = 16;
  public:
    using MeshFunction = std::function<double(double,double)>;
    UnitSquare( unsigned int _MESH_SIZE = N): MESH_SIZE( _MESH_SIZE ) {
      Init( MESH_SIZE );
    }
    unsigned int GetMeshSize(void) const { return MESH_SIZE; }
    void Init( unsigned int );
    void InitializeMesh( MeshFunction );
    UnitSquare( const UnitSquare& ) = delete;
    UnitSquare& operator=( const UnitSquare& ) = delete;
    double ComputeNorm(NORM n=NORM::L2) const { return CommonDataTypes::ComputeNorm( mesh, n ); }
  private:
    unsigned int MESH_SIZE;
    TwoDMesh mesh;
  };
} // end of FEMCommon namespace
    
void FEMCommon::UnitSquare::InitializeMesh( MeshFunction F )
{
  const double Y_STEP = 1.0/MESH_SIZE;
  for( unsigned int i=0; i < MESH_SIZE; ++i ) {
    const double Y = i*Y_STEP;
    const double X_STEP = 1.0/mesh[i].size();
    for( unsigned int j=0; j < mesh[i].size(); ++j ) {
      const double X = j*X_STEP;
      mesh[i][j] = F(X,Y);
    }
  }
}


void FEMCommon::UnitSquare::Init( unsigned int N ) 
{
  mesh.clear();
  MESH_SIZE = N;
  mesh.resize( N );
  for( unsigned int i=0; i < N; ++i ) {
    mesh[i].resize( N );
    for( unsigned int j=0; j < N; ++j ) {
      mesh[i][j] = 0.0;
    }
  }
}

void CommonDataTypes::CheckCommonDataTypes()
{
  std::cout << "OneDMesh = " << typeid( OneDMesh ).name() << std::endl;
  std::cout << "sizeof( OneDMesh ) = " << sizeof( OneDMesh ) << std::endl;
}

template <typename T>
void CodeLogic::RelatedConcepts::InitializeDataStore( T& dataStore )
{
  static_assert( T::is_a_data_store, "Type is not a data-store." );
  dataStore.Init();
}

template <typename T>
auto CheckConformance( const T&  ) -> decltype( T::is_a_data_store )
{
  return ( T::is_a_data_store );
}

template <typename T>
std::string GetTypeString()
{
  return "UNKNOWN";
}

template <> std::string GetTypeString<bool>() { return "bool"; }

static void TestTypeID()
{
  int X = 10;
  std::cout << " X = " << X << " of type " << typeid( X ).name() << std::endl;
}

static void TestFEM()
{
  FEMCommon::UnitSquare SQ( 10 );
  double X = SQ.ComputeNorm();
  std::cout << "|SQ|_2 = " << X << std::endl;
  SQ.InitializeMesh( [](double X, double Y){ return X*Y; } );
  X = SQ.ComputeNorm();
  std::cout << "|SQ|_2 = " << X << std::endl;
}

long double operator "" _g( long double d ) { return d*0.001; }
long double operator "" _kg( long double d ) { return d;       }
long double operator "" _pounds( long double d ) { return 1/2.25*d; }

static void TestUserDefinedLiterals()
{
  long double total_weight = 1.2_g + 1.2_kg;
  std::cout << " Total weight = " << total_weight << std::endl;
}

namespace MemoryPool
{
  class Allocator
  {
  };
}

static void TestByte()
{
  constexpr size_t N = 1024;
  std::byte *p_plain = new std::byte[N];
  std::unique_ptr< std::byte[] > MEM_ARRAY{new std::byte[N]};
  for( size_t i=0; i < N; ++i ) MEM_ARRAY[i] = std::byte{0};
  delete[] p_plain;
}

class DoesFunctionConsumeSpace
{
  using FN = std::function<double(double)>;
private:
  static FN compute;
  double m_input;
public:
  //DoesFunctionConsumeSpace( double input ): compute{[](double x){return x*x;}},m_input{ input }  {}
  DoesFunctionConsumeSpace( double input ): m_input{ input }  {}
  double Run() const { return compute( m_input ); }
};

DoesFunctionConsumeSpace::FN DoesFunctionConsumeSpace::compute = [](double x) { return x*x;};

template <typename T>
class SourceData
{
public:
  static constexpr bool IS_SOURCE = true;
  using DATA_TYPE = T;
  SourceData(T input): m_input{ input } {}
  T get() const { return m_input; }
protected:
  T m_input;
};

template <typename T>
class SinkData
{
public:
  static constexpr bool IS_SINK = true;
  using DATA_TYPE = T;  
  SinkData() = default;
  void put( T data ) const { std::cout << "Recvd. " << data << std::endl; }
  void end() const { std::cout << "End." << std::endl; }
};


template <typename T, typename U, typename V, template<typename=U> class SourceData, template<typename=V> class SinkData>
class ComputeUnit
{
  static_assert( std::is_same<typename SourceData<U>::DATA_TYPE, typename SinkData<V>::DATA_TYPE>::value, "Incompatible source/sink data types." );
  using Source = SourceData<U>;
  using Sink   = SinkData<V>;
  static_assert( Source::IS_SOURCE, "Class is not a source of data." );
  static_assert( Sink::IS_SINK, "Class is not a sink of data." );
public:
  ComputeUnit( const Source& source, const Sink& sink ): m_source{ source }, m_sink{ sink } {}
  void Run() { m_sink.put( m_source.get() ); m_sink.end(); }
protected:
  const Source& m_source;
  const Sink&   m_sink;
};

template <typename U> using S = SourceData<U>;
template <typename V> using D = SinkData<V>;
template <typename T> using COMPUTE = ComputeUnit<T,T,T,S,D>;
static void TestComputeUnit()
{
  //using C = COMPUTE<int>;
  using C = ComputeUnit<int,int,int,S,D>;
  S<int> my_source{10};
  D<int> my_dest;
  C my_compute(my_source, my_dest);
  my_compute.Run();
}



int main()
{
  TestUserDefinedLiterals();
  TestByte();
  TestFEM();
  CommonDataTypes::CheckCommonDataTypes();
  TestTypeID();
  DataStore::IntegerStore IDS;
  DataStore::IntegerStore IDS2 = std::move( IDS );
  IDS = std::move( IDS2 );
  decltype( IDS ) IDS3;
  DataStore::FloatStore   FDS;
  auto X = CheckConformance( IDS3 );
  std::cout << " X = " << X << " of type " << GetTypeString< decltype( X )>() << std::endl;
  std::cout << " X = " << X << " of type " << typeid( X ).name() << std::endl;
  CodeLogic::RelatedConcepts::InitializeDataStore( IDS );
  CodeLogic::RelatedConcepts::InitializeDataStore( FDS );
  std::cout << "sizeof(DoesFunctionConsumeSpace) = " << sizeof(DoesFunctionConsumeSpace) << std::endl;
  DoesFunctionConsumeSpace XDF{10};
  std::cout << "Run := " << XDF.Run() << std::endl;
  TestComputeUnit();
  return( EXIT_SUCCESS );
}
