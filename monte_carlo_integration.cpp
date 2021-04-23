////////////////////////////////////////////////////////////////////////////////
// File   : monte_carlo_integration.cpp
// Author : Sandeep Koranne (C) 2020. All rights reserved.
// Purpose: Use the C++-11 random number functionality to integrate multidim 
//        : functions using sample-reject Monte Carlo algorithm
//
// (1) Compiled on Visual Studio also
////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>
#include <random>
#include <cstdlib> // for atoi
#include <memory>
#include <omp.h>

namespace MonteCarloIntegration {
  using UNIVARIATE_FUNCTION = std::function<double(double)>;
  class MCI
  {
  public:
    MCI( double A, double B, int N, const UNIVARIATE_FUNCTION& F ): m_A{A}, m_B{B}, m_N{N}, m_F{F} { SetupRNG(); }
    double Integral() const;
  private:
    void SetupRNG();
    double getRandomNumberX() const { return (*(URDP_X))(*(DREP)); }
    double getRandomNumberY() const { return (*(URDP_Y))(*(DREP)); }
  protected:
    double m_A, m_B;
    int m_N;
    const UNIVARIATE_FUNCTION& m_F;
    std::unique_ptr<std::mt19937> DREP;
    std::unique_ptr<std::uniform_real_distribution<double>> URDP_X, URDP_Y;
    double MAXIMUM_VALUE_OF_FUNCTION, MINIMUM_VALUE_OF_FUNCTION;
  };
};

void MonteCarloIntegration::MCI::SetupRNG()
{
  // Assumption 1: function is monotonic, so max is either f(A) or f(B)
  MAXIMUM_VALUE_OF_FUNCTION = std::max( m_F(m_A), m_F(m_B) );
  MINIMUM_VALUE_OF_FUNCTION = std::min( m_F(m_A), m_F(m_B) );
  //std::cout << "Function in : " << MINIMUM_VALUE_OF_FUNCTION << "\t" << MAXIMUM_VALUE_OF_FUNCTION << std::endl;
  std::random_device rd;
  std::seed_seq seed{ rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd() };
  DREP.reset(new std::mt19937( seed ));
  URDP_X.reset(new std::uniform_real_distribution<double>(m_A,m_B));
  URDP_Y.reset(new std::uniform_real_distribution<double>(0,MAXIMUM_VALUE_OF_FUNCTION));

  if(false){
    std::ofstream f("a.dat");
    for( int i=0; i < 100; ++i ) { 
      double x = getRandomNumberX();
      f << x << "\t" << m_F(x) << "\t" << getRandomNumberY() << "\n";
    }
  }
  #if 0
  // Generate a normal distribution around that mean
  std::seed_seq seed2{r(), r(), r(), r(), r(), r(), r(), r()}; 
  std::mt19937 e2(seed2);
  std::normal_distribution<> normal_dist(mean, 2);
  #endif
}

////////////////////////////////////////////////////////////////////////////////
// Originally the plan was for multiple random evals in this loop for the
// trials. But there is lot of correlation between the RNG across multiple
// threads. It is thus better if multiple polynomials are computed in
// parallel. The polynomials have the same domain, so actually the X rng
// could in theory be shared, redcuing the work, but in future we want to
// explore multiple domains.
//
////////////////////////////////////////////////////////////////////////////////
double MonteCarloIntegration::MCI::Integral() const
{
  double retval=0;
  for( int i=0; i < m_N; ++i ) {
    double x = getRandomNumberX();
    double y = getRandomNumberY();
    //std::cout << "Trial " << i << "\t X = " << x << "\t f(x) = " << m_F(x) << "\t Y = " << y << "\n";
    if( y < m_F(x) ) retval++;
  }
  return (m_B-m_A)*((MAXIMUM_VALUE_OF_FUNCTION*retval)/(double)m_N);
}

using namespace MonteCarloIntegration;

template <typename T>
class Polynomial: std::vector<T>
{
  using V = std::vector<T>;
public:
  explicit Polynomial(size_t N): V(N) {}
  void RandomCoefficients();    
  Polynomial( const std::initializer_list<T>& input ): V{input} {}
  template <typename U>
  friend std::ostream& operator<<( std::ostream& os, const Polynomial<U>& P);
  double Integral( double A, double B ) const;
  UNIVARIATE_FUNCTION getHorner() const {
    auto constructed_lambda = [=](double x) -> double {
      double b = (*this)[this->size()-1];
      for( int i=static_cast<int>(this->size())-2; i >= 0; --i ) {
	b = (*this)[i] + b*x;
      }
      return b;
    };
    return constructed_lambda;
  }
  UNIVARIATE_FUNCTION getLambda() const {
    auto raise_to_power = [](double x, int n)->double {
      if( n==0 ) return 1.0;
      double retval = 1.0;
      for( int i=0; i < n; ++i ) retval *= x;
      return retval;
    };
    auto constructed_lambda= [=](double x) -> double {
      double retval = 0;
      int counter = 0;
      for( auto coeff : (*this) ) { retval += raise_to_power(x,counter++)*coeff; }
      return retval;
    };
    return constructed_lambda;
  }
};

template <typename T>
double Polynomial<T>::Integral( double A, double B ) const
{
  Polynomial<T> temp( this->size()+1 );
  for( size_t i=1; i < this->size(); ++i ) {
    temp[i+1] = (*this)[i]/(double)(i+1);
  }
  temp[1] = (*this)[0];
  auto L = temp.getHorner();
  return ( L(B)-L(A) );
}

template <typename T>
std::ostream& operator<<( std::ostream& os, const Polynomial<T>& P )
{
  os << "P := [ "; 
  for( auto coeff:P ) { os << coeff << " "; }
  return os << "] ";
}

template <typename T>
void Polynomial<T>::RandomCoefficients()
{
  std::random_device rd;  // produces a seed
  std::mt19937 gen(rd()); 
  //std::uniform_int_distribution<> distribution(0,100);
  std::uniform_real_distribution<>  distribution(0,1.0);
  for( auto& coeff : (*this) ) coeff = distribution(gen);
}

static void TestMonteCarlo( int M, int N )
{
  std::cout << "Running " << __FUNCTION__ << " with N = " << N << std::endl;
  using PX = Polynomial<double>;
  double max_rel_error = 0;
  std::cout << std::setw(8) << "POL INT" << "\t" << std::setw(8) << "MC INT" << "\t" 
	    << std::setw(8) << "ERROR" << "\t" << std::setw(8) << "REL %ERROR" << "\t\t" << "POLYNOMIAL" << std::endl;
  std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
  #pragma omp parallel for
  for( int i=0; i < M; ++i ) {
    PX px(10);
    px.RandomCoefficients();
    double A = 1.15;
    double B = 2.23;
    auto f{px.getHorner()}; // Visual Studio needs this outside
    MCI m(A,B,N,f);
    double mc_int = m.Integral();
    double pc_int = px.Integral(A,B);
    double error  = std::abs(pc_int-mc_int);
    double rel    = 100.0*error/pc_int;
    max_rel_error = std::max( max_rel_error, rel );
    #pragma omp critical
    std::cout << std::setw(8) << pc_int << "\t" << std::setw(8) << mc_int << "\t" 
	      << std::setw(8) << error << "\t" << std::setw(8) << rel << "\t" << px << std::endl;
  }
  std::cout << "------------------------------------------------------------------------------------------------" << std::endl;
  std::cout << "MAX REL ERROR := " << max_rel_error << std::endl;
}

static void Usage( const std::string& programName )
{
  std::cerr << programName << " M (number of polynomials) N (number of trials) T(number threads)" << std::endl;
}
int main(int argc, char* argv[])
{
  if( argc != 4 ) { Usage( argv[0] ); return (-1); }
  const int M = atoi( argv[1] );
  const int N = atoi( argv[2] );
  const int T = atoi( argv[3] );
  omp_set_num_threads( T );
  std::cout << "Running MC Integration with " << omp_get_max_threads() << " OpenMP threads.\n";
  TestMonteCarlo( M, N );
  return ( EXIT_SUCCESS );
}
#if 0
Running TestMonteCarlo with N = 1000000
 POL INT	  MC INT	   ERROR	REL ERROR	POLYNOMIAL
------------------------------------------------------------------------------------------------
 663.734	 663.186	0.548211	0.082595	P := [ 8 7 9 6 5 8 3 10 2 0 ] 
 1191.36	 1188.79	 2.56752	0.215512	P := [ 4 2 6 0 4 8 1 6 6 5 ] 
 974.142	 971.856	  2.2865	0.234719	P := [ 10 0 8 0 1 2 5 6 4 4 ] 
 630.767	 629.217	 1.54948	0.245651	P := [ 4 1 5 6 2 4 0 4 0 4 ] 
 1414.32	 1411.43	 2.88711	0.204134	P := [ 0 1 4 8 10 5 8 5 6 6 ] 
 1496.72	 1491.37	 5.34264	0.356958	P := [ 3 7 0 3 9 5 0 0 6 10 ] 
 1294.96	 1291.57	 3.38848	0.261667	P := [ 8 7 3 4 2 4 4 4 5 7 ] 
   893.3	 892.074	 1.22591	0.137234	P := [ 9 7 8 3 5 9 7 10 3 1 ] 
 1085.21	  1083.4	 1.81161	0.166937	P := [ 5 7 1 9 7 6 4 9 10 0 ] 
 1165.39	 1163.24	 2.15219	0.184676	P := [ 0 5 10 6 0 10 7 8 2 5 ] 
------------------------------------------------------------------------------------------------
MAX REL ERROR := 0.356958
#endif
