#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include <random>
#include <functional>

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
  std::function<double(double)> getHorner() const {
    auto constructed_lambda = [=](double x) -> double {
      double b = (*this)[this->size()-1];
      for( int i=static_cast<int>(this->size())-2; i >= 0; --i ) {
	b = (*this)[i] + b*x;
      }
      return b;
    };
    return constructed_lambda;
  }
  std::function<double(double)> getLambda() const {
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

#endif // POLYNOMIAL_H
