#include "headers.hpp"


double BLAS1(Vector a, Vector b)
{
	double s = 0;
    
  for(int i = 0; i < a.rows(); i++)
  		s += a(i)*b(i);

  return s;
}

Vector BLAS2(Matrix A, Vector x)
{
	Vector s(A.rows());
  	Vector a(x.rows());

  	for (int i = 0; i < A.rows(); ++i)
  	{
  		for (int j = 0; j < x.rows(); ++j)
  			a(j) = A(i,j);
  		s(i) = BLAS1(a, x);
  	}

  return s;
}

Vector BLAS3(Matrix A, Vector x)
{

	Vector s = BLAS2(A, x);
	
  s = BLAS2(A, s) + s + x;

  s.conservativeResize(x.rows());

  return s;
}