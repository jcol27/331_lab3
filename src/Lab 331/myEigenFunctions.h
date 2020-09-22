/*****************************************************************************************************
	myEigenFunctions.h

	This file is where you'll put the header information about the functions you write.  This is just 
	the argument list (as it is in the first line of your function declarations in myEigenFunctions.cpp)
	followed by a semi-colon.  Two example functions have been done for you.

====================================================================================================*/

#include <math.h>
#include <string>
#include <vector>
#include <sstream> //istringstream
#include <iostream> // cout
#include <fstream> // ifstream

using namespace std;

double DotProduct(double **A, double **B, int n, int m);
double DotProduct(double *A, double *B, int n);

// Note that these two functions have the same name, but different argument lists.  This is called 
// "function overloading" and is allowed by C++ compilers.

double* DotProduct(double **A, double *v, int n);

struct Eigenpair {
  double value; //Eigenvalue
  double *vector; //Eigenvector
  int length; //Length of eigenvector
  void normalize() {
    // Set eigenvalue to norm of vector and normalize eigenvector
    value = sqrt(DotProduct(vector, vector, length));
    for (int i=0; i<length; i++)
      vector[i]/= value;
  }; //

  void print() {
    std::cout << value << ": \t";
    for (int i=0; i<length; i++)
      std::cout << vector[i] << "\t";
    std::cout << "\n";
  }
  
  // Constructor
  // Attribute value is set to 0.0 and attribute vector to an array of doubles with lenght n
Eigenpair(const int n) : value(0.0), length(n), vector(new double[n]) {} //Constructor

};

double** ReadData(string inputFileName, int n, int m);

double** CenterMatrix(double **A, int n, int m);

double** CovarianceMatrix(double **A, int n, int m);

Eigenpair power_method(double **A, double *v, int n, double tol);

void deflate(double **A, Eigenpair eigenpair);

void print_matrix(double **A, int n, int m);
