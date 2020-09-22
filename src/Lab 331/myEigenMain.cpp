
/******************************************************************************************************

	This is a template main file for the ENGSCI331 Eigenvectors module.  It demonstrates some new C++ 
	syntax and functions, as described in the accompanying document ENGSCI331_Eigenstuff.pdf.  

	*** There are some examples of "bad" programming in here (bits missing etc) that you will need to 
	find and fix, though this file should compile without errors straight away. ***

	You should use this file to get you started on the assignment.  You're welcome to change whatever
	you'd like to as you go, this is only a starting point.
	
======================================================================================================*/

#define _CRT_SECURE_NO_DEPRECATE
#define COUNT 56 // Number of principle components to calculate

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// This is the header file for your functions.  Usual programming practise would be to use a *.cpp 
// file that has the same name (ie: myEigenFunctions.cpp) and include it as normal in your project.
// Inside this file you'll see some ideas for functions that you could use during this project.  I 
// suggest you plan out your code first to see what kind of functions you'll use repeatedly and then
// write them.  
#include "myEigenFunctions.h"

using namespace std;

#define PI 3.14159265358979323846
#define TOL 1E-6

bool AE(double a, double b) {
	return (abs(a - b) < TOL); // Compare two values with respect to tolerance
}

int main(void) 
{
	// ----------------------------------------------------------------------------------------------
	//
	// PART 1: Initialisation
	//
	// ----------------------------------------------------------------------------------------------
	
	// Defining local variables to be used:

  int n = 2, nc=3;
  double *v, *result = NULL;

  double **A = NULL;
  double **C = NULL;
  // Allocating memory for the 1D arrays - these are the number of masses, n, long:
  v = new double [n];

  // struct comprising eigenvalue, eigenvector pair
  Eigenpair eigenpair(2), expected_eigenpair(2), deflate_eigenpair(3);
    
  // Allocating memory for the 2D arrays - these have dimensions of n*n:
  A = new double* [n];
  C = new double* [nc];
  for(int i = 0; i < n; i++) 
    A[i] = new double [n];
  for(int i = 0; i < nc; i++) 
    C[i] = new double [nc];

  if (eigenpair.value == 0.0)
    cout << "OK: Initialization of eigenpair with eigenvalue 0.\n";
  else
    cout << "ERROR: Initialization of eigenpair with eigenvalue 0.\n";

  eigenpair.vector[0] = 2.0;
  eigenpair.vector[1] = 0.0;
  eigenpair.normalize();
  if (eigenpair.value == 2.0 &&
      eigenpair.vector[0] == 1.0 &&
      eigenpair.vector[1] == 0.0)
    cout << "OK: Normation of eigenpair.\n";
  else
    cout << "ERROR: Normation of eigenpair. ";

  A[0][0] = 4.0;
  A[0][1] = 6.0;
  A[1][0] = 1.0;
  A[1][1] = 3.0;

  v[0] = 3.0;
  v[1] = 1.0;

  result = DotProduct(A, v, 2);
  if (result[0] == 18.0 && result[1] == 6.0)
    cout << "OK: DotProduct of Matrix and Vector\n";
  else
    cout << "ERROR: DotProduct of Matrix and Vector\n";

  v[0] = 1.0;
  v[1] = 1.0;
  eigenpair = power_method(A, v, 2, TOL);

  if (abs(eigenpair.value-6.0)/6.0 < TOL)
    cout << "OK: Power Method computing largest eigenvalue\n";
  else
    cout << "ERROR: Power Method computing largest eigenvalue\n";

  if (abs(eigenpair.vector[0]-0.94868417)/0.94868417 < sqrt(TOL) &&
      abs(eigenpair.vector[1]-0.31622515)/0.31622515 < sqrt(TOL))
    cout << "OK: Power Method computing eigenvector of largest eigenvalue\n";
  else {
    cout << "ERROR: Power Method computing eigenvector of largest eigenvalue\n";
    cout << "(" << eigenpair.vector[0] << "," << eigenpair.vector[1]  << ")\n";
  }

  
  ifstream infile;

  infile.open("C.txt");
  for (int i=0; i < nc; i++)
    for (int j=0; j< nc; j++)
      infile >> C[i][j];
  cout << "Test matrix" << endl;
  print_matrix(C, 3, 3);

    double **cov = CovarianceMatrix(C, 3, 3);
  //  cout << "CovarianceMatrix" << endl;
  // print_matrix(cov, 3, 3);
  if (AE(cov[0][0], 1.0) && AE(cov[0][1], 0.5) && AE(cov[0][2], -1.0) &&
      AE(cov[1][0], 0.5) && AE(cov[1][1], 1.0) && AE(cov[1][2], -1.0) &&
      AE(cov[2][0], -1.0) && AE(cov[2][1],-1.0) && AE(cov[2][2], 4.0/3)
      )
    cout << "OK: covariance matrix computed correctly\n";
 else {
       cout << "ERROR: in computation of covariance matrix\n";
       print_matrix(C, deflate_eigenpair.length, deflate_eigenpair.length);
 }


  deflate_eigenpair.value = 3.0;
  deflate_eigenpair.vector[0] = 1.0;
  deflate_eigenpair.vector[1] = 1.0;
  deflate_eigenpair.vector[2] = 0.0;

  deflate(C, deflate_eigenpair);
  if (AE(C[0][0],-0.5) && AE(C[0][1], 0.5) && AE(C[0][2], 0.0) &&
      AE(C[1][0], 0.5) && AE(C[1][1], -0.5) && AE(C[1][2], 0.0) &&
      AE(C[2][0], 0.0) && AE(C[2][1],  0.0) && AE(C[2][2], 2.0)
      )
    cout << "OK: deflate method computes deflated matrix correctly\n";
 else {
       cout << "ERROR: deflate method does not compute deflated matrix correctly\n";
       print_matrix(C, deflate_eigenpair.length, deflate_eigenpair.length);
 }
  
  // Test reading of CSV file
  double **spectra;

  const int N = 56;   // rows
  const int M = 286;  // columns

  spectra = ReadData("DS19hH2_dk0_FTIR_Spectra_instant_coffee.csv", N, M);
  if (AE(spectra[0][0], 21.227620) && AE(spectra[N-1][M-1], 1.574679))
    cout << "OK: Reading of CSV file" << endl;
  else
    cout << "ERROR: reading of CSV file" << spectra[0][0] << " " << spectra[N-1][M-1] << endl;

  //print_matrix(spectra, N, M);
  cov = CovarianceMatrix(spectra, N, M);

  // Initializations
  double* x = NULL;
  x = new double[M]; // Start vector
  std::vector<Eigenpair> pairs(COUNT, Eigenpair(M));

  for (int i = 0; i < M; i++) {
	  x[i] = 1;
  }
  
  // While loop to solve
  for (int i = 0; i < COUNT; i++) {
	  pairs[i] = power_method(cov, x, M, TOL);
	  deflate(cov, pairs[i]);
  }

  // Order based on decreasing eigenvalue size
  std::vector<Eigenpair> sorted(0, Eigenpair(M));
  while (size(pairs) > 0) {
	  Eigenpair max_eig_pair(M);
	  double max_eig = 0.0;
	  int index = 0;
	  // Find the largest eigenvalue and its index
	  for (int j = 0; j < size(pairs); j++) {
		  if (pairs[j].value > max_eig) {
			  max_eig = pairs[j].value;
			  index = j;
		  }
	  }

	  // Add that pair to sorted, then remove from pairs
	  sorted.push_back(pairs[index]);
	  pairs.erase(pairs.begin()+index);
  }

  // Write components
  std::ofstream file1;
  file1.open("Sorted Eigenvectors.csv", ios::out | ios::trunc);
  for (int i = 0; i < COUNT; i++) {
	  for (int j = 0; j < M - 1; j++) {
		  file1 << sorted[i].vector[j] << ",";
	  }
	  file1 << sorted[i].vector[M - 1] << "\n";
  }
  file1.close();
  
  std::ofstream file2;
  file2.open("Eigenvalues.csv", ios::out | ios::trunc);
  for (int i = 0; i < size(g); i++) {
	  file2 << i+1 << "," << sorted[i].value << "\n";
  }
  file2.close();


  // Memory clean up
  delete[] W;
  delete[] T;
  delete[] spectra;
  delete[] x;
  delete[] cov;
  delete[] A;
  delete[] C;
  delete[] v;



} // end main

