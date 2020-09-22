/*****************************************************************************************************

 This file is where you'll put the source information for the functions you write.  Make sure it's 
 included in your project (shift-alt-A) and that the functions you add are declared in the header
 file, myEigenFunctions.h.  You can add any extra comments or checks that could be relevant to these
 functions as you need to.

====================================================================================================*/

#include <string>

#include "myEigenFunctions.h"
#include <dos.h>

using namespace std;

double DotProduct(double **A, double **B, int n, int m)
{
	//
	//	This is a function that takes two matrices A and B of identical dimensions (n*m) and 
	//  calculates and returns their dot product.
	//
	double dot = 0.0;

	for(int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			dot += A[i][j]*B[i][j];
		}
	}
	return dot;
}

double DotProduct(double *A, double *B, int n) 
{
	//
	//	This is a function that takes two vectors A and B of identical length (n) and 
	//  calculates and returns their dot product.
	//

	double dot = 0.0;

	for(int i = 0; i < n; i++) {
		dot += A[i]*B[i];
	}
	return dot;
}

double* DotProduct(double **A, double *v, int n)
    {
          //
          //  This is a function that takes a nxn-matrix A and an n-dimensional vector v stores
          //  the product A.v at the original location of v
          //

    double *result = new double[n]; // pointer to result vector
  
    for(int i = 0; i < n; i++) {
      result[i] = 0.0; // initialize ith element of result v
      for(int j = 0; j <n; j++){
        result[i] += A[i][j] * v[j]; 
      }
    }

    return result;
  }



double** ReadData(string inputFileName, int n, int m)
    {
          //
          //  This is a function that returns a pointer to a 56x286 matrix, which contains the reflectance spectra.
          //  The first 29 rows contain spectra from Arabica samples, the following 27 rows contain spectra from Robusta samples.
          // The code for parsing the CSV file has been adapted from Bernardo Trindade
          // https://waterprogramming.wordpress.com/2017/08/20/reading-csv-files-in-c

      double **spectra = new double* [n];

     
      for(int sample = 0; sample < n; sample++)
	spectra[sample] = new double [m];

    vector<vector<double> > data;
    cout << inputFileName << endl;
    ifstream inputFile(inputFileName);
    int l = 0;
      
    while (inputFile) {
        string s;
        if (!getline(inputFile, s)) break;
	// cout << s << endl;
        if (l!=0) // ignore first line, which contains wavelengths
	  {
            istringstream ss(s);
            vector<double> record;

	    int column = 0;
            while (ss) {
	      // cout << "Row " << l << " " << ss << endl;
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
		  // cout << "Row " << l << " Column " << column << " line " << line << endl; 
                    spectra[l-1][column] = stod(line);
		    column++;
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
		}

        }
	l++;
    }

    return spectra;
  }

void print_matrix(double **A, int n, int m)
{
  //
  // This procedure prints an nxm matrix.
  //
  for (int i=0; i < n; i++)
    {
    for (int j=0; j < m; j++)
      std::cout << A[i][j] << '\t';
    std::cout << '\n';
    }
}

Eigenpair power_method(double** A, double* v, int n, double tol)
{
	/*
		Definition:
			- This function computes the largest eigenvalue and the corresponding eigenvector
		Inputs:
			- A, a pointer to an array of doubles, representing the initial matrix to solve
			- v, a pointer to a vector of doubles, representing the starting eigenvector
			- n, an int, representing the size of A (nxn) and v (n)
			- tol, a double, representing the tolerance to use in while loop
		Outputs:
			- eigenpair, an Eigenpair, with an eigenvalue, eigenvector, and length
	*/

	Eigenpair eigenpair(n);

	// Your code comes here

	// Initialize
	double lambda = 0.0;
	double temp = 1.0;
	Eigenpair t(n);

	// Initialize eigenvector and normalize
	eigenpair.vector = v;
	eigenpair.length = n;
	eigenpair.value = sqrt(DotProduct(eigenpair.vector, eigenpair.vector, eigenpair.length));
	eigenpair.normalize();
	

	// While loop
	while (fabs(lambda - temp) > tol) {
		// Compute t = A*v
		t.vector = DotProduct(A, eigenpair.vector, n);
		
		// Compute eigenvalue estimate
		temp = lambda;
		lambda = sqrt(DotProduct(t.vector, t.vector, t.length));
		t.value = lambda;

		// Compute normalized eigenvector estimate
		t.normalize();
		eigenpair.vector = t.vector;
		eigenpair.value = lambda;
	}

	cout << "eigenpair.value = " << eigenpair.value << "\n";

  return eigenpair;  
}

void deflate(double **A, Eigenpair eigenpair)
{
	/*
		Definition:
			- This procedure removes eigenpair.vector from the transformation matrix A in place
		Inputs:
			- A, a pointer to an array of doubles, representing the initial matrix to solve
			- eigenpair, an Eigenpair, with an eigenvalue, eigenvector, and length
		Outputs:
			- none
	*/

	// Your code comes here

	// Initializations
	int n = eigenpair.length;
	double *v = NULL;
	v = new double [n];
	double a = 0.0;

	for (int i = 0; i < n; i++) {
		v[i] = eigenpair.vector[i];
	}
	
	// Normalize eigenpair
	for (int i = 0; i < n; i++) {
		a = a + pow(v[i], 2);
	}
	a = pow(a, 0.5);
	for (int i = 0; i < n; i++) {
		v[i] = v[i] / a;
	}
	
	// Remove eigenpair vector from 
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A[i][j] = A[i][j] - eigenpair.value * v[i] * v[j];
		}
	}

	// Memory clean up
	delete[] v;
}

double** CenterMatrix(double **A, int n, int m)
{
	/*
		Definition: 
			- This is a function that takes a nxm-matrix A and subtracts the emperical mean 
			  of each column. The function returns the pointer to the centered matrix
		Inputs:
			- A, a pointer to an array of doubles, representing the initial matrix to solve
			- n and m, ints, representing the size of A, n = rows, m = columns
		Outputs:
			- result, an array of pointers to doubles, representing the centered matrix
    */
      
    double **result = new double* [n];
      
	for (int row = 0; row < n; row++) {
		result[row] = new double[m];
	}

    // Your code comes here

	//Initializations
	double sum = 0.0;

	// Loop over columns of A
	for (int j = 0; j < m; j++) {
		// Sum of column
		sum = 0.0;
		for (int i = 0; i < n; i++) {
			sum = sum + A[i][j];
		}
		// New value
		for (int i = 0; i < n; i++) {
			result[i][j] = A[i][j] - sum / n;
		}
	}
      
    return result;
}


double** CovarianceMatrix(double **A, int n, int m)
{
	/*
		Definition:
			- This is a function that takes a nxm-matrix A and subtracts the emperical mean of 
			  each column. The function returns the pointer to the centered matrix
		Inputs:
			- A, a pointer to an array of doubles, representing the initial matrix to solve
			- n and m, ints, representing the size of A, n = rows, m = columns
		Outputs:
			- cov, a pointer to an array of doubles, representing the covariance matrix
    */
      
	// Initializations
	double **cov = new double* [m];
	double **cA = CenterMatrix(A, n, m);
	double sum = 0.0;
	double *means = NULL;
	means = new double [m];
		
	for (int i = 0; i < m; i++) {
		cov[i] = new double[m];
	}

	// Your code comes here

	// Calculate the means of each column
	for (int j = 0; j < m; j++) {
		sum = 0.0;
		for (int i = 0; i < n; i++) {
			sum = sum + cA[i][j];
		}
		means[j] = sum / n;
	}

	// Loop over cov to find each value
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			// If Variance (forward diagonal)
			if (i == j) {
				cov[i][j] = 0.0;
				for (int u = 0; u < n; u++) {
					cov[i][j] = cov[i][j] + pow((cA[u][i]-means[i]),2);
				}
				cov[i][j] = cov[i][j] / (double(n) - 1);
			}
			// If Covariance (not forward diagonal)
			else {
				cov[i][j] = 0.0;
				for (int u = 0; u < n; u++) {
					cov[i][j] = cov[i][j] + cA[u][i] * cA[u][j];
				}
				cov[i][j] = cov[i][j] / (double(n) - 1);
			}
		}
	}
	

	
	delete[] cA;
	delete[] means;
	return cov;
}
