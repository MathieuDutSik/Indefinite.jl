/*****	This file includes some basic routines 
	for matrix and vector operations	*****/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"

/*****	multiplies 1xn-vector x with nxn-matrix
      A and puts result on y	*****/

void vecmatmul(long long *x, long long **A, long long n, long long *y)
{
  long long	i, j;
  long long *Ai, xi;
  for (j = 0; j < n; ++j)
    y[j] = 0;
  for (i = 0; i < n; ++i)
    {
      if ((xi = x[i]) != 0)
	{
	  Ai = A[i];
	  for (j = 0; j < n; ++j)
	    y[j] += xi * Ai[j];
	}
    }
}

void matmul(long long **A, long long **B, long long n, long long **C)
/*****	multiplies nxn-matrices A and B and puts result
      on C	*****/
{
  long long	i, j, k;
  long long *Ai, *Bk, *Ci, Aik;
  for (i = 0; i < n; ++i)
    {
      Ai = A[i];
      Ci = C[i];
      for (j = 0; j < n; ++j)
	Ci[j] = 0;
      for (k = 0; k < n; ++k)
	{
	  if ((Aik = Ai[k]) != 0)
	    {
	      Bk = B[k];
	      for (j = 0; j < n; ++j)
		Ci[j] += Aik * Bk[j];
	    }
	}
    }
}

long long scp(long long *x, long long **F, long long *y, long long n)
/*****	returns scalar product of 1xn-vectors x and y
      with respect to Gram-matrix F	*****/
{
  long long	i, j;
  long long sum, *Fi, xi;
  sum = 0;
  for (i = 0; i < n; ++i)
    {
      if ((xi = x[i]) != 0)
	{
	  Fi = F[i];
	  for (j = 0; j < n; ++j)
	    sum += xi * Fi[j] * y[j];
	}
    }
  return(sum);
}

long long sscp(long long *x, long long *y, long long n)
/*****	returns standard scalar product of 
      1xn-vectors x and y, heavily used	*****/
{
	long long	i;
	long long sum;
	sum = 0;
	for (i = 0; i < n; ++i)
		sum += *(x++) * *(y++);
	/*	sum += x[i] * y[i];	*/
	return(sum);
}

void psolve(long long **X, long long **A, long long **B, long long n, long long p)
/*****	computes X = B*A^-1 modulo the prime p
      for nxn-matrices A and B, where A is
      invertible,
      A and B are modified !!!	*****/
{
	long long	i, j, k;
	long long *Xi, *Bi, ainv, sum, tmp;
/* convert A to lower triangular matrix and change B accordingly */
	for (i = 0; i < n-1; ++i)
	{
		for (j = i; A[i][j] % p == 0; ++j);
		if (j == n)
		{
		  fprintf(stderr, "Error: matrix is singular modulo ");
		  single2Print(stderr, p);
		  fprintf(stderr, "\n");
		  exit (2);
		}
		if (j != i)
/* interchange columns i and j such that A[i][i] is != 0 */
		{
			for (k = i; k < n; ++k)
			{
				tmp = A[k][i];
				A[k][i] = A[k][j];
				A[k][j] = tmp;
			}
			for (k = 0; k < n; ++k)
			{
				tmp = B[k][i];
				B[k][i] = B[k][j];
				B[k][j] = tmp;
			}
		}
		pgauss(i, A, B, n, p);
	}
/* compute X recursively */
	for (i = 0; i < n; ++i)
	{
		Xi = X[i];
		Bi = B[i];
		for (j = n-1; j >= 0; --j)
		{
			sum = 0;
			for (k = n-1; k > j; --k)
				sum = (sum + Xi[k] * A[k][j]) % p;
/* ainv is the inverse of A[j][j] modulo p */
			for (ainv = 1; abs(A[j][j]*ainv-1) % p != 0; ++ainv);
			Xi[j] = (Bi[j] - sum) * ainv % p;
/* make sure that -p/2 < X[i][j] <= p/2 */
			if (2*Xi[j] > p)
				Xi[j] -= p;
			else if (2*Xi[j] <= -p)
				Xi[j] += p;
		}
	}
}

void pgauss(long long r, long long **A, long long **B, long long n, long long p)
/*****	clears row nr. r of A assuming that the 
      first r-1 rows of A are already cleared 
      and that A[r][r] != 0 modulo p, 
      A and B	are changed !!!	*****/
{
	long long	i, j;
	long long *Ar, f, ainv;

	Ar = A[r];
/* ainv is the inverse of A[r][r] modulo p */
	for (ainv = 1; abs(Ar[r]*ainv-1) % p != 0; ++ainv);
	for (j = r+1; j < n; ++j)
	{
		if (Ar[j] % p != 0)
		{
			f = Ar[j] * ainv % p;
			for (i = r+1; i < n; ++i)
				A[i][j] = (A[i][j] - f * A[i][r]) % p;
			for (i = 0; i < n; ++i)
				B[i][j] = (B[i][j] - f * B[i][r]) % p;
			Ar[j] = 0;
		}
	}
}

long long isprime(long long n)
/*****	checks, whether n is a prime	*****/
{
	long long	i;
	for (i = 2; i <= n/i; ++i)
	{
		if (n % i == 0)
			return 0;
	}
	return 1;
}
