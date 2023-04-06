/*****	This file contains routines for lattice reduction and the
	calculation of shortest vectors	*****/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"

long long lll(long long **F, long long **T, long long dim)
/*****	LLL-reduction of the positive semidefinite 
      Gram-matrix F with transformation matrix T	*****/
{
	long long	i, j, m, l, rank;
	double	**mu, *B;
	long long **gram, *ptmp, tmp;

/* the weights */
	if ((B = (double*)malloc(dim * sizeof(double))) == 0)
		exit(EXIT_FAILURE);
/* the model matrix */
	if ((mu = (double**)malloc(dim * sizeof(double*))) == 0)
		exit(EXIT_FAILURE);
/* F is copied on gram and remains unchanged */
	if ((gram = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit(EXIT_FAILURE);
	for (i = 0; i < dim; ++i)
	{
		if ((mu[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit(EXIT_FAILURE);
		if ((gram[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit(EXIT_FAILURE);
		for (j = 0; j < dim; ++j)
		{
			T[i][j] = 0;
			mu[i][j] = 0.0;
			gram[i][j] = F[i][j];
		}
		T[i][i] = 1;
		mu[i][i] = 1.0;
	}
/* the first basis-vector should not have norm 0 */
	for (i = 0; i < dim  &&  gram[i][i] == 0; ++i);
	if (i == dim)
		return(0);
	else if (i > 0)
	{
		ptmp = gram[i];
		gram[i] = gram[0];
		gram[0] = ptmp;
		for (j = 0; j < dim; ++j)
		{
			tmp = gram[j][i];
			gram[j][i] = gram[j][0];
			gram[j][0] = tmp;
		}
		T[0][0] = 0;
		T[0][i] = 1;
		T[i][i] = 0;
		T[i][0] = 1;
	}
	rank = dim;
	initialize(0, mu, B, gram);
	if (dim > 1)
		initialize(1, mu, B, gram);
/* recursively go through the matrix */
	m = 1;
	l = 0;
	while (m < rank)
	  m = red(&m, &l, gram, T, mu, B, dim, &rank);
	free(B);
	for (i = 0; i < dim; ++i)
	{
		free(mu[i]);
		free(gram[i]);
	}
	free(mu);
	free(gram);
	return(rank);
}

void initialize(long long i, double **mu, double *B, long long **gram)
/*****	initialization of the model	*****/
{
  long long	j, k;
  double	muh, *mui, *muj;
  
  B[i] = (double)gram[i][i];
  mui = mu[i];
  for (j = 0; j < i; ++j)
    {
      muh = (double)gram[i][j];
      muj = mu[j];
      for (k = 0; k < j; ++k)
	muh -= B[k] * mui[k] * muj[k];
      muh /= B[j];
      B[i] -= B[j] * muh * muh;
      mui[j] = muh;
    }
}

long long red(long long *m, long long *l, long long **gram, long long **T, double **mu, double *B, long long dim, long long *rank)
/*****	pair-reduction with vectors 
      m and l, changes Gram-matrix, 
      transformation matrix and 
      model appropriately	*****/
{
	double	r;
	long long	j, jump;
	long long ir;
	long long *gramm, *graml;
	long long *Tm, *Tl, *ptmp;
	jump = 0;
	r = mu[*m][*l];
	if (r > 0.5  ||  r < -0.5)
	{
		ir = MyRound(r);
		Tm = T[*m];
		Tl = T[*l];
		gramm = gram[*m];
		graml = gram[*l];
		for (j = 0; j < dim; ++j)
		{
			Tm[j] -= ir * Tl[j];
			gramm[j] -= ir * graml[j];	
		}
		for (j = 0; j < dim; ++j)
			gram[j][*m] -= ir * gram[j][*l];	
		initialize(*m, mu, B, gram); 
	}
	gramm = gram[*m];
	for (j = 0; j < *rank  &&  gramm[j] == 0; ++j);
	if (j == *rank)
	{
		ptmp = T[*m];
		T[*m] = T[*rank-1];
		T[*rank-1] = ptmp;
		for (j = 0; j < dim; ++j)
		{
			gram[*m][j] = gram[*rank-1][j];
			gram[*rank-1][j] = 0;
		}
		for (j = 0; j < dim; ++j)
		{
			gram[j][*m] = gram[j][*rank-1];
			gram[j][*rank-1] = 0;
		}
		*rank -= 1;
		if (*m < *rank)
		{
			initialize(*m, mu, B, gram); 
			*l = *m-1;
			jump = 1;
		}
		else
			jump = 1;
	}
	if (*l == *m-1  &&  jump == 0)
		check(B, m, l, gram, T, mu, dim, rank);
	else if (jump == 0)
		decrease(m, l, gram, mu, B, rank);
	return(*m);
}

void check(double *B, long long *m, long long *l, long long **gram, long long **T, double **mu, long long dim, long long *rank)
/*****	check the 
      LLL-condition	*****/
{
  if (B[*m] < (LLL_CONST - mu[*m][*m-1]*mu[*m][*m-1]) * B[*m-1])
    interchange(m, l, gram, T, mu, B, dim);
  else
    {
      *l = *m-1;
      decrease(m, l, gram, mu, B, rank);	
    }
}

void decrease(long long *m, long long *l, long long **gram, double **mu, double *B, long long *rank)
/*****	go back in the recursion 
      and change the model	*****/
{
  *l -= 1;
  if (*l < 0)
    {
      *m += 1;
      if (*m < *rank)
	{
	  *l = *m-1;
	  initialize(*m, mu, B, gram);
	}
    }
}

void interchange(long long *m, long long *l, long long **gram, long long **T, double **mu, double *B, long long dim)
/*****	interchange the vectors 
      nr. m and m-1 and 
      change the model 
      appropriately	*****/
{
	long long	i;
	long long tmp, *ptmp;
	ptmp = T[*m];
	T[*m] = T[*m-1];
	T[*m-1] = ptmp;
	ptmp = gram[*m];
	gram[*m] = gram[*m-1];
	gram[*m-1] = ptmp;
	for (i = 0; i < dim; ++i)
	{
		tmp = gram[i][*m];
		gram[i][*m] = gram[i][*m-1];
		gram[i][*m-1] = tmp;
	}
	initialize(*m-1, mu, B, gram);
	if (*m > 1)
	{
		*m -= 1;
		*l = *m-1;
	}
	else
	{
		initialize(*m, mu, B, gram);
		*l = *m-1;
	}
}

long long MyRound(double x)
/*****	rounds x to an integer	*****/
{
	if (x >= 0)
	  {
	    return((long long)(2.0*x + 1.0) / 2);
	  }
	else
	  {
	    return(-((long long)(-2.0*x + 1.0) / 2));
	  }
}

void shortvec(long long **F, long long max, veclist *V)
/*****	calculates the vectors of norm <= max 
      with respect to the Gram-matrix F, 
      writes the vectors on V->v and the 
      length on V->v[i][dim],
      the algorithm uses the model of the 
      LLL-algorithm	*****/
{
	long long     dim, i, j, m, l, rank, con;
	double	**mu, *B;
	long long **gram, **T, *vec;

	dim = V->dim;
	if ((B = (double*)malloc(dim * sizeof(double))) == 0)
		exit(EXIT_FAILURE);
	if ((mu = (double**)malloc(dim * sizeof(double*))) == 0)
		exit(EXIT_FAILURE);
/* F is copied on gram and remains unchanged */
	if ((gram = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit(EXIT_FAILURE);
	if ((T = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit(EXIT_FAILURE);
	if ((vec = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit(EXIT_FAILURE);
	for (i = 0; i < dim; ++i)
	{
		vec[i] = 0;
		if ((mu[i] = (double*)malloc(dim * sizeof(double))) == 0)
			exit(EXIT_FAILURE);
		if ((gram[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit(EXIT_FAILURE);
		if ((T[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit(EXIT_FAILURE);
		for (j = 0; j < dim; ++j)
		{
			T[i][j] = 0;
			mu[i][j] = 0.0;
			gram[i][j] = F[i][j];
		}
		T[i][i] = 1;
		mu[i][i] = 1.0;
	}
/* no basis-vector should have norm 0 */
	for (i = 0; i < dim  &&  gram[i][i] > 0; ++i);
	if (i != dim)
	{
		fprintf(stderr, "Error: matrix is not positive definite\n");
		exit (3);
	}
	rank = dim;
	initialize(0, mu, B, gram);
	if (dim > 1)
		initialize(1, mu, B, gram);
/* do the LLL-reduction */
	m = 1;
	l = 0;
	while (m < rank  &&  rank == dim)
		m = red(&m, &l, gram, T, mu, B, dim, &rank);
	if (rank < dim)
	{
		fprintf(stderr, "Error: matrix is not positive definite\n");
		exit (3);
	}
	if ((V->v = (long long**)malloc(sizeof(long long*))) == 0)
		exit(EXIT_FAILURE);
	if ((V->v[0] = (long long*)malloc(V->len * sizeof(long long))) == 0)
		exit(EXIT_FAILURE);
	V->n = 0;
	con = 0;
/* recursively find the vectors up to length max */
	shrt(dim-1, 0.0, max, V, vec, mu, T, B, &con);
	free(B);
	free(vec);
	for (i = 0; i < dim; ++i)
	{
		free(mu[i]);
		free(T[i]);
		free(gram[i]);
	}
	free(mu);
	free(T);
	free(gram);
}

void shrt(long long k, double norm, long long max, veclist *V, long long *vec, double **mu, long long **T, double *B, long long *con)
/*****	recursion for finding the vectors up 
      to length max, where the components
      vec[k+1],...,vec[dim-1] are already
      chosen and contribute norm to the 
      length of the vector	*****/
{
  double	x, Bk;
  long long	i, j, dim;
  long long *Vvn, Vvni;
  dim = V->dim;
  if (k == -1)
    {
      for (i = 0; i < dim  &&  vec[i] == 0; ++i);
      if (i == dim)
	/* the 0-vector indicates the end of the recursion, since the short vectors
	   are only enumerated up to sign */
	*con = 1;
      else
	/* a short vector is found */
	{
	  ++V->n;
	  if ((V->v = (long long**)realloc(V->v, (V->n+1) * sizeof(long long*))) == 0)
	    exit(EXIT_FAILURE);
	  if ((V->v[V->n] = (long long*)malloc(V->len * sizeof(long long))) == 0)
	    exit(EXIT_FAILURE);
	  Vvn = V->v[V->n];
	  for (i = 0; i < dim; ++i)
	    {
	      Vvni = 0;
	      for (j = 0; j < dim; ++j)
		Vvni += vec[j] * T[j][i];
	      Vvn[i] = Vvni;
	    }
	  Vvn[dim] = MyRound(norm);
	}
    }
  else
    {
      x = 0.0;
      for (j = k+1; j < dim; ++j)
			x += vec[j] * mu[j][k];
      i = -MyRound(x);	
      Bk = B[k];
      if (Bk*(x+i)*(x+i) + norm < max + EPS)
	{
	  while ((Bk*(x+i)*(x+i) + norm <= max + EPS) || 
		 (x+i <= 0))
	    ++i;
	  --i;
	  while ((Bk*(x+i)*(x+i) + norm < max + EPS) && *con == 0)
	    {
	      vec[k] = i;
	      /* go deeper into the recursion */
	      shrt(k-1, Bk*(x+i)*(x+i)+norm, max, V, vec, mu, T, B, con);
	      --i;
	    }
	}
    }
}

void change(long long **v, long long nv, long long **T, long long n)
/*****	changes v according to the 
      transformation matrix T,
      i.e. v := T*v	*****/
{
  long long	i, j, k;
  long long fac, **w, *wi, *Ti, *vi, *vj;
  if ((w = (long long**)malloc(nv * sizeof(long long*))) == 0)
    exit(EXIT_FAILURE);
  for (i = 0; i < nv; ++i)
    {
      if ((w[i] = (long long*)malloc(n * sizeof(long long))) == 0)
	exit(EXIT_FAILURE);
      wi = w[i];
      for (j = 0; j < n; ++j)
	wi[j] = 0;
    }
  for (i = 0; i < nv; ++i)
    {
      wi = w[i];
      Ti = T[i];
      for (j = 0; j < nv; ++j)
	{
	  if ((fac = Ti[j]) != 0)
	    {
	      vj = v[j];
	      for (k = 0; k < n; ++k)
		wi[k] += fac * vj[k];
	    }
	}
    }
  for (i = 0; i < nv; ++i)
    {
      vi = v[i];
      wi = w[i];
      for (j = 0; j < n; ++j)
	vi[j] = wi[j];
      free(wi);
    }
  free(w);
}
