/*****	This file includes some routines for 
	the preprocessing of the algorithm	*****/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"

void checkvecs(veclist *V, invar F, veclist norm)
/*****	removes those vectors from V.v, which 
      have other norm combinations than the 
      base-vectors	*****/
{
	long long	i, j, k, dim;
	long long *Vvi, *normvec;
	dim = F.dim;
	if ((normvec = (long long*)malloc(norm.dim * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
	j = 0;
	for (i = 1; i <= V->n; ++i)
	{
		Vvi = V->v[i];
		for (k = 0; k < F.n; ++k)
		{
			Vvi[dim+k] = scp(Vvi, F.A[k], Vvi, dim);
			normvec[k] = Vvi[dim+k];
		}
		if (abs(numberof(normvec, norm)) > dim)
/* the vector V->V[i] has a different norm combination from those in the list
   norm.v and is therefore deleted */
		{
			free(Vvi);
			++j;
		}
		else
			V->v[i-j] = Vvi;
	}
	V->n -= j;
	free(normvec);
}

long long checkgen(long long **g, invar F)
/*****	checks, whether g fixes the invariant forms	*****/
{
  long long	i, j, k, fix;
  long long **FAk;
  fix = 1;
  for (k = 0; k < F.n  &&  fix == 1; ++k)
    {
      FAk = F.A[k];
      for (i = 0; i < F.dim  &&  fix == 1; ++i)
	{
	  for (j = 0; j <= i  &&  fix == 1; ++j)
	    {
	      if (scp(g[i], FAk, g[j], F.dim) != FAk[i][j])
		fix = 0;
	    }
	}
    }
  return(fix);
}

void fingerprint(fpstruct *fp, invar F, veclist	V)
/*****	a permutation per := fp.per for the
      order of the basis-vectors is chosen
      such that in every step the number of
      possible continuations is minimal,
      for j from per[i] to per[dim-1] the
      value f[i][j] in the fingerprint f is
      the number of vectors, which have the
      same scalar product with the 
      basis-vectors per[0]...per[i-1] as the 
      basis-vector j and the same length as 
      this vector with respect to all 
      invariant forms	*****/
{
	long long	i, j, k, dim, min, tmp, **f;
	long long *Vvj, ***Ftr;

	dim = F.dim;
	if ((Ftr = (long long***)malloc(F.n * sizeof(long long**))) == 0)
	  exit(EXIT_FAILURE);
	for (i = 0; i < F.n; ++i)
	{
		if ((Ftr[i] = (long long**)malloc(dim * sizeof(long long*))) == 0)
			exit (EXIT_FAILURE);
		for (j = 0; j < dim; ++j)
		{
			if ((Ftr[i][j] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
			for (k = 0; k < dim; ++k)
				Ftr[i][j][k] = F.A[i][k][j];
		}
	}
	if ((f = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i < dim; ++i)
	{
		if ((f[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
		for (j = 0; j < dim; ++j)
			f[i][j] = 0;
	}
	for (i = 0; i < dim; ++i)
		fp->per[i] = i;
/* the first row of the fingerprint has as entry nr. i the number of
   vectors, which have the same length as the i-th basis-vector with
   respect to every invariant form */
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		for (i = 0; i < dim; ++i)
		{
			for (k = 0; k < F.n  &&  Vvj[dim+k] == F.A[k][i][i]; ++k);
			if (k == F.n)
				f[0][i] += 2;
		}
	}
	for (i = 0; i < dim-1; ++i)
	{
/* a minimal entry != 0 in the i-th row is chosen */
		min = i;
		for (j = i+1; j < dim; ++j)
		{
			if (f[i][fp->per[j]] < f[i][fp->per[min]])
				min = j;
		}
		tmp = fp->per[i];
		fp->per[i] = fp->per[min];
		fp->per[min] = tmp;
/* the column below the minimal entry is set to 0 */
		for (j = i+1; j < dim; ++j)
			f[j][fp->per[i]] = 0;
/* compute the row i+1 of the fingerprint */
		for (j = i+1; j < dim; ++j)
		  f[i+1][fp->per[j]] = possible(F, Ftr, V, fp->per, i, fp->per[j]);
	}
/* only the diagonal of f will be needed later */
	for (i = 0; i < dim; ++i)
	{
		fp->diag[i] = f[i][fp->per[i]];
		free(f[i]);
	}
	free(f);
	for (i = 0; i < F.n; ++i)
	{
		for (j = 0; j < dim; ++j)
			free(Ftr[i][j]);
		free(Ftr[i]);
	}
	free(Ftr);
}

long long possible(invar F, long long ***Ftr, veclist V, long long *per, long long I, long long J)
/*****	returns the number of vectors, which
      have the same scalar products with the
      basis-vectors nr. per[0]...per[I] as the
      basis-vector nr. J and the same length
      as this vector with respect to all
      invariant forms	*****/
{
	long long	i, j, k, dim, count;
	long long *Vvj;
	dim = F.dim;
	count = 0;
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		i = I+1;
/* check the length of the vector */
		for (k = 0; k < F.n  &&  i > I  &&  Vvj[dim+k] == F.A[k][J][J]; ++k)
/* check the scalar products with the basis-vectors */
		/*	for (i = 0; i <= I  &&  F.v[k][j][per[i]] == F.A[k][per[i]][J]; ++i); */
			for (i = 0; i <= I  &&  sscp(Vvj, Ftr[k][per[i]], dim) == F.A[k][J][per[i]]; ++i);
		if (k == F.n  &&  i > I)
			++count;
/* the same for the negative vector */
		i = I+1;
		for (k = 0; k < F.n  &&  i > I  &&  Vvj[dim+k] == F.A[k][J][J]; ++k)
		/*	for (i = 0; i <= I  &&  F.v[k][j][per[i]] == -F.A[k][per[i]][J]; ++i); */
			for (i = 0; i <= I  &&  sscp(Vvj, Ftr[k][per[i]], dim) == -F.A[k][J][per[i]]; ++i);
		if (k == F.n  &&  i > I)
			++count;
	}
	return(count);
}

void scpvector(long long *scpvec, long long *w, long long *b, long long I, long long dep, invar	F)
/*****	calculates the scalar products
      of the vector w with the base
      vectors v[b[I]] down to
      v[b[I-dep+1]] with respect to
      all invariant forms and puts
      them on scpvec	*****/
{
  long long	i, j, dim, bi;
  
  dim = F.dim;
  for (i = I; i >= 0  &&  i > I-dep; --i)
    {
      if ((bi = b[i]) > 0)
	{
	  for (j = 0; j < F.n; ++j)
	    scpvec[j*dep + I-i] = sscp(w, F.v[j][bi], dim);
	}
      else
	{
	  for (j = 0; j < F.n; ++j)
	    scpvec[j*dep + I-i] = -sscp(w, F.v[j][-bi], dim);
	}
    }
}

void scpvecs(veclist *list, long long ***vec, long long I, long long *b, long long dep, veclist V, invar F)
/*****	computes the list of scalar product 
      combinations of the vectors in V.v with
      the basis-vectors in b	*****/
{
	long long	i, j, dim, len, nr, sign;
	long long *Vvj, *listvn, **listv, *tmp;
	long long *vecnr, *vecn, *scpvec;

	dim = F.dim;
	len = F.n * dep;
/* scpvec is the vector for the scalar product combination */
	if ((scpvec = (long long*)malloc(len * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* the first vector in the list is the 0-vector and is not counted */
	if ((list->v = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	list->dim = len;
	list->len = len;
	if ((list->v[0] = (long long*)malloc(len * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i < len; ++i)
		list->v[0][i] = 0;
	list->n = 0;
	if ((*vec = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	if (((*vec)[0] = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i < dim; ++i)
		(*vec)[0][i] = 0;
	for (j = 1; j <= V.n; ++j)
	{
		Vvj = V.v[j];
		for (i = 0; i < len; ++i)
			scpvec[i] = 0;
		scpvector(scpvec, Vvj, b, I, dep, F);
		for (i = 0; i < len  &&  scpvec[i] == 0; ++i);
		if (i == len)
/* if scpvec is the 0-vector, nr is set to 0, since numberof never returns 0 */
			nr = 0;
		else
		{
			nr = numberof(scpvec, *list);
			sign = nr > 0 ? 1 : -1;
			nr = abs(nr);
		}
/* scpvec is already in list */
		if (nr <= list->n  &&  nr > 0)
		{
			vecnr = (*vec)[nr];
			for (i = 0; i < dim; ++i)
				vecnr[i] += sign * Vvj[i];
		}
/* scpvec is a new scalar product combination */
		else if (nr >= list->n+1)
		{
			++list->n;
			if ((list->v = (long long**)realloc(list->v, (list->n+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			if ((list->v[list->n] = (long long*)malloc(len * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
/* numberof changes the sign of scpvec if the first nonzero entry is < 0,
   hence this is correct */
			listvn = list->v[list->n];
			for (i = 0; i < len; ++i)
				listvn[i] = scpvec[i];
			if ((*vec = (long long**)realloc(*vec, (list->n+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			if (((*vec)[list->n] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
			vecn = (*vec)[list->n];
			for (i = 0; i < dim; ++i)
				vecn[i] = sign * Vvj[i];
/* shuffle the new vector to the correct position, this should be quick enough,
   since the length of the list should not exceed some hundreds */
			listv = list->v;
			for (i = list->n; i > nr+1-list->n; --i)
			{
				tmp = listv[i];
				listv[i] = listv[i-1];
				listv[i-1] = tmp;
				tmp = (*vec)[i];
				(*vec)[i] = (*vec)[i-1];
				(*vec)[i-1] = tmp;
			}
		}
	}
	free(scpvec);
}

void base(scpcomb *com, long long ***b, long long **v, long long **F, long long dim)
/*****	computes a basis b for the lattice 
      generated by the vectors in v, puts a 
      transformation matrix on T, 
      i.e. b = T*v, uses LLL-reduction with 	
      respect to F	*****/
{
	long long	i, j, k, nv, rank, *perm, tmp;
	long long *Fvi, **Fv, *vi, **f, *fj, **tr, *comtr, *norm;
	long long TheMaxNorm, max;
	TheMaxNorm=GetMAXNORM();
	nv = com->list.n + 1;
	max = 1;
/* get the maximal entry in the vector sums */
	for (i = 0; i < nv; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			if (abs(v[i][j]) > max)
				max = abs(v[i][j]);
		}
	}
/* Fv is the list of products of F with the vectors in v, scalar products
   with respect to F can be calculated as standard scalar products of
   vectors in v with vectors in Fv */
	if ((Fv = (long long**)malloc(nv * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
/* the product of the maximal entry in the vector sums with the maximal entry
   in the product of the Gram-matrices with the vector sums should not exceed
   MAXNORM to avoid overflow */
	max = TheMaxNorm / max;
	for (i = 0; i < nv; ++i)
	{
		if ((Fv[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
		Fvi = Fv[i];
		for (j = 0; j < dim; ++j)
		{
			Fvi[j] = sscp(F[j], v[i], dim);
			if (abs(Fvi[j]) > max)
/* an entry in F.v[i] is too large */
			{
			  fprintf(stderr, "Error: Found entry ");
			  single2Print(stderr, Fvi[j]);
			  fprintf(stderr," in vector sum.\nTo avoid overflow, the entries should not exceed ");
			  single2Print(stderr, max);
			  fprintf(stderr, ".\nTry without option -D or use -D with higher value.\n");
			  exit (3);
			}
		}
	}
/* b is the basis for the lattice */
	if ((*b = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	if (((*b)[0] = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* com->trans is the transformation matrix */
	if ((com->trans = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	if ((com->trans[0] = (long long*)malloc(nv * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* f is the Gram-matrix for the basis b with respect to the form F */
	if ((f = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	if ((f[0] = (long long*)malloc(1 * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* tr is the transformation matrix for the LLL-reduced lattice */
	if ((tr = (long long**)malloc(1 * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	if ((tr[0] = (long long*)malloc(1 * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* perm is the new order in which the vectors in v are added to the lattice
   generated by the preceding vectors */
	if ((perm = (long long*)malloc(nv * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
	if ((norm = (long long*)malloc(nv * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
/* bubble sort with respect to the norm of the vectors in v with respect to the
   Gram-matrix F, which is assumed to be positive definite, 
   this should stabilize the LLL-reduction, 
   the bubble sort should be sufficient since the number of vectors here is 
   assumed to be rather small (at most some hundreds) */
	for (i = 0; i < nv; ++i)
	{
		norm[i] = sscp(v[i], Fv[i], dim);
		if (norm[i] > TheMaxNorm)
		{
			exit (3);
		}
		perm[i] = i;
	}
	i = 0;
	while (i < nv-1)
	{
		if (norm[perm[i]] > norm[perm[i+1]])
		{
			tmp = perm[i];
			perm[i] = perm[i+1];
			perm[i+1] = tmp;
			if (i > 0)
				--i;
		}
		else
			++i;
	}
	free(norm);
/* jump over possible 0-vectors */
	i = 0;
	for (j = 0; j < dim  &&  v[perm[i]][j] == 0; ++j);
	while (j == dim  &&  i < nv)
	{
		++i;
		if (i < nv)
			for (j = 0; j < dim  &&  v[perm[i]][j] == 0; ++j);
	}
	rank = 0;
	while (i < nv)
	{
/* v[perm[i]] is the candidate to enlarge the lattice generated by the 
   base vectors in b */
		vi = v[perm[i]];
		Fvi = Fv[perm[i]];
		for (j = 0; j < rank; ++j)
		{
		  f[j][rank] = sscp((*b)[j], Fvi, dim);
		  f[rank][j] = f[j][rank];
		}
		f[rank][rank] = sscp(vi, Fvi, dim);
		if (lll(f, tr, rank+1) > rank)
/* the rank of the lattice generated by v[perm[i]] and the vectors in b is 
   rank+1, i.e. one higher than without v[perm[i]] */
		{
			++rank;
			for (j = 0; j < dim; ++j)
				(*b)[rank-1][j] = vi[j];
			comtr = com->trans[rank-1];
			for (j = 0; j < nv; ++j)
				comtr[j] = 0;
			comtr[perm[i]] = 1;
/* transform basis and transformation matrix */
			change(*b, rank, tr, dim);
			change(com->trans, rank, tr, nv);
			if ((*b = (long long**)realloc(*b, (rank+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			if (((*b)[rank] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
			if ((com->trans = (long long**)realloc(com->trans, (rank+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			if ((com->trans[rank] = (long long*)malloc(nv * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
			if ((f = (long long**)realloc(f, (rank+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			for (j = 0; j < rank; ++j)
			{
				if ((f[j] = (long long*)realloc(f[j], (rank+1) * sizeof(long long))) == 0)
					exit (EXIT_FAILURE);
			}
			if ((f[rank] = (long long*)malloc((rank+1) * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
			for (j = 0; j < rank; ++j)
			{
				fj = f[j];
				for (k = 0; k <= j; ++k)
				{
					fj[k] = scp((*b)[j], F, (*b)[k], dim);
					f[k][j] = fj[k];
				}
			}
			if ((tr = (long long**)realloc(tr, (rank+1) * sizeof(long long*))) == 0)
				exit (EXIT_FAILURE);
			for (j = 0; j < rank; ++j)
			{
				if ((tr[j] = (long long*)realloc(tr[j], (rank+1) * sizeof(long long))) == 0)
					exit (EXIT_FAILURE);
			}
			if ((tr[rank] = (long long*)malloc((rank+1) * sizeof(long long))) == 0)
				exit (EXIT_FAILURE);
		}
		else if (abs(tr[rank][rank]) > 1)
/* v[perm[i]] enlarges the lattice generated by the vectors in b by a finite
   index |tr[rank][rank]| */
		{
			for (j = 0; j < dim; ++j)
				(*b)[rank][j] = vi[j];
			comtr = com->trans[rank];
			for (j = 0; j < nv; ++j)
				comtr[j] = 0;
			comtr[perm[i]] = 1;
/* transform basis and transformation matrix */
			change(*b, rank+1, tr, dim);
			change(com->trans, rank+1, tr, nv);
			for (j = 0; j < rank; ++j)
			{
				fj = f[j];
				for (k = 0; k <= j; ++k)
				{
					fj[k] = scp((*b)[j], F, (*b)[k], dim);
					f[k][j] = fj[k];
				}
			}
		}
		++i;
	}
	com->rank = rank;
	for (i = 0; i < nv; ++i)
		free(Fv[i]);
	free(Fv);
	for (i = 0; i <= rank; ++i)
	{
		free(f[i]);
		free(tr[i]);
	}
	free(f);
	free(tr);
	free(perm);
}

void coef(scpcomb *com, long long **b, long long **v, long long **F, long long dim)
/*****	expresses the vectors in v in the 
      basis b, puts the transformation matrix
      on com->coef, i.e. v = com->coef * b,	
      uses LLL-reduction with respect to F 
      to obtain the coefficients	*****/
{
	long long	i, j, nb, nv, sign;
	long long **Fv, *fnb, **f, *Fvi, *fi, **Fb, **tr, *trnb;
	long long *comci;
	if ((com->coef = (long long**)malloc((com->list.n+1) * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i <= com->list.n; ++i)
	{
		if ((com->coef[i] = (long long*)malloc(com->rank * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
	}
	nb = com->rank;
	nv = com->list.n + 1;
/* Fv is the list of products of F with the vectors in v, scalar products
   with respect to F can be calculated as standard scalar products of
   vectors in v with vectors in Fv */
	if ((Fv = (long long**)malloc(nv * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i < nv; ++i)
	{
	  if ((Fv[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
	    exit (EXIT_FAILURE);
	  for (j = 0; j < dim; ++j)
	    Fv[i][j] = sscp(F[j], v[i], dim);
	}
/* Fb is the list of products of F with the vectors in b */
	if ((Fb = (long long**)malloc(nb * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i < nb; ++i)
	{
		if ((Fb[i] = (long long*)malloc(dim * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
		for (j = 0; j < dim; ++j)
			Fb[i][j] = sscp(F[j], b[i], dim);
	}
	if ((f = (long long**)malloc((nb+1) * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i <= nb; ++i)
	{
		if ((f[i] = (long long*)malloc((nb+1) * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
	}
	for (i = 0; i < nb; ++i)
	{
		fi = f[i];
		for (j = 0; j <= i; ++j)
		{
			fi[j] = sscp(b[i], Fb[j], dim);
			f[j][i] = fi[j];
		}
	}
	if ((tr = (long long**)malloc((nb+1) * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (i = 0; i <= nb; ++i)
	{
		if ((tr[i] = (long long*)malloc((nb+1) * sizeof(long long))) == 0)
			exit (EXIT_FAILURE);
	}
	fnb = f[nb];
	trnb = tr[nb];
	for (i = 0; i < nv; ++i)
	{
		Fvi = Fv[i];
		for (j = 0; j < nb; ++j)
		{
			f[j][nb] = sscp(b[j], Fvi, dim);
			fnb[j] = f[j][nb];
		}
		fnb[nb] = sscp(v[i], Fvi, dim);
/* if the vector v[i] is in the lattice generated by the base b, the rank of the
   Gram-matrix f must be nb and the last row of the transformation matrix tr
   expresses the 0-vector, in particular |tr[nb][nb]| must be <=1, since 
   otherwise the vector v[i] would enlarge the lattice be an index 
   |tr[nb][nb]| */
		if (lll(f, tr, nb+1) > nb)
		{
			fprintf(stderr, "Error: vector lies not in lattice, Case 1\n");
			exit (3);
		}
		else if (abs(trnb[nb]) > 1)
		{
			fprintf(stderr, "Error: vector lies not in lattice, Case 2\n");
			exit (3);
		}
		else
		{
		  comci = com->coef[i];
		  sign = trnb[nb];
		  for (j = 0; j < nb; ++j)
		    comci[j] = -sign * trnb[j];
		}
	}
	for (i = 0; i < nv; ++i)
		free(Fv[i]);
	free(Fv);
	for (i = 0; i < nb; ++i)
		free(Fb[i]);
	free(Fb);
	for (i = 0; i <= nb; ++i)
	{
		free(f[i]);
		free(tr[i]);
	}
	free(f);
	free(tr);
}

void scpforms(scpcomb *com, long long **b, invar F)
/*****	com->F[i] is the Gram-matrix of the
      basis b with respect to F.A[i]	*****/
{
	long long	i, j, k, dim, nb;
	long long **FAi, **Fbi, *comFij;

	dim = F.dim;
	nb = com->rank;
	if ((com->F = (long long***)malloc(F.n * sizeof(long long**))) == 0)
		exit (EXIT_FAILURE);
/* Fbi is the list of products of F.A[i] with the vectors in b */
	if ((Fbi = (long long**)malloc(nb * sizeof(long long*))) == 0)
		exit (EXIT_FAILURE);
	for (j = 0; j < nb; ++j)
	{
	  if ((Fbi[j] = (long long*)malloc(dim * sizeof(long long))) == 0)
	    exit (EXIT_FAILURE);
	}
	for (i = 0; i < F.n; ++i)
	{
	  if ((com->F[i] = (long long**)malloc(nb * sizeof(long long*))) == 0)
	    exit (EXIT_FAILURE);
	  FAi = F.A[i];
	  for (j = 0; j < nb; ++j)
	    {
	      for (k = 0; k < dim; ++k)
		Fbi[j][k] = sscp(FAi[k], b[j], dim);
	    }
	  for (j = 0; j < nb; ++j)
	    {
	      if ((com->F[i][j] = (long long*)malloc(nb * sizeof(long long))) == 0)
		exit (EXIT_FAILURE);
	      comFij = com->F[i][j];
	      for (k = 0; k < nb; ++k)
		comFij[k] = sscp(b[j], Fbi[k], dim);
	    }
	}
	for (j = 0; j < nb; ++j)
	  free(Fbi[j]);
	free(Fbi);
}
