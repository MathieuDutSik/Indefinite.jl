/*****	This file contains routines to calculate 
	and compare Bacher-polynomials	*****/

#include <stdio.h>
#include <stdlib.h>
#include "types.h"

void bacher(bachpol *pol, long long I, long long S, veclist	V, long long **Fv)
/*****	the Bacher-polynomial for v[I] with scalar product S
      is defined as follows: let list be the vectors which
      have the same length as v[I] and with scalar product S
      with v[I], for each vector w in list let n_w be the 
      number of pairs (y,z) of vectors in list, such that all
      scalar products between w,y and z are S, 
      then the Bacher-polynomial is the sum over the w in list
      of the monomials X^n_w	*****/
{
	long long	*list, nlist, *listxy, nxy, *counts;
	long long	i, j, k, dim, sign1, sign2;
	long long *vI, s;
/* the Bacher-polynomials of v[I] and -v[I] are equal */
	I = abs(I);
	dim = V.dim;
	vI = V.v[I];
/* list of vectors that have scalar product S with v[I] */
	if ((list = (long long*)malloc(2*V.n * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < 2*V.n; ++i)
		list[i] = 0;
	nlist = 0;
	for (i = 1; i <= V.n; ++i)
	{
	  if (V.v[i][dim] == vI[dim])
	    {
	      s = sscp(vI, Fv[i], dim);
	      if (s == S)
		{
		  list[nlist] = i;
		  ++nlist;
		}
	      if (-s == S)
		{
		  list[nlist] = -i;
		  ++nlist;
		}
	    }
	}
/* there are nlist vectors that have scalar product S with v[I] */
	pol->sum = nlist;
	if ((counts = (long long*)malloc(nlist * sizeof(long long))) == 0)
		exit (1);
	if ((listxy = (long long*)malloc(nlist * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < nlist; ++i)
	{
/* listxy is the list of the nxy vectors from list that have scalar product S 
   with v[list[i]] */
		for (j = 0; j < nlist; ++j)
			listxy[j] = 0;
		nxy = 0;
		sign1 = list[i] > 0 ? 1 : -1;
		for (j = 0; j < nlist; ++j)
		{
			sign2 = list[j] > 0 ? 1 : -1;
/* note: for i > 0 is v[-i] = -v[i] */
			if (sign1*sign2 * sscp(V.v[sign1*list[i]], Fv[sign2*list[j]], dim) == S)
			{
				listxy[nxy] = list[j];
				nxy += 1;
			}
		}
/* counts[i] is the number of pairs for the vector v[list[i]] */
		counts[i] = 0;
		for (j = 0; j < nxy; ++j)
		{
			sign1 = listxy[j] > 0 ? 1 : -1;
			for (k = j+1; k < nxy; ++k)
			{
				sign2 = listxy[k] > 0 ? 1 : -1;
				if (sign1*sign2 * sscp(V.v[sign1*listxy[j]], Fv[sign2*listxy[k]], dim) == S)
					counts[i] += 1;
			}
		}
	}
/* pol->maxd is the maximal degree of the Bacher-polynomial, 
   pol->mind the minimal degree */
	pol->maxd = counts[0];
	pol->mind = counts[0];
	for (i = 1; i < nlist; ++i)
	{
		if (counts[i] > pol->maxd)
			pol->maxd = counts[i];
		else if (counts[i] < pol->mind)
			pol->mind = counts[i];
	}
	if ((pol->coef = (long long*)malloc((pol->maxd - pol->mind + 1) * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i <= pol->maxd - pol->mind; ++i)
		pol->coef[i] = 0;
	for (i = 0; i < nlist; ++i)
		pol->coef[counts[i] - pol->mind] += 1;
/* the Bacher-polynomial is now: sum from i=pol->mind to pol->maxd over
   pol->coef[i - pol->mind] * X^i */
	free(listxy);
	free(counts);
	free(list);
}

long long bachcomp(bachpol pol, long long I, long long S, veclist V, long long **Fv)
/*****	checks, whether the vector v[I] has 
      the Bacher-polynomial pol	*****/
{
	long long	*co, *list, nlist, *listxy, nxy, count;
	long long	i, j, k, dim, sign1, sign2;
	long long *vI, s;
	I = abs(I);
	dim = V.dim;
	vI = V.v[I];
	if ((co = (long long*)malloc((pol.maxd - pol.mind + 1) * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i <= pol.maxd-pol.mind; ++i)
		co[i] = 0;
	if ((list = (long long*)malloc(pol.sum * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < pol.sum; ++i)
		list[i] = 0;
/* nlist should be equal to pol.sum */
	nlist = 0;
	for (i = 1; i <= V.n  &&  nlist <= pol.sum; ++i)
	{
		if (V.v[i][dim] == vI[dim])
		{
			s = sscp(vI, Fv[i], dim);
			if (s == S)
			{
				if (nlist < pol.sum)
					list[nlist] = i;
				nlist += 1;
			}
			if (-s == S)
			{
				if (nlist < pol.sum)
					list[nlist] = -i;
				nlist += 1;
			}
		}
	}
	if (nlist != pol.sum)
/* the number of vectors with scalar product S is already different */
	{
		free(co);
		free(list);
		return(0);
	}
/* listxy is the list of the nxy vectors from list that have scalar product S 
   with v[list[i]] */
	if ((listxy = (long long*)malloc(nlist * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < nlist; ++i)
	{
		for (j = 0; j < nlist; ++j)
			listxy[j] = 0;
		nxy = 0;
		sign1 = list[i] > 0 ? 1 : -1;
		for (j = 0; j < nlist; ++j)
		{
			sign2 = list[j] > 0 ? 1 : -1;
			if (sign1*sign2 * sscp(V.v[sign1*list[i]], Fv[sign2*list[j]], dim) == S)
			{
				listxy[nxy] = list[j];
				nxy += 1;
			}
		}
/* count is the number of pairs */
		count = 0;
		for (j = 0; j < nxy  &&  count <= pol.maxd; ++j)
		{
			sign1 = listxy[j] > 0 ? 1 : -1;
			for (k = j+1; k < nxy  &&  count <= pol.maxd; ++k)
			{
				sign2 = listxy[k] > 0 ? 1 : -1;
				if (sign1*sign2 * sscp(V.v[sign1*listxy[j]], Fv[sign2*listxy[k]], dim) == S)
					count += 1;
			}
		}
		if (count < pol.mind  ||  count > pol.maxd  || 
		    co[count-pol.mind] >= pol.coef[count-pol.mind])
/* if the number of pairs is smaller than pol.mind or larger than pol.maxd
   or if the coefficient of X^count becomes now larger than the one in pol,
   then the Bacher-polynomials can not be equal */
		{
			free(listxy);
			free(co);
			free(list);
			return(0);
		}
		else
			co[count-pol.mind] += 1;
	}
	free(listxy);
	free(co);
	free(list);
/* the Bacher-polynomials are equal */
	return(1);
}

void fputbach(FILE *outfile, bachpol pol)
/*****	prints a Bacher-polynomial	*****/
{
	long long	i;

	for (i = pol.mind; i <= pol.maxd; ++i)
	{
		if (pol.coef[i-pol.mind] != 0)
			fprintf(outfile, "%lld*x^%lld ", pol.coef[i-pol.mind], i);
	}
	fprintf(outfile, "\n");
}
