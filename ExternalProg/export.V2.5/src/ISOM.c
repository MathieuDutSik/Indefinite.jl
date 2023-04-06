/*****	Main program for the isometry program ISOM	*****/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "types.h"

int main(int argc, char *argv[])
{
	FILE		*infile, *outfile;
	char		*filename;
	bachpol		*bach;
	flagstruct	flags;
	scpcomb		*comb;
	invar		F, FF;
	veclist		V, norm;
	fpstruct	fp;
	long long		dim, fail;
	long long		nG, ngen;
	long long		i, j, k, n, nV1;
	long long      max;
	long long      **FAi, *Vvi, *Vvj, **Fvi, *Fvij;
	long long      **sumveclist, **sumvecbase, *vec, ***G;
	long long      TheMaxNorm, TheMaxEntry;
	TheMaxNorm=GetMAXNORM();
	TheMaxEntry=GetMAXENTRY();

/* get the flags from the command line */
	getflags(&flags, argc, argv);
	if (argc > 1  &&  argv[argc-1][0] != '-')
/* if the last argument is not a flag it is assumed to be the input file */
	{
		filename = argv[argc-1];
		if ((infile = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr, "Error: can not open file %s\n", filename);
			exit (2);
		}
	}
	else
		infile = stdin;
/* F.n is the number of invariant forms */
	F.n = getnumber(infile);
	if ((F.A = (long long***)malloc(F.n * sizeof(long long**))) == 0)
		exit (1);
/* read the invariant forms of the first lattice */
	F.dim = getmat(infile, &F.A[0]);
	dim = F.dim;
	for (i = 1; i < F.n; ++i)
	{
		n = getmat(infile, &F.A[i]);
		if (n != dim)
/* all invariant forms should have the same dimension */
		{
			fprintf(stderr, "Error: dimension %lld should be %lld\n", n, dim);
			exit (2);
		}
	}
	V.dim = dim;
	V.len = dim + F.n;
/* get the short vectors of the first lattice */
	if (flags.VEC == 1  ||  flags.VEC == 3)
/* the short vectors are read from the input stream infile */
		getvecs(infile, &V);
	else
/* the short vectors are calculated (which is the default) */
	{
		max = F.A[0][0][0];
		for (i = 1; i < dim; ++i)
		{
			if (F.A[0][i][i] > max)
				max = F.A[0][i][i];
		}
/* the vectors with length up to the maximum of the diagonal of the first
   form (which is assumed to be positive definite) are required */
		shortvec(F.A[0], max, &V);
	}
/* the first nonzero entry of each vector is made positive and the maximal 
   entry in the vectors is determined */
	max = 1;
	for (i = 1; i <= V.n; ++i)
	{
		Vvi = V.v[i];
		for (j = 0; j < dim  &&  Vvi[j] == 0; ++j);
		if (j < dim  &&  Vvi[j] < 0)
		{
			for (k = j; k < dim; ++k)
				Vvi[k] *= -1;
		}
		for (k = j; k < dim; ++k)
		{
			if (abs(Vvi[k]) > max)
				max = abs(Vvi[k]);
		}
	}
	if (max > TheMaxEntry)
/* there is an entry, which is bigger than MAXENTRY, which might cause overflow,
   since the program doesn't use long integer arithmetic */
	{
	  fprintf(stderr, "Error: Found entry ");
	  single2Print(stderr, max);
	  fprintf(stderr, " in short vectors of the first lattice.\nTo avoid overflow, the entries should not exceed ");
	  single2Print(stderr, TheMaxEntry);
	  fprintf(stderr, ".\n");
	  exit (3);
	}
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base */
	if ((norm.v = (long long**)malloc((dim+1) * sizeof(long long*))) == 0)
		exit (1);
	norm.n = dim;
	norm.dim = F.n;
	norm.len = F.n;
	for (i = 1; i <= norm.n; ++i)
	{
/* norm.v[i] is the norm combination of the (i-1)-th base-vector */
	  if ((norm.v[i] = (long long*)malloc(norm.dim * sizeof(long long))) == 0)
	    exit (1);
	  for (k = 0; k < norm.dim; ++k)
	    norm.v[i][k] = F.A[k][i-1][i-1];
	}
	sortvecs(&norm);
/* delete those vectors, which can not be the image of any vector of the
   standard-base */
	checkvecs(&V, F, norm);
/* F.v[i][j] is the transposed of the product of the Gram-matrix F.A[i] with the
   transposed of the vector v[j], hence the scalar product of v[j] and v[k] with
   respect to the Gram-matrix F.A[i] can be computed as standard scalar product
   of v[j] and F.v[i][k] */
	if ((F.v = (long long***)malloc(F.n * sizeof(long long**))) == 0)
		exit (1);
/* the product of the maximal entry in the short vectors with the maximal entry
   in F.v[i] should not exceed MAXNORM to avoid overflow */
	max = TheMaxNorm / max;
	for (i = 0; i < F.n; ++i)
	{
		FAi = F.A[i];
		if ((F.v[i] = (long long**)malloc((V.n+1) * sizeof(long long*))) == 0)
			exit (1);
		Fvi = F.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (1);
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in F.v[i] is too large */
				{
				  fprintf(stderr, "Error: Found entry ");
				  single2Print(stderr, Fvij[k]);
				  fprintf(stderr, " in F.v.\nTo avoid overflow, the entries should not exceed ");
				  single2Print(stderr, max);
				  fprintf(stderr, ".\n");
				  exit (3);
				}
			}
		}
	}
/* fp.per is the order in which the images of the base vectors are chosen */
	if ((fp.per = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	if ((fp.e = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* fp.diag is the diagonal of the fingerprlong long in the new order */
	if ((fp.diag = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* compute the fingerprint */
	fingerprint(&fp, F, V);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of
   short vectors on ISOM.tmp */
	{
		outfile = fopen("ISOM.tmp", "a");
		fprintf(outfile, "%lld short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %lld", fp.diag[i]);
		fprintf(outfile, "\n");
		fclose(outfile);
	}
/* get the standard basis in the new order */
	if ((vec = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < dim; ++i)
	{
		for (j = 0; j < dim; ++j)
			vec[j] = 0;
		vec[fp.per[i]] = 1;
		fp.e[i] = numberof(vec, V);
		if (abs(fp.e[i]) > V.n)
		{
			fprintf(stderr, "Error: standard basis vector nr. %lld not found\n", i);
			exit (2);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the basis of the first lattice */
	if (flags.DEPTH > 0)
	{
		if ((comb = (scpcomb*)malloc(dim * sizeof(scpcomb))) == 0)
			exit (1);
		for (i = 0; i < dim; ++i)
		{
/* compute the list of scalar product combinations and the corresponding
   vector sums */
			scpvecs(&comb[i].list, &sumveclist, i, fp.e, flags.DEPTH, V, F);
/* compute a basis for the lattice that is generated by the vector sums and
   a transformation matrix that expresses the basis in terms of the 
   vector sums */
			base(&comb[i], &sumvecbase, sumveclist, F.A[0], dim);
			if (flags.PRINT == 1)
/* if the -P option is given, print the rank of the lattice generated by the
   vector sums on level i on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
				fprintf(outfile, "comb[%lld].rank = %lld\n", i, comb[i].rank);
				fclose(outfile);
			}
/* compute the coefficients of the vector sums in terms of the basis */
			coef(&comb[i], sumvecbase, sumveclist, F.A[0], dim);
			for (j = 0; j <= comb[i].list.n; ++j)
				free(sumveclist[j]);
			free(sumveclist);
/* compute the scalar products of the base-vectors */
			scpforms(&comb[i], sumvecbase, F);
			for (j = 0; j < comb[i].rank; ++j)
				free(sumvecbase[j]);
			free(sumvecbase);
		}
	}
/* the Bacher-polynomials for the first BACHDEP base vectors are computed,
   if no scalar product was given as an argument, the scalar product is set
   to 1/2 the norm of the base-vector (with respect to the first form) */
	if (flags.BACH[0] == 1)
	{
		if ((bach = (bachpol*)malloc(flags.BACHDEP * sizeof(bachpol))) == 0)
			exit (1);
		for (i = 0; i < flags.BACHDEP; ++i)
		{
			if (flags.BACH[2] == 0)
/* no scalar product was given as an option */
				flags.BACHSCP = V.v[fp.e[i]][dim] / 2;
/* compute the Bacher-polynomial */
			bacher(&bach[i], fp.e[i], flags.BACHSCP, V, F.v[0]);
			if (flags.PRINT == 1)
/* if the -P option is given, print the Bacher-polynomial on ISOM.tmp */
			{
				outfile = fopen("ISOM.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%lld]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* the vectors of the first lattice are no longer required, only their number */
	nV1 = V.n;
	for (i = 0; i <= V.n; ++i)
		free(V.v[i]);
	free(V.v);
	for (i = 0; i < F.n; ++i)
	{
		for (j = 1; j <= V.n; ++j)
			free(F.v[i][j]);
		free(F.v[i]);
	}
	free(F.v);
/* FF are the invariant forms of the second lattice */
	if ((FF.A = (long long***)malloc(F.n * sizeof(long long**))) == 0)
		exit (1);
	FF.n = F.n;
	FF.dim = dim;
/* read the invariant forms of the second lattice */
	for (i = 0; i < FF.n; ++i)
	{
		n = getmat(infile, &FF.A[i]);
		if (n != dim)
/* all invariant forms should have the same dimension */
		{
			fprintf(stderr, "Error: dimension %lld should be %lld\n", n, dim);
			exit (2);
		}
	}
/* get the short vectors of the second lattice */
	if (flags.VEC == 2  ||  flags.VEC == 3)
/* the short vectors are read from the input stream infile */
		getvecs(infile, &V);
	else
/* the short vectors are calculated (which is the default) */
	{
		max = F.A[0][0][0];
		for (i = 1; i < dim; ++i)
		{
			if (F.A[0][i][i] > max)
				max = F.A[0][i][i];
		}
/* the vectors with length up to the maximum of the diagonal of the first
   form of the first lattice are required */
		shortvec(FF.A[0], max, &V);
	}
/* the first nonzero entry of each vector is made positive and the maximal
   entry of the vectors is determined */
	max = 1;
	for (i = 1; i <= V.n; ++i)
	{
		Vvi = V.v[i];
		for (j = 0; j < dim  &&  Vvi[j] == 0; ++j);
		if (j < dim  &&  Vvi[j] < 0)
		{
			for (k = j; k < dim; ++k)
				Vvi[k] *= -1;
		}
		for (k = j; k < dim; ++k)
		{
			if (abs(Vvi[k]) > max)
				max = abs(Vvi[k]);
		}
	}
	if (max > TheMaxEntry)
/* there is an entry, which is bigger than MAXENTRY, which might cause overflow,
   since the program doesn't use long integer arithmetic */
	{
	  fprintf(stderr, "Error: Found entry ");
	  single2Print(stderr, max);
	  fprintf(stderr, " in short vectors of the second lattice.\nTo avoid overflow, the entries should not exceed ");
	  single2Print(stderr, TheMaxEntry);
	  fprintf(stderr, ".\n");
	  exit (3);
	}
/* V.prime is a prime p, such that the entries of the short vectors remain 
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
	k = F.A[0][0][0];
	for (i = 1; i < dim; ++i)
	{
		if (F.A[0][i][i] > k)
			k = F.A[0][i][i];
	}
	for (i = 0; i < dim; ++i)
	{
		if (FF.A[0][i][i] > k)
			V.prime = DEF_PRIME;
	}
/* sort the vectors and delete doublets */
	sortvecs(&V);
/* the norm-vector (i.e. the vector of the norms with respect to the different
   invariant forms) of each vector must be equal to the norm-vector of one
   of the vectors from the standard-base of the first lattice */
	checkvecs(&V, FF, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
	if (V.n != nV1)
	{
		printf("The lattices are not isomorphic !\nDifferent numbers of short vectors (%lld|%lld)\n", nV1, V.n);
		exit (0);
	}
/* FF.v[i][j] is the transposed of the product of the Gram-matrix FF.A[i] 
   with the transposed of the vector V.v[j], which now is a short vector of
   the second lattice */
	if ((FF.v = (long long***)malloc(FF.n * sizeof(long long**))) == 0)
		exit (1);
/* the product of the maximal entry in the short vectors with the maximal entry
   in FF.v[i] should not exceed MAXNORM to avoid overflow */
	max = TheMaxNorm / max;
	for (i = 0; i < FF.n; ++i)
	{
		FAi = FF.A[i];
		if ((FF.v[i] = (long long**)malloc((V.n+1) * sizeof(long long*))) == 0)
			exit (1);
		Fvi = FF.v[i];
		for (j = 1; j <= V.n; ++j)
		{
			Vvj = V.v[j];
			if ((Fvi[j] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (1);
			Fvij = Fvi[j];
			for (k = 0; k < dim; ++k)
			{
				Fvij[k] = sscp(FAi[k], Vvj, dim);
				if (abs(Fvij[k]) > max)
/* some entry in FF.v[i] is too large */
				{
				  fprintf(stderr, "Error: Found entry ");
				  single2Print(stderr, Fvij[k]);
				  fprintf(stderr, " in FF.v.\nTo avoid overflow, the entries should not exceed ");
				  single2Print(stderr, max);
				  fprintf(stderr, ".\n");
				  exit (3);
				}
			}
		}
	}
/* G are the automorphisms of the second lattice */
	if ((G = (long long***)malloc(1 * sizeof(long long**))) == 0)
		exit (1);
/* nG is the number of automorphisms of the second lattice */
	nG = 0;
	if (flags.GEN == 1)
/* get the automorphisms, which are already known */
	{
		ngen = getnumber(infile);
		if ((G = (long long***)realloc(G, ngen * sizeof(long long**))) == 0)
			exit (1);
		fail = 0;
		while (nG+fail < ngen)
		{
			n = getmat(infile, &G[nG]);
			if (n != dim)
			{
				fprintf(stderr, "Error: dimension %lld should be %lld\n", n, dim);
				exit (2);
			}
/* check whether the matrix is really an automorphism, i.e. fixes the forms */
			if (checkgen(G[nG], FF) == 0)
/* the matrix is not an automorphism */
			{
				++fail;
				for (i = 0; i < dim; ++i)
					free(G[nG][i]);
			}
			else
/* the matrix fixes the forms in FF */
				++nG;
		}
	}
	if (nG == 0)
/* if no automorphisms are known at least take -Id */
	{
		if ((G[0] = (long long**)malloc(dim * sizeof(long long*))) == 0)
			exit (1);
		for (i = 0; i < dim; ++i)
		{
			if ((G[0][i] = (long long*)malloc(dim * sizeof(long long))) == 0)
				exit (1);
			for (j = 0; j < dim; ++j)
				G[0][i][j] = 0;
			G[0][i][i] = -1;
		}
		nG = 1;
	}
	fclose(infile);
/* now search for an isometry */
	isometry(V, F, FF, fp, G, nG, comb, bach, flags);
	return 0;
}
