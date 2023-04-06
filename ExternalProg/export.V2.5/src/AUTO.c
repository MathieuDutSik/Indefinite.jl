/*****	Main program for the automorphism program AUTO	*****/

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
	group		G;
	invar		F;
	veclist		V, norm;
	fpstruct	fp;
	long long		dim, fail;
	long long		nH, ngen;
	long long		i, j, k, l, n;
	long long      max;
	long long      **FAi, ***H, *Vvi, *Vvj, **Fvi, *Fvij;
	long long      **sumveclist, **sumvecbase, *vec;
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
/* read the invariant forms */
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
	if (flags.GEN == 1)
/* some automorphisms are already known */
	{
	  ngen = getnumber(infile);
	  if ((H = (long long***)malloc(ngen * sizeof(long long**))) == 0)
	    exit (1);
	  nH = 0;
	  fail = 0;
	  while (nH+fail < ngen)
	    {
	      n = getmat(infile, &H[nH]);
	      if (n != dim)
		{
		  fprintf(stderr, "Error: dimension %lld should be %lld\n", n, dim);
		  exit (2);
		}
	      /* check whether the matrix is really an automorphism, i.e. fixes the forms */
	      if (checkgen(H[nH], F) == 0)
		/* the matrix is not an automorphism */
		{
		  ++fail;
		  for (i = 0; i < dim; ++i)
		    free(H[nH][i]);
		}
	      else
		/* the matrix fixes all forms in F */
		++nH;
	    }
	}
	else
		nH = 0;
	V.dim = dim;
	V.len = dim + F.n;
/* get the short vectors */
	if (flags.VEC > 0)
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
	fclose(infile);
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
	  fprintf(stderr, " in short vectors.\nTo avoid overflow, the entries should not exceed ");
	  single2Print(stderr, TheMaxEntry);
	  fprintf(stderr, ".\n");
	  exit (3);
	}
/* V.prime is a prime p, such that every entry x in the short vectors remains
   unchanged under reduction mod p in symmetric form, i.e. -p/2 < x <= p/2 */
	for (V.prime = 2*max + 1; isprime(V.prime) == 0; ++V.prime);
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
/* delete those vectors, which can not be the image of any of the standard-base
   vectors */
	checkvecs(&V, F, norm);
	for (i = 1; i <= norm.n; ++i)
		free(norm.v[i]);
	free(norm.v);
/* print the invariant forms, if output is in ASCII-format */
	if (flags.OUTPUT == 0)
	{
		printf("#%lld\n", F.n);
		for (i = 0; i < F.n; ++i)
			putform(F.A[i], dim);
	}
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
/* fp.per is the order in which the images of the base-vectors are set */
	if ((fp.per = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* fp.e[i] is the index in V.v of the i-th vector of the standard-base in the 
   new order */
	if ((fp.e = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* fp.diag is the diagonal of the fingerprint in the new order, fp.diag[i] is
   an upper bound for the length of the orbit of the i-th base-vector under
   the stabilizer of the preceding ones */
	if ((fp.diag = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* compute the fingerprint */
	fingerprint(&fp, F, V);
	if (flags.PRINT == 1)
/* if the -P option is given, print the diagonal fp.diag and the number of short
   vectors on AUTO.tmp */
	{
		outfile = fopen("AUTO.tmp", "a");
		fprintf(outfile, "%lld short vectors\n", 2*V.n);
		fprintf(outfile, "fingerprint diagonal:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %lld", fp.diag[i]);
		fprintf(outfile, "\n");
		fprintf(outfile, "permutation:\n");
		for (i = 0; i < dim; ++i)
			fprintf(outfile, " %lld", fp.per[i]);
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
/* the standard-base must occur in the set of short vectors, since Id is
   certainly an automorphism */
		{
			fprintf(stderr, "Error: standard basis vector nr. %lld not found\n", i);
			exit (2);
		}
	}
	free(vec);
/* if the -D option is given, the scalar product combinations and the 
   corresponding vector sums are computed for the standard-basis */
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
   vector sums on level i on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
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
/* the Bacher-polynomials for the first BACHDEP base-vectors are computed,
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
/* if the -P option is given, print the Bacher-polynomial on AUTO.tmp */
			{
				outfile = fopen("AUTO.tmp", "a");
				fprintf(outfile, "Bacher-polynomial of base vector fp.e[%lld]:\n", i);
				fputbach(outfile, bach[i]);
				fclose(outfile);
			}
		}
	}
/* set up the group: the generators in G.g[i] are matrices that fix the
   first i base-vectors but do not fix fp.e[i] */
	if ((G.g = (long long****)malloc(dim * sizeof(long long***))) == 0)
		exit (1);
/* G.ng is the number of generators in G.g[i] */
	if ((G.ng = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* the first G.nsg[i] generators in G.g[i] are obtained as stabilizer elements
   and are not necessary to generate the group */
	if ((G.nsg = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
/* G.ord[i] is the orbit length of fp.e[i] under the generators in 
   G.g[i] ... G.g[dim-1] */
	if ((G.ord = (long long*)malloc(dim * sizeof(long long))) == 0)
		exit (1);
	for (i = 0; i < dim; ++i)
	{
	  if ((G.g[i] = (long long***)malloc(1 * sizeof(long long**))) == 0)
	    exit (1);
	  G.ng[i] = 0;
	  G.nsg[i] = 0;
	}
	G.dim = dim;
	if ((G.g[0][0] = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit (1);
/* -Id is always an automorphism */
	for (i = 0; i < dim; ++i)
	{
	  if ((G.g[0][0][i] = (long long*)malloc(dim * sizeof(long long))) == 0)
	    exit (1);
	  for (j = 0; j < dim; ++j)
	    G.g[0][0][i][j] = 0;
	  G.g[0][0][i][i] = -1;
	  G.ng[0] = 1;
	}
/* fill in the generators which were given as input */
	for (i = 0; i < nH; ++i)
	{
	  for (j = 0; j < dim  &&  operate(fp.e[j], H[i], V) == fp.e[j]; ++j);
	  if (j < dim)
	    /* the generator is not Id and fixes fp.e[0],...,fp.e[j-1] but not fp.e[j] */
	    {
	      ++G.ng[j];
	      if ((G.g[j] = (long long***)realloc(G.g[j], G.ng[j] * sizeof(long long**))) == 0)
		exit (1);
	      if ((G.g[j][G.ng[j]-1] = (long long**)malloc(dim * sizeof(long long*))) == 0)
		exit (1);
	      for (k = 0; k < dim; ++k)
		{
		  if ((G.g[j][G.ng[j]-1][k] = (long long*)malloc(dim * sizeof(long long))) == 0)
		    exit (1);
		  for (l = 0; l < dim; ++l)
		    G.g[j][G.ng[j]-1][k][l] = H[i][k][l];
		}
	    }
	  for (k = 0; k < dim; ++k)
	    free(H[i][k]);
	  free(H[i]);
	}
	if (nH > 0)
	  free(H);
	nH = 0;
	for (i = 0; i < dim; ++i)
	  nH += G.ng[i];
	if ((H = (long long***)malloc(nH * sizeof(long long**))) == 0)
	  exit (1);
/* calculate the orbit lengths under the automorphisms known so far */
	for (i = 0; i < dim; ++i)
	{
		if (G.ng[i] > 0)
		{
			nH = 0;
			for (j = i; j < dim; ++j)
			{
				for (k = 0; k < G.ng[j]; ++k)
				{
					H[nH] = G.g[j][k];
					++nH;
				}
			}
			G.ord[i] = orbitlen(fp.e[i], fp.diag[i], H, nH, V);
		}
		else
			G.ord[i] = 1;
	}
	free(H);
/* compute the full automorphism group */
	autom(&G, V, F, fp, comb, bach, flags);
/* print the generators of the group, which are necessary for generating */
	putgens(G, flags);
	n = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
		if (G.ord[i] > 1)
			++n;
	}
/* for a base in the sense of permutation groups (a set of vectors the point
   stabilizer of which is trivial) one needs n points (standard-base vectors) */
	k = 0;
	if (flags.OUTPUT >= 1)
/* for GAP or MAGMA print the vectors of such a base */
	{
	  printf("BaseBG := [\n");
	  for (i = flags.STAB; i < dim; ++i)
	    {
	      if (G.ord[i] > 1)
		{
		  printf("[");
		  for (j = 0; j < dim-1; ++j)
		    {
		      fprintf(stdout, " ");
		      single2Print(stdout, V.v[fp.e[i]][j]);
		    }
		  if (k < n-1)
		    {
		      fprintf(stdout, " ");
		      single2Print(stdout, V.v[fp.e[i]][dim-1]);
		      fprintf(stdout, "],\n");
		    }
		  else
		    {
		      fprintf(stdout, " ");
		      single2Print(stdout, V.v[fp.e[i]][dim-1]);
		      fprintf(stdout, "]];\n");
		    }
		  ++k;
		}
	    }
	}
/* print a base to AUTO.tmp if the print-flag is set */
	if (flags.PRINT == 1)
	{
	  outfile = fopen("AUTO.tmp", "a");
	  fprintf(outfile, "%lldx%lld\n", n, dim);
	  for (i = flags.STAB; i < dim; ++i)
	    {
	      if (G.ord[i] > 1)
		{
		  for (j = 0; j < dim; ++j)
		    {
		      fprintf(outfile, " ");
		      single2Print(outfile, V.v[fp.e[i]][j]);
		    }
		  fprintf(outfile, "\n");
		}
	    }
	  fclose(outfile);
	}
/* print the order of the group */
	putord(G, flags);
	return 0;
}
