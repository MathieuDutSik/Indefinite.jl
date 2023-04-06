/*****	This file contains some routines for input/output	*****/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "types.h"

void single2Print(FILE *f, long long TheNb)
{
  fprintf(f, "%2lld", TheNb);
}
void singleScan(FILE *f, long long *TheNb)
{
  fscanf(f, "%lld", TheNb);
}
long long GetMAXNORM()
{
  long long a;
  long long i;
  a=1;
  for (i=1; i<=18; i++)
    {
      a=a*10;
    }
  return a;
}
long long GetMAXENTRY()
{
  long long a;
  long long i;
  a=1;
  for (i=1; i<=8; i++)
    {
      a=a*10;
    }
  return a;
}

void getflags(flagstruct *fl, int argc, char *argv[])
/*****	gets the options from the 
      command line	*****/
{
	long long	i;
	char	*str;

/* depth for the scalar product combinations */
	fl->DEPTH = 0;
/* only the point stabilizer of the first STAB basis-vectors will be computed */
	fl->STAB = 0;
/* flag that Bacher-polynomials will be used */
	fl->BACH[0] = 0;
/* flag that the depth for the Bacher-polynomials is given as an argument, 
   default is 1 */
	fl->BACH[1] = 0;
/* depth for the Bacher-polynomials */
	fl->BACHDEP = 0;
/* flag that the scalar product for the Bacher-polynomials is given as an 
   argument, default is 1/2*norm of the vector */
	fl->BACH[2] = 0;
/* scalar product for the Bacher-polynomials */
	fl->BACHSCP = 0;
/* flag that generators will be read */
	fl->GEN = 0;
/* flag that every new generator is immediately written on the file 
   AUTO.tmp */
	fl->PRINT = 0;
/* flag that the vectors will be read instead of calculated by the program */
	fl->VEC = 0;
/* flag for the output-style: 0 means ASCII, 1 GAP-format, 2 MAGMA-format */
	fl->OUTPUT = 0;
/* scan through the arguments */
	for (i = 1; i < argc; ++i)
	{
/* every option should start with a '-' */
	  if ((str = strchr(argv[i], '-')) != NULL)
	    {
	      if (strlen(str) <= 1)
		fprintf(stderr, "unknown option %s: ignored\n", str);
	      else if (str[1] == 'D')
/* option -Dn where n is some non-negative integer for depth of scalar product
   combinations */
		{
		  if (strlen(str) <= 2  ||  strcspn(str+2, "0123456789") > 0)
		    {
		      fprintf(stderr, "Error: no non-negative integer specified with -D option\n");
		      exit (2);
		    }
		  else
		    fl->DEPTH = atoi(str+2);
		}
	      else if (str[1] == 'S')
		/* option -Sn where n is some non-negative integer for n-point stabilizer */
		{
		  if (strlen(str) <= 2  ||  strcspn(str+2, "0123456789") > 0)
		    {
		      fprintf(stderr, "Error: no non-negative integer specified with -S option\n");
		      exit (2);
		    }
		  else
		    fl->STAB = atoi(str+2);
		}
	      else if (strlen(str) >= 3  &&  strncmp(str, "-BD", 3) == 0)
		/* option -BDn where n is some non-negative integer for depth of 
		   Bacher-polynomials */
		{
		  if (strlen(str) <= 3  ||  strcspn(str+3, "0123456789") > 0)
		    {
		      fprintf(stderr, "Error: no non-negative integer specified with -BD option\n");
		      exit (2);
		    }
		  else
		    {
		      fl->BACHDEP = atoi(str+3);
		      fl->BACH[0] = 1;
		      fl->BACH[1] = 1;
		    }
		}
	      else if (strlen(str) >= 3  &&  strncmp(str, "-BS", 3) == 0)
		/* option -BSn where n is some integer for scalar product of 
		   Bacher-polynomials */
		{
		  if (strlen(str) <= 3  ||  strcspn(str+3, "-0123456789") > 0)
		    {
		      fprintf(stderr, "Error: no integer specified with -BS option\n");
		      exit (2);
		    }
		  else
		    {
		      fl->BACHSCP = atoi(str+3);
		      fl->BACH[0] = 1;
		      fl->BACH[2] = 1;
		    }
		}
	      else if (strlen(str) == 2  &&  str[1] == 'B')
		/* option -B indicates that Bacher-polynomials will be used */
		fl->BACH[0] = 1;
	      else if (strlen(str) == 2  &&  str[1] == 'G')
		/* option -G indicates that some generators can be read from the input stream */
		fl->GEN = 1;
	      else if (strlen(str) == 3  &&  str[1] == 'V')
		/* option -V1 indicates that the short vectors for the first lattice are read 
		   from the input stream, option -V2 that the short vectors for the second
		   lattice are read and option -V3 that the short vectors for both lattices
		   are read (-V2 only makes sense in the isometry-program) */
		{
		  if (strpbrk(str+2, "123") == NULL)
		    {
		      fprintf(stderr, "Error: no integer between 1 and 3 specified with -V option\n");
		      exit (2);
		    }
		  else
		    fl->VEC = atoi(str+2);
		}
	      else if (strlen(str) == 2  &&  str[1] == 'V')
		/* option -V indicates that the short vectors are read from the input stream */
		{
		  if (fl->VEC == 0)
		    fl->VEC = 3;
		}
	      else if (strlen(str) == 2  &&  str[1] == 'P')
		/* option -P indicates that new generators will be written to the file AUTO.tmp
		   immediately */
		fl->PRINT = 1;
	      else if (strlen(str) == 3  &&  str[1] == 'O')
		/* option -OG indicates that the output is converted to GAP-format, -OM that it
		   is converted to MAGMA-format */
		{
		  if (strpbrk(str+2, "GM") == NULL)
		    {
		      fprintf(stderr, "Error: can only convert to GAP (-OG) or MAGMA (-OM)\n");
		      exit (2);
		    }
		  else if (str[2] == 'G')
		    fl->OUTPUT = 1;
		  else if (str[2] == 'M')
		    fl->OUTPUT = 2;
		}
	      else 
		fprintf(stderr, "unknown option %s: ignored\n", str);
	    }
	}
/* if Bacher-polynomials are to be used and no depth is given it is set to
   the default-value 1 */
	if (fl->BACH[0] == 1  &&  fl->BACH[1] == 0)
		fl->BACHDEP = 1;
}

long long getnumber(FILE *infile)
/*****	reads from stream infile, returns n if a line
      with #n exists before a line with n, nxm or ndm,
      puts the stream to the next line,
      otherwise returns 1 and rewinds the stream to 
      its initial position	*****/
{
	char	string[STRLEN], *str;
	long	pos;
	long long	nmat;

	pos = ftell(infile);
	fgets(string, STRLEN, infile);
	str = strtok(string, "%");
	while (strpbrk(str, "#xd0123456789") == NULL)
	{
		fgets(string, STRLEN, infile);
		str = strtok(string, "%");
	}
	if ((str = strpbrk(str, "#")) != NULL)
	{
		sscanf(str+1, "%lld", &nmat);
		return(nmat);
	}
	else
	{
/* rewind the stream to its initial position */
	  fseek(infile, pos, SEEK_SET);
	  return(1);
	}
}

long long getmat(FILE *infile, long long ***A)
/*****	reads a matrix from stream infile and allocates
      its memory, returns the dimension, 
      the input format is:
      n	full nxn-matrix
      nxm	full nxm-matrix
      nx0	symmetric nxn-matrix given as lower
      triangular matrix
      nx-1	skew-symmetric nxn-matrix given as lower
      triangular matrix (with diagonal)
      nd1	diagonal nxn-matrix given as row
      nd0	scalar nxn-matrix given as number
      the rest of the line containing the format 
      string is ignored	*****/
{
	char	string[STRLEN], *str;
	long long	i, j, dim1, dim2, full, symm, skew, diag, scal;

	full = 0;
	symm = 0;
	skew = 0;
	diag = 0;
	scal = 0;
	fgets(string, STRLEN, infile);
	str = strtok(string, "%");
	while (strpbrk(str, "xd0123456789") == NULL)
/* search for the first string containing a format string */
	{
		fgets(string, STRLEN, infile);
		str = strtok(string, "%");
	}
	if ((str = strpbrk(str, "0123456789")) == NULL)
	{
		fprintf(stderr, "Error: cannot interpret format string %s\n", string);
		exit (2);
	}
	sscanf(str, "%lld", &dim1);
	str = strpbrk(str, "dx");
	if (str == NULL)
/* full quadratic matrix of size dim1 */
		full = 1;
	else if (str[0] == 'x')
	{
		sscanf(str+1, "%lld", &dim2);
		if (dim2 == dim1)
/* full quadratic matrix of size dim1 */
			full = 1;
		else if (dim2 == 0)
/* symmetric matrix given as lower triangular */
			symm = 1;
		else if (dim2 < 0)
/* skew-symmetric matrix given as lower triangular */
			skew = 1;
		else	
		{
			fprintf(stderr, "Error: cannot interpret format string %s\n", str);
			exit (2);
		}
	}
	else if (str[0] == 'd')
	{
	  sscanf(str+1, "%lld", &dim2);
	  if (dim2 == 1)
	    /* diagonal matrix given as vector */
	    diag = 1;
	  else if (dim2 == 0)
	    /* scalar matrix given as number */
	    scal = 1;
	  else	
	    {
	      fprintf(stderr, "Error: can not interpret format string %s\n", str);
	      exit (2);
	    }
	}
	else	
	  {
	    fprintf(stderr, "Error: can not interpret format string %s\n", str);
	    exit (2);
	  }
	if ((*A = (long long**)malloc(dim1 * sizeof(long long*))) == 0)
	  exit (1);
	for (i = 0; i < dim1; ++i)
	{
		if (((*A)[i] = (long long*)malloc(dim1 * sizeof(long long))) == 0)
			exit (1);
	}
	if (full == 1)
	{
		for (i = 0; i < dim1; ++i)
		{
			for (j = 0; j < dim1; ++j)
			  {
			    singleScan(infile, &(*A)[i][j]);
			  }
		}
	}
	else if (symm == 1)
	{
		for (i = 0; i < dim1; ++i)
		{
			for (j = 0; j <= i; ++j)
			{
			  singleScan(infile, &(*A)[i][j]);
			  (*A)[j][i] = (*A)[i][j];
			}
		}
	}
	else if (skew == 1)
	{
	  for (i = 0; i < dim1; ++i)
	    {
	      for (j = 0; j <= i; ++j)
		{
		  singleScan(infile, &(*A)[i][j]);
		  (*A)[j][i] = -(*A)[i][j];
		}
	    }
	}
	else if (diag == 1)
	{
		for (i = 0; i < dim1; ++i)
		{
			for (j = 0; j < dim1; ++j)
				(*A)[i][j] = 0;
		}
		for (i = 0; i < dim1; ++i)
		  {
		    singleScan(infile, &(*A)[i][i]);
		  }
	}
	else if (scal == 1)
	{
		for (i = 0; i < dim1; ++i)
		{
			for (j = 0; j < dim1; ++j)
				(*A)[i][j] = 0;
		}
		singleScan(infile, &(*A)[0][0]);
		for (i = 1; i < dim1; ++i)
			(*A)[i][i] = (*A)[0][0];
	}
	return(dim1);
}

void getvecs(FILE *infile, veclist *V)
/*****	reads vectors from stream infile, first the 
      number then the list, every vector together
      with its norm, 
      returns the number of vectors	*****/
{
	char	string[STRLEN], *str;
	long long	i, j;
	long long *Vvi;
	fgets(string, STRLEN, infile);
	while ((str = strpbrk(string, "0123456789")) == NULL)
/* search for a string containing a number, which is interpreted as the number
   of short vectors */
		fgets(string, STRLEN, infile);
	sscanf(str, "%lld", &V->n);
	if ((V->v = (long long**)malloc((V->n+1) * sizeof(long long*))) == 0)
		exit (1);
	for (i = 0; i <= V->n; ++i)
/* allocate vectors of length V.len */
	{
		if ((V->v[i] = (long long*)malloc(V->len * sizeof(long long))) == 0)
			exit (1);
	}
	for (i = 1; i <= V->n; ++i)
/* read vectors of length V.dim: the entries V.v[i][V.dim],...,V.v[i][V.len-1]
   are for the norms of the vector with respect to the invariant forms */
	{
		Vvi = V->v[i];
		for (j = 0; j <= V->dim; ++j)
		  {
		    singleScan(infile, &Vvi[j]);
		  }
	}
}

void putform(long long **A, long long n)
/*****	prints nxn-matrix A with format conventions as 
      getmat()	*****/
{
	long long	i, j, symm, skew, diag, scal;

	symm = 1;
	skew = 1;
	diag = 1;
	scal = 1;
	for (i = 0; i < n  &&  (symm == 1  ||  skew == 1); ++i)
	{
		for (j = 0; j <= i  &&  (symm == 1  ||  skew == 1); ++j)
		{
			if (A[i][j] != A[j][i])
/* the matrix is not symmetric and hence not diagonal or scalar */
			{
				symm = 0;	
				diag = 0;
				scal = 0;
			}
			if (A[i][j] != -A[j][i])
/* the matrix is not skew-symmetric */
				skew = 0;	
			if (j != i  &&  A[i][j] != 0)
/* the matrix is not diagonal and hence not scalar */
			{
				diag = 0;
				scal = 0;
			}
		}
		if (i > 0  &&  A[i][i] != A[i-1][i-1])
/* the matrix is not scalar */
			scal = 0;
	}
	if (scal == 1)
		printf("%lldd0\n", n);
	else if (diag == 1)
		printf("%lldd1\n", n);
	else if (symm == 1)
		printf("%lldx0\n", n);
	else if (skew == 1)
		printf("%lldx-1\n", n);
	else
		printf("%lld\n", n);
	for (i = 0; i < n; ++i)
	{
		if (scal == 1)
		{
			if (i == 0)
			  {
			    printf(" ");
			    single2Print(stdout, A[0][0]);
			    printf("\n");
			  }
		}
		else if (diag == 1)
		{
		  printf(" ");
		  single2Print(stdout, A[i][i]);
		  if (i == n-1)
		    printf("\n");
		}
		else if (symm == 1  ||  skew == 1)
		{
			for (j = 0; j <= i; ++j)
			  {
			    printf(" ");
			    single2Print(stdout, A[i][j]);
			  }
			printf("\n");
		}
		else
		{
			for (j = 0; j < n; ++j)
			  {
			    printf(" ");
			    single2Print(stdout, A[i][j]);
			  }
			printf("\n");
		}
	}
}

void putgens(group G, flagstruct flags)
/*****	prints the generators of the group, which are
      necessary for generating, if flags.OUTPUT = 0
      the output is in ASCII-format, 
      if flags.OUTPUT = 1 in GAP-format,
      if flags.OUTPUT = 2 in MAGMA-format	*****/
{
	long long	i, j, k, l, dim, ngen, nr;

	dim = G.dim;
	ngen = 0;
	for (i = flags.STAB; i < dim; ++i)
		ngen += G.ng[i] - G.nsg[i];
	if (flags.OUTPUT == 0)
		printf("#%lld\n", ngen);
/* in GAP or MAGMA the matrix group is called BG */
	else if (flags.OUTPUT == 1)
		printf("BG := Group(\n");
	else if (flags.OUTPUT == 2)
		printf("BG := MatrixGroup<%lld, Integers() |\n", dim);
	nr = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
	  for (j = G.nsg[i]; j < G.ng[i]; ++j)
	    {
	      if (flags.OUTPUT == 0)
		printf("%lld\n", dim);
	      else if (flags.OUTPUT >= 1)
		printf("[");
	      for (k = 0; k < dim; ++k)
		{
		  if (flags.OUTPUT == 1)
		    printf("[");
		  for (l = 0; l < dim-1; ++l)
		    {
		      if (flags.OUTPUT == 0)
			{
			  printf(" ");
			  single2Print(stdout, G.g[i][j][k][l]);
			}
		      else if (flags.OUTPUT >= 1)
			single2Print(stdout, G.g[i][j][k][l]);
		    }
		  if (flags.OUTPUT == 0)
		    {
		      printf(" ");
		      single2Print(stdout, G.g[i][j][k][dim-1]);
		      printf("\n");
		    }
		  else if (flags.OUTPUT == 1)
		    {
		      if (k < dim-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf("],\n");
			}
		      else if (k == dim-1  &&  nr < ngen-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf("]],\n");
			}
		      else if (k == dim-1  &&  nr == ngen-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf("]]);\n");
			}
		    }
		  else if (flags.OUTPUT == 2)
		    {
		      if (k < dim-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf(",\n");
			}
		      else if (k == dim-1  &&  nr < ngen-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf("],\n");
			}
		      else if (k == dim-1  &&  nr == ngen-1)
			{
			  single2Print(stdout, G.g[i][j][k][dim-1]);
			  printf("]>;\n");
			}
		    }
		}
	      ++nr;
	    }
	}
}

void putord(group G, flagstruct flags)
/*****	prints orbit lengths of the base vectors,
      which are != 1 and the prime power decomposition
      of G.ord[flags.STAB] *...* G.ord[G.dim-1]
      in ASCII-, GAP- or MAGMA-format	*****/
{
	long long	*ex, i, j, dim, nbase, fac;

	dim = G.dim;
/* nbase is the numbers of standard-base vectors needed for a base in the
   sense of permutation groups (i.e. a set of vectors the stabilizer of
   which is trivial) */
	nbase = 0;
	for (i = flags.STAB; i < dim; ++i)
	{
		if (G.ord[i] > 1)
			++nbase;
	}
	j = 0;
	if (flags.OUTPUT == 0)
		printf("\n|Aut| =");
	else if (flags.OUTPUT >= 1)
		printf("OrbitLengthsBG := [");
/* print the orbit lengths of the points of the base under the stabilizer of the
   preceding base points, the product of these lengths is the group order, for 
   GAP or MAGMA these orbit lengths are assigned to a vector OrbitLengthsBG */
	for (i = flags.STAB; i < dim; ++i)
	{
		if (G.ord[i] > 1)
		{
			if (j < nbase-1)
			{
				if (flags.OUTPUT == 0)
					printf(" %lld *", G.ord[i]);
				else if (flags.OUTPUT >= 1)
					printf(" %lld,", G.ord[i]);
			}
			else
			{
				if (flags.OUTPUT == 0)
					printf(" %lld\n", G.ord[i]);
				else if (flags.OUTPUT >= 1)
					printf(" %lld];\n", G.ord[i]);
			}
			++j;
		}
	}
/* the largest possible prime divisor of the order of an integral nxn-matrix 
   is n+1 */
	if ((ex = (long long*)malloc((dim+2) * sizeof(long long))) == 0)
		exit (1);
/* the order of the group is the product over i^ex[i] */
	for (i = 0; i <= dim+1; ++i)
		ex[i] = 0;
	ex[1] = 1;
	for (i = flags.STAB; i < dim; ++i)
	{
		fac = G.ord[i];
		for (j = 2; j <= dim+1  &&  fac != 1; ++j)
		{
			if (fac % j == 0)
			{
				++ex[j];
				fac = fac / j;
				--j;
			}
		}
	}	
	for (j = dim+1; j > 1  &&  ex[j] == 0; --j);
/* now print the group order in factorized form */
	if (flags.OUTPUT == 0)
		printf("\n=");
	else if (flags.OUTPUT >= 1)
/* in GAP or MAGMA the order of the group is assigned to OrderBG */
		printf("OrderBG :=");
	for (i = 2; i < j; ++i)
	{
		if (ex[i] > 0)
		{
			if (ex[i] == 1)
				printf(" %lld *", i);
			else 
				printf(" %lld^%lld *", i, ex[i]);
		}
	}
	if (ex[j] == 1)
		printf(" %lld", j);
	else 
		printf(" %lld^%lld", j, ex[j]);
	if (flags.OUTPUT == 0)
		printf("\n");
	else if (flags.OUTPUT >= 1)
		printf(";\n");
	free(ex);
}

void putiso(long long **X, flagstruct flags, long long dim)
/*****	prints an isometry in ASCII-, GAP- or
      MAGMA-format	*****/
{
	long long	i, j;

	if (flags.OUTPUT == 0)
		printf("%lld\n", dim);
	else if (flags.OUTPUT >= 1)
/* in GAP or MAGMA the isometry is called Iso */
		printf("Iso := [");
	for (i = 0; i < dim; ++i)
	{
		if (flags.OUTPUT == 1)
			printf("[");
		for (j = 0; j < dim-1; ++j)
		{
			if (flags.OUTPUT == 0)
			  {
			    printf(" ");
			    single2Print(stdout, X[i][j]);
			  }
			else if (flags.OUTPUT >= 1)
			  {
			    single2Print(stdout, X[i][j]);
			    printf(",");
			  }
		}
		if (flags.OUTPUT == 0)
		  {
		    printf(" ");
		    single2Print(stdout, X[i][dim-1]);
		    printf("\n");
		  }
		else if (flags.OUTPUT == 1)
		{
			if (i < dim-1)
			  {
			    single2Print(stdout, X[i][dim-1]);
			    printf("],\n");
			  }
			else if (i == dim-1)
			  {
			    single2Print(stdout, X[i][dim-1]);
			    printf("]];\n");
			  }
		}
		else if (flags.OUTPUT == 2)
		{
			if (i < dim-1)
			  {
			    single2Print(stdout, X[i][dim-1]);
			    printf(",\n");
			  }
			else if (i == dim-1)
			  {
			    single2Print(stdout, X[i][dim-1]);
			    printf("];\n");
			  }
		}
	}
}
