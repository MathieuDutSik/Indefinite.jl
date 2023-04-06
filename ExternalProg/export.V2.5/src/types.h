/* some definitions of constants */


/* the maximal entry in the short vectors should not exceed MAXENTRY to
   avoid overflow */
/* it is checked that the product of the maximal entry in the short vectors
   and the maximal entry in the products of Gram-matrices with short vectors
   does not exceed MAXNORM, since the program works just with int's */
#define	STRLEN		80
/* since the LLL-reduction works with a floating-point model, there might be
   some rounding error, which should not reach EPS */
#define	EPS		0.001
/* the constant for the LLL-condition, which should be < 1, but 0.75 is the
   standard value, a higher value might yield a slightly better result but
   increases the running time */
#define	LLL_CONST	0.75
#define	DEF_PRIME	16001

#ifndef SEEK_SET		/* some UNIX-systems do not define SEEK_SET */
#	define SEEK_SET	0	/* Set file pointer to "offset" */
#endif

/* structure to hold the generators of the group according to the stabilizer
   chain */
typedef struct {
	long long	dim;
	long long	*ord;
	long long	*ng;
	long long	*nsg;
	long long	****g;
	} group;

/* structure for a list of vectors */
typedef struct {
  long long	dim;
  long long	len;
  long long	prime;
  long long	n;
  long long	**v;
} veclist;

/* structure for the invariant forms and their products with the list of
   short vectors */
typedef struct {
	long long	dim;
	long long	n;
	long long	***A;
	long long	***v;
	} invar;

/* structure for the diagonal of the fingerprint, the new order and the indices
   of the standard-base in the list of short vectors */
typedef struct {
	long long	*diag;
	long long	*per;
	long long	*e;
	} fpstruct;

/* structure for the scalar product combinations, the transformation matrices
   between the vector sums and that lattice bases and the scalar products of
   the lattice bases */
typedef struct {
	veclist	list;
	long long	rank;
	long long	**trans;
	long long	**coef;
	long long	***F;
	} scpcomb;

/* structure for a Bacher-polynomial */
typedef struct {
	long long	mind;
	long long	maxd;
	long long	sum;
	long long	*coef;
	} bachpol;

/* structure for the option flags */
typedef struct {
	long long	DEPTH;
	long long	STAB;
	long long	BACH[3];
	long long	BACHDEP;
	long long	BACHSCP;
	long long	GEN;
	long long	PRINT;
	long long	VEC;
	long long	OUTPUT;
	} flagstruct;

/* functions in auttools.c */
long long cand(long long *CI, long long I, long long *x, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags);
void autom(group *G, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags);
long long aut(long long step, long long *x, long long **C, group *G, veclist V, invar F, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags);

/* functions in bachtools.c */
void bacher(bachpol *pol, long long I, long long S, veclist	V, long long **Fv);
long long bachcomp(bachpol pol, long long I, long long S, veclist V, long long **Fv);
void fputbach(FILE *outfile, bachpol pol);

/* functions in iotools.c */
void single2Print(FILE *f, long long TheNb);
void getflags(flagstruct *fl, int argc, char *argv[]);
long long getnumber(FILE *infile);
long long getmat(FILE *infile, long long ***A);
void getvecs(FILE *infile, veclist *V);
void putform(long long **A, long long n);
void putgens(group G, flagstruct flags);
long long GetMAXNORM();
long long GetMAXENTRY();
void putord(group G, flagstruct flags);
void putiso(long long **X, flagstruct flags, long long dim);

/* functions in isotools.c */
long long isocand(long long *CI, long long I, long long *x, veclist V, invar F, invar FF, fpstruct fp, scpcomb *comb, bachpol *bach, flagstruct flags);
void isometry(veclist V, invar F, invar FF, fpstruct fp, long long ***G, long long nG, scpcomb *comb, bachpol *bach, flagstruct flags);
long long iso(long long step, long long *x, long long **C, veclist V, invar F, invar FF, fpstruct fp, long long ***G, long long nG, scpcomb *comb, bachpol *bach, flagstruct flags);
long long isostab(long long ****H, long long pt, long long ***G, long long nG, veclist V, long long Maxfail);

/* functions in lattools.c */
long long lll(long long **F, long long **T, long long dim);
void initialize(long long i, double **mu, double *B, long long **gram);
long long red(long long *m, long long *l, long long **gram, long long **T, double **mu, double *B, long long dim, long long *rank);
void check(double *B, long long *m, long long *l, long long **gram, long long **T, double **mu, long long dim, long long *rank);
void decrease(long long *m, long long *l, long long **gram, double **mu, double *B, long long *rank);
void interchange(long long *m, long long *l, long long **gram, long long **T, double **mu, double *B, long long dim);
long long MyRound(double x);
void shortvec(long long **F, long long max, veclist *V);
void shrt(long long k, double norm, long long max, veclist *V, long long *vec, double **mu, long long **T, double *B, long long *con);
void change(long long **v, long long nv, long long **T, long long n);

/* functions in mattools.c */
void vecmatmul(long long *x, long long **A, long long n, long long *y);
void matmul(long long **A, long long **B, long long n, long long **C);
long long scp(long long *x, long long **F, long long *y, long long n);
long long sscp(long long *x, long long *y, long long n);
void psolve(long long **X, long long **A, long long **B, long long n, long long p);
void pgauss(long long r, long long **A, long long **B, long long n, long long p);
long long isprime(long long n);

/* functions in orbtools.c */
long long operate(long long nr, long long **A, veclist V);
long long orbit(long long	*pt, long long npt, long long ***G, long long nG, veclist V, long long **orb);
long long orbitlen(long long pt, long long orblen, long long ***G, long long nG, veclist V);
long long delete(long long *orb1, long long l1, long long *orb2, long long l2);
void stab(long long I, group *G, fpstruct fp, veclist V);
void matgen(long long *x, long long **X, long long dim, long long *per, long long **v);
void stabil(long long **S, long long *x1, long long *x2, long long *per, long long **G, veclist V);

/* functions in preproc.c */
void    checkvecs(veclist *V, invar F, veclist norm);
long long     checkgen(long long **g, invar F);
void    fingerprint(fpstruct *fp, invar F, veclist V);
long long     possible(invar F, long long ***Ftr, veclist V, long long *per, long long I, long long J);
void    scpvector(long long *scpvec, long long *w, long long *b, long long I, long long dep, invar F);
void    scpvecs(veclist *list, long long ***vec, long long I, long long *b, long long dep, veclist V, invar F);
void    base(scpcomb *com, long long ***b, long long **v, long long **F, long long dim);
void    coef(scpcomb *com, long long **b, long long **v, long long **F, long long dim);
void    scpforms(scpcomb *com, long long **b, invar F);

/* functions in sorttools.c */
long long comp(long long *x, long long *y, long long n);
long long numberof(long long *vec, veclist V);
void sortvecs(veclist *V);
void quicksort(long long **v, long long inf, long long sup, long long dim);
