#line 228 "/home/zhiyuan/splash2/codes/null_macros/c.m4.null.POSIX"

#line 1 "lu.C"
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*************************************************************************/
/*                                                                       */
/*  Parallel dense blocked LU factorization (no pivoting)                */
/*                                                                       */
/*  This version contains one dimensional arrays in which the matrix     */
/*  to be factored is stored.                                            */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*  -nN : Decompose NxN matrix.                                          */
/*  -pP : P = number of processors.                                      */
/*  -bB : Use a block size of B. BxB elements should fit in cache for    */
/*        good performance. Small block sizes (B=8, B=16) work well.     */
/*  -s  : Print individual processor timing statistics.                  */
/*  -t  : Test output.                                                   */
/*  -o  : Print out matrix values.                                       */
/*  -h  : Print out command line options.                                */
/*                                                                       */
/*  Note: This version works under both the FORK and SPROC models        */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#line 42
#include <pthread.h>
#line 42
#include <sys/time.h>
#line 42
#include <unistd.h>
#line 42
#include <stdlib.h>
#line 42
#define MAX_THREADS 32
#line 42
pthread_t PThreadTable[MAX_THREADS];
#line 42


#define MAXRAND					32767.0
#define DEFAULT_N				128
#define DEFAULT_P				1
#define DEFAULT_B				16
#define min(a,b) ((a) < (b) ? (a) : (b))
#define PAGE_SIZE				4096

struct GlobalMemory {
  double *t_in_fac;
  double *t_in_solve;
  double *t_in_mod;
  double *t_in_bar;
  double *completion;
  unsigned long starttime;
  unsigned long rf;
  unsigned long rs;
  unsigned long done;
  long id;
  
#line 62
struct {
#line 62
	pthread_mutex_t	mutex;
#line 62
	pthread_cond_t	cv;
#line 62
	unsigned long	counter;
#line 62
	unsigned long	cycle;
#line 62
} (start);
#line 62

  pthread_mutex_t (idlock);
} *Global;

struct LocalCopies {
  double t_in_fac;
  double t_in_solve;
  double t_in_mod;
  double t_in_bar;
};

long n = DEFAULT_N;          /* The size of the matrix */
long P = DEFAULT_P;          /* Number of processors */
long block_size = DEFAULT_B; /* Block dimension */
long nblocks;                /* Number of blocks in each dimension */
long num_rows;               /* Number of processors per row of processor grid */
long num_cols;               /* Number of processors per col of processor grid */
double *a;                   /* a = lu; l and u both placed back in a */
double *rhs;
long *proc_bytes;            /* Bytes to malloc per processor to hold blocks of A*/
long test_result = 0;        /* Test result of factorization? */
long doprint = 0;            /* Print out matrix values? */
long dostats = 0;            /* Print out individual processor statistics? */

void SlaveStart(void);
void OneSolve(long n, long block_size, long MyNum, long dostats);
void lu0(double *a, long n, long stride);
void bdiv(double *a, double *diag, long stride_a, long stride_diag, long dimi, long dimk);
void bmodd(double *a, double *c, long dimi, long dimj, long stride_a, long stride_c);
void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stride);
void daxpy(double *a, double *b, long n, double alpha);
long BlockOwner(long I, long J);
long BlockOwnerColumn(long I, long J);
long BlockOwnerRow(long I, long J);
void lu(long n, long bs, long MyNum, struct LocalCopies *lc, long dostats);
void InitA(double *rhs);
double TouchA(long bs, long MyNum);
void PrintA(void);
void CheckResult(long n, double *a, double *rhs);
void printerr(char *s);

int main(int argc, char *argv[])
{
  long i, ch;
  extern char *optarg;
  double mint, maxt, avgt;
  double min_fac, min_solve, min_mod, min_bar;
  double max_fac, max_solve, max_mod, max_bar;
  double avg_fac, avg_solve, avg_mod, avg_bar;
  unsigned long start;

  {
#line 113
	struct timeval	FullTime;
#line 113

#line 113
	gettimeofday(&FullTime, NULL);
#line 113
	(start) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 113
};

  while ((ch = getopt(argc, argv, "n:p:b:cstoh")) != -1) {
    switch(ch) {
    case 'n': n = atoi(optarg); break;
    case 'p': P = atoi(optarg); break;
    case 'b': block_size = atoi(optarg); break;
    case 's': dostats = 1; break;
    case 't': test_result = !test_result; break;
    case 'o': doprint = !doprint; break;
    case 'h': printf("Usage: LU <options>\n\n");
              printf("options:\n");
              printf("  -nN : Decompose NxN matrix.\n");
              printf("  -pP : P = number of processors.\n");
              printf("  -bB : Use a block size of B. BxB elements should fit in cache for \n");
              printf("        good performance. Small block sizes (B=8, B=16) work well.\n");
              printf("  -c  : Copy non-locally allocated blocks to local memory before use.\n");
              printf("  -s  : Print individual processor timing statistics.\n");
              printf("  -t  : Test output.\n");
              printf("  -o  : Print out matrix values.\n");
              printf("  -h  : Print out command line options.\n\n");
              printf("Default: LU -n%1d -p%1d -b%1d\n",
                     DEFAULT_N,DEFAULT_P,DEFAULT_B);
              exit(0);
              break;
    }
  }

  {;}

  printf("\n");
  printf("Blocked Dense LU Factorization\n");
  printf("     %ld by %ld Matrix\n",n,n);
  printf("     %ld Processors\n",P);
  printf("     %ld by %ld Element Blocks\n",block_size,block_size);
  printf("\n");
  printf("\n");

  num_rows = (long) sqrt((double) P);
  for (;;) {
    num_cols = P/num_rows;
    if (num_rows*num_cols == P)
      break;
    num_rows--;
  }
  nblocks = n/block_size;
  if (block_size * nblocks != n) {
    nblocks++;
  }

  a = (double *) valloc(n*n*sizeof(double));;
  if (a == NULL) {
	  printerr("Could not malloc memory for a.\n");
	  exit(-1);
  }
  rhs = (double *) valloc(n*sizeof(double));;
  if (rhs == NULL) {
	  printerr("Could not malloc memory for rhs.\n");
	  exit(-1);
  }

  Global = (struct GlobalMemory *) valloc(sizeof(struct GlobalMemory));;
  Global->t_in_fac = (double *) valloc(P*sizeof(double));;
  Global->t_in_mod = (double *) valloc(P*sizeof(double));;
  Global->t_in_solve = (double *) valloc(P*sizeof(double));;
  Global->t_in_bar = (double *) valloc(P*sizeof(double));;
  Global->completion = (double *) valloc(P*sizeof(double));;

  if (Global == NULL) {
    printerr("Could not malloc memory for Global\n");
    exit(-1);
  } else if (Global->t_in_fac == NULL) {
    printerr("Could not malloc memory for Global->t_in_fac\n");
    exit(-1);
  } else if (Global->t_in_mod == NULL) {
    printerr("Could not malloc memory for Global->t_in_mod\n");
    exit(-1);
  } else if (Global->t_in_solve == NULL) {
    printerr("Could not malloc memory for Global->t_in_solve\n");
    exit(-1);
  } else if (Global->t_in_bar == NULL) {
    printerr("Could not malloc memory for Global->t_in_bar\n");
    exit(-1);
  } else if (Global->completion == NULL) {
    printerr("Could not malloc memory for Global->completion\n");
    exit(-1);
  }

/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the a
   matrix data across physically distributed memories in a 
   round-robin fashion as desired. */

  {
#line 205
	unsigned long	Error;
#line 205

#line 205
	Error = pthread_mutex_init(&(Global->start).mutex, NULL);
#line 205
	if (Error != 0) {
#line 205
		printf("Error while initializing barrier.\n");
#line 205
		exit(-1);
#line 205
	}
#line 205

#line 205
	Error = pthread_cond_init(&(Global->start).cv, NULL);
#line 205
	if (Error != 0) {
#line 205
		printf("Error while initializing barrier.\n");
#line 205
		pthread_mutex_destroy(&(Global->start).mutex);
#line 205
		exit(-1);
#line 205
	}
#line 205

#line 205
	(Global->start).counter = 0;
#line 205
	(Global->start).cycle = 0;
#line 205
};
  {pthread_mutex_init(&(Global->idlock), NULL);};
  Global->id = 0;

  InitA(rhs);
  if (doprint) {
    printf("Matrix before decomposition:\n");
    PrintA();
  }

  {
#line 215
	long	i, Error;
#line 215

#line 215
	for (i = 0; i < (P) - 1; i++) {
#line 215
		Error = pthread_create(&PThreadTable[i], NULL, (void * (*)(void *))(SlaveStart), NULL);
#line 215
		if (Error != 0) {
#line 215
			printf("Error in pthread_create().\n");
#line 215
			exit(-1);
#line 215
		}
#line 215
	}
#line 215

#line 215
	SlaveStart();
#line 215
};
  {
#line 216
	unsigned long	i, Error;
#line 216
	for (i = 0; i < (P) - 1; i++) {
#line 216
		Error = pthread_join(PThreadTable[i], NULL);
#line 216
		if (Error != 0) {
#line 216
			printf("Error in pthread_join().\n");
#line 216
			exit(-1);
#line 216
		}
#line 216
	}
#line 216
};

  if (doprint) {
    printf("\nMatrix after decomposition:\n");
    PrintA();
  }

  if (dostats) {
    maxt = avgt = mint = Global->completion[0];
    for (i=1; i<P; i++) {
      if (Global->completion[i] > maxt) {
        maxt = Global->completion[i];
      }
      if (Global->completion[i] < mint) {
        mint = Global->completion[i];
      }
      avgt += Global->completion[i];
    }
    avgt = avgt / P;

    min_fac = max_fac = avg_fac = Global->t_in_fac[0];
    min_solve = max_solve = avg_solve = Global->t_in_solve[0];
    min_mod = max_mod = avg_mod = Global->t_in_mod[0];
    min_bar = max_bar = avg_bar = Global->t_in_bar[0];

    for (i=1; i<P; i++) {
      if (Global->t_in_fac[i] > max_fac) {
        max_fac = Global->t_in_fac[i];
      }
      if (Global->t_in_fac[i] < min_fac) {
        min_fac = Global->t_in_fac[i];
      }
      if (Global->t_in_solve[i] > max_solve) {
        max_solve = Global->t_in_solve[i];
      }
      if (Global->t_in_solve[i] < min_solve) {
        min_solve = Global->t_in_solve[i];
      }
      if (Global->t_in_mod[i] > max_mod) {
        max_mod = Global->t_in_mod[i];
      }
      if (Global->t_in_mod[i] < min_mod) {
        min_mod = Global->t_in_mod[i];
      }
      if (Global->t_in_bar[i] > max_bar) {
        max_bar = Global->t_in_bar[i];
      }
      if (Global->t_in_bar[i] < min_bar) {
        min_bar = Global->t_in_bar[i];
      }
      avg_fac += Global->t_in_fac[i];
      avg_solve += Global->t_in_solve[i];
      avg_mod += Global->t_in_mod[i];
      avg_bar += Global->t_in_bar[i];
    }
    avg_fac = avg_fac/P;
    avg_solve = avg_solve/P;
    avg_mod = avg_mod/P;
    avg_bar = avg_bar/P;
  }
  printf("                            PROCESS STATISTICS\n");
  printf("              Total      Diagonal     Perimeter      Interior       Barrier\n");
  printf(" Proc         Time         Time         Time           Time          Time\n");
  printf("    0    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
          Global->completion[0],Global->t_in_fac[0],
          Global->t_in_solve[0],Global->t_in_mod[0],
          Global->t_in_bar[0]);
  if (dostats) {
    for (i=1; i<P; i++) {
      printf("  %3ld    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
              i,Global->completion[i],Global->t_in_fac[i],
              Global->t_in_solve[i],Global->t_in_mod[i],
              Global->t_in_bar[i]);
    }
    printf("  Avg    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           avgt,avg_fac,avg_solve,avg_mod,avg_bar);
    printf("  Min    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           mint,min_fac,min_solve,min_mod,min_bar);
    printf("  Max    %10.0f    %10.0f    %10.0f    %10.0f    %10.0f\n",
           maxt,max_fac,max_solve,max_mod,max_bar);
  }
  printf("\n");
  Global->starttime = start;
  printf("                            TIMING INFORMATION\n");
  printf("Start time                        : %16lu\n", Global->starttime);
  printf("Initialization finish time        : %16lu\n", Global->rs);
  printf("Overall finish time               : %16lu\n", Global->rf);
  printf("Total time with initialization    : %16lu\n", Global->rf-Global->starttime);
  printf("Total time without initialization : %16lu\n", Global->rf-Global->rs);
  printf("\n");

  if (test_result) {
    printf("                             TESTING RESULTS\n");
    CheckResult(n, a, rhs);
  }

  {exit(0);};
}

void SlaveStart()
{
  long MyNum;

  {pthread_mutex_lock(&(Global->idlock));}
    MyNum = Global->id;
    Global->id ++;
  {pthread_mutex_unlock(&(Global->idlock));}

/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration */

  {;};
  OneSolve(n, block_size, MyNum, dostats);
}


void OneSolve(long n, long block_size, long MyNum, long dostats)
{
  unsigned long myrs, myrf, mydone;
  struct LocalCopies *lc;

  lc = (struct LocalCopies *) malloc(sizeof(struct LocalCopies));
  if (lc == NULL) {
    fprintf(stderr,"Proc %ld could not malloc memory for lc\n",MyNum);
    exit(-1);
  }
  lc->t_in_fac = 0.0;
  lc->t_in_solve = 0.0;
  lc->t_in_mod = 0.0;
  lc->t_in_bar = 0.0;

  /* barrier to ensure all initialization is done */
  {
#line 348
	unsigned long	Error, Cycle;
#line 348
	long		Cancel, Temp;
#line 348

#line 348
	Error = pthread_mutex_lock(&(Global->start).mutex);
#line 348
	if (Error != 0) {
#line 348
		printf("Error while trying to get lock in barrier.\n");
#line 348
		exit(-1);
#line 348
	}
#line 348

#line 348
	Cycle = (Global->start).cycle;
#line 348
	if (++(Global->start).counter != (P)) {
#line 348
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 348
		while (Cycle == (Global->start).cycle) {
#line 348
			Error = pthread_cond_wait(&(Global->start).cv, &(Global->start).mutex);
#line 348
			if (Error != 0) {
#line 348
				break;
#line 348
			}
#line 348
		}
#line 348
		pthread_setcancelstate(Cancel, &Temp);
#line 348
	} else {
#line 348
		(Global->start).cycle = !(Global->start).cycle;
#line 348
		(Global->start).counter = 0;
#line 348
		Error = pthread_cond_broadcast(&(Global->start).cv);
#line 348
	}
#line 348
	pthread_mutex_unlock(&(Global->start).mutex);
#line 348
};

  /* to remove cold-start misses, all processors begin by touching a[] */
  TouchA(block_size, MyNum);

  {
#line 353
	unsigned long	Error, Cycle;
#line 353
	long		Cancel, Temp;
#line 353

#line 353
	Error = pthread_mutex_lock(&(Global->start).mutex);
#line 353
	if (Error != 0) {
#line 353
		printf("Error while trying to get lock in barrier.\n");
#line 353
		exit(-1);
#line 353
	}
#line 353

#line 353
	Cycle = (Global->start).cycle;
#line 353
	if (++(Global->start).counter != (P)) {
#line 353
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 353
		while (Cycle == (Global->start).cycle) {
#line 353
			Error = pthread_cond_wait(&(Global->start).cv, &(Global->start).mutex);
#line 353
			if (Error != 0) {
#line 353
				break;
#line 353
			}
#line 353
		}
#line 353
		pthread_setcancelstate(Cancel, &Temp);
#line 353
	} else {
#line 353
		(Global->start).cycle = !(Global->start).cycle;
#line 353
		(Global->start).counter = 0;
#line 353
		Error = pthread_cond_broadcast(&(Global->start).cv);
#line 353
	}
#line 353
	pthread_mutex_unlock(&(Global->start).mutex);
#line 353
};

/* POSSIBLE ENHANCEMENT:  Here is where one might reset the
   statistics that one is measuring about the parallel execution */

  if ((MyNum == 0) || (dostats)) {
    {
#line 359
	struct timeval	FullTime;
#line 359

#line 359
	gettimeofday(&FullTime, NULL);
#line 359
	(myrs) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 359
};
  }

  lu(n, block_size, MyNum, lc, dostats);

  if ((MyNum == 0) || (dostats)) {
    {
#line 365
	struct timeval	FullTime;
#line 365

#line 365
	gettimeofday(&FullTime, NULL);
#line 365
	(mydone) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 365
};
  }

  {
#line 368
	unsigned long	Error, Cycle;
#line 368
	long		Cancel, Temp;
#line 368

#line 368
	Error = pthread_mutex_lock(&(Global->start).mutex);
#line 368
	if (Error != 0) {
#line 368
		printf("Error while trying to get lock in barrier.\n");
#line 368
		exit(-1);
#line 368
	}
#line 368

#line 368
	Cycle = (Global->start).cycle;
#line 368
	if (++(Global->start).counter != (P)) {
#line 368
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 368
		while (Cycle == (Global->start).cycle) {
#line 368
			Error = pthread_cond_wait(&(Global->start).cv, &(Global->start).mutex);
#line 368
			if (Error != 0) {
#line 368
				break;
#line 368
			}
#line 368
		}
#line 368
		pthread_setcancelstate(Cancel, &Temp);
#line 368
	} else {
#line 368
		(Global->start).cycle = !(Global->start).cycle;
#line 368
		(Global->start).counter = 0;
#line 368
		Error = pthread_cond_broadcast(&(Global->start).cv);
#line 368
	}
#line 368
	pthread_mutex_unlock(&(Global->start).mutex);
#line 368
};

  if ((MyNum == 0) || (dostats)) {
    {
#line 371
	struct timeval	FullTime;
#line 371

#line 371
	gettimeofday(&FullTime, NULL);
#line 371
	(myrf) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 371
};
    Global->t_in_fac[MyNum] = lc->t_in_fac;
    Global->t_in_solve[MyNum] = lc->t_in_solve;
    Global->t_in_mod[MyNum] = lc->t_in_mod;
    Global->t_in_bar[MyNum] = lc->t_in_bar;
    Global->completion[MyNum] = mydone-myrs;
  }
  if (MyNum == 0) {
    Global->rs = myrs;
    Global->done = mydone;
    Global->rf = myrf;
  }
}


void lu0(double *a, long n, long stride)
{
  long j, k, length;
  double alpha;

  for (k=0; k<n; k++) {
    /* modify subsequent columns */
    for (j=k+1; j<n; j++) {
      a[k+j*stride] /= a[k+k*stride];
      alpha = -a[k+j*stride];
      length = n-k-1;
      daxpy(&a[k+1+j*stride], &a[k+1+k*stride], n-k-1, alpha);
    }
  }
}


void bdiv(double *a, double *diag, long stride_a, long stride_diag, long dimi, long dimk)
{
  long j, k;
  double alpha;

  for (k=0; k<dimk; k++) {
    for (j=k+1; j<dimk; j++) {
      alpha = -diag[k+j*stride_diag];
      daxpy(&a[j*stride_a], &a[k*stride_a], dimi, alpha);
    }
  }
}


void bmodd(double *a, double *c, long dimi, long dimj, long stride_a, long stride_c)
{
  long j, k, length;
  double alpha;

  for (k=0; k<dimi; k++)
    for (j=0; j<dimj; j++) {
      c[k+j*stride_c] /= a[k+k*stride_a];
      alpha = -c[k+j*stride_c];
      length = dimi - k - 1;
      daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
    }
}


void bmod(double *a, double *b, double *c, long dimi, long dimj, long dimk, long stride)
{
  long j, k;
  double alpha;

  for (k=0; k<dimk; k++) {
    for (j=0; j<dimj; j++) {
      alpha = -b[k+j*stride];
      daxpy(&c[j*stride], &a[k*stride], dimi, alpha);
    }
  }
}


void daxpy(double *a, double *b, long n, double alpha)
{
  long i;

  for (i=0; i<n; i++) {
    a[i] += alpha*b[i];
  }
}


long BlockOwner(long I, long J)
{
//	return((I%num_cols) + (J%num_rows)*num_cols);
	return((I + J) % P);
}

long BlockOwnerColumn(long I, long J)
{
	return(I % P);
}

long BlockOwnerRow(long I, long J)
{
	return(((J % P) + (P / 2)) % P);
}

void lu(long n, long bs, long MyNum, struct LocalCopies *lc, long dostats)
{
  long i, il, j, jl, k, kl, I, J, K;
  double *A, *B, *C, *D;
  long strI;
  unsigned long t1, t2, t3, t4, t11, t22;

  strI = n;
  for (k=0, K=0; k<n; k+=bs, K++) {
    kl = k+bs; 
    if (kl>n) {
      kl = n;
    }

    if ((MyNum == 0) || (dostats)) {
      {
#line 487
	struct timeval	FullTime;
#line 487

#line 487
	gettimeofday(&FullTime, NULL);
#line 487
	(t1) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 487
};
    }

    /* factor diagonal block */
    if (BlockOwner(K, K) == MyNum) {
      A = &(a[k+k*n]); 
      lu0(A, kl-k, strI);
    }

    if ((MyNum == 0) || (dostats)) {
      {
#line 497
	struct timeval	FullTime;
#line 497

#line 497
	gettimeofday(&FullTime, NULL);
#line 497
	(t11) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 497
};
    }

    {
#line 500
	unsigned long	Error, Cycle;
#line 500
	long		Cancel, Temp;
#line 500

#line 500
	Error = pthread_mutex_lock(&(Global->start).mutex);
#line 500
	if (Error != 0) {
#line 500
		printf("Error while trying to get lock in barrier.\n");
#line 500
		exit(-1);
#line 500
	}
#line 500

#line 500
	Cycle = (Global->start).cycle;
#line 500
	if (++(Global->start).counter != (P)) {
#line 500
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 500
		while (Cycle == (Global->start).cycle) {
#line 500
			Error = pthread_cond_wait(&(Global->start).cv, &(Global->start).mutex);
#line 500
			if (Error != 0) {
#line 500
				break;
#line 500
			}
#line 500
		}
#line 500
		pthread_setcancelstate(Cancel, &Temp);
#line 500
	} else {
#line 500
		(Global->start).cycle = !(Global->start).cycle;
#line 500
		(Global->start).counter = 0;
#line 500
		Error = pthread_cond_broadcast(&(Global->start).cv);
#line 500
	}
#line 500
	pthread_mutex_unlock(&(Global->start).mutex);
#line 500
};

    if ((MyNum == 0) || (dostats)) {
      {
#line 503
	struct timeval	FullTime;
#line 503

#line 503
	gettimeofday(&FullTime, NULL);
#line 503
	(t2) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 503
};
    }

    /* divide column k by diagonal block */
    D = &(a[k+k*n]);
    for (i=kl, I=K+1; i<n; i+=bs, I++) {
      if (BlockOwner/*Column*/(I, K) == MyNum) {  /* parcel out blocks */
	      /*if (K == 0) printf("C%lx\n", BlockOwnerColumn(I, K));*/
        il = i + bs;
        if (il > n) {
          il = n;
        }
        A = &(a[i+k*n]);
        bdiv(A, D, strI, n, il-i, kl-k);
      }
    }
    /* modify row k by diagonal block */
    for (j=kl, J=K+1; j<n; j+=bs, J++) {
      if (BlockOwner/*Row*/(K, J) == MyNum) {  /* parcel out blocks */
	      /*if (K == 0) printf("R%lx\n", BlockOwnerRow(K, J));*/
        jl = j+bs;
        if (jl > n) {
          jl = n;
        }
        A = &(a[k+j*n]);
        bmodd(D, A, kl-k, jl-j, n, strI);
      }
    }

    if ((MyNum == 0) || (dostats)) {
      {
#line 533
	struct timeval	FullTime;
#line 533

#line 533
	gettimeofday(&FullTime, NULL);
#line 533
	(t22) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 533
};
    }

    {
#line 536
	unsigned long	Error, Cycle;
#line 536
	long		Cancel, Temp;
#line 536

#line 536
	Error = pthread_mutex_lock(&(Global->start).mutex);
#line 536
	if (Error != 0) {
#line 536
		printf("Error while trying to get lock in barrier.\n");
#line 536
		exit(-1);
#line 536
	}
#line 536

#line 536
	Cycle = (Global->start).cycle;
#line 536
	if (++(Global->start).counter != (P)) {
#line 536
		pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &Cancel);
#line 536
		while (Cycle == (Global->start).cycle) {
#line 536
			Error = pthread_cond_wait(&(Global->start).cv, &(Global->start).mutex);
#line 536
			if (Error != 0) {
#line 536
				break;
#line 536
			}
#line 536
		}
#line 536
		pthread_setcancelstate(Cancel, &Temp);
#line 536
	} else {
#line 536
		(Global->start).cycle = !(Global->start).cycle;
#line 536
		(Global->start).counter = 0;
#line 536
		Error = pthread_cond_broadcast(&(Global->start).cv);
#line 536
	}
#line 536
	pthread_mutex_unlock(&(Global->start).mutex);
#line 536
};

    if ((MyNum == 0) || (dostats)) {
      {
#line 539
	struct timeval	FullTime;
#line 539

#line 539
	gettimeofday(&FullTime, NULL);
#line 539
	(t3) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 539
};
    }

    /* modify subsequent block columns */
    for (i=kl, I=K+1; i<n; i+=bs, I++) {
      il = i+bs;
      if (il > n) {
        il = n;
      }
      A = &(a[i+k*n]);
      for (j=kl, J=K+1; j<n; j+=bs, J++) {
        jl = j + bs;
        if (jl > n) {
          jl = n;
        }
        if (BlockOwner(I, J) == MyNum) {  /* parcel out blocks */
//		if (K == 0) printf("%lx\n", BlockOwner(I, J));
          B = &(a[k+j*n]);
          C = &(a[i+j*n]);
          bmod(A, B, C, il-i, jl-j, kl-k, n);
        }
      }
    }
    if ((MyNum == 0) || (dostats)) {
      {
#line 563
	struct timeval	FullTime;
#line 563

#line 563
	gettimeofday(&FullTime, NULL);
#line 563
	(t4) = (unsigned long)(FullTime.tv_usec + FullTime.tv_sec * 1000000);
#line 563
};
      lc->t_in_fac += (t11-t1);
      lc->t_in_solve += (t22-t2);
      lc->t_in_mod += (t4-t3);
      lc->t_in_bar += (t2-t11) + (t3-t22);
    }
  }
}


void InitA(double *rhs)
{
  long i, j;

  srand48((long) 1);
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      a[i+j*n] = (double) lrand48()/MAXRAND;
      if (i == j) {
	a[i+j*n] *= 10;
      }
    }
  }

  for (j=0; j<n; j++) {
    rhs[j] = 0.0;
  }
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      rhs[i] += a[i+j*n];
    }
  }
}


double TouchA(long bs, long MyNum)
{
  long i, j, I, J;
  double tot = 0.0;

  for (J=0; J*bs<n; J++) {
    for (I=0; I*bs<n; I++) {
      if (BlockOwner(I, J) == MyNum) {
        for (j=J*bs; j<(J+1)*bs && j<n; j++) {
          for (i=I*bs; i<(I+1)*bs && i<n; i++) {
            tot += a[i+j*n];
          }
        }
      }
    }
  }
  return(tot);
}


void PrintA()
{
  long i, j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      printf("%8.1f ", a[i+j*n]);
    }
    printf("\n");
  }
}


void CheckResult(long n, double *a, double *rhs)
{
  long i, j, bogus = 0;
  double *y, diff, max_diff;

  y = (double *) malloc(n*sizeof(double));
  if (y == NULL) {
    printerr("Could not malloc memory for y\n");
    exit(-1);
  }
  for (j=0; j<n; j++) {
    y[j] = rhs[j];
  }
  for (j=0; j<n; j++) {
    y[j] = y[j]/a[j+j*n];
    for (i=j+1; i<n; i++) {
      y[i] -= a[i+j*n]*y[j];
    }
  }

  for (j=n-1; j>=0; j--) {
    for (i=0; i<j; i++) {
      y[i] -= a[i+j*n]*y[j];
    }
  }

  max_diff = 0.0;
  for (j=0; j<n; j++) {
    diff = y[j] - 1.0;
    if (fabs(diff) > 0.00001) {
      bogus = 1;
      max_diff = diff;
    }
  }
  if (bogus) {
    printf("TEST FAILED: (%.5f diff)\n", max_diff);
  } else {
    printf("TEST PASSED\n");
  }
  free(y);
}


void printerr(char *s)
{
  fprintf(stderr,"ERROR: %s\n",s);
}

