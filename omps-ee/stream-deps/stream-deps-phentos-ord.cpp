/*-----------------------------------------------------------------------*/
/* Program: Stream                                                       */
/* Adapted to StarSs by Rosa M. Badia (Barcelona Supercomputing Center)  */
/* Revision: $Id: stream.c,v 5.8 2007/02/19 23:57:39 mccalpin Exp mccalpin $ */
/* Original code developed by John D. McCalpin                           */
/* Programmers: John D. McCalpin                                         */
/*              Joe R. Zagar                                             */
/*                                                                       */
/* This program measures memory transfer rates in MB/s for simple        */
/* computational kernels coded in C.                                     */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2005: John D. McCalpin                                 */
/*-----------------------------------------------------------------------*/
/* License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*         "tuned STREAM benchmark results"                              */
/*         "based on a variant of the STREAM benchmark code"             */
/*         Other comparable, clear and reasonable labelling is           */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.                   */
/*-----------------------------------------------------------------------*/
# include <stdio.h>
# include <math.h>
# include <float.h>
# include "limits.h"
# include <sys/time.h>
# include <stdlib.h>
# include "omp.h"

#include "../../../phentos/phentos/femtos.hpp"

/* INSTRUCTIONS:
 *
 *	1) Stream requires a good bit of memory to run.  Adjust the
 *          value of 'N' (below) to give a 'timing calibration' of 
 *          at least 20 clock-ticks.  This will provide rate estimates
 *          that should be good to about 5% precision.
 */

// # define N	128*1024*1024 
# ifndef N
# define N	16*1024*1024 
# endif
# define NTIMES	3
# define OFFSET	0

/*
 *	3) Compile the code with full optimization.  Many compilers
 *	   generate unreasonably bad code before the optimizer tightens
 *	   things up.  If the results are unreasonably good, on the
 *	   other hand, the optimizer might be too smart for me!
 *
 *         Try compiling with:
 *               cc -O stream_omp.c -o stream_omp
 *
 *         This is known to work on Cray, SGI, IBM, and Sun machines.
 *
 *
 *	4) Mail the results to mccalpin@cs.virginia.edu
 *	   Be sure to include:
 *		a) computer hardware model number and software revision
 *		b) the compiler flags
 *		c) all of the output from the test case.
 * Thanks!
 *
 */

# define HLINE "-------------------------------------------------------------\n"

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif

#define OmpSs
//#define SMPSs
//#define CellSs
//#define CellSs_tracing
#define TUNED
// Definition of BSIZE for CellSs
#ifdef CellSs
#define BSIZE 4000
#endif
// Definition of BSIZE for CellSs when using tracing
#ifdef CellSs_tracing
#define BSIZE 3200
#endif
#ifdef SMPSs
// Definition of BSIZE for SMPSs
#define BSIZE N/64
#endif
#ifdef OmpSs
// Definitions of BSIZE for OmpSs
#define BSIZE N/64
#endif

static double	_a[N+OFFSET],
		_b[N+OFFSET],
		_c[N+OFFSET];

static double	avgtime[4] = {0}, maxtime[4] = {0},
		mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};

static char	*label[4] = {"Copy:      ", "Scale:     ",
    "Add:       ", "Triad:     "};

static double	bytes[4] = {
    2 * sizeof(double) * N,
    2 * sizeof(double) * N,
    3 * sizeof(double) * N,
    3 * sizeof(double) * N
    };

extern double mysecond();
extern void checkSTREAMresults();
#ifdef TUNED
extern void tuned_STREAM_Copy();
extern void tuned_STREAM_Scale(double scalar);
extern void tuned_STREAM_Add();
extern void tuned_STREAM_Triad(double scalar);
#endif

double scalar;

uint64_t num_iterations = 0;
unsigned int numPendingWorkRequests = 0;

// #pragma omp task out ([bs]a, [bs]b, [bs]c)
void init_task(double *a, double *b, double *c, int bs)
{
	int j;	
	for (j=0; j < bs; j++){
	        a[j] = 1.0;
	        b[j] = 2.0;
	        c[j] = 0.0;
		a[j] = 2.0E0 * a[j];
  	}
}

void ORD_init_task(uint64_t swID)
{
	uint64_t i = (uint64_t) metadataArray[swID].depAddresses0[0];
	double *a = &(_a[i]);
	double *b = &(_b[i]);
	double *c = &(_c[i]);

	int j;	
	for (j=0; j < BSIZE; j++){
	        a[j] = 1.0;
	        b[j] = 2.0;
	        c[j] = 0.0;
		a[j] = 2.0E0 * a[j];
  	}
}

void tuned_initialization()
{
	int j;
	uint64_t swID = 0;
        for (j=0; j<N; j+=BSIZE) {
            // init_task (&a[j], &b[j], &c[j], BSIZE); 
			swID = getNewSWID(swID);
			metadataArray[swID].functionAddr = (unsigned long long) ORD_init_task;
			metadataArray[swID].depAddresses0[0] = (unsigned long long) j;

			asm volatile ("fence" ::: "memory");
			num_iterations++;

#ifdef ZERO_PACKETS_V2
			make_submission_request_or_work(48, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
			make_submission_request_or_work(12, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 3, numPendingWorkRequests);
#endif

			submit_three_or_work((unsigned long long) &(_a[j]), 2, numPendingWorkRequests);
			submit_three_or_work((unsigned long long) &(_b[j]), 2, numPendingWorkRequests);
			submit_three_or_work((unsigned long long) &(_c[j]), 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
#endif
		}
}

int
main(int argc, char *argv[])
    {
    int			quantum, checktick();
    int			BytesPerWord;
    register int	j, k;
    double		t, times[4][NTIMES], total_time, total_bytes;



    /* --- SETUP --- determine precision and check timing --- */

    printf(HLINE);
    printf("STREAM version $Revision: 5.8 $\n");
    printf(HLINE);
    BytesPerWord = sizeof(double);
    printf("This system uses %d bytes per DOUBLE PRECISION word.\n",
	BytesPerWord);

    printf(HLINE);
    printf("Array size = %d, Offset = %d\n" , N, OFFSET);
    printf("Total memory required = %.1f MB.\n",
	(3.0 * BytesPerWord) * ( (double) N / 1048576.0));
    printf("Each test is run %d times, but only\n", NTIMES);
    printf("the *best* time for each is used.\n");

#ifdef OmpSs
    printf(HLINE);
    k = omp_get_num_threads();
    printf ("Number of Threads = %i\n",k);
#endif    
#ifdef SMPSs
    printf(HLINE);
    k = atoi (getenv ("CSS_NUM_CPUS"));
    printf ("Number of CSS Threads = %i\n",k);
#endif
#ifdef CellSs
    printf(HLINE);
    k = atoi (getenv ("CSS_NUM_SPUS"));
    printf ("Number of CSS Threads = %i\n",k);
#endif
#ifdef CellSs_tracing
    printf(HLINE);
    k = atoi (getenv ("CSS_NUM_SPUS"));
    printf ("Number of CSS Threads = %i\n",k);
#endif
#ifdef PHENTOS
	femtos_init();
  #ifdef ZERO_PACKETS_V2
		#pragma message ("Benchmark compiled with support for ZERO_PACKETS_V2")
		#pragma message ("Benchmark compiled with support for SW-BASED PADDING")
		#pragma message ("Benchmark compiled with support for RESET COUNT OF SUBMITTED TASKS")
  #endif
#endif

    printf(HLINE);

    printf ("Printing one line per active thread....\n");

    /* Get initial value for system clock. */

total_time = mysecond();
    tuned_initialization();

    /*	--- MAIN LOOP --- repeat test cases NTIMES times --- */

    scalar = 3.0;
    for (k=0; k<NTIMES; k++)
	{
	times[0][k] = mysecond();
#ifdef TUNED
        tuned_STREAM_Copy();
#else
printf("WARNING: This version is a port to SMPSs that only works for TUNED option \n");
	for (j=0; j<N; j++)
	    c[j] = a[j];
#endif
	times[0][k] = mysecond() - times[0][k];
	
	times[1][k] = mysecond();
#ifdef TUNED
        tuned_STREAM_Scale(scalar);
#else
	for (j=0; j<N; j++)
	    b[j] = scalar*c[j];
#endif
	times[1][k] = mysecond() - times[1][k];
	
	times[2][k] = mysecond();
#ifdef TUNED
        tuned_STREAM_Add();
#else
	for (j=0; j<N; j++)
	    c[j] = a[j]+b[j];
#endif
	times[2][k] = mysecond() - times[2][k];
	
	times[3][k] = mysecond();
#ifdef TUNED
        tuned_STREAM_Triad(scalar);
#else
	for (j=0; j<N; j++)
	    a[j] = b[j]+scalar*c[j];
#endif
	times[3][k] = mysecond() - times[3][k];
	}

// #pragma omp taskwait 
task_wait_and_try_executing_tasks(num_iterations);

total_time = mysecond() - total_time;
    /*	--- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
	{
	for (j=0; j<4; j++)
	    {
	    avgtime[j] = avgtime[j] + times[j][k];
	    mintime[j] = MIN(mintime[j], times[j][k]);
	    maxtime[j] = MAX(maxtime[j], times[j][k]);
	    }
	}
    

   total_bytes = bytes[0] + bytes[1] + bytes [2] + bytes [3];
    printf ("Average Rate (MB/s): %11.4f \n", 1.0E-06 * total_bytes*NTIMES/total_time);
    printf ("note: in this version, the average rate per function\n");
    printf ("can not be provided, use tracing to check it\n");

    printf(HLINE);

    printf("TOTAL time (including initialization) =  %11.4f seconds\n", total_time);
    /* --- Check Results --- */
    checkSTREAMresults();
    printf(HLINE);

    return 0;
}

# define	M	20

int
checktick()
    {
    int		i, minDelta, Delta;
    double	t1, t2, timesfound[M];

/*  Collect a sequence of M unique time values from the system. */

    for (i = 0; i < M; i++) {
	t1 = mysecond();
	while( ((t2=mysecond()) - t1) < 1.0E-6 )
	    ;
	timesfound[i] = t1 = t2;
	}

/*
 * Determine the minimum difference between these M values.
 * This result will be our estimate (in microseconds) for the
 * clock granularity.
 */

    minDelta = 1000000;
    for (i = 1; i < M; i++) {
	Delta = (int)( 1.0E6 * (timesfound[i]-timesfound[i-1]));
	minDelta = MIN(minDelta, MAX(Delta,0));
	}

   return(minDelta);
    }



/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond()
{
        struct timeval tp;
        int i;

        i = gettimeofday(&tp,NULL);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void checkSTREAMresults ()
{
	double aj,bj,cj,scalar;
	double asum,bsum,csum;
	double epsilon;
	int	j,k;

    /* reproduce initialization */
	aj = 1.0;
	bj = 2.0;
	cj = 0.0;
    /* a[] is modified during timing check */
	aj = 2.0E0 * aj;
    /* now execute timing loop */
	scalar = 3.0;
	for (k=0; k<NTIMES; k++)
        {
            cj = aj;
            bj = scalar*cj;
            cj = aj+bj;
            aj = bj+scalar*cj;
        }
	aj = aj * (double) (N);
	bj = bj * (double) (N);
	cj = cj * (double) (N);

	asum = 0.0;
	bsum = 0.0;
	csum = 0.0;
	for (j=0; j<N; j++) {
		asum += _a[j];
		bsum += _b[j];
		csum += _c[j];
	}
#ifdef VERBOSE
	printf ("Results Comparison: \n");
	printf ("        Expected  : %f %f %f \n",aj,bj,cj);
	printf ("        Observed  : %f %f %f \n",asum,bsum,csum);
#endif

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
	epsilon = 1.e-8;

	if (abs(aj-asum)/asum > epsilon) {
		printf ("Failed Validation on array a[]\n");
		printf ("        Expected  : %f \n",aj);
		printf ("        Observed  : %f \n",asum);
	}
	else if (abs(bj-bsum)/bsum > epsilon) {
		printf ("Failed Validation on array b[]\n");
		printf ("        Expected  : %f \n",bj);
		printf ("        Observed  : %f \n",bsum);
	}
	else if (abs(cj-csum)/csum > epsilon) {
		printf ("Failed Validation on array c[]\n");
		printf ("        Expected  : %f \n",cj);
		printf ("        Observed  : %f \n",csum);
	}
	else {
		printf ("Solution Validates\n");
	}
}
// #pragma omp task in ([bs]a) out ([bs]c)
void copy_task(double *a, double *c, int bs)
{
	int j;	
	for (j=0; j < bs; j++)
		c[j] = a[j];
}

void ORD_copy_task(uint64_t swID)
{
	uint64_t i = (uint64_t) metadataArray[swID].depAddresses0[0];
	double *a = &(_a[i]);
	double *c = &(_c[i]);

	int j;	
	for (j=0; j < BSIZE; j++)
		c[j] = a[j];
}

void tuned_STREAM_Copy()
{
	uint64_t swID = 0;
	int j;

        for (j=0; j<N; j+=BSIZE) {
            // copy_task (&a[j], &c[j], BSIZE); 
			swID = getNewSWID(swID);
			metadataArray[swID].functionAddr = (unsigned long long) ORD_copy_task;
			metadataArray[swID].depAddresses0[0] = (unsigned long long) j;

			asm volatile ("fence" ::: "memory");
			num_iterations++;

#ifdef ZERO_PACKETS_V2
			make_submission_request_or_work(48, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
			make_submission_request_or_work(9, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 2, numPendingWorkRequests);
#endif

			submit_three_or_work((unsigned long long) &(_a[j]), 0, numPendingWorkRequests);
			submit_three_or_work((unsigned long long) &(_c[j]), 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
#endif
		}
}

// #pragma omp task in ([bs]c) out ([bs]b)
void scale_task (double *b, double *c, double scalar, int bs)
{
	int j;	
	for (j=0; j < bs; j++)
	    b[j] = scalar*c[j];
}

void ORD_scale_task (uint64_t swID)
{
	uint64_t i = (uint64_t) metadataArray[swID].depAddresses0[0];
	double *b = &(_b[i]);
	double *c = &(_c[i]);

	int j;	
	for (j=0; j < BSIZE; j++)
	    b[j] = scalar*c[j];
}

void tuned_STREAM_Scale(double scalar)
{
	int j;
	uint64_t swID = 0;

	for (j=0; j<N; j+=BSIZE) {
			// scale_task (&b[j], &c[j], scalar, BSIZE); 
			swID = getNewSWID(swID);
			metadataArray[swID].functionAddr = (unsigned long long) ORD_scale_task;
			metadataArray[swID].depAddresses0[0] = (unsigned long long) j;

			asm volatile ("fence" ::: "memory");
			num_iterations++;

#ifdef ZERO_PACKETS_V2
			make_submission_request_or_work(48, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
			make_submission_request_or_work(9, 0, numPendingWorkRequests);
			submit_three_or_work(swID, 2, numPendingWorkRequests);
#endif

			submit_three_or_work((unsigned long long) &(_b[j]), 2, numPendingWorkRequests);
			submit_three_or_work((unsigned long long) &(_c[j]), 0, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
			submit_three_or_work(0, 0, numPendingWorkRequests);
#endif
	}
}

// #pragma omp task in ([bs]a, [bs]b) out ([bs]c)
void add_task (double *a, double *b, double *c, int bs)
{
	int j;	
	for (j=0; j < bs; j++)
	   c[j] = a[j]+b[j]; 
}

void ORD_add_task (uint64_t swID)
{
	uint64_t i = (uint64_t) metadataArray[swID].depAddresses0[0];
	double *a = &(_a[i]);
	double *b = &(_b[i]);
	double *c = &(_c[i]);

	int j;	
	for (j=0; j < BSIZE; j++)
	   c[j] = a[j]+b[j]; 
}

void tuned_STREAM_Add()
{
	int j;
	uint64_t swID = 0;
	for (j=0; j<N; j+=BSIZE) {
		// add_task(&a[j], &b[j], &c[j], BSIZE); 
		swID = getNewSWID(swID);
		metadataArray[swID].functionAddr = (unsigned long long) ORD_add_task;
		metadataArray[swID].depAddresses0[0] = (unsigned long long) j;

		asm volatile ("fence" ::: "memory");
		num_iterations++;

#ifdef ZERO_PACKETS_V2
		make_submission_request_or_work(48, 0, numPendingWorkRequests);
		submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
		make_submission_request_or_work(12, 0, numPendingWorkRequests);
		submit_three_or_work(swID, 3, numPendingWorkRequests);
#endif

		submit_three_or_work((unsigned long long) &(_a[j]), 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) &(_b[j]), 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) &(_c[j]), 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
#endif
	}
}

// #pragma omp task in ([bs]b, [bs]c ) out ([bs]a)
void triad_task (double *a, double *b, double *c, double scalar, int bs)
{
	int j;	
	for (j=0; j < bs; j++)
            a[j] = b[j]+scalar*c[j];

}

void ORD_triad_task (uint64_t swID)
{
	uint64_t i = (uint64_t) metadataArray[swID].depAddresses0[0];
	double *a = &(_a[i]);
	double *b = &(_b[i]);
	double *c = &(_c[i]);

	int j;	
	for (j=0; j < BSIZE; j++)
            a[j] = b[j]+scalar*c[j];

}

void tuned_STREAM_Triad(double scalar)
{
	int j;
	uint64_t swID = 0;

	for (j=0; j<N; j+=BSIZE) {
		// triad_task (&a[j], &b[j], &c[j], scalar, BSIZE);
		swID = getNewSWID(swID);
		metadataArray[swID].functionAddr = (unsigned long long) ORD_triad_task;
		metadataArray[swID].depAddresses0[0] = (unsigned long long) j;

		asm volatile ("fence" ::: "memory");
		num_iterations++;

#ifdef ZERO_PACKETS_V2
		make_submission_request_or_work(48, 0, numPendingWorkRequests);
		submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
		make_submission_request_or_work(12, 0, numPendingWorkRequests);
		submit_three_or_work(swID, 3, numPendingWorkRequests);
#endif

		submit_three_or_work((unsigned long long) &(_a[j]), 2, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) &(_b[j]), 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) &(_c[j]), 0, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
		submit_three_or_work(0, 0, numPendingWorkRequests);
#endif
	}
}

