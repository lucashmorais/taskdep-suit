// Copyright (c) 2007 Intel Corp.
// Black-Scholes
// Analytical method for calculating European Options
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice 
// Hall, John C. Hull,

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "../../../c/bench.h"
#include "../../../../phentos/phentos/femtos.hpp"
int allow_out;


// Multi-threaded pthreads header
#ifdef ENABLE_THREADS
#define MAX_THREADS 128
// Add the following line so that icc 9.0 is compatible with pthread lib.
#define __thread __threadp
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K
#endif
#ifndef __USE_UNIX98
#define __USE_UNIX98
#endif
#include <pthread.h>
#include <time.h>

pthread_t _M4_threadsTable[MAX_THREADS];
pthread_mutexattr_t _M4_normalMutexAttr;
int _M4_numThreads = MAX_THREADS;
#undef __thread
#endif

// Multi-threaded OpenMP header
#ifdef ENABLE_OPENMP
#include <omp.h>
#endif

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 100

typedef struct OptionData_ {
        fptype s;          // spot price
        fptype strike;     // strike price
        fptype r;          // risk-free interest rate
        fptype divq;       // dividend rate
        fptype v;          // volatility
        fptype t;          // time to maturity or option expiration in years 
                           //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)  
        char OptionType;   // Option type.  "P"=PUT, "C"=CALL
        fptype divs;       // dividend vals (not used in this test)
        fptype DGrefval;   // DerivaGem Reference Value
} OptionData;

OptionData *data;
fptype *prices;
int numOptions;

int    * otype;
fptype * sptprice;
fptype * strike;
fptype * rate;
fptype * volatility;
fptype * otime;
int numError = 0;
int nThreads;

int BSIZE = 0;

////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

void extra_init() {
	femtos_init();
}

fptype CNDF ( fptype InputX ) 
{
    int sign;

    fptype OutputX;
    fptype xInput;
    fptype xNPrimeofX;
    fptype expValues;
    fptype xK2;
    fptype xK2_2, xK2_3;
    fptype xK2_4, xK2_5;
    fptype xLocal, xLocal_1;
    fptype xLocal_2, xLocal_3;

    // Check for negative value of InputX
    if (InputX < 0.0) {
        InputX = -InputX;
        sign = 1;
    } else 
        sign = 0;

    xInput = InputX;
 
    // Compute NPrimeX term common to both four & six decimal accuracy calcs
    expValues = exp(-0.5f * InputX * InputX);
    xNPrimeofX = expValues;
    xNPrimeofX = xNPrimeofX * inv_sqrt_2xPI;

    xK2 = 0.2316419 * xInput;
    xK2 = 1.0 + xK2;
    xK2 = 1.0 / xK2;
    xK2_2 = xK2 * xK2;
    xK2_3 = xK2_2 * xK2;
    xK2_4 = xK2_3 * xK2;
    xK2_5 = xK2_4 * xK2;
    
    xLocal_1 = xK2 * 0.319381530;
    xLocal_2 = xK2_2 * (-0.356563782);
    xLocal_3 = xK2_3 * 1.781477937;
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_4 * (-1.821255978);
    xLocal_2 = xLocal_2 + xLocal_3;
    xLocal_3 = xK2_5 * 1.330274429;
    xLocal_2 = xLocal_2 + xLocal_3;

    xLocal_1 = xLocal_2 + xLocal_1;
    xLocal   = xLocal_1 * xNPrimeofX;
    xLocal   = 1.0 - xLocal;

    OutputX  = xLocal;
    
    if (sign) {
        OutputX = 1.0 - OutputX;
    }
    
    return OutputX;
} 
void FA_BlkSchlsEqEuroNoDiv(unsigned long long swID) {
  fptype * sptprice = (fptype * ) metadataArray[swID].depAddresses0[0];
  fptype * strike = (fptype * ) metadataArray[swID].depAddresses0[1];
  fptype * rate = (fptype * ) metadataArray[swID].depAddresses0[2];
  fptype * volatility = (fptype * ) metadataArray[swID].depAddresses0[3];
  fptype * time = (fptype * ) metadataArray[swID].depAddresses0[4];
  int * otype = (int * ) metadataArray[swID].depAddresses0[5];
  fptype * OptionPrice = (fptype * ) metadataArray[swID].depAddresses0[6];
  float timet = 0;
  int & bsize = BSIZE;

  int i;

  for (i = 0; i < bsize; i++) {
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1;
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;

    xStockPrice = sptprice[i];
    xStrikePrice = strike[i];
    xRiskFreeRate = rate[i];
    xVolatility = volatility[i];

    xTime = time[i];
    xSqrtTime = sqrt(xTime);

    logValues = log(sptprice[i] / strike[i]);

    xLogTerm = logValues;

    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;

    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 - xDen;

    d1 = xD1;
    d2 = xD2;

    NofXd1 = CNDF(d1);
    NofXd2 = CNDF(d2);

    FutureValueX = strike[i] * (exp(-(rate[i]) * (time[i])));
    if (otype[i] == 0) {
      OptionPrice[i] = (sptprice[i] * NofXd1) - (FutureValueX * NofXd2);
    } else {
      NegNofXd1 = (1.0 - NofXd1);
      NegNofXd2 = (1.0 - NofXd2);
      OptionPrice[i] = (FutureValueX * NegNofXd2) - (sptprice[i] * NegNofXd1);
    }
  }

#ifdef ERR_CHK
	priceDelta = data[i].DGrefval - price;
	if( fabs(priceDelta) >= 1e-4 ){
		if(allow_out) printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
				i, price, data[i].DGrefval, priceDelta);
		numError ++;
	}
#endif
}

//////////////////////////////////////////////////////////////////////////////////////
fptype BlkSchlsEqEuroNoDiv( fptype sptprice,
                            fptype strike, fptype rate, fptype volatility,
                            fptype time, int otype, float timet )
{
    fptype OptionPrice;

    // local private working variables for the calculation
    fptype xStockPrice;
    fptype xStrikePrice;
    fptype xRiskFreeRate;
    fptype xVolatility;
    fptype xTime;
    fptype xSqrtTime;

    fptype logValues;
    fptype xLogTerm;
    fptype xD1; 
    fptype xD2;
    fptype xPowerTerm;
    fptype xDen;
    fptype d1;
    fptype d2;
    fptype FutureValueX;
    fptype NofXd1;
    fptype NofXd2;
    fptype NegNofXd1;
    fptype NegNofXd2;    
    
    xStockPrice = sptprice;
    xStrikePrice = strike;
    xRiskFreeRate = rate;
    xVolatility = volatility;

    xTime = time;
    xSqrtTime = sqrt(xTime);

    logValues = log( sptprice / strike );
        
    xLogTerm = logValues;
        
    
    xPowerTerm = xVolatility * xVolatility;
    xPowerTerm = xPowerTerm * 0.5;
        
    xD1 = xRiskFreeRate + xPowerTerm;
    xD1 = xD1 * xTime;
    xD1 = xD1 + xLogTerm;

    xDen = xVolatility * xSqrtTime;
    xD1 = xD1 / xDen;
    xD2 = xD1 -  xDen;

    d1 = xD1;
    d2 = xD2;
    
    NofXd1 = CNDF( d1 );
    NofXd2 = CNDF( d2 );

    FutureValueX = strike * ( exp( -(rate)*(time) ) );        
    if (otype == 0) {            
        OptionPrice = (sptprice * NofXd1) - (FutureValueX * NofXd2);
    } else { 
        NegNofXd1 = (1.0 - NofXd1);
        NegNofXd2 = (1.0 - NofXd2);
        OptionPrice = (FutureValueX * NegNofXd2) - (sptprice * NegNofXd1);
    }
    
    return OptionPrice;
}


int bs_thread(void *tid_ptr) {
    int i, j;
    fptype price;
    fptype priceDelta;
    int tid = *(int *)tid_ptr;
	BSIZE = numOptions;

	unsigned long long swID = 0;
	unsigned num_iterations = 0;
	unsigned numPendingWorkRequests = 0;


    for (j=0; j<NUM_RUNS; j++) {
		swID = getNewSWID(swID);
		/*
		fptype * sptprice = (fptype * ) metadataArray[swID].depAddresses0[0];
		fptype * strike = (fptype * ) metadataArray[swID].depAddresses0[1];
		fptype * rate = (fptype * ) metadataArray[swID].depAddresses0[2];
		fptype * volatility = (fptype * ) metadataArray[swID].depAddresses0[3];
		fptype * time = (fptype * ) metadataArray[swID].depAddresses0[4];
		int * otype = (int * ) metadataArray[swID].depAddresses0[5];
		fptype * OptionPrice = (fptype * ) metadataArray[swID].depAddresses0[6];
		*/

		metadataArray[swID].functionAddr = (unsigned long long) FA_BlkSchlsEqEuroNoDiv;
		metadataArray[swID].depAddresses0[0] = (unsigned long long) sptprice;
		metadataArray[swID].depAddresses0[1] = (unsigned long long) strike;
		metadataArray[swID].depAddresses0[2] = (unsigned long long) rate;
		metadataArray[swID].depAddresses0[3] = (unsigned long long) volatility;
		metadataArray[swID].depAddresses0[4] = (unsigned long long) time;
		metadataArray[swID].depAddresses0[5] = (unsigned long long) otype;
		metadataArray[swID].depAddresses0[6] = (unsigned long long) prices;

		asm volatile ("fence" ::: "memory");

		num_iterations++;
		make_submission_request_or_work(24, 0, numPendingWorkRequests);
		submit_three_or_work(swID, 7, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) sptprice, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) strike, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) rate, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) volatility, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) time, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) otype, 0, numPendingWorkRequests);
		submit_three_or_work((unsigned long long) prices, 1, numPendingWorkRequests);

		price = BlkSchlsEqEuroNoDiv( sptprice[i], strike[i],
										rate[i], volatility[i], otime[i], 
										otype[i], 0);
    }

	printf("Going to task wait until %d tasks were retired.\n", num_iterations);
	task_wait_and_try_executing_tasks(num_iterations);

    return 0;
}

int main (int argc, char **argv)
{
    allow_out = 1;
    if(getenv("BENCH_SILENT") != NULL) allow_out = 0;

    process_name("parsec-blackscholes");
    process_mode(SEQ);
    process_args(argc, argv);
    process_init();

    FILE *file;
    int i;
    int loopnum;
    fptype * buffer;
    int * buffer2;
    int rv;

   if (argc < 4)
        {
                if(allow_out) printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile>\n", argv[0]);
                exit(1);
        }
    nThreads = atoi(argv[1]);
    char *inputFile = argv[2];
    char *outputFile = argv[3];

    //Read input data from file
    file = fopen(inputFile, "r");
    if(file == NULL) {
      if(allow_out) printf("ERROR: Unable to open file `%s'.\n", inputFile);
      exit(1);
    }
    rv = fscanf(file, "%i", &numOptions);
    if(rv != 1) {
      if(allow_out) printf("ERROR: Unable to read from file `%s'.\n", inputFile);
      fclose(file);
      exit(1);
    }
    if(nThreads > numOptions) {
      /*printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");*/
      nThreads = numOptions;
    }

#if !defined(ENABLE_THREADS) && !defined(ENABLE_OPENMP)
    if(nThreads != 1) {
        if(allow_out) printf("Error: <nthreads> must be 1 (serial version)\n");
        exit(1);
    }
#endif

#ifdef _PHENTOS
	extra_init();
#endif

    // alloc spaces for the option data
    data = (OptionData*)malloc(numOptions*sizeof(OptionData));
    prices = (fptype*)malloc(numOptions*sizeof(fptype));
    for ( loopnum = 0; loopnum < numOptions; ++ loopnum )
    {
        rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
        if(rv != 9) {
          if(allow_out) printf("ERROR: Unable to read from file `%s'.\n", inputFile);
          fclose(file);
          exit(1);
        }
    }
    rv = fclose(file);
    if(rv != 0) {
      if(allow_out) printf("ERROR: Unable to close file `%s'.\n", inputFile);
      exit(1);
    }

#ifdef ENABLE_THREADS
    pthread_mutexattr_init( &_M4_normalMutexAttr);
    //    pthread_mutexattr_settype( &_M4_normalMutexAttr, PTHREAD_MUTEX_NORMAL);
    _M4_numThreads = nThreads;
    {
        int _M4_i;
        for ( _M4_i = 0; _M4_i < MAX_THREADS; _M4_i++) {
            _M4_threadsTable[_M4_i] = -1;
        }
    }
#endif
    if(allow_out) printf("Num of Options: %d\n", numOptions);
    if(allow_out) printf("Num of Runs: %d\n", NUM_RUNS);

#define PAD 256
#define LINESIZE 64

    buffer = (fptype *) malloc(5 * numOptions * sizeof(fptype) + PAD);
    sptprice = (fptype *) (((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
    strike = sptprice + numOptions;
    rate = strike + numOptions;
    volatility = rate + numOptions;
    otime = volatility + numOptions;

    buffer2 = (int *) malloc(numOptions * sizeof(fptype) + PAD);
    otype = (int *) (((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

    for (i=0; i<numOptions; i++) {
        otype[i]      = (data[i].OptionType == 'P') ? 1 : 0;
        sptprice[i]   = data[i].s;
        strike[i]     = data[i].strike;
        rate[i]       = data[i].r;
        volatility[i] = data[i].v;    
        otime[i]      = data[i].t;
    }

    if(allow_out) printf("Size of data: %d\n", numOptions * (sizeof(OptionData) + sizeof(int)));

#ifdef ENABLE_THREADS
    int *tids;
    tids = (int *) malloc (nThreads * sizeof(int));

    for(i=0; i<nThreads; i++) {
        tids[i]=i;
        int _M4_i;
        for ( _M4_i = 0; _M4_i < MAX_THREADS; _M4_i++) {
            if ( _M4_threadsTable[_M4_i] == -1)    break;
        }
        pthread_create(&_M4_threadsTable[_M4_i],NULL,(void *(*)(void *))bs_thread,(void *)&tids[i]);
		// if(allow_out) printf("tid=%d %d\n", _M4_i, _M4_threadsTable[_M4_i]);
    }
    int _M4_i;
    void *_M4_ret;
    for ( _M4_i = 0; _M4_i < MAX_THREADS;_M4_i++) {
        // if(allow_out) printf("tid=%d %d\n", _M4_i, _M4_threadsTable[_M4_i]);
        if ( _M4_threadsTable[_M4_i] == -1)    break;
        pthread_join( _M4_threadsTable[_M4_i], &_M4_ret);
    }
    free(tids);
#else //ENABLE_THREADS
#ifdef ENABLE_OPENMP
    {
        int tid=0;
        omp_set_num_threads(nThreads);
        bs_thread(&tid);
    }
#else //ENABLE_OPENMP
    //serial version
    int tid=0;
    process_start_measure();
    bs_thread(&tid);
    process_stop_measure();
#endif //ENABLE_OPENMP
#endif //ENABLE_THREADS

    //Write prices to output file
    file = fopen(outputFile, "w");
    if(file == NULL) {
      if(allow_out) printf("ERROR: Unable to open file `%s'.\n", outputFile);
      exit(1);
    }
    rv = fprintf(file, "%i\n", numOptions);
    if(rv < 0) {
      if(allow_out) printf("ERROR: Unable to write to file `%s'.\n", outputFile);
      fclose(file);
      exit(1);
    }
    for(i=0; i<numOptions; i++) {
      rv = fprintf(file, "%.18f\n", prices[i]);
      if(rv < 0) {
        if(allow_out) printf("ERROR: Unable to write to file `%s'.\n", outputFile);
        fclose(file);
        exit(1);
      }
    }
    rv = fclose(file);
    if(rv != 0) {
      if(allow_out) printf("ERROR: Unable to close file `%s'.\n", outputFile);
      exit(1);
    }

#ifdef ERR_CHK
    if(allow_out) printf("Num Errors: %d\n", numError);
#endif
    free(data);
    free(prices);

    process_append_file(outputFile);

    dump_csv(stdout);

    return 0;
}

