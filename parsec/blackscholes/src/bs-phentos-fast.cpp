// Copyright (c) 2007 Intel Corp.

// Black-Scholes
// Analytical method for calculating European Options
//
//
// Reference Source: Options, Futures, and Other Derivatives, 3rd Edition, Prentice
// Hall, John C. Hull,
//
// OmpSs/OpenMP 4.0 versions written by Dimitrios Chasapis and Iulian Brumar - Barcelona Supercomputing Center

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>

#include "../../../c/bench.h"
#include "../../../../phentos/phentos/femtos.hpp"
int allow_out;

#define getOneParameter(swID, a)      \
	{                                 \
		a = (swID >> 4) & MASK15BITS; \
	}

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

// Multi-threaded OpenMP header
//#include <omp.h>

#define MASK15BITS 0x7fff
#define BSIZE_UNIT 1024
int BSIZE;

//Precision to use for calculations
#define fptype float

#define NUM_RUNS 100

typedef struct OptionData_
{
	fptype s;		 // spot price
	fptype strike;	 // strike price
	fptype r;		 // risk-free interest rate
	fptype divq;	 // dividend rate
	fptype v;		 // volatility
	fptype t;		 // time to maturity or option expiration in years
					 //     (1yr = 1.0, 6mos = 0.5, 3mos = 0.25, ..., etc)
	char OptionType; // Option type.  "P"=PUT, "C"=CALL
	fptype divs;	 // dividend vals (not used in this test)
	fptype DGrefval; // DerivaGem Reference Value
} OptionData;

OptionData *data;
fptype *prices;
int numOptions;

int *otype;
fptype *sptprice;
fptype *strike;
fptype *rate;
fptype *volatility;
fptype *otime;
int numError = 0;
//int nThreads;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Cumulative Normal Distribution Function
// See Hull, Section 11.8, P.243-244
#define inv_sqrt_2xPI 0.39894228040143270286

void extra_init()
{
	totalNumberOfSubmittedTasks = 0;
	femtos_fast_init();
#ifdef ZERO_PACKETS_V2
	#pragma message ("Benchmark compiled with support for ZERO_PACKETS_V2")
	#pragma message ("Benchmark compiled with support for SW-BASED PADDING")
	#pragma message ("Benchmark compiled with support for RESET COUNT OF SUBMITTED TASKS")
	std::cout << "[feature]: ZERO_PACKETS_V2" << std::endl;
	std::cout << "[feature]: SW-BASED PADDING" << std::endl;
	std::cout << "[feature]: RESET COUNT OF SUBMITTED TASKS" << std::endl;
#endif
}

fptype CNDF(fptype InputX)
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
	if (InputX < 0.0)
	{
		InputX = -InputX;
		sign = 1;
	}
	else
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
	xLocal = xLocal_1 * xNPrimeofX;
	xLocal = 1.0 - xLocal;

	OutputX = xLocal;

	if (sign)
	{
		OutputX = 1.0 - OutputX;
	}

	return OutputX;
}

// For debugging
void print_xmm(fptype in, char *s)
{
	if (allow_out)
		printf("%s: %f\n", s, in);
}

void FA_BlkSchlsEqEuroNoDiv(unsigned long long swID)
{
	int i;
	int pOffset;
	getOneParameter(swID, pOffset);
	
	float timet = 0;
	int &bsize = BSIZE;

	for (i = 0; i < bsize; i++)
	{
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

		xStockPrice = sptprice[i + pOffset];
		xStrikePrice = strike[i + pOffset];
		xRiskFreeRate = rate[i + pOffset];
		xVolatility = volatility[i + pOffset];

		xTime = otime[i + pOffset];
		xSqrtTime = sqrt(xTime);

		logValues = log(sptprice[i + pOffset] / strike[i + pOffset]);

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

		FutureValueX = strike[i + pOffset] * (exp(-(rate[i + pOffset]) * (otime[i + pOffset])));
		if (otype[i + pOffset] == 0)
		{
			prices[i + pOffset] = (sptprice[i + pOffset] * NofXd1) - (FutureValueX * NofXd2);
		}
		else
		{
			NegNofXd1 = (1.0 - NofXd1);
			NegNofXd2 = (1.0 - NofXd2);
			prices[i + pOffset] = (FutureValueX * NegNofXd2) - (sptprice[i + pOffset] * NegNofXd1);
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void bs_thread(void *tid_ptr, fptype *prices)
{
	int i, j;
	fptype priceDelta;
	int tid = *(int *)tid_ptr;

	unsigned long long swID = 0;
	unsigned num_iterations = 0;
	unsigned numPendingWorkRequests = 0;

	functionAddresses[0] = (unsigned long long)FA_BlkSchlsEqEuroNoDiv;
	asm volatile("fence" ::
					 : "memory");

#ifdef DEBUG
	printf("[bs_thread]: We are about to start the first run.\n");
#endif

	for (j = 0; j < NUM_RUNS; j++)
	{
#ifdef DEBUG
		printf("[bs_thread]: Run #%d\n", j);
#endif
		for (i = 0; i <= (numOptions - BSIZE); i += BSIZE)
		{

#ifdef DEBUG
			printf("[bs_thread, run %d, stride %d]: Acquiring swID\n", j, i);
#endif
			swID = ((i & MASK15BITS) << 4) | 0x1;
#ifdef DEBUG
			printf("[bs_thread, run %d, stride %d]: Acquired swID: %llu\n", j, i, swID);
#endif

			num_iterations++;
#ifdef ZERO_PACKETS_V2
				make_submission_request_or_work_fast(48, 0, numPendingWorkRequests);
				submit_three_or_work_fast(swID, 15, numPendingWorkRequests);
#else
			make_submission_request_or_work_fast(24, 0, numPendingWorkRequests);
			submit_three_or_work_fast(swID, 7, numPendingWorkRequests);
#endif

			submit_three_or_work_fast((unsigned long long)(&sptprice[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&strike[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&rate[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&volatility[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&otime[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&otype[i]), 0, numPendingWorkRequests);
			submit_three_or_work_fast((unsigned long long)(&prices[i]), 1, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
			submit_three_or_work_fast(0, 0, numPendingWorkRequests);
#endif

#ifdef DEBUG
			printf("[bs_thread, run %d, stride %d]: We have just sent all dependences\n", j, i);
#endif
		}

		//We put a barrier here to avoid overlapping the execution of
		// tasks in different runs
		printf("Going to task wait until %d tasks were retired.\n", num_iterations);
		fast_task_wait_and_try_executing_tasks(num_iterations);

#ifdef ERR_CHK
		for (i = 0; i < numOptions; i++)
		{
			priceDelta = data[i].DGrefval - prices[i];
			if (fabs(priceDelta) >= 1e-4)
			{
				if (allow_out)
					printf("Error on %d. Computed=%.5f, Ref=%.5f, Delta=%.5f\n",
						   i, prices[i], data[i].DGrefval, priceDelta);
				numError++;
			}
		}
#endif
	}
}

int main(int argc, char **argv)
{
	allow_out = 1;
	if (getenv("BENCH_SILENT") != NULL)
		allow_out = 0;

	process_name("parsec-blackscholes");
	process_mode(OMPSS);
	process_args(argc, argv);
	process_init();

	FILE *file;
	int i;
	int loopnum;
	fptype *buffer;
	int *buffer2;
	int rv;
	struct timeval start;
	struct timeval stop;
	unsigned long elapsed;

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_begin(__parsec_blackscholes);
#endif

	if (argc < 4)
	{
		if (allow_out)
			printf("Usage:\n\t%s <nthreads> <inputFile> <outputFile> [blocksize]\n", argv[0]);
		if (allow_out)
			printf("Warning: nthreads is ignored! Use NX_ARGS=\"--threads=<nthreads>\" instead\n");
		exit(1);
	}

	char *inputFile = argv[2];
	char *outputFile = argv[3];
	int nthreads = atoi(argv[1]);

	if (argc > 4)
	{
		BSIZE = atoi(argv[4]);
	}
	else
	{
		BSIZE = BSIZE_UNIT;
	}

	//Read input data from file
	file = fopen(inputFile, "r");
	if (file == NULL)
	{
		if (allow_out)
			printf("ERROR: Unable to open file `%s'.\n", inputFile);
		exit(1);
	}
	rv = fscanf(file, "%i", &numOptions);
	if (rv != 1)
	{
		if (allow_out)
			printf("ERROR: Unable to read from file `%s'.\n", inputFile);
		fclose(file);
		exit(1);
	}
	if (BSIZE > numOptions)
	{
		if (allow_out)
			printf("ERROR: Block size larger than number of options. Please reduce the block size, or use larger data size.\n");
		exit(1);
		//printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
		//nThreads = numOptions;
	}
	if (numOptions % BSIZE)
	{
		if (allow_out)
			printf("ERROR: Number of options is not a multiple of block size.\n");
		exit(1);
		//printf("WARNING: Not enough work, reducing number of threads to match number of options.\n");
		//nThreads = numOptions;
	}

	// alloc spaces for the option data
	data = (OptionData *)malloc(numOptions * sizeof(OptionData));
	prices = (fptype *)malloc(numOptions * sizeof(fptype));
	for (loopnum = 0; loopnum < numOptions; ++loopnum)
	{
		rv = fscanf(file, "%f %f %f %f %f %f %c %f %f", &data[loopnum].s, &data[loopnum].strike, &data[loopnum].r, &data[loopnum].divq, &data[loopnum].v, &data[loopnum].t, &data[loopnum].OptionType, &data[loopnum].divs, &data[loopnum].DGrefval);
		if (rv != 9)
		{
			if (allow_out)
				printf("ERROR: Unable to read from file `%s'.\n", inputFile);
			fclose(file);
			exit(1);
		}
	}
	rv = fclose(file);
	if (rv != 0)
	{
		if (allow_out)
			printf("ERROR: Unable to close file `%s'.\n", inputFile);
		exit(1);
	}

	if (allow_out)
		printf("Num of Options: %d\n", numOptions);
	if (allow_out)
		printf("Num of Runs: %d\n", NUM_RUNS);

#define PAD 256
#define LINESIZE 64

	buffer = (fptype *)malloc(5 * numOptions * sizeof(fptype) + PAD);
	sptprice = (fptype *)(((unsigned long long)buffer + PAD) & ~(LINESIZE - 1));
	strike = sptprice + numOptions;
	rate = strike + numOptions;
	volatility = rate + numOptions;
	otime = volatility + numOptions;

	buffer2 = (int *)malloc(numOptions * sizeof(fptype) + PAD);
	otype = (int *)(((unsigned long long)buffer2 + PAD) & ~(LINESIZE - 1));

	for (i = 0; i < numOptions; i++)
	{
		otype[i] = (data[i].OptionType == 'P') ? 1 : 0;
		sptprice[i] = data[i].s;
		strike[i] = data[i].strike;
		rate[i] = data[i].r;
		volatility[i] = data[i].v;
		otime[i] = data[i].t;
	}

	if (allow_out)
		printf("Size of data: %d\n", numOptions * (sizeof(OptionData) + sizeof(int)));

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_begin();
#endif

	extra_init();

	//do work
	int tid = 0;
	//omp_set_num_threads(nThreads);
	gettimeofday(&start, NULL);
	process_start_measure();
	bs_thread(&tid, prices);
	process_stop_measure();
	gettimeofday(&stop, NULL);

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_end();
#endif

	//Write prices to output file
	file = fopen(outputFile, "w");
	if (file == NULL)
	{
		if (allow_out)
			printf("ERROR: Unable to open file `%s'.\n", outputFile);
		exit(1);
	}
	rv = fprintf(file, "%i\n", numOptions);
	if (rv < 0)
	{
		if (allow_out)
			printf("ERROR: Unable to write to file `%s'.\n", outputFile);
		fclose(file);
		exit(1);
	}
	for (i = 0; i < numOptions; i++)
	{
		rv = fprintf(file, "%.18f\n", prices[i]);
		if (rv < 0)
		{
			if (allow_out)
				printf("ERROR: Unable to write to file `%s'.\n", outputFile);
			fclose(file);
			exit(1);
		}
	}
	rv = fclose(file);
	if (rv != 0)
	{
		if (allow_out)
			printf("ERROR: Unable to close file `%s'.\n", outputFile);
		exit(1);
	}

#ifdef ERR_CHK
	if (allow_out)
		printf("Num Errors: %d\n", numError);
#endif
	free(data);
	free(prices);

	elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
	elapsed += stop.tv_usec - start.tv_usec;

	if (allow_out)
		printf("par_sec_time_us:%lu\n", elapsed);

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_end();
#endif

	process_append_file(outputFile);

	dump_csv(stdout);

	return 0;
}
