/**********************************************************************************************/
/*  This program is part of the Barcelona OpenMP Tasks Suite                                  */
/*  Copyright (C) 2009 Barcelona Supercomputing Center - Centro Nacional de Supercomputacion  */
/*  Copyright (C) 2009 Universitat Politecnica de Catalunya                                   */
/*                                                                                            */
/*  This program is free software; you can redistribute it and/or modify                      */
/*  it under the terms of the GNU General Public License as published by                      */
/*  the Free Software Foundation; either version 2 of the License, or                         */
/*  (at your option) any later version.                                                       */
/*                                                                                            */
/*  This program is distributed in the hope that it will be useful,                           */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of                            */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             */
/*  GNU General Public License for more details.                                              */
/*                                                                                            */
/*  You should have received a copy of the GNU General Public License                         */
/*  along with this program; if not, write to the Free Software                               */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA            */
/**********************************************************************************************/

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include "sparselu.h"
#include "../../../../phentos/phentos/femtos.hpp"

#define MASK15BITS 0x7fff

#define getOneParameter(swID, a) {\
	a = (swID >> 4) & MASK15BITS;\
}

#define getTwoParameters(swID, a, b) {\
	a = (swID >> 4) & MASK15BITS;\
	b = (swID >> (4 + 15)) & MASK15BITS;\
}

#define getThreeParameters(swID, a, b, c) {\
	a = (swID >> 4) & MASK15BITS;\
	b = (swID >> (4 + 15)) & MASK15BITS;\
	c = (swID >> (4 + 30)) & MASK15BITS;\
}

#define getFourParameters(swID, a, b, c, d) {\
	a = (swID >> 4) & MASK15BITS;\
	b = (swID >> (4 + 15)) & MASK15BITS;\
	c = (swID >> (4 + 30)) & MASK15BITS;\
	d = (swID >> (4 + 45)) & MASK15BITS;\
}

float ** GLOBAL_BENCH;
unsigned int GLOBAL_matrix_size;
unsigned int GLOBAL_submatrix_size;

void resetTaskWaitCounters() {
	numRetiredTasks = 0;
	totalNumberOfSubmittedTasks = 0;
	#pragma message ("Benchmark compiled with support for ZERO_PACKETS_V2")
	#pragma message ("Benchmark compiled with support for SW-BASED PADDING")
	#pragma message ("Benchmark compiled with support for RESET COUNT OF SUBMITTED TASKS")
#ifdef ZERO_PACKETS_V2
	std::cout << "[feature]: ZERO_PACKETS_V2" << std::endl;
	std::cout << "[feature]: SW-BASED PADDING" << std::endl;
	std::cout << "[feature]: RESET COUNT OF SUBMITTED TASKS" << std::endl;
#endif
}

#if 0
void lu0(float *diag, int submatrix_size)
{
    int i, j, k;

    for (k=0; k<submatrix_size; k++)
        for (i=k+1; i<submatrix_size; i++)
        {
            diag[i*submatrix_size+k] = diag[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                diag[i*submatrix_size+j] = diag[i*submatrix_size+j] - diag[i*submatrix_size+k] * diag[k*submatrix_size+j];
        }
}
#endif
void FA_lu0 (uint64_t swID) {
    int i, j, k;
    int kk;

	getOneParameter(swID, kk);

	int matrix_size = GLOBAL_matrix_size;
	int submatrix_size = GLOBAL_submatrix_size;

	float * diag = GLOBAL_BENCH[kk*matrix_size+kk];

    for (k=0; k<submatrix_size; k++)
        for (i=k+1; i<submatrix_size; i++)
        {
            diag[i*submatrix_size+k] = diag[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                diag[i*submatrix_size+j] = diag[i*submatrix_size+j] - diag[i*submatrix_size+k] * diag[k*submatrix_size+j];
        }
}

#if 0
void fwd(float *diag, float *col, int submatrix_size)
{
    int i, j, k;
    for (j=0; j<submatrix_size; j++)
        for (k=0; k<submatrix_size; k++)
            for (i=k+1; i<submatrix_size; i++)
                col[i*submatrix_size+j] = col[i*submatrix_size+j] - diag[i*submatrix_size+k]*col[k*submatrix_size+j];
}
#endif
void FA_fwd (uint64_t swID) {
    int i, j, k;
    int ii, jj, kk;

	getTwoParameters(swID, jj, kk);

	int matrix_size = GLOBAL_matrix_size;
	int submatrix_size = GLOBAL_submatrix_size;

	float * diag = GLOBAL_BENCH[kk*matrix_size+kk];
	float * col = GLOBAL_BENCH[kk*matrix_size+jj];

    for (j=0; j<submatrix_size; j++)
        for (k=0; k<submatrix_size; k++)
            for (i=k+1; i<submatrix_size; i++)
                col[i*submatrix_size+j] = col[i*submatrix_size+j] - diag[i*submatrix_size+k]*col[k*submatrix_size+j];
}

#if 0
void bdiv(float *diag, float *row, int submatrix_size)
{
    int i, j, k;
    for (i=0; i<submatrix_size; i++)
        for (k=0; k<submatrix_size; k++)
        {
            row[i*submatrix_size+k] = row[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                row[i*submatrix_size+j] = row[i*submatrix_size+j] - row[i*submatrix_size+k]*diag[k*submatrix_size+j];
        }
}
#endif
void FA_bdiv (uint64_t swID) {
    int i, j, k;
    int ii, jj, kk;

	getTwoParameters(swID, ii, kk);

	int matrix_size = GLOBAL_matrix_size;
	int submatrix_size = GLOBAL_submatrix_size;

	float * diag = GLOBAL_BENCH[kk*matrix_size+kk];
	float * row = GLOBAL_BENCH[ii*matrix_size+kk];

    for (i=0; i<submatrix_size; i++)
        for (k=0; k<submatrix_size; k++)
        {
            row[i*submatrix_size+k] = row[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                row[i*submatrix_size+j] = row[i*submatrix_size+j] - row[i*submatrix_size+k]*diag[k*submatrix_size+j];
        }
}

#if 0
void bmod(float *row, float *col, float *inner, int submatrix_size)
{
    int i, j, k;
    for (i=0; i<submatrix_size; i++)
        for (j=0; j<submatrix_size; j++)
            for (k=0; k<submatrix_size; k++)
                inner[i*submatrix_size+j] = inner[i*submatrix_size+j] - row[i*submatrix_size+k]*col[k*submatrix_size+j];
}
#endif
void FA_bmod (uint64_t swID) {
    int i, j, k;
	int ii, jj, kk;

	getThreeParameters(swID, jj, ii, kk);

	int matrix_size = GLOBAL_matrix_size;
	int submatrix_size = GLOBAL_submatrix_size;

	float * row = GLOBAL_BENCH[ii*matrix_size+kk];
	float * col = GLOBAL_BENCH[kk*matrix_size+jj];
	float * inner = GLOBAL_BENCH[ii*matrix_size+jj];

    for (i=0; i<submatrix_size; i++)
        for (j=0; j<submatrix_size; j++)
            for (k=0; k<submatrix_size; k++)
                inner[i*submatrix_size+j] = inner[i*submatrix_size+j] - row[i*submatrix_size+k]*col[k*submatrix_size+j];
}

void sparselu_par_call_core(float **BENCH, int matrix_size, int submatrix_size)
{
        unsigned num_iterations = 0;
        unsigned numPendingWorkRequests = 0;

        resetTaskWaitCounters();

        GLOBAL_BENCH = BENCH;
        GLOBAL_matrix_size = matrix_size;
        GLOBAL_submatrix_size = submatrix_size;

        functionAddresses[0] = (unsigned long long) FA_lu0;
        functionAddresses[1] = (unsigned long long) FA_fwd;
        functionAddresses[2] = (unsigned long long) FA_bdiv;
        functionAddresses[3] = (unsigned long long) FA_bmod;

        asm volatile ("fence" ::: "memory");

        unsigned long long ii, jj, kk;
        unsigned long long swID;
        {
                for (kk=0; kk<matrix_size; kk++) {

                        swID = ((kk & MASK15BITS) << 4) | 0x1;

                        num_iterations++;
#ifdef ZERO_PACKETS_V2
												make_submission_request_or_work_fast(48, 0, numPendingWorkRequests);
												submit_three_or_work_fast(swID, 15, numPendingWorkRequests);
#else
                        make_submission_request_or_work_fast(6, 0, numPendingWorkRequests);
                        submit_three_or_work_fast(swID, 1, numPendingWorkRequests);
#endif
                        submit_three_or_work_fast((unsigned long long) BENCH[kk*matrix_size+kk], 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
#endif


                        for (jj=kk+1; jj<matrix_size; jj++)
                                if (BENCH[kk*matrix_size+jj] != NULL) {

                                        swID = ((((kk & MASK15BITS) << 15) | (jj & MASK15BITS)) << 4) | 0x3;

                                        num_iterations++;
#ifdef ZERO_PACKETS_V2
																				make_submission_request_or_work_fast(48, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(swID, 15, numPendingWorkRequests);
#else
                                        make_submission_request_or_work_fast(9, 0, numPendingWorkRequests);
                                        submit_three_or_work_fast(swID, 2, numPendingWorkRequests);
#endif
                                        submit_three_or_work_fast((unsigned long long) BENCH[kk*matrix_size+kk], 0, numPendingWorkRequests);
                                        submit_three_or_work_fast((unsigned long long) BENCH[kk*matrix_size+jj], 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
#endif
                                }

                        for (ii=kk+1; ii<matrix_size; ii++)
                                if (BENCH[ii*matrix_size+kk] != NULL) {
                                        swID = ((((kk & MASK15BITS) << 15) | (ii & MASK15BITS)) << 4) | 0x5;

                                        num_iterations++;
#ifdef ZERO_PACKETS_V2
																				make_submission_request_or_work_fast(48, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(swID, 15, numPendingWorkRequests);
#else
                                        make_submission_request_or_work_fast(9, 0, numPendingWorkRequests);
                                        submit_three_or_work_fast(swID, 2, numPendingWorkRequests);
#endif
                                        submit_three_or_work_fast((unsigned long long) BENCH[kk*matrix_size+kk], 0, numPendingWorkRequests);
                                        submit_three_or_work_fast((unsigned long long) BENCH[ii*matrix_size+kk], 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																				submit_three_or_work_fast(0, 0, numPendingWorkRequests);
#endif
                                }

                        for (ii=kk+1; ii<matrix_size; ii++)
                                if (BENCH[ii*matrix_size+kk] != NULL) {
                                        for (jj=kk+1; jj<matrix_size; jj++)
                                                if (BENCH[kk*matrix_size+jj] != NULL) {
                                                        if (BENCH[ii*matrix_size+jj]==NULL) {
                                                                BENCH[ii*matrix_size+jj] = allocate_clean_block(submatrix_size);
                                                        }

                                                        swID = ((  ((kk & MASK15BITS) << 30) | ((ii & MASK15BITS) << 15) | (jj & MASK15BITS)  ) << 4) | 0x7;

                                                        num_iterations++;
#ifdef ZERO_PACKETS_V2
																												make_submission_request_or_work_fast(48, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(swID, 15, numPendingWorkRequests);
#else
                                                        make_submission_request_or_work_fast(12, 0, numPendingWorkRequests);
                                                        submit_three_or_work_fast(swID, 3, numPendingWorkRequests);
#endif
                                                        submit_three_or_work_fast((unsigned long long) BENCH[ii*matrix_size+kk], 0, numPendingWorkRequests);
                                                        submit_three_or_work_fast((unsigned long long) BENCH[kk*matrix_size+jj], 0, numPendingWorkRequests);

                                                        asm volatile ("fence" ::: "memory");

                                                        submit_three_or_work_fast((unsigned long long) BENCH[ii*matrix_size+jj], 2, numPendingWorkRequests);
#ifdef ZERO_PACKETS_V2
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
																												submit_three_or_work_fast(0, 0, numPendingWorkRequests);
#endif
                                                }
                                }
                }
        }

        printf("Going to task wait until %d tasks were retired.\n", num_iterations);
        fast_task_wait_and_try_executing_tasks(num_iterations);
}

void extra_init() {
	femtos_fast_init();
}

void sparselu_par_call(float **BENCH, int matrix_size, int submatrix_size) {
	sparselu_par_call_core(BENCH, matrix_size, submatrix_size);
}
