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

void resetTaskWaitCounters() {
	numRetiredTasks = 0;
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
void ORD_lu0 (uint64_t swID) {
    int i, j, k;
    int kk;

	asm volatile ("fence" ::: "memory");
	float * diag = (float *) metadataArray[swID].depAddresses0[0];
	int submatrix_size = metadataArray[swID].depAddresses0[1];

    for (k=0; k<submatrix_size; k++)
        for (i=k+1; i<submatrix_size; i++)
        {
            diag[i*submatrix_size+k] = diag[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                diag[i*submatrix_size+j] = diag[i*submatrix_size+j] - diag[i*submatrix_size+k] * diag[k*submatrix_size+j];
        }
	asm volatile ("fence" ::: "memory");
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
void ORD_fwd (uint64_t swID) {
    int i, j, k;
    int ii, jj, kk;

	asm volatile ("fence" ::: "memory");
	float * diag = (float *) metadataArray[swID].depAddresses0[0];
	float * col = (float *) metadataArray[swID].depAddresses0[1];
	int submatrix_size = metadataArray[swID].depAddresses0[2];

    for (j=0; j<submatrix_size; j++)
        for (k=0; k<submatrix_size; k++)
            for (i=k+1; i<submatrix_size; i++)
                col[i*submatrix_size+j] = col[i*submatrix_size+j] - diag[i*submatrix_size+k]*col[k*submatrix_size+j];
	asm volatile ("fence" ::: "memory");
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
void ORD_bdiv (uint64_t swID) {
    int i, j, k;
    int ii, jj, kk;

	asm volatile ("fence" ::: "memory");
	float * diag = (float *) metadataArray[swID].depAddresses0[0];
	float * row = (float *) metadataArray[swID].depAddresses0[1];
	int submatrix_size = metadataArray[swID].depAddresses0[2];

    for (i=0; i<submatrix_size; i++)
        for (k=0; k<submatrix_size; k++)
        {
            row[i*submatrix_size+k] = row[i*submatrix_size+k] / diag[k*submatrix_size+k];
            for (j=k+1; j<submatrix_size; j++)
                row[i*submatrix_size+j] = row[i*submatrix_size+j] - row[i*submatrix_size+k]*diag[k*submatrix_size+j];
        }
	asm volatile ("fence" ::: "memory");
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
void ORD_bmod (uint64_t swID) {
    int i, j, k;
	int ii, jj, kk;

	asm volatile ("fence" ::: "memory");
	float * row = (float *) metadataArray[swID].depAddresses0[0];
	float * col = (float *) metadataArray[swID].depAddresses0[1];
	float * inner = (float *) metadataArray[swID].depAddresses0[2];
	int submatrix_size = metadataArray[swID].depAddresses0[3];

    for (i=0; i<submatrix_size; i++)
        for (j=0; j<submatrix_size; j++)
            for (k=0; k<submatrix_size; k++)
                inner[i*submatrix_size+j] = inner[i*submatrix_size+j] - row[i*submatrix_size+k]*col[k*submatrix_size+j];
	asm volatile ("fence" ::: "memory");
}

void sparselu_par_call_core(float **BENCH, int matrix_size, int submatrix_size)
{
        unsigned num_iterations = 0;
        unsigned numPendingWorkRequests = 0;

        resetTaskWaitCounters();

        functionAddresses[0] = (unsigned long long) ORD_lu0;
        functionAddresses[1] = (unsigned long long) ORD_fwd;
        functionAddresses[2] = (unsigned long long) ORD_bdiv;
        functionAddresses[3] = (unsigned long long) ORD_bmod;
        asm volatile ("fence" ::: "memory");

        unsigned long long ii, jj, kk;
        unsigned long long swID = 0;
        {
                for (kk=0; kk<matrix_size; kk++) {

						swID = getNewSWID(swID);

						metadataArray[swID].functionAddr = (unsigned long long) ORD_lu0;
						metadataArray[swID].depAddresses0[0] = (unsigned long long) BENCH[kk*matrix_size+kk];
						metadataArray[swID].depAddresses0[1] = (unsigned long long) submatrix_size;
						asm volatile ("fence" ::: "memory");

                        num_iterations++;
                        make_submission_request_or_work(6, 0, numPendingWorkRequests);
                        submit_three_or_work(swID, 1, numPendingWorkRequests);
                        submit_three_or_work((unsigned long long) BENCH[kk*matrix_size+kk], 2, numPendingWorkRequests);

                        for (jj=kk+1; jj<matrix_size; jj++)
                                if (BENCH[kk*matrix_size+jj] != NULL) {

										swID = getNewSWID(swID);

										metadataArray[swID].functionAddr = (unsigned long long) ORD_fwd;
										metadataArray[swID].depAddresses0[0] = (unsigned long long) BENCH[kk*matrix_size+kk];
										metadataArray[swID].depAddresses0[1] = (unsigned long long) BENCH[kk*matrix_size+jj];
										metadataArray[swID].depAddresses0[2] = (unsigned long long) submatrix_size;
										asm volatile ("fence" ::: "memory");

                                        num_iterations++;
                                        make_submission_request_or_work(9, 0, numPendingWorkRequests);
                                        submit_three_or_work(swID, 2, numPendingWorkRequests);
                                        submit_three_or_work((unsigned long long) BENCH[kk*matrix_size+kk], 0, numPendingWorkRequests);
                                        submit_three_or_work((unsigned long long) BENCH[kk*matrix_size+jj], 2, numPendingWorkRequests);
                                }

                        for (ii=kk+1; ii<matrix_size; ii++)
                                if (BENCH[ii*matrix_size+kk] != NULL) {
										swID = getNewSWID(swID);

										metadataArray[swID].functionAddr = (unsigned long long) ORD_bdiv;
										metadataArray[swID].depAddresses0[0] = (unsigned long long) BENCH[kk*matrix_size+kk];
										metadataArray[swID].depAddresses0[1] = (unsigned long long) BENCH[ii*matrix_size+kk];
										metadataArray[swID].depAddresses0[2] = (unsigned long long) submatrix_size;
										asm volatile ("fence" ::: "memory");

                                        num_iterations++;
                                        make_submission_request_or_work(9, 0, numPendingWorkRequests);
                                        submit_three_or_work(swID, 2, numPendingWorkRequests);
                                        submit_three_or_work((unsigned long long) BENCH[kk*matrix_size+kk], 0, numPendingWorkRequests);
                                        submit_three_or_work((unsigned long long) BENCH[ii*matrix_size+kk], 2, numPendingWorkRequests);
                                }

                        for (ii=kk+1; ii<matrix_size; ii++)
                                if (BENCH[ii*matrix_size+kk] != NULL) {
                                        for (jj=kk+1; jj<matrix_size; jj++)
                                                if (BENCH[kk*matrix_size+jj] != NULL) {
                                                        if (BENCH[ii*matrix_size+jj]==NULL) {
                                                                BENCH[ii*matrix_size+jj] = allocate_clean_block(submatrix_size);
                                                        }

														swID = getNewSWID(swID);

														metadataArray[swID].functionAddr = (unsigned long long) ORD_bmod;
														metadataArray[swID].depAddresses0[0] = (unsigned long long) BENCH[ii*matrix_size+kk];
														metadataArray[swID].depAddresses0[1] = (unsigned long long) BENCH[kk*matrix_size+jj];
														metadataArray[swID].depAddresses0[1] = (unsigned long long) BENCH[ii*matrix_size+jj];
														metadataArray[swID].depAddresses0[2] = (unsigned long long) submatrix_size;
														asm volatile ("fence" ::: "memory");

                                                        num_iterations++;
                                                        make_submission_request_or_work(12, 0, numPendingWorkRequests);
                                                        submit_three_or_work(swID, 3, numPendingWorkRequests);
                                                        submit_three_or_work((unsigned long long) BENCH[ii*matrix_size+kk], 0, numPendingWorkRequests);
                                                        submit_three_or_work((unsigned long long) BENCH[kk*matrix_size+jj], 0, numPendingWorkRequests);

														asm volatile ("fence" ::: "memory");

                                                        submit_three_or_work((unsigned long long) BENCH[ii*matrix_size+jj], 2, numPendingWorkRequests);
                                                }
                                }
                }
        }

        printf("Going to task wait until %d tasks were retired.\n", num_iterations);
        task_wait_and_try_executing_tasks(num_iterations);
}

void sparselu_par_call(float **BENCH, int matrix_size, int submatrix_size) {
	femtos_init();
	sparselu_par_call_core(BENCH, matrix_size, submatrix_size);
}