# include "poisson.h"
#include "../../../c/bench.h"
#include "../../../../phentos/phentos/femtos.hpp"

double *f_;
double *u_;
double *unew_;

int nx, ny;
double dx, dy;

#define MASK15BITS 0x7fff
#define MASK3BITS 0x7

#define getOneParameter(swID, a) {\
	a = (swID >> 4) & MASK15BITS;\
}

#define buildSWIDWithOneParameter(swID, funcIdx, a) {\
	swID = ((a & MASK15BITS) << 4) | ((funcIdx & MASK3BITS) << 1) | 0x1;\
}


void sweep_partial_a(double *u_array, double *unew_array) {
	for (int j = 0; j < ny; j++) {
		u_array[j] = unew_array[j];
	}
}

void FAST_sweep_partial_a(uint64_t swID) {
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;
	unsigned i;

	getOneParameter(swID, i);

	for (int j = 0; j < ny; j++) {
		(*u)[i][j] = (*unew)[i][j];
	}
}

void sweep_partial_b(int i) {
    double (*f)[nx][ny] = (double (*)[nx][ny])f_;
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;

	for (int j = 0; j < ny; j++) {
		if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
			(*unew)[i][j] = (*f)[i][j];
		} else {
			(*unew)[i][j] = 0.25 * ((*u)[i-1][j] + (*u)[i][j+1]
									+ (*u)[i][j-1] + (*u)[i+1][j]
									+ (*f)[i][j] * dx * dy);
		}
	}
}

void FAST_sweep_partial_b(uint64_t swID) {
    double (*f)[nx][ny] = (double (*)[nx][ny])f_;
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;
	unsigned i;

	getOneParameter(swID, i);

	for (int j = 0; j < ny; j++) {
		if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
			(*unew)[i][j] = (*f)[i][j];
		} else {
			(*unew)[i][j] = 0.25 * ((*u)[i-1][j] + (*u)[i][j+1]
									+ (*u)[i][j-1] + (*u)[i+1][j]
									+ (*f)[i][j] * dx * dy);
		}
	}
}

void extra_init() {
	numRetiredTasks = 0;
	femtos_fast_init();

	functionAddresses[0] = (unsigned long long) FAST_sweep_partial_a;
	functionAddresses[1] = (unsigned long long) FAST_sweep_partial_b;
}

void sweep_old_seq(int nx, int ny, double dx, double dy, double *f_,
        int itold, int itnew, double *u_, double *unew_, int block_size)
{
    int i;
    int it;
    int j;
    double (*f)[nx][ny] = (double (*)[nx][ny])f_;
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;

    for (it = itold + 1; it <= itnew; it++) {
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                (*u)[i][j] = (*unew)[i][j];
            }
        }
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
                    (*unew)[i][j] = (*f)[i][j];
                } else {
                    (*unew)[i][j] = 0.25 * ((*u)[i-1][j] + (*u)[i][j+1]
                            + (*u)[i][j-1] + (*u)[i+1][j]
                            + (*f)[i][j] * dx * dy);
                }
            }
        }
    }
}

/* #pragma omp task depend version of SWEEP. */
void sweep (int nx_, int ny_, double dx_, double dy_, double *f__,
        int itold, int itnew, double *u__, double *unew__, int block_size)
{
    int i;
    int it;
    int j;

	unsigned long long num_iterations = 0;
	unsigned numPendingWorkRequests = 0;
	unsigned long long swID = 0;

	nx = nx_;
	ny = ny_;
	dx = dx_;
	dy = dy_;

	f_ = f__;
	u_ = u__;
	unew_ = unew__;
	
    double (*f)[nx][ny] = (double (*)[nx][ny])f_;
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;
	asm volatile ("fence" ::: "memory");

    {
        for (it = itold + 1; it <= itnew; it++) {
            // Save the current estimate.
            for (i = 0; i < nx; i++) {

				/*
				swID = getNewSWID(swID);

				metadataArray[swID].functionAddr = (unsigned long long) FAST_sweep_partial_a;
				metadataArray[swID].depAddresses0[0] = (unsigned long long) (*u)[i];
				metadataArray[swID].depAddresses0[1] = (unsigned long long) (*unew)[i];
				*/
				buildSWIDWithOneParameter(swID, 0, i)

				asm volatile ("fence" ::: "memory");
				num_iterations++;

				make_submission_request_or_work_fast(9, 0, numPendingWorkRequests);
				submit_three_or_work_fast(swID, 2, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*u)[i][0]), 1, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*unew)[i][0]), 0, numPendingWorkRequests);
            }
            // Compute a new estimate.
            for (i = 0; i < nx; i++) {

				/*
				swID = getNewSWID(swID);

				metadataArray[swID].functionAddr = (unsigned long long) FAST_sweep_partial_b;
				metadataArray[swID].depAddresses0[0] = (unsigned long long) i;
				*/
				buildSWIDWithOneParameter(swID, 1, i)

				asm volatile ("fence" ::: "memory");
				num_iterations++;

				make_submission_request_or_work_fast(18, 0, numPendingWorkRequests);
				submit_three_or_work_fast(swID, 5, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*f)[i][0]), 0, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*u)[i-1][0]), 0, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*u)[i][0]), 0, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*u)[i+1][0]), 0, numPendingWorkRequests);
				submit_three_or_work_fast((unsigned long long) &((*unew)[i][0]), 1, numPendingWorkRequests);
            }
        }
        printf("Going to task wait until %d tasks were retired.\n", num_iterations);
		fast_task_wait_and_try_executing_tasks(num_iterations);
    }
}
