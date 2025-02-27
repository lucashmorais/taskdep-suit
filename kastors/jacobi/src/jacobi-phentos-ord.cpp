# include "poisson.h"
#include <iostream>
#include "../../../c/bench.h"
#include "../../../../phentos/phentos/femtos.hpp"

double *f_;
double *u_;
double *unew_;

int nx, ny;
double dx, dy;

void extra_init() {
	numRetiredTasks = 0;
	totalNumberOfSubmittedTasks = 0;
#ifdef ZERO_PACKETS_V2
	std::cout << "[feature]: ZERO_PACKETS_V2" << std::endl;
	std::cout << "[feature]: SW-BASED PADDING" << std::endl;
	std::cout << "[feature]: RESET COUNT OF SUBMITTED TASKS" << std::endl;
#endif
	femtos_init();
}

void sweep_partial_a(double *u_array, double *unew_array) {
	for (int j = 0; j < ny; j++) {
		u_array[j] = unew_array[j];
	}
}

void ORD_sweep_partial_a(uint64_t swID) {
	double * u_array = (double *) metadataArray[swID].depAddresses0[0];
	double * unew_array = (double *) metadataArray[swID].depAddresses0[1];

	for (int j = 0; j < ny; j++) {
		u_array[j] = unew_array[j];
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

void ORD_sweep_partial_b(uint64_t swID) {
    double (*f)[nx][ny] = (double (*)[nx][ny])f_;
    double (*u)[nx][ny] = (double (*)[nx][ny])u_;
    double (*unew)[nx][ny] = (double (*)[nx][ny])unew_;

	int i = (int) metadataArray[swID].depAddresses0[0];

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
				swID = getNewSWID(swID);

				metadataArray[swID].functionAddr = (unsigned long long) ORD_sweep_partial_a;
				metadataArray[swID].depAddresses0[0] = (unsigned long long) (*u)[i];
				metadataArray[swID].depAddresses0[1] = (unsigned long long) (*unew)[i];

				asm volatile ("fence" ::: "memory");
				num_iterations++;

#ifdef ZERO_PACKETS_V2
				make_submission_request_or_work(48, 0, numPendingWorkRequests);
				submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
				make_submission_request_or_work(9, 0, numPendingWorkRequests);
				submit_three_or_work(swID, 2, numPendingWorkRequests);
#endif
				submit_three_or_work((unsigned long long) &((*u)[i][0]), 1, numPendingWorkRequests);
				submit_three_or_work((unsigned long long) &((*unew)[i][0]), 0, numPendingWorkRequests);

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
            // Compute a new estimate.
            for (i = 0; i < nx; i++) {
				swID = getNewSWID(swID);

				metadataArray[swID].functionAddr = (unsigned long long) ORD_sweep_partial_b;
				metadataArray[swID].depAddresses0[0] = (unsigned long long) i;

				asm volatile ("fence" ::: "memory");
				num_iterations++;

#ifdef ZERO_PACKETS_V2
				make_submission_request_or_work(48, 0, numPendingWorkRequests);
				submit_three_or_work(swID, 15, numPendingWorkRequests);
#else
				make_submission_request_or_work(18, 0, numPendingWorkRequests);
				submit_three_or_work(swID, 5, numPendingWorkRequests);
#endif
				submit_three_or_work((unsigned long long) &((*f)[i][0]), 0, numPendingWorkRequests);
				submit_three_or_work((unsigned long long) &((*u)[i-1][0]), 0, numPendingWorkRequests);
				submit_three_or_work((unsigned long long) &((*u)[i][0]), 0, numPendingWorkRequests);
				submit_three_or_work((unsigned long long) &((*u)[i+1][0]), 0, numPendingWorkRequests);
				submit_three_or_work((unsigned long long) &((*unew)[i][0]), 1, numPendingWorkRequests);

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
#endif
            }
        }
        printf("Going to task wait until %d tasks were retired.\n", num_iterations);
		task_wait_and_try_executing_tasks(num_iterations);
    }
}
