#include <stdio.h>
#include <stdlib.h>
#include "../../common/Utils.h"

void srandom (unsigned seed);
double dboard (int darts, int *pscore);

#define MAX_TASKS    50

int score[MAX_TASKS];

int main (int argc, char *argv[])
{
    int i, n_throws, n_tasks, n_darts;
    int sum_score = 0, sum_darts = 0;
    double pi = 0.0;
    double t_start, t_end;
    if(argc != 3) {
        printf("usage: pi_calc num_throws num_tasks\n");
        exit(1);
    }

    n_throws = atoi(argv[1]);
    n_tasks = atoi(argv[2]);
    n_darts = n_throws / n_tasks;
    if( n_tasks > MAX_TASKS ) {
        printf("too many tasks! maximum tasks will be: %d\n", MAX_TASKS);
        exit(1);
    }


    // Note:
    // - need "parallel" to create multiple worker threads
    // - need "single" to ensure only 1 thread to generate the task
    // - need "task" to generate tasks with delayed execution
    // - need "taskwait": shouldn't there be an implicit barrier at end-of-paralel block?

    t_start = rtclock();
    #pragma omp parallel
    {
    #pragma omp single
    for( i = 0; i < n_tasks; i++) {
        int *pscore = &score[i];
        #pragma omp task firstprivate(i, n_darts, pscore)
        {
            pi = dboard(n_darts, &score[i]);
            printf("%d-th task: pi=%f\n", i, pi);
        }
    }
    }
    #pragma omp taskwait

    sum_score = 0;
    for( i = 0; i < n_tasks; i++ ) {
        {
        sum_score += score[i];
        }
    }
    t_end = rtclock();
    FILE *file;

    sum_darts = n_darts * n_tasks;
    pi = 4.0 * sum_score / sum_darts;
    printf("combine all tasks: pi=%f\n", pi);
    printf("Time: %0.6lf\n", t_end - t_start);
    fprintf(stderr, "Time: %0.6lf\n", t_end - t_start);

    return 0;
}



