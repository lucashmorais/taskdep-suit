#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "../../../c/bench.h"
int allow_out;

#ifdef _OMPSS
    #include <omp.h>
#elif _OPENMP
    #include <omp.h>
#endif

#include "main.h"

#define min(a, b) ((a<b)?a:b)
#define max(a, b) ((a>b)?a:b)

void parse(int argc, char* argv[], struct user_parameters* params)
{
    process_init();
    process_name("karstors-jacobi");
    process_args(argc, argv);
#ifdef _OMPSS
    process_mode(OMPSS);
    task_init_measure();
#elif _OPENMP
    process_mode(OPENMP_TASK);
    task_init_measure();
#else
    process_mode(SEQ);
#endif

    int i;
    params->check = 0;
    for(i=1; i<argc; i++) {
        if(!strcmp(argv[i], "-c"))
            params->check = 1;
        else if(!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) {

            if(allow_out) printf("----------------------------------------------\n");
            if(allow_out) printf("-                KaStORS                     -\n");
            if(allow_out) printf("-   Kaapi Starpu OpenMP Runtime task Suite   -\n");
            if(allow_out) printf("----------------------------------------------\n");
            if(allow_out) printf("-h, --help : Show help information\n");
            if(allow_out) printf("-c : Ask to check result\n");
            if(allow_out) printf("-i : Number of iterations\n");
#ifdef TITER
            if(allow_out) printf("-r : Number ot timestep iteration\n");
#endif
#ifdef MSIZE
            if(allow_out) printf("-n : Matrix size\n");
#endif
#ifdef SMSIZE
            if(allow_out) printf("-m : SubMatrix size\n");
#endif
#ifdef BSIZE
            if(allow_out) printf("-b : Block size\n");
#endif
#ifdef IBSIZE
            if(allow_out) printf("-ib : Internal Block size\n");
#endif
#ifdef CUTOFF_SIZE
            if(allow_out) printf("-s : Cutoff (Size of the matrix)\n");
#endif
#ifdef CUTOFF_DEPTH
            if(allow_out) printf("-d : Cutoff (depth)\n");
#endif
            exit(EXIT_SUCCESS);

        } else if(!strcmp(argv[i], "-i")) {
            if (++i < argc)
                params->niter = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-i requires a number\n");
                exit(EXIT_FAILURE);
            }
#ifdef TITER
        } else if(!strcmp(argv[i], "-r")) {
            if (++i < argc)
                params->titer = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-r requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef MSIZE
        } else if(!strcmp(argv[i], "-n")) {
            if (++i < argc)
                params->matrix_size = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-n requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef SMSIZE
        } else if(!strcmp(argv[i], "-m")) {
            if (++i < argc)
                params->submatrix_size = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-m requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef BSIZE
        } else if(!strcmp(argv[i], "-b")) {
            if (++i < argc)
                params->blocksize = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-b requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef IBSIZE
        } else if(!strcmp(argv[i], "-ib")) {
            if (++i < argc)
                params->iblocksize = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-ib requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef CUTOFF_SIZE
        } else if(!strcmp(argv[i], "-s")) {
            if (++i < argc)
                params->cutoff_size = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-s requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
#ifdef CUTOFF_DEPTH
        } else if(!strcmp(argv[i], "-d")) {
            if (++i < argc)
                params->cutoff_depth = atoi(argv[i]);
            else {
                if(allow_out) fprintf(stderr, "-d requires a number\n");
                exit(EXIT_FAILURE);
            }
#endif
        } else
            if(allow_out) fprintf(stderr, "Unknown parameter : %s\n", argv[i]);
    }
}

int comp (const void * elem1, const void * elem2) 
{
    double f = *((double*)elem1);
    double s = *((double*)elem2);
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;
}

int main(int argc, char* argv[])
{
    allow_out = 1;
    if(getenv("BENCH_SILENT") != NULL) allow_out = 0;

    int num_threads = 1;
    struct user_parameters params;
    memset(&params, 0, sizeof(params));

    /* default value */
    params.niter = 1;

    parse(argc, argv, &params);

// get Number of thread if OpenMP is activated
#ifdef _OMPSS
    num_threads = omp_get_num_threads();
#elif _OPENMP
    #pragma omp parallel
    #pragma omp master
    num_threads = omp_get_num_threads();
#endif

    // warmup
    run(&params);

    double mean = 0.0;
    double meansqr = 0.0;
    double min_ = DBL_MAX;
    double max_ = -1;
    double* all_times = (double*)malloc(sizeof(double) * params.niter); 

    for (int i=0; i<params.niter; ++i)
    {
      double cur_time = run(&params);
      all_times[i] = cur_time;
      mean += cur_time;
      min_ = min(min_, cur_time);
      max_ = max(max_, cur_time);
      meansqr += cur_time * cur_time;
      }
    mean /= params.niter;
    meansqr /= params.niter;
    double stddev = sqrt(meansqr - mean * mean);

    qsort(all_times, params.niter, sizeof(double), comp);
    double median = all_times[params.niter / 2];

    free(all_times);

    if(allow_out) printf("Program : %s\n", argv[0]);
#ifdef MSIZE
    if(allow_out) printf("Size : %d\n", params.matrix_size);
#endif
#ifdef SMSIZE
    if(allow_out) printf("Submatrix size : %d\n", params.submatrix_size);
#endif
#ifdef BSIZE
    if(allow_out) printf("Blocksize : %d\n", params.blocksize);
#endif
#ifdef IBSIZE
    if(allow_out) printf("Internal Blocksize : %d\n", params.iblocksize);
#endif
#ifdef TITER
    if(allow_out) printf("Iteration time : %d\n", params.titer);
#endif
    if(allow_out) printf("Iterations : %d\n", params.niter);
#ifdef CUTOFF_SIZE
    if(allow_out) printf("Cutoff Size : %d\n", params.cutoff_size);
#endif
#ifdef CUTOFF_DEPTH
    if(allow_out) printf("Cutoff depth : %d\n", params.cutoff_depth);
#endif
    if(allow_out) printf("Threads : %d\n", num_threads);
#ifdef GFLOPS
    if(allow_out) printf("Gflops:: ");
#else
    if(allow_out) printf("Time(sec):: ");
#endif
    if(allow_out) printf("avg : %lf :: std : %lf :: min : %lf :: max : %lf :: median : %lf\n",
           mean, stddev, min_, max_, median);
    if(params.check)
        if(allow_out) printf("Check : %s\n", (params.succeed)?
                ((params.succeed > 1)?"not implemented":"success")
                :"fail");
    if (params.string2display !=0)
      if(allow_out) printf("%s", params.string2display);
    if(allow_out) printf("\n");

    if (params.succeed)
	process_append_result("Sucess!", 7);
    else
        process_append_result("Fail", 4);
    dump_csv(stdout);
    return 0;
}
