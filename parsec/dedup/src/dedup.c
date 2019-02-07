#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"
#include "debug.h"
#include "dedupdef.h"
#include "encoder.h"
#include "decoder.h"
#include "config.h"
#include "queue.h"

#include "../../../c/bench.h"

#ifdef ENABLE_DMALLOC
#include <dmalloc.h>
#endif //ENABLE_DMALLOC

#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif //ENABLE_PTHREADS

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif //ENABLE_PARSEC_HOOKS

int allow_out;

config_t * conf;


/*--------------------------------------------------------------------------*/
static void
usage(char* prog)
{
  if(allow_out) printf("usage: %s [-cusfvh] [-w gzip/bzip2/none] [-i file] [-o file] [-t number_of_threads]\n",prog);
  if(allow_out) printf("-c \t\t\tcompress\n");
  if(allow_out) printf("-u \t\t\tuncompress\n");
  if(allow_out) printf("-p \t\t\tpreloading (for benchmarking purposes)\n");
  if(allow_out) printf("-w \t\t\tcompression type: gzip/bzip2/none\n");
  if(allow_out) printf("-i file\t\t\tthe input file\n");
  if(allow_out) printf("-o file\t\t\tthe output file\n");
  if(allow_out) printf("-t \t\t\tnumber of threads per stage \n");
  if(allow_out) printf("-v \t\t\tverbose output\n");
  if(allow_out) printf("-h \t\t\thelp\n");
}
/*--------------------------------------------------------------------------*/
int main(int argc, char** argv) {

  allow_out = 1;
  if(getenv("BENCH_SILENT") != NULL) allow_out = 0;

  process_name("parsec-dedup");
  process_mode(SEQ);
  process_args(argc, argv);
  process_init();

#ifdef ENABLE_OMPSS
  process_mode(OMPSS);
#endif

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
  if(allow_out) printf("PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION)"\n");
#else
  if(allow_out) printf("PARSEC Benchmark Suite\n");
#endif //PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
        __parsec_bench_begin(__parsec_dedup);
#endif //ENABLE_PARSEC_HOOKS

#if defined(ENABLE_OMPSS) || defined(ENABLE_OMP)
  if(allow_out) printf("Warning! Argument -t must be 1, to set number of threads use NX_ARGS for OmpSs or OMP_NUM_THREAD for OpenMP 4.0\n");
  task_init_measure();
#endif

  int32 compress = TRUE;

  //We force the sha1 sum to be integer-aligned, check that the length of a sha1 sum is a multiple of unsigned int
  assert(SHA1_LEN % sizeof(unsigned int) == 0);

  conf = (config_t *) malloc(sizeof(config_t));
  if (conf == NULL) {
    EXIT_TRACE("Memory allocation failed\n");
  }

  strcpy(conf->outfile, "");
  conf->compress_type = COMPRESS_GZIP;
  conf->preloading = 0;
  conf->nthreads = 1;
  conf->verbose = 0;

  //parse the args
  int ch;
  opterr = 0;
  optind = 1;
  while (-1 != (ch = getopt(argc, argv, "cupvo:i:w:t:h"))) {
    switch (ch) {
    case 'c':
      compress = TRUE;
      strcpy(conf->infile, "test.txt");
      strcpy(conf->outfile, "out.ddp");
      break;
    case 'u':
      compress = FALSE;
      strcpy(conf->infile, "out.ddp");
      strcpy(conf->outfile, "new.txt");
      break;
    case 'w':
      if (strcmp(optarg, "gzip") == 0)
        conf->compress_type = COMPRESS_GZIP;
      else if (strcmp(optarg, "bzip2") == 0) 
        conf->compress_type = COMPRESS_BZIP2;
      else if (strcmp(optarg, "none") == 0)
        conf->compress_type = COMPRESS_NONE;
      else {
        fprintf(stdout, "Unknown compression type `%s'.\n", optarg);
        usage(argv[0]);
        return -1;
      }
      break;
    case 'o':
      strcpy(conf->outfile, optarg);
      break;
    case 'i':
      strcpy(conf->infile, optarg);
      break;
    case 'h':
      usage(argv[0]);
      return -1;
    case 'p':
      conf->preloading = TRUE;
      break;
    case 't':
      conf->nthreads = atoi(optarg);
      break;
    case 'v':
      conf->verbose = TRUE;
      break;
    case '?':
      fprintf(stdout, "Unknown option `-%c'.\n", optopt);
      usage(argv[0]);
      return -1;
    }
  }

#ifndef ENABLE_BZIP2_COMPRESSION
 if (conf->compress_type == COMPRESS_BZIP2){
    if(allow_out) printf("Bzip2 compression not supported\n");
    exit(1);
  }
#endif

#ifndef ENABLE_GZIP_COMPRESSION
 if (conf->compress_type == COMPRESS_GZIP){
    if(allow_out) printf("Gzip compression not supported\n");
    exit(1);
  }
#endif

#ifndef ENABLE_STATISTICS
 if (conf->verbose){
    if(allow_out) printf("Statistics collection not supported\n");
    exit(1);
  }
#endif

#ifndef ENABLE_PTHREADS
#ifndef ENABLE_OMPSS
 if (conf->nthreads != 1){
    if(allow_out) printf("Number of threads must be 1 (serial version)\n");
    exit(1);
  }
#endif
#endif

  //process_start_measure();
  if (compress) {
    Encode(conf);
  } else {
    Decode(conf);
  }
  //process_stop_measure();

  free(conf);

#ifdef ENABLE_PARSEC_HOOKS
  __parsec_bench_end();
#endif

  process_append_file(conf->outfile);
  dump_csv(stdout);

  return 0;
}

