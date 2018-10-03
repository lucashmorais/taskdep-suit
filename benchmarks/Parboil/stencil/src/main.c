
/***************************************************************************
 *cr
 *cr            (C) Copyright 2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
#include "../../common/parboil.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "../../common/polybenchUtilFuncts.h"

#include <endian.h>
#include <stdlib.h>
#include <malloc.h>
#include <inttypes.h>

#define Index3D(_nx,_ny,_i,_j,_k) ((_i)+_nx*((_j)+_ny*(_k)))

void
outputData(char* fName, float *h_A0,int nx,int ny,int nz) {
	FILE* fid = fopen(fName, "w");
	uint32_t tmp32;
	if(fid == NULL) {
		fprintf(stderr, "Cannot open output file\n");
		exit(-1);
	}
	tmp32 = nx*ny*nz;
	fwrite(&tmp32, sizeof(uint32_t), 1, fid);
	fwrite(h_A0, sizeof(float), tmp32, fid);

	fclose (fid);
}


void cpu_stencilCPU(float c0,float c1, float *A0,float * Anext,const int nx, const int ny, const int nz)
{
	#pragma omp parallel
	#pragma omp single
	{
		for(int k=1;k<nz-1;k++) {
			#pragma omp task
			for(int j=1;j<ny-1;j++) {
				for(int i=1;i<nx-1;i++) {
					Anext[Index3D (nx, ny, i, j, k)] = 
					(A0[Index3D (nx, ny, i, j, k + 1)] +
					A0[Index3D (nx, ny, i, j, k - 1)] +
					A0[Index3D (nx, ny, i, j + 1, k)] +
					A0[Index3D (nx, ny, i, j - 1, k)] +
					A0[Index3D (nx, ny, i + 1, j, k)] +
					A0[Index3D (nx, ny, i - 1, j, k)])*c1
					- A0[Index3D (nx, ny, i, j, k)]*c0;
				}
			}
		}
	}
}


#define ERROR_THRESHOLD 0.05
#define GPU_DEVICE 1
double t_start, t_end, t_start_GPU, t_end_GPU;

float *h_Anext_CPU;
int NX, NY, NZ;

typedef float DATA_TYPE;
struct pb_Parameters *parameters;

static int
read_data(float *A0, int nx,int ny,int nz,FILE *fp) 
{	
	int s=0;
	int i, j, k;
	for(i=0;i<nz;i++) {
		for(j=0;j<ny;j++) {
			for(k=0;k<nx;k++) {
				fread(A0+s,sizeof(float),1,fp);
				s++;
			}
		}
	}
	return 0;
}

double
stencilCPU(int argc, char** argv) {
	int nx,ny,nz;
	int size;
	int iteration;
	float c0=1.0f/6.0f;
	float c1=1.0f/6.0f/6.0f;

	if (argc<5) 
	{
		printf("Usage: probe nx ny nz tx ty t\n"
			"nx: the grid size x\n"
			"ny: the grid size y\n"
			"nz: the grid size z\n"
			"t: the iteration time\n");
		return -1;
	}

	nx = atoi(argv[1]);
	if (nx<1)
		return -1;
	ny = atoi(argv[2]);
	if (ny<1)
		return -1;
	nz = atoi(argv[3]);
	if (nz<1)
		return -1;
	iteration = atoi(argv[4]);
	if(iteration<1)
		return -1;


	//host data
	float *h_A0;
	float *h_Anext;

	size=nx*ny*nz;

	h_A0=(float*)malloc(sizeof(float)*size);
	h_Anext=(float*)malloc(sizeof(float)*size);
	FILE *fp = fopen(parameters->inpFiles[0], "rb");
	read_data(h_A0, nx,ny,nz,fp);
	fclose(fp);
	memcpy (h_Anext,h_A0 ,sizeof(float)*size);

	int t;
	t_start = rtclock();
	for(t=0;t<iteration;t++) {
		cpu_stencilCPU(c0,c1, h_A0, h_Anext, nx, ny,  nz);
		float *temp=h_A0;
		h_A0 = h_Anext;
		h_Anext = temp;
	}
	t_end = rtclock();

	float *temp=h_A0;
	h_A0 = h_Anext;
	h_Anext_CPU = temp;
	NX = nx;
	NY = ny;
	NZ = nz;

	free (h_A0);

	return t_end - t_start;
}

int
main(int argc, char** argv) {
	double t_CPU;

	parameters = pb_ReadParameters(&argc, argv);

	t_CPU = stencilCPU(argc, argv);
	fprintf(stdout, "CPU Runtime: %0.6lfs\n", t_CPU);

	pb_FreeParameters(parameters);
	free (h_Anext_CPU);

	return 0;
}
