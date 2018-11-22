/*
 * Copyright (C) 2013 Michael Andersch <michael.andersch@mailbox.tu-berlin.de>
 *
 * This file is part of Starbench.
 *
 * Starbench is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Starbench is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Starbench.  If not, see <http://www.gnu.org/licenses/>.
 */

/**********************************************************************************
				INCLUDES & DEFINES
*************************************************************************************/
#include <sys/time.h>
#include <stdlib.h>
#include "rotation_engine.h"
#include "convert_engine.h"
#include "benchmark_engine.h"

#define BAD_EXIT -1;
#define TIME(x) gettimeofday(&x,NULL)

typedef struct timeval timer;
using namespace std;

/**********************************************************************************
				FUNCTION PROTOTYPES
*************************************************************************************/
static long timevaldiff(timer* start, timer* finish);
string* convertToString(char **in, size_t size);
bool parseArgs(string* args, unsigned int &angle, string &inname, string &outname, int &parallel);

/* GLOBAL VARIABLES */
string usage =  "Usage: ./rot-cc <infile> <outfile> <angle> <parallel>\n\n"
                "infile:      input file\n"
                "outfile:     output file\n"
                "angle:       angle to be rotated\n"
                "parallel:       1 if parallel, 0 if sequential\n";
string p_name = "--- StarBENCH - rot-cc Workload ---\n";

timer start, finish, startseq, finishseq;

/*
*	Function: main
*	--------------
*	The program main function.
*/
int main(int argc, char* argv[]) {
    cout << p_name;

    if(argc < 4) {
		cerr << usage;
		return BAD_EXIT;
    }

    string srcfile, destfile;
    unsigned int angle;
    int parallel;
    timer start, finish;
    BenchmarkEngine be, beseq;


    string *args = convertToString(argv, argc);
    if(!parseArgs(args, angle, srcfile, destfile, parallel)) {
        cerr << usage;
        return BAD_EXIT;
    }
	delete [] args;
    if (parallel) {
        printf("Starting parallel code...\n");
        if(!be.init(srcfile, destfile, angle))
            return BAD_EXIT;

    	TIME(start);
        be.run();
    	TIME(finish);

        be.finish();
        cout << "Time: " << (double)timevaldiff(&start, &finish)/1000 << endl;
    } else {
        printf("Starting sequential code...\n");
        if(!beseq.init(srcfile, destfile, angle))
            return BAD_EXIT;

        TIME(startseq);
        beseq.run_seq();
        TIME(finishseq);

        beseq.finish();
        cout << "Time: " << (double)timevaldiff(&startseq, &finishseq)/1000 << endl;
    }



    return 0;
}

/*
*   Function: convertToString
*   -------------------------
*   Converts the c-string program arguments into c++-strings and returns
*   a pointer to an array of such strings.
*/
string* convertToString(char** in, size_t size) {
    string* args = new string[size];
    for(size_t i = 0; i < size; i++) {
       args[i] = in[i];
    }
    return args;
}

/*
*   Function: parseArgs
*   -------------------
*   Extracts the rotation angle as well as the in- and output file names
*   from the string array args, storing them in the specified variables.
*/
bool parseArgs(string* args, unsigned int &angle, string &inname, string &outname, int &parallel) {
    const char *tmp = args[3].c_str();
    angle = atoi(tmp) % 360;
    const char *tmp2 = args[4].c_str();
    parallel = atoi(tmp2) % 360;

    if (parallel < 0 || parallel > 1) {
        cerr << "Parallel flag must be 0 or 1" << endl;
        exit(-1);
    }

    if (angle < 0) {
        cerr << "Bad arguments, exiting" << endl;
        exit(-1);
    }

    inname = args[1];
    outname = args[2];
    return true;
}

/*
*   Function: timevaldiff
*   ---------------------
*   Provides a millisecond-resolution timer, computing the elapsed time
*   in between the two given timeval structures.
*/
static long timevaldiff(timer* start, timer* finish){
	long msec;
	msec = (finish->tv_sec - start->tv_sec)*1000;
	msec += (finish->tv_usec - start->tv_usec)/1000;
	return msec;
}
