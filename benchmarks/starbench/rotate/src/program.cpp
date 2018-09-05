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

#define BAD_EXIT -1;
#define TIME(x) gettimeofday(&x,NULL)

typedef struct timeval timer;
using namespace std;

/**********************************************************************************
				FUNCTION PROTOTYPES
*************************************************************************************/
static long timevaldiff(timer* start, timer* finish);
bool parseArgs(char* argv[], int argc, unsigned int &angle, string &inname, string &outname, int &flag);

/* GLOBAL VARIABLES */
string usage =  "Usage: ./rot <infile> <outfile> <angle> <flag>\n\n"
                "infile:      input file\n"
                "outfile:     output file\n"
                "angle:       angle to be rotated\n"
                "flag:       0 for sequential, 1 for parallel\n";
string p_name = "--- StarBENCH - rotate Kernel ---\n";

timer a,b,c,d;

/*
*	Function: main
*	--------------
*	The program main function.
*/
int main(int argc, char* argv[]) {
    cout << p_name;

    if(argc != 5) {
            printf("ARGC %d\n", argc);
		cerr << usage;
		return BAD_EXIT;
    }

    string srcfile, destfile;
    unsigned int angle;
    int flag;
    timer start, finish, startseq, finishseq;
    RotateEngine re, reseq;

    if(!parseArgs(argv, argc, angle, srcfile, destfile, flag)) {
        cerr << usage;
        return BAD_EXIT;
    }
    if (flag) {
        printf("Starting parallel code...\n");
        if(!re.init(srcfile, destfile, angle)) return BAD_EXIT;

    	//re.printRotationState();

    	TIME(start);
        re.run();
    	TIME(finish);

        re.finish();
        cout << "Conversion timing (parallel):" << endl;
        cout << "Compute: " << (double)timevaldiff(&a, &b)/1000 << "s" << endl;
        cout << "Time: " << (double)timevaldiff(&start, &finish)/1000 << endl;
    } else {
        printf("Starting sequential code...\n");
        if(!reseq.init(srcfile, destfile, angle)) return BAD_EXIT;

        //re.printRotationState();

        TIME(startseq);
        reseq.run_seq();
        TIME(finishseq);

        reseq.finish();
        cout << "Conversion timing (sequential):" << endl;
        cout << "Compute: " << (double)timevaldiff(&c, &d)/1000 << "s" << endl;
        cout << "Time: " << (double)timevaldiff(&startseq, &finishseq)/1000 << endl;
    }
    return 0;
}

/*
*   Function: parseArgs
*   -------------------
*   Extracts the rotation angle as well as the in- and output file names
*   from the string array args, storing them in the specified variables.
*/
bool parseArgs(char* argv[], int argc, unsigned int &angle, string &inname, string &outname, int &flag) {
    if (argc !=5)
        return false;

    const char *tmp = argv[3];
    angle = atoi(tmp) % 360;

    const char *tmp2 = argv[4];
    flag = atoi(tmp2);

    if (flag < 0 || flag > 1) {
        cerr << "Flag must be 0 or 1" << endl;
        exit(-1);
    }
    if (angle < 0) {
        cerr << "Bad arguments, exiting" << endl;
        exit(-1);
    }

    inname = argv[1];
    outname = argv[2];
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
