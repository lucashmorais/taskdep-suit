#include "interface.h"
#include "rules.h"
#include "AI.h"

#include <cstdio>

#include <sys/time.h>
#include <time.h>

using namespace ttt_rules;

void solve_seq(Board a, Turn cur_t) {
/* Check who starts */
    {
        bool finished = false;

        while (!finished) {
            char input[10];

            printf("Would you like to start? [Y/N] \n");

            int r = scanf("%s", input);

            if (input[0] == 'Y' || input[0] == 'y' && r == 1) {
                cur_t = PLAYER;

                finished = true;
            } else if (input[0] == 'N' || input[0] == 'n' && r == 1) {
                cur_t = PC;

                finished = true;
            }
        }
    }

     /* While game isn't over */
    while(!game_over(a)) {
        if (cur_t == PLAYER) {
            bool finished = false; // in order to proceed

            int pos;
            char input[10]; // player input

            while (!finished) {
                a.display_board();

                printf("Please, choose your movement:\n");

                int r = scanf("%3s", input);

                pos = a.read_input(input);

                if (pos != INVALID && r == 1) {
                    a.move(pos, cur_t);

                    finished = true;
                    cur_t = PC;
                } else {
                    printf("Invalid move!\n\n");
                }
            }
        } else {
            a.display_board();
            printf("It's computer's turn:\n\n");

            ttt_AI::minimax_seq(a, cur_t);

            a.move(ttt_AI::next_move, cur_t);

            cur_t = PLAYER;
        }
    }

    /* Finally, get result */
    int result = score(a);

    a.display_board();

    if (result > 0) {
        printf("Congratulations, you won!\n");
    } else if (result == 0) {
        printf("It's a tie!\n");
    } else {
        printf("Computer won, muahaha!\n");
    }
}

void solve_OMP(Board a, Turn cur_t) {
/* Check who starts */
    {
        bool finished = false;

        while (!finished) {
            char input[10];

            printf("Would you like to start? [Y/N] \n");

            int r = scanf("%s", input);

            if (input[0] == 'Y' || input[0] == 'y' && r == 1) {
                cur_t = PLAYER;

                finished = true;
            } else if (input[0] == 'N' || input[0] == 'n' && r == 1) {
                cur_t = PC;

                finished = true;
            }
        }
    }

     /* While game isn't over */
    while(!game_over(a)) {
        if (cur_t == PLAYER) {
            bool finished = false; // in order to proceed

            int pos;
            char input[10]; // player input

            while (!finished) {
                a.display_board();

                printf("Please, choose your movement:\n");

                int r = scanf("%3s", input);

                pos = a.read_input(input);

                if (pos != INVALID && r == 1) {
                    a.move(pos, cur_t);

                    finished = true;
                    cur_t = PC;
                } else {
                    printf("Invalid move!\n\n");
                }
            }
        } else {
            a.display_board();
            printf("It's computer's turn:\n\n");

            #pragma omp parallel
            #pragma omp single
            ttt_AI::minimax_omp(a, cur_t);

            a.move(ttt_AI::next_move, cur_t);

            cur_t = PLAYER;
        }
    }

    /* Finally, get result */
    int result = score(a);

    a.display_board();

    if (result > 0) {
        printf("Congratulations, you won!\n");
    } else if (result == 0) {
        printf("It's a tie!\n");
    } else {
        printf("Computer won, muahaha!\n");
    }
}


int main () {
    Board a;
    Turn cur_t;

    /* Set time measure */
    struct timeval start, stop;
    unsigned long elapsed_parallel, elapsed_seq;
    gettimeofday(&start, NULL);
    solve_OMP(a, cur_t);
    /* Get elapsed time */
    gettimeofday(&stop, NULL);
    elapsed_parallel = 1000000 * (stop.tv_sec - start.tv_sec);
    elapsed_parallel += stop.tv_usec - start.tv_usec;

    gettimeofday(&start, NULL);
    solve_seq(a, cur_t);
    /* Get elapsed time */
    gettimeofday(&stop, NULL);
    elapsed_seq = 1000000 * (stop.tv_sec - start.tv_sec);
    elapsed_seq += stop.tv_usec - start.tv_usec;

    printf("Parallel Time: %lu\n", elapsed_parallel);
    printf("Sequential Time: %lu\n", elapsed_seq);

    return 0;
}
