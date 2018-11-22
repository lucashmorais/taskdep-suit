#include "interface.h"
#include "rules.h"

#include <omp.h>

/* Min. depth to apply task dispatch */
#define MIN_DEPTH   3
#define FIRST_DEPTH 0

using namespace ttt_rules;

/* Tic tac toe namespace */
namespace ttt_AI {
    int next_move = INVALID;

    int minimax_seq(Board b, Turn p, int depth = 0) {
        /* Is it over? */
        if (game_over(b)) {
            return score(b, depth);
        }

        int size = 0;
        int mov[BOARD_SIZE], scores[BOARD_SIZE];

        /* Switch turn */
        Turn next_turn = (p == PLAYER) ? PC : PLAYER;

        /* Get all the possible movements */
        size = possible_mov(b, mov);

        /* For each of the available movements */
        for (int i = 0; i < size; i++) {
            Board t(b); // temp board

            t.move(mov[i], p);

            /* Since we are smaller than min. depth, despatch into tasks.
             * Otherwise, the amount of computation is too small to
             * require a task parallelism. */

            scores[i] = minimax_seq(t, next_turn, depth + 1);

        }

        /* Get index of highest score */
        int index = (p == PLAYER) ? max_index(scores, size) :
                    min_index(scores, size);

        /* If it is in first depth, updates next_move! */
        if (depth == FIRST_DEPTH) {
            next_move = mov[index];
        }

        /* Return fittest score */
        return scores[index];
    }

    int minimax_omp(Board b, Turn p, int depth = 0) {
        /* Is it over? */
        if (game_over(b)) {
            return score(b, depth);
        }

        int size = 0;
        int mov[BOARD_SIZE], scores[BOARD_SIZE];

        /* Switch turn */
        Turn next_turn = (p == PLAYER) ? PC : PLAYER;

        /* Get all the possible movements */
        size = possible_mov(b, mov);

        /* For each of the available movements */
        for (int i = 0; i < size; i++) {
            Board t(b); // temp board

            t.move(mov[i], p);

            /* Since we are smaller than min. depth, despatch into tasks.
             * Otherwise, the amount of computation is too small to
             * require a task parallelism. */
            if (depth < MIN_DEPTH) {
                #pragma omp task firstprivate(t, next_turn, depth) \
                             shared(scores) depend(out: scores[i])
                scores[i] = minimax_omp(t, next_turn, depth + 1);
            } else {
                scores[i] = minimax_omp(t, next_turn, depth + 1);
            }
        }

        /* If we are at min. depth, wait for all the children
         * tasks to execute. */
        if (depth < MIN_DEPTH) {
            #pragma omp taskwait
        }

        /* Get index of highest score */
        int index = (p == PLAYER) ? max_index(scores, size) :
                    min_index(scores, size);

        /* If it is in first depth, updates next_move! */
        if (depth == FIRST_DEPTH) {
            next_move = mov[index];
        }

        /* Return fittest score */
        return scores[index];
    }
}
