#include "interface.h"

#pragma once

/* Definitions for a winning game, given a board b */
#define DIAGONAL_LEFT     b.p[0] != EMPTY && b.p[0] == b.p[4] && b.p[4] == b.p[8]
#define DIAGONAL_RIGHT    b.p[2] != EMPTY && b.p[2] == b.p[4] && b.p[4] == b.p[6]

#define VERTICAL_LEFT     b.p[0] != EMPTY && b.p[0] == b.p[3] && b.p[3] == b.p[6]
#define VERTICAL_CENTER   b.p[1] != EMPTY && b.p[1] == b.p[4] && b.p[4] == b.p[7]
#define VERTICAL_RIGHT    b.p[2] != EMPTY && b.p[2] == b.p[5] && b.p[5] == b.p[8]

#define HORIZONTAL_UP     b.p[0] != EMPTY && b.p[0] == b.p[1] && b.p[1] == b.p[2]
#define HORIZONTAL_CENTER b.p[3] != EMPTY && b.p[3] == b.p[4] && b.p[4] == b.p[5]
#define HORIZONTAL_DOWN   b.p[6] != EMPTY && b.p[6] == b.p[7] && b.p[7] == b.p[8]

#define MAX_VALUE         10000

/* Default score for minimax */
#define DEFAULT_SCORE 10

namespace ttt_rules {
    /* Functions signatures */
    bool game_over(Board b);
    bool full(Board b);
    int score(Board b, int);
    int possible_mov(Board b, int mov[]);

    /* Is the game over already? */
    bool game_over(Board b) {
        if (DIAGONAL_LEFT || DIAGONAL_RIGHT || VERTICAL_LEFT || 
            VERTICAL_CENTER || VERTICAL_RIGHT || HORIZONTAL_UP || 
            HORIZONTAL_CENTER || HORIZONTAL_DOWN || full(b)) {
            return true;
        } else {
            return false;
        }
    }

    /* Returns if the given board is full */
    bool full(Board b) {
        for (int i = 0; i < BOARD_SIZE; i++) {
            if (b.p[i] == EMPTY) {
                return false;
            }
        }

        return true;
    }

    /* Finds out who won the game
     *   depth: level of depth at minimax
     *   returns: score of the round, based on the winner
     * */
    int score(Board b, int depth = 0) {
        Turn winner;

        if (DIAGONAL_LEFT || DIAGONAL_RIGHT || VERTICAL_CENTER || 
            HORIZONTAL_CENTER) {
            winner = (Turn)b.p[4];
        } else if (VERTICAL_LEFT || HORIZONTAL_UP) {
            winner = (Turn)b.p[0];
        } else if (VERTICAL_RIGHT || HORIZONTAL_DOWN) {
            winner = (Turn)b.p[8];
        } else {
            return 0;
        }

        if (winner == PLAYER) {
            return DEFAULT_SCORE - depth;
        } else {
            return depth - DEFAULT_SCORE;
        }
    }

    /* Check all the available movements on present board,
     * storing each of them at mov[]
     *   return: total number of movements 
     * */
    int possible_mov(Board b, int mov[]) {
        int size = 0;

        for (int i = 0; i < BOARD_SIZE; i++) {
            if (b.p[i] == EMPTY) {
                mov[size++] = i;
            }
        }

        return size;
    }

    /* Return index of max. value */
    int max_index(int v[], int size) {
        int max_v = -MAX_VALUE, max_i = 0;

        for (int i = 0; i < size; i++) {
            /* Find max. value */
            if (v[i] > max_v) {
                max_v = v[i];
                max_i = i;
            }
        }

        return max_i;
    }

    /* Return index of min. value */
    int min_index(int v[], int size) {
        int min_v = MAX_VALUE, min_i = 0;

        for (int i = 0; i < size; i++) {
            /* Find min. value */
            if (v[i] < min_v) {
                min_v = v[i];
                min_i = i;
            }
        }

        return min_i;
    }
}
