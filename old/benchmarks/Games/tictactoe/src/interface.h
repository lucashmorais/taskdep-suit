/*
 * interface.h
 *  Responsible for dealing with the interface of tic tac toe game,
 * as well the board structure.
 * */

#pragma once

#include <iostream>     // cin, cout and cerr
#include <string.h>     // strlen

#define INPUT_SIZE 2

#define BOARD_SIZE 9
#define ROW        3

/* Boolean */
#define bool  int
#define false 0
#define true  1

#define INVALID -1

/* Each position corresponds to its ASCII value */
typedef enum Position { EMPTY = 32, O = 111, X = 120 } Position;

/* Each of the possible turns 
 *   PLAYER = 'o'   
 *   PC     = 'x'
 * */
typedef enum Turn { PLAYER = 111, PC = 120 } Turn;

class Board {
public:
    Position p[BOARD_SIZE]; // each position of the board

    /* Default constructor */
    Board() {
        /* Initializes board */
        for (int i = 0; i < BOARD_SIZE; i++) {
            p[i] = EMPTY;
        }
    }

    /* Initialize based on a board b */
    Board (const Board& b) {
        /* Duplicates board */
        for (int i = 0; i < BOARD_SIZE; i++) {
            p[i] = b.p[i];
        }
    }

    /* Displays board at stdout */
    void display_board() {
        char row = 'A';

        for (int i = 0; i < ROW; i++) {
            std::cout << "    " << row << "  ";

            for (int j = 0; j < ROW; j++) {
                std::cout << (char)p[i * ROW + j];

                /* If is not the last of the column */
                if (j + 1 < ROW) {
                    std::cout << " | ";
                } else {
                    std::cout << "\n";
                }
            }

            row++; // increases row count
        }

        std::cout << "\n       1   2   3\n"; // column count
    }

    /* Read move input and translate it into a position
     *   return: where the new position should be in the board,
     *           otherwise, return INVALID
     * */
    int read_input(const char* input) {
        size_t size = strlen(input);
        int i, j;

        /* Check input */
        if (size >= INPUT_SIZE) {
            /* Check row */
            if (input[0] >= 'A' && input[0] <= 'C') {
                i = (int)input[0] - (int)'A';
            } else if (input[0] >= 'a' && input[0] <= 'c') {
                i = (int)input[0] - (int)'a';
            } else {
                return INVALID;
            }

            /* Check column */
            if (input[1] >= '1' && input[1] <= '3') {
                j = input[1] - (int)'1';
            } else {
                return INVALID;
            }
        } else {
            return INVALID;
        }

        /* Finally, check if it is an empty position */
        if (p[i * ROW + j] != EMPTY) {
            return INVALID;
        } else {
            return i * ROW + j;
        }
    }

    /* Make a move based on position and the type of player
     *   return: if it was a valid move 
     * */
    bool move(int position, Turn move) {
        if (p[position] == EMPTY) {
            p[position] = (Position)move;
        } else {
            return false;
        }

        return true;
    }
};
