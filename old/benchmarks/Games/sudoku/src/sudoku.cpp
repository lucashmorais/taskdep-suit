/* Sudoku application
 * Code taken at:
 *                http://stackoverflow.com/a/26195266
 */

#include <iostream>

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <sys/time.h>
#include <time.h>

using namespace std;

#define UNASSIGNED 0

#define N 16
#define BOX 4

#define lvl 20

bool found = false;

typedef struct Board {
    int grid[N][N];
} Board;

bool FindUnassignedLocation(int grid[N][N], int &row, int &col);
bool isSafe(int grid[N][N], int row, int col, int num);
void printGrid(int grid[N][N]);

void SolveSudoku_omp(Board b, int depth = 0)
{
    int row, col;

    if (!FindUnassignedLocation(b.grid, row, col)) {
        /* Print result */
        printGrid(b.grid);

        /* Found! */
        found = true;
    }

    for (int num = 1; num <= N && !found; num++)
    {
        if (isSafe(b.grid, row, col, num))
        {
            b.grid[row][col] = num;

            #pragma omp task firstprivate(b, depth) \
                             final (depth > lvl)
            SolveSudoku_omp(b, depth + 1);

            b.grid[row][col] = UNASSIGNED;
        }
    }

    /* Condition for task wait */
    if (depth <= lvl) {
      #pragma omp taskwait
    }
}

/* Assign values to all unassigned locations for Sudoku solution */
bool SolveSudoku_seq(Board b)
{
    int row, col;

    if (!FindUnassignedLocation(b.grid, row, col)) {
        printGrid(b.grid);

        return true;
    }

    for (int num = 1; num <= N; num++)
    {
      if (isSafe(b.grid, row, col, num))
      {
          b.grid[row][col] = num;

          if (SolveSudoku_seq(b))
              return true;

          b.grid[row][col] = UNASSIGNED;
      }
    }

    return false;
}

/* Searches the grid to find an entry that is still unassigned. */
bool FindUnassignedLocation(int grid[N][N], int &row, int &col)
{
    for (row = 0; row < N; row++)
        for (col = 0; col < N; col++)
            if (grid[row][col] == UNASSIGNED)
                return true;

    return false;
}

/* Returns whether any assigned entry n the specified row matches
   the given number. */
bool UsedInRow(int grid[N][N], int row, int num)
{
    for (int col = 0; col < N; col++)
    if (grid[row][col] == num)
        return true;
    return false;
}

/* Returns whether any assigned entry in the specified column matches
   the given number. */
bool UsedInCol(int grid[N][N], int col, int num)
{
    for (int row = 0; row < N; row++)
        if (grid[row][col] == num)
            return true;
    return false;
}

/* Returns whether any assigned entry within the specified box matches
   the given number. */
bool UsedInBox(int grid[N][N], int boxStartRow, int boxStartCol, int num)
{
    for (int row = 0; row < BOX; row++)
        for (int col = 0; col < BOX; col++)
            if (grid[row+boxStartRow][col+boxStartCol] == num)
                return true;

    return false;
}

/* Returns whether it will be legal to assign num to the given row,col location.
 */
bool isSafe(int grid[N][N], int row, int col, int num)
{
    return !UsedInRow(grid, row, num) && !UsedInCol(grid, col, num) &&
       !UsedInBox(grid, row - row % BOX, col - col % BOX, num);
}

void printGrid(int grid[N][N])
{
    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
            cout << grid[row][col] << "  ";
        cout << endl;
    }
}

int main()
{
    /* To complete! */
    Board b =       {{{6, 1, 0, 16, 11, 14, 0, 0, 0, 0, 0, 12, 15, 8, 4, 2},
                      {0, 12, 7, 0, 4, 0, 0, 0, 0, 0, 6, 0, 0, 0, 1, 16},
                      {10, 0, 0, 4, 0, 0, 0, 9, 3, 0, 16, 7, 0, 0, 11, 5},
                      {11, 8, 9, 14, 16, 12, 0, 1, 0, 4, 0, 0, 0, 0, 0, 13},
                      {2, 0, 0, 5, 0, 0, 0, 12, 0, 14, 11, 0, 8, 15, 0, 0},
                      {0, 6, 0, 7, 0, 10, 0, 13, 1, 15, 4, 0, 0, 0, 0, 0},
                      {0, 0, 15, 0, 0, 0, 0, 5, 0, 0, 0, 0, 6, 3, 2, 9},
                      {0, 9, 0, 0, 0, 15, 0, 0, 16, 0, 0, 0, 0, 0, 0, 14},
                      {4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0, 0, 0, 12, 0},
                      {5, 11, 1, 2, 0, 0, 0, 0, 8, 0, 0, 0, 0, 16, 0, 0},
                      {0, 0, 0, 0, 0, 5, 2, 4, 14, 0, 10, 0, 1, 0, 9, 0},
                      {0, 0, 12, 9, 0, 1, 13, 0, 4, 0, 0, 0, 11, 0, 0, 15},
                      {15, 0, 0, 0, 0, 0, 4, 0, 6, 0, 3, 11, 12, 9, 7, 10},
                      {16, 3, 0, 0, 5, 7, 0, 10, 12, 0, 0, 0, 2, 0, 0, 1},
                      {7, 2, 0, 0, 0, 9, 0, 0, 0, 0, 0, 16, 0, 5, 8, 0},
                      {9, 10, 6, 8, 1, 0, 0, 0, 0, 0, 14, 13, 16, 0, 15, 11}}};

    /* Set time measure */
    struct timeval start, stop;
    unsigned long elapsed_parallel, elapsed_seq;
    gettimeofday(&start, NULL);

    #pragma omp parallel
    #pragma omp single
    SolveSudoku_omp(b);

    if (!found)
        cout << "No solution exists" << endl;

    /* Get elapsed time */
    gettimeofday(&stop, NULL);
    elapsed_parallel = 1000000 * (stop.tv_sec - start.tv_sec);
    elapsed_parallel += stop.tv_usec - start.tv_usec;

    Board a =       {{{6, 1, 0, 16, 11, 14, 0, 0, 0, 0, 0, 12, 15, 8, 4, 2},
                      {0, 12, 7, 0, 4, 0, 0, 0, 0, 0, 6, 0, 0, 0, 1, 16},
                      {10, 0, 0, 4, 0, 0, 0, 9, 3, 0, 16, 7, 0, 0, 11, 5},
                      {11, 8, 9, 14, 16, 12, 0, 1, 0, 4, 0, 0, 0, 0, 0, 13},
                      {2, 0, 0, 5, 0, 0, 0, 12, 0, 14, 11, 0, 8, 15, 0, 0},
                      {0, 6, 0, 7, 0, 10, 0, 13, 1, 15, 4, 0, 0, 0, 0, 0},
                      {0, 0, 15, 0, 0, 0, 0, 5, 0, 0, 0, 0, 6, 3, 2, 9},
                      {0, 9, 0, 0, 0, 15, 0, 0, 16, 0, 0, 0, 0, 0, 0, 14},
                      {4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 2, 0, 0, 0, 12, 0},
                      {5, 11, 1, 2, 0, 0, 0, 0, 8, 0, 0, 0, 0, 16, 0, 0},
                      {0, 0, 0, 0, 0, 5, 2, 4, 14, 0, 10, 0, 1, 0, 9, 0},
                      {0, 0, 12, 9, 0, 1, 13, 0, 4, 0, 0, 0, 11, 0, 0, 15},
                      {15, 0, 0, 0, 0, 0, 4, 0, 6, 0, 3, 11, 12, 9, 7, 10},
                      {16, 3, 0, 0, 5, 7, 0, 10, 12, 0, 0, 0, 2, 0, 0, 1},
                      {7, 2, 0, 0, 0, 9, 0, 0, 0, 0, 0, 16, 0, 5, 8, 0},
                      {9, 10, 6, 8, 1, 0, 0, 0, 0, 0, 14, 13, 16, 0, 15, 11}}};


    gettimeofday(&start, NULL);

   if (SolveSudoku_seq(a) == false)
        cout << "No solution exists" << endl;

    /* Get elapsed time */
    gettimeofday(&stop, NULL);
    elapsed_seq = 1000000 * (stop.tv_sec - start.tv_sec);
    elapsed_seq += stop.tv_usec - start.tv_usec;

    cout << "Parallel Time: " << elapsed_parallel << endl;
    cout << "Sequential Time: " << elapsed_seq << endl;



    return 0;
}

