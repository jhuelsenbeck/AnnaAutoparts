#ifndef Hungarian_H
#define Hungarian_H



class Hungarian {

	public:
                Hungarian(int** m, int nr, int nc);
        int     getNumRows(void) { return numRows; }
        int     getNumCols(void) { return numCols; }
        int     getAssignmentCost(void) { return assignmentCost; }

    private:
        int     hungarianInit(int **costMatrix, int rows, int cols);
        void    hungarianSolve(void);
        int     imax(int a, int b) { return (a < b) ? b : a; }
        int     numRows;
        int     numCols;
        int     assignmentCost;
        int**   cost;
        int**   assignment;
        int**   minimizedMatrix;
        int     matrixSize;
};



#endif
