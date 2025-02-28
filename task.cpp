#include "task.h"
//Result=A
// This function computes the inverse of an n×n matrix A using
// Gauss-Jordan elimination with complete pivoting (global maximal element search).
// Because we perform column swaps on A,
// we track these swaps via columnOrder and then reassemble the inverse at the end.
int gaussJordanInverse(std::vector<std::vector<double>>& A,
                                   std::vector<std::vector<double>>& inv,std::vector<int> &columnOrder, std::vector<int> &undo) {
    int n = A.size();
    if (n == 0 || A[0].size() != n)
        throw std::invalid_argument("Matrix A must be square.");

    // Initialize inv as the identity matrix.
    inv.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
        inv[i][i] = 1.0;

    // Initialize columnOrder
    for (int i = 0; i < n; ++i)
        columnOrder[i] = i;

    const double tol = 1e-19; // Tolerance for singularity

    
    for (int step = 0; step < n; ++step) {
        // --- Pivot Search: find maximal element in submatrix A[step..n-1][step..n-1] ---
        int pivotRow = step, pivotCol = step;
        double maxVal = std::fabs(A[step][step]);
        for (int i = step; i < n; ++i) {
            for (int j = step; j < n; ++j) {
                double current = std::fabs(A[i][j]);
                if (current > maxVal) {
                    maxVal = current;
                    pivotRow = i;
                    pivotCol = j;
                }
            }
        }
        
        if (maxVal < tol){
            throw std::runtime_error("Matrix is close to singular.");
            return 0;
        }

        // --- Swap Rows ---
        if (pivotRow != step) {
            std::swap(A[step], A[pivotRow]);
            std::swap(inv[step], inv[pivotRow]);
        }

        // --- Swap Columns ---
        // (We do not swap columns in inv – instead, we record the swap in columnOrder.)
        if (pivotCol != step) {
            for (int i = 0; i < n; ++i) {
                std::swap(A[i][step], A[i][pivotCol]);
            }
            std::swap(columnOrder[step], columnOrder[pivotCol]);
        }

        // --- Normalize Pivot Row ---
        double pivotVal = A[step][step];
        for (int j = step; j < n; ++j) {
            A[step][j] /= pivotVal;
        }
        for (int j = 0; j < n; ++j) {
            inv[step][columnOrder[j]] /= pivotVal;
        }

        // --- Eliminate all other rows ---
        for (int i = 0; i < n; ++i) {
            if (i == step)
                continue;
            double factor = A[i][step];
            for (int j = step; j < n; ++j) {
                A[i][j] -= factor * A[step][j];
            }
            //Below we should perform on column saved in columnOrder not the actual number
            for (int j = 0; j < n; ++j) {
                inv[i][columnOrder[j]] -= factor * inv[step][columnOrder[j]];
            }
        }

    }

    //find inverse permutations
    for(int q=0;q<n;q++){
        undo[columnOrder[q]]=q;
    }
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++)
            A[i][j]=inv[undo[i]][j];
    }

    return 1;
}


