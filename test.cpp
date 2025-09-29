#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <fstream>

// Function to generate a random matrix with values between -1000 and 1000
std::vector<std::vector<double>> generate_random_matrix(int n) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matrix[i][j] = static_cast<double>(std::rand()) / RAND_MAX * 2000.0 - 1000.0;
        }
    }
    return matrix;
}

// Function to calculate the determinant of a matrix (using recursion for simplicity)
double determinant(const std::vector<std::vector<double>>& matrix, int n) {
    if (n == 1) {
        return matrix[0][0];
    }
    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    double det = 0.0;
    std::vector<std::vector<double>> submatrix(n - 1, std::vector<double>(n - 1));

    for (int i = 0; i < n; ++i) {
        int submatrix_row = 0;
        for (int j = 1; j < n; ++j) {
            int submatrix_col = 0;
            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    submatrix[submatrix_row][submatrix_col] = matrix[j][k];
                    submatrix_col++;
                }
            }
            submatrix_row++;
        }

        // Alternate the sign of the determinant based on column index
        det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * determinant(submatrix, n - 1);
    }

    return det;
}

bool is_invertible(const std::vector<std::vector<double>>& matrix, int n) {
    return determinant(matrix, n) != 0;
}


void write_matrix_to_file(const std::vector<std::vector<double>>& matrix, int n, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        
    }

    // Write matrix to file
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << std::setw(10) << std::fixed << std::setprecision(4) << matrix[i][j] << " "; 
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Matrix written to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    int n;
    n=std::stoi(argv[1]);
    std::srand(std::time(0));
    std::vector<std::vector<double>> matrix;

    // Generate an invertible matrix
    do {
        matrix = generate_random_matrix(n);
    } while (!is_invertible(matrix, n));  // Retry if the matrix is not invertible

    // Write the matrix to a file "input.txt"
    write_matrix_to_file(matrix, n, "input.txt");

    return 0;
}
