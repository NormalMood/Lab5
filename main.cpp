#include <iostream>
#include <vector>
#include <mpi.h>

void printMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}

void gaussElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, int n) {
    for (int i = 0; i < n; ++i) {
        int pivot_row = i;
        double max_val = abs(A[i][i]);
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > max_val) {
                max_val = abs(A[k][i]);
                pivot_row = k;
            }
        }
        if (pivot_row != i) {
            std::swap(A[i], A[pivot_row]);
            std::swap(b[i], b[pivot_row]);
        }
        for (int j = i + 1; j < n; ++j) {
            double ratio = A[j][i] / A[i][i];
            for (int k = i; k < n; ++k) {
                A[j][k] -= ratio * A[i][k];
            }
            b[j] -= ratio * b[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        b[i] /= A[i][i];
        A[i][i] = 1;
        for (int j = i - 1; j >= 0; --j) {
            b[j] -= A[j][i] * b[i];
            A[j][i] = 0;
        }
    }
}

int main(int argc, char* argv[]) {

    int n = 3; // Размерность СЛАУ
    std::vector<std::vector<double>> A(n, std::vector<double>(n)); // Матрица коэффициентов
    std::vector<double> b(n); // Вектор свободных членов

    return 0;
}

