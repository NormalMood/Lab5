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

int main(int argc, char* argv[]) {

    int n = 3; // Размерность СЛАУ
    std::vector<std::vector<double>> A(n, std::vector<double>(n)); // Матрица коэффициентов
    std::vector<double> b(n); // Вектор свободных членов

    return 0;
}

