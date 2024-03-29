#include <iostream>
#include <vector>
#include <mpi.h>

using namespace std;

void printMatrix(const vector<vector<double>>& matrix) {
    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

void gaussElimination(vector<vector<double>>& A, vector<double>& b, int n) {
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
            swap(A[i], A[pivot_row]);
            swap(b[i], b[pivot_row]);
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
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 3; // ����������� ����
    vector<vector<double>> A(n, vector<double>(n)); // ������� �������������
    vector<double> b(n); // ������ ��������� ������

    if (rank == 0) {
        // ������������� ������� A � ������� b
        // ������ �������������:
        A = {{2, -1, 1},
             {3, 3, 9},
             {3, 3, 5}};
        b = {2, 15, 19};
        // ������������� ������ ����� ����������
        for (int i = 1; i < size; ++i) {
            MPI_Send(A[i].data(), n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            MPI_Send(b.data(), n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(A[rank].data(), n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b.data(), n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    cout<<"Matrix A:"<<endl;
    printMatrix(A);

    cout<<"array b:"<<endl;
    for (int i = 0; i < b.size(); i++)
        cout<<b[i]<<endl;

    gaussElimination(A, b, n);

    if (rank == 0) {
        cout << "Solution: ";
        for (double val : b) {
            cout << val << " ";
        }
        cout << endl;
    }

    MPI_Finalize();
    return 0;
}

