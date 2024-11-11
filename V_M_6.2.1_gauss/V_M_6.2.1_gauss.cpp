#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double EPSILON = 1e-9;

void print_matrix(const std::vector<std::vector<double>>& A, const std::vector<double>& b) {
    for (size_t i = 0; i < A.size(); i++) {
        for (size_t j = 0; j < A[i].size(); j++) {
            std::cout << std::setw(10) << A[i][j] << " ";
        }
        std::cout << "| " << std::setw(10) << b[i] << std::endl;
    }
}

std::vector<double> gaussian_elimination(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {

        std::cout << "After pivot step" << i + 1 << ":" << std::endl;
        print_matrix(A, b);
        std::cout << std::endl;

        for (int j = i + 1; j < n; j++) {
            double coef = A[j][i] / A[i][i];
            b[j] -= coef * b[i];
            for (int k = i; k < n; k++) {
                A[j][k] -= coef * A[i][k];
            }
        }

        std::cout << "After elimination step " << i + 1 << ":" << std::endl;
        print_matrix(A, b);
        std::cout << std::endl;
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

int main() {
    setlocale(LC_ALL, "");
    int n = 3;
    std::vector<std::vector<double>> A = {
        {3.1, 2.8, 1.9},
        {1.9, 3.1, 2.1},
        {7.5, 3.8, 4.8}
    };

    std::vector<double> b = { 0.2, 2.1, 5.6 };

    std::cout << "Initial Matrix:" << std::endl;
    print_matrix(A, b);
    std::cout << std::endl;

    std::vector<double> x = gaussian_elimination(A, b);

    std::cout << "Solution:" << std::scientific << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "x" << (i + 1) << " = " << x[i] << std::scientific << std::endl;
    }
    std::cout << std::endl;
    std::vector<double> res(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += A[i][j] * x[j];
        }
    }

    std::cout << "Residual:" << std::scientific << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "r" << (i + 1) << " = " << std::scientific << std::abs(res[i] - b[i]) << std::endl;
    }
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    std::cout << "Press Enter to exit...";
    std::cin.get();

    return 0;
}