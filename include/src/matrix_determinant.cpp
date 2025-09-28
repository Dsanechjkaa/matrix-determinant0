#include "matrix_determinant/matrix_determinant.hpp"
#include <iostream>
#include <cmath>

namespace matrix_determinant {

bool isSquare(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    if (n == 0) return true;
    
    for (int i = 0; i < n; i++) {
        if (matrix[i].size() != n) {
            return false;
        }
    }
    return true;
}

double determinant(const std::vector<std::vector<double>>& matrix) {
    if (!isSquare(matrix)) {
        return 0;
    }
    
    int n = matrix.size();
    if (n == 0) return 1;
    if (n == 1) return matrix[0][0];
    if (n == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    
    std::vector<std::vector<double>> a = matrix;
    double det = 1.0;
    
    for (int k = 0; k < n; k++) {
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > abs(a[max_row][k])) {
                max_row = i;
            }
        }
        
        if (max_row != k) {
            std::swap(a[k], a[max_row]);
            det = -det;
        }
        
        if (abs(a[k][k]) < 1e-10) {
            return 0;
        }
        
        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k] / a[k][k];
            for (int j = k; j < n; j++) {
                a[i][j] -= factor * a[k][j];
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        det *= a[i][i];
    }
    
    return det;
}

std::vector<double> gaussSolve(const std::vector<std::vector<double>>& matrix, 
                              const std::vector<double>& vector) {
    int n = matrix.size();
    std::vector<std::vector<double>> a = matrix;
    std::vector<double> y = vector;
    std::vector<double> x(n);
    
    for (int k = 0; k < n; k++) {
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > abs(a[max_row][k])) {
                max_row = i;
            }
        }
        
        if (max_row != k) {
            std::swap(a[k], a[max_row]);
            std::swap(y[k], y[max_row]);
        }
        
        if (abs(a[k][k]) < 1e-10) {
            std::cout << "Система не имеет решения!" << std::endl;
            return std::vector<double>(n, 0);
        }
        
        double temp = a[k][k];
        for (int j = k; j < n; j++) {
            a[k][j] /= temp;
        }
        y[k] /= temp;
        
        for (int i = k + 1; i < n; i++) {
            double factor = a[i][k];
            for (int j = k; j < n; j++) {
                a[i][j] -= factor * a[k][j];
            }
            y[i] -= factor * y[k];
        }
    }
    
    for (int k = n - 1; k >= 0; k--) {
        x[k] = y[k];
        for (int i = 0; i < k; i++) {
            y[i] -= a[i][k] * x[k];
        }
    }
    
    return x;
}

}
