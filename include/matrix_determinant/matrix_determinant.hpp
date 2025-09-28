#ifndef MATRIX_DETERMINANT_HPP
#define MATRIX_DETERMINANT_HPP

#include <vector>

namespace matrix_determinant {

bool isSquare(const std::vector<std::vector<double>>& matrix);

double determinant(const std::vector<std::vector<double>>& matrix);

std::vector<double> gaussSolve(const std::vector<std::vector<double>>& matrix, 
                              const std::vector<double>& vector);

}

#endif
