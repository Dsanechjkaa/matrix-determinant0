cmake_minimum_required(VERSION 3.10)
project(matrix_determinant)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Включение заголовочных файлов
include_directories(include)

# Создание библиотеки
add_library(matrix_determinant_lib src/matrix_determinant.cpp)

# Создание исполняемого файла для тестов
add_executable(test_matrix_determinant tests/test_matrix_determinant.cpp)
target_link_libraries(test_matrix_determinant matrix_determinant_lib)

# Включение тестов
enable_testing()
add_test(NAME MatrixDeterminantTests COMMAND test_matrix_determinant)
