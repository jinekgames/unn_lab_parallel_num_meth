#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#define OUTPUT_CELL_WIDTH -1
#include "./helpers.h"
#include "./math_types.h"


#define USE_MP        false
#define PRINT_CONTENT false
#define MATRIX_SIZE   5000
#define ITERATIONS    200000
#define ACCURASY      0.1


extern "C" {
int SolveBCGF64(double** A, double* b, double* x, const size_t size,
                const size_t iterations, const double accuracy);
int MP_SolveBCGF64(double** A, double* b, double* x, const size_t size,
                   const size_t iterations, const double accuracy);
}


int main() {

    constexpr size_t size = MATRIX_SIZE;
    Matrix<double, size> matrixA{"A"};
    Vector<double, size> vectorB{"b"};
    Vector<double, size> vectorX{"x"};

#if MATRIX_SIZE == 3
    matrixA = std::vector<std::vector<double>>{
        { 1, 1, 1, },
        { 1, 1, 0, },
        { 0, 1, 1, },
    };
    vectorB = std::vector<double>{ 1, 2, 3, };
#else
    matrixA.FillDiagonal(5);
    // matrixA.FillNumbers();
    vectorB.FillNumbers();
#endif
#if PRINT_CONTENT
    std::cout << matrixA << vectorB << std::endl;
#endif

    Timer timer{"Execution"};
    auto iterations =
#if USE_MP
    MP_SolveBCGF64(matrixA, vectorB, vectorX, matrixA.GetSize(), ITERATIONS, ACCURASY);
#else
    SolveBCGF64(matrixA, vectorB, vectorX, matrixA.GetSize(), ITERATIONS, ACCURASY);
#endif
    timer.Finish();
    if (iterations < 0) {
        std::cout << std::endl << "Not solved in " << ITERATIONS <<
                    " iteration(s) (accur.: "<< ACCURASY << ")" << std::endl;
    } else {
        std::cout << std::endl << "Finished in " << iterations << " of " << ITERATIONS <<
                    " iteration(s) (accur.: "<< ACCURASY << ")" << std::endl;
    }

#if PRINT_CONTENT
    std::cout << std::endl << vectorX << std::endl;
    std::cout << matrixA * vectorX << std::endl;
#endif

    return 0;
}