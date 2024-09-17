#include <cstdint>
#include <vector>
#include <thread>
#include <iomanip>
#include <iostream>

#include <omp.h>

namespace {

template<class T>
T** Convert(T* src, size_t n) {
    auto ret = new T*[n];
    for(size_t i = 0; i < n; ++i) {
        ret[i] = &src[i * n];
    }
    return ret;
}

template<class T>
void Copy(T* dst, T** src, size_t n) {
    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            dst[i * n + j] = src[i][j];
        }
    }
}

template<class T>
T* Convert(T** src, size_t n) {
    auto ret = new T[n * n];

    Copy(ret, src, n);

    return ret;
}

template<class T>
void DeallocInterpret(T** interpret)
{
    delete[] interpret;
}
template<class T>
void Dealloc(T** matrix, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}
template<class T>
void Dealloc(T* matrix)
{
    delete[] matrix;
}

template<class T>
T** Allocate(size_t size) {
    T** ret;
    ret = new T*[size];
    for (size_t i = 0; i < size; ++i) {
        ret[i] = new T[size];
    }
    return ret;
}

template<class T>
T* AllocateSingle(size_t sizeX, size_t sizeY = sizeX) {
    return new T[sizeX * sizeY];
}

template<class T>
void Copy(T** dst, T** src, size_t size,
          size_t srcOffsetX = 0, size_t srcOffsetY = 0,
          size_t dstOffsetX = 0, size_t dstOffsetY = 0) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            dst[i + dstOffsetY][j + dstOffsetX] = src[i + srcOffsetY][j + srcOffsetX];
        }
    }
}

template<class T>
T** Copy(T** src, size_t size,
         size_t srcOffsetX = 0, size_t srcOffsetY = 0) {
    auto out = Allocate<T>(size);
    Copy(out, src, size, srcOffsetX, srcOffsetY);
    return out;
}

template<class T>
void FillZeros(T** matrix, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix[i][j] = 0;
        }
    }
}

template<class T>
void Print(T** matrix, size_t n) {
    std::cout << "[" << n << " x " << n << "]:" << std::endl;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cout << std::setw(5) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

inline uint32_t GetLogicalThreadsCount() {
    return std::thread::hardware_concurrency();
}

size_t CalculateBlocksSize(size_t mtxsize)
{
    return 1;
    const auto orig = mtxsize;
    for (;;) {
        const auto threads = orig / mtxsize;
        if (threads >= GetLogicalThreadsCount() || mtxsize <= 3) {
            return threads;
        }

        if (mtxsize % 2 == 0) {
            mtxsize /= 2;
        } else if (mtxsize % 3 == 0) {
            mtxsize /= 3;
        } else if (mtxsize % 5 == 0) {
            mtxsize /= 5;
        } else if (mtxsize % 7 == 0) {
            mtxsize /= 7;
        } else if (mtxsize % 11 == 0) {
            mtxsize /= 11;
        } else if (mtxsize % 13 == 0) {
            mtxsize /= 13;
        } else {
            break;
        }
    }
    return orig / mtxsize;
}

} // anonimous namespace

template<class T>
std::tuple<T**, T**> ExtractDecomposition(T** matrix, size_t n) {
    auto L = Allocate<T>(n);
    auto U = Allocate<T>(n);
    FillZeros(L, n);
    FillZeros(U, n);
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            L[i][j] = matrix[i][j];
        }
    }
    for (size_t i = 0; i < n; ++i) {
        L[i][i] = 1;
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            U[i][j] = matrix[i][j];
        }
    }
    return std::make_tuple(L, U);
}

template<class T>
std::tuple<T*, T*> ExtractDecomposition(T* matrix, T* L, T* U, size_t n) {

    std::fill(L, L + n * n, 0);
    std::fill(U, U + n * n, 0);

#pragma omp parallel for
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            L[i * n + j] = matrix[i * n + j];
        }
    }

#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        L[i * n + i] = 1;
    }

#pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            U[i * n + j] = matrix[i * n + j];
        }
    }

    return std::make_tuple(L, U);
}

template<class T>
void RectangleLUijParallel(T* A, const size_t m, const size_t n,
                           const size_t i0, const size_t j0, size_t size)
{
    size_t min = std::min(m - 1, n);

    for (size_t i = 0; i < min; ++i) {
#pragma omp parallel for default(none) shared(A, i, min, m, n, i0, j0)
        for (size_t j = i + 1; j < m; ++j)
            A[(j + i0) * size + i + j0] /= A[(i + i0) * size + i + j0];
        if (i < n - 1) {
#pragma omp parallel for default(none) shared(A, i, min, m, n, i0, j0)
            for (size_t j = i + 1; j < m; ++j)
                for (size_t k = i + 1; k < n; ++k)
                    A[(j + i0) * size + k + j0] -=
                        A[(j + i0) * size + i + j0] *
                        A[(i + i0) * size + k + j0];
        }
    }
}

template<class T>
void LU_Decomposition_Block(T* A, size_t size, size_t blockSize) {

    T* U23 = AllocateSingle<double>(blockSize, size - blockSize);
    T* L32 = AllocateSingle<double>(size - blockSize, blockSize);

    for (size_t i = 0; i < size - 1; i += blockSize) {

        RectangleLUijParallel(A, size - i, blockSize, i, i, size);

#pragma omp parallel default(none) shared(A,L32,U23,i)
        {
#pragma omp for
            for (size_t j = 0; j < size - blockSize - i; ++j) {
                for (size_t k = 0; k < blockSize; ++k) {
                    T sum = 0;
                    for (size_t q = 0; q < k; ++q) {
                        sum += A[(i + k) * size + (i + q)] * A[(q + i) * size + blockSize + i + j];
                    }
                    A[(k + i) * size + blockSize + i + j] -= sum;
                }

            }

#pragma omp for
            for (size_t j = 0; j < size - blockSize - i; ++j) {
                for (size_t k = 0; k < blockSize; ++k) {
                    U23[k * (size - i - blockSize) + j] = A[(i + k) * size + j + i + blockSize];
                    L32[j * blockSize + k] = A[(i + j + blockSize) * size + k + i];
                }
            }


#pragma omp for
            for (size_t j = blockSize + i; j < size; ++j) {
                for (size_t c = 0; c < blockSize; ++c) {
                    for (size_t k = blockSize + i; k < size; ++k) {
                        A[j * size + k] -= L32[(j - i - blockSize) * blockSize + c] * U23[c * (size - i - blockSize) + k - (i + blockSize)];
                    }

                }
            }
        }


    }

    Dealloc(U23);
    Dealloc(L32);
}

void LU_Decomposition(double* A, double* L, double* U, int size) {

    const auto optimalSize = CalculateBlocksSize(size);
    // std::cout << "optimalSize - " << optimalSize << std::endl;
    const auto blockSize = size / optimalSize;
    LU_Decomposition_Block(A, size, blockSize);
    auto ret = ExtractDecomposition(A, L, U, size);
}

#if defined(WIN32)
#define EXPORT_API __declspec(dllexport)
#define IMPORT_API __declspec(dllimport)
#define CC __cdecl
#else
#define EXPORT_API __attribute__((visibility("default"), used))
#define IMPORT_API
#define CC
#endif

#define USE_BLOCK true //false //true

extern "C" {

void EXPORT_LU_Decomposition(double** A, double** L, double** U, int n)
{
    auto Aconv = Convert<double>(A, n);
    auto Lconv = Convert<double>(L, n);
    auto Uconv = Convert<double>(U, n);
    LU_Decomposition(Aconv, Lconv, Uconv, n);
    auto Lorig = Convert<double>(Lconv, n);
    auto Uorig = Convert<double>(Uconv, n);
    Copy<double>(L, Lorig, n);
    Copy<double>(U, Uorig, n);
// #if USE_BLOCK
//     MP_Decompose_Block(A, L, U, n);
// #else
//     Decompose(A, n);
//     auto ret = ExtractDecomposition(A, n);
//     Copy(L, std::get<0>(ret), n);
//     Copy(U, std::get<1>(ret), n);
// #endif
}

}