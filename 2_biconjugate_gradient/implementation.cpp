#include <algorithm>
#include <cstdint>
#include <functional>
#include <optional>
#include <thread>
#include <type_traits>
#include <vector>


namespace Math {
using namespace std;
#include <cmath>
} // namespace Math


// namespace {

template<class T>
concept Number = std::is_same_v<T, float> || std::is_same_v<T, double>;

template<Number T>
void FillFunc(T* vector, const size_t size, std::function<void(T&)> allocator) {
    std::for_each(vector, vector + size, allocator);
}

template<Number T>
void Fill(T* vector, const size_t size, const T value) {
    FillFunc<T>(vector, size, [=](T& el) { el = value; });
}

template<Number T>
T* AllocateVector(const size_t size, std::optional<T> initValue = {}) {
    auto ret = new T[size];
    if (initValue) {
        Fill(ret, size, initValue.value());
    }
    return ret;
}

template<Number T>
T** AllocateMatrix(const size_t size) {
    T** ret = nullptr;
    ret = new T*[size];
    for (size_t i = 0; i < size; ++i) {
        ret[i] = AllocateVector<T>(size);
    }
    return ret;
}

template<Number T>
void Destroy(T* vector) {
    if (!vector) {
        return;
    }
    delete[] vector;
    vector = nullptr;
}

template<Number T>
void Destroy(T** matrix, size_t size) {
    if (!matrix) {
        return;
    }
    for (size_t i = 0; i < size; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
    matrix = nullptr;
}

template<Number T>
void Copy(T** dst, T** src, const size_t size) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

template<Number T>
T** Copy(T** src, const size_t size) {
    auto ret = AllocateMatrix<T>(size);
    Copy(ret, src);
    return ret;
}

template<Number T>
void Transpose(T** src, T** dst, const size_t size) {
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            dst[i][j] = src[j][i];
        }
    }
}

template<Number T>
T** Transpose(T** matrix, const size_t size) {
    auto ret = AllocateMatrix<T>(size);
    Transpose<T>(matrix, ret, size);
    return ret;
}

template<Number T>
T Mult(T* left, T* right, const size_t size) {
    T ret = 0;
    for (int i = 0; i < size; i++) {
        ret += left[i] * right[i];
    }
    return ret;
}

template<Number T>
void Mult(T** matrix, T* vector, T* res, const size_t size) {
    for (int i = 0; i < size; i++) {
        res[i] = Mult<T>(matrix[i], vector, size);
    }
}

template<Number T>
T* Mult(T** matrix, T* vector, const size_t size) {
    auto ret = AllocateVector<T>(size, static_cast<T>(0));
    Mult<T>(matrix, vector, ret, size);
    return ret;
}

inline uint32_t GetLogicalThreadsCount() {
    return std::thread::hardware_concurrency();
}

void SmartMPCycle(size_t start, size_t end, const std::function<void(size_t)>& op, size_t TargetItPerThrd = 200) {
    std::vector<std::thread> procs;
    static const size_t threadsCountSystem   = GetLogicalThreadsCount();
    const size_t        numIterations  = end - start;
    static const size_t threadsCount   = threadsCountSystem * Math::fminl(1, static_cast<double>(numIterations) / TargetItPerThrd);
    const size_t        itersPerThread = numIterations / threadsCount + 1;
    for (size_t tid = 0; tid < threadsCount; ++tid) {
        size_t newStart = start + itersPerThread * tid;
        size_t newEnd   = Math::min(newStart + itersPerThread, end);
        procs.push_back(std::thread{[newStart, newEnd, &op]() {
            for (size_t i = newStart; i < newEnd; ++i) {
                op(i);
            }
        }});
    }
    for (auto& proc: procs) {
        proc.join();
    }
}

template<Number T, bool MP = false>
int SolveBCG(T** A, T* b, T* x, const size_t size, const size_t iterations,
             const T accuracy) {

    // Get transposed A
    auto At = Transpose<T>(A, size);

    T norm = Math::sqrt(Mult(b, b, size));

    // first approach
    Fill<T>(x, size, 1);

    // method init
    auto R      = AllocateVector<T>(size);
    auto biR    = AllocateVector<T>(size);
    auto nR     = AllocateVector<T>(size);
    auto nbiR   = AllocateVector<T>(size);
    auto P      = AllocateVector<T>(size);
    auto nP     = AllocateVector<T>(size);
    auto biP    = AllocateVector<T>(size);
    auto nbiP   = AllocateVector<T>(size);
    auto mAtbiP = AllocateVector<T>(size);
    auto mAP    = Mult<T>(A, x, size);
    if constexpr (MP) {
        SmartMPCycle(0, size, [&](size_t i) {
            R[i] = biR[i] = P[i] = biP[i] = b[i] - mAP[i];
        });
    } else {
        for (size_t i = 0; i < size; ++i) {
            R[i] = biR[i] = P[i] = biP[i] = b[i] - mAP[i];
        }
    }

    // method realization
    bool solved = false;
    size_t iter;
    for (iter = 0; iter < iterations; ++iter) {
        Mult<T>(A,  P,   mAP,    size);
        Mult<T>(At, biP, mAtbiP, size);

        T numerator   = Mult<T>(biR, R,   size);
        T denominator = Mult<T>(biP, mAP, size);
        T alfa = numerator / denominator;

        if constexpr (MP) {
            SmartMPCycle(0, size, [&](size_t i) {
                nR[i] = R[i] - alfa * mAP[i];
                nbiR[i] = biR[i] - alfa * mAtbiP[i];
            });
        } else {
            for(size_t i = 0; i < size; ++i) {
                nR[i] = R[i] - alfa * mAP[i];
                nbiR[i] = biR[i] - alfa * mAtbiP[i];
            }
        }

        denominator = numerator;
        numerator = Mult<T>(nbiR, nR, size);
        T beta = numerator / denominator;

        if constexpr (MP) {
            SmartMPCycle(0, size, [&](size_t i) {
                nP[i] = nR[i] + beta * P[i];
                nbiP[i] = nbiR[i] + beta * biP[i];
            });
        } else {
            for(size_t i = 0; i < size; ++i) {
                nP[i] = nR[i] + beta * P[i];
                nbiP[i] = nbiR[i] + beta * biP[i];
            }
        }

        // accuracy control
        T curAccuracy = Math::sqrt(Mult<T>(R, R, size) / norm);
        if (curAccuracy < accuracy) {
            solved = true;
            break;
        }

        for(size_t i = 0; i < size; ++i)
        {
            x[i] += alfa * P[i];
        }

        // swap current and next-iteration arrays
        std::swap(R, nR);
        std::swap(P, nP);
        std::swap(biR, nbiR);
        std::swap(biP, nbiP);
    }

    // clearing
    Destroy(R);
    Destroy(biR);
    Destroy(nR);
    Destroy(nbiR);
    Destroy(P);
    Destroy(nP);
    Destroy(biP);
    Destroy(nbiP);
    Destroy(mAtbiP);
    Destroy(mAP);
    Destroy(At, size);

    return (solved) ? static_cast<int>(iter) : -1;
}

// } // anonimous namespace

#if defined(WIN32)
#define EXPORT_API __declspec(dllexport)
#define IMPORT_API __declspec(dllimport)
#define CC __cdecl
#else
#define EXPORT_API __attribute__((visibility("default"), used))
#define IMPORT_API
#define CC
#endif

extern "C" {

EXPORT_API int CC SolveBCGF32(float** A, float* b, float* x, const size_t size,
                                 const size_t iterations, const float accuracy) {
    return SolveBCG<float, false>(A, b, x, size, iterations, accuracy);
}
EXPORT_API int CC SolveBCGF64(double** A, double* b, double* x, const size_t size,
                                 const size_t iterations, const double accuracy) {
    return SolveBCG<double, false>(A, b, x, size, iterations, accuracy);
}

EXPORT_API int CC MP_SolveBCGF32(float** A, float* b, float* x, const size_t size,
                                 const size_t iterations, const float accuracy) {
    return SolveBCG<float, true>(A, b, x, size, iterations, accuracy);
}
EXPORT_API int CC MP_SolveBCGF64(double** A, double* b, double* x, const size_t size,
                                 const size_t iterations, const double accuracy) {
    return SolveBCG<double, true>(A, b, x, size, iterations, accuracy);
}

}