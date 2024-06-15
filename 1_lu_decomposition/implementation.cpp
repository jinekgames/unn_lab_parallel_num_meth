#include <cstdint>
#include <vector>
#include <thread>


namespace {

template<class T>
T** Allocate(size_t size) {
    T** out;
    out = new T*[size];
    for (size_t i = 0; i < size; ++i) {
        out[i] = new T[size];
    }
}

template<class T>
void Copy(T** dst, T** src, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            dst[i][j] = src[i][j];
        }
    }
}

template<class T>
T** Copy(T** src, size_t size) {
    auto out = Allocate(size);
    Copy(out, src);
    return out;
}

template<class T>
void Decompose(T** matrix, size_t size) {
    for (size_t i = 1; i < size; ++i) {
        for (size_t k = 0; k <= i - 1; ++k) {
            matrix[i][k] = matrix[i][k] / matrix[k][k];
            for (size_t j = k + 1; j < size; ++j) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
        }
    }
}

template<class T>
void MP_Decompose(T** matrix, size_t size) {
    std::vector<std::thread> procs;
    for (size_t i = 1; i < size; ++i) {
        procs.push_back(std::thread{[size, i, &matrix]() {
            for (size_t k = 0; k <= i - 1; ++k) {
                matrix[i][k] = matrix[i][k] / matrix[k][k];
                for (size_t j = k + 1; j < size; ++j) {
                    matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
                }
            }
        }});
    }
    for (auto& proc: procs) {
        proc.join();
    }
}

} // anonimous namespace

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

EXPORT_API void CC DecomposeF32(float** matrix, size_t size) {
    return Decompose<float>(matrix, size);
}
EXPORT_API void CC DecomposeF64(double** matrix, size_t size) {
    return Decompose<double>(matrix, size);
}

EXPORT_API void CC MP_DecomposeF32(float** matrix, size_t size) {
    return MP_Decompose<float>(matrix, size);
}
EXPORT_API void CC MP_DecomposeF64(double** matrix, size_t size) {
    return MP_Decompose<double>(matrix, size);
}

}