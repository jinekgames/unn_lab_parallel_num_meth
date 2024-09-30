#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

namespace {

template<class T>
void Fill(T* vector, const size_t size, const T value) {
    std::for_each(vector, vector + size, [&](T& el) { el = value; });
}

template<class T>
T* AllocateVector(const size_t size, bool init = false, T initValue = {}) {
    auto ret = new T[size];
    if (init) {
        Fill(ret, size, initValue);
    }
    return ret;
}

template<class T>
void Destroy(T* vector) {
    if (!vector) {
        return;
    }
    delete[] vector;
    vector = nullptr;
}

template<class T>
T Mult(T* left, T* right, const size_t size) {
    T ret = 0;
    for (int i = 0; i < size; i++) {
        ret += left[i] * right[i];
    }
    return ret;
}

} // anonimous namespace

#ifdef JNK_TEST_BUILD
struct CRSMatrix
{
    int n  = 0;                   // rows count
    int m  = 0;                   // columns count
    int nz = 0;                   // non-zero elements count
    std::vector<double> val;      // matrix values by rows
    std::vector<int>    colIndex; // columns indeces
    std::vector<int>    rowPtr;   // lines start indices
};
#endif

void Multiplicate(CRSMatrix A, double *x, double* b)
{
    for (int i = 0; i < A.n; i++) {
        b[i] = 0.0;
        for (int j = A.rowPtr[i]; j < A.rowPtr[i + 1]; j++) {
            b[i] += A.val[j] * x[A.colIndex[j]];
        }
    }
}

double* Multiplicate(CRSMatrix A, double *x)
{
    double* b = AllocateVector<double>(A.n);
    Multiplicate(A, x, b);
    return b;
}

CRSMatrix sparse_transpose(const CRSMatrix& input) {
    CRSMatrix res;
    res.n = input.n;
    res.m = input.m;
    res.nz = input.nz;
    res.val = std::vector<double>(input.nz, 0.0);
    res.colIndex = std::vector<int>(input.nz, 0);
    res.rowPtr = std::vector<int>(input.m + 2, 0); // one extra

    // count per column
    for (int i = 0; i < input.nz; ++i) {
        ++res.rowPtr[input.colIndex[i] + 2];
    }

    // from count per column generate new rowPtr (but shifted)
    for (int i = 2; i < res.rowPtr.size(); ++i) {
        // create incremental sum
        res.rowPtr[i] += res.rowPtr[i - 1];
    }

    // perform the main part
    for (int i = 0; i < input.n; ++i) {
        for (int j = input.rowPtr[i]; j < input.rowPtr[i + 1]; ++j) {
            // calculate index to transposed matrix at which we should place current element, and at the same time build final rowPtr
            const int new_index = res.rowPtr[input.colIndex[j] + 1]++;
            res.val[new_index] = input.val[j];
            res.colIndex[new_index] = i;
        }
    }
    res.rowPtr.pop_back(); // pop that one extra

    return res;
}

void SLE_Solver_CRS_BICG(CRSMatrix& A, double* b, double eps, int max_iter,
                         double* x, int& count) {

    using T = double;

    auto size = A.n;

    // Get transposed A
    auto At = sparse_transpose(A);

    T norm = std::sqrt(Mult(b, b, size));

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
    auto mAP    = Multiplicate(A, x);

    // init vectors
// #pragma omp parallel for
    for (size_t i = 0; i < size; ++i) {
        R[i] = biR[i] = P[i] = biP[i] = b[i] - mAP[i];
    }

    // method realization
    bool solved = false;
    size_t iter;
    for (iter = 0; iter < max_iter; ++iter) {
        Multiplicate(A, P, mAP);
        Multiplicate(At, biP, mAtbiP);

        T numerator   = Mult<T>(biR, R,   size);
        T denominator = Mult<T>(biP, mAP, size);
        T alfa = numerator / denominator;

// #pragma omp parallel for default(none) shared(size,R,nR,biR,nbiR,mAP,mAtbiP,alfa)
        for(size_t i = 0; i < size; ++i) {
            nR[i] = R[i] - alfa * mAP[i];
            nbiR[i] = biR[i] - alfa * mAtbiP[i];
        }

        denominator = numerator;
        numerator = Mult<T>(nbiR, nR, size);
        T beta = numerator / denominator;

// #pragma omp parallel for default(none) shared(size,nR,nbiR,P,nP,biP,nbiP,beta)
        for(size_t i = 0; i < size; ++i) {
            nP[i] = nR[i] + beta * P[i];
            nbiP[i] = nbiR[i] + beta * biP[i];
        }

        // accuracy control
        T curAccuracy = std::sqrt(Mult<T>(R, R, size) / norm);
        if (curAccuracy < eps) {
            solved = true;
            break;
        }

// #pragma omp parallel for default(none) shared(size,x,alfa,P)
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

    count = static_cast<int>(iter);

    return;
}
