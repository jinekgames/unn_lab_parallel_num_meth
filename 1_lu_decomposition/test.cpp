#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <optional>
#include <string>
#include <sstream>
#include <tuple>
#include <vector>


extern "C" {
void DecomposeF64(double** matrix, size_t size);
void MP_DecomposeF64(double** matrix, size_t size);
}

#ifndef MATRIX_OUTPUT_CELL_WIDTH
#define MATRIX_OUTPUT_CELL_WIDTH 4
#endif

#define USE_MP        0
#define PRINT_CONTENT 0
#define MATRIX_SIZE   5000


template<class T, size_t Size>
class Matrix {
public:

    Matrix(std::optional<std::string> name = {})
        : m_Name{name} {
        Allocate();
    }
    Matrix(const Matrix& right)
        : m_Name{right.m_Name} {
        Allocate();
        Copy(right);
    }
    Matrix(Matrix&& right)
        : m_Name{right.m_Name} {
        Allocate();
        Copy(right);
        right.Destroy();
    }
    ~Matrix() {
        Destroy();
    }
    Matrix& operator = (const Matrix& right) {
        m_Name = right.m_Name;
        Allocate();
        Copy(right);
    }
    Matrix& operator = (Matrix&& right) {
        m_Name = std::move(right.m_Name);
        Allocate();
        Copy(right);
        right.Destroy();
    }

    void Copy(const Matrix& right) {
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                m_Data[i][j] = right[i][j];
            }
        }
    }

    void FillNumbers() {
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                m_Data[i][j] = 1 + static_cast<double>(i) * Size + j;
            }
        }
    }

    void FillZeros() {
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                m_Data[i][j] = 0;
            }
        }
    }

    T** Get() { return m_Data; }

    size_t GetSize() const {
        return Size;
    }

    std::string GetName() const {
        return m_Name.value_or("");
    }

    void SetName(const std::string& name) {
        m_Name.emplace(name);
    }

    void UnsetName() {
        m_Name.reset();
    }

    std::string ToStr() const {
        std::stringstream out;
        if (m_Name) {
            out << m_Name.value() << " ";
        }
        out << "[" << Size << " x " << Size << "]:" << std::endl;
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                out << std::setw(MATRIX_OUTPUT_CELL_WIDTH) << m_Data[i][j] << " ";
            }
            out << std::endl;
        }
        return out.str();
    }

    friend std::ostream& operator << (std::ostream& stream, const Matrix& obj) {
        return stream << obj.ToStr();
    }

    operator T**() {
        return Get();
    }

    template<class I>
    void operator = (const I& input) {
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                m_Data[i][j] = static_cast<T>(input[i][j]);
            }
        }
    }

    T*& operator [] (size_t i) {
        return m_Data[i];
    }

    const T* operator [] (size_t i) const {
        return m_Data[i];
    }

    Matrix operator * (const Matrix& right) {
        Matrix<T, Size> ret;
        for (int i = 0; i < Size; i++) {
            for (int j = 0; j < Size; j++) {
                ret[i][j] = 0;
                for (int k = 0; k < Size; k++) {
                    ret[i][j] += m_Data[i][k] * right[k][j];
                }
            }
        }
        if (m_Name && right.m_Name) {
            ret.SetName(GetName() + " * " + right.GetName());
        }
        return ret;
    }

private:

    void Allocate() {
        if (m_Data) {
            Destroy();
        }
        m_Data = new T*[Size];
        for (size_t i = 0; i < Size; ++i) {
            m_Data[i] = new T[Size];
        }
    }

    void Destroy() {
        if (!m_Data) {
            return;
        }
        for (size_t i = 0; i < Size; ++i) {
            delete[] m_Data[i];
        }
        delete[] m_Data;
        m_Data = nullptr;
    }

private:

    T** m_Data = nullptr;
    std::optional<std::string> m_Name;
};


template<class T, size_t Size>
std::tuple<Matrix<T, Size>, Matrix<T, Size>> ExtractDecoposition(const Matrix<T, Size>& matrix) {
    Matrix<T, Size> L{"L"};
    Matrix<T, Size> U{"U"};
    L.FillZeros();
    U.FillZeros();
    for (size_t i = 1; i < Size; ++i) {
        for (size_t j = 0; j < i; ++j) {
            L[i][j] = matrix[i][j];
        }
    }
    for (size_t i = 0; i < Size; ++i) {
        L[i][i] = 1;
    }
    for (size_t i = 0; i < Size; ++i) {
        for (size_t j = i; j < Size; ++j) {
            U[i][j] = matrix[i][j];
        }
    }
    return std::make_tuple(L, U);
}

inline uint64_t GetTick() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

class Timer final {

public:

    Timer(const std::string& name = "Unnamed")
        : m_Name(name)
        , m_StartTime(GetTick()) {}

    Timer(const Timer&) = delete;
    Timer(Timer&&)      = delete;

    ~Timer() {
        if (!m_IsFinished)
            Finish();
    }

    Timer& operator = (const Timer&) = delete;
    Timer& operator = (Timer&&)      = delete;

    bool IsFinished() { return m_IsFinished; }

    void Finish() {
        auto res = GetResult();
        std::cout << m_Name << " timer: " << res << " us" << std::endl;;
        m_IsFinished = true;
    }

    uint64_t GetResult() {
        return GetTick() - m_StartTime;
    }

private:

    bool m_IsFinished = false;
    uint64_t m_StartTime;
    std::string m_Name;

};


int main() {

    constexpr size_t size = MATRIX_SIZE;
    Matrix<double, size> matrix{"Sample"};

    // matrix = std::vector<std::vector<double>>{
    //     { 1, 1, 1, 0, 0 },
    //     { 0, 2, 0, 0, 0 },
    //     { 1, 0, 2, 0, 1 },
    //     { 0, 0, 0, 4, 0 },
    //     { 0, 1, 0, 1, 4 }
    // };
    matrix.FillNumbers();
#if PRINT_CONTENT
    std::cout << matrix << std::endl;
#endif

    Timer timer{"Execution"};
#if USE_MP
    MP_DecomposeF64(matrix, matrix.GetSize());
#else
    DecomposeF64(matrix, matrix.GetSize());
#endif
    timer.Finish();

#if PRINT_CONTENT
    matrix.SetName("Decomposition result");
    std::cout << std::endl << matrix << std::endl;

    auto [L, U] = ExtractDecoposition(matrix);
    std::cout << L << U << std::endl;
    std::cout << L * U << std::endl;
#endif

    return 0;
}