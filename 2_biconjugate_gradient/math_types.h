#pragma once


#include <optional>
#include <string>
#include <sstream>


namespace Math {
using namespace std;
#include <cmath>
} // namespace Math


#ifndef OUTPUT_CELL_WIDTH
#define OUTPUT_CELL_WIDTH 4
#endif
#define OUTPUT_CELL_WIDTH_ADJUST 1


template<class T, size_t Size>
class Matrix;


template<class T, size_t Size>
class Vector {
public:

    Vector(std::optional<std::string> name = {})
        : m_Name{name} {
        Allocate();
    }
    Vector(const Vector& right)
        : m_Name{right.m_Name} {
        Allocate();
        Copy(right);
    }
    Vector(Vector&& right)
        : m_Name{right.m_Name} {
        Allocate();
        Copy(right);
        right.Destroy();
    }
    ~Vector() {
        Destroy();
    }
    Vector& operator = (const Vector& right) {
        m_Name = right.m_Name;
        Allocate();
        Copy(right);
    }
    Vector& operator = (Vector&& right) {
        m_Name = std::move(right.m_Name);
        Allocate();
        Copy(right);
        right.Destroy();
    }

    void Copy(const Vector& right) {
        for (size_t i = 0; i < Size; ++i) {
            m_Data[i] = right[i];
        }
    }

    void FillNumbers() {
        for (size_t i = 0; i < Size; ++i) {
            m_Data[i] = 1 + static_cast<double>(i);
        }
    }

    void FillZeros() {
        for (size_t i = 0; i < Size; ++i) {
            m_Data[i] = 0;
        }
    }

    T* Get() { return m_Data; }

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
        out << "[" << Size << "]:" << std::endl;
        int width = OUTPUT_CELL_WIDTH;
        if (width < 0) {
            width = GetMaxWidth() + OUTPUT_CELL_WIDTH_ADJUST;
        }
        for (size_t i = 0; i < Size; ++i) {
            out << std::setw(width) << static_cast<T>(m_Data[i]) << " ";
        }
        out << std::endl;
        return out.str();
    }

    friend std::ostream& operator << (std::ostream& stream, const Vector& obj) {
        return stream << obj.ToStr();
    }

    operator T*() {
        return m_Data;
    }

    template<class I>
    void operator = (const I& input) {
        for (size_t i = 0; i < Size; ++i) {
            m_Data[i] = static_cast<T>(input[i]);
        }
    }

    T& operator [] (size_t i) {
        return m_Data[i];
    }

    const T& operator [] (size_t i) const {
        return m_Data[i];
    }

    T operator * (const Vector& right) {
        T ret = 0;
        for (int i = 0; i < Size; i++) {
            ret += m_Data[i] * right[i];
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
        m_Data = new T[Size];
    }

    void Destroy() {
        if (!m_Data) {
            return;
        }
        delete[] m_Data;
        m_Data = nullptr;
    }

    size_t GetMaxWidth() const {
        size_t ret = 0;
        for (size_t i = 0; i < Size; ++i) {
            std::stringstream out;
            out << static_cast<T>(m_Data[i]);
            size_t len = out.str().length();
            if (len > ret) {
                ret = len;
            }
        }
        return ret;
    }

private:

    T* m_Data = nullptr;
    std::optional<std::string> m_Name;

    friend class Matrix<T, Size>;
};


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

    void FillDiagonal(size_t width = 1) {
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                if (Math::llabs(static_cast<int64_t>(j) - static_cast<int64_t>(i)) <= width / 2) {
                    m_Data[i][j] = 1;
                } else {
                    m_Data[i][j] = 0;
                }
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
        int width = OUTPUT_CELL_WIDTH;
        if (width < 0) {
            width = GetMaxWidth() + OUTPUT_CELL_WIDTH_ADJUST;
        }
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                out << std::setw(width) << static_cast<T>(m_Data[i][j]) << " ";
            }
            out << std::endl;
        }
        return out.str();
    }

    friend std::ostream& operator << (std::ostream& stream, const Matrix& obj) {
        return stream << obj.ToStr();
    }

    operator T**() {
        return m_Data;
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

    const T*& operator [] (size_t i) const {
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

    Vector<T, Size> operator * (const Vector<T, Size>& right) {
        Vector<T, Size> ret;
        for (int i = 0; i < Size; i++) {
            ret[i] = 0;
            for (size_t j = 0; j < Size; ++j){
                ret[i] += m_Data[i][j] * right[j];
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

    size_t GetMaxWidth() const {
        size_t ret = 0;
        for (size_t i = 0; i < Size; ++i) {
            for (size_t j = 0; j < Size; ++j) {
                std::stringstream out;
                out << static_cast<T>(m_Data[i][j]);
                size_t len = out.str().length();
                if (len > ret) {
                    ret = len;
                }
            }
        }
        return ret;
    }

private:

    T** m_Data = nullptr;
    std::optional<std::string> m_Name;
};
