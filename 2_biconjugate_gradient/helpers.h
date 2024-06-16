#pragma once


#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>


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
        std::cout << m_Name << " timer: " << res << " us" << std::endl;
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
