#include <iostream>
#include <chrono>
#include <string>

// Simple timer to clock the implementations

#define PRINT_TIME_AS_MS 1 // In case measures are more familiar in milliseconds

class Timer {
    public:
        Timer() : _start(std::chrono::high_resolution_clock::now()) {}

        std::chrono::microseconds clock() {
            auto end = std::chrono::high_resolution_clock::now();
            auto diff = end - _start;
            _start = end;
            return std::chrono::duration_cast<std::chrono::microseconds>(diff);
        }
private:
    std::chrono::high_resolution_clock::time_point _start;    
};

std::ostream& operator<<(std::ostream& os, std::chrono::microseconds t) {
#if PRINT_TIME_AS_MS
    os << std::chrono::duration_cast<std::chrono::milliseconds>(t).count() << "ms";
#else
    os << t.count() << "\u00B5s";
#endif
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::plus<T>) {
    os << "x";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::multiplies<T>) {
    os << "^";
    return os;
}
