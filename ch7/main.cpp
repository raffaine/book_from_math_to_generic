#include <iostream>
#include <chrono>
#include <random>
#include <functional>
#include <concepts>
#include <string>

// Simple timer to clock the implementations

#define PRINT_TIME_AS_MS 0 // In case measures are more familiar in milliseconds

using namespace std::chrono;

class Timer {
    public:
        Timer() : _start(high_resolution_clock::now()) {}

        microseconds clock() {
            auto end = high_resolution_clock::now();
            auto diff = end - _start;
            _start = end;
            return duration_cast<microseconds>(diff);
        }
private:
    high_resolution_clock::time_point _start;    
};

std::ostream& operator<<(std::ostream& os, microseconds t) {
#if PRINT_TIME_AS_MS
    os << duration_cast<milliseconds>(t).count() << "ms";
#else
    os << t.count() << " microseconds";
#endif
    return os;
}


////// Math functions

// Concepts

template <typename T>
concept Additive = requires (T a, T b) {
    a + b;
};

template <typename T>
concept NoncommutativeAdditiveSemigroup = std::regular<T> && Additive<T>;

template <typename T>
concept Integer = std::integral<T>;


template <Integer T>
bool odd(T n) {
    return n & 1;
}

template <Integer T>
T half(T n) {
    return n >> 1;
}

/*
    Given r == 0, calculates n * a through sequential accumulation
    i.e. result == a + a + a ... (repeated n times)
    although uses an approach that does O(log(n)) operations (instead of n)
*/
template <NoncommutativeAdditiveSemigroup A, Integer N>
A multiply_accumulate(A r, N  n, A a) {
    while(true) {
        if (odd(n)) {
            r = r + a;
            if (n == 1) return r;
        }
        n = half(n);
        a = a + a;
    }
}






int main(int argc, char** argv) {
    // Number of inputs to be considered
    const int N = 1000000;

    
    // Seed with a real random value, if available
    std::random_device r;
 
    // Choose a random mean between 1000 and 10 000
    std::default_random_engine e(r());
    std::uniform_int_distribution<int> udist(1000, 10000);

    // Seed the input vector with a fresh set of random numbers
    std::vector<int> input(N);
    for (int i=0; i < N; i++) {
        input[i] = udist(e);
    }

    // Create 2 output vectors, one for each approach
    std::vector<int> output(N), output2(N);

    // Time N multiplications using our custom approach
    Timer t;
    for (int i=0; i < N; i++) {
        output2[i] = multiply_accumulate(0, input[i], input[N - i - 1]);
    }
    auto custom_duration = t.clock();
    
    // Time N multiplications using the hardware implementation
    for (int i=0; i < N; i++) {
        output[i] = input[i] * input[N - i - 1];
    }
    auto base_duration = t.clock();
    
    std::cout << "Custom Mult took " << custom_duration << std::endl;    
    std::cout << "Multiplication took " << base_duration << std::endl;
    std::cout << "Custom runs " << ((double) base_duration.count()) / custom_duration.count() << " times faster." << std::endl;
    if ( auto failed_ind = std::invoke([&] {
        for (int i=0; i < N; i++) {
            if (output[i] != output2[i]) return i;
        }
        return -1;
    }) < 0 ) {
        std::cout << "Results match" << std::endl;
    } else {
        std::cout << "Results don't match " << output[failed_ind] << " != " << output2[failed_ind] << std::endl;
        std::cout << "Values multiplied: " << input[failed_ind] << " x " << input[N - failed_ind - 1] << std::endl;
    }

    return 0;
}
