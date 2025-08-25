#include <iostream>
#include <chrono>
#include <random>
#include <functional>
#include <concepts>
#include <string>

// Simple timer to clock the implementations

#define PRINT_TIME_AS_MS 1 // In case measures are more familiar in milliseconds

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
concept Multiplicative = requires (T a, T b) {
    a * b;
};

template <typename T>
concept NoncommutativeAdditiveSemigroup = std::regular<T> && Additive<T>;

template <typename T>
concept NoncommutativeMultiplicativeSemigroup = std::regular<T> && Multiplicative<T>;

template <typename T>
concept Integer = std::integral<T>;

template <typename T>
concept NoncommutativeAdditiveMonoid = NoncommutativeAdditiveSemigroup<T> && requires(T a) {
    a + T(0) == a;
    T(0) + a == a;
};

template <typename T>
concept NoncommutativeMultiplicativeMonoid = NoncommutativeMultiplicativeSemigroup<T> && requires(T a) {
    a * T(1) == a;
    T(1) * a == a;
};

template <typename T>
T reciprocal(T a) {
    return -a;
}

template <typename T>
concept NoncommutativeAdditiveGroup = NoncommutativeAdditiveMonoid<T> && requires(T a) {
    a + reciprocal(a) == T(0);
    reciprocal(a) + a == T(0);
};

template <typename T>
T inverse(T a) {
    return T(1) / a;
}

template <typename T>
concept NoncommutativeMultiplicativeGroup = NoncommutativeMultiplicativeMonoid<T> && requires(T a) {
    a * inverse(a) == T(1);
    inverse(a) * a == T(1);
};

// Operations

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
A multiply_accumulate_semigroup(A r, N  n, A a) {
    // A precondition for correctness is that n >= 0
    if (n == 0) return r;
    while(true) {
        if (odd(n)) {
            r = r + a;
            if (n == 1) return r;
        }
        n = half(n);
        a = a + a;
    }
}

template <NoncommutativeAdditiveSemigroup A, Integer N>
A multiply_semigroup(N n, A a) {
    // A precondition is that n > 0
    while (!odd(n)) {
        a = a + a;
        n = half(n);
    }
    if (n == 1) return a;
    return multiply_accumulate_semigroup(a, half(n - 1), a + a);
}

template <NoncommutativeAdditiveMonoid A, Integer N>
A multiply_monoid(N n, A a) {
    if (n == 0) return A(0);
    return multiply_semigroup(n, a);
}

template <NoncommutativeAdditiveGroup A, Integer N>
A multiply_group(N n, A a) {
    if (n < 0) {
        n = -n;
        a = -a;
    }
    return multiply_monoid(n, a);
}

template <NoncommutativeMultiplicativeSemigroup A, Integer N>
A power_accumulate_semigroup(A r, A a, N n) {
    // A precondition for correctness is that n >= 0
    if (n == 0) return r;
    while(true) {
        if (odd(n)) {
            r = r * a;
            if (n == 1) return r;
        }
        n = half(n);
        a = a * a;
    }
}

template <NoncommutativeMultiplicativeSemigroup A, Integer N>
A power_semigroup(A a, N n) {
    while (!odd(n)) {
        a = a * a;
        n = half(n);
    }
    if (n == 1) return a;
    return power_accumulate_semigroup(a, a * a, half(n - 1));
}

template <NoncommutativeMultiplicativeMonoid A, Integer N>
A power_monoid(A a, N n) {
    if (n == 0) return A(0);
    return power_semigroup(a, n);
}

template <NoncommutativeMultiplicativeGroup A, Integer N>
A power_group(A a, N n) {
    if (n < 0) {
        n = -n;
        a = inverse(a);
    }
    return power_monoid(a, n);
}

template <typename T, Integer N>
void run_operation(T min_bound, T max_bound, N min_bound2, N max_bound2, std::function<T(N,T)> base_op, std::function<T(N,T)> custom_op) {
    // Number of inputs to be considered
    const int NUM = 1000000;

    // Seed with a real random value, if available
    std::random_device r;
 
    // Choose a random mean between min_bound and max_bound
    std::default_random_engine e(r());
    std::uniform_int_distribution<T> u1dist(min_bound, max_bound);
    std::uniform_int_distribution<T> u2dist(min_bound2, max_bound2);

    // Seed the input vectors with a fresh set of random numbers
    std::vector<T> input(NUM);
    std::vector<N> input2(NUM);
    for (int i=0; i < NUM; i++) {
        input[i] = u1dist(e);
        input2[i] = u2dist(e);
    }

    // Create 2 output vectors, one for each approach
    std::vector<T> output(NUM), output2(NUM);

    // Time N operations using our custom approach
    Timer t;
    for (int i=0; i < NUM; i++) {
        output2[i] = custom_op(input2[i], input[i]);
    }
    auto custom_duration = t.clock();
    
    // Time N multiplications using the hardware implementation
    for (int i=0; i < NUM; i++) {
        output[i] = base_op(input2[i], input[i]);
    }
    auto base_duration = t.clock();
    
    std::cout << "Custom Operation took " << custom_duration << std::endl;    
    std::cout << "Base Operation took " << base_duration << std::endl;
    std::cout << "Custom runs " << ((double) base_duration.count()) / custom_duration.count() << " times faster." << std::endl;
    if ( auto failed_ind = std::invoke([&] {
        for (int i=0; i < NUM; i++) {
            if (output[i] != output2[i]) return i;
        }
        return -1;
    }) < 0 ) {
        std::cout << "Results match" << std::endl;
    } else {
        std::cout << "Results don't match " << +output[failed_ind] << " != " << +output2[failed_ind] << std::endl;
        std::cout << "Values multiplied: " << +input[failed_ind] << " x " << +input2[failed_ind] << std::endl;
    }
}

int8_t my_multiply(int a, int8_t b) {
    return static_cast<int8_t>(a) * b;
}


int main(int argc, char** argv) {

    std::cout << "Results for 64bit signed integer" << std::endl;
    run_operation<int64_t, int64_t>(-500000, 500000, -500000, 500000, std::multiplies<int64_t>(), multiply_group<int64_t, int64_t>);
    std::cout << std::endl;
    std::cout << "Results for 32bit signed integer" << std::endl;
    run_operation<int, int>(-5000, 5000, -5000, 5000, std::multiplies<int>(), multiply_group<int, int>);
    std::cout << std::endl;
    std::cout << "Results for 32bit unsigned integer" << std::endl;
    run_operation<uint32_t, uint32_t>(1000U, 10000U, 1000U, 10000U, std::multiplies<uint32_t>(), multiply_group<uint32_t, uint32_t>);

    std::cout << std::endl;
    std::cout << "Results for 8bit signed integer" << std::endl;
    run_operation<int8_t, int>(-10, 10, -10, 10, my_multiply, multiply_group<int, int8_t>);

    return 0;
}
