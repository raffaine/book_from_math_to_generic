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

template <typename T>
concept Regular = std::regular<T>;

template <typename Op, typename T>
concept BinaryOperation = Regular<T> && std::invocable<Op, T, T>;

template <typename Op, typename T>
concept SemigroupOperation = BinaryOperation<Op, T> && requires (Op op, T a, T b, T c) {
    op(a, op(b, c)) == op(op(a, b), b);
};

template <Regular A, Integer N, typename Op>
requires SemigroupOperation<Op, A>
A power_accumulate_semigroup(A r, A a, N n, Op op) {
    // A precondition for correctness is that n >= 0
    if (n == 0) return r;
    while(true) {
        if (odd(n)) {
            r = op(r, a);
            if (n == 1) return r;
        }
        n = half(n);
        a = op(a, a);
    }
}

template <Regular A, Integer N, typename Op>
requires SemigroupOperation<Op, A>
A power_semigroup(A a, N n, Op op) {
    // Precondition n > 0
    while (!odd(n)) {
        a = op(a, a);
        n = half(n);
    }
    if (n == 1) return a;
    return power_accumulate_semigroup(a, op(a, a), half(n - 1), op);
}

template <typename T>
T identity_element(std::plus<T>) {
    return T(0);
}
template <typename T>
T identity_element(std::multiplies<T>) {
    return T(1);
}

template <typename Op, typename T>
concept MonoidOperation = SemigroupOperation<Op, T> && requires (Op op, T a) {
    op(a, identity_element(op)) == a;
    op(identity_element(op), a) == a;
};

template <Regular A, Integer N, typename Op>
requires MonoidOperation<Op, A>
A power_monoid(A a, N n, Op op) {
    // Precondition n > 0
    if (n == 0) return identity_element(op);
    return power_semigroup(a, n, op);
}

template <typename T>
requires MonoidOperation<typename std::plus<T>, T>
T inverse(const T& a, std::plus<T>) {
    return -a;
}

template <typename T>
requires MonoidOperation<typename std::multiplies<T>, T>
T inverse(const T& x, std::multiplies<T>) {
    return T(1) / x;
}

template <typename Op, typename T>
concept GroupOperation = MonoidOperation<Op, T> && requires (Op op, T a) {
    op(a, inverse(a, op)) == identity_element(op);
    op(inverse(a, op), a) == identity_element(op);
};

template <Regular A, Integer N, typename Op>
requires GroupOperation<Op, A>
A power_group(A a, N n, Op op) {
    if (n < 0) {
        n = -n;
        a = inverse(a, op);
    }
    return power_monoid(a, n, op);
}

template <typename T, Integer N, typename Op, typename R = T>
void run_operation(T min_bound, T max_bound, N min_bound2, N max_bound2, std::function<T(N,T)> base_op, const int NUM = 1000000) {
    // Seed with a real random value, if available
    std::random_device r;
 
    // Choose a random mean between min_bound and max_bound
    std::default_random_engine e(r());
    std::uniform_int_distribution<R> u1dist(min_bound, max_bound);
    std::uniform_int_distribution<N> u2dist(min_bound2, max_bound2);

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
        output2[i] = power_group(input[i], input2[i], Op());
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

template <typename T, Integer N>
T my_power(T a, N n) {
    if (n < 0) {
        a = inverse(a, std::multiplies<T>());
        n = -n;
    }
    T res = a;
    while(--n > 0) {
        res = std::multiplies<T>()(res, a);
    }
    return res;
}

int main(int argc, char** argv) {

    std::cout << "Results for 64bit signed integer" << std::endl;
    run_operation<int64_t, int64_t, std::plus<int64_t>>(-500000, 500000, -500000, 500000, std::multiplies<int64_t>());
    std::cout << std::endl;
    std::cout << "Results for 32bit signed integer" << std::endl;
    run_operation<int, int, std::plus<int>>(-5000, 5000, -5000, 5000, std::multiplies<int>());
    std::cout << std::endl;
    std::cout << "Results for 32bit unsigned integer" << std::endl;
    run_operation<uint32_t, uint32_t, std::plus<uint32_t>>(1000U, 10000U, 1000U, 10000U, std::multiplies<uint32_t>());
    std::cout << std::endl;
    std::cout << "Results for 8bit signed integer" << std::endl;
    run_operation<int8_t, int, std::plus<int8_t>>(-10, 10, -10, 10, my_multiply);

    std::cout << std::endl;
    std::cout << "Results for 32bit signed float with 32bit integer exponents" << std::endl;
    run_operation<float, int, std::multiplies<float>, int>(-15, 15, -10, 10, my_power<float, int>);

    return 0;
}
